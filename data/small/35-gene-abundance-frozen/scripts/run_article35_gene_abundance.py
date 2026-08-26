#!/usr/bin/env python3
"""Run Article 35 historical and audited read-to-gene quantification branches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import subprocess
from pathlib import Path


SEED = 20260735


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class Runner:
    def __init__(self, work: Path, env_prefix: Path, threads: int, force: bool):
        self.work = work
        self.env_prefix = env_prefix
        self.threads = threads
        self.force = force
        self.logs = work / "logs"
        self.markers = work / "markers"
        self.logs.mkdir(parents=True, exist_ok=True)
        self.markers.mkdir(parents=True, exist_ok=True)
        self.records: list[dict] = []
        self.environment = os.environ.copy()
        self.environment.update(
            {
                "LC_ALL": "C",
                "LANG": "C",
                "OMP_NUM_THREADS": str(threads),
                "PYTHONHASHSEED": "0",
                "PYTHONNOUSERSITE": "1",
                "PYTHONPATH": "",
                "PATH": f"{env_prefix / 'bin'}:{self.environment.get('PATH', '')}",
            }
        )

    def binary(self, name: str) -> Path:
        path = self.env_prefix / "bin" / name
        if not path.is_file():
            raise FileNotFoundError(path)
        return path

    def _contract(self, command: list[str], inputs: list[Path]) -> dict:
        return {
            "command": command,
            "inputs": {
                str(path): {
                    "resolved": str(path.resolve()),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path.resolve()),
                }
                for path in inputs
                if path.is_file()
            },
        }

    def _cached(self, marker: Path, contract: dict, outputs: list[Path]) -> bool:
        if self.force or not marker.is_file() or not all(path.is_file() and path.stat().st_size > 0 for path in outputs):
            return False
        return json.loads(marker.read_text(encoding="utf-8")).get("contract") == contract

    def _finish(self, marker: Path, contract: dict, outputs: list[Path]) -> None:
        marker.write_text(
            json.dumps(
                {
                    "contract": contract,
                    "outputs": {str(path): {"bytes": path.stat().st_size, "sha256": sha256(path)} for path in outputs},
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    def run(
        self,
        label: str,
        command: list[str],
        outputs: list[Path],
        inputs: list[Path],
        stdout_output: Path | None = None,
    ) -> None:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "-", label)
        marker = self.markers / f"{safe}.json"
        contract = self._contract(command, inputs)
        cached = self._cached(marker, contract, outputs)
        if not cached:
            for output in outputs:
                output.parent.mkdir(parents=True, exist_ok=True)
            stdout_log = stdout_output or self.logs / f"{safe}.stdout.log"
            stderr_log = self.logs / f"{safe}.stderr.log"
            time_log = self.logs / f"{safe}.time.txt"
            timed = ["/usr/bin/time", "-v", "-o", str(time_log), *command]
            with stdout_log.open("w", encoding="utf-8") as stdout, stderr_log.open("w", encoding="utf-8") as stderr:
                subprocess.run(timed, check=True, env=self.environment, stdout=stdout, stderr=stderr)
            missing = [str(path) for path in outputs if not path.is_file() or path.stat().st_size == 0]
            if missing:
                raise RuntimeError(f"Step {label} did not create nonempty outputs: {missing}")
            self._finish(marker, contract, outputs)
        self.records.append(
            {
                "Step": label,
                "Cached": "yes" if cached else "no",
                "Command": shlex.join(command),
                "Outputs": ";".join(str(path.relative_to(self.work)) for path in outputs),
            }
        )

    def mapping_pipeline(
        self,
        label: str,
        bowtie_command: list[str],
        parser_command: list[str],
        outputs: list[Path],
        inputs: list[Path],
    ) -> None:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "-", label)
        marker = self.markers / f"{safe}.json"
        combined = [*bowtie_command, "|", *parser_command]
        contract = self._contract(combined, inputs)
        cached = self._cached(marker, contract, outputs)
        if not cached:
            for output in outputs:
                output.parent.mkdir(parents=True, exist_ok=True)
            bowtie_stderr = self.logs / f"{safe}.bowtie2.stderr.log"
            parser_stdout = self.logs / f"{safe}.parser.stdout.log"
            parser_stderr = self.logs / f"{safe}.parser.stderr.log"
            time_log = self.logs / f"{safe}.time.txt"
            timed_bowtie = ["/usr/bin/time", "-v", "-o", str(time_log), *bowtie_command]
            with bowtie_stderr.open("w", encoding="utf-8") as bt_err, parser_stdout.open(
                "w", encoding="utf-8"
            ) as parse_out, parser_stderr.open("w", encoding="utf-8") as parse_err:
                bowtie = subprocess.Popen(timed_bowtie, env=self.environment, stdout=subprocess.PIPE, stderr=bt_err, text=True)
                assert bowtie.stdout is not None
                parser = subprocess.Popen(
                    parser_command,
                    env=self.environment,
                    stdin=bowtie.stdout,
                    stdout=parse_out,
                    stderr=parse_err,
                    text=True,
                )
                bowtie.stdout.close()
                parser_status = parser.wait()
                bowtie_status = bowtie.wait()
            if bowtie_status != 0 or parser_status != 0:
                raise RuntimeError(f"Pipeline {label} failed: bowtie2={bowtie_status}, parser={parser_status}")
            missing = [str(path) for path in outputs if not path.is_file() or path.stat().st_size == 0]
            if missing:
                raise RuntimeError(f"Pipeline {label} did not create nonempty outputs: {missing}")
            self._finish(marker, contract, outputs)
        self.records.append(
            {
                "Step": label,
                "Cached": "yes" if cached else "no",
                "Command": f"{shlex.join(bowtie_command)} | {shlex.join(parser_command)}",
                "Outputs": ";".join(str(path.relative_to(self.work)) for path in outputs),
            }
        )


def conda_package_version(prefix: Path, package: str) -> str:
    records = sorted((prefix / "conda-meta").glob(f"{package}-*.json"))
    if len(records) != 1:
        raise ValueError(f"Expected one conda record for {package}, found {len(records)}")
    payload = json.loads(records[0].read_text(encoding="utf-8"))
    return f"{payload['version']} ({payload['build']})"


def tool_versions(runner: Runner) -> list[dict]:
    probes = (
        ("Bowtie2", [str(runner.binary("bowtie2")), "--version"], r"version\s+([^\s]+)"),
        ("SAMtools", [str(runner.binary("samtools")), "--version"], r"samtools\s+([^\s]+)"),
        ("HTSeq", [str(runner.binary("htseq-count")), "--version"], r"([^\s]+)"),
        ("DIAMOND", [str(runner.binary("diamond")), "version"], r"diamond version\s+([^\s]+)"),
        ("Python", [str(runner.binary("python")), "--version"], r"Python\s+([^\s]+)"),
    )
    rows: list[dict] = []
    for tool, command, pattern in probes:
        result = subprocess.run(command, check=True, env=runner.environment, capture_output=True, text=True)
        combined = result.stdout + "\n" + result.stderr
        match = re.search(pattern, combined)
        if not match:
            raise ValueError(f"Unable to parse {tool} version from {combined!r}")
        rows.append({"Tool": tool, "Version": match.group(1), "Command": shlex.join(command), "ExitCode": 0})
    rows.append(
        {
            "Tool": "seqtk",
            "Version": conda_package_version(runner.env_prefix, "seqtk"),
            "Command": "conda-meta/seqtk-*.json",
            "ExitCode": 0,
        }
    )
    return rows


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    env_prefix = args.env_prefix.resolve()
    if not (work / ".article35-inputs-complete").is_file():
        raise SystemExit("Run prepare_article35_gene_abundance_inputs.py first")
    contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))
    runner = Runner(work, env_prefix, args.threads, args.force)
    write_tsv(work / "tool-versions.tsv", tool_versions(runner))

    reference = work / "reference/megahit-mix-primary.fna"
    protein = work / "reference/megahit-mix-primary.faa"
    gff = work / "reference/megahit-mix-primary.gff"
    index_prefix = work / "index/megahit-mix-primary"
    index_outputs = [Path(f"{index_prefix}.{suffix}.bt2") for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")]
    runner.run(
        "build-bowtie2-index",
        [str(runner.binary("bowtie2-build")), "--threads", str(args.threads), str(reference), str(index_prefix)],
        index_outputs,
        [reference],
    )

    for sample in ("MOCK1", "MOCK2"):
        r1 = work / "inputs" / f"{sample}_R1.fastq.gz"
        r2 = work / "inputs" / f"{sample}_R2.fastq.gz"
        legacy_fastq = work / "legacy" / f"{sample}.R1.seed100.n10000.fastq"
        runner.run(
            f"legacy-select-{sample}",
            [str(runner.binary("seqtk")), "sample", "-s100", str(r1), "10000"],
            [legacy_fastq],
            [r1],
            stdout_output=legacy_fastq,
        )
        legacy_sam = work / "legacy" / f"{sample}.local.sam"
        runner.run(
            f"legacy-map-{sample}",
            [
                str(runner.binary("bowtie2")),
                "--local",
                "--seed",
                "100",
                "-p",
                str(args.threads),
                "-x",
                str(index_prefix),
                "-U",
                str(legacy_fastq),
                "-S",
                str(legacy_sam),
            ],
            [legacy_sam],
            [legacy_fastq, *index_outputs],
        )
        legacy_bam = work / "legacy" / f"{sample}.local.sorted.bam"
        runner.run(
            f"legacy-bam-{sample}",
            [
                str(runner.binary("samtools")),
                "sort",
                "-@",
                str(args.threads),
                "-o",
                str(legacy_bam),
                str(legacy_sam),
            ],
            [legacy_bam],
            [legacy_sam],
        )
        legacy_counts = work / "legacy" / f"{sample}.htseq-counts.tsv"
        runner.run(
            f"legacy-htseq-{sample}",
            [
                str(runner.binary("htseq-count")),
                "-f",
                "bam",
                "-r",
                "pos",
                "-t",
                "CDS",
                "-i",
                "ID",
                "-s",
                "no",
                "-a",
                "0",
                "--secondary-alignments",
                "ignore",
                "--supplementary-alignments",
                "ignore",
                str(legacy_bam),
                str(gff),
            ],
            [legacy_counts],
            [legacy_bam, gff],
            stdout_output=legacy_counts,
        )

        total_reads = int(contract["reads"][sample]["R1"]["records"]) + int(contract["reads"][sample]["R2"]["records"])
        main_counts = work / "mapping" / f"{sample}.policy-counts.tsv"
        main_summary = work / "mapping" / f"{sample}.mapping-summary.json"
        quality_hist = work / "mapping" / f"{sample}.quality-histogram.tsv"
        bowtie_command = [
            str(runner.binary("bowtie2")),
            "--very-sensitive-local",
            "-k",
            "2",
            "--seed",
            str(SEED),
            "--no-unal",
            "-p",
            str(args.threads),
            "-x",
            str(index_prefix),
            "-U",
            f"{r1},{r2}",
        ]
        parser_command = [
            str(runner.binary("python")),
            str(root / "scripts/parse_article35_sam.py"),
            "--sample",
            sample,
            "--total-reads",
            str(total_reads),
            "--counts",
            str(main_counts),
            "--summary",
            str(main_summary),
            "--quality-histogram",
            str(quality_hist),
        ]
        runner.mapping_pipeline(
            f"audited-map-{sample}",
            bowtie_command,
            parser_command,
            [main_counts, main_summary, quality_hist],
            [r1, r2, *index_outputs, root / "scripts/parse_article35_sam.py"],
        )

    annotation = contract["annotation_assets"]
    hits = work / "annotation/uniref90-top5.tsv"
    tmp = work / "tmp/diamond"
    tmp.mkdir(parents=True, exist_ok=True)
    runner.run(
        "annotate-uniref90",
        [
            str(runner.binary("diamond")),
            "blastp",
            "--db",
            annotation["uniref90_db"],
            "--query",
            str(protein),
            "--out",
            str(hits),
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qlen",
            "slen",
            "evalue",
            "bitscore",
            "qcovhsp",
            "--threads",
            str(args.threads),
            "--iterate",
            "faster",
            "--sensitive",
            "--id",
            "50",
            "--query-cover",
            "80",
            "--evalue",
            "1e-5",
            "--max-target-seqs",
            "5",
            "--masking",
            "1",
            "--tmpdir",
            str(tmp),
        ],
        [hits],
        [protein, Path(annotation["uniref90_db"])],
    )

    write_tsv(work / "command-log.tsv", runner.records)
    (work / ".article35-run-complete").write_text("completed\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "steps": len(runner.records), "samples": 2, "seed": SEED}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
