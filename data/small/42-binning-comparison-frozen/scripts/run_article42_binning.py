#!/usr/bin/env python3
"""Run the five deterministic Article 42 binner branches."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--assembly-env", type=Path, required=True)
    parser.add_argument("--binning-env", type=Path, required=True)
    parser.add_argument("--kraken-db", type=Path, required=True)
    parser.add_argument("--article41-work-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=24)
    return parser.parse_args()


def run_timed(
    label: str,
    command: list[str],
    work: Path,
    runtime: dict[str, str],
    rows: list[dict[str, object]],
) -> None:
    stdout_path = work / "logs" / f"{label}.stdout.log"
    stderr_path = work / "logs" / f"{label}.stderr.log"
    resource_path = work / "logs" / f"{label}.time.txt"
    started = datetime.now(timezone.utc)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        result = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(resource_path), *command],
            stdout=stdout,
            stderr=stderr,
            env=runtime,
            text=True,
        )
    ended = datetime.now(timezone.utc)
    rows.append(
        {
            "Label": label,
            "Command": shlex.join(command),
            "ReturnCode": result.returncode,
            "StartedUTC": started.isoformat(),
            "EndedUTC": ended.isoformat(),
            "Stdout": str(stdout_path),
            "Stderr": str(stderr_path),
            "ResourceLog": str(resource_path),
        }
    )
    if result.returncode:
        tail = stderr_path.read_text(encoding="utf-8", errors="replace")[-8000:]
        raise RuntimeError(f"{label} failed ({result.returncode})\n{tail}")


def record_reuse(label: str, command: list[str], work: Path, rows: list[dict[str, object]]) -> None:
    rows.append(
        {
            "Label": label,
            "Command": shlex.join(command),
            "ReturnCode": 0,
            "StartedUTC": "reused-complete-output",
            "EndedUTC": "reused-complete-output",
            "Stdout": str(work / "logs" / f"{label}.stdout.log"),
            "Stderr": str(work / "logs" / f"{label}.stderr.log"),
            "ResourceLog": str(work / "logs" / f"{label}.time.txt"),
        }
    )


def count_bins(path: Path) -> int:
    return sum(
        1
        for candidate in path.rglob("*")
        if candidate.is_file() and candidate.suffix in {".fa", ".fna", ".fasta"} and candidate.stat().st_size > 0
    )


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    assembly_env = args.assembly_env.resolve()
    binning_env = args.binning_env.resolve()
    kraken_db = args.kraken_db.resolve()
    article41_work = args.article41_work_dir.resolve()
    if not (work / ".article42-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article42_binning.py first")
    if (work / ".article42-run-complete").is_file():
        print(f"Article 42 binning already complete: {work}")
        return 0
    runtime = os.environ.copy()
    runtime.update(
        {
            "PATH": f"{binning_env / 'bin'}:{assembly_env / 'bin'}:/usr/bin:/bin",
            "LC_ALL": "C",
            "TZ": "UTC",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "PYTHONPATH": "",
            "TMPDIR": str(work / "tmp"),
            "OMP_NUM_THREADS": str(args.threads),
            "MKL_NUM_THREADS": str(args.threads),
            "OPENBLAS_NUM_THREADS": str(args.threads),
        }
    )
    fasta = work / "inputs/megahit-coassembly.ge1500.fna"
    depth_multi = work / "inputs/jgi-depth.ge1500.tsv"
    depth_single = work / "inputs/jgi-depth.MOCK1-only.tsv"
    abundance = work / "inputs/vamb-abundance.tsv"
    metabat2 = assembly_env / "bin/metabat2"
    semibin = binning_env / "bin/SemiBin2"
    vamb = binning_env / "bin/vamb"
    kraken2 = binning_env / "bin/kraken2"
    python = binning_env / "bin/python"
    rows: list[dict[str, object]] = []

    kraken_output = work / "taxonomy/contigs.kraken.tsv"
    kraken_command = [
            str(kraken2),
            "--db",
            str(kraken_db),
            "--threads",
            str(args.threads),
            "--confidence",
            "0.05",
            "--minimum-hit-groups",
            "2",
            "--report",
            str(work / "taxonomy/contigs.kreport.tsv"),
            "--output",
            str(kraken_output),
            str(fasta),
        ]
    if kraken_output.is_file() and (work / "taxonomy/contigs.kreport.tsv").is_file():
        record_reuse("kraken2-contig-taxonomy", kraken_command, work, rows)
    else:
        run_timed("kraken2-contig-taxonomy", kraken_command, work, runtime, rows)
    convert_command = [
            str(python),
            str(root / "scripts/convert_article42_kraken_taxonomy.py"),
            "--kraken-output",
            str(kraken_output),
            "--nodes",
            str(kraken_db / "nodes.dmp"),
            "--names",
            str(kraken_db / "names.dmp"),
            "--output",
            str(work / "taxonomy/taxvamb-taxonomy.tsv"),
            "--summary",
            str(work / "taxonomy/taxonomy-summary.tsv"),
        ]
    if (work / "taxonomy/taxvamb-taxonomy.tsv").is_file() and (work / "taxonomy/taxonomy-summary.tsv").is_file():
        record_reuse("convert-taxvamb-taxonomy", convert_command, work, rows)
    else:
        run_timed("convert-taxvamb-taxonomy", convert_command, work, runtime, rows)

    metabat_specs = (
        ("metabat2-MOCK1-only", depth_single),
        ("metabat2-multisample", depth_multi),
    )
    for label, depth in metabat_specs:
        output_dir = work / "bins" / label
        command = [
                str(metabat2),
                "-i",
                str(fasta),
                "-a",
                str(depth),
                "-o",
                str(output_dir / "bin"),
                "-m",
                "1500",
                "-s",
                "200000",
                "-t",
                str(args.threads),
                "--seed",
                "20260742",
                "--saveCls",
            ]
        if count_bins(output_dir) > 0:
            record_reuse(label, command, work, rows)
        else:
            output_dir.mkdir(exist_ok=True)
            run_timed(label, command, work, runtime, rows)
        if count_bins(output_dir) == 0:
            raise RuntimeError(f"{label} produced no bins")

    semibin_out = work / "bins/semibin2-self-supervised"
    bam_paths = [
        article41_work / "bam/MOCK1.primary.sorted.bam",
        article41_work / "bam/MOCK2.primary.sorted.bam",
    ]
    for bam in bam_paths:
        if not bam.is_file() or not bam.with_suffix(bam.suffix + ".bai").is_file():
            raise FileNotFoundError(f"Article 41 BAM/index is missing: {bam}")
    semibin_command = [
            str(semibin),
            "single_easy_bin",
            "--input-fasta",
            str(fasta),
            "--input-bam",
            *map(str, bam_paths),
            "--output",
            str(semibin_out),
            "--threads",
            str(args.threads),
            "--min-len",
            "1500",
            "--minfasta-kbs",
            "200",
            "--self-supervised",
            "--engine",
            "cpu",
            "--random-seed",
            "20260742",
            "--compression",
            "none",
            "--tag-output",
            "SemiBin2",
        ]
    if count_bins(semibin_out / "output_bins") > 0:
        record_reuse("semibin2-self-supervised", semibin_command, work, rows)
    else:
        if semibin_out.exists():
            shutil.rmtree(semibin_out)
        run_timed("semibin2-self-supervised", semibin_command, work, runtime, rows)
    if count_bins(semibin_out / "output_bins") == 0:
        raise RuntimeError("SemiBin2 produced no bins")

    vamb_out = work / "bins/vamb-taxonomy-free"
    vamb_command = [
            str(vamb),
            "bin",
            "default",
            "--outdir",
            str(vamb_out),
            "--fasta",
            str(fasta),
            "--abundance_tsv",
            str(abundance),
            "-m",
            "1500",
            "-p",
            str(args.threads),
            "--seed",
            "20260742",
            "--minfasta",
            "200000",
            "-o",
        ]
    if count_bins(vamb_out / "bins") > 0:
        record_reuse("vamb-taxonomy-free", vamb_command, work, rows)
    else:
        if vamb_out.exists():
            shutil.rmtree(vamb_out)
        run_timed("vamb-taxonomy-free", vamb_command, work, runtime, rows)
    if count_bins(vamb_out / "bins") == 0:
        raise RuntimeError("VAMB produced no bins")

    taxvamb_out = work / "bins/taxvamb-kraken2"
    taxvamb_command = [
            str(vamb),
            "bin",
            "taxvamb",
            "--outdir",
            str(taxvamb_out),
            "--fasta",
            str(fasta),
            "--abundance_tsv",
            str(abundance),
            "--taxonomy",
            str(work / "taxonomy/taxvamb-taxonomy.tsv"),
            "-m",
            "1500",
            "-p",
            str(args.threads),
            "--seed",
            "20260742",
            "--minfasta",
            "200000",
            "-o",
        ]
    if count_bins(taxvamb_out / "bins") > 0:
        record_reuse("taxvamb-kraken2", taxvamb_command, work, rows)
    else:
        if taxvamb_out.exists():
            shutil.rmtree(taxvamb_out)
        run_timed("taxvamb-kraken2", taxvamb_command, work, runtime, rows)
    if count_bins(taxvamb_out / "bins") == 0:
        raise RuntimeError("TaxVAMB produced no bins")

    with (work / "command-log.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "article": 42,
        "commands": len(rows),
        "threads": args.threads,
        "bin_counts": {
            "MetaBAT2-MOCK1-only": count_bins(work / "bins/metabat2-MOCK1-only"),
            "MetaBAT2-multisample": count_bins(work / "bins/metabat2-multisample"),
            "SemiBin2-self-supervised": count_bins(semibin_out / "output_bins"),
            "VAMB-taxonomy-free": count_bins(vamb_out / "bins"),
            "TaxVAMB-Kraken2": count_bins(taxvamb_out / "bins"),
        },
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    (work / ".article42-run-complete").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
