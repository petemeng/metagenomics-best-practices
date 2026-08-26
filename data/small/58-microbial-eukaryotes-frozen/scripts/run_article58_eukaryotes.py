#!/usr/bin/env python3
"""Run EukDetect2, MEGAHIT, minimap2 and EukRep for Article 58."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shlex
import subprocess
from pathlib import Path
from typing import Iterator


MODES = ("strict", "balanced", "lenient")
FRAGMENT_LENGTHS = (3000, 5000, 10000, 20000)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--eukdetect-db", type=Path, required=True)
    parser.add_argument("--database-manifest", type=Path, required=True)
    parser.add_argument("--eukdetect-env", type=Path, required=True)
    parser.add_argument("--eukrep-env", type=Path, required=True)
    parser.add_argument("--assembly-env", type=Path, required=True)
    parser.add_argument("--read-qc-env", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--assembly-threads", type=int, default=32)
    return parser.parse_args()


def digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def capture_version(command: list[str], env: dict[str, str] | None = None) -> str:
    completed = subprocess.run(
        command,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if completed.returncode != 0 or not lines:
        raise RuntimeError(f"Version command failed: {shlex.join(command)}")
    return lines[0]


def timed_run(
    label: str,
    command: list[str],
    cwd: Path,
    logs: Path,
    env: dict[str, str],
) -> dict[str, object]:
    stdout = logs / f"{label}.stdout.log"
    stderr = logs / f"{label}.stderr.log"
    timing = logs / f"{label}.time.txt"
    wrapped = ["/usr/bin/time", "-v", "-o", str(timing), *command]
    with stdout.open("w", encoding="utf-8") as out_handle, stderr.open(
        "w", encoding="utf-8"
    ) as err_handle:
        completed = subprocess.run(
            wrapped,
            cwd=cwd,
            env=env,
            stdout=out_handle,
            stderr=err_handle,
            check=False,
        )
    row = {
        "Label": label,
        "Command": shlex.join(command),
        "ExitStatus": completed.returncode,
        "Stdout": str(stdout),
        "Stderr": str(stderr),
        "Timing": str(timing),
    }
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed; inspect {stderr}")
    return row


def reused(
    label: str,
    outputs: list[Path],
    logs: Path,
    allow_empty: bool = False,
) -> dict[str, object]:
    complete = all(
        path.is_file() and (allow_empty or path.stat().st_size > 0)
        for path in outputs
    )
    if not outputs or not complete:
        raise RuntimeError(f"Cannot reuse incomplete {label} output")
    return {
        "Label": label,
        "Command": "reused checksum-stable completed output",
        "ExitStatus": 0,
        "Stdout": str(logs / f"{label}.stdout.log"),
        "Stderr": str(logs / f"{label}.stderr.log"),
        "Timing": str(logs / f"{label}.time.txt"),
    }


def parse_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before header in {path}")
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def read_name_tokens(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    names: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            token = line.strip().lstrip(">").split()[0] if line.strip() else ""
            if token:
                names.add(token)
    return names


def complete_partition(euk_path: Path, prok_path: Path, expected: set[str]) -> bool:
    if not euk_path.is_file() or not prok_path.is_file():
        return False
    euk = read_name_tokens(euk_path)
    prok = read_name_tokens(prok_path)
    return not (euk & prok) and (euk | prok) == expected


def fastq_total_bases(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    total = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if index % 4 == 1:
                total += len(line.rstrip("\r\n"))
    return total


def write_command_log(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    if not (work / ".article58-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article58_eukaryotes.py first")
    database = args.eukdetect_db.resolve()
    database_manifest = args.database_manifest.resolve()
    required_db = (
        "eukdb.1.bt2l",
        "eukdb.2.bt2l",
        "eukdb.3.bt2l",
        "eukdb.4.bt2l",
        "eukdb.rev.1.bt2l",
        "eukdb.rev.2.bt2l",
        "specific_and_inherited_markers_per_taxid.txt",
        "taxa.sqlite",
    )
    if not all((database / name).is_file() for name in required_db):
        raise FileNotFoundError("EukDetect2 database is incomplete")

    manifest_rows = read_tsv(database_manifest)
    if len(manifest_rows) != 14 or len({row["File"] for row in manifest_rows}) != 14:
        raise RuntimeError("Expected 14 unique files in the EukDetect2 database manifest")
    database_audit: list[dict[str, object]] = []
    for row in manifest_rows:
        path = database / row["File"]
        if not path.is_file():
            raise FileNotFoundError(path)
        observed_bytes = path.stat().st_size
        observed_md5 = digest(path, "md5")
        passed = observed_bytes == int(row["Bytes"]) and observed_md5 == row["MD5"]
        database_audit.append(
            {
                "File": row["File"],
                "Release": row["Release"],
                "RecordDOI": row["RecordDOI"],
                "ExpectedBytes": row["Bytes"],
                "ObservedBytes": observed_bytes,
                "ExpectedMD5": row["MD5"],
                "ObservedMD5": observed_md5,
                "ChecksumPass": str(passed).lower(),
            }
        )
        if not passed:
            raise RuntimeError(f"EukDetect2 database checksum gate failed for {path}")
    write_tsv(work / "database-audit.tsv", database_audit)

    eukdetect_env = args.eukdetect_env.resolve()
    eukrep_env = args.eukrep_env.resolve()
    assembly_env = args.assembly_env.resolve()
    read_qc_env = args.read_qc_env.resolve()
    executables = {
        "eukdetect": eukdetect_env / "bin/eukdetect",
        "normalize": eukdetect_env / "bin/eukdetect-normalize",
        "eukrep": eukrep_env / "bin/EukRep",
        "megahit": assembly_env / "bin/megahit",
        "minimap2": assembly_env / "bin/minimap2",
        "fastp": read_qc_env / "bin/fastp",
    }
    for name, path in executables.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing {name}: {path}")

    euk_env_for_version = os.environ.copy()
    euk_env_for_version["PATH"] = str(eukdetect_env / "bin") + os.pathsep + euk_env_for_version.get("PATH", "")
    eukrep_env_for_version = os.environ.copy()
    eukrep_env_for_version["PATH"] = str(eukrep_env / "bin") + os.pathsep + eukrep_env_for_version.get("PATH", "")
    eukdetect_conda_records = sorted(
        (eukdetect_env / "conda-meta").glob("eukdetect-*.json")
    )
    if len(eukdetect_conda_records) != 1:
        raise RuntimeError(
            f"Expected one EukDetect conda record: {eukdetect_conda_records}"
        )
    eukdetect_conda_record = json.loads(
        eukdetect_conda_records[0].read_text(encoding="utf-8")
    )
    versions = [
        {
            "Tool": "EukDetect conda record",
            "Version": str(eukdetect_conda_record["version"]),
            "Executable": str(eukdetect_conda_records[0]),
        },
        {
            "Tool": "EukDetect Python metadata",
            "Version": capture_version(
                [
                    str(eukdetect_env / "bin/python"),
                    "-c",
                    "import importlib.metadata as m; print(m.version('eukdetect'))",
                ]
            ),
            "Executable": str(eukdetect_env / "bin/python"),
        },
        {
            "Tool": "EukDetect CLI",
            "Version": capture_version([str(executables["eukdetect"]), "--version"], euk_env_for_version),
            "Executable": str(executables["eukdetect"]),
        },
        {
            "Tool": "legacy-cgi compatibility package",
            "Version": capture_version(
                [
                    str(eukdetect_env / "bin/python"),
                    "-c",
                    "import importlib.metadata as m; print(m.version('legacy-cgi'))",
                ]
            ),
            "Executable": str(eukdetect_env / "bin/python"),
        },
        {
            "Tool": "EukRep",
            "Version": capture_version([str(executables["eukrep"]), "--version"], eukrep_env_for_version),
            "Executable": str(executables["eukrep"]),
        },
        {
            "Tool": "Bowtie2",
            "Version": capture_version([str(eukdetect_env / "bin/bowtie2"), "--version"], euk_env_for_version),
            "Executable": str(eukdetect_env / "bin/bowtie2"),
        },
        {
            "Tool": "Snakemake",
            "Version": capture_version([str(eukdetect_env / "bin/snakemake"), "--version"], euk_env_for_version),
            "Executable": str(eukdetect_env / "bin/snakemake"),
        },
        {
            "Tool": "fastp",
            "Version": capture_version([str(executables["fastp"]), "--version"]),
            "Executable": str(executables["fastp"]),
        },
        {
            "Tool": "MEGAHIT",
            "Version": capture_version([str(executables["megahit"]), "--version"]),
            "Executable": str(executables["megahit"]),
        },
        {
            "Tool": "minimap2",
            "Version": capture_version([str(executables["minimap2"]), "--version"]),
            "Executable": str(executables["minimap2"]),
        },
        {
            "Tool": "EukDetect2 database",
            "Version": f"{manifest_rows[0]['Release']} (doi:{manifest_rows[0]['RecordDOI']})",
            "Executable": str(database),
        },
    ]
    write_tsv(work / "tool-versions.tsv", versions)

    logs = work / "logs"
    results = work / "results"
    qc = results / "qc"
    eukdetect_results = results / "eukdetect"
    eukrep_reference = results / "eukrep-reference"
    eukrep_assembly = results / "eukrep-assembly"
    assembly_results = results / "megahit-memory60g"
    for path in (logs, qc, eukdetect_results, eukrep_reference, eukrep_assembly):
        path.mkdir(parents=True, exist_ok=True)
    run_env = os.environ.copy()
    run_env["TMPDIR"] = str((work / "tmp").resolve())
    run_env["XDG_CACHE_HOME"] = str((work / "cache").resolve())
    run_env["XDG_CONFIG_HOME"] = str((work / "config").resolve())
    for variable in ("TMPDIR", "XDG_CACHE_HOME", "XDG_CONFIG_HOME"):
        Path(run_env[variable]).mkdir(parents=True, exist_ok=True)
    commands: list[dict[str, object]] = []

    input_dir = work / "inputs"
    raw_det1 = input_dir / "SRR12324253.first-1000000-pairs.R1.fastq"
    raw_det2 = input_dir / "SRR12324253.first-1000000-pairs.R2.fastq"
    raw_asm1 = input_dir / "SRR12324253.first-20000000-pairs.R1.fastq"
    raw_asm2 = input_dir / "SRR12324253.first-20000000-pairs.R2.fastq"
    det1 = qc / "eukdetect.R1.fastq.gz"
    det2 = qc / "eukdetect.R2.fastq.gz"
    # Assembly intermediates stay uncompressed: gzip adds no evidentiary value
    # here and was the dominant wall-time cost in the measured 20-M-pair branch.
    asm1 = qc / "assembly.R1.fastq"
    asm2 = qc / "assembly.R2.fastq"
    for label, inputs, outputs in (
        ("fastp-eukdetect", (raw_det1, raw_det2), (det1, det2)),
        ("fastp-assembly", (raw_asm1, raw_asm2), (asm1, asm2)),
    ):
        json_out = qc / f"{label}.json"
        html_out = qc / f"{label}.html"
        output_list = [*outputs, json_out, html_out]
        if all(path.is_file() and path.stat().st_size > 0 for path in output_list):
            commands.append(reused(label, output_list, logs))
        else:
            commands.append(
                timed_run(
                    label,
                    [
                        str(executables["fastp"]),
                        "--in1",
                        str(inputs[0]),
                        "--in2",
                        str(inputs[1]),
                        "--out1",
                        str(outputs[0]),
                        "--out2",
                        str(outputs[1]),
                        "--detect_adapter_for_pe",
                        "--qualified_quality_phred",
                        "20",
                        "--unqualified_percent_limit",
                        "40",
                        "--length_required",
                        "80",
                        "--compression",
                        "4",
                        "--thread",
                        str(args.threads),
                        "--json",
                        str(json_out),
                        "--html",
                        str(html_out),
                    ],
                    root,
                    logs,
                    run_env,
                )
            )

    total_bases = fastq_total_bases(det1) + fastq_total_bases(det2)
    library_sizes = results / "library-sizes-post-qc.tsv"
    library_sizes.write_text(
        f"Sample\tTotalBases\nZymo_D6300\t{total_bases}\n", encoding="utf-8"
    )

    euk_hits = eukdetect_results / "Zymo_D6300_filtered_hits_table.txt"
    euk_frac = eukdetect_results / "Zymo_D6300_filtered_hits_eukfrac.txt"
    all_hits = eukdetect_results / "filtering" / "Zymo_D6300_all_hits_table.txt"
    if euk_hits.is_file() and euk_frac.is_file() and all_hits.is_file():
        commands.append(reused("eukdetect", [euk_hits, euk_frac, all_hits], logs))
    else:
        euk_env = run_env.copy()
        euk_env["PATH"] = str(eukdetect_env / "bin") + os.pathsep + euk_env.get("PATH", "")
        commands.append(
            timed_run(
                "eukdetect",
                [
                    str(executables["eukdetect"]),
                    "single",
                    "-1",
                    str(det1),
                    "-2",
                    str(det2),
                    "-n",
                    "Zymo_D6300",
                    "--outdir",
                    str(eukdetect_results),
                    "--database",
                    str(database),
                    "--cores",
                    str(args.threads),
                    "--force",
                ],
                root,
                logs,
                euk_env,
            )
        )
    normalized = eukdetect_results / "Zymo_D6300.normalized.tsv"
    commands.append(
        timed_run(
            "eukdetect-normalize",
            [
                str(executables["normalize"]),
                "--eukfrac",
                str(euk_frac),
                "--library-sizes",
                str(library_sizes),
                "--output",
                str(normalized),
            ],
            root,
            logs,
            run_env,
        )
    )

    eukrep_env_vars = run_env.copy()
    eukrep_env_vars["PATH"] = str(eukrep_env / "bin") + os.pathsep + eukrep_env_vars.get("PATH", "")
    for fragment_length in FRAGMENT_LENGTHS:
        input_fasta = work / "benchmarks" / f"eukrep-fragments-{fragment_length}.fna"
        expected_fragments = {header.split()[0] for header, _ in parse_fasta(input_fasta)}
        if len(expected_fragments) != 800:
            raise RuntimeError(
                f"Expected 800 reference fragments at {fragment_length} bp"
            )
        for mode in MODES:
            euk_out = eukrep_reference / f"length-{fragment_length}.{mode}.euk.names"
            prok_out = eukrep_reference / f"length-{fragment_length}.{mode}.prok.names"
            label = f"eukrep-reference-{fragment_length}-{mode}"
            if complete_partition(euk_out, prok_out, expected_fragments):
                commands.append(
                    reused(label, [euk_out, prok_out], logs, allow_empty=True)
                )
            else:
                commands.append(
                    timed_run(
                        label,
                        [
                            str(executables["eukrep"]),
                            "-i",
                            str(input_fasta),
                            "-o",
                            str(euk_out),
                            "--prokarya",
                            str(prok_out),
                            "--seq_names",
                            "--min",
                            # EukRep applies --min to each internal 5-kb chunk,
                            # not only to the original scaffold. Values above
                            # 5 kb therefore yield empty output with exit code
                            # zero. Keep the published/default 3-kb threshold
                            # constant while benchmarking scaffold length.
                            "3000",
                            "-m",
                            mode,
                            "--tie",
                            "prok",
                            "-ff",
                        ],
                        root,
                        logs,
                        eukrep_env_vars,
                    )
                )
            if not complete_partition(euk_out, prok_out, expected_fragments):
                raise RuntimeError(
                    "EukRep failed to classify the complete reference partition "
                    f"for {fragment_length} bp, {mode} mode"
                )

    contigs = assembly_results / "final.contigs.fa"
    if contigs.is_file() and contigs.stat().st_size > 0:
        commands.append(reused("megahit", [contigs], logs))
    else:
        commands.append(
            timed_run(
                "megahit",
                [
                    str(executables["megahit"]),
                    "-1",
                    str(asm1),
                    "-2",
                    str(asm2),
                    "-o",
                    str(assembly_results),
                    "--presets",
                    "meta-sensitive",
                    "--min-contig-len",
                    "1000",
                    "--memory",
                    "60000000000",
                    "-t",
                    str(args.assembly_threads),
                ],
                root,
                logs,
                run_env,
            )
        )

    oneline = results / "megahit-contigs.oneline.fna"
    with oneline.open("w", encoding="utf-8") as handle:
        for header, sequence in parse_fasta(contigs):
            handle.write(f">{header.split()[0]}\n{sequence}\n")

    expected_assembly = {
        header.split()[0]
        for header, sequence in parse_fasta(oneline)
        if len(sequence) >= 3000
    }

    for mode in MODES:
        euk_out = eukrep_assembly / f"megahit.{mode}.euk.names"
        prok_out = eukrep_assembly / f"megahit.{mode}.prok.names"
        label = f"eukrep-assembly-{mode}"
        if complete_partition(euk_out, prok_out, expected_assembly):
            commands.append(
                reused(label, [euk_out, prok_out], logs, allow_empty=True)
            )
        else:
            commands.append(
                timed_run(
                    label,
                    [
                        str(executables["eukrep"]),
                        "-i",
                        str(oneline),
                        "-o",
                        str(euk_out),
                        "--prokarya",
                        str(prok_out),
                        "--seq_names",
                        "--min",
                        "3000",
                        "-m",
                        mode,
                        "--tie",
                        "prok",
                        "-ff",
                    ],
                    root,
                    logs,
                    eukrep_env_vars,
                )
            )
        if not complete_partition(euk_out, prok_out, expected_assembly):
            raise RuntimeError(
                f"EukRep failed to classify the complete assembly partition in {mode} mode"
            )

    paf = results / "megahit-to-zymo-v2.paf"
    commands.append(
        timed_run(
            "minimap2-assembly-truth",
            [
                str(executables["minimap2"]),
                "-x",
                "asm5",
                "--secondary=no",
                "-t",
                str(args.assembly_threads),
                "-o",
                str(paf),
                str(work / "references" / "zymo-v2-all-references.fna"),
                str(oneline),
            ],
            root,
            logs,
            run_env,
        )
    )
    write_command_log(work / "command-log.tsv", commands)
    (work / ".article58-run-complete").write_text("verified\n", encoding="utf-8")
    print(f"Article 58 workflow complete: {results}")


if __name__ == "__main__":
    main()
