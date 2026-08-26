#!/usr/bin/env python3
"""Run Bowtie2/SAMtools/JGI depth workflow for Article 41."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--assembly-env", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=24)
    return parser.parse_args()


def q(value: object) -> str:
    return shlex.quote(str(value))


def run_timed(
    label: str,
    command: list[str],
    work: Path,
    runtime: dict[str, str],
    rows: list[dict[str, object]],
    *,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
) -> None:
    resource = work / "logs" / f"{label}.time.txt"
    stdout_path = stdout_path or work / "logs" / f"{label}.stdout.log"
    stderr_path = stderr_path or work / "logs" / f"{label}.stderr.log"
    started = datetime.now(timezone.utc)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        result = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(resource), *command],
            env=runtime,
            stdout=stdout,
            stderr=stderr,
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
            "ResourceLog": str(resource),
            "Stdout": str(stdout_path),
            "Stderr": str(stderr_path),
        }
    )
    if result.returncode != 0:
        tail = stderr_path.read_text(encoding="utf-8", errors="replace")[-5000:]
        raise RuntimeError(f"{label} failed ({result.returncode})\n{tail}")


def record_reuse(
    label: str,
    command: list[str],
    work: Path,
    rows: list[dict[str, object]],
    *,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
) -> None:
    rows.append(
        {
            "Label": label,
            "Command": shlex.join(command),
            "ReturnCode": 0,
            "StartedUTC": "reused-complete-output",
            "EndedUTC": "reused-complete-output",
            "ResourceLog": str(work / "logs" / f"{label}.time.txt"),
            "Stdout": str(stdout_path or work / "logs" / f"{label}.stdout.log"),
            "Stderr": str(stderr_path or work / "logs" / f"{label}.stderr.log"),
        }
    )


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    prefix = args.assembly_env.resolve()
    if not (work / ".article41-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article41_mapping_depth.py first")
    if (work / ".article41-run-complete").is_file():
        print(f"Article 41 mapping already complete: {work}")
        return 0
    runtime = os.environ.copy()
    runtime.update(
        {
            "PATH": f"{prefix / 'bin'}:/usr/bin:/bin",
            "LC_ALL": "C",
            "TZ": "UTC",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "TMPDIR": str(work / "tmp"),
            "OMP_NUM_THREADS": str(args.threads),
        }
    )
    bowtie2 = prefix / "bin/bowtie2"
    bowtie2_build = prefix / "bin/bowtie2-build"
    samtools = prefix / "bin/samtools"
    jgi = prefix / "bin/jgi_summarize_bam_contig_depths"
    for path in (bowtie2, bowtie2_build, samtools, jgi):
        if not path.is_file():
            raise FileNotFoundError(path)
    assembly = work / "inputs/megahit-coassembly.ge1000.fna"
    index = work / "index/megahit-coassembly.ge1000"
    rows: list[dict[str, object]] = []

    index_command = [str(bowtie2_build), "--threads", str(args.threads), str(assembly), str(index)]
    index_files = list((work / "index").glob("megahit-coassembly.ge1000*.bt2"))
    if len(index_files) == 6 and all(path.stat().st_size > 0 for path in index_files) and (work / "logs/bowtie2-build.time.txt").stat().st_size > 0:
        record_reuse("bowtie2-build", index_command, work, rows)
    else:
        run_timed("bowtie2-build", index_command, work, runtime, rows)
    faidx_command = [str(samtools), "faidx", str(assembly)]
    if (assembly.with_suffix(assembly.suffix + ".fai")).is_file() and (work / "logs/samtools-faidx.time.txt").stat().st_size > 0:
        record_reuse("samtools-faidx", faidx_command, work, rows)
    else:
        run_timed("samtools-faidx", faidx_command, work, runtime, rows)

    sample_rows = list(csv.DictReader((work / "inputs/samples.tsv").open(encoding="utf-8"), delimiter="\t"))
    bams = []
    for row in sample_rows:
        sample = row["Sample"]
        bam = work / f"bam/{sample}.primary.sorted.bam"
        bams.append(bam)
        bowtie_log = work / f"logs/map-{sample}.bowtie2.log"
        pipeline_log = work / f"logs/map-{sample}.pipeline.stderr.log"
        pipeline = (
            f"{q(bowtie2)} --very-sensitive --seed 20260741 -p 16 "
            f"-x {q(index)} -1 {q(row['R1'])} -2 {q(row['R2'])} "
            f"2> {q(bowtie_log)} | "
            f"{q(samtools)} view -@ 4 -b -F 3588 -q 0 - | "
            f"{q(samtools)} sort -@ 8 -m 2G -T {q(work / f'tmp/{sample}.sort')} -o {q(bam)} -"
        )
        run_timed(
            f"map-{sample}",
            ["/bin/bash", "-o", "pipefail", "-c", pipeline],
            work,
            runtime,
            rows,
            stderr_path=pipeline_log,
        )
        subprocess.run([str(samtools), "quickcheck", "-v", str(bam)], check=True, env=runtime)
        run_timed(
            f"index-{sample}",
            [str(samtools), "index", "-@", "8", str(bam)],
            work,
            runtime,
            rows,
        )
        for operation in ("flagstat", "idxstats", "stats", "coverage"):
            output = work / f"depth/{sample}.{operation}.tsv"
            command = [str(samtools), operation, "-@", "8", str(bam)]
            if operation == "coverage":
                command = [str(samtools), "coverage", "-o", str(output), str(bam)]
                run_timed(f"coverage-{sample}", command, work, runtime, rows)
            else:
                run_timed(
                    f"{operation}-{sample}",
                    command,
                    work,
                    runtime,
                    rows,
                    stdout_path=output,
                )

    run_timed(
        "jgi-depth",
        [
            str(jgi),
            "--outputDepth",
            str(work / "depth/jgi-depth.tsv"),
            "--pairedContigs",
            str(work / "depth/paired-contigs.tsv"),
            "--percentIdentity",
            "97",
            "--minMapQual",
            "0",
            "--referenceFasta",
            str(assembly),
            *map(str, bams),
        ],
        work,
        runtime,
        rows,
    )
    with (work / "command-log.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "article": 41,
        "commands": len(rows),
        "samples": [row["Sample"] for row in sample_rows],
        "threads": args.threads,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    (work / ".article41-run-complete").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
