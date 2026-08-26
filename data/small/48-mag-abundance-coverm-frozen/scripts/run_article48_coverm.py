#!/usr/bin/env python3
"""Run CoverM with a primary and stricter read-identity branch."""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
from pathlib import Path

from article41_44_utils import read_tsv, write_tsv


def run_timed(label: str, command: list[str], work: Path, env: dict[str, str]) -> dict[str, object]:
    logs = work / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    stdout_path = logs / f"{label}.stdout.log"
    stderr_path = logs / f"{label}.stderr.log"
    time_path = logs / f"{label}.time.txt"
    timed = ["/usr/bin/time", "-v", "-o", str(time_path), *command]
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(timed, stdout=stdout, stderr=stderr, env=env, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed; inspect {stderr_path}")
    return {
        "Label": label,
        "ExitStatus": completed.returncode,
        "Command": shlex.join(command),
        "Stdout": str(stdout_path),
        "Stderr": str(stderr_path),
        "TimeLog": str(time_path),
    }


def version(environment: str, command: list[str], env: dict[str, str]) -> str:
    completed = subprocess.run(
        ["conda", "run", "-n", environment, *command],
        text=True, capture_output=True, env=env, check=False,
    )
    return (completed.stdout + " " + completed.stderr).strip().replace("\n", " ")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--environment", default="metagenome-assembly-2026.07")
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article48-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article48_coverm.py first")
    samples = read_tsv(work / "samples.tsv")
    if [row["Sample"] for row in samples] != ["MOCK1", "MOCK2"]:
        raise ValueError("Article 48 sample order changed")

    env = os.environ.copy()
    env.update({"TMPDIR": str(work / "tmp"), "RUST_BACKTRACE": "1"})
    Path(env["TMPDIR"]).mkdir(parents=True, exist_ok=True)
    commands: list[dict[str, object]] = []
    raw = work / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    for branch, identity in (("identity95", "95"), ("identity97", "97")):
        bam_dir = work / "bam" / branch
        if bam_dir.exists():
            shutil.rmtree(bam_dir)
        bam_dir.parent.mkdir(parents=True, exist_ok=True)
        output = raw / f"coverm-{branch}.tsv"
        if output.exists():
            output.unlink()
        command = [
            "conda", "run", "-n", args.environment,
            "coverm", "genome",
            "--coupled",
            *[path for row in samples for path in (row["R1"], row["R2"])],
            "--genome-fasta-directory", str(work / "inputs/genomes"),
            "--genome-fasta-extension", "fna",
            "--mapper", "strobealign",
            "--min-read-percent-identity", identity,
            "--min-read-aligned-percent", "75",
            "--proper-pairs-only",
            "--exclude-supplementary",
            "--min-covered-fraction", "0",
            "--contig-end-exclusion", "75",
            "--trim-min", "5", "--trim-max", "95",
            "--methods", "mean", "trimmed_mean", "covered_fraction",
            "relative_abundance", "count", "anir", "length",
            "--output-format", "sparse",
            "--cache-unfiltered-bam-directory", str(bam_dir),
            "--threads", str(args.threads),
            "--output-file", str(output),
        ]
        commands.append(run_timed(f"coverm-{branch}", command, work, env))

    tools = [
        {"Tool": "CoverM", "Version": version(args.environment, ["coverm", "--version"], env), "Role": "mapping filters and genome-level metrics"},
        {"Tool": "Strobealign", "Version": version(args.environment, ["strobealign", "--version"], env), "Role": "short-read mapper"},
        {"Tool": "SAMtools", "Version": version(args.environment, ["samtools", "--version"], env).split(" ")[1], "Role": "BAM backend"},
    ]
    write_tsv(work / "command-log.tsv", commands)
    write_tsv(work / "tool-versions.tsv", tools)
    (work / ".article48-run-complete").write_text("complete\n", encoding="utf-8")
    print("Article 48 CoverM branches completed")


if __name__ == "__main__":
    main()
