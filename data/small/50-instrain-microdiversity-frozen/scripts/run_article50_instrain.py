#!/usr/bin/env python3
"""Run explicit inStrain profiles and a two-sample popANI comparison for Article 50."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path

from article41_44_utils import parse_time, write_tsv


def run_timed(
    label: str, command: list[str], work: Path, env: dict[str, str]
) -> dict[str, object]:
    stdout_path = work / "logs" / f"{label}.stdout.log"
    stderr_path = work / "logs" / f"{label}.stderr.log"
    time_path = work / "logs" / f"{label}.time.txt"
    timed = ["/usr/bin/time", "-v", "-o", str(time_path), *command]
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        completed = subprocess.run(
            timed, stdout=stdout, stderr=stderr, env=env, check=False
        )
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
        text=True,
        capture_output=True,
        env=env,
        check=False,
    )
    lines = [
        line.strip()
        for line in (completed.stdout + "\n" + completed.stderr).splitlines()
        if line.strip() and not line.startswith("Matplotlib created")
    ]
    return lines[0] if lines else "VERSION_NOT_REPORTED"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--environment", default="metagenome-instrain-2026.07")
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article50-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article50_instrain.py first")
    contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))

    for path in (work / "profiles", work / "comparison"):
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)
    for sample in contract["samples"]:
        index = work / "inputs" / f"{sample}.identity95.unfiltered.bam.bai"
        if index.exists():
            index.unlink()

    env = os.environ.copy()
    cache = work / "matplotlib-cache"
    temporary = work / "tmp"
    cache.mkdir(exist_ok=True)
    temporary.mkdir(exist_ok=True)
    env.update(
        {
            "MPLCONFIGDIR": str(cache),
            "TMPDIR": str(temporary),
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
    )
    commands: list[dict[str, object]] = []
    for sample in contract["samples"]:
        bam = work / "inputs" / f"{sample}.identity95.unfiltered.bam"
        index = work / "inputs" / f"{sample}.identity95.unfiltered.bam.bai"
        command = [
            "conda",
            "run",
            "-n",
            args.environment,
            "samtools",
            "index",
            "-@",
            str(args.threads),
            "-o",
            str(index),
            str(bam),
        ]
        commands.append(run_timed(f"samtools-index-{sample.lower()}", command, work, env))

    catalog = work / "inputs/article50-sgb-catalog.fna"
    stb = work / "inputs/article50-scaffold-to-bin.tsv"
    for sample in contract["samples"]:
        bam = work / "inputs" / f"{sample}.identity95.unfiltered.bam"
        output = work / "profiles" / sample
        command = [
            "conda",
            "run",
            "-n",
            args.environment,
            "inStrain",
            "profile",
            str(bam),
            str(catalog),
            "-o",
            str(output),
            "-p",
            str(args.threads),
            "-s",
            str(stb),
            "--min_read_ani",
            str(contract["min_read_ani"]),
            "--min_mapq",
            str(contract["min_mapq"]),
            "--pairing_filter",
            contract["pairing_filter"],
            "--min_cov",
            str(contract["min_cov"]),
            "--min_freq",
            str(contract["min_freq"]),
            "--fdr",
            str(contract["fdr"]),
            "--rarefied_coverage",
            str(contract["rarefied_coverage"]),
            "--skip_mm_profiling",
            "--skip_plot_generation",
            "--force_compress",
        ]
        commands.append(run_timed(f"instrain-profile-{sample.lower()}", command, work, env))

    compare = work / "comparison/MOCK1-vs-MOCK2"
    command = [
        "conda",
        "run",
        "-n",
        args.environment,
        "inStrain",
        "compare",
        "-i",
        str(work / "profiles/MOCK1"),
        str(work / "profiles/MOCK2"),
        "-o",
        str(compare),
        "-p",
        str(args.threads),
        "-s",
        str(stb),
        "--min_cov",
        str(contract["min_cov"]),
        "--min_freq",
        str(contract["min_freq"]),
        "--fdr",
        str(contract["fdr"]),
        "--database_mode",
        "--breadth",
        str(contract["compare_presence_breadth"]),
        "--ani_threshold",
        "0.99999",
        "--coverage_treshold",
        "0.5",
        "--store_mismatch_locations",
        "--skip_plot_generation",
        "--force_compress",
    ]
    commands.append(run_timed("instrain-compare", command, work, env))
    write_tsv(work / "command-log.tsv", commands)
    write_tsv(
        work / "resource-summary.tsv",
        [parse_time(Path(row["TimeLog"])) for row in commands],
    )
    write_tsv(
        work / "tool-versions.tsv",
        [
            {
                "Tool": "inStrain",
                "Version": version(args.environment, ["inStrain", "profile", "--version"], env),
                "Role": "SNV, nucleotide-diversity and popANI profiling",
            },
            {
                "Tool": "SAMtools",
                "Version": version(args.environment, ["samtools", "--version"], env),
                "Role": "BAM indexing and input audit",
            },
            {
                "Tool": "Python",
                "Version": version(args.environment, ["python", "--version"], env),
                "Role": "inStrain runtime",
            },
        ],
    )
    (work / ".article50-run-complete").write_text("complete\n", encoding="utf-8")
    print("Article 50 inStrain profiles and comparison completed")


if __name__ == "__main__":
    main()
