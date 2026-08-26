#!/usr/bin/env python3
"""Run fail-closed GTDB-Tk R232 installation audit and classify_wf."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path

from article41_44_utils import write_tsv


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
        ["conda", "run", "-n", environment, *command], text=True,
        capture_output=True, env=env, check=False,
    )
    return (completed.stdout + " " + completed.stderr).strip().replace("\n", " ")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--environment", default="metagenome-gtdbtk-2026.07")
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article46-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article46_gtdbtk.py first")
    contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))
    database = Path(contract["database_path"])
    env = os.environ.copy()
    env.update({
        "GTDBTK_DATA_PATH": str(database),
        "TMPDIR": str(work / "tmp"),
        "PYTHONHASHSEED": "0",
    })
    Path(env["TMPDIR"]).mkdir(parents=True, exist_ok=True)
    output = work / "gtdbtk"
    if output.exists():
        shutil.rmtree(output)
    commands = []
    commands.append(run_timed(
        "gtdbtk-check-install",
        ["conda", "run", "-n", args.environment, "gtdbtk", "check_install"],
        work, env,
    ))
    commands.append(run_timed(
        "gtdbtk-classify-wf",
        [
            "conda", "run", "-n", args.environment,
            "gtdbtk", "classify_wf",
            "--genome_dir", str(work / "inputs/genomes"),
            "--out_dir", str(output),
            "--extension", "fna",
            "--prefix", "article46",
            "--cpus", str(args.threads),
            "--pplacer_cpus", "1",
            "--min_perc_aa", "10",
            "--min_af", "0.5",
            "--write_single_copy_genes",
            "--keep_intermediates",
            "--tmpdir", str(work / "tmp"),
        ],
        work, env,
    ))
    tools = [
        {"Tool": "GTDB-Tk", "Version": version(args.environment, ["gtdbtk", "--version"], env), "Role": "taxonomy workflow"},
        {"Tool": "Prodigal", "Version": version(args.environment, ["prodigal", "-v"], env), "Role": "gene prediction"},
        {"Tool": "pplacer", "Version": version(args.environment, ["pplacer", "--version"], env), "Role": "marker-tree placement"},
        {"Tool": "skani", "Version": version(args.environment, ["skani", "--version"], env), "Role": "ANI screening and classification"},
    ]
    write_tsv(work / "command-log.tsv", commands)
    write_tsv(work / "tool-versions.tsv", tools)
    (work / ".article46-run-complete").write_text("complete\n", encoding="utf-8")
    print("Article 46 GTDB-Tk classification completed")


if __name__ == "__main__":
    main()
