#!/usr/bin/env python3
"""Run dRep species-level dereplication plus a 99.9% ANI sensitivity branch."""

from __future__ import annotations

import argparse
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--environment", default="metagenome-drep-2026.07")
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()

    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    genomes = work / "genomes.txt"
    genome_info = work / "genome-info.csv"
    if not genomes.is_file() or not genome_info.is_file():
        raise FileNotFoundError("Run prepare_article45_drep.py first")

    env = os.environ.copy()
    env.update(
        {
            "PYTHONHASHSEED": "0",
            "MPLCONFIGDIR": str(work / "tmp" / "matplotlib"),
            "XDG_CACHE_HOME": str(work / "tmp" / "xdg"),
        }
    )
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

    commands: list[dict[str, object]] = []
    stale = work / "drep" / "strain99"
    if stale.exists():
        shutil.rmtree(stale)
    for label, ani in (("species95", "0.95"), ("nearclone999", "0.999")):
        output = work / "drep" / label
        if output.exists():
            shutil.rmtree(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        command = [
            "conda",
            "run",
            "-n",
            args.environment,
            "dRep",
            "dereplicate",
            str(output),
            "-g",
            str(genomes),
            "--genomeInfo",
            str(genome_info),
            "-comp",
            "50",
            "-con",
            "10",
            "--S_algorithm",
            "fastANI",
            "-pa",
            "0.90",
            "-sa",
            ani,
            "-nc",
            "0.30",
            "-p",
            str(args.threads),
            "--gen_warnings",
        ]
        commands.append(run_timed(f"drep-{label}", command, work, env))

    versions: list[dict[str, str]] = []
    for tool, command, role in (
        (
            "dRep",
            ["python", "-c", "import importlib.metadata as m; print(m.version('drep'))"],
            "dereplication workflow",
        ),
        ("Mash", ["mash", "--version"], "primary ANI screen"),
        ("fastANI", ["fastANI", "--version"], "secondary ANI"),
        ("Prodigal", ["prodigal", "-v"], "available dRep dependency"),
    ):
        completed = subprocess.run(
            ["conda", "run", "-n", args.environment, *command],
            text=True,
            capture_output=True,
            env=env,
            check=False,
        )
        value = (completed.stdout + " " + completed.stderr).strip().replace("\n", " ")
        versions.append({"Tool": tool, "Version": value, "Role": role})
    write_tsv(work / "command-log.tsv", commands)
    write_tsv(work / "tool-versions.tsv", versions)
    print("Article 45 dRep branches completed")


if __name__ == "__main__":
    main()
