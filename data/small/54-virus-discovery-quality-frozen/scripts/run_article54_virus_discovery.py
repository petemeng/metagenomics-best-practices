#!/usr/bin/env python3
"""Run the checksum-locked Article 54 virus discovery benchmark."""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
from pathlib import Path

from article41_44_utils import write_tsv


def timed_run(
    *, label: str, command: list[str], cwd: Path, logs: Path, env: dict[str, str]
) -> dict[str, object]:
    stdout = logs / f"{label}.stdout.log"
    stderr = logs / f"{label}.stderr.log"
    timing = logs / f"{label}.time.txt"
    wrapped = ["/usr/bin/time", "-v", "-o", str(timing), *command]
    with stdout.open("w", encoding="utf-8") as out_handle:
        with stderr.open("w", encoding="utf-8") as err_handle:
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--genomad-db", type=Path, required=True)
    parser.add_argument("--checkv-db", type=Path, required=True)
    parser.add_argument("--virsorter2-db", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    results = args.results_dir.resolve()
    input_fna = work / "input/checkv-test-sequences.fna"
    if not input_fna.is_file():
        raise FileNotFoundError(
            "Run prepare_article54_virus_discovery.py before this command"
        )
    for database in (args.genomad_db, args.checkv_db, args.virsorter2_db):
        if not database.resolve().is_dir():
            raise FileNotFoundError(f"Database directory is missing: {database}")

    logs = work / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    run_env = os.environ.copy()
    for key in (
        "http_proxy",
        "https_proxy",
        "HTTP_PROXY",
        "HTTPS_PROXY",
        "ALL_PROXY",
        "all_proxy",
    ):
        run_env.pop(key, None)
    home = work / "virsorter2-home"
    home.mkdir(parents=True, exist_ok=True)
    run_env["HOME"] = str(home)
    threads = str(args.threads)
    command_rows = []

    command_rows.append(
        timed_run(
            label="genomad",
            command=[
                "conda",
                "run",
                "-n",
                "metagenome-virus-discovery-2026.07",
                "genomad",
                "end-to-end",
                "--cleanup",
                "--restart",
                "--threads",
                threads,
                "--splits",
                "8",
                str(input_fna),
                str(results / "genomad"),
                str(args.genomad_db.resolve()),
            ],
            cwd=root,
            logs=logs,
            env=run_env,
        )
    )
    command_rows.append(
        timed_run(
            label="checkv",
            command=[
                "conda",
                "run",
                "-n",
                "metagenome-virus-discovery-2026.07",
                "checkv",
                "end_to_end",
                str(input_fna),
                str(results / "checkv"),
                "-d",
                str(args.checkv_db.resolve()),
                "-t",
                threads,
                "--restart",
            ],
            cwd=root,
            logs=logs,
            env=run_env,
        )
    )

    config_command = [
        "conda",
        "run",
        "-n",
        "metagenome-virsorter2-2026.07",
        "virsorter",
        "config",
        "--init-source",
        "--db-dir",
        str(args.virsorter2_db.resolve()),
    ]
    subprocess.run(config_command, cwd=root, env=run_env, check=True)
    command_rows.append(
        timed_run(
            label="virsorter2",
            command=[
                "conda",
                "run",
                "-n",
                "metagenome-virsorter2-2026.07",
                "virsorter",
                "run",
                "-w",
                str(results / "virsorter2"),
                "-d",
                str(args.virsorter2_db.resolve()),
                "-i",
                str(input_fna),
                "--min-length",
                "1500",
                "--keep-original-seq",
                "--include-groups",
                "dsDNAphage,ssDNA",
                "--min-score",
                "0.5",
                "-j",
                threads,
                "all",
            ],
            cwd=root,
            logs=logs,
            env=run_env,
        )
    )

    votu = results / "votu"
    votu.mkdir(parents=True, exist_ok=True)
    blast_db = votu / "checkv-test"
    commands = [
        [
            "makeblastdb",
            "-in",
            str(input_fna),
            "-dbtype",
            "nucl",
            "-out",
            str(blast_db),
        ],
        [
            "blastn",
            "-query",
            str(input_fna),
            "-db",
            str(blast_db),
            "-out",
            str(votu / "all-vs-all-blast.tsv"),
            "-outfmt",
            "6 std qlen slen",
            "-max_target_seqs",
            "10000",
            "-num_threads",
            threads,
        ],
        [
            "anicalc",
            "-i",
            str(votu / "all-vs-all-blast.tsv"),
            "-o",
            str(votu / "ani.tsv"),
        ],
        [
            "aniclust",
            "--fna",
            str(input_fna),
            "--ani",
            str(votu / "ani.tsv"),
            "--out",
            str(votu / "votu-clusters.tsv"),
            "--min_ani",
            "95",
            "--min_tcov",
            "85",
        ],
    ]
    for index, command in enumerate(commands, start=1):
        wrapped = [
            "conda",
            "run",
            "-n",
            "metagenome-virus-discovery-2026.07",
            *command,
        ]
        completed = subprocess.run(wrapped, cwd=root, env=run_env, check=False)
        command_rows.append(
            {
                "Label": f"votu-{index}",
                "Command": shlex.join(wrapped),
                "ExitStatus": completed.returncode,
                "Stdout": "terminal",
                "Stderr": "terminal",
                "Timing": "not separately measured",
            }
        )
        if completed.returncode != 0:
            raise RuntimeError(f"vOTU command failed: {shlex.join(wrapped)}")

    write_tsv(work / "command-log.tsv", command_rows)
    print(f"Article 54 run complete: {results}")


if __name__ == "__main__":
    main()
