#!/usr/bin/env python3
"""Run geNomad taxonomy and 95/85 vOTU clustering for Article 55."""

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
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    results = args.results_dir.resolve()
    database = args.genomad_db.resolve()
    input_fna = work / "input/cook15-phage-reference.fna"
    if not input_fna.is_file():
        raise FileNotFoundError(
            "Run prepare_article55_virus_abundance.py before this command"
        )
    if not database.is_dir():
        raise FileNotFoundError(f"geNomad database directory is missing: {database}")

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
    run_env["HOME"] = str(work / "runtime-home")
    Path(run_env["HOME"]).mkdir(parents=True, exist_ok=True)
    env_prefix = ["conda", "run", "-n", "metagenome-virus-discovery-2026.07"]
    threads = str(args.threads)
    command_rows = []

    command_rows.append(
        timed_run(
            label="genomad",
            command=[
                *env_prefix,
                "genomad",
                "end-to-end",
                "--cleanup",
                "--restart",
                "--threads",
                threads,
                "--splits",
                "16",
                str(input_fna),
                str(results / "genomad"),
                str(database),
            ],
            cwd=root,
            logs=logs,
            env=run_env,
        )
    )

    votu = results / "votu"
    votu.mkdir(parents=True, exist_ok=True)
    blast_db = votu / "cook15"
    commands = (
        (
            "makeblastdb",
            [
                *env_prefix,
                "makeblastdb",
                "-in",
                str(input_fna),
                "-dbtype",
                "nucl",
                "-out",
                str(blast_db),
            ],
        ),
        (
            "blastn",
            [
                *env_prefix,
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
        ),
        (
            "anicalc",
            [
                *env_prefix,
                "anicalc",
                "-i",
                str(votu / "all-vs-all-blast.tsv"),
                "-o",
                str(votu / "ani.tsv"),
            ],
        ),
        (
            "aniclust",
            [
                *env_prefix,
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
        ),
    )
    for label, command in commands:
        command_rows.append(
            timed_run(
                label=label,
                command=command,
                cwd=root,
                logs=logs,
                env=run_env,
            )
        )

    write_tsv(work / "command-log.tsv", command_rows)
    print(f"Article 55 taxonomy and vOTU run complete: {results}")


if __name__ == "__main__":
    main()
