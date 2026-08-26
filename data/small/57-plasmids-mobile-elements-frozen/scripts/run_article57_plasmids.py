#!/usr/bin/env python3
"""Run geNomad and truth alignment for Article 57."""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--genomad-db", type=Path, required=True)
    parser.add_argument("--genomad-env", type=Path, required=True)
    parser.add_argument("--minimap2-exe", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--splits", type=int, default=16)
    return parser.parse_args()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Empty command ledger")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


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
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    database = args.genomad_db.resolve()
    env_prefix = args.genomad_env.resolve()
    if not (work / ".article57-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article57_plasmids.py first")
    if not (database / "version.txt").is_file():
        raise FileNotFoundError(database / "version.txt")
    if (database / "version.txt").read_text(encoding="utf-8").strip() != "1.9":
        raise RuntimeError("Article 57 requires geNomad database v1.9")
    genomad = env_prefix / "bin/genomad"
    minimap2 = args.minimap2_exe.resolve()
    if not genomad.is_file():
        raise FileNotFoundError("geNomad executable missing from locked environment")
    if not minimap2.is_file():
        raise FileNotFoundError("minimap2 executable missing from locked assembly environment")

    logs = work / "logs"
    results = work / "results"
    logs.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    run_env = os.environ.copy()
    run_env["PATH"] = str(env_prefix / "bin") + os.pathsep + run_env.get("PATH", "")
    run_env["TMPDIR"] = str((work / "tmp").resolve())
    Path(run_env["TMPDIR"]).mkdir(parents=True, exist_ok=True)
    threads = str(args.threads)
    splits = str(args.splits)
    commands: list[dict[str, object]] = []

    branches = (
        ("reference-benchmark", work / "inputs/matched-reference-replicons.fna.gz"),
        ("coassembly", work / "inputs/coassembly.fna.gz"),
        ("staphylococcus", work / "inputs/staphylococcus-USA300.fna"),
    )
    for label, input_fasta in branches:
        completed_summaries = sorted((results / label).glob("*_summary/*_summary.json"))
        completed_plasmids = sorted(
            (results / label).glob("*_summary/*_plasmid_summary.tsv")
        )
        if len(completed_summaries) == 1 and len(completed_plasmids) == 1:
            commands.append(
                {
                    "Label": f"genomad-{label}",
                    "Command": "reused checksum-stable completed output",
                    "ExitStatus": 0,
                    "Stdout": str(logs / f"genomad-{label}.stdout.log"),
                    "Stderr": str(logs / f"genomad-{label}.stderr.log"),
                    "Timing": str(logs / f"genomad-{label}.time.txt"),
                }
            )
            continue
        commands.append(
            timed_run(
                f"genomad-{label}",
                [
                    str(genomad),
                    "end-to-end",
                    "--cleanup",
                    "--restart",
                    "--disable-find-proviruses",
                    "--threads",
                    threads,
                    "--splits",
                    splits,
                    str(input_fasta.resolve()),
                    str((results / label).resolve()),
                    str(database),
                ],
                root,
                logs,
                run_env,
            )
        )

    plasmid_fastas = sorted(
        (results / "coassembly").glob("*_summary/*_plasmid.fna")
    )
    if len(plasmid_fastas) != 1:
        raise RuntimeError(
            f"Expected one co-assembly plasmid FASTA, observed {plasmid_fastas}"
        )
    plasmid_fasta = plasmid_fastas[0]
    paf = results / "coassembly-plasmid-to-reference.paf"
    commands.append(
        timed_run(
            "minimap2-plasmid-truth",
            [
                str(minimap2),
                "-x",
                "asm5",
                "--secondary=no",
                "-t",
                threads,
                "-o",
                str(paf.resolve()),
                str((work / "inputs/mock2-all-references.fna.gz").resolve()),
                str(plasmid_fasta.resolve()),
            ],
            root,
            logs,
            run_env,
        )
    )
    write_tsv(work / "command-log.tsv", commands)
    (work / ".article57-run-complete").write_text("verified\n", encoding="utf-8")
    print(f"Article 57 geNomad run complete: {results}")


if __name__ == "__main__":
    main()
