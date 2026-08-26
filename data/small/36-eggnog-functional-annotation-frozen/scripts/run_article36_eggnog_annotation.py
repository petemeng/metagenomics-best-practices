#!/usr/bin/env python3
"""Run the one-time full Article 36 eggNOG-mapper annotation and GO audit."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import platform
import shlex
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--database-dir", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=32)
    parser.add_argument("--skip-main", action="store_true")
    return parser.parse_args()


def run_timed(
    command: list[str],
    env: dict[str, str],
    stdout_path: Path,
    stderr_path: Path,
    time_path: Path,
) -> None:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(time_path), *command],
            env=env,
            stdout=stdout,
            stderr=stderr,
            text=True,
            check=False,
        )
    if completed.returncode != 0:
        tail = "\n".join(stderr_path.read_text(encoding="utf-8", errors="replace").splitlines()[-40:])
        raise RuntimeError(f"Command failed ({completed.returncode}): {shlex.join(command)}\n{tail}")


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    db_dir = args.database_dir.resolve()
    env_prefix = args.env_prefix.resolve()
    if not (work / ".article36-inputs-ready").is_file():
        raise FileNotFoundError("Run prepare_article36_eggnog_inputs.py first")
    if args.cpu < 1:
        raise ValueError("--cpu must be positive")

    emapper = env_prefix / "bin/emapper.py"
    query = work / "inputs/catalog.faa"
    main_dir = work / "annotation/main"
    go_all_dir = work / "annotation/go-all"
    main_dir.mkdir(parents=True, exist_ok=True)
    go_all_dir.mkdir(parents=True, exist_ok=True)
    runtime_env = os.environ.copy()
    runtime_env.update(
        {
            "PATH": f"{env_prefix / 'bin'}:/usr/bin:/bin",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "PYTHONPATH": "",
            "TMPDIR": str(work / "tmp"),
        }
    )

    common_annotation = [
        "--seed_ortholog_evalue", "0.001",
        "--tax_scope", "auto",
        "--tax_scope_mode", "inner_narrowest",
        "--target_orthologs", "all",
        "--pfam_realign", "none",
        "--data_dir", str(db_dir),
    ]
    main_command = [
        str(emapper),
        "-i", str(query),
        "--itype", "proteins",
        "-m", "diamond",
        "--sensmode", "sensitive",
        "--dmnd_iterate", "yes",
        "--evalue", "0.001",
        "--outfmt_short",
        *common_annotation,
        "--go_evidence", "non-electronic",
        "--md5",
        "--cpu", str(args.cpu),
        "-o", "catalog-main",
        "--output_dir", str(main_dir),
        "--override",
    ]
    sensitivity_command = [
        str(emapper),
        "-m", "no_search",
        "--annotate_hits_table", str(main_dir / "catalog-main.emapper.seed_orthologs"),
        *common_annotation,
        "--go_evidence", "all",
        "--cpu", str(min(args.cpu, 8)),
        "-o", "catalog-go-all",
        "--output_dir", str(go_all_dir),
        "--override",
    ]

    started = dt.datetime.now(dt.timezone.utc)
    commands: list[dict[str, str]] = []
    main_output = main_dir / "catalog-main.emapper.annotations"
    if not args.skip_main:
        run_timed(
            main_command,
            runtime_env,
            work / "logs/emapper-main.stdout.log",
            work / "logs/emapper-main.stderr.log",
            work / "logs/emapper-main.time.txt",
        )
        commands.append({"Step": "main-non-electronic", "Command": shlex.join(main_command), "Status": "PASS"})
    elif not main_output.is_file():
        raise FileNotFoundError(f"--skip-main requested but output is missing: {main_output}")
    else:
        commands.append({"Step": "main-non-electronic", "Command": shlex.join(main_command), "Status": "REUSED"})

    for required in (
        main_output,
        main_dir / "catalog-main.emapper.hits",
        main_dir / "catalog-main.emapper.seed_orthologs",
    ):
        if not required.is_file() or required.stat().st_size == 0:
            raise FileNotFoundError(required)

    run_timed(
        sensitivity_command,
        runtime_env,
        work / "logs/emapper-go-all.stdout.log",
        work / "logs/emapper-go-all.stderr.log",
        work / "logs/emapper-go-all.time.txt",
    )
    commands.append({"Step": "go-all-sensitivity", "Command": shlex.join(sensitivity_command), "Status": "PASS"})
    sensitivity_output = go_all_dir / "catalog-go-all.emapper.annotations"
    if not sensitivity_output.is_file() or sensitivity_output.stat().st_size == 0:
        raise FileNotFoundError(sensitivity_output)

    with (work / "command-log.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Step", "Command", "Status"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(commands)
    ended = dt.datetime.now(dt.timezone.utc)
    metadata = {
        "article": 36,
        "started_utc": started.isoformat(),
        "ended_utc": ended.isoformat(),
        "elapsed_seconds": (ended - started).total_seconds(),
        "hostname": platform.node(),
        "platform": platform.platform(),
        "cpu_requested": args.cpu,
        "main_output_bytes": main_output.stat().st_size,
        "go_all_output_bytes": sensitivity_output.stat().st_size,
        "database_dir": str(db_dir),
        "environment_prefix": str(env_prefix),
    }
    (work / "run-metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (work / ".article36-run-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
