#!/usr/bin/env python3
"""Run and replay the deterministic Article 53 supplementary-data analysis."""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
from pathlib import Path

from article41_44_utils import parse_time, read_tsv, sha256, write_tsv


def timed(label: str, command: list[str], work: Path, env: dict[str, str]) -> dict[str, object]:
    stdout_path = work / "logs" / f"{label}.stdout.log"
    stderr_path = work / "logs" / f"{label}.stderr.log"
    time_path = work / "logs" / f"{label}.time.txt"
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        completed = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(time_path), *command],
            stdout=stdout,
            stderr=stderr,
            env=env,
            check=False,
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


def version(environment: str, expression: str) -> str:
    completed = subprocess.run(
        ["conda", "run", "-n", environment, "python", "-c", expression],
        text=True,
        capture_output=True,
        check=True,
    )
    return completed.stdout.strip()


def verify_assets(work: Path) -> None:
    rows = []
    for row in read_tsv(work / "asset-manifest.tsv"):
        path = work / row["SourceFile"]
        exists = path.is_file()
        observed_bytes = path.stat().st_size if exists else -1
        observed_hash = sha256(path) if exists else "MISSING"
        passed = (
            exists
            and observed_bytes == int(row["Bytes"])
            and observed_hash == row["SHA256"]
        )
        rows.append(
            {
                "Asset": row["Asset"],
                "SourceFile": row["SourceFile"],
                "ExpectedBytes": row["Bytes"],
                "ObservedBytes": observed_bytes,
                "ExpectedSHA256": row["SHA256"],
                "ObservedSHA256": observed_hash,
                "ChecksumPass": passed,
            }
        )
    write_tsv(work / "asset-check-audit.tsv", rows)
    if len(rows) != 5 or not all(row["ChecksumPass"] for row in rows):
        raise ValueError("Refusing to parse failed Article 53 assets")


def compare_directories(primary: Path, replay: Path, work: Path) -> None:
    primary_files = sorted(
        path.relative_to(primary) for path in primary.rglob("*") if path.is_file()
    )
    replay_files = sorted(
        path.relative_to(replay) for path in replay.rglob("*") if path.is_file()
    )
    if primary_files != replay_files:
        raise ValueError("Primary and replay output file sets differ")
    rows = []
    for relative in primary_files:
        left = primary / relative
        right = replay / relative
        same = left.stat().st_size == right.stat().st_size and sha256(left) == sha256(right)
        rows.append(
            {
                "File": relative.as_posix(),
                "PrimaryBytes": left.stat().st_size,
                "ReplayBytes": right.stat().st_size,
                "PrimarySHA256": sha256(left),
                "ReplaySHA256": sha256(right),
                "ByteIdentical": same,
                "Seed": 20260753,
                "RandomOutputRequested": False,
            }
        )
    write_tsv(work / "determinism-audit.tsv", rows)
    if not rows or not all(row["ByteIdentical"] for row in rows):
        raise ValueError("Article 53 deterministic replay mismatch")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument(
        "--environment", default="metagenome-strain-transmission-2026.07"
    )
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    if not (work / ".article53-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article53_transmission.py first")
    verify_assets(work)
    (work / "logs").mkdir(exist_ok=True)
    replay = work / "replay-summary"
    if replay.exists():
        shutil.rmtree(replay)

    env = os.environ.copy()
    env.update(
        {
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
    )
    script = root / "scripts/summarize_article53_transmission.py"
    commands = []
    for label, output in (
        ("transmission-primary", work / "summary"),
        ("transmission-replay", replay),
    ):
        command = [
            "conda",
            "run",
            "-n",
            args.environment,
            "python",
            str(script),
            "--project-root",
            str(root),
            "--work-dir",
            str(work),
            "--output-dir",
            str(output),
        ]
        commands.append(timed(label, command, work, env))
    compare_directories(work / "summary", replay, work)
    write_tsv(work / "command-log.tsv", commands)
    write_tsv(
        work / "resource-summary.tsv",
        [parse_time(work / "logs" / f"{row['Label']}.time.txt") for row in commands],
    )
    write_tsv(
        work / "tool-versions.tsv",
        [
            {
                "Software": "Python",
                "Version": version(
                    args.environment,
                    "import sys; print('.'.join(map(str, sys.version_info[:3])))",
                ),
            },
            {
                "Software": "openpyxl",
                "Version": version(
                    args.environment, "import openpyxl; print(openpyxl.__version__)"
                ),
            },
            {
                "Software": "lxml",
                "Version": version(args.environment, "import lxml; print(lxml.__version__)"),
            },
        ],
    )
    (work / ".article53-run-complete").write_text("complete\n", encoding="utf-8")
    print(f"Article 53 replay verified across {len(list((work / 'summary').iterdir()))} outputs")


if __name__ == "__main__":
    main()

