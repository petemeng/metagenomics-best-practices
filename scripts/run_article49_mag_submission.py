#!/usr/bin/env python3
"""Time the deterministic Article 49 review-package build and retain its logs."""

from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from article41_44_utils import parse_time, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    command = [
        sys.executable,
        str(root / "scripts/prepare_article49_mag_submission.py"),
        "--project-root", str(root),
        "--work-dir", str(work),
    ]
    with tempfile.TemporaryDirectory(prefix="article49-") as temporary:
        temporary_path = Path(temporary)
        stdout_path = temporary_path / "article49-prep.stdout.log"
        stderr_path = temporary_path / "article49-prep.stderr.log"
        time_path = temporary_path / "article49-prep.time.txt"
        timed = ["/usr/bin/time", "-v", "-o", str(time_path), *command]
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
            completed = subprocess.run(timed, stdout=stdout, stderr=stderr, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"Article 49 preparation failed; temporary logs: {temporary}")
        logs = work / "logs"
        logs.mkdir(exist_ok=True)
        for path in (stdout_path, stderr_path, time_path):
            shutil.copy2(path, logs / path.name)
    write_tsv(work / "command-log.tsv", [{
        "Label": "article49-prep",
        "ExitStatus": completed.returncode,
        "Command": shlex.join(command),
        "Stdout": str(logs / "article49-prep.stdout.log"),
        "Stderr": str(logs / "article49-prep.stderr.log"),
        "TimeLog": str(logs / "article49-prep.time.txt"),
    }])
    write_tsv(work / "tool-versions.tsv", [{
        "Tool": "Python",
        "Version": sys.version.split()[0],
        "Role": "deterministic evidence assembly and robust anomaly screening",
    }])
    write_tsv(work / "summary/resource-summary.tsv", [parse_time(logs / "article49-prep.time.txt")])
    (work / ".article49-run-complete").write_text("complete\n", encoding="utf-8")
    print("Article 49 review package build completed")


if __name__ == "__main__":
    main()
