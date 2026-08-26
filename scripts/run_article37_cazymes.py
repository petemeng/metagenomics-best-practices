#!/usr/bin/env python3
"""Run dbCAN 5.2.9 on the full real catalog and a real CGC substrate control."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--database-dir", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=32)
    return parser.parse_args()


def run(label: str, command: list[str], work: Path, env: dict[str, str]) -> dict[str, object]:
    stdout = work / "logs" / f"{label}.stdout.log"
    stderr = work / "logs" / f"{label}.stderr.log"
    resource = work / "logs" / f"{label}.time.txt"
    full = ["/usr/bin/time", "-v", "-o", str(resource), *command]
    started = datetime.now(timezone.utc)
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        result = subprocess.run(full, env=env, stdout=out, stderr=err, text=True, check=False)
    ended = datetime.now(timezone.utc)
    if result.returncode != 0:
        tail = stderr.read_text(encoding="utf-8", errors="replace")[-4000:]
        raise RuntimeError(f"{label} failed ({result.returncode})\n{tail}")
    return {
        "Label": label, "Command": " ".join(command), "ReturnCode": result.returncode,
        "StartedUTC": started.isoformat(), "EndedUTC": ended.isoformat(),
        "Stdout": str(stdout), "Stderr": str(stderr), "ResourceLog": str(resource),
    }


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    db = args.database_dir.resolve()
    env_prefix = args.env_prefix.resolve()
    if not (work / ".article37-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article37_cazymes.py first")
    for name in ("catalog", "btheta"):
        target = work / name
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True)
    runtime = os.environ.copy()
    runtime.update({
        "PATH": f"{env_prefix / 'bin'}:/usr/bin:/bin", "PYTHONHASHSEED": "0",
        "PYTHONNOUSERSITE": "1", "PYTHONPATH": "", "OMP_NUM_THREADS": str(args.threads),
        "TMPDIR": str(work / "tmp"),
    })
    binary = env_prefix / "bin/run_dbcan"
    rows = []
    rows.append(run("dbcan-catalog", [
        str(binary), "CAZyme_annotation", "--mode", "protein",
        "--input_raw_data", str(work / "inputs/catalog.faa"),
        "--output_dir", str(work / "catalog"), "--db_dir", str(db),
        "--methods", "diamond,hmm,dbCANsub", "--threads", str(args.threads),
        "--e_value_threshold", "1e-102", "--coverage_threshold_dbcan", "0.35",
        "--e_value_threshold_dbcan", "1e-15", "--coverage_threshold_dbsub", "0.35",
        "--e_value_threshold_dbsub", "1e-15", "--batch_size", "20000",
        "--batch_size_dbsub", "5000", "--log-level", "INFO",
        "--log-file", str(work / "logs/dbcan-catalog.internal.log"),
    ], work, runtime))
    rows.append(run("dbcan-btheta-substrate", [
        str(binary), "easy_substrate", "--mode", "prok",
        "--input_raw_data", str(work / "inputs/btheta.fna"),
        "--output_dir", str(work / "btheta"), "--db_dir", str(db),
        "--methods", "diamond,hmm,dbCANsub", "--threads", str(min(args.threads, 16)),
        "--e_value_threshold", "1e-102", "--coverage_threshold_dbcan", "0.35",
        "--e_value_threshold_dbcan", "1e-15", "--coverage_threshold_dbsub", "0.35",
        "--e_value_threshold_dbsub", "1e-15", "--num_null_gene", "2",
        "--additional_genes", "TC", "--log-level", "INFO",
        "--log-file", str(work / "logs/dbcan-btheta-substrate.internal.log"),
    ], work, runtime))
    with (work / "command-log.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    payload = {"article": 37, "commands": len(rows), "threads": args.threads, "completed_utc": datetime.now(timezone.utc).isoformat()}
    (work / ".article37-run-complete").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
