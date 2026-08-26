#!/usr/bin/env python3
"""Run ABRicate on the real gene catalog, co-assembly, and two controls."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--database-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=32)
    return parser.parse_args()


def run(label: str, command: list[str], output: Path, work: Path, env: dict[str, str]) -> dict[str, object]:
    stderr = work / "logs" / f"{label}.stderr.log"
    resource = work / "logs" / f"{label}.time.txt"
    started = datetime.now(timezone.utc)
    with output.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        result = subprocess.run(["/usr/bin/time", "-v", "-o", str(resource), *command], env=env, stdout=out, stderr=err, text=True)
    ended = datetime.now(timezone.utc)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed ({result.returncode})\n{stderr.read_text(encoding='utf-8', errors='replace')[-4000:]}")
    return {"Label": label, "Command": " ".join(command), "ReturnCode": result.returncode, "StartedUTC": started.isoformat(), "EndedUTC": ended.isoformat(), "Output": str(output), "Stderr": str(stderr), "ResourceLog": str(resource)}


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    database = args.database_root.resolve()
    env_prefix = args.env_prefix.resolve()
    if not (work / ".article39-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article39_virulome.py first")
    runtime = os.environ.copy()
    runtime.update({"PATH": f"{env_prefix / 'bin'}:/usr/bin:/bin", "LC_ALL": "C", "PYTHONHASHSEED": "0"})
    abricate = str(env_prefix / "bin/abricate")
    datadir = str(database / "abricate")
    common = [abricate, "--datadir", datadir, "--threads", str(args.threads), "--nopath"]
    specifications = (
        ("abricate-catalog-core-primary", "catalog.fna", "catalog-core-primary/hits.tsv", "vfdb-core", "90", "80"),
        ("abricate-catalog-core-sensitive", "catalog.fna", "catalog-core-sensitive/hits.tsv", "vfdb-core", "80", "80"),
        ("abricate-catalog-full-primary", "catalog.fna", "catalog-full-primary/hits.tsv", "vfdb-full", "90", "80"),
        ("abricate-coassembly-core-primary", "coassembly.fna", "coassembly/hits.tsv", "vfdb-core", "90", "80"),
        ("abricate-pseudomonas-core-primary", "pseudomonas.fna", "pseudomonas/hits.tsv", "vfdb-core", "90", "80"),
        ("abricate-staphylococcus-core-primary", "staphylococcus.fna", "staphylococcus/hits.tsv", "vfdb-core", "90", "80"),
    )
    rows = []
    for label, input_name, output_name, db, minid, mincov in specifications:
        output = work / output_name
        output.parent.mkdir(parents=True, exist_ok=True)
        rows.append(run(label, [*common, "--db", db, "--minid", minid, "--mincov", mincov, str(work / "inputs" / input_name)], output, work, runtime))
    with (work / "command-log.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    payload = {"article": 39, "commands": len(rows), "threads": args.threads, "completed_utc": datetime.now(timezone.utc).isoformat()}
    (work / ".article39-run-complete").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
