#!/usr/bin/env python3
"""Run RGI 6.0.8 on the real protein catalog, co-assembly, and two controls."""

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
    parser.add_argument("--database-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=32)
    return parser.parse_args()


def run(label: str, command: list[str], work: Path, env: dict[str, str]) -> dict[str, object]:
    stdout = work / "logs" / f"{label}.stdout.log"
    stderr = work / "logs" / f"{label}.stderr.log"
    resource = work / "logs" / f"{label}.time.txt"
    started = datetime.now(timezone.utc)
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        result = subprocess.run(["/usr/bin/time", "-v", "-o", str(resource), *command], env=env, stdout=out, stderr=err, text=True)
    ended = datetime.now(timezone.utc)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed ({result.returncode})\n{stderr.read_text(encoding='utf-8', errors='replace')[-4000:]}")
    return {
        "Label": label, "Command": " ".join(command), "ReturnCode": result.returncode,
        "StartedUTC": started.isoformat(), "EndedUTC": ended.isoformat(),
        "Stdout": str(stdout), "Stderr": str(stderr), "ResourceLog": str(resource),
    }


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    database = args.database_root.resolve()
    env_prefix = args.env_prefix.resolve()
    if not (work / ".article38-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article38_resistome.py first")
    for name in ("catalog", "coassembly", "pseudomonas", "staphylococcus"):
        target = work / name
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True)
    runtime = os.environ.copy()
    runtime.update({
        "PATH": f"{env_prefix / 'bin'}:/usr/bin:/bin", "DATA_PATH": str(database / "rgi-db"),
        "MPLCONFIGDIR": str(work / "tmp/matplotlib"), "PYTHONHASHSEED": "0",
        "PYTHONNOUSERSITE": "1", "PYTHONPATH": "", "OMP_NUM_THREADS": str(args.threads),
        "TMPDIR": str(work / "tmp"),
    })
    rgi = str(env_prefix / "bin/rgi")
    shared = ["-a", "DIAMOND", "-n", str(args.threads), "--include_loose", "--clean"]
    rows = []
    rows.append(run("rgi-catalog-protein", [rgi, "main", "-i", str(work / "inputs/catalog.faa"), "-o", str(work / "catalog/catalog"), "-t", "protein", *shared], work, runtime))
    for label, input_name, output_dir, dtype in (
        ("rgi-coassembly-contig", "coassembly.fna", "coassembly", "wgs"),
        ("rgi-pseudomonas-control", "pseudomonas.fna", "pseudomonas", "chromosome"),
        ("rgi-staphylococcus-control", "staphylococcus.fna", "staphylococcus", "chromosome"),
    ):
        rows.append(run(label, [
            rgi, "main", "-i", str(work / "inputs" / input_name),
            "-o", str(work / output_dir / output_dir), "-t", "contig",
            *shared, "--low_quality", "-g", "PYRODIGAL", "-d", dtype,
        ], work, runtime))
    with (work / "command-log.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    payload = {"article": 38, "commands": len(rows), "threads": args.threads, "completed_utc": datetime.now(timezone.utc).isoformat()}
    (work / ".article38-run-complete").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
