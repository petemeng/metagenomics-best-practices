#!/usr/bin/env python3
"""Run antiSMASH 8.0.4 and GECCO 0.10.3 on Article 40 real inputs."""

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
    parser.add_argument("--antismash-database", type=Path, required=True)
    parser.add_argument("--antismash-env", type=Path, required=True)
    parser.add_argument("--gecco-env", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=32)
    return parser.parse_args()


def run(label: str, command: list[str], work: Path, runtime: dict[str, str]) -> dict[str, object]:
    stdout = work / "logs" / f"{label}.stdout.log"
    stderr = work / "logs" / f"{label}.stderr.log"
    resource = work / "logs" / f"{label}.time.txt"
    started = datetime.now(timezone.utc)
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        result = subprocess.run(["/usr/bin/time", "-v", "-o", str(resource), *command], env=runtime, stdout=out, stderr=err, text=True)
    ended = datetime.now(timezone.utc)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed ({result.returncode})\n{stderr.read_text(encoding='utf-8', errors='replace')[-5000:]}")
    return {
        "Label": label, "Command": " ".join(command), "ReturnCode": result.returncode,
        "StartedUTC": started.isoformat(), "EndedUTC": ended.isoformat(),
        "Stdout": str(stdout), "Stderr": str(stderr), "ResourceLog": str(resource),
    }


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    database = args.antismash_database.resolve()
    antismash_env = args.antismash_env.resolve()
    gecco_env = args.gecco_env.resolve()
    if not (work / ".article40-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article40_bgc.py first")
    for tool in ("antismash", "gecco"):
        target = work / tool
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True)
    antismash_runtime = os.environ.copy()
    antismash_runtime.update({
        "PATH": f"{antismash_env / 'bin'}:/usr/bin:/bin", "MPLCONFIGDIR": str(work / "tmp/antismash-mpl"),
        "PYTHONHASHSEED": "0", "PYTHONNOUSERSITE": "1", "PYTHONPATH": "", "OMP_NUM_THREADS": str(args.threads),
        "TMPDIR": str(work / "tmp"),
    })
    gecco_runtime = os.environ.copy()
    gecco_runtime.update({
        "PATH": f"{gecco_env / 'bin'}:/usr/bin:/bin", "PYTHONHASHSEED": "0", "PYTHONNOUSERSITE": "1",
        "PYTHONPATH": "", "OMP_NUM_THREADS": str(args.threads), "TMPDIR": str(work / "tmp"),
    })
    antismash = str(antismash_env / "bin/antismash")
    gecco = str(gecco_env / "bin/gecco")
    specifications = (
        ("salinispora-full", "salinispora-full.fna", "prodigal"),
        ("salinispora-fragmented", "salinispora-fragmented.fna", "prodigal-m"),
        ("nostoc", "nostoc.fna", "prodigal"),
        ("coassembly-ge20kb", "coassembly-ge20kb.fna", "prodigal-m"),
    )
    rows: list[dict[str, object]] = []
    for label, input_name, gene_finder in specifications:
        anti_out = work / "antismash" / label
        rows.append(run(f"antismash-{label}", [
            antismash, "--taxon", "bacteria", "--cpus", str(args.threads), "--databases", str(database),
            "--output-dir", str(anti_out), "--output-basename", label,
            "--genefinding-tool", gene_finder, "--minimal", "--cc-mibig", "--cb-knownclusters",
            "--no-enable-html", "--no-zip-output", "--no-region-gbks", str(work / "inputs" / input_name),
        ], work, antismash_runtime))
        gecco_out = work / "gecco" / label
        rows.append(run(f"gecco-{label}", [
            gecco, "run", "--genome", str(work / "inputs" / input_name), "--jobs", str(args.threads),
            "--threshold", "0.8", "--cds", "3", "--output-dir", str(gecco_out),
            "--force-tsv", "--merge-gbk",
        ], work, gecco_runtime))
    with (work / "command-log.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    payload = {"article": 40, "commands": len(rows), "threads": args.threads, "completed_utc": datetime.now(timezone.utc).isoformat()}
    (work / ".article40-run-complete").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
