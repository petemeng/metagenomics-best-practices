#!/usr/bin/env python3
"""Summarize one real Article 08 read-prefix stratum deterministically."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--run-accession", required=True)
    parser.add_argument("--expected-rows", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def stable_float(value: float) -> float:
    return round(value, 6)


def main() -> int:
    args = parse_args()
    random.seed(args.seed)
    if args.seed != 20260719:
        raise SystemExit("Article 10 smoke seed must be 20260719")
    if not args.metrics.is_file():
        raise SystemExit(f"missing metrics table: {args.metrics}")

    row_count = 0
    total_bases = 0
    total_expected_errors = 0.0
    read_id_digest = hashlib.sha256()
    min_length: int | None = None
    max_length: int | None = None

    with args.metrics.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "PlatformKey",
            "RunAccession",
            "ReadID",
            "ReadLength",
            "ExpectedErrors",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"metrics table missing columns: {sorted(missing)}")
        for row in reader:
            if row["PlatformKey"] != args.platform:
                continue
            if row["RunAccession"] != args.run_accession:
                raise SystemExit(
                    f"unexpected run for {args.platform}: {row['RunAccession']}"
                )
            length = int(row["ReadLength"])
            row_count += 1
            total_bases += length
            total_expected_errors += float(row["ExpectedErrors"])
            read_id_digest.update(row["ReadID"].encode("utf-8"))
            read_id_digest.update(b"\n")
            min_length = length if min_length is None else min(min_length, length)
            max_length = length if max_length is None else max(max_length, length)

    if row_count != args.expected_rows:
        raise SystemExit(
            f"{args.platform}: expected {args.expected_rows} rows, observed {row_count}"
        )

    result = {
        "status": "passed",
        "scope": "workflow-plumbing-smoke-not-assembler-benchmark",
        "platform": args.platform,
        "run_accession": args.run_accession,
        "seed": args.seed,
        "rows": row_count,
        "total_bases": total_bases,
        "mean_read_length": stable_float(total_bases / row_count),
        "minimum_read_length": min_length,
        "maximum_read_length": max_length,
        "mean_expected_errors": stable_float(total_expected_errors / row_count),
        "read_id_sha256": read_id_digest.hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
