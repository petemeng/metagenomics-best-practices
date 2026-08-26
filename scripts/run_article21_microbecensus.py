#!/usr/bin/env python3
"""Run an unmodified MicrobeCensus model with a Python 3 preflight shim.

MicrobeCensus v1.1.1 supports Python 3 for its analysis code, but its RAPsearch2
version preflight compares subprocess bytes with text.  This wrapper replaces
only that preflight.  Marker data, alignment commands, thresholds, weights,
average-genome-size estimation, base counting, and output writing remain the
upstream implementation from the pinned source checkout.
"""

from __future__ import annotations

import argparse
import importlib
import json
import subprocess
import sys
from pathlib import Path


EXPECTED_COMMIT = "dfc42d356bfd7943633cde6c0fbfc0b116f29ae2"
EXPECTED_RAPSEARCH_BANNER = (
    "rapsearch v2.15: Fast protein similarity search tool for short reads"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--input-r1", type=Path, required=True)
    parser.add_argument("--input-r2", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    parser.add_argument("--read-length", type=int, choices=(100, 150), required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--nreads", type=int, default=100_000_000)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def git_value(root: Path, *args: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(root), *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def main() -> int:
    args = parse_args()
    source_root = args.source_root.resolve()
    commit = git_value(source_root, "rev-parse", "HEAD")
    if commit != EXPECTED_COMMIT:
        raise SystemExit(
            f"MicrobeCensus source commit mismatch: {commit} != {EXPECTED_COMMIT}"
        )

    sys.path.insert(0, str(source_root))
    module = importlib.import_module("microbe_census.microbe_census")
    upstream_check = module.check_rapsearch

    def check_rapsearch_python3(rapsearch: str) -> None:
        completed = subprocess.run(
            [rapsearch, "-h"],
            check=False,
            capture_output=True,
            text=True,
            errors="replace",
        )
        combined = "\n".join((completed.stdout, completed.stderr))
        if EXPECTED_RAPSEARCH_BANNER not in combined:
            raise RuntimeError(
                "MicrobeCensus requires the bundled RAPsearch2 v2.15 binary; "
                f"banner not found for {rapsearch}"
            )

    module.check_rapsearch = check_rapsearch_python3
    try:
        pipeline_args = {
            "seqfiles": [str(args.input_r1.resolve()), str(args.input_r2.resolve())],
            "outfile": str(args.output.resolve()),
            "nreads": args.nreads,
            "threads": args.threads,
            "read_length": args.read_length,
            "min_quality": -5,
            "mean_quality": -5,
            "filter_dups": False,
            "max_unknown": 100,
            "keep_tmp": False,
            "verbose": args.verbose,
            "no_equivs": False,
        }
        result = module.run_pipeline(pipeline_args)
        if result is None:
            raise RuntimeError("MicrobeCensus returned no result")
        average_genome_size, resolved_args = result
        total_bases = module.count_bases(resolved_args)
        module.report_results(resolved_args, average_genome_size, total_bases)
    finally:
        module.check_rapsearch = upstream_check

    metadata = {
        "status": "passed",
        "source_tag": "v1.1.1",
        "source_commit": commit,
        "internal_version": str(module.__version__),
        "compatibility_shim": "RAPsearch2 preflight bytes-to-text only",
        "read_length": int(args.read_length),
        "nreads_ceiling": int(args.nreads),
        "threads": int(args.threads),
        "reads_sampled": int(resolved_args["sampled_reads"]),
        "average_genome_size": float(average_genome_size),
        "total_bases": int(total_bases),
        "genome_equivalents": float(total_bases / average_genome_size),
        "quality_filters": False,
        "duplicate_filter": False,
    }
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_output.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metadata, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
