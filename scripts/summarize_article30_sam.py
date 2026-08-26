#!/usr/bin/env python3
"""Summarize a Bowtie2 SAM stream without retaining BAM/SAM files."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assembly-branch", required=True)
    parser.add_argument("--expected-pairs", type=int, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    total_reads = 0
    mapped_reads = 0
    total_pairs = 0
    both_mapped_pairs = 0
    proper_pairs = 0
    singleton_pairs = 0
    discordant_pairs = 0

    for raw in sys.stdin:
        if not raw or raw.startswith("@"):
            continue
        fields = raw.split("\t", 3)
        if len(fields) < 2:
            raise SystemExit("Malformed SAM alignment line")
        flag = int(fields[1])
        if flag & (0x100 | 0x800):
            continue
        total_reads += 1
        if not flag & 0x4:
            mapped_reads += 1
        if flag & 0x40:
            total_pairs += 1
            this_mapped = not flag & 0x4
            mate_mapped = not flag & 0x8
            if this_mapped and mate_mapped:
                both_mapped_pairs += 1
                if flag & 0x2:
                    proper_pairs += 1
                else:
                    discordant_pairs += 1
            elif this_mapped != mate_mapped:
                singleton_pairs += 1

    if total_pairs != args.expected_pairs or total_reads != 2 * args.expected_pairs:
        raise SystemExit(
            f"SAM primary-record invariant failed: reads={total_reads}, "
            f"pairs={total_pairs}, expected_pairs={args.expected_pairs}"
        )
    payload = {
        "sample": args.sample,
        "assembly_branch": args.assembly_branch,
        "total_primary_reads": total_reads,
        "mapped_primary_reads": mapped_reads,
        "mapped_read_fraction": mapped_reads / total_reads,
        "total_pairs": total_pairs,
        "both_mapped_pairs": both_mapped_pairs,
        "both_mapped_pair_fraction": both_mapped_pairs / total_pairs,
        "proper_pairs": proper_pairs,
        "proper_pair_fraction": proper_pairs / total_pairs,
        "discordant_pairs": discordant_pairs,
        "singleton_pairs": singleton_pairs,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
