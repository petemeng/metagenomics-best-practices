#!/usr/bin/env python3
"""Summarize long-read recruitment and split-alignment alarms from PAF."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import statistics
from pathlib import Path
from typing import TextIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--branch", required=True)
    parser.add_argument("--assembler", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--expected-reads", type=int, required=True)
    parser.add_argument("--expected-bases", type=int, required=True)
    parser.add_argument("--reference-threshold-bp", type=int, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def union_length(intervals: list[tuple[int, int]]) -> tuple[int, int, int, int]:
    if not intervals:
        return 0, 0, 0, 0
    ordered = sorted(intervals)
    merged: list[list[int]] = []
    for start, end in ordered:
        if end < start:
            raise ValueError(f"Invalid query interval: {start}-{end}")
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    covered = sum(end - start for start, end in merged)
    maximum_gap = max(
        (merged[index + 1][0] - merged[index][1] for index in range(len(merged) - 1)),
        default=0,
    )
    return covered, merged[0][0], merged[-1][1], maximum_gap


def main() -> int:
    args = parse_args()
    reads: dict[str, dict[str, object]] = {}
    paf_records = 0
    base_level_records = 0
    secondary_records = 0
    with open_text(args.paf) as handle:
        for number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {number}: {args.paf}")
            tags = {
                tag.split(":", 2)[0]: tag.split(":", 2)[2]
                for tag in fields[12:]
                if tag.count(":") >= 2
            }
            if tags.get("tp", "P") != "P":
                secondary_records += 1
                continue
            if "cg" not in tags or "NM" not in tags:
                raise ValueError(
                    f"PAF line {number} lacks base-level cg/NM tags; rerun minimap2 with -c: {args.paf}"
                )
            qname = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            target = fields[5]
            matches = int(fields[9])
            block = int(fields[10])
            mapq = int(fields[11])
            if (
                qlen <= 0
                or block <= 0
                or not (0 <= matches <= block)
                or not (0 <= qstart <= qend <= qlen)
            ):
                raise ValueError(f"Invalid PAF geometry on line {number}: {args.paf}")
            row = reads.setdefault(
                qname,
                {
                    "qlen": qlen,
                    "intervals": [],
                    "targets": set(),
                    "segments": 0,
                    "matches": 0,
                    "blocks": 0,
                    "best_identity": 0.0,
                    "best_mapq": 0,
                },
            )
            if row["qlen"] != qlen:
                raise ValueError(f"Inconsistent query length for {qname}")
            row["intervals"].append((qstart, qend))
            row["targets"].add(target)
            row["segments"] += 1
            row["matches"] += matches
            row["blocks"] += block
            row["best_identity"] = max(float(row["best_identity"]), matches / block)
            row["best_mapq"] = max(int(row["best_mapq"]), mapq)
            paf_records += 1
            base_level_records += 1

    if len(reads) > args.expected_reads:
        raise ValueError("More mapped query names than expected FASTQ records")
    aligned_union = 0
    mapq20 = 0
    multi_segment = 0
    multi_contig = 0
    end_clipped = 0
    internal_gap = 0
    read_fractions: list[float] = []
    best_identities: list[float] = []
    total_matches = 0
    total_blocks = 0
    for row in reads.values():
        qlen = int(row["qlen"])
        covered, left, right, maximum_gap = union_length(row["intervals"])
        aligned_union += covered
        read_fractions.append(covered / qlen)
        best_identities.append(float(row["best_identity"]))
        mapq20 += int(row["best_mapq"] >= 20)
        multi_segment += int(row["segments"] >= 2)
        multi_contig += int(len(row["targets"]) >= 2)
        end_clipped += int(left >= 1_000 or qlen - right >= 1_000)
        internal_gap += int(maximum_gap >= 1_000)
        total_matches += int(row["matches"])
        total_blocks += int(row["blocks"])

    mapped = len(reads)
    result = {
        "Branch": args.branch,
        "Assembler": args.assembler,
        "Platform": args.platform,
        "ExpectedReads": args.expected_reads,
        "ExpectedBases": args.expected_bases,
        "ReferenceThresholdBp": args.reference_threshold_bp,
        "ReferenceBytes": args.reference_fasta.stat().st_size,
        "ReferenceSHA256": sha256(args.reference_fasta),
        "MappedReads": mapped,
        "MappedReadPct": 100 * mapped / args.expected_reads,
        "MapQ20Reads": mapq20,
        "MapQ20ReadPct": 100 * mapq20 / args.expected_reads,
        "AlignedQueryBasesUnion": aligned_union,
        "AlignedQueryBasePct": 100 * aligned_union / args.expected_bases,
        "IdentityComputation": "minimap2-c-paf-matches-over-block",
        "ExactMatchingBases": total_matches,
        "AlignmentBlockBases": total_blocks,
        "WeightedAlignmentIdentityPct": 100 * total_matches / total_blocks
        if total_blocks
        else 0.0,
        "MedianBestSegmentIdentityPct": 100 * statistics.median(best_identities)
        if best_identities
        else 0.0,
        "MedianReadAlignedFractionPct": 100 * statistics.median(read_fractions)
        if read_fractions
        else 0.0,
        "MultiSegmentReads": multi_segment,
        "MultiSegmentPctOfMapped": 100 * multi_segment / mapped if mapped else 0.0,
        "MultiContigReads": multi_contig,
        "MultiContigPctOfMapped": 100 * multi_contig / mapped if mapped else 0.0,
        "EndClippedGE1kbReads": end_clipped,
        "EndClippedGE1kbPctOfMapped": 100 * end_clipped / mapped if mapped else 0.0,
        "InternalQueryGapGE1kbReads": internal_gap,
        "InternalQueryGapGE1kbPctOfMapped": 100 * internal_gap / mapped if mapped else 0.0,
        "PrimaryPAFRecords": paf_records,
        "BaseLevelCigarPAFRecords": base_level_records,
        "DiscardedSecondaryPAFRecords": secondary_records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
