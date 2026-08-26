#!/usr/bin/env python3
"""Stream Bowtie2 SAM and emit compact, denominator-aware Article 35 counts."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
POLICIES = ("AllPrimary", "IdentityQcov", "Main", "Strict")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--total-reads", type=int, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--quality-histogram", type=Path, required=True)
    return parser.parse_args()


def alignment_metrics(cigar: str, sequence: str, fields: list[str]) -> tuple[int, int, float, float]:
    operations = CIGAR.findall(cigar)
    if not operations or "".join(f"{n}{op}" for n, op in operations) != cigar:
        raise ValueError(f"Unsupported CIGAR: {cigar}")
    query_aligned = sum(int(n) for n, op in operations if op in {"M", "I", "=", "X"})
    alignment_columns = sum(int(n) for n, op in operations if op in {"M", "I", "D", "=", "X"})
    query_length = len(sequence) if sequence != "*" else sum(int(n) for n, op in operations if op in {"M", "I", "S", "=", "X"})
    tags = {item.split(":", 2)[0]: item.split(":", 2)[2] for item in fields if item.count(":") >= 2}
    if "NM" not in tags:
        raise ValueError("Bowtie2 alignment lacks NM tag")
    edit_distance = int(tags["NM"])
    if query_aligned <= 0 or alignment_columns <= 0 or query_length <= 0:
        raise ValueError(f"Non-positive aligned/read length: {cigar}")
    identity = max(0.0, (alignment_columns - edit_distance) / alignment_columns)
    query_coverage = query_aligned / query_length
    return query_aligned, query_length, identity, query_coverage


def mapq_bin(value: int) -> str:
    if value == 0:
        return "0"
    if value < 10:
        return "1-9"
    if value < 20:
        return "10-19"
    if value < 40:
        return "20-39"
    return "40+"


def fraction_bin(value: float) -> str:
    if value < 0.80:
        return "<0.80"
    if value < 0.90:
        return "0.80-0.89"
    if value < 0.95:
        return "0.90-0.94"
    if value < 0.97:
        return "0.95-0.96"
    if value < 0.99:
        return "0.97-0.98"
    return "0.99-1.00"


def main() -> int:
    args = parse_args()
    if args.total_reads <= 0:
        raise SystemExit("--total-reads must be positive")
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    histograms: dict[str, Counter[str]] = {
        "MAPQ": Counter(),
        "Identity": Counter(),
        "QueryCoverage": Counter(),
    }
    primary = 0
    secondary = 0
    supplementary = 0
    malformed = 0
    policy_totals = Counter()

    for line_number, raw in enumerate(sys.stdin, start=1):
        if not raw or raw.startswith("@"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 11:
            malformed += 1
            continue
        flag = int(fields[1])
        if flag & 0x4:
            continue
        if flag & 0x800:
            supplementary += 1
            continue
        if flag & 0x100:
            secondary += 1
            continue
        gene = fields[2]
        mapq = int(fields[4])
        _, _, identity, query_coverage = alignment_metrics(fields[5], fields[9], fields[11:])
        primary += 1
        histograms["MAPQ"][mapq_bin(mapq)] += 1
        histograms["Identity"][fraction_bin(identity)] += 1
        histograms["QueryCoverage"][fraction_bin(query_coverage)] += 1

        passed = {
            "AllPrimary": True,
            "IdentityQcov": identity >= 0.95 and query_coverage >= 0.80,
            "Main": mapq >= 10 and identity >= 0.95 and query_coverage >= 0.80,
            "Strict": mapq >= 20 and identity >= 0.97 and query_coverage >= 0.90,
        }
        for policy, keep in passed.items():
            if keep:
                counts[gene][policy] += 1
                policy_totals[policy] += 1

    if malformed:
        raise ValueError(f"Malformed SAM records: {malformed}")
    if primary > args.total_reads:
        raise ValueError(f"Primary mapped records exceed input reads: {primary} > {args.total_reads}")
    if not (policy_totals["Strict"] <= policy_totals["Main"] <= policy_totals["IdentityQcov"] <= primary):
        raise ValueError("Mapping-policy nesting failed")

    args.counts.parent.mkdir(parents=True, exist_ok=True)
    with args.counts.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["GeneID", *POLICIES], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in sorted(counts):
            writer.writerow({"GeneID": gene, **{policy: counts[gene][policy] for policy in POLICIES}})

    metric_order = {
        "MAPQ": ["0", "1-9", "10-19", "20-39", "40+"],
        "Identity": ["<0.80", "0.80-0.89", "0.90-0.94", "0.95-0.96", "0.97-0.98", "0.99-1.00"],
        "QueryCoverage": ["<0.80", "0.80-0.89", "0.90-0.94", "0.95-0.96", "0.97-0.98", "0.99-1.00"],
    }
    with args.quality_histogram.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Sample", "Metric", "Bin", "Count"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for metric in ("MAPQ", "Identity", "QueryCoverage"):
            for label in metric_order[metric]:
                writer.writerow({"Sample": args.sample, "Metric": metric, "Bin": label, "Count": histograms[metric][label]})

    summary = {
        "sample": args.sample,
        "total_input_reads": args.total_reads,
        "primary_mapped_reads": primary,
        "unmapped_reads": args.total_reads - primary,
        "secondary_alignments": secondary,
        "supplementary_alignments": supplementary,
        "policy_assigned_reads": {policy: policy_totals[policy] for policy in POLICIES},
        "policy_detected_genes": {policy: sum(1 for values in counts.values() if values[policy] > 0) for policy in POLICIES},
        "main_filter": {"minimum_mapq": 10, "minimum_identity": 0.95, "minimum_query_coverage": 0.80},
        "strict_filter": {"minimum_mapq": 20, "minimum_identity": 0.97, "minimum_query_coverage": 0.90},
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
