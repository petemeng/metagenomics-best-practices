#!/usr/bin/env python3
"""Count conservative read support across software-declared circular junctions."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    normalized = work / "normalized"
    candidates = read_tsv(normalized / "circular-candidates.tsv")
    by_junction = {row["JunctionID"]: row for row in candidates}
    support: dict[str, set[str]] = {key: set() for key in by_junction}
    audited_records: dict[str, int] = {key: 0 for key in by_junction}
    base_level_records: dict[str, int] = {key: 0 for key in by_junction}

    branches = sorted({row["Branch"] for row in candidates})
    for branch in branches:
        paf = work / "junction-mapping" / f"{branch}.paf"
        eligible = {
            row["JunctionID"]
            for row in candidates
            if row["Branch"] == branch and row["JunctionEligible"] == "TRUE"
        }
        if not eligible:
            if paf.exists() and paf.stat().st_size:
                raise ValueError(f"Unexpected junction alignments for {branch}")
            continue
        if not paf.is_file():
            raise FileNotFoundError(paf)
        identity_floor = 0.75 if by_junction[next(iter(eligible))]["Platform"] == "ONT R9" else 0.95
        with paf.open(encoding="utf-8") as handle:
            for number, raw in enumerate(handle, start=1):
                if not raw.strip():
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 12:
                    raise ValueError(f"Malformed junction PAF line {number}: {paf}")
                tags = {
                    tag.split(":", 2)[0]: tag.split(":", 2)[2]
                    for tag in fields[12:]
                    if tag.count(":") >= 2
                }
                if tags.get("tp", "P") != "P":
                    continue
                if "cg" not in tags or "NM" not in tags:
                    raise ValueError(
                        f"Junction PAF line {number} lacks base-level cg/NM tags; rerun minimap2 with -c: {paf}"
                    )
                query = fields[0]
                target = fields[5]
                target_length = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                matches = int(fields[9])
                block = int(fields[10])
                mapq = int(fields[11])
                if target not in eligible:
                    raise ValueError(f"Unknown junction target in {paf}: {target}")
                if block <= 0 or not (0 <= matches <= block):
                    raise ValueError(f"Invalid base-level identity fields on line {number}: {paf}")
                audited_records[target] += 1
                base_level_records[target] += 1
                center = target_length // 2
                crosses = target_start <= center - 1_000 and target_end >= center + 1_000
                identity = matches / block if block else 0.0
                if crosses and mapq >= 20 and identity >= identity_floor:
                    support[target].add(query)

    rows: list[dict[str, object]] = []
    for row in candidates:
        target = row["JunctionID"]
        count = len(support[target])
        rows.append(
            {
                **row,
                "IdentityFloor": "0.75" if row["Platform"] == "ONT R9" else "0.95",
                "IdentityComputation": "minimap2-c-paf-matches-over-block",
                "AuditedPrimaryPAFRecords": audited_records[target],
                "BaseLevelCigarPAFRecords": base_level_records[target],
                "JunctionSpanningReads": count,
                "SupportedByGE3Reads": "TRUE" if count >= 3 else "FALSE",
            }
        )

    fields = list(candidates[0].keys()) if candidates else [
        "Branch",
        "Assembler",
        "Platform",
        "JunctionID",
        "ContigID",
        "LengthBp",
        "JunctionEligible",
        "JunctionFlankBp",
        "SequenceSHA256",
        "OriginalHeader",
    ]
    write_tsv(
        normalized / "junction-support.tsv",
        rows,
        fields
        + [
            "IdentityFloor",
            "IdentityComputation",
            "AuditedPrimaryPAFRecords",
            "BaseLevelCigarPAFRecords",
            "JunctionSpanningReads",
            "SupportedByGE3Reads",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
