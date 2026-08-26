#!/usr/bin/env python3
"""Summarize CoverM breadth, depth, abundance, capture, and stringency evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path

from article41_44_utils import dump_json, parse_time, read_tsv, write_tsv


def number(value: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2
        for position in order[start:end]:
            result[position] = rank
        start = end
    return result


def pearson(first: list[float], second: list[float]) -> float:
    mean_first, mean_second = statistics.mean(first), statistics.mean(second)
    numerator = sum((x - mean_first) * (y - mean_second) for x, y in zip(first, second))
    denominator = math.sqrt(sum((x - mean_first) ** 2 for x in first) * sum((y - mean_second) ** 2 for y in second))
    return numerator / denominator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article48-run-complete").is_file():
        raise FileNotFoundError("Run run_article48_coverm.py first")
    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    genomes = {row["SGB"]: row for row in read_tsv(work / "genome-ledger.tsv")}
    samples = read_tsv(work / "samples.tsv")
    accession_to_sample = {row["RunAccession"]: row["Sample"] for row in samples}
    sample_to_accession = {row["Sample"]: row["RunAccession"] for row in samples}

    long_rows: list[dict[str, object]] = []
    capture_rows: list[dict[str, object]] = []
    by_key: dict[tuple[str, str, str], dict[str, object]] = {}
    branch_specs = (("Primary 95% identity", "identity95", 95), ("Strict 97% identity", "identity97", 97))
    for branch_label, branch_file, identity in branch_specs:
        path = work / "raw" / f"coverm-{branch_file}.tsv"
        with path.open(encoding="utf-8", newline="") as handle:
            raw_rows = list(csv.DictReader(handle, delimiter="\t"))
        for sample in ("MOCK1", "MOCK2"):
            sample_rows = [row for row in raw_rows if row["Sample"].startswith(sample_to_accession[sample])]
            if len(sample_rows) != 25:
                raise ValueError(f"{branch_label}/{sample}: expected 24 genomes plus unmapped")
            unmapped = next(row for row in sample_rows if row["Genome"] == "unmapped")
            genome_rows = [row for row in sample_rows if row["Genome"] != "unmapped"]
            observed = {row["Genome"] for row in genome_rows}
            if observed != set(genomes):
                raise ValueError(f"{branch_label}/{sample}: genome coordinate mismatch")
            for row in sorted(genome_rows, key=lambda value: value["Genome"]):
                ledger = genomes[row["Genome"]]
                mean_depth = number(row["Mean"])
                trimmed = number(row["Trimmed Mean"])
                breadth = 100 * number(row["Covered Fraction"])
                abundance = number(row["Relative Abundance (%)"])
                read_count = int(float(row["Read Count"]))
                anir = 100 * number(row["ANIr"])
                result = {
                    "Branch": branch_label,
                    "ReadIdentityThresholdPct": identity,
                    "Sample": sample,
                    "RunAccession": sample_to_accession[sample],
                    "SGB": row["Genome"],
                    "Representative": ledger["Representative"],
                    "MeanDepth": mean_depth,
                    "TrimmedMeanDepth": trimmed,
                    "CoverageUniformityPct": 100 * trimmed / mean_depth if mean_depth else math.nan,
                    "CoveredFractionPct": breadth,
                    "RelativeAbundancePct": abundance,
                    "ReadCount": read_count,
                    "MeanReadIdentityPct": anir,
                    "GenomeBp": int(float(row["Length"])),
                    "DetectedBreadth50Depth1": breadth >= 50 and mean_depth >= 1,
                    "HighBreadth90": breadth >= 90,
                }
                long_rows.append(result)
                by_key[(branch_label, sample, row["Genome"])] = result
            assigned = sum(float(row["Relative Abundance (%)"]) for row in genome_rows)
            unmapped_pct = float(unmapped["Relative Abundance (%)"])
            capture_rows.append({
                "Branch": branch_label,
                "ReadIdentityThresholdPct": identity,
                "Sample": sample,
                "RunAccession": sample_to_accession[sample],
                "CatalogRelativeAbundancePct": assigned,
                "UnmappedPct": unmapped_pct,
                "CompositionSumPct": assigned + unmapped_pct,
                "AssignedReadCount": sum(int(float(row["Read Count"])) for row in genome_rows),
                "DetectedSGBsBreadth50Depth1": sum(bool(by_key[(branch_label, sample, sgb)]["DetectedBreadth50Depth1"]) for sgb in genomes),
                "HighBreadthSGBs": sum(bool(by_key[(branch_label, sample, sgb)]["HighBreadth90"]) for sgb in genomes),
                "MedianMeanDepth": statistics.median(float(by_key[(branch_label, sample, sgb)]["MeanDepth"]) for sgb in genomes),
                "MedianCoveredFractionPct": statistics.median(float(by_key[(branch_label, sample, sgb)]["CoveredFractionPct"]) for sgb in genomes),
            })

    capture_lookup = {(row["Branch"], row["Sample"]): float(row["CatalogRelativeAbundancePct"]) for row in capture_rows}
    for row in long_rows:
        row["CatalogNormalizedAbundancePct"] = 100 * float(row["RelativeAbundancePct"]) / capture_lookup[(str(row["Branch"]), str(row["Sample"]))]

    detection_rows: list[dict[str, object]] = []
    for branch_label, _, identity in branch_specs:
        for sample in ("MOCK1", "MOCK2"):
            rows = [row for row in long_rows if row["Branch"] == branch_label and row["Sample"] == sample]
            for breadth_cutoff in (10, 50, 90):
                detection_rows.append({
                    "Branch": branch_label,
                    "ReadIdentityThresholdPct": identity,
                    "Sample": sample,
                    "BreadthCutoffPct": breadth_cutoff,
                    "MeanDepthCutoff": 1,
                    "DetectedSGBs": sum(float(row["CoveredFractionPct"]) >= breadth_cutoff and float(row["MeanDepth"]) >= 1 for row in rows),
                    "CatalogSGBs": len(rows),
                })

    sensitivity_rows: list[dict[str, object]] = []
    for sample in ("MOCK1", "MOCK2"):
        for sgb in sorted(genomes):
            primary = by_key[("Primary 95% identity", sample, sgb)]
            strict = by_key[("Strict 97% identity", sample, sgb)]
            sensitivity_rows.append({
                "Sample": sample,
                "SGB": sgb,
                "PrimaryRelativeAbundancePct": primary["RelativeAbundancePct"],
                "StrictRelativeAbundancePct": strict["RelativeAbundancePct"],
                "DeltaRelativeAbundancePctPoints": float(strict["RelativeAbundancePct"]) - float(primary["RelativeAbundancePct"]),
                "PrimaryCoveredFractionPct": primary["CoveredFractionPct"],
                "StrictCoveredFractionPct": strict["CoveredFractionPct"],
                "DeltaCoveredFractionPctPoints": float(strict["CoveredFractionPct"]) - float(primary["CoveredFractionPct"]),
                "PrimaryReadCount": primary["ReadCount"],
                "StrictReadCount": strict["ReadCount"],
                "ReadCountRetainedPct": 100 * int(strict["ReadCount"]) / int(primary["ReadCount"]),
            })

    abundance_wide: list[dict[str, object]] = []
    for sgb in sorted(genomes):
        abundance_wide.append({
            "SGB": sgb,
            "MOCK1RelativeAbundancePct": by_key[("Primary 95% identity", "MOCK1", sgb)]["RelativeAbundancePct"],
            "MOCK2RelativeAbundancePct": by_key[("Primary 95% identity", "MOCK2", sgb)]["RelativeAbundancePct"],
            "MOCK1MeanDepth": by_key[("Primary 95% identity", "MOCK1", sgb)]["MeanDepth"],
            "MOCK2MeanDepth": by_key[("Primary 95% identity", "MOCK2", sgb)]["MeanDepth"],
            "MOCK1CoveredFractionPct": by_key[("Primary 95% identity", "MOCK1", sgb)]["CoveredFractionPct"],
            "MOCK2CoveredFractionPct": by_key[("Primary 95% identity", "MOCK2", sgb)]["CoveredFractionPct"],
        })
    first = [float(row["MOCK1RelativeAbundancePct"]) for row in abundance_wide]
    second = [float(row["MOCK2RelativeAbundancePct"]) for row in abundance_wide]
    spearman = pearson(ranks(first), ranks(second))
    resources = [
        parse_time(work / "logs/coverm-identity95.time.txt"),
        parse_time(work / "logs/coverm-identity97.time.txt"),
    ]
    write_tsv(summary_dir / "coverm-long.tsv.gz", long_rows)
    write_tsv(summary_dir / "sample-capture-summary.tsv", capture_rows)
    write_tsv(summary_dir / "detection-threshold-summary.tsv", detection_rows)
    write_tsv(summary_dir / "stringency-sensitivity.tsv", sensitivity_rows)
    write_tsv(summary_dir / "abundance-wide.tsv", abundance_wide)
    write_tsv(summary_dir / "resource-summary.tsv", resources)
    main_capture = {row["Sample"]: round(float(row["CatalogRelativeAbundancePct"]), 6) for row in capture_rows if row["Branch"] == "Primary 95% identity"}
    strict_capture = {row["Sample"]: round(float(row["CatalogRelativeAbundancePct"]), 6) for row in capture_rows if row["Branch"] == "Strict 97% identity"}
    result = {
        "article": 48,
        "catalog_sgbs": 24,
        "samples": 2,
        "main_catalog_capture_pct": main_capture,
        "strict_catalog_capture_pct": strict_capture,
        "detected_sgbs_main_breadth50_depth1": {
            row["Sample"]: row["DetectedSGBsBreadth50Depth1"]
            for row in capture_rows if row["Branch"] == "Primary 95% identity"
        },
        "high_breadth_sgbs_main": {
            row["Sample"]: row["HighBreadthSGBs"]
            for row in capture_rows if row["Branch"] == "Primary 95% identity"
        },
        "mock1_mock2_relative_abundance_spearman": round(spearman, 6),
        "maximum_absolute_stringency_delta_pct_points": round(max(abs(float(row["DeltaRelativeAbundancePctPoints"])) for row in sensitivity_rows), 6),
        "truth_used_for_mapping_or_detection": False,
    }
    dump_json(summary_dir / "run-summary.json", result)
    (work / ".article48-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
