#!/usr/bin/env python3
"""Summarize genome-level inStrain microdiversity and two-sample popANI evidence."""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics
from pathlib import Path

from article41_44_utils import dump_json, read_tsv, write_tsv


def number(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def integer(value: object) -> int:
    parsed = number(value)
    return int(round(parsed)) if math.isfinite(parsed) else 0


def one_file(pattern: str, base: Path) -> Path:
    matches = sorted(base.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"Expected one {pattern} under {base}, observed {matches}")
    return matches[0]


def parse_read_filter(log: Path, sample: str) -> dict[str, object]:
    text = log.read_text(encoding="utf-8", errors="replace")
    removed = re.search(
        r"([0-9.]+)% of pairs and ([0-9.]+)% of singletons were removed during filtering",
        text,
    )
    retained = re.search(r"([0-9,]+) read pairs remain \(([0-9.]+) Gbp\)", text)
    if removed is None or retained is None:
        raise ValueError(f"Could not parse the inStrain read-filter audit from {log}")
    return {
        "Sample": sample,
        "ReadPairsRemovedPct": float(removed.group(1)),
        "SingletonsRemovedPct": float(removed.group(2)),
        "ReadPairsRetained": int(retained.group(1).replace(",", "")),
        "RetainedGbp": float(retained.group(2)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article50-run-complete").is_file():
        raise FileNotFoundError("Run run_article50_instrain.py first")

    summary = work / "summary"
    summary.mkdir(parents=True, exist_ok=True)
    contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))
    genome_ledger = {row["SGB"]: row for row in read_tsv(work / "genome-ledger.tsv")}
    expected = set(genome_ledger)

    profile_rows: list[dict[str, object]] = []
    profile_lookup: dict[tuple[str, str], dict[str, object]] = {}
    filter_rows: list[dict[str, object]] = []
    for sample in contract["samples"]:
        profile_file = one_file("*_genome_info.tsv*", work / "profiles" / sample / "output")
        rows = read_tsv(profile_file)
        observed = {row["genome"] for row in rows}
        if observed != expected:
            raise ValueError(
                f"{sample}: genome_info coordinate mismatch; missing={sorted(expected-observed)}, "
                f"extra={sorted(observed-expected)}"
            )
        for row in sorted(rows, key=lambda item: item["genome"]):
            sgb = row["genome"]
            length = integer(row.get("length"))
            breadth_min_cov = number(row.get("breadth_minCov"))
            considered = length * breadth_min_cov if math.isfinite(breadth_min_cov) else 0.0
            snv_count = integer(row.get("SNV_count"))
            raw_rarefied = number(row.get("nucl_diversity_rarefied"))
            rarefied_enabled = bool(contract["rarefied_metric_enabled"])
            entry = {
                "Sample": sample,
                "RunAccession": {
                    "MOCK1": "ERR9765746",
                    "MOCK2": "ERR9765747",
                }[sample],
                "SGB": sgb,
                "GTDBTaxonomy": genome_ledger[sgb]["GTDBTaxonomy"],
                "GenomeBp": length,
                "Coverage": number(row.get("coverage")),
                "BreadthPct": 100 * number(row.get("breadth")),
                "BreadthAtMinCovPct": 100 * breadth_min_cov,
                "MinCov": int(contract["min_cov"]),
                "NucleotideDiversity": number(row.get("nucl_diversity")),
                "RarefiedNucleotideDiversity": (
                    raw_rarefied if rarefied_enabled else math.nan
                ),
                "RawRarefiedNucleotideDiversity": raw_rarefied,
                "RarefiedMetricStatus": (
                    "ENABLED"
                    if rarefied_enabled
                    else "DISABLED_UNREACHED_COVERAGE_SENTINEL_NORMALIZED_TO_NA"
                ),
                "ConANIToReferencePct": 100 * number(row.get("conANI_reference")),
                "PopANIToReferencePct": 100 * number(row.get("popANI_reference")),
                "SNVCount": snv_count,
                "SNSCount": integer(row.get("SNS_count")),
                "ConsensusDivergentSites": integer(row.get("consensus_divergent_sites")),
                "PopulationDivergentSites": integer(row.get("population_divergent_sites")),
                "ConsideredBpAtMinCov": round(considered),
                "SNVsPerMbpConsidered": 1e6 * snv_count / considered if considered > 0 else math.nan,
                "PresentBreadth50AtMinCov": breadth_min_cov >= float(contract["compare_presence_breadth"]),
                "FilteredReadPairCount": integer(row.get("filtered_read_pair_count")),
                "ReadsMeanPIDPct": 100 * number(row.get("reads_mean_PID")),
            }
            profile_rows.append(entry)
            profile_lookup[(sample, sgb)] = entry
        filter_rows.append(
            parse_read_filter(work / "logs" / f"instrain-profile-{sample.lower()}.stderr.log", sample)
        )

    compare_file = one_file(
        "*_genomeWide_compare.tsv*", work / "comparison/MOCK1-vs-MOCK2/output"
    )
    observed_compare: dict[str, dict[str, object]] = {}
    for row in read_tsv(compare_file):
        if row.get("name1") == row.get("name2"):
            continue
        sgb = row["genome"]
        if sgb in observed_compare:
            raise ValueError(f"Duplicate genome-wide comparison for {sgb}")
        observed_compare[sgb] = row

    popani_min = float(contract["same_strain_reporting_rule"]["popANI_min_pct"])
    overlap_min = float(
        contract["same_strain_reporting_rule"]["percent_genome_compared_min_pct"]
    )
    comparison_rows: list[dict[str, object]] = []
    for sgb in sorted(expected):
        row = observed_compare.get(sgb)
        present_first = bool(profile_lookup[("MOCK1", sgb)]["PresentBreadth50AtMinCov"])
        present_second = bool(profile_lookup[("MOCK2", sgb)]["PresentBreadth50AtMinCov"])
        popani = 100 * number(row.get("popANI")) if row else math.nan
        conani = 100 * number(row.get("conANI")) if row else math.nan
        overlap = 100 * number(row.get("percent_compared")) if row else math.nan
        comparison_rows.append(
            {
                "SGB": sgb,
                "GTDBTaxonomy": genome_ledger[sgb]["GTDBTaxonomy"],
                "Sample1": "MOCK1",
                "Sample2": "MOCK2",
                "NativeSample1": row.get("name1", "") if row else "",
                "NativeSample2": row.get("name2", "") if row else "",
                "PresentInMOCK1": present_first,
                "PresentInMOCK2": present_second,
                "Compared": row is not None,
                "ComparedBases": integer(row.get("compared_bases_count")) if row else 0,
                "PercentGenomeCompared": overlap,
                "PopulationSNPs": integer(row.get("population_SNPs")) if row else 0,
                "ConsensusSNPs": integer(row.get("consensus_SNPs")) if row else 0,
                "PopANIPct": popani,
                "ConANIPct": conani,
                "SameStrainRule": bool(
                    row is not None
                    and math.isfinite(popani)
                    and math.isfinite(overlap)
                    and popani >= popani_min
                    and overlap >= overlap_min
                ),
                "ReportingPopANIThresholdPct": popani_min,
                "ReportingGenomeComparedThresholdPct": overlap_min,
            }
        )

    write_tsv(summary / "profile-genome-long.tsv.gz", profile_rows)
    write_tsv(summary / "pairwise-genome-comparison.tsv", comparison_rows)
    write_tsv(summary / "read-filter-audit.tsv", filter_rows)

    compared = [row for row in comparison_rows if row["Compared"]]
    same = [row for row in compared if row["SameStrainRule"]]
    diversity_summary: dict[str, dict[str, float]] = {}
    for sample in contract["samples"]:
        values = [
            float(profile_lookup[(sample, sgb)]["NucleotideDiversity"])
            for sgb in expected
            if profile_lookup[(sample, sgb)]["PresentBreadth50AtMinCov"]
        ]
        diversity_summary[sample] = {
            "min": min(values),
            "median": statistics.median(values),
            "max": max(values),
        }
    run_summary = {
        "article": 50,
        "catalog_sgbs": len(expected),
        "samples": len(contract["samples"]),
        "min_read_ani": contract["min_read_ani"],
        "min_cov": contract["min_cov"],
        "min_freq": contract["min_freq"],
        "presence_breadth_min_pct": 100 * contract["compare_presence_breadth"],
        "present_sgbs": {
            sample: sum(
                bool(profile_lookup[(sample, sgb)]["PresentBreadth50AtMinCov"])
                for sgb in expected
            )
            for sample in contract["samples"]
        },
        "nucleotide_diversity_present_sgbs": diversity_summary,
        "read_pairs_retained": {
            row["Sample"]: row["ReadPairsRetained"] for row in filter_rows
        },
        "genomes_compared": len(compared),
        "same_strain_calls": len(same),
        "same_strain_sgbs": [row["SGB"] for row in same],
        "popani_range_pct": [
            min(float(row["PopANIPct"]) for row in compared),
            max(float(row["PopANIPct"]) for row in compared),
        ] if compared else [],
        "percent_genome_compared_range": [
            min(float(row["PercentGenomeCompared"]) for row in compared),
            max(float(row["PercentGenomeCompared"]) for row in compared),
        ] if compared else [],
        "truth_used_for_profiling_or_comparison": False,
    }
    dump_json(summary / "run-summary.json", run_summary)
    (work / ".article50-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(run_summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
