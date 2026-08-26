#!/usr/bin/env python3
"""Summarize reference-free QUAST and truth-aware MetaQUAST evidence for Article 33."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path
from typing import Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict], fields: list[str] | None = None) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    fields = fields or list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value: str | int | float | None):
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text in {"", "-", "NA", "None"}:
        return None
    parsed = float(text)
    return int(parsed) if parsed.is_integer() else parsed


def first_integer(value: str | None):
    if value is None:
        return None
    match = re.match(r"\s*([0-9,]+)", value)
    return int(match.group(1).replace(",", "")) if match else None


def parse_elapsed(value: str) -> float:
    parts = [float(part) for part in value.split(":" )]
    if len(parts) == 3:
        return parts[0] * 3600 + parts[1] * 60 + parts[2]
    if len(parts) == 2:
        return parts[0] * 60 + parts[1]
    return parts[0]


def parse_time_file(path: Path) -> dict:
    values = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ": " in raw:
            key, value = raw.strip().split(": ", 1)
            values[key] = value
    elapsed_key = next((key for key in values if key.startswith("Elapsed (wall clock) time")), None)
    elapsed = parse_elapsed(values[elapsed_key]) if elapsed_key else None
    rss = numeric(values.get("Maximum resident set size (kbytes)"))
    return {
        "Step": path.stem,
        "ElapsedSeconds": f"{elapsed:.3f}" if elapsed is not None else "NA",
        "PeakRSSGiB": f"{rss / 1024 / 1024:.3f}" if rss is not None else "NA",
        "CPUPercent": values.get("Percent of CPU this job got", "NA").rstrip("%"),
        "FileSystemInputs": numeric(values.get("File system inputs")) or 0,
        "FileSystemOutputs": numeric(values.get("File system outputs")) or 0,
        "ExitStatus": numeric(values.get("Exit status")),
    }


def rankdata(values: list[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][1] == ordered[start][1]:
            end += 1
        rank = (start + 1 + end) / 2
        for index in range(start, end):
            ranks[ordered[index][0]] = rank
        start = end
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    x_mean = statistics.mean(x)
    y_mean = statistics.mean(y)
    numerator = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y))
    denominator = math.sqrt(sum((a - x_mean) ** 2 for a in x) * sum((b - y_mean) ** 2 for b in y))
    return numerator / denominator if denominator else float("nan")


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(rankdata(x), rankdata(y))


def quast_label_map(truth_rows: list[dict[str, str]]) -> dict[str, str]:
    result = {}
    for row in truth_rows:
        label = row["Reference"]
        result[label] = label
        result[re.sub(r"[^A-Za-z0-9_.-]", "_", label)] = label
    return result


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    lineage = read_tsv(work / "summary/input-lineage.tsv")
    if len(lineage) != 17 or len({row["Branch"] for row in lineage}) != 17:
        raise ValueError("Article 33 requires 17 unique evaluation branches")
    lineage_by_branch = {row["Branch"]: row for row in lineage}
    truth = read_tsv(work / "summary/truth-manifest.tsv")
    truth_by_set = {
        set_name: [row for row in truth if row["EvaluationSet"] == set_name]
        for set_name in ("MOCK1", "MOCK2", "MOCK1+MOCK2")
    }
    expected_truth = {"MOCK1": 71, "MOCK2": 87, "MOCK1+MOCK2": 87}
    for set_name, count in expected_truth.items():
        if len(truth_by_set[set_name]) != count:
            raise ValueError(f"Unexpected truth count for {set_name}")

    reference_free_rows = read_tsv(work / "quast/transposed_report.tsv")
    reference_free = {row["Assembly"]: row for row in reference_free_rows if row["Assembly"] in lineage_by_branch}
    if set(reference_free) != set(lineage_by_branch):
        raise ValueError(f"Reference-free QUAST branches differ: {sorted(set(lineage_by_branch) - set(reference_free))}")

    meta_by_branch: dict[str, dict[str, str]] = {}
    split_rows = []
    meta_dirs = {"MOCK1": "MOCK1", "MOCK2": "MOCK2", "MOCK1+MOCK2": "MOCK1_MOCK2"}
    for set_name, safe_name in meta_dirs.items():
        group_branches = [row["Branch"] for row in lineage if row["EvaluationSet"] == set_name]
        report = read_tsv(work / f"metaquast/{safe_name}/combined_reference/transposed_report.tsv")
        report_by_name = {row["Assembly"]: row for row in report}
        missing = set(group_branches) - set(report_by_name)
        if missing:
            raise ValueError(f"Missing combined-reference rows for {set_name}: {sorted(missing)}")
        for branch in group_branches:
            meta_by_branch[branch] = report_by_name[branch]
        allowed = set(group_branches) | {f"{branch}_broken" for branch in group_branches}
        unexpected = set(report_by_name) - allowed
        if unexpected:
            raise ValueError(f"Unexpected MetaQUAST branches for {set_name}: {sorted(unexpected)}")
        for name, row in report_by_name.items():
            if name.endswith("_broken"):
                split_rows.append(
                    {
                        "EvaluationSet": set_name,
                        "Assembly": name,
                        "ParentBranch": name.removesuffix("_broken"),
                        "Contigs": numeric(row.get("# contigs")),
                        "TotalLengthBp": numeric(row.get("Total length")),
                        "N50Bp": numeric(row.get("N50")),
                        "NA50Bp": numeric(row.get("NA50")),
                        "Misassemblies": numeric(row.get("# misassemblies")),
                        "GenomeFractionPct": numeric(row.get("Genome fraction (%)")),
                    }
                )
    if split_rows:
        write_tsv(output / "split-scaffold-sensitivity.tsv", split_rows)

    per_genome = []
    physical_report_counts: dict[str, int] = {}
    for set_name, safe_name in meta_dirs.items():
        branches = [row["Branch"] for row in lineage if row["EvaluationSet"] == set_name]
        truth_rows = truth_by_set[set_name]
        label_map = quast_label_map(truth_rows)
        reports: dict[str, Path] = {}
        reports_root = work / f"metaquast/{safe_name}/runs_per_reference"
        for directory in sorted(reports_root.iterdir()):
            report = directory / "transposed_report.tsv"
            if not directory.is_dir() or not report.is_file():
                continue
            label = label_map.get(directory.name)
            if label is None:
                raise ValueError(f"Unknown per-reference report for {set_name}: {directory.name}")
            if label in reports:
                raise ValueError(f"Duplicate per-reference report for {set_name}/{label}")
            reports[label] = report
        physical_report_counts[set_name] = len(reports)
        for truth_row in truth_rows:
            label = truth_row["Reference"]
            report_rows = {row["Assembly"]: row for row in read_tsv(reports[label])} if label in reports else {}
            unexpected = set(report_rows) - set(branches) - {f"{branch}_broken" for branch in branches}
            if unexpected:
                raise ValueError(f"Unexpected branch in {set_name}/{label}: {sorted(unexpected)}")
            abundance = float(truth_row["ExpectedAbundancePct"])
            abundance_bin = "<0.1%" if abundance < 0.1 else ("0.1-<1%" if abundance < 1 else ">=1%")
            for branch in branches:
                row = report_rows.get(branch)
                fraction = numeric(row.get("Genome fraction (%)")) if row else None
                per_genome.append(
                    {
                        "EvaluationSet": set_name,
                        "Reference": label,
                        "GenBankAssembly": truth_row["GenBankAssembly"],
                        "ExpectedAbundancePct": truth_row["ExpectedAbundancePct"],
                        "AbundanceBin": abundance_bin,
                        "Branch": branch,
                        "EvidenceClass": lineage_by_branch[branch]["EvidenceClass"],
                        "GenomeFractionPct": fraction if fraction is not None else 0,
                        "RecoveredGe90Pct": "yes" if fraction is not None and fraction >= 90 else "no",
                        "FullGenomeGe99Pct": "yes" if fraction is not None and fraction >= 99 else "no",
                        "MismatchesPer100Kbp": numeric(row.get("# mismatches per 100 kbp")) if row and numeric(row.get("# mismatches per 100 kbp")) is not None else "NA",
                        "IndelsPer100Kbp": numeric(row.get("# indels per 100 kbp")) if row and numeric(row.get("# indels per 100 kbp")) is not None else "NA",
                        "ReferenceReportPresent": "yes" if label in reports else "no",
                        "BranchReportPresent": "yes" if row is not None else "no",
                    }
                )
    expected_per_genome = 13 * 71 + 2 * 87 + 2 * 87
    if len(per_genome) != expected_per_genome:
        raise ValueError(f"Expected {expected_per_genome} per-genome rows, found {len(per_genome)}")
    write_tsv(output / "per-genome-metaquast.tsv", per_genome)

    branch_metrics = []
    for local in lineage:
        branch = local["Branch"]
        row = meta_by_branch[branch]
        genome_rows = [item for item in per_genome if item["Branch"] == branch]
        truth_count = len(genome_rows)
        recovered = sum(item["RecoveredGe90Pct"] == "yes" for item in genome_rows)
        full = sum(item["FullGenomeGe99Pct"] == "yes" for item in genome_rows)
        branch_metrics.append(
            {
                "Branch": branch,
                "Display": local["Display"],
                "EvaluationSet": local["EvaluationSet"],
                "EvidenceClass": local["EvidenceClass"],
                "Family": local["Family"],
                "Assembler": local["Assembler"],
                "Strategy": local["Strategy"],
                "TruthGenomes": truth_count,
                "ContigsGe1kb": numeric(local["Contigs"]),
                "TotalLengthBp": numeric(local["TotalBp"]),
                "LargestBp": numeric(local["LargestBp"]),
                "N50Bp": numeric(local["N50Bp"]),
                "L50": numeric(local["L50"]),
                "N90Bp": numeric(local["N90Bp"]),
                "L90": numeric(local["L90"]),
                "GCPercent": numeric(local["GCPercent"]),
                "NBases": numeric(local["NBases"]),
                "Misassemblies": numeric(row.get("# misassemblies")),
                "MisassembledContigs": numeric(row.get("# misassembled contigs")),
                "LocalMisassemblies": numeric(row.get("# local misassemblies")),
                "FullyUnalignedContigs": first_integer(row.get("# unaligned contigs")),
                "UnalignedLengthBp": numeric(row.get("Unaligned length")),
                "GenomeFractionPct": numeric(row.get("Genome fraction (%)")),
                "DuplicationRatio": numeric(row.get("Duplication ratio")),
                "NsPer100Kbp": numeric(row.get("# N's per 100 kbp")),
                "MismatchesPer100Kbp": numeric(row.get("# mismatches per 100 kbp")),
                "IndelsPer100Kbp": numeric(row.get("# indels per 100 kbp")),
                "LargestAlignmentBp": numeric(row.get("Largest alignment")),
                "TotalAlignedLengthBp": numeric(row.get("Total aligned length")),
                "NA50Bp": numeric(row.get("NA50")),
                "LA50": numeric(row.get("LA50")),
                "NGA50Bp": numeric(row.get("NGA50")),
                "LGA50": numeric(row.get("LGA50")),
                "RecoveredGenomesGe90Pct": recovered,
                "RecoveredFractionGe90Pct": f"{recovered / truth_count:.9f}",
                "FullGenomesGe99Pct": full,
                "FullGenomeFractionGe99Pct": f"{full / truth_count:.9f}",
                "CanonicalSHA256": local["CanonicalSHA256"],
            }
        )
    write_tsv(output / "branch-metrics.tsv", branch_metrics)

    abundance_rows = []
    for branch in lineage_by_branch:
        for abundance_bin in ("<0.1%", "0.1-<1%", ">=1%"):
            rows = [row for row in per_genome if row["Branch"] == branch and row["AbundanceBin"] == abundance_bin]
            if not rows:
                continue
            abundance_rows.append(
                {
                    "Branch": branch,
                    "EvaluationSet": lineage_by_branch[branch]["EvaluationSet"],
                    "AbundanceBin": abundance_bin,
                    "TruthGenomes": len(rows),
                    "MedianGenomeFractionPct": f"{statistics.median(float(row['GenomeFractionPct']) for row in rows):.6f}",
                    "RecoveredGenomesGe90Pct": sum(row["RecoveredGe90Pct"] == "yes" for row in rows),
                    "FullGenomesGe99Pct": sum(row["FullGenomeGe99Pct"] == "yes" for row in rows),
                }
            )
    write_tsv(output / "abundance-bin-recovery.tsv", abundance_rows)

    threshold_rows = []
    for row in lineage:
        total = float(row["TotalBp"])
        threshold_rows.append(
            {
                "Branch": row["Branch"],
                "EvidenceClass": row["EvidenceClass"],
                "BasesGe1kb": row["TotalBp"],
                "BasesGe10kb": row["BasesGe10kb"],
                "BasesGe100kb": row["BasesGe100kb"],
                "RetainedAt10kbPct": f"{100 * float(row['BasesGe10kb']) / total:.6f}",
                "RetainedAt100kbPct": f"{100 * float(row['BasesGe100kb']) / total:.6f}",
                "ContigsGe1kb": row["Contigs"],
                "ContigsGe10kb": row["ContigsGe10kb"],
                "ContigsGe100kb": row["ContigsGe100kb"],
            }
        )
    write_tsv(output / "length-threshold-sensitivity.tsv", threshold_rows)

    m1_biological = [row for row in branch_metrics if row["EvaluationSet"] == "MOCK1" and row["EvidenceClass"] == "biological"]
    correlation_specs = (
        ("N50Bp", "NA50Bp"),
        ("N50Bp", "GenomeFractionPct"),
        ("N50Bp", "RecoveredFractionGe90Pct"),
        ("N50Bp", "Misassemblies"),
        ("N50Bp", "MismatchesPer100Kbp"),
        ("N50Bp", "IndelsPer100Kbp"),
    )
    correlations = []
    for x_name, y_name in correlation_specs:
        pairs = [
            (float(row[x_name]), float(row[y_name]))
            for row in m1_biological
            if row[x_name] not in (None, "NA") and row[y_name] not in (None, "NA")
        ]
        correlations.append(
            {
                "Scope": "MOCK1 biological workflows",
                "X": x_name,
                "Y": y_name,
                "N": len(pairs),
                "SpearmanRho": f"{spearman([x for x, _ in pairs], [y for _, y in pairs]):.9f}",
                "Interpretation": "descriptive cross-workflow association; not a causal platform comparison",
            }
        )
    write_tsv(output / "metric-correlation-audit.tsv", correlations)

    by_branch = {row["Branch"]: row for row in branch_metrics}
    source = by_branch["lr-flye-hifi"]
    control_rows = []
    for branch in ("diagnostic-fragmented-50kb", "diagnostic-chimeric-rotation"):
        row = by_branch[branch]
        control_rows.append(
            {
                "Control": branch,
                "Source": "lr-flye-hifi",
                "TotalLengthDeltaBp": int(row["TotalLengthBp"]) - int(source["TotalLengthBp"]),
                "N50SourceBp": source["N50Bp"],
                "N50ControlBp": row["N50Bp"],
                "N50DeltaBp": int(row["N50Bp"]) - int(source["N50Bp"]),
                "NA50SourceBp": source["NA50Bp"],
                "NA50ControlBp": row["NA50Bp"],
                "MisassembliesSource": source["Misassemblies"],
                "MisassembliesControl": row["Misassemblies"],
                "MisassemblyDelta": int(row["Misassemblies"]) - int(source["Misassemblies"]),
                "GenomeFractionSourcePct": source["GenomeFractionPct"],
                "GenomeFractionControlPct": row["GenomeFractionPct"],
                "MismatchesSourcePer100Kbp": source["MismatchesPer100Kbp"],
                "MismatchesControlPer100Kbp": row["MismatchesPer100Kbp"],
                "Purpose": "diagnostic transformation of a real assembly; not a biological workflow",
            }
        )
    write_tsv(output / "diagnostic-control-effects.tsv", control_rows)

    resources = [parse_time_file(path) for path in sorted((work / "resources").glob("*.txt"))]
    if resources:
        write_tsv(output / "resource-usage.tsv", resources)

    prepare_summary = json.loads((work / "summary/prepare-summary.json").read_text())
    summary = {
        **prepare_summary,
        "per_genome_rows": len(per_genome),
        "physical_reference_reports": physical_report_counts,
        "physical_reference_reports_total": sum(physical_report_counts.values()),
        "branch_metrics_rows": len(branch_metrics),
        "mock1_biological_branches": len(m1_biological),
        "quast_version": "5.3.0",
        "minimum_contig_bp": 1000,
        "minimum_alignment_bp": 500,
        "minimum_identity_pct": 97,
        "fragmented_reference_mode": True,
        "split_scaffolds_sensitivity": True,
        "n50_is_correctness": False,
        "universal_assembly_threshold_claimed": False,
        "diagnostic_controls_are_biological_results": False,
        "chimeric_control_n50_invariant": by_branch["diagnostic-chimeric-rotation"]["N50Bp"] == source["N50Bp"],
        "chimeric_control_misassembly_delta": control_rows[1]["MisassemblyDelta"],
        "fragmented_control_n50_delta_bp": control_rows[0]["N50DeltaBp"],
    }
    (output / "run-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
