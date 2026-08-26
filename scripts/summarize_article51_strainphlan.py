#!/usr/bin/env python3
"""Summarize StrainPhlAn filtering, diversity, topology and baseline concordance."""

from __future__ import annotations

import argparse
import copy
import json
import math
import re
import statistics
from pathlib import Path

from Bio import Phylo, SeqIO

from article41_44_utils import dump_json, read_tsv, sha256, write_tsv


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def info_integer(text: str, patterns: tuple[str, ...], label: str) -> int:
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.MULTILINE)
        if match:
            return int(match.group(1))
    raise ValueError(f"Could not parse {label} from StrainPhlAn info file")


def parse_info_values(text: str) -> dict[str, int]:
    return {
        "InputSamples": info_integer(
            text,
            (r"^Number of samples:\s*(\d+)$", r"^Number of main samples:\s*(\d+)$"),
            "input samples",
        ),
        "InputReferences": info_integer(
            text,
            (r"^Number of references:\s*(\d+)$", r"^Number of main references:\s*(\d+)$"),
            "input references",
        ),
        "AvailableMarkers": info_integer(
            text,
            (r"^Number of available markers for the clade:\s*(\d+)$",),
            "available markers",
        ),
        "SelectedMarkers": info_integer(
            text,
            (r"^Number of markers selected after filtering:\s*(\d+)$",),
            "selected markers",
        ),
        "RetainedSamples": info_integer(
            text,
            (r"^Number of samples after filtering:\s*(\d+)$", r"^Number of main samples after filtering:\s*(\d+)$"),
            "retained samples",
        ),
        "RetainedReferences": info_integer(
            text,
            (r"^Number of references after filtering:\s*(\d+)$", r"^Number of main references after filtering:\s*(\d+)$"),
            "retained references",
        ),
    }


def read_alignment(path: Path) -> dict[str, str]:
    sequences = {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}
    if not sequences:
        raise ValueError(f"Empty alignment: {path}")
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(f"Non-rectangular alignment: {sorted(lengths)}")
    return sequences


def pair_distance(first: str, second: str) -> tuple[int, int, float]:
    allowed = set("ACGT")
    comparable = 0
    differences = 0
    for left, right in zip(first, second):
        if left in allowed and right in allowed:
            comparable += 1
            differences += left != right
    distance = differences / comparable if comparable else math.nan
    return comparable, differences, distance


def prune_to(tree: Phylo.BaseTree.Tree, keep: set[str]) -> Phylo.BaseTree.Tree:
    pruned = copy.deepcopy(tree)
    for terminal in list(pruned.get_terminals()):
        if terminal.name not in keep:
            pruned.prune(terminal)
    return pruned


def unrooted_splits(tree: Phylo.BaseTree.Tree) -> set[tuple[str, ...]]:
    all_tips = {terminal.name for terminal in tree.get_terminals()}
    splits: set[tuple[str, ...]] = set()
    for clade in tree.get_nonterminals(order="preorder"):
        side = {terminal.name for terminal in clade.get_terminals()}
        other = all_tips - side
        if min(len(side), len(other)) < 2:
            continue
        left = tuple(sorted(side))
        right = tuple(sorted(other))
        splits.add(left if (len(left), left) <= (len(right), right) else right)
    return splits


def median(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return statistics.median(finite) if finite else math.nan


def topology_row(
    comparison: str,
    branch: str,
    actual_tree: Phylo.BaseTree.Tree,
    actual_path: Path,
    baseline_tree: Phylo.BaseTree.Tree,
    baseline_path: Path,
) -> dict[str, object]:
    actual_tips = {terminal.name for terminal in actual_tree.get_terminals()}
    baseline_tips = {terminal.name for terminal in baseline_tree.get_terminals()}
    common_tips = actual_tips & baseline_tips
    if len(common_tips) < 4:
        raise ValueError(f"Too few shared tips for topology comparison: {branch}")
    actual_splits = unrooted_splits(prune_to(actual_tree, common_tips))
    baseline_splits = unrooted_splits(prune_to(baseline_tree, common_tips))
    rf_distance = len(actual_splits.symmetric_difference(baseline_splits))
    rf_denominator = len(actual_splits) + len(baseline_splits)
    return {
        "Comparison": comparison,
        "Branch": branch,
        "ActualTips": len(actual_tips),
        "BaselineTips": len(baseline_tips),
        "CommonTips": len(common_tips),
        "ActualSplits": len(actual_splits),
        "BaselineSplits": len(baseline_splits),
        "RobinsonFouldsDistance": rf_distance,
        "NormalizedRF": round(rf_distance / rf_denominator, 6) if rf_denominator else math.nan,
        "ExactTopologyMatch": rf_distance == 0,
        "BaselineTreeSHA256": sha256(baseline_path),
        "ActualTreeSHA256": sha256(actual_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article51-run-complete").is_file():
        raise FileNotFoundError("Run run_article51_strainphlan.py first")
    contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))
    output_paths = read_tsv(work / "output-paths.tsv")[0]
    info_path = Path(output_paths["Info"])
    polymorphic_path = Path(output_paths["Polymorphic"])
    alignment_path = Path(output_paths["Alignment"])
    tree_path = Path(output_paths["Tree"])
    official_threshold_info_path = Path(output_paths["OfficialThresholdInfo"])
    official_threshold_alignment_path = Path(output_paths["OfficialThresholdAlignment"])
    official_threshold_tree_path = Path(output_paths["OfficialThresholdTree"])
    metadata = {row["Sample"]: row for row in read_tsv(work / "sample-metadata.tsv")}

    info_values = parse_info_values(info_path.read_text(encoding="utf-8", errors="replace"))
    official_threshold_info_values = parse_info_values(
        official_threshold_info_path.read_text(encoding="utf-8", errors="replace")
    )

    sequences = read_alignment(alignment_path)
    alignment_sites = len(next(iter(sequences.values())))
    tree = Phylo.read(tree_path, "newick")
    tree_tips = {terminal.name for terminal in tree.get_terminals()}
    if tree_tips != set(sequences):
        raise ValueError("Tree and concatenated alignment contain different tips")
    sample_tips = sorted(set(metadata) & tree_tips)
    reference_tips = sorted(tree_tips - set(metadata))
    if len(sample_tips) != info_values["RetainedSamples"]:
        raise ValueError("Retained sample count disagrees with alignment")
    if len(reference_tips) != info_values["RetainedReferences"]:
        raise ValueError("Retained reference count disagrees with alignment")

    official_threshold_sequences = read_alignment(official_threshold_alignment_path)
    official_threshold_alignment_sites = len(next(iter(official_threshold_sequences.values())))
    official_threshold_tree = Phylo.read(official_threshold_tree_path, "newick")
    official_threshold_tree_tips = {
        terminal.name for terminal in official_threshold_tree.get_terminals()
    }
    if official_threshold_tree_tips != set(official_threshold_sequences):
        raise ValueError("Official-threshold tree and alignment contain different tips")
    official_threshold_sample_tips = set(metadata) & official_threshold_tree_tips
    official_threshold_reference_tips = official_threshold_tree_tips - set(metadata)
    if len(official_threshold_sample_tips) != official_threshold_info_values["RetainedSamples"]:
        raise ValueError("Official-threshold sample count disagrees with alignment")
    if len(official_threshold_reference_tips) != official_threshold_info_values["RetainedReferences"]:
        raise ValueError("Official-threshold reference count disagrees with alignment")

    tip_rows: list[dict[str, object]] = []
    for tip in sorted(tree_tips):
        observed = sum(base in "ACGT" for base in sequences[tip])
        meta = metadata.get(tip, {})
        tip_rows.append(
            {
                "Tip": tip,
                "Type": "Metagenome sample" if tip in metadata else "Reference genome",
                "Study": meta.get("Study", "Reference"),
                "Country": meta.get("Country", "Reference"),
                "AlignmentSites": alignment_sites,
                "ObservedACGTSites": observed,
                "OccupancyPct": round(100 * observed / alignment_sites, 6),
            }
        )

    distance_rows: list[dict[str, object]] = []
    for index, first in enumerate(sample_tips):
        for second in sample_tips[index + 1 :]:
            comparable, differences, distance = pair_distance(sequences[first], sequences[second])
            first_meta = metadata[first]
            second_meta = metadata[second]
            if first_meta["Study"] == second_meta["Study"]:
                stratum = "Same study"
            elif first_meta["Country"] == second_meta["Country"]:
                stratum = "Same country, different study"
            else:
                stratum = "Different country"
            distance_rows.append(
                {
                    "Sample1": first,
                    "Sample2": second,
                    "Study1": first_meta["Study"],
                    "Study2": second_meta["Study"],
                    "Country1": first_meta["Country"],
                    "Country2": second_meta["Country"],
                    "PairStratum": stratum,
                    "ComparableSites": comparable,
                    "SNPs": differences,
                    "PDistance": round(distance, 10),
                    "DifferencesPer10kb": round(distance * 10_000, 6),
                }
            )

    nearest_rows: list[dict[str, object]] = []
    for sample in sample_tips:
        candidates = []
        for row in distance_rows:
            if row["Sample1"] == sample:
                candidates.append((float(row["PDistance"]), str(row["Sample2"])))
            elif row["Sample2"] == sample:
                candidates.append((float(row["PDistance"]), str(row["Sample1"])))
        distance, neighbor = min(candidates, key=lambda value: (value[0], value[1]))
        nearest_rows.append(
            {
                "Sample": sample,
                "NearestSample": neighbor,
                "PDistance": distance,
                "SameStudy": metadata[sample]["Study"] == metadata[neighbor]["Study"],
                "SameCountry": metadata[sample]["Country"] == metadata[neighbor]["Country"],
                "Study": metadata[sample]["Study"],
                "Country": metadata[sample]["Country"],
            }
        )

    polymorphic_rows: list[dict[str, object]] = []
    for row in read_tsv(polymorphic_path):
        sample = row["sample"]
        if sample not in metadata:
            raise ValueError(f"Unknown sample in polymorphism table: {sample}")
        polymorphic_rows.append(
            {
                **row,
                "Study": metadata[sample]["Study"],
                "Country": metadata[sample]["Country"],
            }
        )

    baseline_tree_path = work / "official_baseline/RAxML_bestTree.t__SGB4933_group.StrainPhlAn4.tre"
    baseline_tree = Phylo.read(baseline_tree_path, "newick")
    topology_rows = [
        topology_row(
            "Current explicit thresholds vs official tutorial baseline",
            "Current explicit thresholds",
            tree,
            tree_path,
            baseline_tree,
            baseline_tree_path,
        ),
        topology_row(
            "Official 2022 thresholds rerun vs official tutorial baseline",
            "Official 2022 thresholds",
            official_threshold_tree,
            official_threshold_tree_path,
            baseline_tree,
            baseline_tree_path,
        ),
    ]

    current_thresholds = {
        "sample_with_n_markers": contract["sample_with_n_markers"],
        "sample_with_n_markers_perc": contract["sample_with_n_markers_perc"],
        "marker_in_n_samples_perc": contract["marker_in_n_samples_perc"],
        "sample_with_n_markers_after_filt": contract["sample_with_n_markers_after_filt"],
        "sample_with_n_markers_after_filt_perc": contract["sample_with_n_markers_after_filt_perc"],
        "breadth_thres": contract["breadth_thres"],
    }
    official_thresholds = contract["official_baseline_threshold_sensitivity"]
    threshold_rows = []
    for branch, thresholds, values, sites, tips in (
        ("Current explicit thresholds", current_thresholds, info_values, alignment_sites, len(tree_tips)),
        (
            "Official 2022 thresholds",
            official_thresholds,
            official_threshold_info_values,
            official_threshold_alignment_sites,
            len(official_threshold_tree_tips),
        ),
    ):
        threshold_rows.append(
            {
                "Branch": branch,
                **thresholds,
                **values,
                "AlignmentSites": sites,
                "TreeTips": tips,
            }
        )

    filtering_rows = [
        {"Metric": key, "Value": value}
        for key, value in info_values.items()
    ] + [
        {"Metric": "AlignmentSites", "Value": alignment_sites},
        {"Metric": "TreeTips", "Value": len(tree_tips)},
        {"Metric": "RAxMLSeed", "Value": int(contract["seed"])},
    ]
    stratum_medians = {
        stratum: median(
            [float(row["PDistance"]) for row in distance_rows if row["PairStratum"] == stratum]
        )
        for stratum in ("Same study", "Same country, different study", "Different country")
    }
    summary = {
        "article": 51,
        "input_samples": info_values["InputSamples"],
        "input_references": info_values["InputReferences"],
        "available_markers": info_values["AvailableMarkers"],
        "selected_markers": info_values["SelectedMarkers"],
        "retained_samples": info_values["RetainedSamples"],
        "retained_references": info_values["RetainedReferences"],
        "alignment_sites": alignment_sites,
        "sample_pairs": len(distance_rows),
        "median_pairwise_p_distance": median([float(row["PDistance"]) for row in distance_rows]),
        "median_pairwise_p_distance_by_stratum": stratum_medians,
        "median_polymorphic_sites_pct": median(
            [as_float(str(row["percentage_of_polymorphic_sites"])) for row in polymorphic_rows]
        ),
        "nearest_neighbor_same_study": sum(str(row["SameStudy"]).lower() == "true" for row in nearest_rows),
        "nearest_neighbor_same_country": sum(str(row["SameCountry"]).lower() == "true" for row in nearest_rows),
        "topology_common_tips": topology_rows[0]["CommonTips"],
        "normalized_rf_vs_official_baseline": topology_rows[0]["NormalizedRF"],
        "exact_topology_match_vs_official_baseline": topology_rows[0]["ExactTopologyMatch"],
        "official_thresholds_selected_markers": official_threshold_info_values["SelectedMarkers"],
        "official_thresholds_alignment_sites": official_threshold_alignment_sites,
        "official_thresholds_topology_common_tips": topology_rows[1]["CommonTips"],
        "official_thresholds_normalized_rf_vs_baseline": topology_rows[1]["NormalizedRF"],
        "official_thresholds_exact_topology_match_vs_baseline": topology_rows[1]["ExactTopologyMatch"],
        "raxml_seed": int(contract["seed"]),
        "geography_used_for_tree_inference": False,
    }

    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(summary_dir / "filtering-summary.tsv", filtering_rows)
    write_tsv(summary_dir / "threshold-branch-summary.tsv", threshold_rows)
    write_tsv(summary_dir / "tip-metadata.tsv", tip_rows)
    write_tsv(summary_dir / "pairwise-p-distance.tsv", distance_rows)
    write_tsv(summary_dir / "nearest-neighbor-audit.tsv", nearest_rows)
    write_tsv(summary_dir / "polymorphism-by-sample.tsv", polymorphic_rows)
    write_tsv(summary_dir / "topology-baseline-audit.tsv", topology_rows)
    dump_json(summary_dir / "run-summary.json", summary)
    (work / ".article51-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
