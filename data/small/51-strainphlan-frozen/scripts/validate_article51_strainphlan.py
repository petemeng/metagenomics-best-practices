#!/usr/bin/env python3
"""Fail-closed validation for Article 51 StrainPhlAn evidence."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

from article42_44_validation_utils import (
    Audit,
    as_bool,
    audit_checksums,
    audit_figures,
    finish,
    read_tsv,
)


FIGURES = (
    "51-strainphlan-tree",
    "51-marker-distance-pcoa",
    "51-pairwise-distance-strata",
    "51-polymorphism-audit",
)


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot = (frozen / "scripts/plot_article51_strainphlan.R").read_text(encoding="utf-8")
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260751" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "耗时")),
        "official-real-data": all(token in text for token in ("25", "13 studies", "9 countries", "SGB4933")),
        "frozen-input": "51-strainphlan-frozen" in text,
        "methods-results": all(token in text for token in ("Methods template", "Results template")),
        "versions": all(token in text for token in ("StrainPhlAn 4.2.5", "MetaPhlAn 4.2.5", "RAxML binary 8.2.12", "Bioconda package 8.2.13")),
        "database": "mpa_vJan21_CHOCOPhlAnSGB_202103" in text,
        "filters": all(token in text for token in (
            "--trim_sequences 50", "--sample_with_n_markers 20",
            "--sample_with_n_markers_perc 25", "--marker_in_n_samples_perc 50",
            "--sample_with_n_markers_after_filt 20",
            "--sample_with_n_markers_after_filt_perc 25", "--breadth_thres 80",
        )),
        "official-threshold-sensitivity": all(token in text for token in (
            "--sample_with_n_markers_perc 80", "--marker_in_n_samples_perc 80",
            "--sample_with_n_markers_after_filt_perc 80",
        )),
        "seeded-raxml": "-p 20260751" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot) is None,
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-meta-prose": not any(token in text for token in (
            "本篇可独立", "本文可独立", "全系列约定", "接口只学一次",
            "作者代码通常长这样", "（即本文）", "无头服务器",
        )),
    }
    sections = (
        "这一步对应论文里的哪张图", "理论：", "准备工作", "可复制代码", "审计与升级",
        "出版级美化", "常见坑", "这段 Methods 怎么写", "换成你自己的数据怎么做", "参考",
    )
    for section in sections:
        checks[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def finite(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    args = parser.parse_args()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_checksums(frozen, audit)

    summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    contract = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    assets = read_tsv(frozen / "asset-manifest.tsv")
    asset_checks = read_tsv(frozen / "asset-check-audit.tsv")
    metadata = read_tsv(frozen / "sample-metadata.tsv")
    filtering = {row["Metric"]: int(row["Value"]) for row in read_tsv(frozen / "filtering-summary.tsv")}
    threshold_rows = read_tsv(frozen / "threshold-branch-summary.tsv")
    tips = read_tsv(frozen / "tip-metadata.tsv")
    pairs = read_tsv(frozen / "pairwise-p-distance.tsv")
    nearest = read_tsv(frozen / "nearest-neighbor-audit.tsv")
    polymorphism = read_tsv(frozen / "polymorphism-by-sample.tsv")
    topology = read_tsv(frozen / "topology-baseline-audit.tsv")
    determinism = read_tsv(frozen / "determinism-audit.tsv")
    commands = read_tsv(frozen / "command-log.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    versions = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}

    audit.add("Identity", "article", summary.get("article") == 51 and contract.get("article") == 51, summary)
    audit.add("Identity", "seed", contract.get("seed") == 20260751 and summary.get("raxml_seed") == 20260751, contract)
    audit.add("Identity", "target", contract.get("target_clade") == "t__SGB4933_group", contract.get("target_clade"))
    audit.add("Identity", "database", contract.get("database_release") == "mpa_vJan21_CHOCOPhlAnSGB_202103", contract.get("database_release"))
    audit.add("Identity", "truth-blind", contract.get("truth_used_for_tree") is False and summary.get("geography_used_for_tree_inference") is False, summary)

    counts = Counter(row["AssetType"] for row in assets)
    audit.add("Input", "37-assets", len(assets) == 37 and counts == Counter({
        "consensus marker profile": 25,
        "clade marker FASTA": 1,
        "reference genome": 6,
        "official precomputed output": 4,
        "MetaPhlAn metadata extracted by HTTP range": 1,
    }), counts)
    audit.add("Input", "asset-sha256", all(re.fullmatch(r"[0-9a-f]{64}", row["SHA256"]) for row in assets), len(assets))
    audit.add("Input", "asset-checks", len(asset_checks) == 37 and all(as_bool(row["ChecksumPass"]) and row["ExpectedSHA256"] == row["ObservedSHA256"] and row["ExpectedBytes"] == row["ObservedBytes"] for row in asset_checks), len(asset_checks))
    audit.add("Input", "metadata", len(metadata) == 25 and len({row["Study"] for row in metadata}) == 13 and len({row["Country"] for row in metadata}) == 9, len(metadata))
    audit.add("Input", "pickles-excluded", not any(frozen.rglob("*.pkl")), "serialized inputs remain outside frozen bundle")

    audit.add("Contract", "filters", contract.get("trim_sequences") == 50 and contract.get("sample_with_n_markers") == 20 and contract.get("sample_with_n_markers_perc") == 25 and contract.get("marker_in_n_samples_perc") == 50 and contract.get("sample_with_n_markers_after_filt") == 20 and contract.get("sample_with_n_markers_after_filt_perc") == 25 and contract.get("breadth_thres") == 80 and contract.get("phylophlan_mode") == "fast", contract)
    sensitivity = contract.get("official_baseline_threshold_sensitivity", {})
    audit.add("Contract", "official-threshold-sensitivity", sensitivity == {
        "breadth_thres": 80,
        "marker_in_n_samples_perc": 80,
        "sample_with_n_markers": 20,
        "sample_with_n_markers_after_filt": 20,
        "sample_with_n_markers_after_filt_perc": 80,
        "sample_with_n_markers_perc": 80,
    }, sensitivity)
    audit.add("Contract", "determinism", len(determinism) == 3 and all(as_bool(row["Pass"]) for row in determinism) and any(row["Control"] == "-p 20260751" for row in determinism), determinism)
    config = (frozen / "raw/phylophlan-seeded.cfg").read_text(encoding="utf-8")
    audit.add("Contract", "seeded-config", "params = -p 20260751 -m GTRCAT" in config and "--thread 1" in config and "-p 1989" not in config, config)

    audit.add("Filtering", "input-counts", filtering.get("InputSamples") == 25 and filtering.get("InputReferences") == 6, filtering)
    audit.add("Filtering", "marker-counts", filtering.get("AvailableMarkers", 0) > 0 and 0 < filtering.get("SelectedMarkers", 0) <= filtering.get("AvailableMarkers", 0), filtering)
    audit.add("Filtering", "retained-samples", filtering.get("RetainedSamples") == 25 and summary.get("retained_samples") == 25, filtering)
    audit.add("Filtering", "retained-references", 1 <= filtering.get("RetainedReferences", 0) <= 6 and summary.get("retained_references") == filtering.get("RetainedReferences"), filtering)
    audit.add("Filtering", "alignment", filtering.get("AlignmentSites", 0) > 1000 and filtering.get("TreeTips") == filtering.get("RetainedSamples") + filtering.get("RetainedReferences"), filtering)
    threshold_map = {row["Branch"]: row for row in threshold_rows}
    current_threshold_row = threshold_map.get("Current explicit thresholds", {})
    official_threshold_row = threshold_map.get("Official 2022 thresholds", {})
    audit.add("Filtering", "two-threshold-branches", len(threshold_rows) == 2 and set(threshold_map) == {"Current explicit thresholds", "Official 2022 thresholds"}, threshold_map)
    audit.add("Filtering", "current-threshold-ledger", current_threshold_row.get("sample_with_n_markers_perc") == "25" and current_threshold_row.get("marker_in_n_samples_perc") == "50" and current_threshold_row.get("SelectedMarkers") == str(filtering.get("SelectedMarkers")), current_threshold_row)
    audit.add("Filtering", "official-threshold-ledger", official_threshold_row.get("sample_with_n_markers_perc") == "80" and official_threshold_row.get("marker_in_n_samples_perc") == "80" and int(official_threshold_row.get("SelectedMarkers", 0)) > 0 and int(official_threshold_row.get("AlignmentSites", 0)) > 1000 and int(official_threshold_row.get("TreeTips", 0)) == int(official_threshold_row.get("RetainedSamples", 0)) + int(official_threshold_row.get("RetainedReferences", 0)), official_threshold_row)

    sample_tips = [row for row in tips if row["Type"] == "Metagenome sample"]
    reference_tips = [row for row in tips if row["Type"] == "Reference genome"]
    audit.add("Output", "tip-coordinate", len(sample_tips) == 25 and {row["Tip"] for row in sample_tips} == {row["Sample"] for row in metadata} and len(reference_tips) == filtering["RetainedReferences"], len(tips))
    audit.add("Output", "tip-occupancy", all(0 < float(row["OccupancyPct"]) <= 100 and int(row["ObservedACGTSites"]) <= int(row["AlignmentSites"]) for row in tips), len(tips))
    audit.add("Output", "pair-coordinate", len(pairs) == 300 and len({tuple(sorted((row["Sample1"], row["Sample2"]))) for row in pairs}) == 300, len(pairs))
    audit.add("Output", "pair-bounds", all(int(row["ComparableSites"]) > 0 and int(row["SNPs"]) >= 0 and 0 <= float(row["PDistance"]) <= 1 and abs(float(row["DifferencesPer10kb"]) - float(row["PDistance"]) * 10000) < 0.01 for row in pairs), len(pairs))
    audit.add("Output", "pair-strata", {row["PairStratum"] for row in pairs} == {"Same study", "Same country, different study", "Different country"}, Counter(row["PairStratum"] for row in pairs))
    audit.add("Output", "nearest-neighbors", len(nearest) == 25 and all(row["Sample"] != row["NearestSample"] and finite(row["PDistance"]) for row in nearest), len(nearest))
    audit.add("Output", "polymorphism", len(polymorphism) == 25 and {row["sample"] for row in polymorphism} == {row["Sample"] for row in metadata} and all(0 <= float(row["percentage_of_polymorphic_sites"]) <= 100 for row in polymorphism), len(polymorphism))
    topology_map = {row["Branch"]: row for row in topology}
    current_topology = topology_map.get("Current explicit thresholds", {})
    official_topology = topology_map.get("Official 2022 thresholds", {})
    audit.add("Output", "two-topology-audits", len(topology) == 2 and set(topology_map) == {"Current explicit thresholds", "Official 2022 thresholds"}, topology_map)
    audit.add("Output", "current-topology-audit", int(current_topology.get("CommonTips", 0)) >= 25 and int(current_topology.get("RobinsonFouldsDistance", -1)) >= 0 and 0 <= float(current_topology.get("NormalizedRF", -1)) <= 1 and summary.get("normalized_rf_vs_official_baseline") == float(current_topology.get("NormalizedRF", "nan")), current_topology)
    audit.add("Output", "official-threshold-topology-audit", int(official_topology.get("CommonTips", 0)) >= 25 and int(official_topology.get("RobinsonFouldsDistance", -1)) >= 0 and 0 <= float(official_topology.get("NormalizedRF", -1)) <= 1 and summary.get("official_thresholds_normalized_rf_vs_baseline") == float(official_topology.get("NormalizedRF", "nan")) and summary.get("official_thresholds_selected_markers") == int(official_threshold_row.get("SelectedMarkers", -1)), official_topology)

    command_map = {row["Label"]: row["Command"] for row in commands}
    audit.add("Execution", "three-commands", len(commands) == 3 and all(int(row["ExitStatus"]) == 0 for row in commands), len(commands))
    analysis_command = command_map.get("strainphlan-analysis", "")
    audit.add("Execution", "analysis-command", all(token in analysis_command for token in (
        "--sample_list", "--reference_list", "t__SGB4933_group", "mpa_vJan21_CHOCOPhlAnSGB_202103.pkl",
        "--trim_sequences 50", "--sample_with_n_markers 20", "--sample_with_n_markers_perc 25",
        "--marker_in_n_samples_perc 50", "--sample_with_n_markers_after_filt 20",
        "--sample_with_n_markers_after_filt_perc 25", "--breadth_thres 80",
        "--phylophlan_mode fast", "--phylophlan_configuration",
    )), analysis_command)
    official_threshold_command = command_map.get("strainphlan-official-thresholds", "")
    audit.add("Execution", "official-threshold-command", all(token in official_threshold_command for token in (
        "--sample_with_n_markers 20", "--sample_with_n_markers_perc 80",
        "--marker_in_n_samples_perc 80", "--sample_with_n_markers_after_filt 20",
        "--sample_with_n_markers_after_filt_perc 80", "--breadth_thres 80",
    )), official_threshold_command)
    run_source = (frozen / "scripts/run_article51_strainphlan.py").read_text(encoding="utf-8")
    audit.add("Execution", "checksum-before-pickle", run_source.index("verify_assets(work)") < run_source.index('commands.append(run_timed("strainphlan-analysis"'), "source order")
    audit.add("Execution", "resources", len(resources) == 3 and all(int(row["ExitStatus"]) == 0 and float(row["WallSeconds"]) > 0 and float(row["PeakRAMGiB"]) > 0 for row in resources), resources)
    audit.add("Version", "strainphlan", "4.2.5" in versions.get("StrainPhlAn", ""), versions)
    audit.add("Version", "metaphlan", "4.2.5" in versions.get("MetaPhlAn", ""), versions)
    audit.add("Version", "raxml", "8.2.12" in versions.get("RAxML", ""), versions)
    audit.add("Native", "required-outputs", all((frozen / path).is_file() and (frozen / path).stat().st_size > 0 for path in (
        "raw/strainphlan.info", "raw/strainphlan.polymorphic.tsv",
        "raw/strainphlan-concatenated.aln", "raw/strainphlan-raxml.tree",
        "raw/official-thresholds-strainphlan.info",
        "raw/official-thresholds-strainphlan.polymorphic.tsv",
        "raw/official-thresholds-strainphlan-concatenated.aln",
        "raw/official-thresholds-strainphlan-raxml.tree",
        "raw/official-baseline.info", "raw/official-baseline.tree",
    )), "raw outputs")

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=51,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "samples": summary["retained_samples"],
            "references": summary["retained_references"],
            "selected_markers": summary["selected_markers"],
            "alignment_sites": summary["alignment_sites"],
            "official_thresholds_selected_markers": summary["official_thresholds_selected_markers"],
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
