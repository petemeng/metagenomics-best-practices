#!/usr/bin/env python3
"""Fail-closed validation for the Article 52 PanPhlAn tutorial and frozen evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
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
    sha256,
)


FIGURES = (
    "52-pangenome-prevalence",
    "52-gene-content-pcoa",
    "52-gene-content-dendrogram",
    "52-accessory-gene-heatmap",
    "52-plateau-sensitivity",
)
CATEGORIES = {"Core >=95%", "Accessory 5-<95%", "Rare >0-<5%", "Undetected"}


def compressed_shape(path: Path) -> tuple[int, int, list[str]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = sum(1 for row in reader if row)
    return rows, len(header) - 1, header[1:]


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot = (frozen / "scripts/plot_article52_panphlan.R").read_text(encoding="utf-8")
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260752" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "scale_color_pub <-", "scale_fill_pub", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "耗时")),
        "official-real-data": all(token in text for token in ("25 个真实", "13 studies", "9 countries", "15 个参考基因组", "11,069")),
        "frozen-input": "data/small/52-panphlan-pangenome-frozen" in text,
        "methods-results": all(token in text for token in ("Methods template", "Results template")),
        "software-versions": all(token in text for token in ("PanPhlAn profiling source 3.1", "Python 3.12.13", "NumPy 2.3.5", "Bowtie2 2.5.5", "SAMtools 1.23.1", "R 4.4.1")),
        "source-provenance": all(token in text for token in ("4294c3f5c92be9b9ef7d61b69df43e7f27c51601", "328be1e6074015c5183668913ed5e0cf0f879fe5", "411d8736790ef591e8131d58e096d68a15c477bc21eb39b41c337af950d4dfc6", "5f9bfecbe3e5459e06f48677076b49a4a37263935ac8cff8422cf1ef6bc9bc50")),
        "primary-thresholds": all(token in text for token in ("--min_coverage 2", "--left_max 1.25", "--right_min 0.75", "--th_non_present 0.25", "--th_present 0.5", "--th_multicopy 1.5")),
        "sensitivity-thresholds": all(token in text for token in ("--min_coverage 1", "--left_max 1.70", "--right_min 0.30", "--i_covmat")),
        "result-counts": all(token in text for token in ("1,591", "4,542", "1,613", "3,323", "1,571", "3,191")),
        "interpretation-boundaries": all(token in text for token in ("potential multi-strain", "不是 species phylogeny", "CARD", "VFDB", "antiSMASH", "伪重复")),
        "citations": all(token in text for token in ("@scholz2016panphlan", "@segatalab2022panphlan3tutorial", "@beghini2021biobakery", "@suzek2007uniref")),
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot) is None,
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-meta-prose": not any(token in text for token in ("本篇可独立", "本文可独立", "全系列约定", "接口只学一次", "作者代码通常长这样", "（即本文）", "无头服务器")),
    }
    sections = (
        "这一步对应论文里的哪张图", "理论：", "准备工作", "可复制代码", "审计与升级",
        "出版级美化", "常见坑", "这段 Methods 怎么写", "换成你自己的数据怎么做", "参考",
    )
    for section in sections:
        checks[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


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

    summary = json.loads((frozen / "summary.json").read_text(encoding="utf-8"))
    expected_summary = {
        "official_samples": 25,
        "studies": 13,
        "countries": 9,
        "reference_genomes": 15,
        "pangenome_gene_families": 11069,
        "coverage_matrix_gene_families": 10703,
        "primary_retained_samples": 22,
        "sensitivity_retained_samples": 25,
        "primary_multistrain_warnings": 13,
        "primary_operational_core": 1591,
        "primary_accessory": 4542,
        "primary_rare": 1613,
        "primary_undetected": 3323,
        "sensitivity_operational_core": 1571,
        "reference_strict_core": 1335,
        "annotation_unique_gene_families": 11061,
        "annotation_missing_gene_families": 8,
        "annotation_format_anomalies": 5,
    }
    for key, expected in expected_summary.items():
        audit.add("Summary", key, summary.get(key) == expected, {"expected": expected, "observed": summary.get(key)})
    audit.add("Summary", "excluded-samples", set(summary["primary_excluded_samples"]) == {"RA023", "SRS011302", "W1.27.ST"}, summary["primary_excluded_samples"])
    audit.add("Summary", "shared-calls-identical", summary["common_sample_calls_byte_equivalent"] is True, summary["common_sample_calls_byte_equivalent"])
    audit.add("Summary", "nearest-pair", set(summary["nearest_pair"]) == {"SZAXPI003417-4", "T2D-025"} and math.isclose(float(summary["nearest_pair_jaccard"]), 0.24074074, abs_tol=1e-8), summary["nearest_pair"])
    audit.add("Summary", "ordination", math.isclose(float(summary["pcoa_axis1_pct"]), 21.6106, abs_tol=1e-4) and math.isclose(float(summary["pcoa_axis2_pct"]), 15.1056, abs_tol=1e-4), [summary["pcoa_axis1_pct"], summary["pcoa_axis2_pct"]])

    assets = read_tsv(frozen / "asset-check-audit.tsv")
    audit.add("Input", "asset-count", len(assets) == 38, len(assets))
    audit.add("Input", "all-assets-checksum-pass", all(as_bool(row["ChecksumPass"]) for row in assets), Counter(row["AssetType"] for row in assets))
    metadata = read_tsv(frozen / "sample-metadata.tsv")
    audit.add("Input", "metadata-coordinate", len(metadata) == 25 and len({row["Sample"] for row in metadata}) == 25 and len({row["Study"] for row in metadata}) == 13 and len({row["Country"] for row in metadata}) == 9, len(metadata))
    audit.add("Input", "lossless-decode", all("no value transformation" in row["CompatibilityTransformation"] for row in metadata), "25 profiles")
    pangenome_manifest = read_tsv(frozen / "pangenome-file-manifest.tsv")
    audit.add("Input", "pangenome-members", len(pangenome_manifest) == 9 and all(int(row["Bytes"]) > 0 and len(row["SHA256"]) == 64 for row in pangenome_manifest), len(pangenome_manifest))
    audit.add("Input", "pinned-source", sha256(frozen / "pinned-source/panphlan_profiling.py") == "901b4dc3710a26145dca064216d62769772b743db55683a956b66253565e7480", "PanPhlAn source")

    branch_rows = {row["Branch"]: row for row in read_tsv(frozen / "pangenome-summary.tsv")}
    primary = branch_rows.get("Primary", {})
    sensitivity = branch_rows.get("Sensitivity", {})
    references = branch_rows.get("Reference genomes", {})
    audit.add("Output", "three-summary-branches", set(branch_rows) == {"Primary", "Sensitivity", "Reference genomes"}, list(branch_rows))
    audit.add("Output", "primary-partition", [int(primary.get(key, -1)) for key in ("OperationalCore95Pct", "Accessory5To95Pct", "RareBelow5Pct", "Undetected")] == [1591, 4542, 1613, 3323], primary)
    audit.add("Output", "sensitivity-partition", [int(sensitivity.get(key, -1)) for key in ("OperationalCore95Pct", "Accessory5To95Pct", "RareBelow5Pct", "Undetected")] == [1571, 4914, 1720, 2864], sensitivity)
    audit.add("Output", "reference-partition", int(references.get("StrictCore100Pct", -1)) == 1335 and int(references.get("Undetected", -1)) == 0 and float(references.get("MedianGeneFamiliesPerSample", -1)) == 3042.0, references)

    prevalence = read_tsv(frozen / "gene-family-prevalence.tsv.gz")
    audit.add("Output", "prevalence-coordinate", len(prevalence) == 11069 and len({row["GeneFamily"] for row in prevalence}) == 11069, len(prevalence))
    audit.add("Output", "prevalence-bounds", all(0 <= float(row["PrimaryPrevalence"]) <= 1 and 0 <= float(row["SensitivityPrevalence"]) <= 1 and 0 <= float(row["ReferencePrevalence"]) <= 1 for row in prevalence), "bounded")
    audit.add("Output", "prevalence-categories", set(row["PrimaryCategory"] for row in prevalence) == CATEGORIES and Counter(row["PrimaryCategory"] for row in prevalence) == Counter({"Core >=95%": 1591, "Accessory 5-<95%": 4542, "Rare >0-<5%": 1613, "Undetected": 3323}), Counter(row["PrimaryCategory"] for row in prevalence))

    samples = read_tsv(frozen / "sample-filter-audit.tsv")
    primary_retained = [row for row in samples if as_bool(row["PrimaryRetained"])]
    primary_excluded = [row for row in samples if not as_bool(row["PrimaryRetained"])]
    audit.add("Output", "sample-filter", len(samples) == 25 and len(primary_retained) == 22 and all(as_bool(row["SensitivityRetained"]) for row in samples), [len(primary_retained), len(primary_excluded)])
    audit.add("Output", "sample-exclusions", {row["Sample"] for row in primary_excluded} == {"RA023", "SRS011302", "W1.27.ST"} and all(float(row["LeftCoverage"]) > 1.25 and row["PrimaryDecision"] == "Failed primary left-side plateau gate" for row in primary_excluded), primary_excluded)
    audit.add("Output", "multistrain-warnings", sum(as_bool(row["MultiStrainWarning"]) for row in primary_retained) == 13, 13)

    pairs = read_tsv(frozen / "pairwise-jaccard.tsv")
    pair_counts = Counter(row["PairStratum"] for row in pairs)
    unique_pairs = {tuple(sorted((row["Sample1"], row["Sample2"]))) for row in pairs}
    audit.add("Output", "pair-coordinate", len(pairs) == 231 and len(unique_pairs) == 231, len(pairs))
    audit.add("Output", "pair-strata", pair_counts == Counter({"Different country": 170, "Same country, different study": 45, "Same study": 16}), pair_counts)
    audit.add("Output", "jaccard-bounds", all(0 <= float(row["JaccardDistance"]) <= 1 for row in pairs), "bounded")
    pair_summary = {row["PairStratum"]: row for row in read_tsv(frozen / "pairwise-jaccard-summary.tsv")}
    medians = {key: float(row["MedianJaccard"]) for key, row in pair_summary.items()}
    audit.add("Output", "jaccard-medians", math.isclose(medians["Same study"], 0.36702577, abs_tol=1e-8) and math.isclose(medians["Same country, different study"], 0.36640798, abs_tol=1e-8) and math.isclose(medians["Different country"], 0.43573732, abs_tol=1e-8), medians)

    pcoa = read_tsv(frozen / "pcoa-jaccard.tsv")
    audit.add("Output", "pcoa-coordinate", len(pcoa) == 22 and {row["Sample"] for row in pcoa} == {row["Sample"] for row in primary_retained}, len(pcoa))
    heatmap = read_tsv(frozen / "accessory-feature-heatmap.tsv")
    audit.add("Output", "heatmap-coordinate", len(heatmap) == 440 and len({row["GeneFamily"] for row in heatmap}) == 20 and len({row["Sample"] for row in heatmap}) == 22 and {row["Present"] for row in heatmap} == {"0", "1"}, len(heatmap))
    anomalies = read_tsv(frozen / "annotation-format-audit.tsv")
    audit.add("Output", "annotation-anomalies-audited", len(anomalies) == 5 and {int(row["Line"]) for row in anomalies} == {8286, 9449, 12381, 16671, 16672}, anomalies)
    transitions = read_tsv(frozen / "category-transition.tsv")
    transition_map = {(row["PrimaryCategory"], row["SensitivityCategory"]): int(row["GeneFamilies"]) for row in transitions}
    audit.add("Output", "category-transition", len(transitions) == 16 and sum(transition_map.values()) == 11069 and transition_map[("Core >=95%", "Accessory 5-<95%")] == 20 and transition_map[("Rare >0-<5%", "Accessory 5-<95%")] == 335 and transition_map[("Undetected", "Rare >0-<5%")] == 442, transitions)

    native = read_tsv(frozen / "native-output-manifest.tsv")
    audit.add("Native", "compressed-output-count", len(native) == 4, [row["File"] for row in native])
    audit.add("Native", "compressed-output-checksums", all((frozen / row["File"]).is_file() and (frozen / row["File"]).stat().st_size == int(row["Bytes"]) and sha256(frozen / row["File"]) == row["SHA256"] for row in native), "four compressed tables")
    primary_shape = compressed_shape(frozen / "primary-presence-absence.tsv.gz")
    sensitive_shape = compressed_shape(frozen / "sensitive-presence-absence.tsv.gz")
    coverage_shape = compressed_shape(frozen / "coverage-matrix.tsv.gz")
    audit.add("Native", "matrix-shapes", primary_shape[:2] == (11069, 37) and sensitive_shape[:2] == (11069, 40) and coverage_shape[:2] == (10703, 25), [primary_shape[:2], sensitive_shape[:2], coverage_shape[:2]])

    determinism = read_tsv(frozen / "determinism-audit.tsv")
    audit.add("Execution", "byte-identical-replay", len(determinism) == 3 and all(as_bool(row["ByteIdentical"]) and row["Seed"] == "20260752" and not as_bool(row["RandomCoveragePlotRequested"]) for row in determinism), determinism)
    commands = read_tsv(frozen / "command-log.tsv")
    command_map = {row["Label"]: row["Command"] for row in commands}
    audit.add("Execution", "three-commands", len(commands) == 3 and all(int(row["ExitStatus"]) == 0 for row in commands), len(commands))
    primary_command = command_map.get("panphlan-primary", "")
    sensitive_command = command_map.get("panphlan-sensitivity", "")
    audit.add("Execution", "primary-command", all(token in primary_command for token in ("--min_coverage 2", "--left_max 1.25", "--right_min 0.75", "--th_non_present 0.25", "--th_present 0.5", "--th_multicopy 1.5", "--add_ref")) and "--o_covplot_normed" not in primary_command, primary_command)
    audit.add("Execution", "sensitivity-command", all(token in sensitive_command for token in ("--i_covmat", "--min_coverage 1", "--left_max 1.70", "--right_min 0.30", "--add_ref")) and "--o_covplot_normed" not in sensitive_command, sensitive_command)
    resources = read_tsv(frozen / "resource-summary.tsv")
    audit.add("Execution", "resources", len(resources) == 3 and all(int(row["ExitStatus"]) == 0 and float(row["WallSeconds"]) > 0 and float(row["PeakRAMGiB"]) > 0 for row in resources), resources)
    versions = {row["Software"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    audit.add("Version", "analysis-stack", versions.get("PanPhlAn profiling source") == "3.1" and versions.get("Python") == "3.12.13" and versions.get("NumPy") == "2.3.5" and "2.5.5" in versions.get("Bowtie2", "") and "1.23.1" in versions.get("SAMtools", ""), versions)
    plot_versions = {row["Software"]: row["Version"] for row in read_tsv(frozen / "plot-software-versions.tsv")}
    audit.add("Version", "plot-stack", plot_versions == {"R": "4.4.1", "ggplot2": "3.5.2", "dplyr": "1.1.4", "readr": "2.1.5", "ggrepel": "0.9.5", "patchwork": "1.3.2", "ggdendro": "0.2.0"}, plot_versions)

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=52,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "samples": summary["primary_retained_samples"],
            "sensitivity_samples": summary["sensitivity_retained_samples"],
            "gene_families": summary["pangenome_gene_families"],
            "operational_core": summary["primary_operational_core"],
            "accessory": summary["primary_accessory"],
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
