#!/usr/bin/env python3
"""Fail-closed validation for Article 48 CoverM MAG abundance evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

from article42_44_validation_utils import Audit, as_bool, audit_checksums, audit_figures, finish, read_tsv


FIGURES = (
    "48-mag-abundance-heatmap",
    "48-breadth-depth-audit",
    "48-identity-sensitivity",
    "48-catalog-capture",
)


def raw_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot_script = (frozen / "scripts/plot_article48_coverm.R").read_text(encoding="utf-8")
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260748" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "耗时")),
        "real-study": all(token in text for token in ("PRJEB52977", "ERR9765746", "ERR9765747")),
        "frozen-inputs": all(token in text for token in ("45-drep-dereplication-frozen", "48-mag-abundance-coverm-frozen")),
        "methods-results": all(token in text for token in ("Methods template", "Results template")),
        "versions": all(token in text for token in ("CoverM 0.8.0", "Strobealign 0.17.0", "SAMtools 1.23.1")),
        "mapping-contract": all(token in text for token in ("--min-read-percent-identity 95", "--min-read-aligned-percent 75", "--proper-pairs-only", "--min-covered-fraction 0")),
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot_script) is None,
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-meta-prose": not any(token in text for token in (
            "本篇可独立", "本文可独立", "全系列约定", "接口只学一次",
            "作者代码通常长这样", "（即本文）", "无头服务器",
        )),
    }
    sections = (
        "对应论文里的哪张图", "理论：", "准备工作", "可复制代码", "审计与升级",
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
    summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    contract = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    genomes = read_tsv(frozen / "genome-ledger.tsv")
    samples = read_tsv(frozen / "samples.tsv")
    coverage = read_tsv(frozen / "coverm-long.tsv.gz")
    capture = read_tsv(frozen / "sample-capture-summary.tsv")
    thresholds = read_tsv(frozen / "detection-threshold-summary.tsv")
    sensitivity = read_tsv(frozen / "stringency-sensitivity.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    commands = read_tsv(frozen / "command-log.tsv")
    tools = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}

    audit.add("Identity", "article", summary.get("article") == 48 and contract.get("article") == 48, summary)
    audit.add("Identity", "seed", contract.get("seed") == 20260748, contract.get("seed"))
    audit.add("Identity", "truth-blind", summary.get("truth_used_for_mapping_or_detection") is False and contract.get("truth_used_for_mapping_or_detection") is False, summary)
    audit.add("Contract", "deterministic", contract.get("random_process") is False and contract.get("genome_order") == [f"SGB_{i:03d}" for i in range(1, 25)], contract)
    audit.add("Contract", "filters", contract.get("main_read_identity_pct") == 95 and contract.get("strict_read_identity_pct") == 97 and contract.get("minimum_read_aligned_pct") == 75 and contract.get("proper_pairs_only") is True and contract.get("minimum_covered_fraction_pct") == 0, contract)
    audit.add("Input", "catalog-coordinate", len(genomes) == 24 and [row["SGB"] for row in genomes] == [f"SGB_{i:03d}" for i in range(1, 25)], len(genomes))
    audit.add("Input", "catalog-hashes", all(re.fullmatch(r"[0-9a-f]{64}", row["RepresentativeSHA256"]) for row in genomes) and len({row["RepresentativeSHA256"] for row in genomes}) == 24, len(genomes))
    audit.add("Input", "samples", len(samples) == 2 and [row["Sample"] for row in samples] == ["MOCK1", "MOCK2"] and [row["RunAccession"] for row in samples] == ["ERR9765746", "ERR9765747"], samples)
    audit.add("Input", "read-hashes", all(re.fullmatch(r"[0-9a-f]{64}", row[key]) for row in samples for key in ("R1SHA256", "R2SHA256")), len(samples))

    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in coverage:
        groups[(row["Branch"], row["Sample"])].append(row)
    expected_groups = {
        ("Primary 95% identity", "MOCK1"), ("Primary 95% identity", "MOCK2"),
        ("Strict 97% identity", "MOCK1"), ("Strict 97% identity", "MOCK2"),
    }
    audit.add("Coverage", "96-rows", len(coverage) == 96 and set(groups) == expected_groups, len(coverage))
    genome_ids = {row["SGB"] for row in genomes}
    for key in sorted(expected_groups):
        rows = groups[key]
        audit.add("Coverage", f"{key}-coordinate", len(rows) == 24 and {row["SGB"] for row in rows} == genome_ids, len(rows))
        audit.add("Coverage", f"{key}-bounds", all(float(row["MeanDepth"]) >= 0 and 0 <= float(row["CoveredFractionPct"]) <= 100 and 0 <= float(row["RelativeAbundancePct"]) <= 100 and 0 <= float(row["MeanReadIdentityPct"]) <= 100 for row in rows), key)
        audit.add("Coverage", f"{key}-detection-formula", all(as_bool(row["DetectedBreadth50Depth1"]) == (float(row["CoveredFractionPct"]) >= 50 and float(row["MeanDepth"]) >= 1) for row in rows), key)
        audit.add("Coverage", f"{key}-catalog-normalization", abs(sum(float(row["CatalogNormalizedAbundancePct"]) for row in rows) - 100) < 1e-6, sum(float(row["CatalogNormalizedAbundancePct"]) for row in rows))

    capture_map = {(row["Branch"], row["Sample"]): row for row in capture}
    audit.add("Capture", "four-groups", len(capture) == 4 and set(capture_map) == expected_groups, len(capture))
    for key, row in capture_map.items():
        cover_rows = groups[key]
        audit.add("Capture", f"{key}-sum", abs(float(row["CatalogRelativeAbundancePct"]) - sum(float(item["RelativeAbundancePct"]) for item in cover_rows)) < 1e-6 and abs(float(row["CompositionSumPct"]) - 100) < 2e-5, row)
        audit.add("Capture", f"{key}-detected", int(row["DetectedSGBsBreadth50Depth1"]) == sum(as_bool(item["DetectedBreadth50Depth1"]) for item in cover_rows), row["DetectedSGBsBreadth50Depth1"])
    audit.add("Capture", "primary-exact", abs(float(capture_map[("Primary 95% identity", "MOCK1")]["CatalogRelativeAbundancePct"]) - 70.1154315) < 1e-6 and abs(float(capture_map[("Primary 95% identity", "MOCK2")]["CatalogRelativeAbundancePct"]) - 69.86899505) < 1e-6, summary)
    audit.add("Detection", "threshold-table", len(thresholds) == 12 and all(int(row["CatalogSGBs"]) == 24 for row in thresholds), len(thresholds))
    audit.add("Detection", "primary-50pct-all", all(int(row["DetectedSGBs"]) == 24 for row in thresholds if row["Branch"] == "Primary 95% identity" and int(row["BreadthCutoffPct"]) == 50), summary)

    sensitivity_map = {(row["Sample"], row["SGB"]): row for row in sensitivity}
    audit.add("Sensitivity", "coordinate", len(sensitivity) == 48 and set(sensitivity_map) == {(sample, sgb) for sample in ("MOCK1", "MOCK2") for sgb in genome_ids}, len(sensitivity))
    audit.add("Sensitivity", "strict-count-monotonic", all(int(row["StrictReadCount"]) <= int(row["PrimaryReadCount"]) for row in sensitivity), len(sensitivity))
    max_delta = max(abs(float(row["DeltaRelativeAbundancePctPoints"])) for row in sensitivity)
    audit.add("Sensitivity", "max-delta", abs(max_delta - float(summary["maximum_absolute_stringency_delta_pct_points"])) < 1e-6 and max_delta < 0.30, max_delta)
    audit.add("Result", "summary-coordinate", summary.get("catalog_sgbs") == 24 and summary.get("samples") == 2 and summary.get("detected_sgbs_main_breadth50_depth1") == {"MOCK1": 24, "MOCK2": 24} and summary.get("high_breadth_sgbs_main") == {"MOCK1": 22, "MOCK2": 21}, summary)
    audit.add("Result", "sample-correlation", abs(float(summary["mock1_mock2_relative_abundance_spearman"]) - 0.966957) < 1e-6, summary["mock1_mock2_relative_abundance_spearman"])

    for branch, identity in (("identity95", "95"), ("identity97", "97")):
        raw = raw_tsv(frozen / f"raw/coverm-{branch}.tsv")
        audit.add("Native", f"{branch}-rows", len(raw) == 50 and sum(row["Genome"] == "unmapped" for row in raw) == 2, len(raw))
        command = next(row["Command"] for row in commands if row["Label"] == f"coverm-{branch}")
        audit.add("Native", f"{branch}-command", f"--min-read-percent-identity {identity}" in command and "--min-read-aligned-percent 75" in command and "--proper-pairs-only" in command and "--exclude-supplementary" in command and "--min-covered-fraction 0" in command, command)
    audit.add("Version", "tools", tools == {"CoverM": "coverm 0.8.0", "Strobealign": "0.17.0", "SAMtools": "1.23.1"}, tools)
    audit.add("Execution", "commands", len(commands) == 2 and all(int(row["ExitStatus"]) == 0 for row in commands), len(commands))
    audit.add("Execution", "resources", len(resources) == 2 and all(int(row["ExitStatus"]) == 0 and float(row["WallSeconds"]) > 0 and float(row["PeakRAMGiB"]) > 0 for row in resources), resources)
    audit.add("Execution", "bam-excluded", not (frozen / "bam").exists(), "cached BAMs remain outside frozen evidence")

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=48, audit=audit, output=args.output_dir.resolve(),
        payload={"catalog_sgbs": 24, "samples": 2, "main_capture_pct": summary["main_catalog_capture_pct"]},
    )


if __name__ == "__main__":
    raise SystemExit(main())
