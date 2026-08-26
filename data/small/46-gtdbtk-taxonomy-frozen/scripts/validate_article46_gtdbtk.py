#!/usr/bin/env python3
"""Fail-closed validation for Article 46 GTDB-Tk R232 evidence."""

from __future__ import annotations

import argparse
import csv
import json
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
    "46-rank-resolution",
    "46-phylum-composition",
    "46-classification-route",
    "46-ani-af-audit",
)


def raw_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot_script = (frozen / "scripts/plot_article46_gtdbtk.R").read_text(encoding="utf-8")
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260746" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "耗时")),
        "real-catalog": all(token in text for token in ("PRJEB52977", "24", "SGB")),
        "frozen-evidence": "data/small/46-gtdbtk-taxonomy-frozen" in text,
        "methods-results": all(token in text for token in ("Methods template", "Results template")),
        "versions": all(token in text for token in ("GTDB-Tk 2.7.2", "R232", "skani 0.3.2")),
        "database-contract": all(token in text for token in ("25a59e0352b1fd150c589f56559767d4", "60,806,405,195", "100,482,230,310")),
        "classification-contract": all(token in text for token in ("--min_perc_aa 10", "--min_af 0.5", "--pplacer_cpus 1")),
        "result-coordinate": all(token in text for token in ("16 个 Bacteria", "8 个 Archaea", "24/24", "ani_screen")),
        "truth-blind": "不参与分类、参考基因组选择或物种判定" in text,
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
    taxonomy = read_tsv(frozen / "taxonomy-summary.tsv")
    ranks = read_tsv(frozen / "rank-resolution.tsv")
    phyla = read_tsv(frozen / "phylum-summary.tsv")
    methods = read_tsv(frozen / "classification-method-summary.tsv")
    ani = read_tsv(frozen / "fastani-reference-audit.tsv")
    markers = read_tsv(frozen / "marker-file-inventory.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    commands = read_tsv(frozen / "command-log.tsv")
    tools = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    database = read_tsv(frozen / "database-audit.tsv")

    expected_summary = {
        "article": 46,
        "classification_methods": {"ani_screen": 24},
        "classified_sgbs": 24,
        "domain_counts": {"Archaea": 8, "Bacteria": 16},
        "fastani_reference_hits": 24,
        "gtdb_release": "R232",
        "gtdbtk_version": "2.7.2",
        "input_sgbs": 24,
        "single_copy_or_msa_files": 0,
        "species_assigned": 24,
        "species_unresolved_candidates": 0,
        "truth_used_for_taxonomy": False,
        "warnings_nonempty": 0,
    }
    audit.add("Identity", "summary", summary == expected_summary, summary)
    audit.add("Identity", "contract", contract.get("article") == 46 and contract.get("seed") == 20260746 and contract.get("random_process") is False and contract.get("truth_used_for_taxonomy") is False, contract)
    audit.add("Contract", "release-and-tool", contract.get("gtdbtk_version") == "2.7.2" and contract.get("gtdb_release") == "R232" and contract.get("full_tree") is False, contract)
    audit.add("Contract", "thresholds", contract.get("minimum_percent_aa") == 10 and contract.get("minimum_alignment_fraction") == 0.5 and contract.get("pplacer_cpus") == 1, contract)
    audit.add("Input", "catalog-coordinate", len(taxonomy) == 24 and [row["SGB"] for row in taxonomy] == [f"SGB_{i:03d}" for i in range(1, 25)], len(taxonomy))
    audit.add("Taxonomy", "domain-counts", Counter(row["Domain"] for row in taxonomy) == {"d__Bacteria": 16, "d__Archaea": 8}, Counter(row["Domain"] for row in taxonomy))
    audit.add("Taxonomy", "seven-ranks-resolved", all(all(len(row[rank]) > 3 for rank in ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) for row in taxonomy), len(taxonomy))
    audit.add("Taxonomy", "species-resolved", all(as_bool(row["SpeciesAssigned"]) and not as_bool(row["GTDBSpeciesUnresolvedCandidate"]) for row in taxonomy), len(taxonomy))
    audit.add("Taxonomy", "ani-only", all(row["ClassificationMethod"] == "ani_screen" for row in taxonomy) and methods == [{"ClassificationMethod": "ani_screen", "SGBs": "24"}], methods)
    audit.add("Taxonomy", "no-warnings", all(not row["Warnings"] for row in taxonomy), len(taxonomy))
    audit.add("Ranks", "complete", len(ranks) == 14 and all(int(row["ResolvedSGBs"]) == int(row["SGBs"]) and int(row["UnresolvedSGBs"]) == 0 for row in ranks), len(ranks))
    audit.add("Ranks", "phylum-counts", len(phyla) == 16 and sum(int(row["SGBs"]) for row in phyla) == 24, len(phyla))
    audit.add("ANI", "coordinate", len(ani) == 24 and {row["SGB"] for row in ani} == {row["SGB"] for row in taxonomy}, len(ani))
    audit.add("ANI", "versioned-references", all(re.fullmatch(r"GC[AF]_\d+\.\d+", row["FastANIReference"]) for row in ani), len(ani))
    audit.add("ANI", "radius-and-af", all(float(row["ANIMarginToReferenceRadiusPctPoints"]) > 0 and float(row["FastANIAFPct"]) >= 50 for row in ani), min(float(row["ANIMarginToReferenceRadiusPctPoints"]) for row in ani))
    audit.add("ANI", "observed-bounds", min(float(row["FastANIANI"]) for row in ani) == 97.21 and max(float(row["FastANIANI"]) for row in ani) == 100.0 and min(float(row["FastANIAFPct"]) for row in ani) > 82.39, {"ani_min": min(float(row["FastANIANI"]) for row in ani), "af_min": min(float(row["FastANIAFPct"]) for row in ani)})
    audit.add("Markers", "not-invoked-after-ani-resolution", len(markers) == 0, len(markers))
    audit.add("Database", "one-row", len(database) == 1, database)
    if database:
        row = database[0]
        audit.add("Database", "release", row["ToolVersion"] == "2.7.2" and row["Release"] == "R232" and row["ArchiveChecksum"] == "25a59e0352b1fd150c589f56559767d4", row)
        audit.add("Database", "sizes", int(row["ArchiveBytes"]) == 60_806_405_195 and int(row["InstalledFiles"]) == 263 and int(row["InstalledBytes"]) == 100_482_230_310, row)
        audit.add("Database", "verified", row["ManifestStatus"] == "VERIFIED_LOCAL_CHECK_INSTALL_PASS" and row["LocalStatus"] == "ARCHIVE_MD5_BYTES_EXTRACTION_AND_GTDBTK_CHECK_INSTALL_PASS", row)
    audit.add("Version", "tools", tools == {
        "GTDB-Tk": "gtdbtk: version 2.7.2 Copyright 2017 Pierre-Alain Chaumeil, Aaron Mussig and Donovan Parks",
        "Prodigal": "Prodigal V2.6.3: February, 2016",
        "pplacer": "v1.1.alpha19-0-g807f6f3",
        "skani": "skani 0.3.2",
    }, tools)
    audit.add("Execution", "commands", len(commands) == 2 and all(int(row["ExitStatus"]) == 0 for row in commands) and "--min_perc_aa 10" in commands[1]["Command"] and "--min_af 0.5" in commands[1]["Command"] and "--pplacer_cpus 1" in commands[1]["Command"], commands)
    audit.add("Execution", "resources", len(resources) == 2 and all(int(row["ExitStatus"]) == 0 and float(row["WallSeconds"]) > 0 and float(row["PeakRAMGiB"]) > 0 for row in resources) and 14 < float(resources[1]["PeakRAMGiB"]) < 16, resources)

    bac = raw_tsv(frozen / "raw/gtdbtk/classify/article46.bac120.summary.tsv")
    arc = raw_tsv(frozen / "raw/gtdbtk/classify/article46.ar53.summary.tsv")
    raw_ani_bac = raw_tsv(frozen / "raw/gtdbtk/classify/ani_screen/article46.bac120.ani_summary.tsv")
    raw_ani_arc = raw_tsv(frozen / "raw/gtdbtk/classify/ani_screen/article46.ar53.ani_summary.tsv")
    audit.add("Native", "summary-counts", len(bac) == 16 and len(arc) == 8 and {row["user_genome"] for row in bac + arc} == {f"SGB_{i:03d}" for i in range(1, 25)}, {"bac": len(bac), "arc": len(arc)})
    audit.add("Native", "ani-counts", len(raw_ani_bac) == 16 and len(raw_ani_arc) == 8, {"bac": len(raw_ani_bac), "arc": len(raw_ani_arc)})
    audit.add("Native", "raw-method", all(row["classification_method"] == "ani_screen" for row in bac + arc), len(bac + arc))

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=46,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={"classified_sgbs": 24, "bacteria": 16, "archaea": 8, "species_assigned": 24},
    )


if __name__ == "__main__":
    raise SystemExit(main())
