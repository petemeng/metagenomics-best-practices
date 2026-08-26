#!/usr/bin/env python3
"""Fail-closed validation for Article 47 novelty and phylogenomics evidence."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

from Bio import AlignIO, Phylo

from article42_44_validation_utils import (
    Audit,
    as_bool,
    audit_checksums,
    audit_figures,
    finish,
    read_tsv,
)


FIGURES = (
    "47-bacteria-phylogenomics",
    "47-archaea-phylogenomics",
    "47-novelty-gates",
    "47-marker-recovery-audit",
)


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot = (frozen / "scripts/plot_article47_phylogenomics.R").read_text(encoding="utf-8")
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seeds": all(token in text for token in ("20260747", "20260748")),
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "耗时")),
        "real-catalog": all(token in text for token in ("PRJEB52977", "24 SGBs")),
        "frozen-evidence": "data/small/47-novel-mag-phylogenomics-frozen" in text,
        "methods-results": all(token in text for token in ("Methods template", "Results template")),
        "versions": all(token in text for token in ("GToTree 1.8.17", "IQ-TREE 3.1.3", "GTDB release R232")),
        "tree-contract": all(token in text for token in (
            "-H Bacteria", "-H Archaea", "-m MFP", "-B 1000", "--alrt 1000",
            "--seed 20260747", "--seed 20260748",
        )),
        "result-coordinate": all(token in text for token in ("23/24", "SGB_015", "0 novel species")),
        "truth-blind": "不参与 novelty 判定或参考基因组选择" in text,
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


def close(left: object, right: float, tolerance: float = 1e-6) -> bool:
    try:
        return math.isclose(float(left), right, abs_tol=tolerance)
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
    genomes = read_tsv(frozen / "genome-ledger.tsv")
    references = read_tsv(frozen / "reference-request-ledger.tsv")
    pretree = read_tsv(frozen / "novelty-pretree-audit.tsv")
    novelty = read_tsv(frozen / "novelty-decision-audit.tsv")
    domains = {row["Domain"]: row for row in read_tsv(frozen / "domain-tree-summary.tsv")}
    occupancy = read_tsv(frozen / "alignment-occupancy.tsv")
    tree_tips = read_tsv(frozen / "tree-tip-ledger.tsv")
    gtotree = read_tsv(frozen / "gtotree-genome-audit.tsv")
    inventory = read_tsv(frozen / "reference-download-inventory.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    commands = read_tsv(frozen / "command-log.tsv")
    tools = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}

    audit.add("Identity", "article", summary.get("article") == 47 and contract.get("article") == 47, summary)
    audit.add("Identity", "seed", contract.get("seed") == 20260747, contract)
    audit.add("Identity", "truth-blind", contract.get("truth_used_for_novelty_or_reference_selection") is False and summary.get("truth_used_for_novelty_or_reference_selection") is False, summary)
    audit.add("Identity", "claim-policy", contract.get("novel_species_claim_requires_manual_taxonomic_review") is True and summary.get("novel_species_claims") == 0 and summary.get("candidatus_names_proposed") == 0, summary)

    audit.add("Input", "24-sgbs", len(genomes) == 24 and [row["SGB"] for row in genomes] == [f"SGB_{index:03d}" for index in range(1, 25)], len(genomes))
    audit.add("Input", "domains", Counter(row["Domain"] for row in genomes) == {"Bacteria": 16, "Archaea": 8}, Counter(row["Domain"] for row in genomes))
    audit.add("Input", "genome-checksums", all(re.fullmatch(r"[0-9a-f]{64}", row["RepresentativeSHA256"]) for row in genomes), len(genomes))
    audit.add("Input", "24-versioned-references", len(references) == 24 and all(as_bool(row["ImmutableAccessionVersion"]) and re.fullmatch(r"GC[AF]_\d+\.\d+", row["ReferenceAccession"]) for row in references), len(references))
    audit.add("Input", "reference-inventory-coordinate", len(inventory) == 24 and {row["ReferenceAccession"] for row in inventory} == {row["ReferenceAccession"] for row in references}, len(inventory))
    audit.add("Input", "upstream-checksums", all((frozen / path).is_file() for path in ("upstream/article45-file-checksums.sha256", "upstream/article46-file-checksums.sha256")), "Article 45/46 manifests")

    audit.add("Novelty", "coordinate", len(pretree) == 24 and len(novelty) == 24 and {row["SGB"] for row in novelty} == {row["SGB"] for row in genomes}, len(novelty))
    audit.add("Novelty", "species-assigned", all(as_bool(row["SpeciesAssigned"]) and as_bool(row["SpeciesRadiusAndAFCleared"]) for row in novelty) and summary.get("species_fields_assigned") == 24, len(novelty))
    audit.add("Novelty", "ani-af-bounds", min(float(row["FastANIANI"]) - float(row["SpeciesRadiusPct"]) for row in novelty) > 2.20 and min(float(row["AlignmentFractionPct"]) for row in novelty) > 82.39, "all clear reference radius and AF")
    audit.add("Novelty", "no-candidates-or-claims", all(not as_bool(row["PreTreeNovelSpeciesCandidate"]) and not as_bool(row["NovelSpeciesClaim"]) and not as_bool(row["CandidatusNameProposed"]) for row in novelty) and summary.get("pretree_novelty_candidates") == 0, summary)
    audit.add("Novelty", "quality-eligibility", sum(as_bool(row["QualityEligibleForNoveltyReview"]) for row in novelty) == 17, Counter(row["QualityEligibleForNoveltyReview"] for row in novelty))
    audit.add("Novelty", "tree-inclusion", sum(as_bool(row["PhylogenomicTreeIncluded"]) for row in novelty) == 23 and [row["SGB"] for row in novelty if not as_bool(row["PhylogenomicTreeIncluded"])] == ["SGB_015"], summary)

    bacteria = domains.get("Bacteria", {})
    archaea = domains.get("Archaea", {})
    audit.add("Tree", "two-domains", set(domains) == {"Bacteria", "Archaea"}, domains)
    audit.add("Tree", "bacteria-coordinate", bacteria.get("InputQueryMAGs") == "16" and bacteria.get("QueryMAGs") == "15" and bacteria.get("ExcludedQueryMAGs") == "1" and bacteria.get("ExcludedQueryIDs") == "SGB_015" and bacteria.get("UniqueReferenceTips") == "16" and bacteria.get("TreeTips") == "31" and bacteria.get("AlignmentSites") == "10733", bacteria)
    audit.add("Tree", "archaea-coordinate", archaea.get("InputQueryMAGs") == "8" and archaea.get("QueryMAGs") == "8" and archaea.get("ExcludedQueryMAGs") == "0" and archaea.get("UniqueReferenceTips") == "8" and archaea.get("TreeTips") == "16" and archaea.get("AlignmentSites") == "14804", archaea)
    audit.add("Tree", "models-and-seeds", bacteria.get("BestFitModel") == "LG+I+R4" and bacteria.get("Seed") == "20260747" and archaea.get("BestFitModel") == "LG+F+I+R3" and archaea.get("Seed") == "20260748", domains)
    audit.add("Tree", "supports", bacteria.get("SupportedInternalBranches") == "28" and close(bacteria.get("MedianUFBoot"), 100) and close(bacteria.get("PctUFBootGe95"), 89.285714) and archaea.get("SupportedInternalBranches") == "13" and close(archaea.get("MedianUFBoot"), 100) and close(archaea.get("PctUFBootGe95"), 100), domains)
    audit.add("Tree", "occupancy-coordinate", len(occupancy) == 47 and Counter(row["Domain"] for row in occupancy) == {"Bacteria": 31, "Archaea": 16} and all(0 < float(row["OccupancyPct"]) <= 100 for row in occupancy), len(occupancy))
    audit.add("Tree", "tip-coordinate", len(tree_tips) == 47 and len({(row["Domain"], row["Tip"]) for row in tree_tips}) == 47 and all(as_bool(row["InAlignment"]) for row in tree_tips), len(tree_tips))
    excluded = [row for row in gtotree if not as_bool(row["InFinalTree"])]
    audit.add("Tree", "gtotree-filter", len(gtotree) == 48 and len(excluded) == 1 and excluded[0]["AssemblyID"] == "SGB_015" and excluded[0]["UniqueSCGHits"] == "38" and excluded[0]["SCGsAfterLengthFilter"] == "33", excluded)

    command_map = {row["Label"]: row["Command"] for row in commands}
    audit.add("Execution", "four-commands", len(commands) == 4 and all(int(row["ExitStatus"]) == 0 for row in commands), len(commands))
    audit.add("Execution", "gtotree-domains", "-H Bacteria" in command_map.get("gtotree-bacteria", "") and "-H Archaea" in command_map.get("gtotree-archaea", ""), command_map)
    audit.add("Execution", "iqtree-contract", all(token in command_map.get("iqtree-bacteria", "") for token in ("-m MFP", "-B 1000", "--alrt 1000", "--seed 20260747")) and all(token in command_map.get("iqtree-archaea", "") for token in ("-m MFP", "-B 1000", "--alrt 1000", "--seed 20260748")), command_map)
    audit.add("Execution", "resources", len(resources) == 4 and all(int(row["ExitStatus"]) == 0 and float(row["WallSeconds"]) > 0 and float(row["PeakRAMGiB"]) > 0 for row in resources) and max(float(row["PeakRAMGiB"]) for row in resources) < 0.4, resources)
    audit.add("Version", "tools", "1.8.17" in tools.get("GToTree", "") and "3.1.3" in tools.get("IQ-TREE", "") and "5.1" in tools.get("MUSCLE", "") and "1.5" in tools.get("trimAl", "") and "2.6.3" in tools.get("Prodigal", ""), tools)

    for domain, expected_tips, expected_sites in (("bacteria", 31, 10733), ("archaea", 16, 14804)):
        alignment = AlignIO.read(frozen / f"raw/{domain}/alignment.faa", "fasta")
        tree = Phylo.read(frozen / f"raw/{domain}/iqtree.treefile", "newick")
        audit.add("Native", f"{domain}-alignment-tree", len(alignment) == expected_tips and alignment.get_alignment_length() == expected_sites and {record.id for record in alignment} == {tip.name for tip in tree.get_terminals()}, {"tips": len(alignment), "sites": alignment.get_alignment_length()})
        audit.add("Native", f"{domain}-reports", all((frozen / f"raw/{domain}/{name}").is_file() and (frozen / f"raw/{domain}/{name}").stat().st_size > 0 for name in ("gtotree-genomes.tsv", "gtotree-scg-hits.tsv", "gtotree-runlog.txt", "iqtree-report.txt", "iqtree.log")), domain)

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=47,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "input_sgbs": 24,
            "tree_included_sgbs": 23,
            "tree_excluded_sgbs": 1,
            "novel_species_claims": 0,
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
