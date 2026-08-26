#!/usr/bin/env python3
"""Fail-closed validation for Article 50 inStrain microdiversity evidence."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
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
    "50-nucleotide-diversity",
    "50-snv-density",
    "50-profile-detection-audit",
    "50-popani-overlap",
)


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot = (frozen / "scripts/plot_article50_instrain.R").read_text(encoding="utf-8")
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260750" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "耗时")),
        "real-study": all(token in text for token in ("PRJEB52977", "ERR9765746", "ERR9765747")),
        "frozen-input": "50-instrain-microdiversity-frozen" in text,
        "methods-results": all(token in text for token in ("Methods template", "Results template")),
        "versions": all(token in text for token in ("inStrain 1.10.0", "SAMtools 1.23.1", "Python 3.8.20")),
        "profile-contract": all(token in text for token in (
            "--min_read_ani 0.95", "--min_mapq 2", "--pairing_filter paired_only",
            "--min_cov 5", "--min_freq 0.05", "--fdr 1e-06",
            "--rarefied_coverage 1000000000",
        )),
        "comparison-rule": all(token in text for token in ("99.999%", "50%", "popANI")),
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
        "对应论文里的哪张图", "理论：", "准备工作", "可复制代码", "审计与升级",
        "出版级美化", "常见坑", "这段 Methods 怎么写", "换成你自己的数据怎么做", "参考",
    )
    for section in sections:
        checks[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def finite(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except ValueError:
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
    samples = read_tsv(frozen / "sample-ledger.tsv")
    profiles = read_tsv(frozen / "profile-genome-long.tsv.gz")
    comparisons = read_tsv(frozen / "pairwise-genome-comparison.tsv")
    filtering = read_tsv(frozen / "read-filter-audit.tsv")
    commands = read_tsv(frozen / "command-log.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    versions = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}

    expected_sgbs = {f"SGB_{index:03d}" for index in range(1, 25)}
    audit.add("Identity", "article", summary.get("article") == 50 and contract.get("article") == 50, summary)
    audit.add("Identity", "seed", contract.get("seed") == 20260750, contract.get("seed"))
    audit.add("Identity", "truth-blind", summary.get("truth_used_for_profiling_or_comparison") is False and contract.get("truth_used_for_profiling_or_comparison") is False, summary)
    audit.add("Contract", "nonrandom", contract.get("random_process") is False, contract)
    audit.add("Contract", "profile-filters", contract.get("min_read_ani") == 0.95 and contract.get("min_mapq") == 2 and contract.get("pairing_filter") == "paired_only" and contract.get("min_cov") == 5 and contract.get("min_freq") == 0.05 and contract.get("fdr") == 1e-6, contract)
    audit.add("Contract", "stochastic-rarefaction-disabled", contract.get("rarefied_metric_enabled") is False and contract.get("rarefied_coverage") == 1000000000 and "unseeded" in contract.get("rarefaction_audit", ""), contract)
    audit.add("Contract", "comparison-gates", contract.get("compare_database_mode") is True and contract.get("compare_presence_breadth") == 0.5 and contract.get("same_strain_reporting_rule") == {"popANI_min_pct": 99.999, "percent_genome_compared_min_pct": 50}, contract)

    audit.add("Input", "catalog-coordinate", len(genomes) == 24 and {row["SGB"] for row in genomes} == expected_sgbs, len(genomes))
    audit.add("Input", "catalog-hashes", all(re.fullmatch(r"[0-9a-f]{64}", row["RepresentativeSequenceSHA256"]) for row in genomes) and len({row["RepresentativeSequenceSHA256"] for row in genomes}) == 24, len(genomes))
    audit.add("Input", "samples", len(samples) == 2 and [row["Sample"] for row in samples] == ["MOCK1", "MOCK2"] and [row["RunAccession"] for row in samples] == ["ERR9765746", "ERR9765747"], samples)
    audit.add("Input", "bam-hashes", all(re.fullmatch(r"[0-9a-f]{64}", row["BAMSHA256"]) for row in samples), samples)
    audit.add("Input", "large-bam-excluded", not any(frozen.rglob("*.bam")), "BAMs excluded")

    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in profiles:
        by_sample[row["Sample"]].append(row)
    audit.add("Profile", "48-rows", len(profiles) == 48 and set(by_sample) == {"MOCK1", "MOCK2"}, len(profiles))
    for sample, rows in sorted(by_sample.items()):
        audit.add("Profile", f"{sample}-coordinate", len(rows) == 24 and {row["SGB"] for row in rows} == expected_sgbs, len(rows))
        audit.add("Profile", f"{sample}-bounds", all(float(row["Coverage"]) >= 0 and 0 <= float(row["BreadthPct"]) <= 100 and 0 <= float(row["BreadthAtMinCovPct"]) <= 100 and int(row["SNVCount"]) >= 0 for row in rows), sample)
        audit.add("Profile", f"{sample}-presence-formula", all(as_bool(row["PresentBreadth50AtMinCov"]) == (float(row["BreadthAtMinCovPct"]) >= 50) for row in rows), sample)
        audit.add("Profile", f"{sample}-finite-present", all(finite(row["NucleotideDiversity"]) and finite(row["SNVsPerMbpConsidered"]) for row in rows if as_bool(row["PresentBreadth50AtMinCov"])), sample)
        audit.add("Profile", f"{sample}-rarefaction-disabled", all(not finite(row["RarefiedNucleotideDiversity"]) and float(row["RawRarefiedNucleotideDiversity"]) == 0 and row["RarefiedMetricStatus"] == "DISABLED_UNREACHED_COVERAGE_SENTINEL_NORMALIZED_TO_NA" for row in rows), sample)
        audit.add("Profile", f"{sample}-summary-count", int(summary["present_sgbs"][sample]) == sum(as_bool(row["PresentBreadth50AtMinCov"]) for row in rows), summary["present_sgbs"])

    audit.add("Compare", "catalog-coordinate", len(comparisons) == 24 and {row["SGB"] for row in comparisons} == expected_sgbs, len(comparisons))
    compared = [row for row in comparisons if as_bool(row["Compared"])]
    same = [row for row in compared if as_bool(row["SameStrainRule"])]
    audit.add("Compare", "summary-count", len(compared) == int(summary["genomes_compared"]) and len(same) == int(summary["same_strain_calls"]), summary)
    audit.add("Compare", "presence-required", all(as_bool(row["PresentInMOCK1"]) and as_bool(row["PresentInMOCK2"]) for row in compared), len(compared))
    audit.add("Compare", "bounds", all(int(row["ComparedBases"]) > 0 and 0 <= float(row["PercentGenomeCompared"]) <= 100 and 0 <= float(row["PopANIPct"]) <= 100 and 0 <= float(row["ConANIPct"]) <= 100 for row in compared), len(compared))
    audit.add("Compare", "same-strain-formula", all(as_bool(row["SameStrainRule"]) == (float(row["PopANIPct"]) >= 99.999 and float(row["PercentGenomeCompared"]) >= 50) for row in compared), len(same))
    audit.add("Compare", "same-strain-list", [row["SGB"] for row in same] == summary["same_strain_sgbs"], summary["same_strain_sgbs"])

    audit.add("Filtering", "two-samples", len(filtering) == 2 and {row["Sample"] for row in filtering} == {"MOCK1", "MOCK2"}, filtering)
    audit.add("Filtering", "retained-pairs", all(int(row["ReadPairsRetained"]) > 0 and 0 <= float(row["ReadPairsRemovedPct"]) <= 100 and float(row["RetainedGbp"]) > 0 for row in filtering), filtering)
    audit.add("Native", "profile-tables", all(len(read_tsv(frozen / f"raw/{sample}-genome_info.tsv.gz")) == 24 for sample in ("MOCK1", "MOCK2")), "24 rows each")
    audit.add("Native", "compare-table", len(read_tsv(frozen / "raw/MOCK1-vs-MOCK2-genomeWide_compare.tsv.gz")) == len(compared), len(compared))

    command_map = {row["Label"]: row["Command"] for row in commands}
    audit.add("Execution", "five-commands", len(commands) == 5 and all(int(row["ExitStatus"]) == 0 for row in commands), len(commands))
    for sample in ("mock1", "mock2"):
        command = command_map.get(f"instrain-profile-{sample}", "")
        audit.add("Execution", f"{sample}-profile-command", all(token in command for token in ("--min_read_ani 0.95", "--min_mapq 2", "--pairing_filter paired_only", "--min_cov 5", "--min_freq 0.05", "--fdr 1e-06", "--rarefied_coverage 1000000000")), command)
    compare_command = command_map.get("instrain-compare", "")
    audit.add("Execution", "compare-command", all(token in compare_command for token in ("--database_mode", "--breadth 0.5", "--ani_threshold 0.99999", "--coverage_treshold 0.5")), compare_command)
    audit.add("Execution", "resources", len(resources) == 5 and all(int(row["ExitStatus"]) == 0 and float(row["WallSeconds"]) > 0 and float(row["PeakRAMGiB"]) > 0 for row in resources), resources)
    audit.add("Version", "instrain", "1.10.0" in versions.get("inStrain", ""), versions)
    audit.add("Version", "samtools", "1.23.1" in versions.get("SAMtools", ""), versions)
    audit.add("Version", "python", "3.8.20" in versions.get("Python", ""), versions)
    audit.add("Execution", "no-socket-error", all("PermissionError" not in (frozen / "logs" / f"instrain-profile-{sample}.stderr.log").read_text(encoding="utf-8", errors="replace") for sample in ("mock1", "mock2")), "clean elevated run")

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=50,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "catalog_sgbs": 24,
            "samples": 2,
            "genomes_compared": summary["genomes_compared"],
            "same_strain_calls": summary["same_strain_calls"],
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
