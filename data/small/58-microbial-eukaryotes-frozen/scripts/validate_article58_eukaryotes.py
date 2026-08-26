#!/usr/bin/env python3
"""Fail-closed validation for Article 58 microbial-eukaryote evidence."""

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
    "58-eukdetect-marker-evidence",
    "58-eukdetect-abundance-denominator",
    "58-eukrep-reference-benchmark",
    "58-eukrep-assembly-audit",
    "58-eukaryote-evidence-ladder",
)
LENGTHS = {"3000", "5000", "10000", "20000"}
MODES = {"Strict", "Balanced", "Lenient"}


def near(value: object, expected: float, tolerance: float = 1e-6) -> bool:
    try:
        return math.isclose(float(value), expected, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def printable_counter(counter: Counter[object]) -> dict[str, int]:
    return {
        "|".join(str(part) for part in key) if isinstance(key, tuple) else str(key): value
        for key, value in counter.items()
    }


def valid_tilde_fences(text: str) -> bool:
    opened = False
    for line in text.splitlines():
        if not line.startswith("~~~"):
            continue
        if opened:
            if line.strip() != "~~~":
                return False
            opened = False
        else:
            if line.strip() not in {"~~~bash", "~~~r", "~~~text", "~~~python"}:
                return False
            opened = True
    return not opened


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot = (frozen / "scripts/plot_article58_eukaryotes.R").read_text(encoding="utf-8")
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260758" in text,
        "inline-theme": all(
            token in text
            for token in (
                "pal_pub <-", "scale_color_pub <-", "scale_fill_pub",
                "theme_pub <-", "save_pub <-",
            )
        ),
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "耗时")),
        "real-data": all(
            token in text
            for token in (
                "SRR12324253", "PRJNA648136", "10.5281/zenodo.3935737",
                "1485049c9792c7e43d267fe7bb84a2bd",
                "fbd8494da79fc796d6725a4e242a9b9c",
            )
        ),
        "frozen-input": "data/small/58-microbial-eukaryotes-frozen" in text,
        "version-contract": all(
            token in text
            for token in (
                "2.0.2", "2.0.0", "v2.0.1", "legacy-cgi 2.6.4",
                "EukRep 0.6.7", "MEGAHIT 1.2.9", "minimap2 2.31-r1302",
                "2026-03-16", "10.5281/zenodo.19056625",
            )
        ),
        "eukdetect-boundary": all(
            token in text
            for token in (
                "RPKS", "RPKSB", "RelEuk", "2 个 marker", "4 条 reads",
                "不是细胞丰度", "不能证明活性", "不能证明定植",
            )
        ),
        "eukrep-boundary": all(
            token in text
            for token in (
                "--min 3000", "--tie prok", "5 kb chunk", "exit code 0",
                "不提供物种分类", "不提供丰度",
            )
        ),
        "methods-results": all(token in text for token in ("Methods template", "Results template")),
        "citations": all(
            token in text
            for token in (
                "@lind2021eukdetect", "@shih2026eukdetect2",
                "@west2018eukrep", "@zymo2020refgenomes",
                "@zymo2026d6300", "@crouch2024eukaryome",
                "@manni2021busco", "@levykarin2020metaeuk",
                "@saary2020eukcc",
            )
        ),
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(
            token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')
        ),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot) is None,
        "code-fences": valid_tilde_fences(text),
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z][A-Z0-9_]*__|\bTODO\b|\bTBD\b|\bNNN\b", text) is None,
        "no-meta-prose": not any(
            token in text
            for token in (
                "本篇可独立", "本文可独立", "全系列约定", "接口只学一次",
                "作者代码通常长这样", "（即本文）", "无头服务器",
            )
        ),
    }
    sections = (
        "这一步对应论文里的哪张图",
        "理论：",
        "准备工作",
        "可复制代码",
        "审计与升级",
        "出版级美化",
        "常见坑",
        "这段 Methods 怎么写",
        "换成你自己的数据怎么做",
        "参考",
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

    run = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    exact_run = {
        "article": 58,
        "seed": 20260758,
        "run_accession": "SRR12324253",
        "bioproject": "PRJNA648136",
        "eukdetect_pairs": 1000000,
        "assembly_pairs": 20000000,
        "eukdetect_total_bases_pre_qc": 299110793,
        "assembly_total_bases_pre_qc": 5982154152,
        "reference_doi": "10.5281/zenodo.3935737",
        "reference_genomes": 10,
        "fragments_per_species_length": 80,
        "expected_eukaryotes": 2,
        "random_output_requested": False,
    }
    for key, expected in exact_run.items():
        observed = run.get(key)
        audit.add("Run contract", key, observed == expected, {"expected": expected, "observed": observed})
    audit.add("Run contract", "fragment-lengths", set(run["fragment_lengths"]) == {3000, 5000, 10000, 20000}, run["fragment_lengths"])

    assets = read_tsv(frozen / "asset-check-audit.tsv")
    audit.add("Input", "asset-count", len(assets) == 3, len(assets))
    audit.add("Input", "asset-checksums", all(as_bool(row["ChecksumPass"]) for row in assets), assets)
    audit.add(
        "Input", "asset-roles",
        {row["Role"] for row in assets} == {"zymo-reference-v2", "srr12324253-r1-fastq", "srr12324253-r2-fastq"},
        sorted(row["Role"] for row in assets),
    )

    database = read_tsv(frozen / "database-audit.tsv")
    audit.add("Database", "file-count", len(database) == 14, len(database))
    audit.add("Database", "byte-total", sum(int(row["ObservedBytes"]) for row in database) == 7123869920, sum(int(row["ObservedBytes"]) for row in database))
    audit.add("Database", "checksums", all(as_bool(row["ChecksumPass"]) for row in database), database)
    audit.add("Database", "release", {row["Release"] for row in database} == {"2026-03-16"} and {row["RecordDOI"] for row in database} == {"10.5281/zenodo.19056625"}, sorted({row["Release"] for row in database}))

    prepared = read_tsv(frozen / "prepared-fastq-audit.tsv")
    prepared_pairs = {row["Role"]: int(row["Pairs"]) for row in prepared}
    audit.add("Input", "prepared-fastq-count", len(prepared) == 4, len(prepared))
    audit.add(
        "Input", "prepared-pair-counts",
        prepared_pairs == {
            "assembly_r1": 20000000, "assembly_r2": 20000000,
            "eukdetect_r1": 1000000, "eukdetect_r2": 1000000,
        }, prepared_pairs,
    )

    references = read_tsv(frozen / "reference-sequence-ledger.tsv")
    reference_domains = Counter(row["DomainTruth"] for row in references)
    audit.add("Input", "reference-sequence-count", len(references) == 7218, len(references))
    audit.add("Input", "reference-domain-counts", reference_domains == Counter({"Eukaryote": 7205, "Prokaryote": 13}), reference_domains)
    audit.add("Input", "reference-species-count", len({row["Species"] for row in references}) == 10, sorted({row["Species"] for row in references}))

    fragments = read_tsv(frozen / "eukrep-fragment-ledger.tsv")
    fragment_domains = Counter(row["DomainTruth"] for row in fragments)
    length_counts = Counter(row["FragmentLength"] for row in fragments)
    species_length = Counter((row["Species"], row["FragmentLength"]) for row in fragments)
    audit.add("Input", "fragment-count", len(fragments) == 3200, len(fragments))
    audit.add("Input", "fragment-domains", fragment_domains == Counter({"Prokaryote": 2560, "Eukaryote": 640}), fragment_domains)
    audit.add("Input", "fragment-length-balance", length_counts == Counter({length: 800 for length in LENGTHS}), length_counts)
    audit.add("Input", "fragment-species-balance", len(species_length) == 40 and set(species_length.values()) == {80}, printable_counter(species_length))

    summary = json.loads((frozen / "summary.json").read_text(encoding="utf-8"))
    exact_summary = {
        "run_accession": "SRR12324253",
        "eukdetect_species_calls": 3,
        "eukdetect_expected_species_detected": 2,
        "eukdetect_unexpected_species_calls": 0,
        "reference_benchmark_fragments": 3200,
        "assembly_contigs_ge1000": 4222,
        "assembly_total_bases": 37166382,
        "assembly_n50": 122804,
        "assembly_contigs_ge3000": 691,
        "assembly_truth_resolved": 609,
        "assembly_truth_eukaryote": 221,
        "assembly_truth_prokaryote": 388,
        "assembly_truth_unresolved": 82,
        "measured_commands": 21,
    }
    for key, expected in exact_summary.items():
        audit.add("Summary", key, summary.get(key) == expected, {"expected": expected, "observed": summary.get(key)})

    evidence = read_tsv(frozen / "eukdetect-species-evidence.tsv")
    by_name = {row["Name"]: row for row in evidence}
    exact_markers = {
        "Saccharomyces cerevisiae": (60, 235.0),
        "Cryptococcus deneoformans": (33, 77.0),
        "Cryptococcus neoformans": (22, 53.0),
    }
    audit.add("EukDetect", "species-set", set(by_name) == set(exact_markers), sorted(by_name))
    audit.add(
        "EukDetect", "marker-evidence",
        all(int(by_name[name]["ObservedMarkers"]) == expected[0] and near(by_name[name]["TotalMarkerReads"], expected[1]) for name, expected in exact_markers.items()),
        {name: (row["ObservedMarkers"], row["TotalMarkerReads"]) for name, row in by_name.items()},
    )
    audit.add("EukDetect", "expected-only", all(as_bool(row["ExpectedInMock"]) for row in evidence), evidence)
    audit.add("EukDetect", "releuk-sum", near(sum(float(row["RelEukPercent"]) for row in evidence), 99.9999, 1e-4), sum(float(row["RelEukPercent"]) for row in evidence))
    detection = read_tsv(frozen / "eukdetect-detection-audit.tsv")
    audit.add("EukDetect", "two-mock-targets", len(detection) == 2 and all(as_bool(row["DetectedAtSpeciesRank"]) for row in detection), detection)

    reference_calls = read_tsv(frozen / "eukrep-reference-calls.tsv")
    reference_metrics = read_tsv(frozen / "eukrep-reference-metrics.tsv")
    call_cells = Counter((row["FragmentLength"], row["Mode"]) for row in reference_calls)
    audit.add("EukRep", "reference-call-count", len(reference_calls) == 9600, len(reference_calls))
    audit.add("EukRep", "reference-complete-partitions", len(call_cells) == 12 and set(call_cells.values()) == {800}, printable_counter(call_cells))
    audit.add("EukRep", "reference-no-unclassified", {row["Prediction"] for row in reference_calls} <= {"Eukaryote", "Prokaryote"}, Counter(row["Prediction"] for row in reference_calls))
    audit.add("EukRep", "reference-metric-count", len(reference_metrics) == 48, len(reference_metrics))
    audit.add("EukRep", "reference-metric-grid", {row["FragmentLength"] for row in reference_metrics} == LENGTHS and {row["Mode"] for row in reference_metrics} == MODES and {row["Metric"] for row in reference_metrics} == {"Sensitivity", "Specificity", "Precision", "Classification rate"}, "4 lengths x 3 modes x 4 metrics")

    assembly_calls = read_tsv(frozen / "eukrep-assembly-calls.tsv")
    assembly_metrics = read_tsv(frozen / "eukrep-assembly-metrics.tsv")
    mode_counts = Counter(row["Mode"] for row in assembly_calls)
    truth_counts = Counter((row["Mode"], row["Truth"]) for row in assembly_calls)
    audit.add("Assembly", "call-count", len(assembly_calls) == 2073 and mode_counts == Counter({mode: 691 for mode in MODES}), mode_counts)
    audit.add("Assembly", "truth-counts", all(truth_counts[(mode, "Eukaryote")] == 221 and truth_counts[(mode, "Prokaryote")] == 388 and truth_counts[(mode, "Unresolved")] == 82 for mode in MODES), printable_counter(truth_counts))
    audit.add("Assembly", "query-coverage-domain", all(0 <= float(row["QueryCoveragePercent"]) <= 100 for row in assembly_calls), max(float(row["QueryCoveragePercent"]) for row in assembly_calls))
    audit.add("Assembly", "prediction-domain", {row["Prediction"] for row in assembly_calls} <= {"Eukaryote", "Prokaryote"}, Counter(row["Prediction"] for row in assembly_calls))
    audit.add("Assembly", "metric-grid", len(assembly_metrics) == 9 and {row["Metric"] for row in assembly_metrics} == {"Sensitivity", "Specificity", "Precision"}, len(assembly_metrics))
    selected = {(row["Mode"], row["Metric"]): float(row["Estimate"]) for row in assembly_metrics}
    audit.add("Assembly", "balanced-performance", near(selected[("Balanced", "Sensitivity")], 163 / 221) and near(selected[("Balanced", "Specificity")], 386 / 388), {"|".join(key): value for key, value in selected.items()})

    ladder = read_tsv(frozen / "evidence-ladder.tsv")
    audit.add("Evidence", "ladder", len(ladder) == 5 and [int(row["Rank"]) for row in ladder] == list(range(1, 6)), ladder)
    audit.add("Evidence", "claim-ceiling", ladder[0]["ClaimCeiling"] == "Taxon detected" and ladder[-1]["ClaimCeiling"] == "Persistence or colonization supported", ladder)

    resources = read_tsv(frozen / "resource-usage.tsv")
    labels = {row["Label"] for row in resources}
    audit.add("Resource", "measured-command-count", len(resources) == 21 and len(labels) == 21, sorted(labels))
    audit.add("Resource", "positive-runtime-memory", all(float(row["WallSeconds"]) > 0 and float(row["PeakRSSGiB"]) > 0 for row in resources), resources)
    megahit = next(row for row in resources if row["Label"] == "megahit")
    audit.add("Resource", "megahit-measured", near(megahit["WallSeconds"], 1599.23, 0.1) and near(megahit["PeakRSSGiB"], 7.0805, 1e-4), megahit)

    versions = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    exact_versions = {
        "EukDetect conda record": "2.0.2",
        "EukDetect Python metadata": "2.0.0",
        "EukDetect CLI": "EukDetect v2.0.1",
        "legacy-cgi compatibility package": "2.6.4",
        "EukRep": "EukRep 0.6.7",
        "Snakemake": "9.14.5",
        "fastp": "fastp 1.3.6",
        "MEGAHIT": "MEGAHIT v1.2.9",
        "minimap2": "2.31-r1302",
        "EukDetect2 database": "2026-03-16 (doi:10.5281/zenodo.19056625)",
    }
    audit.add("Environment", "exact-versions", all(versions.get(tool) == version for tool, version in exact_versions.items()), versions)
    audit.add("Environment", "bowtie-version", versions.get("Bowtie2", "").endswith("version 2.5.4"), versions.get("Bowtie2"))
    locks = [
        (frozen / "env/microbial-eukaryotes-linux-64.lock").read_text(encoding="utf-8"),
        (frozen / "env/eukrep-legacy-linux-64.lock").read_text(encoding="utf-8"),
        (frozen / "env/assembly-linux-64.lock").read_text(encoding="utf-8"),
        (frozen / "env/read-qc-linux-64.lock").read_text(encoding="utf-8"),
    ]
    audit.add("Environment", "explicit-locks", all("@EXPLICIT" in lock for lock in locks), "four Linux exact locks")
    pip_lock = (frozen / "env/microbial-eukaryotes-pip-lock.txt").read_text(encoding="utf-8")
    audit.add("Environment", "hashed-legacy-cgi", "legacy-cgi==2.6.4" in pip_lock and "sha256:7e235ce58bf1e25d1fc9b2d299015e4e2cd37305eccafec1e6bac3fc04b878cd" in pip_lock, pip_lock.strip())

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=58,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "frozen_dir": str(frozen),
            "chapter": str(args.chapter.resolve()),
            "figures": list(FIGURES),
            "read_branch": "EukDetect2 on first 1,000,000 read pairs",
            "assembly_branch": "MEGAHIT + EukRep on first 20,000,000 read pairs",
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
