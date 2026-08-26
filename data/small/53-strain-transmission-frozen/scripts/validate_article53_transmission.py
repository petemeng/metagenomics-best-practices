#!/usr/bin/env python3
"""Fail-closed validation for Article 53 strain transmission evidence."""

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
    sha256,
)


FIGURES = (
    "53-species-sharing-negative-control",
    "53-mother-infant-strain-evidence",
    "53-relatedness-classifier",
    "53-fmt-casewise-sharing",
    "53-fmt-source-events",
)
EXPECTED_RAW = {
    "sys001172080st8.xlsx": (
        10_992,
        "8aa0eb1e386eb78b983fac927a1ba1ab2099956570688704d3f9a8a6ca7299a3",
    ),
    "sys001172080st9.xlsx": (
        118_131,
        "12b5dc05b3237b65512c3221812dffef9b749287f46af32ca9d983835c0e9144",
    ),
    "PMC5264247.xml": (
        201_829,
        "1891a7c22bcb3616e77dd94757f1ef6be6b2a04bd7271a0deea1c4a87a8d17bf",
    ),
    "40168_2022_1251_MOESM7_ESM.xlsx": (
        168_598,
        "0f51d05906a2a574970070cab4302204c0f04c2a040bce5d923cd921d91ce757",
    ),
    "PMC8951724.xml": (
        141_847,
        "382a06a553b33b4e6edbcdba900f476986794d699eced2d00ad96caa62219411",
    ),
}


def near(value: object, expected: float, tolerance: float = 1e-8) -> bool:
    return math.isclose(float(value), expected, abs_tol=tolerance)


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot = (frozen / "scripts/plot_article53_transmission.R").read_text(
        encoding="utf-8"
    )
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260753" in text,
        "inline-theme": all(
            token in text
            for token in (
                "pal_pub <-",
                "scale_color_pub <-",
                "scale_fill_pub",
                "theme_pub <-",
                "save_pub <-",
            )
        ),
        "resource-contract": all(
            token in text for token in ("CPU", "RAM", "磁盘", "秒")
        ),
        "real-data-coordinate": all(
            token in text
            for token in ("24 个真实", "27 cases", "408", "294 profiles")
        ),
        "frozen-input": "data/small/53-strain-transmission-frozen" in text,
        "methods-results": all(
            token in text
            for token in (
                "Published-output reanalysis template",
                "Results template",
                "Raw-data SameStr insert",
            )
        ),
        "versions": all(
            token in text
            for token in (
                "Python 3.12.13",
                "openpyxl 3.1.5",
                "lxml 6.0.2",
                "MetaPhlAn2 2.6.0",
                "SAMtools 0.1.19",
                "SameStr 1.2025.111",
            )
        ),
        "database-identity": all(
            token in text
            for token in (
                "db_v20 / mpa_v20_m200",
                "mpa_vJun23_CHOCOPhlAnSGB_202307",
                "ada3ca5a30ce3e0d869fb304182916e302f7ee56",
            )
        ),
        "hard-gates": all(
            token in text
            for token in (
                "MVS ≥ 0.999",
                "≥5,000 bp",
                "--aln-pair-min-overlap 5000",
                "--aln-pair-min-similarity 0.999",
                "--sample-var-min-f-vcov 0.10",
            )
        ),
        "result-counts": all(
            token in text
            for token in (
                "3/7",
                "0.0471",
                "0.0351",
                "0.9337/0.9380",
                "207 donor",
                "119 self",
                "25 both",
                "57 unique",
            )
        ),
        "interpretation-boundaries": all(
            token in text
            for token in (
                "共同环境来源未被排除",
                "transmission-compatible",
                "低覆盖阴性",
                "不是独立患者",
                "隐私",
            )
        ),
        "citations": all(
            token in text
            for token in (
                "@asnicar2017transmission",
                "@podlesny2022samestr",
                "@truong2017strainphlan",
                "@smillie2018engraftment",
                "@samestr2025software",
            )
        ),
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(
            token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')
        ),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot) is None,
        "no-source-theme": 'source("R/theme_pub.R")' not in text
        and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN", text)
        is None,
        "no-meta-prose": not any(
            token in text
            for token in (
                "本篇可独立",
                "本文可独立",
                "全系列约定",
                "接口只学一次",
                "作者代码通常长这样",
                "（即本文）",
                "无头服务器",
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
        checks[f"section-{section}"] = (
            re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
        )
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
    exact = {
        "mother_infant_samples": 24,
        "mother_infant_timepoint_sets": 8,
        "mother_infant_families": 5,
        "mother_infant_named_species": 247,
        "primary_infants_with_negative_controls": 7,
        "primary_matched_ranked_first": 3,
        "quantitative_mother_infant_strain_events": 3,
        "fmt_workbook_samples": 294,
        "fmt_rcdi_samples": 92,
        "fmt_control_samples": 202,
        "fmt_cases": 27,
        "fmt_complete_paired_cases": 25,
        "fmt_competing_strain_events": 408,
        "same_str_overlap_bp": 5000,
        "seed": 20260753,
    }
    for key, expected in exact.items():
        audit.add(
            "Summary",
            key,
            summary.get(key) == expected,
            {"expected": expected, "observed": summary.get(key)},
        )
    numeric = {
        "primary_presence_threshold_pct_exclusive": 0.1,
        "primary_matched_median_shared_species": 2.0,
        "primary_unrelated_median_shared_species": 2.0,
        "primary_matched_median_jaccard": 0.04714913,
        "primary_unrelated_median_jaccard": 0.03513098,
        "fmt_shared_strain_prepost_median": 4.0,
        "fmt_shared_strain_donorpost_median": 15.0,
        "fmt_shared_species_prepost_median": 45.0,
        "fmt_shared_species_donorpost_median": 67.0,
        "same_str_mvs_threshold": 0.999,
    }
    for key, expected in numeric.items():
        audit.add("Summary", key, near(summary.get(key), expected), summary.get(key))
    audit.add(
        "Summary",
        "source-counts",
        summary.get("fmt_source_counts")
        == {"donor": 207, "self": 119, "both": 25, "unique": 57},
        summary.get("fmt_source_counts"),
    )
    audit.add(
        "Summary",
        "no-random-output",
        summary.get("random_output_requested") is False,
        summary.get("random_output_requested"),
    )

    assets = read_tsv(frozen / "asset-check-audit.tsv")
    audit.add(
        "Input",
        "five-checksum-gated-assets",
        len(assets) == 5 and all(as_bool(row["ChecksumPass"]) for row in assets),
        len(assets),
    )
    raw_results = []
    for name, (expected_bytes, expected_hash) in EXPECTED_RAW.items():
        path = frozen / "raw" / name
        raw_results.append(
            path.is_file()
            and path.stat().st_size == expected_bytes
            and sha256(path) == expected_hash
        )
    audit.add("Input", "raw-source-snapshots", all(raw_results), raw_results)
    manifest = read_tsv(frozen / "asset-manifest.tsv")
    audit.add(
        "Input",
        "open-licenses-and-dois",
        len(manifest) == 5
        and all(row["License"] == "CC BY 4.0" for row in manifest)
        and {row["DOI"] for row in manifest}
        == {"10.1128/mSystems.00164-16", "10.1186/s40168-022-01251-w"},
        Counter(row["DOI"] for row in manifest),
    )

    design = read_tsv(frozen / "mother-infant-design.tsv")
    audit.add(
        "MotherInfant",
        "design-coordinate",
        len(design) == 8
        and {int(row["Pair"]) for row in design} == {1, 2, 3, 4, 5}
        and Counter(row["TimePoint"] for row in design)
        == Counter({"T1": 5, "T2": 2, "T3": 1}),
        len(design),
    )
    sample_meta = read_tsv(frozen / "mother-infant-sample-metadata.tsv")
    audit.add(
        "MotherInfant",
        "sample-coordinate",
        len(sample_meta) == 24
        and Counter(row["SampleType"] for row in sample_meta)
        == Counter({"Infant": 8, "Mother": 8, "Milk": 8}),
        Counter(row["SampleType"] for row in sample_meta),
    )
    species_profile = read_tsv(frozen / "mother-infant-species-profile.tsv.gz")
    audit.add(
        "MotherInfant",
        "profile-coordinate",
        len(species_profile) == 247 * 24
        and len({row["Species"] for row in species_profile}) == 247
        and len({row["Sample"] for row in species_profile}) == 24,
        len(species_profile),
    )
    audit.add(
        "MotherInfant",
        "profile-values",
        all(float(row["RelativeAbundancePct"]) >= 0 for row in species_profile)
        and not any(row["Species"].endswith("_unclassified") for row in species_profile),
        "non-negative named species",
    )
    pairwise = read_tsv(frozen / "mother-infant-pairwise-sharing.tsv")
    primary_pairs = [
        row for row in pairwise if near(row["ThresholdPctExclusive"], 0.1)
    ]
    audit.add(
        "MotherInfant",
        "pairwise-coordinate",
        len(pairwise) == 120
        and len(primary_pairs) == 30
        and sum(as_bool(row["MatchedPair"]) for row in primary_pairs) == 8
        and sum(not as_bool(row["MatchedPair"]) for row in primary_pairs) == 22,
        [len(pairwise), len(primary_pairs)],
    )
    ranks = read_tsv(frozen / "mother-infant-rank-sensitivity.tsv")
    primary_ranks = [
        row
        for row in ranks
        if near(row["ThresholdPctExclusive"], 0.1)
        and int(row["UnrelatedNegativeControls"]) > 0
    ]
    audit.add(
        "MotherInfant",
        "rank-negative-control",
        len(ranks) == 32
        and len(primary_ranks) == 7
        and sum(as_bool(row["UniquelyRankedFirst"]) for row in primary_ranks) == 3,
        len(primary_ranks),
    )
    sharing = {
        float(row["ThresholdPctExclusive"]): row
        for row in read_tsv(frozen / "mother-infant-sharing-summary.tsv")
    }
    audit.add(
        "MotherInfant",
        "four-thresholds",
        set(sharing) == {0.0, 0.01, 0.1, 1.0}
        and int(sharing[0.1]["MatchedRankedFirstIncludingTies"]) == 3
        and int(sharing[1.0]["MatchedUniquelyRankedFirst"]) == 2,
        sharing,
    )
    strain = read_tsv(frozen / "published-mother-infant-strain-evidence.tsv")
    strain_map = {row["Species"]: row for row in strain}
    audit.add(
        "MotherInfant",
        "published-strain-evidence",
        len(strain) == 3
        and near(strain_map["Bifidobacterium bifidum"]["IntraPairDivergencePct"], 0.04)
        and near(strain_map["Coprococcus comes"]["ClosestOtherDivergencePct"], 1.60)
        and near(strain_map["Ruminococcus bromii"]["ClosestOtherDivergencePct"], 1.53)
        and all(row["IndependentGeneContentEvidence"] == "PanPhlAn Figure S5" for row in strain),
        strain_map,
    )
    ledger = read_tsv(frozen / "mother-infant-evidence-ledger.tsv")
    audit.add(
        "MotherInfant",
        "evidence-boundary",
        len(ledger) == 6
        and all(as_bool(row["MarkerSNVTree"]) for row in ledger)
        and all(not as_bool(row["SharedEnvironmentExcluded"]) for row in ledger),
        len(ledger),
    )

    fmt_samples = read_tsv(frozen / "fmt-sample-metadata.tsv.gz")
    audit.add(
        "FMT",
        "sample-coordinate",
        len(fmt_samples) == 294
        and Counter(row["StudyType"] for row in fmt_samples)
        == Counter({"Control": 202, "rCDI": 92}),
        Counter(row["StudyType"] for row in fmt_samples),
    )
    cases = read_tsv(frozen / "fmt-casewise-sharing.tsv")
    audit.add(
        "FMT",
        "case-coordinate",
        len(cases) == 27
        and Counter(row["FMTOutcome"] for row in cases)
        == Counter({"Resolved": 21, "Failed": 6}),
        Counter(row["FMTOutcome"] for row in cases),
    )
    long_cases = read_tsv(frozen / "fmt-casewise-sharing-long.tsv")
    audit.add("FMT", "long-coordinate", len(long_cases) == 104, len(long_cases))
    paired = {
        row["Resolution"]: row
        for row in read_tsv(frozen / "fmt-casewise-sharing-summary.tsv")
    }
    audit.add(
        "FMT",
        "paired-strain-summary",
        int(paired["Strain"]["CompleteCases"]) == 25
        and near(paired["Strain"]["PrePostMedian"], 4)
        and near(paired["Strain"]["DonorPostMedian"], 15)
        and [int(paired["Strain"][key]) for key in ("DonorGreater", "Equal", "DonorLess")]
        == [19, 1, 5]
        and near(paired["Strain"]["ExactTwoSidedSignP"], 0.006610751152, 1e-12),
        paired["Strain"],
    )
    audit.add(
        "FMT",
        "paired-species-summary",
        near(paired["Species"]["PrePostMedian"], 45)
        and near(paired["Species"]["DonorPostMedian"], 67)
        and [int(paired["Species"][key]) for key in ("DonorGreater", "Equal", "DonorLess")]
        == [22, 0, 3],
        paired["Species"],
    )
    classifier = {
        (row["TaxonomicLevel"], row["TestSet"]): row
        for row in read_tsv(frozen / "relatedness-classifier-performance.tsv")
    }
    audit.add(
        "FMT",
        "classifier-performance",
        len(classifier) == 8
        and near(classifier[("Species", "Control hold-out")]["AUPR"], 0.4745)
        and near(classifier[("Strain", "Control hold-out")]["AUROC"], 1.0)
        and near(classifier[("Strain", "rCDI / FMT")]["AUROC"], 0.9337)
        and near(classifier[("Strain", "rCDI / FMT")]["AUPR"], 0.938),
        list(classifier.values()),
    )
    events = read_tsv(frozen / "fmt-competing-strain-events.tsv.gz")
    audit.add(
        "FMT",
        "event-coordinate",
        len(events) == 408
        and Counter(row["Source"] for row in events)
        == Counter({"donor": 207, "self": 119, "unique": 57, "both": 25}),
        Counter(row["Source"] for row in events),
    )
    source_summary = read_tsv(frozen / "fmt-source-event-summary.tsv")
    audit.add(
        "FMT",
        "outcome-event-denominators",
        len(source_summary) == 8
        and {row["FMTOutcome"]: int(row["TotalEvents"]) for row in source_summary}
        == {"Resolved": 314, "Failed": 94},
        source_summary,
    )
    pure = read_tsv(frozen / "fmt-pure-source-gate-audit.tsv.gz")
    audit.add(
        "FMT",
        "pure-source-hard-gates",
        len(pure) == 383
        and Counter(row["Source"] for row in pure)
        == Counter({"donor": 207, "self": 119, "unique": 57})
        and all(as_bool(row["ExpectedPureCategory"]) for row in pure),
        Counter(row["Source"] for row in pure),
    )
    ladder = read_tsv(frozen / "transmission-evidence-ladder.tsv")
    audit.add(
        "Output",
        "five-level-evidence-ladder",
        len(ladder) == 5
        and [int(row["Level"]) for row in ladder] == [1, 2, 3, 4, 5],
        len(ladder),
    )

    determinism = read_tsv(frozen / "determinism-audit.tsv")
    audit.add(
        "Execution",
        "byte-identical-replay",
        len(determinism) == 18
        and all(
            as_bool(row["ByteIdentical"])
            and row["Seed"] == "20260753"
            and not as_bool(row["RandomOutputRequested"])
            for row in determinism
        ),
        len(determinism),
    )
    commands = read_tsv(frozen / "command-log.tsv")
    audit.add(
        "Execution",
        "two-successful-runs",
        len(commands) == 2 and all(int(row["ExitStatus"]) == 0 for row in commands),
        len(commands),
    )
    resources = read_tsv(frozen / "resource-summary.tsv")
    audit.add(
        "Execution",
        "measured-resources",
        len(resources) == 2
        and all(
            int(row["ExitStatus"]) == 0
            and float(row["WallSeconds"]) > 0
            and float(row["PeakRAMGiB"]) > 0
            for row in resources
        ),
        resources,
    )
    versions = {
        row["Software"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")
    }
    audit.add(
        "Version",
        "parser-stack",
        versions
        == {"Python": "3.12.13", "openpyxl": "3.1.5", "lxml": "6.0.2"},
        versions,
    )
    plot_versions = {
        row["Software"]: row["Version"]
        for row in read_tsv(frozen / "plot-software-versions.tsv")
    }
    audit.add(
        "Version",
        "plot-stack",
        plot_versions
        == {
            "R": "4.4.1",
            "ggplot2": "3.5.2",
            "dplyr": "1.1.4",
            "readr": "2.1.5",
            "tidyr": "1.3.1",
            "scales": "1.3.0",
            "patchwork": "1.3.2",
        },
        plot_versions,
    )

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=53,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "mother_infant_profiles": summary["mother_infant_samples"],
            "matched_mother_ranked_first": summary["primary_matched_ranked_first"],
            "fmt_cases": summary["fmt_cases"],
            "fmt_events": summary["fmt_competing_strain_events"],
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
