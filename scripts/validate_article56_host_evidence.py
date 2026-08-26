#!/usr/bin/env python3
"""Fail-closed validation for Article 56 virus-host evidence."""

from __future__ import annotations

import argparse
import json
import math
import re
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
    "56-evidence-ladder",
    "56-iphop-benchmark-scope",
    "56-confidence-contract",
    "56-negative-control",
    "56-claim-ceiling",
)
EXPECTED_RAW = {
    "PMC10155999.xml": (
        260_284,
        "c2194569972cc1572a0d709e58916c69fd7f1603c7f4645098194d641319e378",
    ),
    "PMC6871006.xml": (
        183_966,
        "289c221ece18f17a0a8930a354cfc3e01b712751abfdd52ff7b054c6f31c36fa",
    ),
}


def near(value: object, expected: float, tolerance: float = 1e-6) -> bool:
    try:
        return math.isclose(float(value), expected, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def valid_tilde_fences(text: str) -> bool:
    open_fence = False
    for line in text.splitlines():
        if not line.startswith("~~~"):
            continue
        if open_fence:
            if line.strip() != "~~~":
                return False
            open_fence = False
        else:
            if line.strip() not in {"~~~bash", "~~~r", "~~~text"}:
                return False
            open_fence = True
    return not open_fence


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot = (frozen / "scripts/plot_article56_host_evidence.R").read_text(
        encoding="utf-8"
    )
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260756" in text,
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
            token in text for token in ("CPU", "RAM", "磁盘", "耗时", "300 GB")
        ),
        "real-data": all(
            token in text
            for token in (
                "PMC10155999.xml",
                "PMC6871006.xml",
                "1,870",
                "216,015",
                "8,128",
                "1,018",
            )
        ),
        "frozen-input": "data/small/56-virus-host-evidence-frozen" in text,
        "historical-current-separation": all(
            token in text
            for token in (
                "iPHoP v1.0",
                "iPHoP v1.4.2",
                "iPHoP_db_Sept21",
                "Sept_2021_pub",
                "env/virus-host-linux-64.lock",
            )
        ),
        "evidence-contract": all(
            token in text
            for token in (
                "prophage",
                "CRISPR spacer",
                "long-read/Hi-C",
                "k-mer/composition",
                "co-abundance",
                "same MAG/bin",
                "List of methods",
                "conflict_flag",
            )
        ),
        "score-contract": all(
            token in text
            for token in (
                "score ≥90",
                "score ≥95",
                "score ≥75",
                "约 ≤10% FDR",
                "经验 PPV",
            )
        ),
        "negative-control-contract": all(
            token in text
            for token in (
                "12.5%",
                "85%",
                "90%",
                "Riboviria",
                "Monodnaviria",
                "域错配负对照",
            )
        ),
        "methods-results": all(
            token in text for token in ("Methods template", "Results template")
        ),
        "citations": all(
            token in text
            for token in (
                "@roux2023iphop",
                "@roux2019miuvig",
                "@edwards2016phagehost",
                "@zhang2021spacepharer",
                "@galiez2017wish",
                "@coutinho2021rafah",
                "@marbouty2017proximity",
                "@bickhart2019assignment",
                "@coenen2018limitations",
                "@iphop2026software",
            )
        ),
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(
            token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')
        ),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot) is None,
        "code-fences": valid_tilde_fences(text),
        "no-source-theme": 'source("R/theme_pub.R")' not in text
        and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(
            r"__[A-Z][A-Z0-9_]*__|\bTODO\b|\bTBD\b|\bNNN\b", text
        )
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
        "article": 56,
        "seed": 20260756,
        "test_viral_genomes": 1870,
        "test_host_genera": 170,
        "temperate_stratum": 949,
        "virulent_stratum": 663,
        "neither_displayed_lifestyle_stratum": 258,
        "imgvr_hq_genomes": 216015,
        "eukaryotic_negative_controls": 8128,
        "eukaryotic_false_host_calls": 1018,
        "false_calls_from_kmer_pct": 85.0,
        "false_calls_below_score90_pct": 90.0,
        "riboviria_false_calls": 640,
        "monodnaviria_false_calls": 155,
        "default_minimum_score": 90,
        "default_estimated_max_fdr_pct": 10,
        "published_runtime_minutes": 12,
        "published_runtime_genomes": 5,
        "published_runtime_cpus": 6,
        "published_iphop_version": "1.0",
        "published_database": "iPHoP_db_Sept21",
        "current_locked_iphop_version": "1.4.2",
        "evidence_tiers": 6,
        "source_assertions": 10,
        "supplementary_spreadsheets_used": False,
        "random_output_requested": False,
    }
    for key, expected in exact.items():
        observed = summary.get(key)
        status = near(observed, expected) if isinstance(expected, float) else observed == expected
        audit.add("Summary", key, status, {"expected": expected, "observed": observed})
    audit.add(
        "Summary",
        "eukaryotic-false-host-rate",
        near(summary.get("eukaryotic_false_host_call_pct"), 12.524606),
        summary.get("eukaryotic_false_host_call_pct"),
    )

    for name, (expected_bytes, expected_hash) in EXPECTED_RAW.items():
        path = frozen / "raw" / name
        exists = path.is_file()
        audit.add("Raw", f"{name}-exists", exists, str(path))
        if exists:
            audit.add(
                "Raw", f"{name}-bytes", path.stat().st_size == expected_bytes,
                {"expected": expected_bytes, "observed": path.stat().st_size},
            )
            audit.add(
                "Raw", f"{name}-sha256", sha256(path) == expected_hash,
                sha256(path),
            )

    assets = read_tsv(frozen / "asset-check-audit.tsv")
    audit.add("Source", "asset-count", len(assets) == 2, len(assets))
    audit.add(
        "Source",
        "all-assets-pass",
        all(as_bool(row["ChecksumPass"]) for row in assets),
        assets,
    )
    assertions = read_tsv(frozen / "source-assertions.tsv")
    audit.add("Source", "assertion-count", len(assertions) == 10, len(assertions))
    audit.add(
        "Source",
        "assertion-ids-unique",
        len({row["AssertionID"] for row in assertions}) == 10,
        [row["AssertionID"] for row in assertions],
    )
    audit.add(
        "Source",
        "all-assertions-pass",
        all(as_bool(row["SourceCheckPass"]) for row in assertions),
        assertions,
    )

    scope = read_tsv(frozen / "benchmark-scope.tsv")
    scope_map = {row["Category"]: int(row["Count"]) for row in scope}
    expected_scope = {
        "All test viruses": 1870,
        "Host genera": 170,
        "Temperate stratum": 949,
        "Virulent stratum": 663,
        "Not in either displayed lifestyle stratum": 258,
        "High-quality prokaryotic-virus genomes": 216015,
        "Eukaryotic-virus genomes": 8128,
    }
    audit.add("Table", "benchmark-scope", scope_map == expected_scope, scope_map)

    hierarchy = read_tsv(frozen / "evidence-hierarchy.tsv")
    audit.add("Table", "evidence-tier-count", len(hierarchy) == 6, len(hierarchy))
    audit.add(
        "Table",
        "evidence-rank-order",
        [int(row["Rank"]) for row in hierarchy] == list(range(1, 7)),
        [row["Evidence"] for row in hierarchy],
    )
    audit.add(
        "Table",
        "prophage-flank-gate",
        hierarchy[0]["Evidence"] == "Integrated prophage + host flanks"
        and "co-binned" in hierarchy[0]["MainFalsePositive"],
        hierarchy[0],
    )

    confidence = read_tsv(frozen / "confidence-contract.tsv")
    confidence_map = {
        int(row["MinimumScore"]): int(row["NominalMaximumFDRPct"])
        for row in confidence
    }
    audit.add(
        "Table", "confidence-contract", confidence_map == {75: 25, 90: 10, 95: 5},
        confidence_map,
    )

    negative = read_tsv(frozen / "negative-control.tsv")
    negative_map = {row["Metric"]: row for row in negative}
    audit.add(
        "Table",
        "negative-control-counts",
        int(negative_map["Eukaryotic viruses tested"]["Count"]) == 8128
        and int(negative_map["Erroneous prokaryotic-host predictions"]["Count"]) == 1018,
        negative_map,
    )
    audit.add(
        "Table",
        "negative-control-percentages",
        near(negative_map["Erroneous prokaryotic-host predictions"]["Percent"], 12.524606)
        and near(negative_map["Errors originating from k-mer comparison"]["Percent"], 85)
        and near(negative_map["Errors with iPHoP score below 90"]["Percent"], 90),
        negative_map,
    )

    components = read_tsv(frozen / "iphop-component-ledger.tsv")
    audit.add("Table", "component-count", len(components) == 6, len(components))
    audit.add(
        "Table",
        "component-membership",
        sum(as_bool(row["IncludedInComposite"]) for row in components) == 4
        and sum(not as_bool(row["IncludedInComposite"]) for row in components) == 2,
        components,
    )
    claims = read_tsv(frozen / "claim-ceiling.tsv")
    audit.add("Table", "claim-matrix-size", len(claims) == 36, len(claims))
    audit.add(
        "Table",
        "claim-status-domain",
        {row["Status"] for row in claims} == {"Allowed", "Conditional", "Avoid"},
        sorted({row["Status"] for row in claims}),
    )

    lock = (frozen / "env/virus-host-linux-64.lock").read_text(encoding="utf-8")
    audit.add(
        "Environment", "iphop-lock", "iphop-1.4.2-pyhdfd78af_0" in lock, "iPHoP 1.4.2"
    )
    audit.add(
        "Environment", "explicit-lock", "@EXPLICIT" in lock and len(lock.splitlines()) == 347,
        len(lock.splitlines()),
    )

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=56,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "frozen_dir": str(frozen),
            "chapter": str(args.chapter.resolve()),
            "figures": list(FIGURES),
            "source_assertions": len(assertions),
            "supplementary_spreadsheets_used": False,
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
