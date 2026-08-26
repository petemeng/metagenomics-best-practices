#!/usr/bin/env python3
"""Fail-closed validation for Article 61 community-metabolism evidence."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image

from article42_44_validation_utils import (
    Audit,
    audit_checksums,
    audit_figures,
    finish,
    read_tsv,
    sha256,
)


FIGURES = (
    "61-model-coverage",
    "61-tradeoff-curves",
    "61-constraint-sensitivity",
    "61-net-community-flux",
    "61-micom-crossfeeding",
    "61-smetana-concordance",
)

SUMMARY_FILES = (
    "sample-selection-audit.tsv",
    "selected-samples.tsv",
    "model-coverage.tsv",
    "taxon-match-audit.tsv",
    "tradeoff-summary.tsv",
    "primary-growth.tsv",
    "primary-growth-summary.tsv",
    "medium-sensitivity-summary.tsv",
    "abundance-sensitivity-summary.tsv",
    "primary-exchanges.tsv",
    "net-community-flux.tsv",
    "micom-crossfeeding-potential.tsv",
    "focal-micom-fluxes.tsv",
    "focal-micom-potential-edges.tsv",
    "smetana-global.tsv",
    "smetana-detailed.tsv",
    "smetana-compatibility-audit.tsv",
    "smetana-component-summary.tsv",
    "smetana-pair-summary.tsv",
    "cross-method-concordance.tsv",
    "run-ledger.tsv",
    "analysis-metrics.json",
)

RAW_FOR_REANALYSIS = (
    "selected-samples.tsv",
    "model-coverage.tsv",
    "sample-selection-audit.tsv",
    "taxon-match-audit.tsv",
    "micom-tradeoff-growth.tsv",
    "micom-primary-growth.tsv",
    "micom-primary-exchanges.tsv",
    "micom-medium-sensitivity.tsv",
    "micom-equal-abundance-growth.tsv",
    "smetana-subcommunity.tsv",
    "smetana-medium-compatibility.tsv",
    "smetana-compatibility-audit.tsv",
    "smetana-compatibility-summary.json",
    "smetana-global.tsv",
    "smetana-detailed.tsv",
    "run-ledger.tsv",
)

EXPECTED_RESOURCES = {
    "AGORA201RefSeq216Species": (
        "264422408",
        "8fb2ade5b970aadefd7e3c381b2fca3854e32082361507cf0e6f8dc32a3394dc",
    ),
    "AGORA2RefSeqSpeciesManifest": (
        "4775059",
        "fe1c7fd3ede1c2f96097c311294ecc674d7ca9320c0f806b619a4b83337905c3",
    ),
    "WesternDietGutAGORA": (
        "9047",
        "a53d83b9c892e6c11dbb6ec22ab46fda136bebe7d0ecd14d4c5b716e1c5574fe",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stage-dir", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    return parser.parse_args()


def near(value: object, expected: float, tolerance: float = 1e-8) -> bool:
    try:
        return math.isclose(float(value), expected, abs_tol=tolerance, rel_tol=tolerance)
    except (TypeError, ValueError):
        return False


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def truthy(value: object) -> bool:
    return str(value).lower() == "true"


def image_pixel_sha(path: Path) -> str:
    with Image.open(path) as image:
        normalized = image.convert("RGBA")
        payload = (
            normalized.size[0].to_bytes(8, "little")
            + normalized.size[1].to_bytes(8, "little")
            + normalized.tobytes()
        )
    return hashlib.sha256(payload).hexdigest()


def audit_chapter_61(chapter: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    sections = (
        "这一步对应论文里的哪张图",
        "理论",
        "准备工作",
        "可复制代码",
        "审计与升级",
        "出版级美化",
        "常见坑",
        "这段 Methods 怎么写",
        "换成你自己的数据怎么做",
        "参考",
    )
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "analysis-seed": "61001" in text,
        "plot-seed": "20260761" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("RAM", "CPU", "磁盘", "小时")),
        "real-data": all(token in text for token in ("AsnicarF_2017", "SRR4052043", "curatedMetagenomicData")),
        "database-lock": all(token in text for token in ("AGORA 2.01", "8fb2ade5", "1,746")),
        "methods-template": "**Personalized community metabolic modeling.**" in text,
        "frozen-input": "data/small/61-community-metabolism-frozen" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "required-science": all(token in text for token in ("cooperative trade-off", "pFBA", "MRO", "MIP", "SCS", "MUS", "MPS")),
        "smetana-fail-closed": all(token in text for token in (
            "not estimable", "interacting=False", "7,853", "not evidence of biological absence",
        )),
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-meta-prose": not any(token in text for token in (
            "本篇可独立", "本文可独立", "全系列约定", "接口只学一次",
            "作者代码通常长这样", "（即本文）", "无头服务器",
        )),
    }
    for section in sections:
        checks[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def stage_reanalysis(
    root: Path,
    frozen: Path,
    stage: Path,
    python: Path,
    audit: Audit,
) -> None:
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    for name in RAW_FOR_REANALYSIS:
        shutil.copy2(frozen / "raw" / name, stage / name)
    summary = stage / "summary"
    environment = os.environ.copy()
    environment["MPLCONFIGDIR"] = str(stage / "matplotlib")

    summarized = subprocess.run(
        [
            str(python),
            str(root / "scripts/summarize_article61_community.py"),
            "--work-dir",
            str(stage),
            "--summary-dir",
            str(summary),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    audit.add(
        "Reanalysis",
        "summary-exit",
        summarized.returncode == 0,
        (summarized.stdout + summarized.stderr)[-4000:],
    )
    if summarized.returncode != 0:
        return
    for name in SUMMARY_FILES:
        audit.add(
            "Reanalysis",
            f"summary-hash-{name}",
            sha256(summary / name) == sha256(frozen / name),
            sha256(summary / name),
        )

    staged_figures = stage / "figures"
    plotted = subprocess.run(
        [
            str(python),
            str(root / "scripts/plot_article61_community.py"),
            "--summary-dir",
            str(summary),
            "--figure-dir",
            str(staged_figures),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    audit.add(
        "Reanalysis",
        "plot-exit",
        plotted.returncode == 0,
        (plotted.stdout + plotted.stderr)[-4000:],
    )
    if plotted.returncode != 0:
        return
    audit_figures(staged_figures, audit, FIGURES)
    for stem in FIGURES:
        staged = staged_figures / f"{stem}.png"
        published = root / "figures" / f"{stem}.png"
        audit.add(
            "Reanalysis",
            f"pixel-match-{stem}",
            staged.is_file()
            and published.is_file()
            and image_pixel_sha(staged) == image_pixel_sha(published),
            image_pixel_sha(staged) if staged.is_file() else "MISSING",
        )


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    figures = args.figure_dir.resolve()
    audit = Audit()

    audit_checksums(frozen, audit)
    contract = json.loads((frozen / "frozen-contract.json").read_text(encoding="utf-8"))
    audit.add("Contract", "article", contract.get("article") == 61, contract)
    audit.add("Contract", "sample-count", contract.get("selected_samples") == 6, contract)
    audit.add("Contract", "subject-count", contract.get("independent_subjects") == 6, contract)
    audit.add("Contract", "smetana-size", contract.get("smetana_subcommunity_size") == 6, contract)
    audit.add("Contract", "smetana-global-limited", contract.get("smetana_global_estimable") is False, contract)
    audit.add("Contract", "smetana-comparison-limited", contract.get("smetana_cross_method_comparison_estimable") is False, contract)
    audit.add("Contract", "seed", contract.get("seed") == 61001, contract)
    audit.add("Contract", "tradeoff", near(contract.get("primary_tradeoff"), 0.5), contract)
    audit.add("Contract", "large-db-excluded", contract.get("agora_qza_included") is False and contract.get("micom_pickles_included") is False, contract)
    audit.add("Contract", "source-profile-included", contract.get("source_profile_included") is True, contract)

    resources = read_tsv(frozen / "database-resource-manifest.tsv")
    by_asset = {row["Asset"]: row for row in resources}
    audit.add("Database", "resource-count", len(resources) == 3 and len(by_asset) == 3, len(resources))
    for asset, (size, digest) in EXPECTED_RESOURCES.items():
        row = by_asset.get(asset, {})
        audit.add(
            "Database",
            f"locked-{asset}",
            row.get("ExpectedBytes") == size and row.get("SHA256") == digest,
            row,
        )

    internal_manifest = read_tsv(frozen / "raw/database-manifest.tsv")
    audit.add("Database", "species-model-count", len(internal_manifest) == 1746, len(internal_manifest))
    audit.add("Database", "unique-species", len({row["species"] for row in internal_manifest}) == 1746, len(internal_manifest))
    audit.add("Database", "summary-rank", {row["summary_rank"] for row in internal_manifest} == {"species"}, sorted({row["summary_rank"] for row in internal_manifest}))

    selected = read_tsv(frozen / "selected-samples.tsv")
    audit.add("Input", "selected-six", len(selected) == 6, selected)
    audit.add("Input", "independent-subjects", len({row["SubjectID"] for row in selected}) == 6, selected)
    audit.add("Input", "read-gate", all(int(row["Reads"]) >= 1_000_000 for row in selected), selected)
    audit.add("Input", "adult-only", {row["age_category"] for row in selected} == {"adult"}, selected)
    audit.add("Input", "accessions", all(row["NCBI_accession"].startswith("SRR") for row in selected), selected)

    coverage = read_tsv(frozen / "model-coverage.tsv")
    modeled = [float(row["ModeledAbundance"]) for row in coverage]
    audit.add("Coverage", "six-rows", len(coverage) == 6, len(coverage))
    audit.add("Coverage", "minimum-gate", min(modeled) >= 0.5, min(modeled))
    audit.add("Coverage", "locked-minimum", near(min(modeled), 0.7479334), min(modeled))
    audit.add("Coverage", "locked-maximum", near(max(modeled), 0.9853743), max(modeled))
    audit.add("Coverage", "original-denominator", any(float(row["SpeciesResolvedAbundance"]) < 0.95 for row in coverage), coverage)

    medium = read_tsv(frozen / "raw/medium.tsv")
    positive = [row for row in medium if str(row["PositiveFlux"]).lower() == "true"]
    oxygen = [row for row in medium if row["Compound"] == "o2"]
    audit.add("Medium", "reaction-count", len(medium) == 171, len(medium))
    audit.add("Medium", "positive-count", len(positive) == 160, len(positive))
    audit.add("Medium", "trace-oxygen", len(oxygen) == 1 and near(oxygen[0]["flux"], 0.001), oxygen)
    smetana_medium = read_tsv(frozen / "raw/smetana-media.tsv")
    audit.add("Medium", "smetana-membership", len(smetana_medium) == 159, len(smetana_medium))
    audit.add("Medium", "smetana-anoxic", all(row["compound"] != "o2" for row in smetana_medium), smetana_medium[:5])
    medium_compatibility = read_tsv(frozen / "raw/smetana-medium-compatibility.tsv")
    audit.add("Medium", "pool-audit-size", len(medium_compatibility) == 159, len(medium_compatibility))
    audit.add(
        "Medium",
        "pool-matches",
        sum(truthy(row["PresentInSubcommunity"]) for row in medium_compatibility) == 120,
        sum(truthy(row["PresentInSubcommunity"]) for row in medium_compatibility),
    )

    model_audit = read_tsv(frozen / "raw/smetana-model-audit.tsv")
    audit.add("Model", "smetana-six", len(model_audit) == 6, len(model_audit))
    audit.add("Model", "unique-ids", len({row["ModelID"] for row in model_audit}) == 6, model_audit)
    audit.add("Model", "positive-structure", all(int(row["Reactions"]) > 0 and int(row["Metabolites"]) > 0 and int(row["Genes"]) > 0 for row in model_audit), model_audit)
    audit.add("Model", "smetana-id-convention", {row["SMETANAIDConvention"] for row in model_audit} == {"metabolite_compartment"}, model_audit)
    archives = sorted((frozen / "models/smetana").glob("*.xml.gz"))
    audit.add("Model", "bundled-six", len(archives) == 6, len(archives))
    for path in archives:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            prefix = handle.read(2000).lower()
        audit.add("Model", f"sbml-{path.name}", "<sbml" in prefix, path.stat().st_size)

    ledger = read_tsv(frozen / "run-ledger.tsv")
    audit.add("Run", "ten-steps", len(ledger) == 10, len(ledger))
    ledger_status = {row["Step"]: row["Status"] for row in ledger}
    audit.add("Run", "micom-eight-passed", sum(row["Status"] == "passed" for row in ledger) == 8, ledger)
    audit.add("Run", "smetana-global-not-estimable", ledger_status.get("SMETANA global MIP/MRO") == "not_estimable", ledger_status)
    audit.add("Run", "smetana-detailed-limited", ledger_status.get("SMETANA detailed SCS/MUS/MPS") == "passed_with_limitation", ledger_status)
    audit.add("Run", "thread-lock", json.loads((frozen / "raw/software-versions.json").read_text(encoding="utf-8")).get("threads") == 6, "threads=6")

    tradeoff = read_tsv(frozen / "tradeoff-summary.tsv")
    audit.add("Result", "tradeoff-grid", len(tradeoff) == 54 and {row["tradeoff"] for row in tradeoff} == {f"{x/10:.1f}" for x in range(1, 10)}, len(tradeoff))
    primary = read_tsv(frozen / "primary-growth-summary.tsv")
    audit.add("Result", "primary-six", len(primary) == 6, primary)
    sensitivity = read_tsv(frozen / "medium-sensitivity-summary.tsv")
    audit.add("Result", "medium-grid", len(sensitivity) == 18 and {row["MediumScale"] for row in sensitivity} == {"0.5", "1.0", "2.0"}, len(sensitivity))
    abundance = read_tsv(frozen / "abundance-sensitivity-summary.tsv")
    audit.add("Result", "abundance-six", len(abundance) == 6, abundance)

    exchanges = read_tsv(frozen / "raw/micom-primary-exchanges.tsv")
    direction_ok = all(
        (float(row["flux"]) > 0 and row["direction"] == "export")
        or (float(row["flux"]) < 0 and row["direction"] == "import")
        for row in exchanges
    )
    audit.add("Result", "exchange-rows", len(exchanges) > 0, len(exchanges))
    audit.add("Result", "exchange-directions", direction_ok, len(exchanges))
    audit.add("Result", "net-flux", len(read_tsv(frozen / "net-community-flux.tsv")) > 0, "non-empty")
    audit.add("Result", "crossfeeding-table", (frozen / "micom-crossfeeding-potential.tsv").stat().st_size > 0, "present")

    smetana_global = read_tsv(frozen / "smetana-global.tsv")
    smetana_detailed = read_tsv(frozen / "smetana-detailed.tsv")
    audit.add("Result", "smetana-global-row", len(smetana_global) == 1, smetana_global)
    audit.add("Result", "smetana-global-not-finite", len(smetana_global) == 1 and not finite(smetana_global[0]["mip"]) and not finite(smetana_global[0]["mro"]), smetana_global)
    audit.add("Result", "smetana-detailed", len(smetana_detailed) == 7853, len(smetana_detailed))
    audit.add("Result", "smetana-score-columns", {"scs", "mus", "mps", "smetana"} <= set(smetana_detailed[0]), sorted(smetana_detailed[0]) if smetana_detailed else [])
    score_counts = {
        column: sum(float(row[column]) > 0 for row in smetana_detailed)
        for column in ("scs", "mus", "mps", "smetana")
    }
    audit.add("Result", "smetana-component-counts", score_counts == {"scs": 0, "mus": 1330, "mps": 395, "smetana": 0}, score_counts)
    compatibility_audit = read_tsv(frozen / "smetana-compatibility-audit.tsv")
    audit.add("Result", "compatibility-six", len(compatibility_audit) == 6, len(compatibility_audit))
    audit.add("Result", "standalone-growth-control", all(truthy(row["StandalonePositive"]) for row in compatibility_audit), compatibility_audit)
    audit.add("Result", "interacting-growth-control", all(truthy(row["InteractingPositive"]) for row in compatibility_audit), compatibility_audit)
    audit.add("Result", "legacy-noninteracting-zero", not any(truthy(row["NoninteractingPositive"]) for row in compatibility_audit), compatibility_audit)
    compatibility_summary = json.loads((frozen / "raw/smetana-compatibility-summary.json").read_text(encoding="utf-8"))
    audit.add("Result", "compatibility-interpretation", compatibility_summary.get("interpretation") == "software_model_interface_limitation_not_biological_absence", compatibility_summary)
    concordance = read_tsv(frozen / "cross-method-concordance.tsv")
    audit.add("Result", "concordance-withheld", {row["EvidenceClass"] for row in concordance} == {"MICOM candidate; SMETANA unavailable"}, sorted({row["EvidenceClass"] for row in concordance}))

    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    audit.add("Summary", "article", metrics.get("article") == 61, metrics)
    audit.add("Summary", "counts", metrics.get("selected_samples") == 6 and metrics.get("independent_subjects") == 6 and metrics.get("smetana_subcommunity_size") == 6, metrics)
    audit.add("Summary", "completed-steps", metrics.get("completed_steps") == 8 and metrics.get("passed_with_limitation_steps") == 1 and metrics.get("not_estimable_steps") == 1, metrics)
    audit.add("Summary", "smetana-limitation", metrics.get("smetana_global_estimable") is False and metrics.get("smetana_cross_method_comparison_estimable") is False, metrics)

    audit_figures(figures, audit, FIGURES)
    plot_source = (root / "scripts/plot_article61_community.py").read_text(encoding="utf-8")
    audit.add("Figure", "english-only-source", re.search(r"[\u4e00-\u9fff]", plot_source) is None, "plot source contains no CJK text")
    audit_chapter_61(args.chapter.resolve(), audit)
    stage_reanalysis(root, frozen, args.stage_dir.resolve(), args.python.resolve(), audit)

    return finish(
        article=61,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "samples": len(selected),
            "figures": len(FIGURES),
            "smetana_rows": len(smetana_detailed),
            "exchange_rows": len(exchanges),
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
