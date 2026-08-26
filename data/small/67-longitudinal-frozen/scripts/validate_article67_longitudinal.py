#!/usr/bin/env python3
"""Fail-closed validation for Article 67 longitudinal analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from article42_44_validation_utils import Audit, audit_checksums, audit_figures, finish


FIGURES = (
    "67-sampling-design",
    "67-dysbiosis-trajectories",
    "67-lag-stability",
    "67-short-interval-shifts",
    "67-species-retention",
    "67-prevotella-trajectories",
    "67-mixed-model-effects",
)
ANCHOR_SHA = "b5ea6e10c6036945c97377705a3128c2d64164de0e84dc9d3874c2d242e5eed7"
EXPECTED_SOURCES = {
    "hmp2-metadata.csv": (
        9_074_342,
        "656b7bd97660ddb875548805e30bede31f2d1208293f7170d2d5755e33862ec9",
        "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/metadata/hmp2_metadata_2018-08-20.csv",
    ),
    "taxonomic_profiles.tsv.gz": (
        666_634,
        "5728531a6a1236371be6795b1da84ff5f6dd029035179d24d3c3be72d814e72c",
        "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles.tsv.gz",
    ),
    "lloyd-price-fig3-original.png": (
        349_963,
        ANCHOR_SHA,
        "https://media.springernature.com/lw1200/springer-static/image/art%3A10.1038%2Fs41586-019-1237-9/MediaObjects/41586_2019_1237_Fig3_HTML.png",
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


def near(value: object, expected: float, tolerance: float = 1e-9) -> bool:
    try:
        return math.isclose(float(value), expected, rel_tol=tolerance, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pixel_sha(path: Path) -> str:
    with Image.open(path) as image:
        normalized = image.convert("RGBA")
        payload = (
            normalized.size[0].to_bytes(8, "little")
            + normalized.size[1].to_bytes(8, "little")
            + normalized.tobytes()
        )
    return hashlib.sha256(payload).hexdigest()


def audit_chapter(chapter: Path, audit: Audit) -> None:
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
        "eval-true": re.search(r"execute:\s*\n\s+eval:\s+true", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "analysis-seed": "67001" in text,
        "plot-seed": "20260767" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("RAM", "CPU", "磁盘", "耗时")),
        "real-study": all(token in text for token in ("Lloyd-Price", "10.1038/s41586-019-1237-9", "PMC6650278", "HMP2")),
        "sample-contract": all(token in text for token in ("1,638", "1,595", "1,523", "107", "43")),
        "eligibility-contract": all(token in text for token in ("1,000,000", "≥4", "filtered reads")),
        "model-contract": all(token in text for token in ("random intercept", "random slope", "Kenward--Roger", "Satterthwaite", "0.0514", "0.0322", "0.0424")),
        "uncertainty-contract": all(token in text for token in ("2,000", "subject bootstrap", "主体级 bootstrap")),
        "shift-boundary": all(token in text for token in ("0.54", "sensitivity analysis", "1--3 周", "lag-dependent")),
        "colonization-boundary": all(token in text for token in ("Detection stability", "permanent colonization", "strain-resolved")),
        "composition-contract": all(token in text for token in ("63", "315", "CLR", "Benjamini--Hochberg", "10^{-6}")),
        "methods-template": "**Longitudinal metagenomic analysis.**" in text,
        "results-template": "**Results template.**" in text,
        "frozen-input": "data/small/67-longitudinal-frozen" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "anchor": "../figures/67-lloyd-price-fig3-original.png" in text,
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-meta-prose": not any(token in text for token in ("本篇可独立", "本文可独立", "全系列约定", "接口只学一次", "作者代码通常长这样", "（即本文）", "无头服务器")),
    }
    for section in sections:
        checks[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def stage_figures(
    root: Path,
    frozen: Path,
    stage: Path,
    python: Path,
    final: Path,
    audit: Audit,
) -> None:
    if stage.exists():
        shutil.rmtree(stage)
    figures = stage / "figures"
    environment = os.environ.copy()
    environment.update(
        {
            "MPLCONFIGDIR": str(stage / "matplotlib"),
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    result = subprocess.run(
        [
            str(python),
            str(frozen / "scripts/plot_article67_longitudinal.py"),
            "--input-dir",
            str(frozen),
            "--figure-dir",
            str(figures),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        timeout=240,
    )
    audit.add("Reanalysis", "plot-script-exit", result.returncode == 0, result.stdout + result.stderr)
    for stem in FIGURES:
        staged = figures / f"{stem}.png"
        published = final / f"{stem}.png"
        status = staged.is_file() and published.is_file() and pixel_sha(staged) == pixel_sha(published)
        audit.add(
            "Reanalysis",
            f"{stem}-pixel-identical",
            status,
            pixel_sha(staged) if staged.is_file() else "MISSING",
        )
    staged_anchor = figures / "67-lloyd-price-fig3-original.png"
    final_anchor = final / "67-lloyd-price-fig3-original.png"
    audit.add(
        "Reanalysis",
        "anchor-byte-identical",
        staged_anchor.is_file()
        and final_anchor.is_file()
        and sha256(staged_anchor) == sha256(final_anchor),
        sha256(staged_anchor) if staged_anchor.is_file() else "MISSING",
    )


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_checksums(frozen, audit)

    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    expected_metrics = {
        "article": 67,
        "bootstrap_iterations": 2000,
        "published_profiles": 1638,
        "published_subjects": 130,
        "primary_biospecimens": 1595,
        "excluded_duplicate_profiles": 43,
        "eligible_samples": 1523,
        "eligible_subjects": 107,
        "terminal_species_rows": 572,
        "ever_observed_species": 568,
        "selected_species": 63,
        "dysbiosis_reference_samples": 207,
        "dysbiosis_reference_subjects": 26,
        "short_intervals": 1005,
        "shift_events": 123,
        "pcopri_selected_subjects": 12,
        "read_gate": 1_000_000,
        "minimum_visits": 4,
        "seed": 67_001,
        "plot_seed": 20_260_767,
    }
    for key, expected in expected_metrics.items():
        audit.add("Metric", key, metrics.get(key) == expected, metrics.get(key))
    audit.add("Metric", "dysbiosis-threshold", near(metrics.get("dysbiosis_threshold"), 0.8323583208493152), metrics.get("dysbiosis_threshold"))
    audit.add("Metric", "shift-threshold", near(metrics.get("shift_threshold"), 0.54), metrics.get("shift_threshold"))

    manifest = json.loads((frozen / "source-manifest.json").read_text(encoding="utf-8"))
    audit.add("Source", "doi", manifest.get("doi") == "10.1038/s41586-019-1237-9", manifest.get("doi"))
    audit.add("Source", "pmcid", manifest.get("pmcid") == "PMC6650278", manifest.get("pmcid"))
    audit.add("Source", "profile-release", manifest.get("profile_release") == "2018-05-04", manifest.get("profile_release"))
    audit.add("Source", "metadata-release", manifest.get("metadata_release") == "2018-08-20", manifest.get("metadata_release"))
    for name, (size, digest, url) in EXPECTED_SOURCES.items():
        record = manifest.get("resources", {}).get(name, {})
        audit.add("Source", f"{name}-bytes", record.get("bytes") == size, record.get("bytes"))
        audit.add("Source", f"{name}-sha256", record.get("sha256") == digest, record.get("sha256"))
        audit.add("Source", f"{name}-official-url", record.get("url") == url, record.get("url"))

    selection = pd.read_csv(frozen / "profile-selection-ledger.tsv", sep="\t")
    audit.add("Selection", "rows", len(selection) == 1638, len(selection))
    audit.add("Selection", "biospecimens", selection["site_sub_coll"].nunique() == 1595, selection["site_sub_coll"].nunique())
    primary = selection["PrimaryProfile"].astype(bool)
    audit.add("Selection", "one-primary-each", primary.sum() == 1595 and selection.loc[primary, "site_sub_coll"].is_unique, int(primary.sum()))
    audit.add("Selection", "duplicates-retained-in-ledger", (~primary).sum() == 43, int((~primary).sum()))
    audit.add("Selection", "rank-contract", selection.loc[primary, "SelectionRank"].eq(1).all(), selection.loc[primary, "SelectionRank"].value_counts().to_dict())

    attrition = pd.read_csv(frozen / "sample-attrition.tsv", sep="\t")
    audit.add("Attrition", "profiles", attrition["Profiles"].tolist() == [1638, 1595, 1553, 1523], attrition["Profiles"].tolist())
    audit.add("Attrition", "subjects", attrition["Subjects"].tolist() == [130, 130, 130, 107], attrition["Subjects"].tolist())
    audit.add("Attrition", "eligible-groups", tuple(attrition.iloc[-1][["Control", "CD", "UC"]].astype(int)) == (411, 699, 413), attrition.iloc[-1].to_dict())

    normalization = pd.read_csv(frozen / "normalization-audit.tsv", sep="\t")
    audit.add("Composition", "normalization-rows", len(normalization) == 1523, len(normalization))
    audit.add("Composition", "positive-source-sums", normalization["TerminalSpeciesSumBeforeRenormalization"].gt(0).all(), float(normalization["TerminalSpeciesSumBeforeRenormalization"].min()))
    audit.add("Composition", "closed-sums", np.allclose(normalization["TerminalSpeciesSumAfterRenormalization"], 1.0), float(np.abs(normalization["TerminalSpeciesSumAfterRenormalization"] - 1).max()))
    feature = pd.read_csv(frozen / "species-feature-audit.tsv", sep="\t")
    audit.add("Composition", "features", len(feature) == 568, len(feature))
    audit.add("Composition", "feature-ids-unique", feature["FeatureID"].is_unique, int(feature["FeatureID"].duplicated().sum()))
    audit.add("Composition", "screened-features", feature["SelectedForMixedModel"].astype(bool).sum() == 63, int(feature["SelectedForMixedModel"].astype(bool).sum()))
    audit.add("Composition", "pcopri-one-row", feature["Species"].eq("s__Prevotella_copri").sum() == 1, int(feature["Species"].eq("s__Prevotella_copri").sum()))

    sample = pd.read_csv(frozen / "sample-ledger.tsv", sep="\t")
    audit.add("Sample", "rows-subjects", len(sample) == 1523 and sample["SubjectID"].nunique() == 107, {"rows": len(sample), "subjects": sample["SubjectID"].nunique()})
    observed_groups = sample.groupby("Diagnosis").agg(Samples=("SampleID", "size"), Subjects=("SubjectID", "nunique")).to_dict("index")
    expected_groups = {"Control": {"Samples": 411, "Subjects": 27}, "CD": {"Samples": 699, "Subjects": 50}, "UC": {"Samples": 413, "Subjects": 30}}
    audit.add("Sample", "diagnosis-groups", observed_groups == expected_groups, observed_groups)
    audit.add("Sample", "read-gate", sample["FilteredReads"].ge(1_000_000).all(), float(sample["FilteredReads"].min()))
    audit.add("Sample", "minimum-visits", sample.groupby("SubjectID").size().min() >= 4, int(sample.groupby("SubjectID").size().min()))

    dysbiosis = pd.read_csv(frozen / "dysbiosis-summary.tsv", sep="\t").set_index("Diagnosis")
    expected_fractions = {"Control": 41 / 411 * 100, "CD": 178 / 699 * 100, "UC": 53 / 413 * 100}
    for diagnosis, expected in expected_fractions.items():
        value = dysbiosis.loc[diagnosis, "DysbioticSamplePercent"]
        audit.add("Dysbiosis", diagnosis, near(value, expected), value)

    pairs = pd.read_csv(frozen / "within-subject-pairs.tsv.gz", sep="\t")
    audit.add("Distance", "pair-rows", len(pairs) == 11275, len(pairs))
    audit.add("Distance", "positive-lags", pairs["LagWeeks"].gt(0).all(), float(pairs["LagWeeks"].min()))
    audit.add("Distance", "bounded-bray", pairs["BrayCurtis"].between(0, 1).all(), [float(pairs["BrayCurtis"].min()), float(pairs["BrayCurtis"].max())])
    lag = pd.read_csv(frozen / "lag-summary.tsv", sep="\t")
    audit.add("Distance", "lag-grid", len(lag) == 18 and lag["Diagnosis"].nunique() == 3 and lag["LagBin"].nunique() == 6, lag.shape)

    consecutive = pd.read_csv(frozen / "consecutive-intervals.tsv.gz", sep="\t")
    audit.add("Shift", "consecutive-rows", len(consecutive) == 1416, len(consecutive))
    audit.add("Shift", "short-rows", consecutive["ShortInterval"].astype(bool).sum() == 1005, int(consecutive["ShortInterval"].astype(bool).sum()))
    audit.add("Shift", "events", consecutive["Shift054"].astype(bool).sum() == 123, int(consecutive["Shift054"].astype(bool).sum()))
    shift = pd.read_csv(frozen / "shift-summary.tsv", sep="\t").set_index("Diagnosis")
    expected_events = {"Control": (281, 34), "CD": (454, 63), "UC": (270, 26)}
    for diagnosis, expected in expected_events.items():
        observed = tuple(shift.loc[diagnosis, ["ShortIntervals", "ShiftEvents"]].astype(int))
        audit.add("Shift", diagnosis, observed == expected, observed)

    retention = pd.read_csv(frozen / "retention-summary.tsv", sep="\t")
    audit.add("Retention", "grid", len(retention) == 9 and retention["Subjects"].min() == 27, retention.shape)
    audit.add("Retention", "bounded", retention[["SubjectMedianRetention", "CILower", "CIUpper"]].apply(lambda x: x.between(0, 1)).all().all(), "bounded")
    pcopri = pd.read_csv(frozen / "prevotella-selected-trajectories.tsv", sep="\t")
    selected_counts = pcopri.groupby("Diagnosis")["SubjectID"].nunique().to_dict()
    audit.add("Trajectory", "pcopri-four-per-group", selected_counts == {"CD": 4, "Control": 4, "UC": 4}, selected_counts)

    fixed = pd.read_csv(frozen / "primary-fixed-effects.tsv", sep="\t").set_index("Term")
    expected_fixed = {
        "DiagnosisCD": (0.0514220832237571, 0.00304603064567568),
        "DiagnosisUC": (0.0184522377672547, 0.329476617026153),
        "AntibioticsYes": (0.0322190173016928, 1.97094896054332e-09),
        "DiagnosisCD:WeekYearCentered": (0.0423905165522957, 0.0228204520755285),
    }
    for term, (estimate, pvalue) in expected_fixed.items():
        observed = fixed.loc[term]
        audit.add("Model", f"{term}-estimate", near(observed["Estimate"], estimate), observed["Estimate"])
        audit.add("Model", f"{term}-p", near(observed["PValue"], pvalue), observed["PValue"])
    model_audit = pd.read_csv(frozen / "primary-model-audit.tsv", sep="\t")
    audit.add("Model", "three-primary-fits", len(model_audit) == 3, len(model_audit))
    audit.add("Model", "non-singular", (~model_audit["Singular"].astype(bool)).all(), model_audit["Singular"].tolist())
    audit.add("Model", "converged", model_audit["ConvergenceMessage"].eq("none").all() and model_audit["MaxAbsGradient"].max() < 1e-3, model_audit.to_dict("records"))
    comparison = pd.read_csv(frozen / "random-slope-model-comparison.tsv", sep="\t")
    audit.add("Model", "random-slope-lrt", near(comparison.iloc[-1]["Pr(>Chisq)"], 6.6768371248067e-12), comparison.iloc[-1]["Pr(>Chisq)"])

    secondary = pd.read_csv(frozen / "species-mixed-model-results.tsv", sep="\t")
    secondary_audit = pd.read_csv(frozen / "species-model-audit.tsv", sep="\t")
    audit.add("Model", "secondary-tests", len(secondary) == 315 and secondary["FeatureID"].nunique() == 63 and secondary["Term"].nunique() == 5, secondary.shape)
    audit.add("Model", "secondary-global-bh", int((secondary["QValueGlobal"] <= 0.05).sum()) == 36, int((secondary["QValueGlobal"] <= 0.05).sum()))
    audit.add("Model", "secondary-convergence", len(secondary_audit) == 63 and secondary_audit["Success"].astype(bool).all() and (~secondary_audit["Singular"].astype(bool)).all() and secondary_audit["MaxAbsGradient"].max() < 1e-3, secondary_audit.to_dict("records")[:3])

    r_versions = pd.read_csv(frozen / "software-versions-r.tsv", sep="\t").set_index("Package")["Version"].astype(str).to_dict()
    expected_versions = {"R": "4.4.1", "Matrix": "1.6.5", "lme4": "1.1.35.4", "lmerTest": "3.1.3", "pbkrtest": "0.5.2"}
    audit.add("Environment", "r-versions", r_versions == expected_versions, r_versions)

    anchor = frozen / "lloyd-price-fig3-original.png"
    audit.add("Anchor", "source-checksum", anchor.is_file() and sha256(anchor) == ANCHOR_SHA, sha256(anchor) if anchor.is_file() else "MISSING")
    audit_chapter(args.chapter.resolve(), audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    stage_figures(root, frozen, args.stage_dir.resolve(), args.python.resolve(), args.figure_dir.resolve(), audit)
    return_code = finish(
        article=67,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "study": "Lloyd-Price et al. Nature 2019",
            "eligible_samples": len(sample),
            "eligible_subjects": int(sample["SubjectID"].nunique()),
            "selected_species": int(feature["SelectedForMixedModel"].astype(bool).sum()),
            "figure_sets": len(FIGURES),
        },
    )
    raise SystemExit(return_code)


if __name__ == "__main__":
    main()
