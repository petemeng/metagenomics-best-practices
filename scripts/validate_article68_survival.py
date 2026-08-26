#!/usr/bin/env python3
"""Offline acceptance tests for the frozen Article 68 evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


FIGURES = (
    "68-cohort-design",
    "68-cox-forest",
    "68-assumption-and-spline",
    "68-cutoff-leakage",
    "68-time-dependent-roc",
    "68-crossvalidated-calibration",
    "68-influence-sensitivity",
)
EXPECTED_SOURCES = {
    "spencer-human-wgs.xlsx": (
        22_718_390,
        "4fa937513f33a9fe2e1127554caf89c4287698a422123213428c73d0f2bb0968",
    ),
    "spencer-data-dictionary.xlsx": (
        101_198,
        "dc200443076f6fc252eed6b9239446712cfbf65aded581896a7498c0827de15c",
    ),
    "spencer-paper.pdf": (
        3_071_441,
        "414c9350d6fef5c1916b464c43b271992d14d6c814ae9817dfa8dd51b642c707",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--qa-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pixel_sha(path: Path) -> str:
    with Image.open(path) as image:
        normalized = image.convert("RGBA")
        digest = hashlib.sha256()
        digest.update(f"{normalized.width}x{normalized.height}".encode())
        digest.update(normalized.tobytes())
        return digest.hexdigest()


def near(value: object, expected: float, tolerance: float = 1e-9) -> bool:
    try:
        return bool(np.isclose(float(value), expected, rtol=tolerance, atol=tolerance))
    except (TypeError, ValueError):
        return False


@dataclass
class Check:
    category: str
    check: str
    status: bool
    detail: str


class Audit:
    def __init__(self) -> None:
        self.rows: list[Check] = []

    def add(self, category: str, check: str, status: bool, detail: object = "") -> None:
        self.rows.append(Check(category, check, bool(status), str(detail)))

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "Category": row.category,
                    "Check": row.check,
                    "Status": "PASS" if row.status else "FAIL",
                    "Detail": row.detail,
                }
                for row in self.rows
            ]
        )


def audit_checksums(frozen: Path, audit: Audit) -> None:
    checksum_file = frozen / "file-checksums.sha256"
    audit.add("Bundle", "checksum-file", checksum_file.is_file(), checksum_file)
    if not checksum_file.is_file():
        return
    lines = [line for line in checksum_file.read_text(encoding="utf-8").splitlines() if line]
    audit.add("Bundle", "checksum-count", len(lines) == 48, len(lines))
    for line in lines:
        digest, relative = line.split("  ", 1)
        path = frozen / relative
        audit.add("Checksum", relative, path.is_file() and sha256(path) == digest, digest)


def audit_sources(frozen: Path, audit: Audit) -> None:
    manifest = json.loads((frozen / "source-manifest.json").read_text(encoding="utf-8"))
    audit.add("Source", "article", manifest.get("article") == 68, manifest.get("article"))
    audit.add("Source", "doi", manifest.get("doi") == "10.1126/science.aaz7015", manifest.get("doi"))
    audit.add("Source", "pmcid", manifest.get("pmcid") == "PMC8970537", manifest.get("pmcid"))
    audit.add("Source", "bioproject", manifest.get("bioproject") == "PRJNA770295", manifest.get("bioproject"))
    audit.add(
        "Source",
        "repository-commit",
        manifest.get("repository_commit") == "a95b1a020b890dbe93a960dec946e197b63b7d15",
        manifest.get("repository_commit"),
    )
    for name, (size, digest) in EXPECTED_SOURCES.items():
        record = manifest.get("resources", {}).get(name, {})
        audit.add("Source", f"{name}-bytes", record.get("bytes") == size, record.get("bytes"))
        audit.add("Source", f"{name}-sha256", record.get("sha256") == digest, record.get("sha256"))
        audit.add("Source", f"{name}-https", str(record.get("url", "")).startswith("https://"), record.get("url"))
    dictionary = frozen / "spencer-data-dictionary.xlsx"
    audit.add("Source", "frozen-data-dictionary", dictionary.stat().st_size == 101_198 and sha256(dictionary) == EXPECTED_SOURCES["spencer-data-dictionary.xlsx"][1], sha256(dictionary))


def audit_cohort(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    expected = {
        "article": 68,
        "public_wgs_samples": 167,
        "complete_pfs_samples": 110,
        "pfs_events": 61,
        "pfs_censored": 49,
        "reported_lkt_features": 225,
        "faecalibacterium_lkt_rows": 3,
        "faecalibacterium_detected_samples": 109,
        "quality_sensitivity_samples": 109,
        "quality_sensitivity_events": 60,
        "analysis_seed": 68_001,
        "plot_seed": 20_260_768,
    }
    for key, value in expected.items():
        audit.add("Metric", key, metrics.get(key) == value, metrics.get(key))
    audit.add("Metric", "pseudocount", near(metrics.get("pseudocount_ppm"), 25.0), metrics.get("pseudocount_ppm"))

    cohort = pd.read_csv(frozen / "survival-cohort.tsv", sep="\t")
    audit.add("Cohort", "rows", len(cohort) == 110, len(cohort))
    audit.add("Cohort", "unique-patient-samples", cohort["SampleID"].nunique() == 110, cohort["SampleID"].nunique())
    audit.add("Cohort", "event-count", int(cohort["Event"].sum()) == 61 and set(cohort["Event"]) == {0, 1}, cohort["Event"].value_counts().to_dict())
    audit.add("Cohort", "positive-times", cohort["PFS_days"].gt(0).all(), cohort["PFS_days"].min())
    audit.add("Cohort", "time-conversion", np.allclose(cohort["PFS_months"], cohort["PFS_days"] / 30.4375), float(np.max(np.abs(cohort["PFS_months"] - cohort["PFS_days"] / 30.4375))))
    audit.add("Cohort", "primary-covariates-complete", not cohort[["BMI", "PrimarySubtype", "AdvancedSubstage", "LDH"]].isna().any().any(), cohort[["BMI", "PrimarySubtype", "AdvancedSubstage", "LDH"]].isna().sum().to_dict())
    audit.add("Cohort", "subtype-levels", set(cohort["PrimarySubtype"]) == {"Cutaneous_or_unknown", "Mucosal_or_acral"}, cohort["PrimarySubtype"].value_counts().to_dict())
    audit.add("Cohort", "substage-levels", set(cohort["AdvancedSubstage"]) == {"Stage_M1C", "Stage_M1D"}, cohort["AdvancedSubstage"].value_counts().to_dict())
    audit.add("Cohort", "ldh-levels", set(cohort["LDH"]) == {"No", "Yes"}, cohort["LDH"].value_counts().to_dict())
    quality = cohort["QualitySensitivityPass"].astype(str).str.lower().eq("true")
    audit.add("Cohort", "quality-sensitivity-one-exclusion", int((~quality).sum()) == 1 and int(cohort.loc[quality, "Event"].sum()) == 60, {"pass": int(quality.sum()), "events": int(cohort.loc[quality, "Event"].sum())})
    audit.add("Cohort", "faec-prevalence", int(cohort["FaecalibacteriumPPM"].gt(0).sum()) == 109, int(cohort["FaecalibacteriumPPM"].gt(0).sum()))
    audit.add("Cohort", "faec-transform", np.allclose(cohort["FaecalibacteriumLog2"], np.log2(cohort["FaecalibacteriumPPM"] + 25)), float(np.max(np.abs(cohort["FaecalibacteriumLog2"] - np.log2(cohort["FaecalibacteriumPPM"] + 25)))))

    attrition = pd.read_csv(frozen / "cohort-attrition.tsv", sep="\t")
    audit.add("Attrition", "samples", attrition["Samples"].tolist() == [167, 110, 110, 109], attrition["Samples"].tolist())
    audit.add("Attrition", "events", attrition["Events"].tolist() == [61, 61, 61, 60], attrition["Events"].tolist())

    matrix = pd.read_csv(frozen / "lkt-survival-ppm.tsv.gz", sep="\t").set_index("LKTFeature")
    feature = pd.read_csv(frozen / "lkt-feature-map.tsv", sep="\t").set_index("LKTFeature")
    audit.add("Composition", "matrix-shape", matrix.shape == (225, 110), matrix.shape)
    audit.add("Composition", "feature-map-shape", len(feature) == 225 and feature.index.is_unique, {"rows": len(feature), "unique": feature.index.is_unique})
    faec_rows = feature.index[feature["Genus"].eq("g__Faecalibacterium")]
    audit.add("Composition", "three-faec-rows", len(faec_rows) == 3, faec_rows.tolist())
    aggregate = matrix.loc[faec_rows].sum(axis=0)
    aggregate = aggregate.loc[cohort["SampleID"]].to_numpy(float)
    audit.add("Composition", "faec-aggregation", np.allclose(aggregate, cohort["FaecalibacteriumPPM"]), float(np.max(np.abs(aggregate - cohort["FaecalibacteriumPPM"]))))
    totals = matrix.sum(axis=0).loc[cohort["SampleID"]].to_numpy(float)
    audit.add("Composition", "reported-ppm-sums", np.allclose(totals, cohort["ReportedLKTPPM"]), float(np.max(np.abs(totals - cohort["ReportedLKTPPM"]))))


def audit_models(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "model-metrics.json").read_text(encoding="utf-8"))
    exact = {
        "article": 68,
        "analysis_seed": 68_001,
        "bootstrap_iterations": 1000,
        "cv_repeats": 20,
        "cv_folds": 5,
        "train_samples": 77,
        "test_samples": 33,
        "train_events": 42,
        "test_events": 19,
    }
    for key, value in exact.items():
        audit.add("Model metric", key, metrics.get(key) == value, metrics.get(key))
    expected_float = {
        "unadjusted_faecalibacterium_hr_per_doubling": 0.902673869660564,
        "unadjusted_faecalibacterium_p": 0.0194706089190765,
        "adjusted_faecalibacterium_hr_per_doubling": 0.932942365742919,
        "adjusted_faecalibacterium_p": 0.0959456260571256,
        "quality_sensitivity_faecalibacterium_hr_per_doubling": 0.916467448717361,
        "quality_sensitivity_faecalibacterium_p": 0.0587956005825838,
        "global_ph_p": 0.767439649263038,
        "nonlinearity_p": 0.0622156483326983,
        "training_cutoff_log2": 14.8838842416909,
        "training_cutoff_ppm": 30209.0,
        "leaky_full_cutoff_log2": 14.8838842416909,
        "heldout_cutoff_hr": 0.630383181627237,
        "heldout_cutoff_p": 0.312338763242907,
        "cv_cindex_clinical": 0.63634422378817,
        "cv_cindex_microbiome": 0.664104206705104,
        "cv_auc365_clinical": 0.629307099758868,
        "cv_auc365_microbiome": 0.676617003633348,
    }
    for key, value in expected_float.items():
        audit.add("Model metric", key, near(metrics.get(key), value, 1e-8), metrics.get(key))

    estimates = pd.read_csv(frozen / "cox-model-estimates.tsv", sep="\t")
    faec = estimates.loc[estimates["Term"].eq("FaecalibacteriumLog2")]
    audit.add("Cox", "four-predeclared-faec-models", len(faec) == 4 and set(faec["Model"]) == {"Unadjusted", "Adjusted primary", "Sequencing-QC sensitivity", "Anti-PD1 sensitivity"}, faec["Model"].tolist())
    audit.add("Cox", "hazard-ratios-positive", estimates["HazardRatio"].gt(0).all() and estimates["CILower"].gt(0).all(), float(estimates["CILower"].min()))
    audit.add("Cox", "ci-order", (estimates["CILower"] <= estimates["HazardRatio"]).all() and (estimates["HazardRatio"] <= estimates["CIUpper"]).all(), "all terms")

    ph = pd.read_csv(frozen / "proportional-hazards-tests.tsv", sep="\t")
    audit.add("Assumption", "ph-terms", set(ph["Term"]) == {"FaecalibacteriumLog2", "PrimarySubtype", "AdvancedSubstage", "LDH", "BMIz", "GLOBAL"}, ph["Term"].tolist())
    audit.add("Assumption", "global-ph", near(ph.loc[ph["Term"].eq("GLOBAL"), "PValue"].iloc[0], 0.767439649263038), ph.loc[ph["Term"].eq("GLOBAL"), "PValue"].iloc[0])
    residual = pd.read_csv(frozen / "faecalibacterium-schoenfeld.tsv", sep="\t")
    audit.add("Assumption", "event-residual-count", len(residual) == 61 and residual["EventTimeDays"].gt(0).all(), len(residual))
    nonlinear = pd.read_csv(frozen / "nonlinearity-test.tsv", sep="\t")
    audit.add("Assumption", "nonlinearity-test", len(nonlinear) == 1 and near(nonlinear["PValue"].iloc[0], 0.0622156483326983), nonlinear.to_dict("records"))
    spline = pd.read_csv(frozen / "spline-effect-curve.tsv", sep="\t")
    audit.add("Assumption", "spline-grid", len(spline) == 121 and spline["HazardRatio"].gt(0).all(), len(spline))
    influence = pd.read_csv(frozen / "leave-one-out-influence.tsv", sep="\t")
    audit.add("Influence", "one-refit-per-patient", len(influence) == 110 and influence["OmittedSampleID"].nunique() == 110, len(influence))
    audit.add("Influence", "stable-direction", influence["HazardRatio"].lt(1).all(), [float(influence["HazardRatio"].min()), float(influence["HazardRatio"].max())])


def audit_leakage_and_prediction(frozen: Path, audit: Audit) -> None:
    split = pd.read_csv(frozen / "cutoff-split-ledger.tsv", sep="\t")
    training_ids = set(split.loc[split["Split"].eq("Training"), "SampleID"])
    test_ids = set(split.loc[split["Split"].eq("Held-out test"), "SampleID"])
    audit.add("Cutoff", "split-counts", len(training_ids) == 77 and len(test_ids) == 33, {"train": len(training_ids), "test": len(test_ids)})
    audit.add("Cutoff", "split-disjoint", training_ids.isdisjoint(test_ids) and len(training_ids | test_ids) == 110, len(training_ids & test_ids))
    audit.add("Cutoff", "event-balance", int(split.loc[split["Split"].eq("Training"), "Event"].sum()) == 42 and int(split.loc[split["Split"].eq("Held-out test"), "Event"].sum()) == 19, split.groupby("Split")["Event"].sum().to_dict())

    search_training = pd.read_csv(frozen / "cutoff-search-training.tsv", sep="\t")
    search_full = pd.read_csv(frozen / "cutoff-search-full-leaky.tsv", sep="\t")
    selected_training = search_training["Selected"].astype(str).str.lower().eq("true")
    selected_full = search_full["Selected"].astype(str).str.lower().eq("true")
    audit.add("Cutoff", "one-training-selection", selected_training.sum() == 1, int(selected_training.sum()))
    audit.add("Cutoff", "one-leaky-selection", selected_full.sum() == 1, int(selected_full.sum()))
    train_cutoff = search_training.loc[selected_training, "CutoffLog2"].iloc[0]
    full_cutoff = search_full.loc[selected_full, "CutoffLog2"].iloc[0]
    audit.add("Cutoff", "training-optimum", near(train_cutoff, search_training.loc[search_training["PValue"].idxmin(), "CutoffLog2"]), train_cutoff)
    audit.add("Cutoff", "full-optimum", near(full_cutoff, search_full.loc[search_full["PValue"].idxmin(), "CutoffLog2"]), full_cutoff)
    cutoff_eval = pd.read_csv(frozen / "cutoff-evaluation.tsv", sep="\t")
    heldout = cutoff_eval.loc[cutoff_eval["EvaluationData"].eq("Held-out test")].iloc[0]
    leaky = cutoff_eval.loc[cutoff_eval["EvaluationData"].eq("Full data")].iloc[0]
    audit.add("Cutoff", "heldout-uses-training-cutoff", heldout["CutoffSource"] == "Training" and near(heldout["CutoffLog2"], train_cutoff), heldout.to_dict())
    audit.add("Cutoff", "leaky-labeled", leaky["CutoffSource"] == "Full data (leaky)", leaky["CutoffSource"])
    audit.add("Cutoff", "optimism-visible", heldout["LogRankPValue"] > 0.30 and leaky["LogRankPValue"] < 0.003, {"heldout": heldout["LogRankPValue"], "leaky": leaky["LogRankPValue"]})
    audit.add("Cutoff", "group-count-conservation", (cutoff_eval["LowN"] + cutoff_eval["HighN"]).tolist() == [77, 33, 110], (cutoff_eval["LowN"] + cutoff_eval["HighN"]).tolist())

    cv_long = pd.read_csv(frozen / "repeated-cv-predictions-long.tsv", sep="\t")
    audit.add("Prediction", "cv-long-rows", len(cv_long) == 2200, len(cv_long))
    audit.add("Prediction", "twenty-oof-per-patient", cv_long.groupby("SampleID").size().eq(20).all() and cv_long["SampleID"].nunique() == 110, cv_long.groupby("SampleID").size().value_counts().to_dict())
    audit.add("Prediction", "one-fold-per-repeat", cv_long.groupby(["Repeat", "SampleID"]).size().eq(1).all(), int(cv_long.groupby(["Repeat", "SampleID"]).size().max()))
    audit.add("Prediction", "five-folds-each-repeat", cv_long.groupby("Repeat")["Fold"].nunique().eq(5).all() and cv_long["Repeat"].nunique() == 20, cv_long.groupby("Repeat")["Fold"].nunique().to_dict())
    audit.add("Prediction", "finite-oof-predictions", np.isfinite(cv_long.filter(regex="LP|Risk").to_numpy(float)).all(), cv_long.filter(regex="LP|Risk").isna().sum().sum())
    risk_columns = [column for column in cv_long if "Risk" in column]
    audit.add("Prediction", "risks-in-unit-interval", ((cv_long[risk_columns] >= 0) & (cv_long[risk_columns] <= 1)).all().all(), {"min": float(cv_long[risk_columns].min().min()), "max": float(cv_long[risk_columns].max().max())})

    fold = pd.read_csv(frozen / "repeated-cv-fold-audit.tsv", sep="\t")
    audit.add("Prediction", "fold-audit-rows", len(fold) == 100, len(fold))
    for _, row in fold.iterrows():
        train = set(str(row["TrainingSampleHash"]).split(";"))
        test = set(str(row["TestSampleHash"]).split(";"))
        status = train.isdisjoint(test) and len(train | test) == 110 and len(train) == int(row["TrainingN"]) and len(test) == int(row["TestN"])
        audit.add("Fold leakage", f"repeat-{int(row['Repeat']):02d}-fold-{int(row['Fold'])}", status, {"train": len(train), "test": len(test), "overlap": len(train & test)})

    cv = pd.read_csv(frozen / "repeated-cv-predictions.tsv", sep="\t")
    audit.add("Prediction", "one-averaged-row-per-patient", len(cv) == 110 and cv["SampleID"].nunique() == 110, len(cv))
    perf = pd.read_csv(frozen / "cv-performance-summary.tsv", sep="\t")
    audit.add("Prediction", "performance-rows", len(perf) == 12, len(perf))
    audit.add("Prediction", "bootstrap-complete", perf["ValidBootstrap"].eq(1000).all(), perf["ValidBootstrap"].value_counts().to_dict())
    bootstrap = pd.read_csv(frozen / "performance-bootstrap.tsv", sep="\t")
    audit.add("Prediction", "bootstrap-rows", len(bootstrap) == 1000 and bootstrap["Iteration"].is_unique, len(bootstrap))
    audit.add("Prediction", "bootstrap-finite", np.isfinite(bootstrap.drop(columns="Iteration").to_numpy(float)).all(), int(bootstrap.isna().sum().sum()))
    calibration = pd.read_csv(frozen / "calibration-365d.tsv", sep="\t")
    audit.add("Prediction", "calibration-ten-bins", len(calibration) == 10 and set(calibration["RiskQuintile"]) == {1, 2, 3, 4, 5}, len(calibration))
    audit.add("Prediction", "calibration-counts", calibration.groupby("Model")["N"].sum().eq(110).all(), calibration.groupby("Model")["N"].sum().to_dict())
    roc = pd.read_csv(frozen / "time-dependent-roc-365d.tsv", sep="\t")
    audit.add("Prediction", "roc-models", set(roc["Model"]) == {"Clinical", "Clinical + Faecalibacterium"}, roc["Model"].value_counts().to_dict())
    audit.add("Prediction", "roc-unit-square", ((roc[["FalsePositiveRate", "TruePositiveRate"]] >= 0) & (roc[["FalsePositiveRate", "TruePositiveRate"]] <= 1)).all().all(), "all points")


def audit_chapter(root: Path, audit: Audit) -> None:
    chapter = root / "chapters" / "68-survival-analysis.qmd"
    audit.add("Chapter", "exists", chapter.is_file(), chapter)
    if not chapter.is_file():
        return
    text = chapter.read_text(encoding="utf-8")
    audit.add("Chapter", "published", "draft: false" in text, "draft: false")
    audit.add("Chapter", "eval-true", "eval: true" in text, "eval: true")
    audit.add("Chapter", "freeze-auto", "freeze: auto" in text, "freeze: auto")
    required = (
        "对应论文里的哪张图",
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
    for heading in required:
        audit.add("Chapter structure", heading, heading in text, heading)
    audit.add("Chapter", "extra-follow-up-data-stated", "PFS/OS 随访时间" in text and "结局指示变量" in text, "follow-up time + event indicator")
    audit.add("Chapter", "cutoff-leakage-explicit", "cutoff 数据泄漏" in text and "训练集" in text and "测试集" in text, "training-derived cutoff")
    audit.add("Chapter", "random-seeds", "set.seed(68001)" in text and "20260768" in text, "68001 / 20260768")
    audit.add("Chapter", "inline-theme", "theme_pub <- function" in text and "save_pub <- function" in text, "inline functions")
    audit.add("Chapter", "no-source-theme", 'source("R/theme_pub.R")' not in text, "no source() dependency")
    forbidden = ("本篇可独立跑通", "这体现全系列", "接口只学一次", "作者代码通常长这样", "（即本文）")
    for phrase in forbidden:
        audit.add("Chapter prose", f"forbidden-{phrase}", phrase not in text, phrase)
    for stem in FIGURES:
        audit.add("Chapter figure", stem, f"../figures/{stem}.png" in text, stem)
    audit.add("Chapter figure", "anchor", "../figures/68-spencer-fig3a-original.png" in text, "anchor")


def audit_figures(root: Path, frozen: Path, audit: Audit) -> None:
    figures = root / "figures"
    for stem in FIGURES:
        for suffix in ("pdf", "png", "tiff"):
            path = figures / f"{stem}.{suffix}"
            audit.add("Figure file", f"{stem}.{suffix}", path.is_file() and path.stat().st_size > 10_000, path.stat().st_size if path.is_file() else "MISSING")
        png = figures / f"{stem}.png"
        tiff = figures / f"{stem}.tiff"
        if png.is_file():
            with Image.open(png) as image:
                dpi = image.info.get("dpi", (0, 0))
                audit.add("Figure raster", f"{stem}-png-size", image.width >= 1800 and image.height >= 1200, (image.width, image.height))
                audit.add("Figure raster", f"{stem}-png-dpi", min(dpi) >= 300, dpi)
        if tiff.is_file():
            with Image.open(tiff) as image:
                dpi = image.info.get("dpi", (0, 0))
                compression = image.tag_v2.get(259)
                audit.add("Figure raster", f"{stem}-tiff-dpi", min(dpi) >= 300, dpi)
                audit.add("Figure raster", f"{stem}-tiff-lzw", compression == 5, compression)
        pdf = figures / f"{stem}.pdf"
        if pdf.is_file():
            result = subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, text=True, timeout=30)
            audit.add("Figure text", f"{stem}-pdf-text", result.returncode == 0 and len(result.stdout.strip()) > 10, len(result.stdout))
            audit.add("Figure text", f"{stem}-english-only", not bool(re.search(r"[\u3400-\u9fff]", result.stdout)), "no CJK glyphs")

    anchor = figures / "68-spencer-fig3a-original.png"
    frozen_anchor = frozen / "spencer-fig3a-original.png"
    audit.add("Figure file", "anchor", anchor.is_file() and frozen_anchor.is_file() and sha256(anchor) == sha256(frozen_anchor), sha256(anchor) if anchor.is_file() else "MISSING")

    with tempfile.TemporaryDirectory(prefix="article68-replot-") as temporary:
        staged = Path(temporary) / "figures"
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(Path(temporary) / "mpl")
        result = subprocess.run(
            [
                sys.executable,
                str(frozen / "scripts" / "plot_article68_survival.py"),
                "--input-dir",
                str(frozen),
                "--figure-dir",
                str(staged),
            ],
            capture_output=True,
            text=True,
            timeout=240,
            env=environment,
        )
        audit.add("Reanalysis", "plot-script-exit", result.returncode == 0, result.stdout + result.stderr)
        for stem in FIGURES:
            staged_png = staged / f"{stem}.png"
            published_png = figures / f"{stem}.png"
            status = staged_png.is_file() and published_png.is_file() and pixel_sha(staged_png) == pixel_sha(published_png)
            audit.add("Reanalysis", f"{stem}-pixel-identical", status, pixel_sha(staged_png) if staged_png.is_file() else "MISSING")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    qa = args.qa_dir.resolve()
    qa.mkdir(parents=True, exist_ok=True)
    audit = Audit()
    audit_checksums(frozen, audit)
    audit_sources(frozen, audit)
    audit_cohort(frozen, audit)
    audit_models(frozen, audit)
    audit_leakage_and_prediction(frozen, audit)
    audit_chapter(root, audit)
    audit_figures(root, frozen, audit)
    report = audit.frame()
    report.to_csv(qa / "qa_report.tsv", sep="\t", index=False)
    passed = int(report["Status"].eq("PASS").sum())
    failed = int(report["Status"].eq("FAIL").sum())
    payload = {
        "article": 68,
        "status": "passed" if failed == 0 else "failed",
        "checks": len(report),
        "passed": passed,
        "failed": failed,
        "failed_checks": report.loc[report["Status"].eq("FAIL"), ["Category", "Check", "Detail"]].to_dict("records"),
    }
    (qa / "qa_report.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
