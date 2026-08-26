#!/usr/bin/env python3
"""Offline acceptance tests for the frozen Article 69 evidence bundle."""

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
    "69-cohort-observed-path",
    "69-dag-identification",
    "69-path-models",
    "69-effect-decomposition",
    "69-bootstrap-distributions",
    "69-positivity-balance",
    "69-sensitivity-analysis",
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
    audit.add("Bundle", "checksum-count", len(lines) == 38, len(lines))
    for line in lines:
        digest, relative = line.split("  ", 1)
        path = frozen / relative
        audit.add(
            "Checksum",
            relative,
            path.is_file() and sha256(path) == digest,
            digest,
        )
    manifest = json.loads((frozen / "bundle-manifest.json").read_text(encoding="utf-8"))
    audit.add("Bundle", "article", manifest.get("article") == 69, manifest.get("article"))
    audit.add("Bundle", "payload-count", manifest.get("payload_files") == 29, manifest.get("payload_files"))
    audit.add("Bundle", "script-count", manifest.get("script_files") == 6, manifest.get("script_files"))
    audit.add("Bundle", "environment-count", manifest.get("environment_files") == 2, manifest.get("environment_files"))


def audit_sources(frozen: Path, audit: Audit) -> None:
    manifest = json.loads((frozen / "source-manifest.json").read_text(encoding="utf-8"))
    expected = {
        "article": 69,
        "doi": "10.1126/science.aaz7015",
        "pmcid": "PMC8970537",
        "bioproject": "PRJNA770295",
        "repository_commit": "a95b1a020b890dbe93a960dec946e197b63b7d15",
    }
    for key, value in expected.items():
        audit.add("Source", key, manifest.get(key) == value, manifest.get(key))
    audit.add("Source", "profile-type", "last-known-taxon" in manifest.get("profile_type", ""), manifest.get("profile_type"))
    for name, (size, digest) in EXPECTED_SOURCES.items():
        record = manifest.get("resources", {}).get(name, {})
        audit.add("Source", f"{name}-bytes", record.get("bytes") == size, record.get("bytes"))
        audit.add("Source", f"{name}-sha256", record.get("sha256") == digest, record.get("sha256"))
        audit.add("Source", f"{name}-https", str(record.get("url", "")).startswith("https://"), record.get("url"))
    dictionary = frozen / "spencer-data-dictionary.xlsx"
    audit.add(
        "Source",
        "frozen-data-dictionary",
        dictionary.stat().st_size == 101_198
        and sha256(dictionary) == EXPECTED_SOURCES["spencer-data-dictionary.xlsx"][1],
        sha256(dictionary),
    )

    contract = json.loads((frozen / "methods-contract.json").read_text(encoding="utf-8"))
    audit.add("Contract", "article", contract.get("article") == 69, contract.get("article"))
    audit.add("Contract", "independent-unit", "patient" in contract.get("independent_unit", ""), contract.get("independent_unit"))
    audit.add("Contract", "fiber-threshold", ">=20" in contract.get("exposure", ""), contract.get("exposure"))
    audit.add("Contract", "mediator-transform", "log2" in contract.get("mediator", "") and "+ 25" in contract.get("mediator", ""), contract.get("mediator"))
    audit.add("Contract", "binary-outcome", "binary" in contract.get("outcome", ""), contract.get("outcome"))
    audit.add("Contract", "timing-limit", "contemporaneous" in contract.get("timing_limit", ""), contract.get("timing_limit"))


def audit_cohort(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    expected = {
        "article": 69,
        "public_wgs_samples": 167,
        "mediation_samples": 94,
        "sufficient_fiber": 23,
        "insufficient_fiber": 71,
        "responders": 60,
        "nonresponders": 34,
        "faecalibacterium_lkt_rows": 3,
        "ruminococcaceae_lkt_rows": 26,
        "faecalibacterium_detected": 93,
        "quality_sensitivity_samples": 93,
        "analysis_seed": 69_001,
        "plot_seed": 20_260_769,
    }
    for key, value in expected.items():
        audit.add("Metric", key, metrics.get(key) == value, metrics.get(key))
    audit.add("Metric", "fiber-threshold", near(metrics.get("fiber_threshold_g_day"), 20.0), metrics.get("fiber_threshold_g_day"))
    audit.add("Metric", "pseudocount", near(metrics.get("pseudocount_ppm"), 25.0), metrics.get("pseudocount_ppm"))

    cohort = pd.read_csv(frozen / "mediation-cohort.tsv", sep="\t")
    audit.add("Cohort", "rows", len(cohort) == 94, len(cohort))
    audit.add("Cohort", "unique-patient-samples", cohort["SampleID"].nunique() == 94, cohort["SampleID"].nunique())
    audit.add("Cohort", "binary-exposure", set(cohort["ExposureSufficient"]) == {0, 1}, cohort["ExposureSufficient"].value_counts().to_dict())
    audit.add("Cohort", "exposure-count", int(cohort["ExposureSufficient"].sum()) == 23, int(cohort["ExposureSufficient"].sum()))
    audit.add("Cohort", "binary-outcome", set(cohort["OutcomeResponder"]) == {0, 1}, cohort["OutcomeResponder"].value_counts().to_dict())
    audit.add("Cohort", "response-count", int(cohort["OutcomeResponder"].sum()) == 60, int(cohort["OutcomeResponder"].sum()))
    audit.add("Cohort", "fiber-complete", cohort["FiberGrams"].notna().all(), int(cohort["FiberGrams"].isna().sum()))
    audit.add("Cohort", "threshold-exposed", cohort.loc[cohort["ExposureSufficient"].eq(1), "FiberGrams"].ge(20).all(), cohort.loc[cohort["ExposureSufficient"].eq(1), "FiberGrams"].min())
    audit.add("Cohort", "threshold-unexposed", cohort.loc[cohort["ExposureSufficient"].eq(0), "FiberGrams"].lt(20).all(), cohort.loc[cohort["ExposureSufficient"].eq(0), "FiberGrams"].max())
    audit.add("Cohort", "response-labels", set(cohort["Response"]) == {"Responder", "Non_Responder"}, cohort["Response"].value_counts().to_dict())
    audit.add("Cohort", "covariates-complete", not cohort[["BMI", "PrimarySubtype", "AdvancedSubstage", "LDH"]].isna().any().any(), cohort[["BMI", "PrimarySubtype", "AdvancedSubstage", "LDH"]].isna().sum().to_dict())
    audit.add("Cohort", "subtype-levels", set(cohort["PrimarySubtype"]) == {"Cutaneous_or_unknown", "Mucosal_or_acral"}, cohort["PrimarySubtype"].value_counts().to_dict())
    audit.add("Cohort", "substage-levels", set(cohort["AdvancedSubstage"]) == {"Stage_M1C", "Stage_M1D"}, cohort["AdvancedSubstage"].value_counts().to_dict())
    audit.add("Cohort", "ldh-levels", set(cohort["LDH"]) == {"No", "Yes"}, cohort["LDH"].value_counts().to_dict())
    quality = cohort["QualitySensitivityPass"].astype(str).str.lower().eq("true")
    audit.add("Cohort", "quality-subset", int(quality.sum()) == 93 and int(cohort.loc[quality, "OutcomeResponder"].sum()) == 59, {"n": int(quality.sum()), "responders": int(cohort.loc[quality, "OutcomeResponder"].sum())})
    audit.add("Cohort", "faec-prevalence", int(cohort["FaecalibacteriumPPM"].gt(0).sum()) == 93, int(cohort["FaecalibacteriumPPM"].gt(0).sum()))
    audit.add("Cohort", "faec-transform", np.allclose(cohort["FaecalibacteriumLog2"], np.log2(cohort["FaecalibacteriumPPM"] + 25)), float(np.max(np.abs(cohort["FaecalibacteriumLog2"] - np.log2(cohort["FaecalibacteriumPPM"] + 25)))))
    audit.add("Cohort", "rumi-transform", np.allclose(cohort["RuminococcaceaeLog2"], np.log2(cohort["RuminococcaceaePPM"] + 25)), float(np.max(np.abs(cohort["RuminococcaceaeLog2"] - np.log2(cohort["RuminococcaceaePPM"] + 25)))))

    attrition = pd.read_csv(frozen / "cohort-attrition.tsv", sep="\t")
    audit.add("Attrition", "stage-order", attrition["StageOrder"].tolist() == [1, 2, 3, 4, 5], attrition["StageOrder"].tolist())
    audit.add("Attrition", "samples", attrition["Samples"].tolist() == [167, 167, 94, 94, 93], attrition["Samples"].tolist())
    audit.add("Attrition", "responders", attrition["Responders"].tolist() == [106, 106, 60, 60, 59], attrition["Responders"].tolist())

    outcome = pd.read_csv(frozen / "exposure-outcome-summary.tsv", sep="\t").set_index("ExposureSufficient")
    audit.add("Observed", "group-counts", outcome["Patients"].to_dict() == {0: 71, 1: 23}, outcome["Patients"].to_dict())
    audit.add("Observed", "responder-counts", outcome["Responders"].to_dict() == {0: 41, 1: 19}, outcome["Responders"].to_dict())
    audit.add("Observed", "response-rate-consistency", np.allclose(outcome["ResponseRate"], outcome["Responders"] / outcome["Patients"]), outcome["ResponseRate"].to_dict())

    matrix = pd.read_csv(frozen / "lkt-mediation-ppm.tsv.gz", sep="\t").set_index("LKTFeature")
    feature = pd.read_csv(frozen / "lkt-feature-map.tsv", sep="\t").set_index("LKTFeature")
    audit.add("Composition", "matrix-shape", matrix.shape == (225, 94), matrix.shape)
    audit.add("Composition", "feature-map", len(feature) == 225 and feature.index.is_unique, {"rows": len(feature), "unique": feature.index.is_unique})
    audit.add("Composition", "same-features", set(matrix.index) == set(feature.index), len(set(matrix.index) ^ set(feature.index)))
    audit.add("Composition", "same-samples", set(matrix.columns) == set(cohort["SampleID"]), len(set(matrix.columns) ^ set(cohort["SampleID"])))
    audit.add("Composition", "nonnegative-ppm", matrix.ge(0).all().all(), float(matrix.min().min()))
    faec_rows = feature.index[feature["Genus"].eq("g__Faecalibacterium")]
    rumi_rows = feature.index[feature["Family"].eq("f__Ruminococcaceae")]
    audit.add("Composition", "three-faec-rows", len(faec_rows) == 3, faec_rows.tolist())
    audit.add("Composition", "twenty-six-rumi-rows", len(rumi_rows) == 26, len(rumi_rows))
    faec = matrix.loc[faec_rows].sum(axis=0).loc[cohort["SampleID"]].to_numpy(float)
    rumi = matrix.loc[rumi_rows].sum(axis=0).loc[cohort["SampleID"]].to_numpy(float)
    audit.add("Composition", "faec-aggregation", np.allclose(faec, cohort["FaecalibacteriumPPM"]), float(np.max(np.abs(faec - cohort["FaecalibacteriumPPM"]))))
    audit.add("Composition", "rumi-aggregation", np.allclose(rumi, cohort["RuminococcaceaePPM"]), float(np.max(np.abs(rumi - cohort["RuminococcaceaePPM"]))))
    feature_audit = pd.read_csv(frozen / "mediator-feature-audit.tsv", sep="\t")
    audit.add("Composition", "feature-audit-rows", len(feature_audit) == 29, len(feature_audit))
    audit.add("Composition", "feature-audit-definitions", feature_audit["MediatorDefinition"].value_counts().to_dict() == {"Sensitivity: family Ruminococcaceae": 26, "Primary: genus Faecalibacterium": 3}, feature_audit["MediatorDefinition"].value_counts().to_dict())


def audit_models(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "model-metrics.json").read_text(encoding="utf-8"))
    exact = {
        "article": 69,
        "analysis_seed": 69_001,
        "primary_bootstrap": 5000,
        "sensitivity_bootstrap": 2000,
        "quadrature_nodes": 20,
    }
    for key, value in exact.items():
        audit.add("Model metric", key, metrics.get(key) == value, metrics.get(key))
    expected_float = {
        "mediator_exposure_beta": 1.1531102873736607,
        "mediator_exposure_p": 0.07812040445454206,
        "outcome_mediator_or": 1.0617037910997646,
        "outcome_mediator_p": 0.4979495300835234,
        "outcome_direct_exposure_or": 3.209621891895187,
        "outcome_direct_exposure_p": 0.07566091019075551,
        "total_exposure_or": 3.4841541340516535,
        "total_exposure_p": 0.05322512630126602,
        "primary_direct_rd": 0.20597590434195256,
        "primary_direct_ci_lower": -0.002766445992316389,
        "primary_direct_ci_upper": 0.40146477485097304,
        "primary_indirect_rd": 0.009872572305464455,
        "primary_indirect_ci_lower": -0.01848100422688635,
        "primary_indirect_ci_upper": 0.05676491420821034,
        "primary_total_rd": 0.21584847664741702,
        "primary_total_ci_lower": 0.017850237094758668,
        "primary_total_ci_upper": 0.40701301654231947,
        "primary_proportion_mediated": 0.0457384386436558,
        "propensity_min": 0.06438357158504718,
        "propensity_max": 0.4183825789403116,
        "maximum_iptw": 9.470251586380279,
        "acme_zero_rho": 0.05,
    }
    for key, value in expected_float.items():
        audit.add("Model metric", key, near(metrics.get(key), value, 1e-8), metrics.get(key))

    paths = pd.read_csv(frozen / "path-model-estimates.tsv", sep="\t")
    expected_models = {
        "Mediator model",
        "Outcome model",
        "Total-association model",
        "Outcome interaction sensitivity",
        "Continuous-fiber mediator model",
        "Continuous-fiber outcome model",
    }
    audit.add("Path model", "six-models", set(paths["Model"]) == expected_models, paths["Model"].value_counts().to_dict())
    audit.add("Path model", "rows", len(paths) == 40, len(paths))
    audit.add("Path model", "ci-order", (paths["CILower"] <= paths["Estimate"]).all() and (paths["Estimate"] <= paths["CIUpper"]).all(), "all coefficients")
    odds = paths.dropna(subset=["OddsRatio"])
    audit.add("Path model", "positive-odds-ratios", odds[["OddsRatio", "OddsRatioLower", "OddsRatioUpper"]].gt(0).all().all(), float(odds["OddsRatioLower"].min()))
    key = paths.set_index(["Model", "Term"])
    audit.add("Path model", "a-to-m", near(key.loc[("Mediator model", "A"), "Estimate"], expected_float["mediator_exposure_beta"], 1e-8), key.loc[("Mediator model", "A"), "Estimate"])
    audit.add("Path model", "m-to-y", near(key.loc[("Outcome model", "M"), "OddsRatio"], expected_float["outcome_mediator_or"], 1e-8), key.loc[("Outcome model", "M"), "OddsRatio"])
    audit.add("Path model", "direct-model-a", near(key.loc[("Outcome model", "A"), "OddsRatio"], expected_float["outcome_direct_exposure_or"], 1e-8), key.loc[("Outcome model", "A"), "OddsRatio"])
    audit.add("Path model", "total-model-a", near(key.loc[("Total-association model", "A"), "OddsRatio"], expected_float["total_exposure_or"], 1e-8), key.loc[("Total-association model", "A"), "OddsRatio"])

    summary = pd.read_csv(frozen / "gformula-effect-summary.tsv", sep="\t")
    variants = {"Primary", "Exposure-mediator interaction", "Sequencing-QC subset", "Ruminococcaceae mediator"}
    effects = {"P00", "P10", "P11", "P01", "Direct", "Indirect", "Total", "ProportionMediated"}
    audit.add("G-formula", "summary-rows", len(summary) == 32, len(summary))
    audit.add("G-formula", "variants", set(summary["Variant"]) == variants, summary["Variant"].value_counts().to_dict())
    audit.add("G-formula", "effects-per-variant", summary.groupby("Variant")["Effect"].apply(set).eq(effects).all(), summary.groupby("Variant")["Effect"].nunique().to_dict())
    audit.add("G-formula", "probabilities-in-range", summary.loc[summary["Effect"].isin({"P00", "P10", "P11", "P01"}), "Estimate"].between(0, 1).all(), "all standardized probabilities")
    for variant, group in summary.groupby("Variant"):
        values = group.set_index("Effect")["Estimate"]
        audit.add("G-formula identity", f"{variant}-direct", near(values["Direct"], values["P10"] - values["P00"], 1e-8), values["Direct"])
        audit.add("G-formula identity", f"{variant}-indirect", near(values["Indirect"], values["P11"] - values["P10"], 1e-8), values["Indirect"])
        audit.add("G-formula identity", f"{variant}-total", near(values["Total"], values["P11"] - values["P00"], 1e-8), values["Total"])
        audit.add("G-formula identity", f"{variant}-sum", near(values["Total"], values["Direct"] + values["Indirect"], 1e-8), values["Total"])

    primary = pd.read_csv(frozen / "primary-gformula-bootstrap.tsv.gz", sep="\t")
    audit.add("Bootstrap", "primary-rows", len(primary) == 5000, len(primary))
    audit.add("Bootstrap", "primary-iterations", primary["Iteration"].tolist() == list(range(1, 5001)), [primary["Iteration"].min(), primary["Iteration"].max()])
    audit.add("Bootstrap", "primary-variant", set(primary["Variant"]) == {"Primary"}, primary["Variant"].value_counts().to_dict())
    numeric = ["P00", "P10", "P11", "P01", "Direct", "Indirect", "Total", "ProportionMediated"]
    audit.add("Bootstrap", "primary-finite", np.isfinite(primary[numeric].to_numpy(float)).all(), int(primary[numeric].isna().sum().sum()))
    audit.add("Bootstrap", "primary-probability-range", ((primary[["P00", "P10", "P11", "P01"]] >= 0) & (primary[["P00", "P10", "P11", "P01"]] <= 1)).all().all(), {"min": float(primary[["P00", "P10", "P11", "P01"]].min().min()), "max": float(primary[["P00", "P10", "P11", "P01"]].max().max())})
    audit.add("Bootstrap", "primary-direct-identity", np.allclose(primary["Direct"], primary["P10"] - primary["P00"]), float(np.max(np.abs(primary["Direct"] - (primary["P10"] - primary["P00"])))))
    audit.add("Bootstrap", "primary-indirect-identity", np.allclose(primary["Indirect"], primary["P11"] - primary["P10"]), float(np.max(np.abs(primary["Indirect"] - (primary["P11"] - primary["P10"])))))
    audit.add("Bootstrap", "primary-total-identity", np.allclose(primary["Total"], primary["Direct"] + primary["Indirect"]), float(np.max(np.abs(primary["Total"] - primary["Direct"] - primary["Indirect"]))))

    sensitivity = pd.read_csv(frozen / "sensitivity-gformula-bootstrap.tsv.gz", sep="\t")
    audit.add("Bootstrap", "sensitivity-rows", len(sensitivity) == 6000, len(sensitivity))
    expected_valid = {"Exposure-mediator interaction": 1892, "Sequencing-QC subset": 2000, "Ruminococcaceae mediator": 2000}
    for variant, valid in expected_valid.items():
        group = sensitivity.loc[sensitivity["Variant"].eq(variant)]
        complete = group[numeric].notna().all(axis=1)
        audit.add("Bootstrap sensitivity", f"{variant}-rows", len(group) == 2000 and group["Iteration"].tolist() == list(range(1, 2001)), len(group))
        audit.add("Bootstrap sensitivity", f"{variant}-valid", int(complete.sum()) == valid, int(complete.sum()))
        audit.add("Bootstrap sensitivity", f"{variant}-all-or-none", (group[numeric].isna().nunique(axis=1) == 1).all(), group[numeric].isna().sum().to_dict())
        good = group.loc[complete]
        audit.add("Bootstrap sensitivity", f"{variant}-finite", np.isfinite(good[numeric].to_numpy(float)).all(), int(good[numeric].isna().sum().sum()))
        audit.add("Bootstrap sensitivity", f"{variant}-identity", np.allclose(good["Total"], good["Direct"] + good["Indirect"]), float(np.max(np.abs(good["Total"] - good["Direct"] - good["Indirect"]))))

    bootstrap_by_variant = {"Primary": primary}
    bootstrap_by_variant.update({variant: sensitivity.loc[sensitivity["Variant"].eq(variant)] for variant in expected_valid})
    for _, row in summary.iterrows():
        draws = bootstrap_by_variant[row["Variant"]][row["Effect"]]
        draws = draws[np.isfinite(draws)].to_numpy(float)
        audit.add("Summary reproduction", f"{row['Variant']}-{row['Effect']}-valid", int(row["ValidBootstrap"]) == len(draws), {"reported": row["ValidBootstrap"], "recomputed": len(draws)})
        audit.add("Summary reproduction", f"{row['Variant']}-{row['Effect']}-ci", near(row["CILower"], np.quantile(draws, 0.025), 1e-8) and near(row["CIUpper"], np.quantile(draws, 0.975), 1e-8), [row["CILower"], row["CIUpper"]])
        if row["Effect"] in {"Direct", "Indirect", "Total"}:
            pvalue = min(1.0, 2 * min(np.mean(draws <= 0), np.mean(draws >= 0)))
            audit.add("Summary reproduction", f"{row['Variant']}-{row['Effect']}-p", near(row["BootstrapP"], pvalue, 1e-8), {"reported": row["BootstrapP"], "recomputed": pvalue})

    individual = pd.read_csv(frozen / "individual-standardized-risks.tsv", sep="\t")
    audit.add("Standardization", "individual-rows", len(individual) == 94 and individual["SampleID"].nunique() == 94, len(individual))
    audit.add("Standardization", "individual-probabilities", ((individual[["P00", "P10", "P11", "P01"]] >= 0) & (individual[["P00", "P10", "P11", "P01"]] <= 1)).all().all(), "all patients")
    audit.add("Standardization", "individual-direct", np.allclose(individual["Direct"], individual["P10"] - individual["P00"]), float(np.max(np.abs(individual["Direct"] - individual["P10"] + individual["P00"]))))
    audit.add("Standardization", "individual-indirect", np.allclose(individual["Indirect"], individual["P11"] - individual["P10"]), float(np.max(np.abs(individual["Indirect"] - individual["P11"] + individual["P10"]))))
    audit.add("Standardization", "individual-total", np.allclose(individual["Total"], individual["Direct"] + individual["Indirect"]), float(np.max(np.abs(individual["Total"] - individual["Direct"] - individual["Indirect"]))))
    primary_point = summary.loc[summary["Variant"].eq("Primary")].set_index("Effect")["Estimate"]
    for effect in ("P00", "P10", "P11", "P01", "Direct", "Indirect", "Total"):
        audit.add("Standardization", f"mean-{effect}", near(individual[effect].mean(), primary_point[effect], 1e-8), {"mean": individual[effect].mean(), "summary": primary_point[effect]})

    mediate = pd.read_csv(frozen / "probit-mediate-summary.tsv", sep="\t")
    audit.add("Sensitivity", "mediate-effects", mediate["Effect"].tolist() == ["ACME control", "ACME treated", "ACME average", "ADE average", "Total effect", "Proportion mediated average"], mediate["Effect"].tolist())
    acme = mediate.set_index("Effect").loc["ACME average"]
    audit.add("Sensitivity", "acme-average", near(acme["Estimate"], 0.011567553, 1e-7) and near(acme["PValue"], 0.5232, 1e-8), acme.to_dict())
    rho = pd.read_csv(frozen / "residual-correlation-sensitivity.tsv", sep="\t")
    audit.add("Sensitivity", "rho-grid", len(rho) == 39 and near(rho["Rho"].min(), -0.95) and near(rho["Rho"].max(), 0.95) and np.allclose(np.diff(rho["Rho"]), 0.05), [rho["Rho"].min(), rho["Rho"].max()])
    audit.add("Sensitivity", "rho-average", np.allclose(rho["ACMEAverage"], (rho["ACMEControl"] + rho["ACMETreated"]) / 2), float(np.max(np.abs(rho["ACMEAverage"] - (rho["ACMEControl"] + rho["ACMETreated"]) / 2))))
    sensitivity_audit = pd.read_csv(frozen / "unmeasured-confounding-sensitivity-audit.tsv", sep="\t").set_index("Quantity")
    audit.add("Sensitivity", "control-zero-crossing", near(sensitivity_audit.loc["Residual correlation where ACME crosses zero (control)", "Value"], 0.05), sensitivity_audit.iloc[:, 0].to_dict())
    audit.add("Sensitivity", "treated-zero-crossing", near(sensitivity_audit.loc["Residual correlation where ACME crosses zero (treated)", "Value"], 0.10), sensitivity_audit.iloc[:, 0].to_dict())


def audit_overlap(frozen: Path, audit: Audit) -> None:
    cohort = pd.read_csv(frozen / "mediation-cohort.tsv", sep="\t").set_index("SampleID")
    overlap = pd.read_csv(frozen / "exposure-overlap.tsv", sep="\t").set_index("SampleID")
    audit.add("Overlap", "rows", len(overlap) == 94 and overlap.index.is_unique, len(overlap))
    audit.add("Overlap", "same-samples", set(overlap.index) == set(cohort.index), len(set(overlap.index) ^ set(cohort.index)))
    audit.add("Overlap", "propensity-range", overlap["Propensity"].between(0, 1, inclusive="neither").all(), [float(overlap["Propensity"].min()), float(overlap["Propensity"].max())])
    expected_weights = np.where(overlap["A"].eq(1), 1 / overlap["Propensity"], 1 / (1 - overlap["Propensity"]))
    audit.add("Overlap", "weight-formula", np.allclose(overlap["IPTW"], expected_weights), float(np.max(np.abs(overlap["IPTW"] - expected_weights))))
    audit.add("Overlap", "max-weight", near(overlap["IPTW"].max(), 9.470251586380279, 1e-8), overlap["IPTW"].max())
    summary = pd.read_csv(frozen / "exposure-overlap-summary.tsv", sep="\t").set_index("Exposure")
    audit.add("Overlap", "summary-groups", summary["N"].to_dict() == {0: 71, 1: 23}, summary["N"].to_dict())
    for exposure, group in overlap.groupby("A"):
        row = summary.loc[exposure]
        audit.add("Overlap summary", f"group-{exposure}-minimum", near(row["Minimum"], group["Propensity"].min(), 1e-8), row["Minimum"])
        audit.add("Overlap summary", f"group-{exposure}-median", near(row["Median"], group["Propensity"].median(), 1e-8), row["Median"])
        audit.add("Overlap summary", f"group-{exposure}-maximum", near(row["Maximum"], group["Propensity"].max(), 1e-8), row["Maximum"])
        audit.add("Overlap summary", f"group-{exposure}-weight", near(row["MaximumIPTW"], group["IPTW"].max(), 1e-8), row["MaximumIPTW"])
    balance = pd.read_csv(frozen / "exposure-balance.tsv", sep="\t")
    audit.add("Balance", "four-covariates", balance["Covariate"].tolist() == ["BMIz", "Mucosal", "StageM1D", "LDHHigh"], balance["Covariate"].tolist())
    audit.add("Balance", "finite", np.isfinite(balance[["SMDUnweighted", "SMDIPTW"]].to_numpy(float)).all(), balance.to_dict("records"))
    expected = {
        "BMIz": (-0.273439, -0.067677),
        "Mucosal": (-0.059887, 0.120715),
        "StageM1D": (0.258130, -0.030873),
        "LDHHigh": (-0.046841, -0.014183),
    }
    for covariate, (before, after) in expected.items():
        row = balance.set_index("Covariate").loc[covariate]
        audit.add("Balance", f"{covariate}-unweighted", near(row["SMDUnweighted"], before, 1e-5), row["SMDUnweighted"])
        audit.add("Balance", f"{covariate}-iptw", near(row["SMDIPTW"], after, 1e-5), row["SMDIPTW"])


def audit_chapter(root: Path, audit: Audit) -> None:
    chapter = root / "chapters" / "69-mediation-analysis.qmd"
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
    audit.add("Chapter", "timing-limit-explicit", "同一个基线窗口" in text and "时间顺序" in text, "same baseline window")
    audit.add("Chapter", "identification-explicit", "sequential ignorability" in text and "Positivity" in text, "identification assumptions")
    audit.add("Chapter", "estimand-explicit", "interventional-analog" in text and "risk-difference" in text, "estimand and scale")
    audit.add("Chapter", "bootstrap-refits", "5,000" in text and "重新拟合" in text, "patient bootstrap")
    audit.add("Chapter", "random-seeds", "set.seed(69001)" in text and "20260769" in text, "69001 / 20260769")
    audit.add("Chapter", "inline-theme", "theme_pub <- function" in text and "save_pub <- function" in text, "inline functions")
    audit.add("Chapter", "no-source-theme", 'source("R/theme_pub.R")' not in text, "no source() dependency")
    forbidden = ("本篇可独立跑通", "这体现全系列", "接口只学一次", "作者代码通常长这样", "（即本文）")
    for phrase in forbidden:
        audit.add("Chapter prose", f"forbidden-{phrase}", phrase not in text, phrase)
    for stem in FIGURES:
        audit.add("Chapter figure", stem, f"../figures/{stem}.png" in text, stem)
    audit.add("Chapter figure", "anchor", "../figures/69-spencer-fig3ab-original.png" in text, "anchor")


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

    anchor = figures / "69-spencer-fig3ab-original.png"
    frozen_anchor = frozen / "spencer-fig3ab-original.png"
    audit.add("Figure file", "anchor", anchor.is_file() and frozen_anchor.is_file() and sha256(anchor) == sha256(frozen_anchor), sha256(anchor) if anchor.is_file() else "MISSING")
    if anchor.is_file():
        with Image.open(anchor) as image:
            audit.add("Figure raster", "anchor-dimensions", image.width >= 1700 and image.height >= 600, (image.width, image.height))

    with tempfile.TemporaryDirectory(prefix="article69-replot-") as temporary:
        staged = Path(temporary) / "figures"
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(Path(temporary) / "mpl")
        result = subprocess.run(
            [
                sys.executable,
                str(frozen / "scripts" / "plot_article69_mediation.py"),
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
    audit_overlap(frozen, audit)
    audit_chapter(root, audit)
    audit_figures(root, frozen, audit)
    report = audit.frame()
    report.to_csv(qa / "qa_report.tsv", sep="\t", index=False)
    passed = int(report["Status"].eq("PASS").sum())
    failed = int(report["Status"].eq("FAIL").sum())
    payload = {
        "article": 69,
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
