#!/usr/bin/env python3
"""Offline acceptance tests for Article 71's frozen piecewise-SEM analysis."""

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
    "71-data-positivity",
    "71-prespecified-dag",
    "71-local-paths",
    "71-path-decomposition",
    "71-model-fit-dsep",
    "71-overlap-influence",
    "71-transport-sensitivity",
)

SOURCE_FILES = {
    "genera.tsv": (
        18_101_016,
        "c4a541fe198a147beccd72d52fb2ebbf75a8cdf75cb3df75f823290971409d3f",
    ),
    "metadata.tsv": (
        39_838,
        "f7396e3d6838b3b30f78b02bd568753757f84c956cd351966dbe654d50285376",
    ),
    "franzosa-fig1-original.png": (
        170_306,
        "7b81b865ae65659ad476d6f5210a3bc383b4eadad50d2ad7793a0b99df2450eb",
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
        return bool(
            np.isclose(float(value), expected, rtol=tolerance, atol=tolerance)
        )
    except (TypeError, ValueError):
        return False


def logical(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map({"true": True, "false": False})


@dataclass
class Check:
    category: str
    check: str
    status: bool
    detail: str


class Audit:
    def __init__(self) -> None:
        self.rows: list[Check] = []

    def add(
        self, category: str, check: str, status: bool, detail: object = ""
    ) -> None:
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
    lines = [
        line
        for line in checksum_file.read_text(encoding="utf-8").splitlines()
        if line
    ]
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

    manifest = json.loads(
        (frozen / "bundle-manifest.json").read_text(encoding="utf-8")
    )
    audit.add("Bundle", "article", manifest.get("article") == 71, manifest.get("article"))
    audit.add(
        "Bundle",
        "payload-count",
        manifest.get("payload_files") == 29,
        manifest.get("payload_files"),
    )
    audit.add(
        "Bundle",
        "script-count",
        manifest.get("script_files") == 6,
        manifest.get("script_files"),
    )
    audit.add(
        "Bundle",
        "environment-count",
        manifest.get("environment_files") == 2,
        manifest.get("environment_files"),
    )


def audit_sources(frozen: Path, audit: Audit) -> None:
    manifest = json.loads((frozen / "source-manifest.json").read_text())
    expected = {
        "article": 71,
        "dataset": "FRANZOSA_IBD_2019",
        "repository": "borenstein-lab/microbiome-metabolome-curated-data",
        "repository_commit": "89a519d8c832008fbc6e650453e83e2f04858d02",
        "paper_doi": "10.1038/s41564-018-0306-4",
        "resource_doi": "10.1038/s41522-022-00345-5",
        "genus_features": 11720,
        "profile_rows": 220,
        "independent_subjects": 220,
        "primary_subjects": 90,
        "validation_strict_subjects": 38,
        "validation_broad_subjects": 60,
    }
    for key, value in expected.items():
        audit.add("Source", key, manifest.get(key) == value, manifest.get(key))
    audit.add(
        "Source",
        "anchor-identity",
        manifest.get("anchor_figure") == "Franzosa et al. 2019 Figure 1",
        manifest.get("anchor_figure"),
    )
    audit.add(
        "Source",
        "faecalibacterium-feature",
        str(manifest.get("faecalibacterium_feature", "")).endswith(
            "g__Faecalibacterium"
        ),
        manifest.get("faecalibacterium_feature"),
    )
    for name, (size, digest) in SOURCE_FILES.items():
        record = manifest["resources"][name]
        audit.add("Source", f"{name}-bytes", record["bytes"] == size, record["bytes"])
        audit.add(
            "Source", f"{name}-sha256", record["sha256"] == digest, record["sha256"]
        )
        audit.add(
            "Source",
            f"{name}-https",
            record["url"].startswith("https://"),
            record["url"],
        )
    anchor = frozen / "franzosa-fig1-original.png"
    audit.add(
        "Source",
        "frozen-anchor",
        anchor.stat().st_size == SOURCE_FILES["franzosa-fig1-original.png"][0]
        and sha256(anchor) == SOURCE_FILES["franzosa-fig1-original.png"][1],
        sha256(anchor),
    )

    contract = json.loads((frozen / "methods-contract.json").read_text())
    audit.add("Contract", "article", contract["article"] == 71, contract["article"])
    audit.add(
        "Contract",
        "cross-sectional",
        "cross-sectional" in contract["design"],
        contract["design"],
    )
    audit.add(
        "Contract",
        "measured-nodes",
        "antibiotic" in contract["environment"]
        and "Shannon" in contract["microbiome"]
        and "calprotectin" in contract["phenotype"],
        contract,
    )
    audit.add(
        "Contract",
        "bootstrap",
        contract["bootstrap"] == 5000,
        contract["bootstrap"],
    )
    audit.add(
        "Contract",
        "causal-limit",
        "not identified causal mediation effects" in contract["interpretation_limit"],
        contract["interpretation_limit"],
    )


def audit_cohorts(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "analysis-metrics.json").read_text())
    expected = {
        "article": 71,
        "analysis_seed": 71001,
        "plot_seed": 20260771,
        "public_subjects": 220,
        "genus_features": 11720,
        "calprotectin_available": 153,
        "primary_subjects": 90,
        "primary_antibiotic_exposed": 13,
        "primary_controls": 20,
        "primary_cd": 46,
        "primary_uc": 24,
        "primary_exposed_controls": 0,
        "validation_strict_subjects": 38,
        "validation_antibiotic_exposed": 0,
    }
    for key, value in expected.items():
        audit.add("Input metric", key, metrics.get(key) == value, metrics.get(key))
    audit.add(
        "Input metric",
        "pseudocount",
        near(metrics["pseudocount"], 1e-6),
        metrics["pseudocount"],
    )

    all_samples = pd.read_csv(frozen / "all-sample-metrics.tsv", sep="\t")
    primary = pd.read_csv(frozen / "sem-primary-cohort.tsv", sep="\t")
    validation = pd.read_csv(frozen / "sem-validation-cohort.tsv", sep="\t")
    audit.add("Cohort", "all-rows", len(all_samples) == 220, len(all_samples))
    audit.add(
        "Cohort", "all-unique-subjects", all_samples["Subject"].is_unique, all_samples["Subject"].nunique()
    )
    audit.add("Cohort", "primary-rows", len(primary) == 90, len(primary))
    audit.add("Cohort", "validation-rows", len(validation) == 38, len(validation))
    audit.add("Cohort", "primary-subjects", primary["Subject"].is_unique, primary["Subject"].nunique())
    audit.add("Cohort", "validation-subjects", validation["Subject"].is_unique, validation["Subject"].nunique())
    audit.add(
        "Cohort",
        "disjoint-cohorts",
        set(primary["Subject"]).isdisjoint(validation["Subject"]),
        len(set(primary["Subject"]) & set(validation["Subject"])),
    )
    audit.add("Cohort", "primary-label", primary["Cohort"].eq("PRISM").all(), primary["Cohort"].unique())
    audit.add("Cohort", "validation-label", validation["Cohort"].eq("Validation").all(), validation["Cohort"].unique())
    audit.add(
        "Cohort",
        "primary-diagnoses",
        primary["Diagnosis"].value_counts().to_dict()
        == {"CD": 46, "UC": 24, "Control": 20},
        primary["Diagnosis"].value_counts().to_dict(),
    )
    audit.add(
        "Cohort",
        "primary-exposure",
        int(primary["Antibiotic"].sum()) == 13,
        int(primary["Antibiotic"].sum()),
    )
    audit.add(
        "Cohort",
        "no-exposed-controls",
        int(primary.loc[primary["Diagnosis"].eq("Control"), "Antibiotic"].sum()) == 0,
        primary.groupby("Diagnosis")["Antibiotic"].sum().to_dict(),
    )
    audit.add(
        "Cohort",
        "validation-no-exposure",
        int(validation["Antibiotic"].sum()) == 0,
        int(validation["Antibiotic"].sum()),
    )
    required = [
        "AgeZ",
        "Antibiotic",
        "Immunosuppressant",
        "Mesalamine",
        "Steroids",
        "CD",
        "UC",
        "LogCalprotectinZ",
        "ShannonZ",
        "LogFaecalibacteriumZ",
    ]
    audit.add(
        "Cohort",
        "complete-model-fields",
        primary[required].notna().all().all() and validation[required].notna().all().all(),
        required,
    )
    binary = ["Antibiotic", "Immunosuppressant", "Mesalamine", "Steroids", "CD", "UC"]
    audit.add(
        "Cohort",
        "binary-fields",
        primary[binary].isin([0, 1]).all().all()
        and validation[binary].isin([0, 1]).all().all(),
        binary,
    )
    audit.add(
        "Cohort",
        "diagnosis-indicators",
        ((primary["Diagnosis"].eq("CD")).astype(int) == primary["CD"]).all()
        and ((primary["Diagnosis"].eq("UC")).astype(int) == primary["UC"]).all(),
        "CD and UC indicators",
    )
    audit.add(
        "Cohort",
        "profile-sums",
        np.allclose(primary["ProfileSum"], 1, atol=1e-12)
        and np.allclose(validation["ProfileSum"], 1, atol=1e-12),
        (primary["ProfileSum"].min(), primary["ProfileSum"].max()),
    )
    audit.add(
        "Cohort",
        "calprotectin-transform",
        np.allclose(primary["LogCalprotectin"], np.log1p(primary["Fecal.Calprotectin"]), atol=1e-12)
        and np.allclose(validation["LogCalprotectin"], np.log1p(validation["Fecal.Calprotectin"]), atol=1e-12),
        "log1p",
    )
    audit.add(
        "Cohort",
        "faecalibacterium-transform",
        np.allclose(
            primary["Log10Faecalibacterium"],
            np.log10(primary["FaecalibacteriumRelative"] + 1e-6),
            atol=1e-12,
        ),
        "log10(relative abundance + 1e-6)",
    )

    standards = pd.read_csv(frozen / "standardization-parameters.tsv", sep="\t")
    audit.add("Scale", "four-parameters", len(standards) == 4, len(standards))
    standard_map = standards.set_index("Variable")
    for raw, zed in {
        "Age": "AgeZ",
        "Shannon": "ShannonZ",
        "LogCalprotectin": "LogCalprotectinZ",
        "Log10Faecalibacterium": "LogFaecalibacteriumZ",
    }.items():
        mean = float(standard_map.loc[raw, "Mean"])
        sd = float(standard_map.loc[raw, "SD"])
        audit.add(
            "Scale",
            f"{raw}-primary-parameters",
            near(primary[raw].mean(), mean) and near(primary[raw].std(ddof=1), sd),
            (mean, sd),
        )
        audit.add(
            "Scale",
            f"{raw}-primary-z",
            np.allclose(primary[zed], (primary[raw] - mean) / sd, atol=1e-12),
            zed,
        )
        audit.add(
            "Scale",
            f"{raw}-validation-reference",
            np.allclose(validation[zed], (validation[raw] - mean) / sd, atol=1e-12),
            zed,
        )

    attrition = pd.read_csv(frozen / "sample-attrition.tsv", sep="\t")
    audit.add(
        "Cohort",
        "attrition-counts",
        attrition["Subjects"].tolist() == [220, 153, 150, 128, 90, 38],
        attrition["Subjects"].tolist(),
    )
    overlap = pd.read_csv(frozen / "antibiotic-overlap-by-diagnosis.tsv", sep="\t")
    audit.add(
        "Cohort",
        "overlap-ledger",
        overlap[["Diagnosis", "Unexposed", "Exposed"]].to_dict("records")
        == [
            {"Diagnosis": "Control", "Unexposed": 20, "Exposed": 0},
            {"Diagnosis": "CD", "Unexposed": 38, "Exposed": 8},
            {"Diagnosis": "UC", "Unexposed": 19, "Exposed": 5},
        ],
        overlap.to_dict("records"),
    )


def audit_models(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "model-metrics.json").read_text())
    expected = {
        "article": 71,
        "analysis_seed": 71001,
        "plot_seed": 20260771,
        "primary_subjects": 90,
        "primary_antibiotic_exposed": 13,
        "bootstrap": 5000,
        "bootstrap_valid": 5000,
        "constrained_fisher_df": 2,
        "validation_antibiotic_exposed": 0,
    }
    for key, value in expected.items():
        audit.add("Model metric", key, metrics.get(key) == value, metrics.get(key))
    numeric_expected = {
        "shannon_a": -0.9037712391973121,
        "shannon_b": -0.1467041786134764,
        "shannon_direct": -0.5139175462673824,
        "shannon_indirect": 0.13258701730092537,
        "shannon_total": -0.3813305289664570,
        "shannon_indirect_ci_lower": -0.10929280094374684,
        "shannon_indirect_ci_upper": 0.36392085814760944,
        "shannon_indirect_p": 0.2836,
        "mediator_r2": 0.4443137130154484,
        "outcome_r2": 0.2981634594244468,
        "constrained_fisher_c": 4.655,
        "constrained_fisher_p": 0.098,
        "primary_aic": 462.062,
        "constrained_aic": 463.131,
        "reverse_aic": 462.062,
        "propensity_max_abs_coefficient": 19.379301618712365,
        "leave_one_out_indirect_min": 0.07455974069105198,
        "leave_one_out_indirect_max": 0.1636875092898683,
    }
    for key, value in numeric_expected.items():
        audit.add("Model metric", key, near(metrics[key], value, 1e-8), metrics[key])
    audit.add(
        "Model identity",
        "indirect-a-times-b",
        near(metrics["shannon_indirect"], metrics["shannon_a"] * metrics["shannon_b"]),
        metrics["shannon_indirect"],
    )
    audit.add(
        "Model identity",
        "total-direct-plus-indirect",
        near(metrics["shannon_total"], metrics["shannon_direct"] + metrics["shannon_indirect"]),
        metrics["shannon_total"],
    )
    audit.add(
        "Model identity",
        "reverse-same-aic",
        near(metrics["primary_aic"], metrics["reverse_aic"]),
        (metrics["primary_aic"], metrics["reverse_aic"]),
    )
    audit.add(
        "Model identity",
        "propensity-converged-but-separated",
        metrics["propensity_converged"] is True
        and metrics["propensity_max_abs_coefficient"] > 19,
        metrics["propensity_max_abs_coefficient"],
    )

    coefficients = pd.read_csv(frozen / "local-path-coefficients-hc3.tsv", sep="\t")
    audit.add("Local model", "coefficient-rows", len(coefficients) == 42, len(coefficients))
    key_rows = {
        ("Microbiome node", "Antibiotic"): (-0.9037712391973121, 0.295636625933503, 0.00301744610214472),
        ("Phenotype node", "ShannonZ"): (-0.1467041786134764, 0.134390519843115, 0.27823395591125),
        ("Phenotype node", "Antibiotic"): (-0.5139175462673824, 0.352423774704048, 0.148641768283997),
        ("Total-association node", "Antibiotic"): (-0.3813305289664570, 0.337737369997513, 0.262158853355584),
        ("Constrained phenotype node", "ShannonZ"): (-0.0696552533023216, 0.136446472153862, 0.61107620588278),
    }
    indexed = coefficients.set_index(["Model", "Term"])
    for key, (estimate, se, pvalue) in key_rows.items():
        row = indexed.loc[key]
        audit.add("Local model", f"{key[0]}-{key[1]}-estimate", near(row["Estimate"], estimate), row["Estimate"])
        audit.add("Local model", f"{key[0]}-{key[1]}-hc3-se", near(row["RobustSE"], se), row["RobustSE"])
        audit.add("Local model", f"{key[0]}-{key[1]}-p", near(row["PValue"], pvalue), row["PValue"])
        audit.add(
            "Local model",
            f"{key[0]}-{key[1]}-ci",
            near(row["CILower"], estimate - 1.96 * se, 1e-8)
            and near(row["CIUpper"], estimate + 1.96 * se, 1e-8),
            (row["CILower"], row["CIUpper"]),
        )

    effects = pd.read_csv(frozen / "path-effect-summary.tsv", sep="\t")
    audit.add("Bootstrap", "summary-rows", len(effects) == 10, len(effects))
    primary = effects.loc[effects["Model"].eq("Primary Shannon path")].set_index("Effect")
    audit.add(
        "Bootstrap",
        "primary-effects",
        set(primary.index) == {"A", "B", "Direct", "Indirect", "Total"},
        primary.index.tolist(),
    )
    audit.add(
        "Bootstrap",
        "primary-valid",
        primary["ValidBootstrap"].eq(5000).all()
        and primary["RequestedBootstrap"].eq(5000).all(),
        primary[["ValidBootstrap", "RequestedBootstrap"]].to_dict("index"),
    )
    audit.add(
        "Bootstrap",
        "indirect-summary",
        near(primary.loc["Indirect", "Estimate"], 0.132587017300925)
        and near(primary.loc["Indirect", "CILower"], -0.109292800943747)
        and near(primary.loc["Indirect", "CIUpper"], 0.363920858147609)
        and near(primary.loc["Indirect", "BootstrapP"], 0.2836),
        primary.loc["Indirect"].to_dict(),
    )
    faec = effects.loc[
        effects["Model"].eq("Faecalibacterium sensitivity")
    ].set_index("Effect")
    audit.add(
        "Bootstrap",
        "faecalibacterium-valid",
        faec["ValidBootstrap"].eq(2000).all()
        and faec["RequestedBootstrap"].eq(2000).all(),
        faec[["ValidBootstrap", "RequestedBootstrap"]].to_dict("index"),
    )
    audit.add(
        "Bootstrap",
        "faecalibacterium-indirect",
        near(faec.loc["Indirect", "Estimate"], 0.098302990157214)
        and near(faec.loc["Indirect", "CILower"], -0.0717121247594993)
        and near(faec.loc["Indirect", "CIUpper"], 0.273407434577941),
        faec.loc["Indirect"].to_dict(),
    )

    draws = pd.read_csv(frozen / "sem-path-bootstrap.tsv.gz", sep="\t")
    audit.add("Bootstrap", "draw-count", len(draws) == 5000, len(draws))
    audit.add(
        "Bootstrap",
        "iterations",
        draws["Iteration"].tolist() == list(range(1, 5001)),
        (draws["Iteration"].min(), draws["Iteration"].max()),
    )
    audit.add(
        "Bootstrap",
        "finite-draws",
        np.isfinite(draws[["A", "B", "Direct", "Indirect", "Total"]].to_numpy()).all(),
        "finite",
    )
    audit.add(
        "Bootstrap",
        "draw-path-identities",
        np.allclose(draws["Indirect"], draws["A"] * draws["B"], atol=1e-12)
        and np.allclose(draws["Total"], draws["Direct"] + draws["Indirect"], atol=1e-12)
        and draws["IdentityError"].abs().max() < 1e-12,
        draws["IdentityError"].abs().max(),
    )
    q025, q975 = draws["Indirect"].quantile([0.025, 0.975])
    empirical_p = 2 * min((draws["Indirect"] <= 0).mean(), (draws["Indirect"] >= 0).mean())
    audit.add(
        "Bootstrap",
        "summary-recomputed",
        near(q025, primary.loc["Indirect", "CILower"])
        and near(q975, primary.loc["Indirect", "CIUpper"])
        and near(empirical_p, primary.loc["Indirect", "BootstrapP"]),
        (q025, q975, empirical_p),
    )
    faec_draws = pd.read_csv(
        frozen / "faecalibacterium-path-bootstrap.tsv.gz", sep="\t"
    )
    audit.add("Bootstrap", "faec-draw-count", len(faec_draws) == 2000, len(faec_draws))
    audit.add(
        "Bootstrap",
        "faec-path-identities",
        np.allclose(faec_draws["Indirect"], faec_draws["A"] * faec_draws["B"], atol=1e-12)
        and np.allclose(faec_draws["Total"], faec_draws["Direct"] + faec_draws["Indirect"], atol=1e-12),
        "a*b and direct+indirect",
    )

    fits = pd.read_csv(frozen / "sem-fit-comparison.tsv", sep="\t")
    audit.add("Graph fit", "three-graphs", len(fits) == 3, len(fits))
    audit.add(
        "Graph fit",
        "aic-ledger",
        np.allclose(fits["AIC"], [462.062, 463.131, 462.062], atol=1e-9),
        fits["AIC"].tolist(),
    )
    audit.add(
        "Graph fit",
        "saturated-ledger",
        logical(fits["Saturated"]).tolist() == [True, False, True]
        and fits["IndependenceClaims"].tolist() == [0, 1, 0],
        fits[["Saturated", "IndependenceClaims"]].to_dict("records"),
    )
    constrained = fits.iloc[1]
    audit.add(
        "Graph fit",
        "fisher-c",
        near(constrained["FisherC"], 4.655)
        and int(constrained["FisherDF"]) == 2
        and near(constrained["FisherP"], 0.098),
        constrained.to_dict(),
    )
    dsep = pd.read_csv(frozen / "directed-separation-claims.tsv", sep="\t")
    audit.add("Graph fit", "one-dsep-claim", len(dsep) == 1, len(dsep))
    audit.add(
        "Graph fit",
        "omitted-direct-path",
        dsep.iloc[0]["IndependenceClaim"].startswith("LogCalprotectinZ ~ Antibiotic")
        and near(dsep.iloc[0]["PValue"], 0.0975441197705989),
        dsep.iloc[0].to_dict(),
    )

    leave = pd.read_csv(frozen / "leave-one-out-paths.tsv", sep="\t")
    audit.add("Sensitivity", "leave-one-out-count", len(leave) == 90, len(leave))
    audit.add(
        "Sensitivity",
        "leave-one-out-unique",
        leave["OmittedSample"].is_unique,
        leave["OmittedSample"].nunique(),
    )
    audit.add(
        "Sensitivity",
        "leave-one-out-range",
        near(leave["Indirect"].min(), 0.074559740691052)
        and near(leave["Indirect"].max(), 0.163687509289868),
        (leave["Indirect"].min(), leave["Indirect"].max()),
    )
    audit.add(
        "Sensitivity",
        "leave-one-out-identities",
        leave["IdentityError"].abs().max() < 1e-12,
        leave["IdentityError"].abs().max(),
    )

    transport = pd.read_csv(frozen / "outcome-path-transport.tsv", sep="\t")
    audit.add("Transport", "two-cohorts", len(transport) == 2, len(transport))
    audit.add(
        "Transport",
        "same-variable-counts",
        transport["N"].tolist() == [90, 38]
        and transport["AntibioticExposed"].tolist() == [13, 0],
        transport[["N", "AntibioticExposed"]].to_dict("records"),
    )
    audit.add(
        "Transport",
        "direction-change",
        near(transport.iloc[0]["Estimate"], -0.0696552533023216)
        and near(transport.iloc[1]["Estimate"], 0.0523204795602166)
        and np.sign(transport.iloc[0]["Estimate"]) != np.sign(transport.iloc[1]["Estimate"]),
        transport["Estimate"].tolist(),
    )

    positivity = pd.read_csv(frozen / "propensity-positivity-audit.tsv", sep="\t")
    pos = positivity.set_index("Quantity")["Value"]
    audit.add("Positivity", "control-zero", near(pos["Control antibiotic exposed"], 0), pos["Control antibiotic exposed"])
    audit.add("Positivity", "validation-zero", near(pos["Validation antibiotic exposed"], 0), pos["Validation antibiotic exposed"])
    audit.add("Positivity", "max-coefficient", near(pos["Maximum absolute coefficient"], 19.3793016187124), pos["Maximum absolute coefficient"])
    audit.add("Positivity", "weights-not-used", near(pos["Weights used for inference"], 0), pos["Weights used for inference"])

    diagnostics = pd.read_csv(frozen / "local-model-diagnostics.tsv", sep="\t")
    audit.add("Diagnostics", "three-local-models", len(diagnostics) == 3, len(diagnostics))
    audit.add("Diagnostics", "condition-numbers", diagnostics["ConditionNumber"].lt(8).all(), diagnostics["ConditionNumber"].tolist())
    audit.add("Diagnostics", "cook-ledger", diagnostics["CookAbove4OverN"].tolist() == [7, 7, 8], diagnostics["CookAbove4OverN"].tolist())
    vif = pd.read_csv(frozen / "variance-inflation.tsv", sep="\t")
    audit.add("Diagnostics", "vif-count", len(vif) == 15, len(vif))
    audit.add("Diagnostics", "vif-below-3.1", vif["VIF"].max() < 3.1, vif["VIF"].max())


def audit_software(frozen: Path, audit: Audit) -> None:
    versions = pd.read_csv(
        frozen / "software-versions-r.tsv", sep="\t"
    ).set_index("Package")["Version"].astype(str)
    expected = {
        "R": "4.4.1",
        "piecewiseSEM": "2.3.0.1",
        "sandwich": "3.1.0",
        "lmtest": "0.9.40",
        "car": "3.1.2",
        "jsonlite": "1.8.8",
    }
    for package, version in expected.items():
        audit.add("Software", package, versions.get(package) == version, versions.get(package))
    python_versions = json.loads((frozen / "software-versions-python.json").read_text())
    for package, key in (("python", "Python"), ("pandas", "pandas"), ("numpy", "numpy")):
        audit.add("Software", f"python-{package}", bool(python_versions.get(key)), python_versions.get(key))

    packages = pd.read_csv(
        frozen / "env" / "multiomics-r-packages.tsv", sep="\t"
    ).set_index("Package")["Version"].astype(str)
    environment_expected = {
        "R": "4.4.1",
        "piecewiseSEM": "2.3.0.1",
        "sandwich": "3.1.0",
        "lmtest": "0.9-40",
        "car": "3.1-2",
        "jsonlite": "1.8.8",
    }
    for package, version in environment_expected.items():
        audit.add("Environment", package, packages.get(package) == version, packages.get(package))


def audit_chapter(root: Path, audit: Audit) -> None:
    chapter = root / "chapters" / "71-structural-equation-model.qmd"
    audit.add("Chapter", "exists", chapter.is_file(), chapter)
    if not chapter.is_file():
        return
    text = chapter.read_text(encoding="utf-8")
    lower_text = text.lower()
    audit.add("Chapter", "published", "draft: false" in text, "draft: false")
    audit.add("Chapter", "eval-true", "eval: true" in text, "eval: true")
    audit.add("Chapter", "freeze-auto", "freeze: auto" in text, "freeze: auto")
    audit.add(
        "Chapter",
        "native-quarto-fences",
        text.count("```{r}") >= 9
        and text.count("```{bash}") >= 2
        and "~~~{" not in text,
        "executable code-cell fences",
    )
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
    audit.add(
        "Chapter",
        "extra-data-explicit",
        "你还需要什么数据" in text
        and "实测环境暴露" in text
        and "独立测量的表型" in text
        and "混杂因素" in text,
        "environment, phenotype, confounders",
    )
    audit.add(
        "Chapter",
        "cross-sectional-boundary",
        "不能识别因果中介" in text
        and "同一次横断面访视" in text
        and "治疗启动时间未知" in text,
        "temporal ordering boundary",
    )
    audit.add(
        "Chapter",
        "piecewise-definition",
        "局部条件模型" in text
        and "directed separation" in lower_text
        and "Fisher's C" in text,
        "piecewise SEM definition",
    )
    audit.add(
        "Chapter",
        "saturated-warning",
        "df=0" in text
        and "饱和图" in text
        and "不能验证模型正确" in text,
        "saturated model limitation",
    )
    audit.add(
        "Chapter",
        "positivity-warning",
        "没有 antibiotic-exposed control" in text
        and "validation 队列 0/38" in text
        and "不使用该权重" in text,
        "structural positivity",
    )
    audit.add(
        "Chapter",
        "direction-not-learned",
        "AIC 都是 462.062" in text
        and "不能从 AIC 学出方向" in text,
        "reverse orientation",
    )
    audit.add(
        "Chapter",
        "indirect-boundary",
        "indirect=0.133" in text
        and "−0.109–0.364" in text
        and "不报告 proportion mediated" in text,
        "uncertain inconsistent mediation",
    )
    audit.add(
        "Chapter",
        "sample-size-boundary",
        "每个参数 10 个样本" in text
        and "Monte Carlo" in text,
        "simulation-based sample size",
    )
    audit.add(
        "Chapter",
        "pls-pm-boundary",
        "PLS-PM" in text
        and "不会赋予因果识别" in text,
        "prediction is not identification",
    )
    audit.add(
        "Chapter",
        "composition-boundary",
        "compositional" in text
        and "relative abundance" in text
        and "Faecalibacterium" in text,
        "microbiome scale",
    )
    audit.add(
        "Chapter",
        "seeds",
        "set.seed(71001)" in text and "20260771" in text,
        "71001 / 20260771",
    )
    audit.add(
        "Chapter",
        "inline-theme",
        "theme_pub <- function" in text and "save_pub <- function" in text,
        "inline functions",
    )
    audit.add(
        "Chapter",
        "no-source-theme",
        'source("R/theme_pub.R")' not in text,
        "no source() dependency",
    )
    audit.add(
        "Chapter",
        "versions",
        all(
            version in text
            for version in (
                "piecewiseSEM 2.3.0.1",
                "sandwich 3.1.0",
                "lmtest 0.9-40",
                "car 3.1-2",
            )
        ),
        "software versions",
    )
    forbidden = (
        "本篇可独立跑通",
        "这体现全系列",
        "接口只学一次",
        "作者代码通常长这样",
        "（即本文）",
        "/media/desk16",
        "/tmp/article71",
    )
    for phrase in forbidden:
        audit.add("Chapter prose", f"forbidden-{phrase}", phrase not in text, phrase)
    for stem in FIGURES:
        audit.add(
            "Chapter figure",
            stem,
            f"../figures/{stem}.png" in text,
            stem,
        )
    audit.add(
        "Chapter figure",
        "anchor",
        "../figures/71-franzosa-fig1-original.png" in text,
        "anchor",
    )
    citation_keys = (
        "franzosa2019ibd",
        "muller2022curatedmultiomics",
        "lefcheck2016piecewisesem",
        "shipley2009confirmatory",
        "wolf2013samplesize",
        "gloor2017compositional",
        "imai2010general",
    )
    for key in citation_keys:
        audit.add("Chapter citation", key, f"@{key}" in text, key)


def audit_figures(root: Path, frozen: Path, audit: Audit) -> None:
    figures = root / "figures"
    for stem in FIGURES:
        for suffix in ("pdf", "png", "tiff"):
            path = figures / f"{stem}.{suffix}"
            audit.add(
                "Figure file",
                f"{stem}.{suffix}",
                path.is_file() and path.stat().st_size > 10_000,
                path.stat().st_size if path.is_file() else "MISSING",
            )
        png = figures / f"{stem}.png"
        tiff = figures / f"{stem}.tiff"
        if png.is_file():
            with Image.open(png) as image:
                dpi = image.info.get("dpi", (0, 0))
                audit.add(
                    "Figure raster",
                    f"{stem}-png-size",
                    image.width >= 1800 and image.height >= 1100,
                    (image.width, image.height),
                )
                audit.add("Figure raster", f"{stem}-png-dpi", min(dpi) >= 300, dpi)
        if tiff.is_file():
            with Image.open(tiff) as image:
                dpi = image.info.get("dpi", (0, 0))
                compression = image.tag_v2.get(259)
                audit.add("Figure raster", f"{stem}-tiff-dpi", min(dpi) >= 300, dpi)
                audit.add("Figure raster", f"{stem}-tiff-lzw", compression == 5, compression)
        pdf = figures / f"{stem}.pdf"
        if pdf.is_file():
            result = subprocess.run(
                ["pdftotext", str(pdf), "-"],
                capture_output=True,
                text=True,
                timeout=30,
            )
            audit.add(
                "Figure text",
                f"{stem}-pdf-text",
                result.returncode == 0 and len(result.stdout.strip()) > 10,
                len(result.stdout),
            )
            audit.add(
                "Figure text",
                f"{stem}-english-only",
                not bool(re.search(r"[\u3400-\u9fff]", result.stdout)),
                "no CJK glyphs",
            )

    anchor = figures / "71-franzosa-fig1-original.png"
    frozen_anchor = frozen / "franzosa-fig1-original.png"
    audit.add(
        "Figure file",
        "anchor",
        anchor.is_file()
        and frozen_anchor.is_file()
        and sha256(anchor) == sha256(frozen_anchor),
        sha256(anchor) if anchor.is_file() else "MISSING",
    )
    if anchor.is_file():
        with Image.open(anchor) as image:
            audit.add(
                "Figure raster",
                "anchor-dimensions",
                (image.width, image.height) == (2123, 710),
                (image.width, image.height),
            )

    with tempfile.TemporaryDirectory(prefix="article71-replot-") as temporary:
        staged = Path(temporary) / "figures"
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(Path(temporary) / "mpl")
        result = subprocess.run(
            [
                sys.executable,
                str(frozen / "scripts" / "plot_article71_sem.py"),
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
        audit.add(
            "Reanalysis",
            "plot-script-exit",
            result.returncode == 0,
            result.stdout + result.stderr,
        )
        for stem in FIGURES:
            staged_png = staged / f"{stem}.png"
            published_png = figures / f"{stem}.png"
            status = (
                staged_png.is_file()
                and published_png.is_file()
                and pixel_sha(staged_png) == pixel_sha(published_png)
            )
            audit.add(
                "Reanalysis",
                f"{stem}-pixel-identical",
                status,
                pixel_sha(staged_png) if staged_png.is_file() else "MISSING",
            )


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    qa = args.qa_dir.resolve()
    qa.mkdir(parents=True, exist_ok=True)
    audit = Audit()
    audit_checksums(frozen, audit)
    audit_sources(frozen, audit)
    audit_cohorts(frozen, audit)
    audit_models(frozen, audit)
    audit_software(frozen, audit)
    audit_chapter(root, audit)
    audit_figures(root, frozen, audit)
    report = audit.frame()
    report.to_csv(qa / "qa_report.tsv", sep="\t", index=False)
    passed = int(report["Status"].eq("PASS").sum())
    failed = int(report["Status"].eq("FAIL").sum())
    payload = {
        "article": 71,
        "status": "passed" if failed == 0 else "failed",
        "checks": len(report),
        "passed": passed,
        "failed": failed,
        "failed_checks": report.loc[
            report["Status"].eq("FAIL"), ["Category", "Check", "Detail"]
        ].to_dict("records"),
    }
    (qa / "qa_report.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
