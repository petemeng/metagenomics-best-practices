#!/usr/bin/env python3
"""Offline acceptance tests for the frozen Article 70 MR evidence bundle."""

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
    "70-study-design-assumptions",
    "70-instrument-harmonisation",
    "70-mr-method-comparison",
    "70-mr-scatter",
    "70-heterogeneity-pleiotropy",
    "70-leave-one-out",
    "70-outlier-sensitivity",
    "70-steiger-directionality",
)

SOURCE_RDATA = (
    900_980,
    "7a13b142efeafc0fee9b80888d60d1db4b06f175afe02e3a3c45c0bc11d63502",
)
SOURCE_PDF = (
    2_256_456,
    "f4db594e534dc54417755c5368a782ff139d429e6b1f17b64af97c77d074876e",
)


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
    audit.add("Bundle", "checksum-count", len(lines) == 41, len(lines))
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
    audit.add("Bundle", "article", manifest.get("article") == 70, manifest.get("article"))
    audit.add(
        "Bundle",
        "payload-count",
        manifest.get("payload_files") == 32,
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
        "article": 70,
        "package": "TwoSampleMR",
        "package_version": "0.7.9",
        "repository_tag": "v0.7.9",
        "repository_commit": "3d119f20d6fc164b0c7f710f5590fee9580f2c7b",
    }
    for key, value in expected.items():
        audit.add("Source", key, manifest.get(key) == value, manifest.get(key))
    exposure = manifest["exposure"]
    outcome = manifest["outcome"]
    audit.add("Source", "exposure-id", exposure["opengwas_id"] == "ieu-a-2", exposure)
    audit.add(
        "Source",
        "exposure-study",
        exposure["doi"] == "10.1038/nature14177"
        and exposure["maximum_sample_size"] == 339224,
        exposure,
    )
    audit.add("Source", "outcome-id", outcome["opengwas_id"] == "ieu-a-7", outcome)
    audit.add(
        "Source",
        "outcome-study",
        outcome["doi"] == "10.1038/ng.3396"
        and outcome["cases"] == 60801
        and outcome["controls"] == 123504
        and outcome["sample_size"] == 184305,
        outcome,
    )
    audit.add(
        "Source",
        "anchor",
        manifest["anchor"]["doi"] == "10.7554/eLife.34408"
        and manifest["anchor"]["figure"] == "Figure 1",
        manifest["anchor"],
    )

    resources = manifest["resources"]
    for name, (size, digest) in {
        "twosamplemr-vig-perform-mr.RData": SOURCE_RDATA,
        "hemani-mrbase-paper.pdf": SOURCE_PDF,
    }.items():
        record = resources[name]
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
    rdata = frozen / "twosamplemr-vig-perform-mr.RData"
    audit.add(
        "Source",
        "frozen-rdata",
        rdata.stat().st_size == SOURCE_RDATA[0] and sha256(rdata) == SOURCE_RDATA[1],
        sha256(rdata),
    )

    contract = json.loads((frozen / "methods-contract.json").read_text())
    audit.add("Contract", "article", contract["article"] == 70, contract["article"])
    audit.add(
        "Contract",
        "two-sample-design",
        "two-sample" in contract["design"],
        contract["design"],
    )
    audit.add(
        "Contract",
        "79-instruments",
        "79 genome-wide" in contract["instruments"],
        contract["instruments"],
    )
    audit.add(
        "Contract",
        "ld-limit-recorded",
        "does not preserve the LD reference-panel" in contract["instruments"],
        contract["instruments"],
    )
    audit.add(
        "Contract",
        "presso-2000",
        contract["mr_presso_null_distributions"] == 2000,
        contract["mr_presso_null_distributions"],
    )
    audit.add(
        "Contract",
        "prevalence-grid",
        np.allclose(contract["steiger_prevalence_grid"], [0.03, 0.06, 0.10, 0.20]),
        contract["steiger_prevalence_grid"],
    )
    audit.add(
        "Contract",
        "microbiome-limit",
        "microbiome GWAS" in contract["interpretation_limit"],
        contract["interpretation_limit"],
    )


def audit_inputs(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "analysis-metrics.json").read_text())
    expected = {
        "article": 70,
        "analysis_seed": 70001,
        "plot_seed": 20260770,
        "exposure_id": "ieu-a-2",
        "outcome_id": "ieu-a-7",
        "exposure_instruments": 79,
        "outcome_associations": 79,
        "snp_overlap": 79,
        "exposure_sample_size_min": 218359,
        "exposure_sample_size_max": 339205,
        "exposure_study_maximum_n": 339224,
        "outcome_sample_size": 184305,
        "outcome_cases": 60801,
        "outcome_controls": 123504,
        "two_sample_mr_version": "0.7.9",
    }
    for key, value in expected.items():
        audit.add("Input metric", key, metrics.get(key) == value, metrics.get(key))
    audit.add(
        "Input metric",
        "max-exposure-p",
        near(metrics["exposure_pvalue_max"], 4.97794e-08, 1e-10),
        metrics["exposure_pvalue_max"],
    )

    exposure = pd.read_csv(
        frozen / "bmi-instruments-raw.tsv.gz", sep="\t"
    )
    outcome = pd.read_csv(
        frozen / "chd-associations-raw.tsv.gz", sep="\t"
    )
    audit.add(
        "Input",
        "row-counts",
        len(exposure) == 79 and len(outcome) == 79,
        (len(exposure), len(outcome)),
    )
    audit.add(
        "Input",
        "unique-snps",
        exposure["SNP"].is_unique and outcome["SNP"].is_unique,
        (exposure["SNP"].nunique(), outcome["SNP"].nunique()),
    )
    audit.add(
        "Input",
        "same-snps",
        set(exposure["SNP"]) == set(outcome["SNP"]),
        len(set(exposure["SNP"]) ^ set(outcome["SNP"])),
    )
    audit.add(
        "Input",
        "genome-wide-exposure-p",
        exposure["pval.exposure"].lt(5e-8).all(),
        exposure["pval.exposure"].max(),
    )
    audit.add(
        "Input",
        "per-snp-exposure-n",
        exposure["samplesize.exposure"].between(218359, 339205).all(),
        (
            exposure["samplesize.exposure"].min(),
            exposure["samplesize.exposure"].max(),
        ),
    )
    audit.add(
        "Input",
        "outcome-n",
        outcome["samplesize.outcome"].eq(184305).all(),
        outcome["samplesize.outcome"].unique().tolist(),
    )
    allele_columns = [
        "effect_allele.exposure",
        "other_allele.exposure",
    ]
    audit.add(
        "Input",
        "exposure-alleles",
        exposure[allele_columns].isin(list("ACGT")).all().all(),
        exposure[allele_columns].stack().unique().tolist(),
    )
    audit.add(
        "Input",
        "exposure-eaf",
        exposure["eaf.exposure"].between(0, 1, inclusive="neither").all(),
        (
            exposure["eaf.exposure"].min(),
            exposure["eaf.exposure"].max(),
        ),
    )
    audit.add(
        "Input",
        "finite-effects",
        np.isfinite(
            exposure[["beta.exposure", "se.exposure"]].to_numpy(float)
        ).all()
        and np.isfinite(
            outcome[["beta.outcome", "se.outcome"]].to_numpy(float)
        ).all(),
        "exposure and outcome beta/se",
    )
    audit.add(
        "Input",
        "positive-se",
        exposure["se.exposure"].gt(0).all() and outcome["se.outcome"].gt(0).all(),
        (
            exposure["se.exposure"].min(),
            outcome["se.outcome"].min(),
        ),
    )

    attrition = pd.read_csv(frozen / "input-attrition.tsv", sep="\t")
    audit.add(
        "Attrition",
        "stage-order",
        attrition["StageOrder"].tolist() == [1, 2, 3, 4],
        attrition["StageOrder"].tolist(),
    )
    audit.add(
        "Attrition",
        "no-hidden-loss",
        attrition["SNPs"].tolist() == [79, 79, 79, 79],
        attrition["SNPs"].tolist(),
    )


def audit_harmonisation(frozen: Path, audit: Audit) -> pd.DataFrame:
    harmonised = pd.read_csv(
        frozen / "harmonised-instruments.tsv.gz", sep="\t"
    )
    audit.add(
        "Harmonisation",
        "rows-unique",
        len(harmonised) == 79 and harmonised["SNP"].is_unique,
        (len(harmonised), harmonised["SNP"].nunique()),
    )
    for column in (
        "remove",
        "palindromic",
        "ambiguous",
        "mr_keep",
        "OutcomeBetaFlipped",
    ):
        audit.add(
            "Harmonisation",
            f"{column}-logical",
            logical(harmonised[column]).notna().all(),
            logical(harmonised[column]).value_counts().to_dict(),
        )
    audit.add(
        "Harmonisation",
        "all-retained",
        logical(harmonised["mr_keep"]).sum() == 79,
        logical(harmonised["mr_keep"]).sum(),
    )
    audit.add(
        "Harmonisation",
        "four-palindromes",
        logical(harmonised["palindromic"]).sum() == 4,
        logical(harmonised["palindromic"]).sum(),
    )
    audit.add(
        "Harmonisation",
        "no-ambiguous",
        logical(harmonised["ambiguous"]).sum() == 0,
        logical(harmonised["ambiguous"]).sum(),
    )
    audit.add(
        "Harmonisation",
        "no-remove",
        logical(harmonised["remove"]).sum() == 0,
        logical(harmonised["remove"]).sum(),
    )
    audit.add(
        "Harmonisation",
        "no-outcome-flips",
        logical(harmonised["OutcomeBetaFlipped"]).sum() == 0,
        logical(harmonised["OutcomeBetaFlipped"]).sum(),
    )
    audit.add(
        "Harmonisation",
        "action-two",
        harmonised["action"].eq(2).all(),
        harmonised["action"].unique().tolist(),
    )

    ledger = pd.read_csv(frozen / "harmonisation-audit.tsv", sep="\t")
    counts = ledger.set_index("Quantity")["Count"].to_dict()
    expected_counts = {
        "Exposure instruments": 79,
        "Outcome associations": 79,
        "Matched SNPs": 79,
        "Retained for MR": 79,
        "Palindromic SNPs": 4,
        "Ambiguous palindromic SNPs": 0,
        "Incompatible alleles": 0,
        "Outcome beta sign flips": 0,
        "Proxy variants": 0,
    }
    audit.add(
        "Harmonisation",
        "ledger",
        counts == expected_counts,
        counts,
    )

    f_recomputed = (
        harmonised["beta.exposure"] / harmonised["se.exposure"]
    ) ** 2
    audit.add(
        "Strength",
        "f-recomputed",
        np.allclose(f_recomputed, harmonised["FStatistic"]),
        np.max(np.abs(f_recomputed - harmonised["FStatistic"])),
    )
    audit.add(
        "Strength",
        "minimum-f",
        near(harmonised["FStatistic"].min(), 29.9072265625, 1e-8),
        harmonised["FStatistic"].min(),
    )
    audit.add(
        "Strength",
        "no-f-below-ten",
        harmonised["FStatistic"].ge(10).all(),
        int(harmonised["FStatistic"].lt(10).sum()),
    )
    audit.add(
        "Strength",
        "r2-sum",
        near(harmonised["R2Exposure"].sum(), 0.015807160388406353, 1e-8),
        harmonised["R2Exposure"].sum(),
    )
    audit.add(
        "Strength",
        "oriented-exposure-positive",
        harmonised["BetaExposureOriented"].ge(0).all(),
        harmonised["BetaExposureOriented"].min(),
    )

    strength = pd.read_csv(frozen / "instrument-strength-summary.tsv", sep="\t")
    values = strength.set_index("Quantity")["Value"]
    expected_strength = {
        "Instrument count": 79,
        "Minimum F": 29.9072265625,
        "Median F": 40.6193777777778,
        "Mean F": 65.5650265385058,
        "Maximum F": 716.454444444444,
        "Instruments with F < 10": 0,
        "I2GX": 0.983490318679763,
        "Sum of approximate exposure R2": 0.0158071603884064,
        "Approximate exposure variance explained (%)": 1.58071603884064,
    }
    for key, value in expected_strength.items():
        audit.add("Strength summary", key, near(values[key], value, 1e-8), values[key])
    return harmonised


def audit_models(frozen: Path, harmonised: pd.DataFrame, audit: Audit) -> None:
    metrics = json.loads((frozen / "model-metrics.json").read_text())
    exact = {
        "article": 70,
        "analysis_seed": 70001,
        "plot_seed": 20260770,
        "instruments": 79,
        "retained_instruments": 79,
        "palindromic_instruments": 4,
        "ambiguous_instruments": 0,
        "outcome_beta_flips": 0,
        "mr_presso_distributions": 2000,
        "mr_presso_outlier_count": 2,
        "radial_nominal_outliers": 10,
        "steiger_per_snp_forward": 78,
    }
    for key, value in exact.items():
        audit.add("Model metric", key, metrics.get(key) == value, metrics.get(key))
    expected_float = {
        "minimum_f": 29.9072265625,
        "mean_f": 65.565026538505805,
        "i2gx": 0.98349031867976344,
        "exposure_r2_sum": 0.015807160388406353,
        "ivw_beta": 0.44590909695339842,
        "ivw_se": 0.058983018762925624,
        "ivw_p": 4.0320203273486104e-14,
        "ivw_or": 1.5619094777501918,
        "ivw_or_ci_lower": 1.3913887932512883,
        "ivw_or_ci_upper": 1.7533282059756294,
        "weighted_median_beta": 0.38700648154444744,
        "weighted_median_p": 9.231687155680036e-08,
        "egger_beta": 0.50249350972972906,
        "egger_p": 0.0008012589918497686,
        "ivw_q": 143.65084092503372,
        "ivw_q_p": 8.728419952475621e-06,
        "egger_intercept": -0.0017193040849560093,
        "egger_intercept_p": 0.6674265921126499,
        "mr_presso_corrected_beta": 0.4797865522442949,
        "mr_presso_distortion_p": 0.566,
        "leave_one_out_beta_min": 0.41303728251147015,
        "leave_one_out_beta_max": 0.4640005643577653,
    }
    for key, value in expected_float.items():
        audit.add(
            "Model metric", key, near(metrics.get(key), value, 1e-8), metrics.get(key)
        )
    audit.add(
        "Model metric",
        "presso-global",
        metrics["mr_presso_global_p"] == "<5e-04",
        metrics["mr_presso_global_p"],
    )
    audit.add(
        "Model metric",
        "presso-snps",
        metrics["mr_presso_outliers"] == ["rs6713510", "rs7903146"],
        metrics["mr_presso_outliers"],
    )
    audit.add(
        "Model metric",
        "steiger-grid-forward",
        metrics["steiger_all_prevalences_forward"] is True,
        metrics["steiger_all_prevalences_forward"],
    )

    estimates = pd.read_csv(frozen / "mr-estimates.tsv", sep="\t")
    expected_methods = {
        "Inverse variance weighted",
        "Weighted median",
        "MR Egger",
        "Simple mode",
        "Weighted mode",
    }
    audit.add(
        "MR estimate",
        "methods",
        set(estimates["method"]) == expected_methods and len(estimates) == 5,
        estimates["method"].tolist(),
    )
    audit.add(
        "MR estimate",
        "79-snps",
        estimates["nsnp"].eq(79).all(),
        estimates["nsnp"].tolist(),
    )
    audit.add(
        "MR estimate",
        "ci-identity",
        np.allclose(estimates["lo_ci"], estimates["b"] - 1.96 * estimates["se"])
        and np.allclose(estimates["up_ci"], estimates["b"] + 1.96 * estimates["se"]),
        "beta +/- 1.96 SE",
    )
    audit.add(
        "MR estimate",
        "or-identity",
        np.allclose(estimates["or"], np.exp(estimates["b"]))
        and np.allclose(estimates["or_lci95"], np.exp(estimates["lo_ci"]))
        and np.allclose(estimates["or_uci95"], np.exp(estimates["up_ci"])),
        "exp(beta and CI)",
    )
    audit.add(
        "MR estimate",
        "all-positive",
        estimates["lo_ci"].gt(0).all(),
        estimates[["method", "lo_ci"]].to_dict("records"),
    )
    ivw = estimates.set_index("method").loc["Inverse variance weighted"]
    audit.add(
        "MR estimate",
        "ivw-exact",
        near(ivw["b"], expected_float["ivw_beta"], 1e-8)
        and near(ivw["or"], expected_float["ivw_or"], 1e-8),
        ivw.to_dict(),
    )

    heterogeneity = pd.read_csv(frozen / "mr-heterogeneity.tsv", sep="\t")
    audit.add(
        "Heterogeneity",
        "two-estimators",
        set(heterogeneity["method"])
        == {"MR Egger", "Inverse variance weighted"},
        heterogeneity["method"].tolist(),
    )
    ivw_q = heterogeneity.set_index("method").loc["Inverse variance weighted"]
    audit.add(
        "Heterogeneity",
        "ivw-q",
        near(ivw_q["Q"], 143.650840925034, 1e-8)
        and int(ivw_q["Q_df"]) == 78
        and near(ivw_q["Q_pval"], 8.72841995247562e-06, 1e-8),
        ivw_q.to_dict(),
    )
    audit.add(
        "Heterogeneity",
        "i2-recomputed",
        np.allclose(
            heterogeneity["I2Percent"],
            np.maximum(
                0,
                100
                * (heterogeneity["Q"] - heterogeneity["Q_df"])
                / heterogeneity["Q"],
            ),
        ),
        heterogeneity["I2Percent"].tolist(),
    )

    egger = pd.read_csv(frozen / "egger-intercept.tsv", sep="\t").iloc[0]
    audit.add(
        "Pleiotropy",
        "egger-intercept",
        near(egger["egger_intercept"], -0.00171930408495601, 1e-8)
        and near(egger["pval"], 0.66742659211265, 1e-8),
        egger.to_dict(),
    )

    tests = pd.read_csv(frozen / "mr-presso-tests.tsv", sep="\t").set_index("Test")
    audit.add(
        "MR-PRESSO",
        "global",
        tests.loc["Global", "PValueText"] == "<5e-04"
        and int(tests.loc["Global", "OutlierCount"]) == 2,
        tests.loc["Global"].to_dict(),
    )
    audit.add(
        "MR-PRESSO",
        "distortion",
        near(tests.loc["Distortion", "Statistic"], -7.06094306570872, 1e-8)
        and near(tests.loc["Distortion", "PValueText"], 0.566, 1e-8),
        tests.loc["Distortion"].to_dict(),
    )
    presso_outliers = pd.read_csv(frozen / "mr-presso-outliers.tsv", sep="\t")
    flagged = presso_outliers.loc[
        logical(presso_outliers["DistortionOutlier"]), "SNP"
    ].tolist()
    audit.add(
        "MR-PRESSO",
        "two-flagged-snps",
        flagged == ["rs6713510", "rs7903146"],
        flagged,
    )
    presso = pd.read_csv(frozen / "mr-presso-estimates.tsv", sep="\t")
    audit.add(
        "MR-PRESSO",
        "raw-corrected",
        presso["Analysis"].tolist() == ["Raw", "Outlier-corrected"],
        presso["Analysis"].tolist(),
    )
    audit.add(
        "MR-PRESSO",
        "ci-or-identity",
        np.allclose(
            presso["CILower"], presso["Estimate"] - 1.96 * presso["StandardError"]
        )
        and np.allclose(
            presso["CIUpper"], presso["Estimate"] + 1.96 * presso["StandardError"]
        )
        and np.allclose(presso["OddsRatio"], np.exp(presso["Estimate"])),
        "beta/OR identities",
    )
    corrected = presso.set_index("Analysis").loc["Outlier-corrected"]
    audit.add(
        "MR-PRESSO",
        "corrected-exact",
        near(corrected["Estimate"], 0.479786552244295, 1e-8)
        and near(corrected["OddsRatio"], 1.6157294915502, 1e-8),
        corrected.to_dict(),
    )

    radial = pd.read_csv(frozen / "radial-ivw-outliers.tsv", sep="\t")
    expected_radial = {
        "rs1000940",
        "rs11727676",
        "rs13078960",
        "rs3849570",
        "rs6457796",
        "rs6567160",
        "rs6713510",
        "rs7550711",
        "rs7903146",
        "rs9304665",
    }
    audit.add(
        "Radial IVW",
        "ten-nominal-outliers",
        len(radial) == 10 and set(radial["SNP"]) == expected_radial,
        radial["SNP"].tolist(),
    )
    audit.add(
        "Radial IVW",
        "nominal-p",
        radial["p.value"].lt(0.05).all(),
        radial["p.value"].max(),
    )

    leave = pd.read_csv(frozen / "leave-one-out.tsv", sep="\t")
    leave_snp = leave.loc[leave["SNP"].ne("All")]
    audit.add(
        "Leave-one-out",
        "79-omissions-plus-all",
        len(leave) == 80 and len(leave_snp) == 79,
        (len(leave), len(leave_snp)),
    )
    audit.add(
        "Leave-one-out",
        "all-snps-once",
        set(leave_snp["SNP"]) == set(harmonised["SNP"])
        and leave_snp["SNP"].is_unique,
        leave_snp["SNP"].nunique(),
    )
    audit.add(
        "Leave-one-out",
        "ci-identity",
        np.allclose(leave["CILower"], leave["b"] - 1.96 * leave["se"])
        and np.allclose(leave["OR"], np.exp(leave["b"])),
        "all rows",
    )
    audit.add(
        "Leave-one-out",
        "positive-range",
        near(leave_snp["b"].min(), 0.41303728251147015, 1e-8)
        and near(leave_snp["b"].max(), 0.4640005643577653, 1e-8)
        and leave_snp["CILower"].gt(0).all(),
        (leave_snp["b"].min(), leave_snp["b"].max()),
    )

    single = pd.read_csv(frozen / "single-snp-estimates.tsv.gz", sep="\t")
    audit.add(
        "Single SNP",
        "wald-plus-two-all",
        len(single) == 81
        and single["SNP"].str.startswith("All -").sum() == 2,
        (len(single), single["SNP"].str.startswith("All -").sum()),
    )

    steiger = pd.read_csv(frozen / "steiger-directionality.tsv", sep="\t")
    audit.add(
        "Steiger",
        "prevalence-grid",
        np.allclose(steiger["AssumedCHDPrevalence"], [0.03, 0.06, 0.10, 0.20]),
        steiger["AssumedCHDPrevalence"].tolist(),
    )
    audit.add(
        "Steiger",
        "aggregate-forward",
        logical(steiger["CorrectCausalDirection"]).all(),
        steiger["CorrectCausalDirection"].tolist(),
    )
    audit.add(
        "Steiger",
        "r2-values",
        np.allclose(steiger["ExposureR2"], 0.0158071603884064)
        and steiger["OutcomeLiabilityR2"].between(
            0.0016984, 0.0016992
        ).all(),
        steiger[["ExposureR2", "OutcomeLiabilityR2"]].to_dict("records"),
    )
    per_snp = pd.read_csv(frozen / "steiger-per-snp.tsv", sep="\t")
    audit.add(
        "Steiger",
        "per-snp-79",
        len(per_snp) == 79
        and set(per_snp["SNP"]) == set(harmonised["SNP"]),
        len(per_snp),
    )
    audit.add(
        "Steiger",
        "per-snp-78-forward",
        logical(per_snp["ExposureExplainsMore"]).sum() == 78,
        logical(per_snp["ExposureExplainsMore"]).value_counts().to_dict(),
    )

    design = pd.read_csv(frozen / "design-audit.tsv", sep="\t")
    unresolved = design["Status"].str.contains(
        "Not encoded|not encoded", regex=True
    )
    audit.add(
        "Design audit",
        "ten-items",
        len(design) == 10 and design["Item"].is_unique,
        design["Item"].tolist(),
    )
    audit.add(
        "Design audit",
        "unresolved-visible",
        unresolved.sum() >= 4,
        design.loc[unresolved, "Item"].tolist(),
    )
    audit.add(
        "Design audit",
        "ld-reference-unresolved",
        design.loc[
            design["Item"].eq("LD clumping reference and release"), "Status"
        ].str.contains("Not encoded").all(),
        design.loc[
            design["Item"].eq("LD clumping reference and release"), "Status"
        ].tolist(),
    )

    versions = pd.read_csv(
        frozen / "software-versions-analysis.tsv", sep="\t"
    ).set_index("Package")["Version"].astype(str)
    expected_versions = {
        "R": "4.4.1",
        "TwoSampleMR": "0.7.9",
        "ieugwasr": "1.1.0",
        "MRPRESSO": "1.0",
        "RadialMR": "1.1",
        "MendelianRandomization": "0.10.0",
        "jsonlite": "1.8.8",
    }
    for package, version in expected_versions.items():
        audit.add(
            "Software",
            package,
            versions.get(package) == version,
            versions.get(package),
        )


def audit_chapter(root: Path, audit: Audit) -> None:
    chapter = root / "chapters" / "70-mendelian-randomization.qmd"
    audit.add("Chapter", "exists", chapter.is_file(), chapter)
    if not chapter.is_file():
        return
    text = chapter.read_text(encoding="utf-8")
    audit.add("Chapter", "published", "draft: false" in text, "draft: false")
    audit.add("Chapter", "eval-true", "eval: true" in text, "eval: true")
    audit.add("Chapter", "freeze-auto", "freeze: auto" in text, "freeze: auto")
    audit.add(
        "Chapter",
        "native-quarto-fences",
        text.count("```{r}") >= 10
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
        "外部宿主遗传 mGWAS" in text
        and "结局 GWAS" in text
        and "你还需要什么数据" in text,
        "external mGWAS and outcome GWAS",
    )
    audit.add(
        "Chapter",
        "teaching-example-boundary",
        "不是微生物导致冠心病的证据" in text,
        "BMI example boundary",
    )
    audit.add(
        "Chapter",
        "three-assumptions",
        all(
            phrase in text
            for phrase in ("Relevance（相关性）", "Independence（独立性）", "Exclusion restriction（排除限制）")
        ),
        "IV assumptions",
    )
    audit.add(
        "Chapter",
        "ld-limit",
        "LD reference/release" in text
        and "不能声称已重建 independence" in text,
        "LD provenance limitation",
    )
    audit.add(
        "Chapter",
        "sample-overlap",
        "participant overlap" in text and "样本完全独立时" in text,
        "two-sample compatibility",
    )
    audit.add(
        "Chapter",
        "microbiome-specific-audit",
        "profiler/database" in text
        and "reverse MR" in text
        and "colocalization" in text,
        "microbiome validity",
    )
    audit.add(
        "Chapter",
        "multiple-testing",
        "FDR/Bonferroni" in text
        and "relaxed P<1e-5" in text,
        "screen-wide inference",
    )
    audit.add(
        "Chapter",
        "seeds",
        "set.seed(70001)" in text and "20260770" in text,
        "70001 / 20260770",
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
                "TwoSampleMR 0.7.9",
                "ieugwasr 1.1.0",
                "MRPRESSO 1.0",
                "RadialMR 1.1",
            )
        ),
        "software versions",
    )
    audit.add(
        "Chapter",
        "methods-limit",
        "treated as unresolved rather than reconstructed retrospectively"
        in text,
        "unresolved design fields",
    )
    forbidden = (
        "本篇可独立跑通",
        "这体现全系列",
        "接口只学一次",
        "作者代码通常长这样",
        "（即本文）",
        "/media/desk16",
        "/tmp/article70",
    )
    for phrase in forbidden:
        audit.add(
            "Chapter prose",
            f"forbidden-{phrase}",
            phrase not in text,
            phrase,
        )
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
        "../figures/70-hemani-figure1-original.png" in text,
        "anchor",
    )
    citation_keys = (
        "hemani2018mrbase",
        "locke2015bmi",
        "nikpay2015cad",
        "sanderson2022mr",
        "bowden2015egger",
        "bowden2016median",
        "hartwig2017mode",
        "verbanck2018presso",
        "hemani2017steiger",
        "burgess2016overlap",
        "skrivankova2021strobemr",
        "wade2020microbiomemr",
        "kurilshikov2021mibiogen",
        "giambartolomei2014coloc",
    )
    for key in citation_keys:
        audit.add("Chapter citation", key, f"@{key}" in text, key)


def audit_environment(frozen: Path, audit: Audit) -> None:
    packages = pd.read_csv(
        frozen / "env" / "multiomics-r-packages.tsv", sep="\t"
    ).set_index("Package")["Version"].astype(str)
    expected = {
        "TwoSampleMR": "0.7.9",
        "ieugwasr": "1.1.0",
        "MRPRESSO": "1.0",
        "RadialMR": "1.1",
        "MendelianRandomization": "0.10.0",
        "png": "0.1-8",
    }
    for package, version in expected.items():
        audit.add(
            "Environment",
            package,
            packages.get(package) == version,
            packages.get(package),
        )


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
                audit.add(
                    "Figure raster",
                    f"{stem}-png-dpi",
                    min(dpi) >= 300,
                    dpi,
                )
        if tiff.is_file():
            with Image.open(tiff) as image:
                dpi = image.info.get("dpi", (0, 0))
                compression = image.tag_v2.get(259)
                audit.add(
                    "Figure raster",
                    f"{stem}-tiff-dpi",
                    min(dpi) >= 300,
                    dpi,
                )
                audit.add(
                    "Figure raster",
                    f"{stem}-tiff-lzw",
                    compression == 5,
                    compression,
                )
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
            audit.add(
                "Figure text",
                f"{stem}-no-superscript-glyphs",
                not bool(re.search(r"[⁻⁸]", result.stdout)),
                "portable Latin/math text",
            )

    anchor = figures / "70-hemani-figure1-original.png"
    frozen_anchor = frozen / "hemani-figure1-original.png"
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
                image.width >= 1500 and image.height >= 1500,
                (image.width, image.height),
            )

    with tempfile.TemporaryDirectory(prefix="article70-replot-") as temporary:
        staged = Path(temporary) / "figures"
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(Path(temporary) / "mpl")
        result = subprocess.run(
            [
                sys.executable,
                str(frozen / "scripts" / "plot_article70_mr.py"),
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
    audit_inputs(frozen, audit)
    harmonised = audit_harmonisation(frozen, audit)
    audit_models(frozen, harmonised, audit)
    audit_chapter(root, audit)
    audit_environment(frozen, audit)
    audit_figures(root, frozen, audit)
    report = audit.frame()
    report.to_csv(qa / "qa_report.tsv", sep="\t", index=False)
    passed = int(report["Status"].eq("PASS").sum())
    failed = int(report["Status"].eq("FAIL").sum())
    payload = {
        "article": 70,
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
