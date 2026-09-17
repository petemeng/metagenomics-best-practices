#!/usr/bin/env python3
"""Offline, evidence-level acceptance tests for Article 71."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


FIGURES = (
    "71-study-data-anchor",
    "71-data-positivity",
    "71-prespecified-dag",
    "71-local-paths",
    "71-path-decomposition",
    "71-model-fit-dsep",
    "71-overlap-influence",
    "71-transport-sensitivity",
    "71-power-positive-control",
    "71-mediation-rho-sensitivity",
)

PLOT_FUNCTIONS = (
    "study_anchor",
    "data_positivity",
    "prespecified_dag",
    "local_paths",
    "path_decomposition",
    "model_fit_dsep",
    "overlap_influence",
    "transport_sensitivity",
    "power_positive_control",
    "mediation_rho",
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
        rgba = image.convert("RGBA")
        digest = hashlib.sha256()
        digest.update(f"{rgba.width}x{rgba.height}".encode())
        digest.update(rgba.tobytes())
        return digest.hexdigest()


def near(value: object, expected: float, tolerance: float = 1e-8) -> bool:
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
            {
                "Category": [row.category for row in self.rows],
                "Check": [row.check for row in self.rows],
                "Status": ["PASS" if row.status else "FAIL" for row in self.rows],
                "Detail": [row.detail for row in self.rows],
            }
        )


def read_tsv(frozen: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(frozen / name, sep="\t")


def audit_bundle(frozen: Path, audit: Audit) -> None:
    checksum_file = frozen / "file-checksums.sha256"
    audit.add("Bundle", "checksum manifest exists", checksum_file.is_file(), checksum_file)
    if not checksum_file.is_file():
        return

    lines = [line for line in checksum_file.read_text().splitlines() if line]
    manifest = json.loads((frozen / "bundle-manifest.json").read_text())
    expected_count = (
        manifest["payload_files"]
        + manifest["script_files"]
        + manifest["environment_files"]
        + 1
    )
    audit.add("Bundle", "checksum entry count", len(lines) == expected_count, len(lines))
    for line in lines:
        digest, relative = line.split("  ", 1)
        path = frozen / relative
        audit.add(
            "Checksum",
            relative,
            path.is_file() and sha256(path) == digest,
            digest,
        )
    audit.add("Bundle", "article id", manifest.get("article") == 71, manifest.get("article"))
    source_work_dir = str(manifest.get("source_work_dir", ""))
    audit.add(
        "Bundle",
        "local build path omitted",
        not source_work_dir.startswith(("/", "~"))
        and "/media/" not in source_work_dir
        and "\\Users\\" not in source_work_dir,
        source_work_dir,
    )
    audit.add(
        "Bundle",
        "expanded evidence contract",
        all(
            phrase in manifest.get("contract", "")
            for phrase in ("power simulation", "positive control", "No publisher artwork")
        ),
        manifest.get("contract"),
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
        "genus_features": 11_720,
        "profile_rows": 220,
        "independent_subjects": 220,
        "primary_subjects": 90,
        "validation_strict_subjects": 38,
    }
    for key, expected_value in expected.items():
        audit.add("Source", key, manifest.get(key) == expected_value, manifest.get(key))

    resources = manifest.get("resources", {})
    audit.add(
        "Source",
        "only public data tables are downloaded",
        set(resources) == set(SOURCE_FILES),
        sorted(resources),
    )
    for name, (expected_bytes, expected_sha) in SOURCE_FILES.items():
        record = resources.get(name, {})
        audit.add("Source", f"{name} bytes", record.get("bytes") == expected_bytes, record)
        audit.add("Source", f"{name} sha256", record.get("sha256") == expected_sha, record)
        audit.add(
            "Source", f"{name} pinned https", str(record.get("url", "")).startswith("https://"), record
        )
    audit.add(
        "Rights",
        "publisher artwork exclusion",
        "downloaded or redistributed" in manifest.get("publisher_figure_policy", "")
        and manifest.get("publisher_figure_policy", "").startswith("No "),
        manifest.get("publisher_figure_policy"),
    )
    audit.add(
        "Rights",
        "no publisher image in frozen root",
        not any("fig1-original" in path.name.lower() for path in frozen.rglob("*")),
        "frozen bundle scan",
    )

    contract = json.loads((frozen / "methods-contract.json").read_text())
    for key, expected_value in {
        "article": 71,
        "bootstrap": 5000,
        "power_simulation_repetitions": 1000,
        "positive_control_subjects": 500,
        "positive_control_bootstrap": 2000,
    }.items():
        audit.add("Contract", key, contract.get(key) == expected_value, contract.get(key))
    audit.add(
        "Contract",
        "cross-sectional interpretation limit",
        "not identified causal mediation" in contract.get("interpretation_limit", ""),
        contract.get("interpretation_limit"),
    )


def audit_cohorts(frozen: Path, audit: Audit) -> None:
    primary = read_tsv(frozen, "sem-primary-cohort.tsv")
    validation = read_tsv(frozen, "sem-validation-cohort.tsv")
    checks = {
        "primary n=90": len(primary) == 90,
        "primary 13 exposed": int(primary["Antibiotic"].sum()) == 13,
        "primary sample unique": primary["Sample"].is_unique,
        "primary subject unique": primary["Subject"].is_unique,
        "primary 20 controls": int(primary["Diagnosis"].eq("Control").sum()) == 20,
        "primary zero exposed controls": int(
            (primary["Diagnosis"].eq("Control") & primary["Antibiotic"].eq(1)).sum()
        )
        == 0,
        "validation n=38": len(validation) == 38,
        "validation zero exposed": int(validation["Antibiotic"].sum()) == 0,
        "profile closure": np.allclose(primary["ProfileSum"], 1.0, atol=1e-10),
        "primary Shannon z mean": near(primary["ShannonZ"].mean(), 0),
        "primary Shannon z sd": near(primary["ShannonZ"].std(ddof=1), 1),
    }
    for name, status in checks.items():
        audit.add("Cohort", name, status, name)

    attrition = read_tsv(frozen, "sample-attrition.tsv")
    audit.add(
        "Cohort",
        "attrition sequence",
        attrition["Subjects"].tolist() == [220, 153, 150, 128, 90, 38],
        attrition["Subjects"].tolist(),
    )
    overlap = read_tsv(frozen, "antibiotic-overlap-by-diagnosis.tsv")
    audit.add(
        "Cohort",
        "diagnosis overlap table",
        overlap[["Unexposed", "Exposed"]].to_numpy().tolist()
        == [[20, 0], [38, 8], [19, 5]],
        overlap.to_dict("records"),
    )
    comparison = read_tsv(frozen, "prism-complete-case-comparison.tsv")
    audit.add("Missingness", "nine comparison rows", len(comparison) == 9, len(comparison))
    cal = comparison.loc[comparison["Variable"].eq("log1p fecal calprotectin")].iloc[0]
    audit.add(
        "Missingness",
        "excluded calprotectin observed 3/65",
        int(cal["ExcludedObserved"]) == 3 and near(cal["ExcludedMissingPct"], 95.3846153846),
        cal.to_dict(),
    )


def audit_models(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "model-metrics.json").read_text())
    expected = {
        "primary_subjects": 90,
        "primary_antibiotic_exposed": 13,
        "bootstrap": 5000,
        "bootstrap_valid": 5000,
        "shannon_a": -0.903771239197312,
        "shannon_b": -0.146704178613476,
        "shannon_direct": -0.513917546267382,
        "shannon_indirect": 0.132587017300925,
        "shannon_total": -0.381330528966457,
        "shannon_indirect_ci_lower": -0.109292800943747,
        "shannon_indirect_ci_upper": 0.363920858147609,
        "shannon_indirect_p": 0.2836,
        "constrained_fisher_p": 0.0975441197705992,
        "primary_aic": 462.062126482154,
        "reverse_aic": 462.062126482154,
        "validation_antibiotic_exposed": 0,
        "first_tested_n_with_b_power_80": 500,
        "positive_control_n": 500,
        "medsens_rho_zero": -0.15,
    }
    for key, expected_value in expected.items():
        actual = metrics.get(key)
        status = actual == expected_value if isinstance(expected_value, int) else near(actual, expected_value)
        audit.add("Model", key, status, actual)

    effects = read_tsv(frozen, "path-effect-summary.tsv")
    primary = effects.loc[effects["Model"].eq("Primary Shannon path")]
    audit.add(
        "Model",
        "five primary path summaries",
        primary["Effect"].tolist() == ["A", "B", "Direct", "Indirect", "Total"],
        primary["Effect"].tolist(),
    )
    identity = metrics["shannon_total"] - metrics["shannon_direct"] - metrics["shannon_indirect"]
    audit.add("Model", "linear path identity", abs(identity) < 1e-10, identity)

    fit = read_tsv(frozen, "sem-fit-comparison.tsv")
    forward = fit.loc[fit["Model"].eq("Primary partial-path model")].iloc[0]
    reverse = fit.loc[fit["Model"].eq("Reverse cross-sectional orientation")].iloc[0]
    constrained = fit.loc[fit["Model"].eq("Constrained microbiome-only path")].iloc[0]
    audit.add("Graph", "forward/reverse AIC identical", near(forward.AIC, reverse.AIC), (forward.AIC, reverse.AIC))
    audit.add("Graph", "saturated graphs have zero claims", int(forward.IndependenceClaims) == 0 and int(reverse.IndependenceClaims) == 0, fit.to_dict("records"))
    audit.add("Graph", "constrained graph one claim", int(constrained.IndependenceClaims) == 1, constrained.to_dict())

    positivity = read_tsv(frozen, "propensity-positivity-audit.tsv")
    pos = dict(zip(positivity["Quantity"], positivity["Value"]))
    audit.add("Positivity", "GLM converged", near(pos["GLM converged"], 1), pos["GLM converged"])
    audit.add("Positivity", "large separated coefficient", near(pos["Maximum absolute coefficient"], 19.3793016187124), pos["Maximum absolute coefficient"])
    audit.add("Positivity", "minimum propensity", near(pos["Minimum fitted propensity"], 2.25594627914396e-09), pos["Minimum fitted propensity"])
    audit.add("Positivity", "weights not used", near(pos["Weights used for inference"], 0), pos["Weights used for inference"])


def audit_simulations(frozen: Path, audit: Audit) -> None:
    power = read_tsv(frozen, "path-power-simulation.tsv")
    audit.add("Simulation", "power grid has 16 rows", len(power) == 16, len(power))
    audit.add("Simulation", "1,000 repetitions per scenario", power["Repetitions"].eq(1000).all(), power["Repetitions"].unique())
    b_path = power.loc[power["Target"].eq("B path (HC3)")]
    n90 = b_path.loc[b_path["SampleSize"].eq(90)].iloc[0]
    n500 = b_path.loc[b_path["SampleSize"].eq(500)].iloc[0]
    audit.add("Simulation", "n=90 power 0.180", near(n90.Power, 0.18), n90.to_dict())
    audit.add("Simulation", "n=500 power 0.825", near(n500.Power, 0.825), n500.to_dict())
    first = b_path.loc[b_path["Power"].ge(0.8), "SampleSize"].min()
    audit.add("Simulation", "first tested n above 80% is 500", int(first) == 500, first)

    positive = read_tsv(frozen, "positive-control-paths.tsv")
    overlap = read_tsv(frozen, "positive-control-overlap.tsv")
    audit.add("Positive control", "overlap in all strata", (overlap[["Unexposed", "Exposed"]] > 0).all().all(), overlap.to_dict("records"))
    for effect in ("A", "B", "Indirect"):
        row = positive.loc[positive["Effect"].eq(effect)].iloc[0]
        excludes_zero = row.CIUpper < 0 if effect in {"A", "B"} else row.CILower > 0
        covers_truth = row.CILower <= row.TrueValue <= row.CIUpper
        audit.add("Positive control", f"{effect} excludes zero", excludes_zero, row.to_dict())
        audit.add("Positive control", f"{effect} covers truth", covers_truth, row.to_dict())
    positive_audit = read_tsv(frozen, "positive-control-audit.tsv")
    aic = positive_audit.loc[positive_audit["Criterion"].eq("Forward versus reverse AIC")].iloc[0]
    audit.add("Positive control", "AIC still identical", aic.Result == "Identical", aic.to_dict())

    rho = read_tsv(frozen, "mediation-rho-summary.tsv")
    rho_map = dict(zip(rho["Quantity"], rho["Value"]))
    audit.add(
        "Sensitivity",
        "point indirect crosses at rho -0.15",
        near(rho_map["Residual rho where point indirect effect is zero"], -0.15),
        rho_map,
    )
    audit.add(
        "Sensitivity",
        "interval crosses at rho zero",
        near(rho_map["Indirect interval at rho=0 crosses zero"], 1),
        rho_map,
    )


def audit_figures(root: Path, frozen: Path, audit: Audit) -> None:
    figure_dir = root / "figures"
    for stem in FIGURES:
        for suffix in ("png", "pdf", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            audit.add("Figure", f"{stem}.{suffix} exists", path.is_file(), path)
        png = figure_dir / f"{stem}.png"
        if png.is_file():
            with Image.open(png) as image:
                width, height = image.size
                dpi = image.info.get("dpi", (0, 0))[0]
            audit.add("Figure", f"{stem} dimensions", width >= 2800 and height >= 1600, (width, height))
            audit.add("Figure", f"{stem} dpi", dpi >= 350, dpi)

    plot_script = frozen / "scripts" / "plot_article71_sem.R"
    if not plot_script.is_file():
        audit.add("Figure", "frozen plot script", False, plot_script)
        return
    with tempfile.TemporaryDirectory(prefix="article71-figures-") as temp_dir:
        command = [
            "Rscript",
            str(plot_script),
            "--prepared-dir",
            str(frozen),
            "--analysis-dir",
            str(frozen),
            "--figure-dir",
            temp_dir,
        ]
        environment = os.environ.copy()
        result = subprocess.run(command, text=True, capture_output=True, env=environment)
        audit.add("Figure", "offline ggplot rerender", result.returncode == 0, result.stderr[-2000:])
        if result.returncode == 0:
            for stem in FIGURES:
                expected = figure_dir / f"{stem}.png"
                observed = Path(temp_dir) / f"{stem}.png"
                audit.add(
                    "Figure reproducibility",
                    stem,
                    expected.is_file() and observed.is_file() and pixel_sha(expected) == pixel_sha(observed),
                    "pixel SHA-256",
                )


def audit_article(root: Path, frozen: Path, audit: Audit) -> None:
    article = root / "chapters" / "71-structural-equation-model.qmd"
    text = article.read_text(encoding="utf-8")
    required_sections = (
        "## 这一章用什么数据",
        "## Shannon 到底算的是什么",
        "## 先看谁能和谁比较",
        "## 两个队列共用一把尺子",
        "## 把三条箭头写成回归",
        "## 间接效应要整条路径一起重拟合",
        "## 删掉 direct path 试试",
        "## 敏感性与外部队列能覆盖到哪里",
        "## 当流程通过时，输出应该长什么样",
        "## Key Takeaways",
        "## 附录 A：Methods / Results 模板",
        "## 附录 B：换成你自己的数据",
        "## 附录 C：十张图如何复现",
    )
    positions = [text.find(section) for section in required_sections]
    audit.add("Article", "linear pipeline sections present", all(position >= 0 for position in positions), positions)
    audit.add("Article", "linear pipeline section order", positions == sorted(positions), positions)
    audit.add("Article", "expected_images is 10", "expected_images: 10" in text, "front matter")

    for forbidden in (
        "piecewise SEM 的核心理论",
        "## 常见坑",
        "## 出版级美化",
        "卡在三道门",
        "节点合同",
        "尺度合同",
        "71-franzosa-fig1-original",
        "franzosa-fig1-original.png",
    ):
        audit.add("Article", f"forbidden phrase absent: {forbidden}", forbidden not in text, forbidden)
    audit.add("Article", "no source() dependency", "source(" not in text, "source(")
    audit.add("Article", "limited 不等于 construction", text.count("不等于") <= 2, text.count("不等于"))
    audit.add("Article", "ten unique figure references", all(f"../figures/{stem}.png" in text for stem in FIGURES), FIGURES)
    audit.add("Article", "no publisher figure URL", "nature-assets" not in text and "springernature" not in text, "rights surface")

    required_citations = (
        "franzosa2019ibd",
        "muller2022curatedmultiomics",
        "gloor2017compositional",
        "petersen2012positivity",
        "imai2010general",
        "tingley2014mediation",
        "shipley2009confirmatory",
        "verma1990equivalence",
        "shi2021cmaverse",
        "sohn2019compositional",
        "zhang2021microhima",
        "yue2022ldmmed",
        "wolf2013samplesize",
    )
    for key in required_citations:
        audit.add("Citation", key, f"@{key}" in text, key)

    plot_script_text = (root / "scripts" / "plot_article71_sem.R").read_text()
    for stem, function in zip(FIGURES, PLOT_FUNCTIONS, strict=True):
        audit.add("Plot code", f"{function} function", f"{function} <- function" in plot_script_text, function)
        audit.add("Plot code", f"{stem} export", f'"{stem}"' in plot_script_text, stem)
        audit.add("Plot code", f"{function} documented", f"`{function}()`" in text, function)

    downloader = (root / "scripts" / "download_article71_sem_data.py").read_text()
    freezer = (root / "scripts" / "freeze_article71_sem.py").read_text()
    model_script = (root / "scripts" / "run_article71_sem_models.R").read_text()
    audit.add("Code", "downloader excludes publisher artwork", "franzosa-fig1-original" not in downloader, "download script")
    audit.add("Code", "freezer excludes publisher artwork", "franzosa-fig1-original" not in freezer, "freeze script")
    audit.add("Code", "analysis seed fixed", "SEED <- 71001L" in model_script, "R model script")
    audit.add("Code", "power repetitions fixed", "POWER_REPETITIONS <- 1000L" in model_script, "R model script")
    audit.add("Code", "positive control included", "Positive control" in model_script, "R model script")
    audit.add("Code", "medsens included", "medsens(" in model_script, "R model script")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    qa_dir = args.qa_dir.resolve()
    qa_dir.mkdir(parents=True, exist_ok=True)

    audit = Audit()
    audit_bundle(frozen, audit)
    audit_sources(frozen, audit)
    audit_cohorts(frozen, audit)
    audit_models(frozen, audit)
    audit_simulations(frozen, audit)
    audit_figures(root, frozen, audit)
    audit_article(root, frozen, audit)

    frame = audit.frame()
    frame.to_csv(qa_dir / "qa-report.tsv", sep="\t", index=False)
    failed = frame.loc[frame["Status"].eq("FAIL")]
    summary = {
        "article": 71,
        "checks": int(len(frame)),
        "passed": int(frame["Status"].eq("PASS").sum()),
        "failed": int(len(failed)),
        "status": "PASS" if failed.empty else "FAIL",
        "failures": [asdict(row) for row in audit.rows if not row.status],
    }
    (qa_dir / "qa-report.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if not failed.empty:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
