#!/usr/bin/env python3
"""Fail-closed validation for Article 63 paired multiomics evidence."""

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

from article42_44_validation_utils import (
    Audit,
    audit_checksums,
    audit_figures,
    finish,
    sha256,
)


FIGURES = (
    "63-paired-data-audit",
    "63-procrustes-concordance",
    "63-global-concordance",
    "63-halla-discovery-replication",
    "63-diablo-external-validation",
    "63-diablo-feature-stability",
)

SUMMARY_FILES = (
    "analysis-metrics.json",
    "diablo-selected-stability.tsv",
    "halla-branch-overlap.tsv",
    "halla-branch-summary.tsv",
    "halla-pair-validation.tsv.gz",
    "halla-replication-summary.tsv",
    "top-replicated-pairs.tsv",
)

EXPECTED_RESOURCES = {
    "genera.tsv": (
        18_101_016,
        "c4a541fe198a147beccd72d52fb2ebbf75a8cdf75cb3df75f823290971409d3f",
    ),
    "mtb.tsv": (
        12_485_385,
        "528b5e5953bd3697dd1ecf551d810d536c0679bd922e3fa3a6956c1412c6288c",
    ),
    "mtb.map.tsv": (
        695_256,
        "0dcdcce04a4e9b2b9b1632a410959baa4802ea9e14fc7c44f63bc17f699e5c65",
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
    parser.add_argument("--chapter", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stage-dir", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    return parser.parse_args()


def near(value: object, expected: float, tolerance: float = 1e-9) -> bool:
    try:
        return math.isclose(float(value), expected, abs_tol=tolerance, rel_tol=tolerance)
    except (TypeError, ValueError):
        return False


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def image_pixel_sha(path: Path) -> str:
    with Image.open(path) as image:
        normalized = image.convert("RGBA")
        payload = (
            normalized.size[0].to_bytes(8, "little")
            + normalized.size[1].to_bytes(8, "little")
            + normalized.tobytes()
        )
    return hashlib.sha256(payload).hexdigest()


def audit_chapter_63(chapter: Path, audit: Audit) -> None:
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
        "analysis-seed": "63001" in text,
        "plot-seed": "20260763" in text,
        "inline-theme": all(
            token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")
        ),
        "real-data": all(
            token in text
            for token in ("FRANZOSA_IBD_2019", "89a519d8", "155", "65", "220")
        ),
        "feature-contract": all(
            token in text
            for token in ("166", "153", "20%", "0.05%", "10^-6", "log1p")
        ),
        "global-contract": all(
            token in text for token in ("10 维", "9,999", "2,000", "diagnosis")
        ),
        "halla-contract": all(
            token in text
            for token in ("HAllA 0.8.40", "analytic Spearman", "FNR", "953")
        ),
        "diablo-contract": all(
            token in text
            for token in ("mixOmics 6.26.0", "5-fold", "100", "200", "0.635")
        ),
        "interpretation-boundary": all(
            token in text for token in ("production", "consumption", "causality")
        ),
        "methods-template": "**Paired microbiome–metabolome integration.**" in text,
        "results-template": "**Results template.**" in text,
        "frozen-input": "data/small/63-metabolomics-integration-frozen" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(
            token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')
        ),
        "no-source-theme": 'source("R/theme_pub.R")' not in text
        and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text)
        is None,
        "no-planning-comment": "Planned chapter" not in text
        and "Do not publish" not in text,
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
    for section in sections:
        checks[f"section-{section}"] = (
            re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
        )
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def audit_matrix(
    frozen: Path,
    name: str,
    expected_features: int,
    samples: list[str],
    audit: Audit,
    *,
    clr: bool = False,
) -> None:
    frame = read_tsv(frozen / name)
    values = frame.iloc[:, 1:].to_numpy(float)
    status = (
        frame.shape == (220, expected_features + 1)
        and frame.iloc[:, 0].astype(str).tolist() == samples
        and np.isfinite(values).all()
    )
    audit.add("Matrix", f"{name}-shape-alignment-finite", status, frame.shape)
    if clr:
        audit.add(
            "Matrix",
            f"{name}-clr-zero-sum",
            np.allclose(values.sum(axis=1), 0.0, atol=1e-10, rtol=0),
            float(np.abs(values.sum(axis=1)).max()),
        )


def stage_reanalysis(
    root: Path,
    frozen: Path,
    stage: Path,
    python: Path,
    audit: Audit,
) -> None:
    if stage.exists():
        shutil.rmtree(stage)
    summary = stage / "summary"
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
    summarized = subprocess.run(
        [
            str(python),
            str(frozen / "scripts/summarize_article63_multiomics.py"),
            "--input-dir",
            str(frozen),
            "--output-dir",
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
        observed = sha256(summary / name) if (summary / name).is_file() else "MISSING"
        expected = sha256(frozen / "summary" / name)
        audit.add("Reanalysis", f"summary-hash-{name}", observed == expected, observed)

    plotted = subprocess.run(
        [
            str(python),
            str(frozen / "scripts/plot_article63_multiomics.py"),
            "--input-dir",
            str(frozen),
            "--summary-dir",
            str(summary),
            "--figure-dir",
            str(figures),
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
    audit_figures(figures, audit, FIGURES)
    for stem in FIGURES:
        staged = figures / f"{stem}.png"
        published = root / "figures" / f"{stem}.png"
        status = staged.is_file() and published.is_file()
        observed = image_pixel_sha(staged) if staged.is_file() else "MISSING"
        expected = image_pixel_sha(published) if published.is_file() else "MISSING"
        audit.add(
            "Reanalysis",
            f"pixel-match-{stem}",
            status and observed == expected,
            observed,
        )


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    figures = args.figure_dir.resolve()
    audit = Audit()

    audit_checksums(frozen, audit)
    contract = json.loads((frozen / "frozen-contract.json").read_text(encoding="utf-8"))
    audit.add("Contract", "article", contract.get("article") == 63, contract)
    audit.add(
        "Contract",
        "dimensions",
        contract.get("samples") == 220
        and contract.get("independent_subjects") == 220
        and contract.get("discovery_samples") == 155
        and contract.get("external_validation_samples") == 65,
        contract,
    )
    audit.add(
        "Contract",
        "features",
        contract.get("selected_microbes") == 166
        and contract.get("selected_metabolites") == 153,
        contract,
    )
    audit.add("Contract", "seed", contract.get("seed") == 63001, contract)
    audit.add(
        "Contract",
        "resampling",
        contract.get("global_permutations") == 9999
        and contract.get("global_bootstraps") == 2000
        and contract.get("diablo_training_bootstraps") == 100
        and contract.get("diablo_label_permutations") == 200,
        contract,
    )
    audit.add(
        "Contract",
        "source-policy",
        contract.get("source_tables_included") is False
        and contract.get("source_manifest_and_checksums_included") is True
        and contract.get("transformed_analysis_matrices_included") is True,
        contract,
    )

    source = json.loads(
        (frozen / "source/download-manifest.json").read_text(encoding="utf-8")
    )
    audit.add(
        "Source",
        "identity",
        source.get("dataset") == "FRANZOSA_IBD_2019"
        and source.get("commit")
        == "89a519d8c832008fbc6e650453e83e2f04858d02"
        and source.get("paper_doi") == "10.1038/s41564-018-0306-4",
        source,
    )
    resources = {row["name"]: row for row in source["resources"]}
    for name, (size, digest) in EXPECTED_RESOURCES.items():
        row = resources.get(name, {})
        audit.add(
            "Source",
            f"locked-{name}",
            row.get("bytes") == size and row.get("sha256") == digest,
            row,
        )

    metadata = read_tsv(frozen / "sample-metadata.tsv")
    samples = metadata["Sample"].astype(str).tolist()
    observed_groups = (
        metadata.groupby(["Cohort", "Study.Group"], observed=True).size().to_dict()
    )
    expected_groups = {
        ("PRISM", "CD"): 68,
        ("PRISM", "UC"): 53,
        ("PRISM", "Control"): 34,
        ("Validation", "CD"): 20,
        ("Validation", "UC"): 23,
        ("Validation", "Control"): 22,
    }
    audit.add(
        "Input",
        "paired-independent-samples",
        len(metadata) == 220
        and metadata.Sample.nunique() == 220
        and metadata.Subject.astype(str).nunique() == 220,
        metadata.shape,
    )
    audit.add(
        "Input",
        "group-counts",
        observed_groups == expected_groups,
        {f"{cohort}|{group}": int(count) for (cohort, group), count in observed_groups.items()},
    )

    audit_matrix(frozen, "microbiome-relative.tsv.gz", 166, samples, audit)
    audit_matrix(frozen, "microbiome-clr.tsv.gz", 166, samples, audit, clr=True)
    audit_matrix(frozen, "metabolome-intensity.tsv.gz", 153, samples, audit)
    audit_matrix(frozen, "metabolome-log1p.tsv.gz", 153, samples, audit)
    feature = read_tsv(frozen / "feature-audit.tsv")
    audit.add(
        "Input",
        "feature-ledger",
        len(feature) == 319
        and feature.FeatureID.nunique() == 319
        and (feature.Modality == "Microbiome").sum() == 166
        and (feature.Modality == "Metabolome").sum() == 153,
        feature.shape,
    )

    global_tests = read_tsv(frozen / "global/global-concordance.tsv")
    prism = global_tests.loc[
        global_tests.Cohort.eq("PRISM")
        & global_tests.Restriction.eq("Within diagnosis")
    ].iloc[0]
    validation = global_tests.loc[
        global_tests.Cohort.eq("Validation")
        & global_tests.Restriction.eq("Within diagnosis")
    ].iloc[0]
    audit.add(
        "Global",
        "four-tests",
        len(global_tests) == 4
        and set(global_tests.Restriction) == {"Unrestricted", "Within diagnosis"}
        and set(global_tests.Permutations) == {9999}
        and set(global_tests.Bootstraps) == {2000},
        global_tests.shape,
    )
    audit.add(
        "Global",
        "prism-lock",
        near(prism.ProcrustesR, 0.619492546198694)
        and near(prism.MantelRho, 0.550298081490437)
        and near(prism.ProcrustesP, 0.0001)
        and near(prism.MantelP, 0.0001),
        prism.to_json(),
    )
    audit.add(
        "Global",
        "validation-lock",
        near(validation.ProcrustesR, 0.699745168346495)
        and near(validation.MantelRho, 0.674658006577757)
        and near(validation.ProcrustesP, 0.0001)
        and near(validation.MantelP, 0.0001),
        validation.to_json(),
    )

    halla_manifest = json.loads(
        (frozen / "halla-results/run-manifest.json").read_text(encoding="utf-8")
    )
    audit.add(
        "HAllA",
        "method-contract",
        halla_manifest.get("association") == "spearman"
        and halla_manifest.get("linkage") == "average"
        and near(halla_manifest.get("fdr_alpha"), 0.05)
        and near(halla_manifest.get("fnr_threshold"), 0.2)
        and "analytic Spearman" in halla_manifest.get("pvalue_mode", ""),
        halla_manifest,
    )
    for branch, expected_q, expected_blocks in (
        ("adjusted", 4416, 1970),
        ("raw", 7845, 2387),
    ):
        associations = read_tsv(
            frozen / f"halla-results/{branch}/all_associations.txt"
        )
        clusters = read_tsv(frozen / f"halla-results/{branch}/sig_clusters.txt")
        audit.add(
            "HAllA",
            f"{branch}-pair-universe",
            len(associations) == 25_398
            and not associations[["X_features", "Y_features"]].duplicated().any(),
            associations.shape,
        )
        audit.add(
            "HAllA",
            f"{branch}-discoveries",
            int((associations["q-values"] < 0.05).sum()) == expected_q
            and len(clusters) == expected_blocks,
            {"pairs": int((associations["q-values"] < 0.05).sum()), "blocks": len(clusters)},
        )

    replication = read_tsv(frozen / "summary/halla-replication-summary.tsv")
    replication_by_stage = dict(zip(replication.Stage, replication.Pairs))
    audit.add(
        "HAllA",
        "external-replication",
        replication_by_stage.get("PRISM discovery BH q < 0.05") == 4416
        and replication_by_stage.get("Same direction in Validation") == 3812
        and replication_by_stage.get("Validation BH q < 0.05") == 959
        and replication_by_stage.get("Same direction + Validation BH") == 953,
        replication_by_stage,
    )

    tuning = read_tsv(frozen / "diablo/tuning-summary.tsv")
    selected_keep = {
        (row.Block, int(row.Component)): int(row.KeepX)
        for row in tuning.itertuples(index=False)
    }
    audit.add(
        "DIABLO",
        "training-only-tuning",
        selected_keep
        == {
            ("microbiome", 1): 20,
            ("microbiome", 2): 20,
            ("metabolome", 1): 20,
            ("metabolome", 2): 5,
        }
        and set(tuning.Folds) == {5}
        and set(tuning.Repeats) == {5}
        and set(tuning.DesignWeight) == {0.1},
        str(selected_keep),
    )
    external = read_tsv(frozen / "diablo/external-metrics.tsv").set_index("Metric")
    null = read_tsv(frozen / "diablo/label-permutation-summary.tsv").iloc[0]
    confusion = read_tsv(frozen / "diablo/external-confusion.tsv")
    audit.add(
        "DIABLO",
        "external-metrics",
        near(external.loc["BalancedAccuracy", "Estimate"], 0.63498023715415)
        and near(external.loc["BalancedAccuracy", "Low"], 0.527865612648221)
        and near(external.loc["BalancedAccuracy", "High"], 0.740793807641634)
        and near(external.loc["MacroF1", "Estimate"], 0.625394088006252)
        and int(confusion.Samples.sum()) == 65,
        external.reset_index().to_json(orient="records"),
    )
    audit.add(
        "DIABLO",
        "label-null",
        int(null.Permutations) == 200
        and near(null.EmpiricalP, 0.00497512437810945),
        null.to_json(),
    )
    stable = read_tsv(frozen / "summary/diablo-selected-stability.tsv")
    audit.add(
        "DIABLO",
        "feature-stability",
        len(stable) == 65
        and int(stable.Stable70.astype(str).str.lower().eq("true").sum()) == 32
        and set(stable.Bootstraps) == {100},
        stable.shape,
    )

    metrics = json.loads(
        (frozen / "summary/analysis-metrics.json").read_text(encoding="utf-8")
    )
    audit.add(
        "Summary",
        "locked-counts",
        metrics.get("halla_adjusted_significant_pairs") == 4416
        and metrics.get("halla_raw_significant_pairs") == 7845
        and metrics.get("halla_replicated_pairs") == 953
        and metrics.get("diablo_stable_selected_features") == 32,
        metrics,
    )

    audit_figures(figures, audit, FIGURES)
    plot_source = (root / "scripts/plot_article63_multiomics.py").read_text(
        encoding="utf-8"
    )
    audit.add(
        "Figure",
        "english-only-source",
        re.search(r"[\u4e00-\u9fff]", plot_source) is None,
        "plot source contains no CJK text",
    )
    audit_chapter_63(args.chapter.resolve(), audit)
    stage_reanalysis(
        root,
        frozen,
        args.stage_dir.resolve(),
        args.python.resolve(),
        audit,
    )

    return finish(
        article=63,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "samples": 220,
            "discovery_samples": 155,
            "external_validation_samples": 65,
            "figures": len(FIGURES),
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
