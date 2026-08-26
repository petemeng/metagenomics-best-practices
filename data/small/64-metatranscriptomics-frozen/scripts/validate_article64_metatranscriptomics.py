#!/usr/bin/env python3
"""Fail-closed validation for Article 64 paired metatranscriptomics evidence."""

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
    "64-pairing-audit",
    "64-dna-rna-concordance",
    "64-sample-concordance",
    "64-relative-activity",
    "64-diagnosis-audit",
    "64-ratio-sensitivity",
)

EXPECTED_RESOURCES = {
    "hmp2-metadata.csv": (9_074_342, "656b7bd97660ddb875548805e30bede31f2d1208293f7170d2d5755e33862ec9"),
    "mgx-pathabundance-rela.tsv.gz": (9_037_246, "e840fb86ed8049bc697a20a3904da11d29878ed02f43f50764a28c93d2111216"),
    "mtx-pathabundance-rela.tsv.gz": (1_715_915, "ee2a1afb69b66bbdac014b48db6a692dc09146d539e4a93d3fea0c2d9903ac08"),
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


def image_pixel_sha(path: Path) -> str:
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
        "analysis-seed": "64001" in text,
        "plot-seed": "20260764" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "real-data": all(token in text for token in ("IBDMDB", "Lloyd-Price", "10.1038/s41586-019-1237-9")),
        "sample-contract": all(token in text for token in ("734", "19", "715", "711", "104")),
        "feature-contract": all(token in text for token in ("418", "211", "20%", "co-detection")),
        "results-contract": all(token in text for token in ("0.793", "173", "202", "0.342")),
        "interpretation-boundary": all(token in text for token in ("absolute transcription", "per-cell", "flux", "causality")),
        "methods-template": "**Metagenome–metatranscriptome integration.**" in text,
        "results-template": "**Results template.**" in text,
        "frozen-input": "data/small/64-metatranscriptomics-frozen" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
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
        checks[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def audit_matrix(path: Path, rows: int, columns: int, first_column: str, audit: Audit) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    audit.add("Matrix", f"{path.name}-shape", frame.shape == (rows, columns), frame.shape)
    audit.add("Matrix", f"{path.name}-first-column", frame.columns[0] == first_column, frame.columns[0])
    values = frame.iloc[:, 1:].to_numpy(float)
    audit.add("Matrix", f"{path.name}-finite-or-na", not np.isinf(values).any(), "no infinity")
    return frame


def stage_figures(root: Path, frozen: Path, stage: Path, python: Path, final_figures: Path, audit: Audit) -> None:
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
            str(frozen / "scripts/plot_article64_metatranscriptomics.py"),
            "--input-dir",
            str(frozen),
            "--figure-dir",
            str(figures),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
    )
    audit.add("Reanalysis", "plot-script-exit", result.returncode == 0, result.stdout + result.stderr)
    for stem in FIGURES:
        staged = figures / f"{stem}.png"
        final = final_figures / f"{stem}.png"
        status = staged.is_file() and final.is_file() and image_pixel_sha(staged) == image_pixel_sha(final)
        audit.add("Reanalysis", f"{stem}-pixel-identical", status, image_pixel_sha(staged) if staged.is_file() else "MISSING")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_checksums(frozen, audit)

    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    expected_metrics = {
        "article": 64,
        "source_mgx_profiles": 1638,
        "source_mtx_profiles": 735,
        "raw_paired_columns": 734,
        "technical_replicates_excluded": 19,
        "unique_paired_biospecimens": 715,
        "zero_layer_samples_excluded": 4,
        "analysis_samples": 711,
        "independent_subjects": 104,
        "shared_unstratified_pathways": 418,
        "selected_pathways": 211,
        "relative_activity_bh05": 173,
        "relative_activity_bh05_positive": 80,
        "relative_activity_bh05_negative": 93,
        "diagnosis_tests": 202,
        "diagnosis_bh05": 0,
    }
    for key, expected in expected_metrics.items():
        audit.add("Metric", key, metrics.get(key) == expected, metrics.get(key))
    audit.add("Metric", "median-subject-spearman", near(metrics.get("median_subject_spearman"), 0.7931967768969652), metrics.get("median_subject_spearman"))
    audit.add("Metric", "minimum-diagnosis-q", near(metrics.get("minimum_diagnosis_q"), 0.3416596554463312), metrics.get("minimum_diagnosis_q"))

    manifest = json.loads((frozen / "source-manifest.json").read_text(encoding="utf-8"))
    resources = manifest.get("resources", {})
    for name, (size, digest) in EXPECTED_RESOURCES.items():
        record = resources.get(name, {})
        audit.add("Source", f"{name}-bytes", record.get("bytes") == size, record.get("bytes"))
        audit.add("Source", f"{name}-sha256", record.get("sha256") == digest, record.get("sha256"))
        audit.add("Source", f"{name}-official-url", str(record.get("url", "")).startswith("https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/"), record.get("url"))

    metadata = pd.read_csv(frozen / "sample-metadata.tsv", sep="\t")
    audit.add("Sample", "rows-unique-samples", len(metadata) == metadata["SampleID"].nunique() == 711, len(metadata))
    audit.add("Sample", "subjects", metadata["SubjectID"].nunique() == 104, metadata["SubjectID"].nunique())
    audit.add("Sample", "no-technical-replicates", not metadata["SampleID"].str.endswith("_TR").any(), "no _TR")
    audit.add("Sample", "positive-layer-sums", metadata[["DNAPathwaySum", "RNAPathwaySum"]].gt(0).all().all(), metadata[["DNAPathwaySum", "RNAPathwaySum"]].min().to_dict())
    subject_counts = metadata.drop_duplicates("SubjectID")["Diagnosis"].value_counts().to_dict()
    audit.add("Sample", "diagnosis-subject-counts", subject_counts == {"CD": 50, "UC": 28, "Control": 26}, subject_counts)

    feature = pd.read_csv(frozen / "feature-audit.tsv", sep="\t")
    audit.add("Feature", "shared-pathways", len(feature) == 418, len(feature))
    audit.add("Feature", "selected-pathways", int(feature["Selected"].sum()) == 211, int(feature["Selected"].sum()))
    audit.add("Feature", "selected-gate", feature.loc[feature["Selected"], "CoDetection"].ge(0.20).all(), feature.loc[feature["Selected"], "CoDetection"].min())
    audit.add("Feature", "unselected-below-gate", feature.loc[~feature["Selected"], "CoDetection"].lt(0.20).all(), feature.loc[~feature["Selected"], "CoDetection"].max())

    dna = audit_matrix(frozen / "dna-relative.tsv.gz", 711, 212, "SampleID", audit)
    rna = audit_matrix(frozen / "rna-relative.tsv.gz", 711, 212, "SampleID", audit)
    activity = audit_matrix(frozen / "activity-log2-rna-dna.tsv.gz", 711, 212, "SampleID", audit)
    subject_activity = audit_matrix(frozen / "subject-activity.tsv.gz", 104, 212, "SubjectID", audit)
    audit.add("Matrix", "sample-order-dna-rna", dna["SampleID"].tolist() == rna["SampleID"].tolist() == activity["SampleID"].tolist(), "aligned")
    audit.add("Matrix", "activity-has-defined-values", activity.iloc[:, 1:].notna().sum().min() >= 142, int(activity.iloc[:, 1:].notna().sum().min()))
    audit.add("Matrix", "subject-activity-rows", subject_activity["SubjectID"].nunique() == 104, subject_activity["SubjectID"].nunique())

    activity_results = pd.read_csv(frozen / "activity-results.tsv", sep="\t")
    diagnosis = pd.read_csv(frozen / "diagnosis-results.tsv", sep="\t")
    sensitivity = pd.read_csv(frozen / "sensitivity-summary.tsv", sep="\t")
    audit.add("Statistics", "activity-tests", len(activity_results) == 211, len(activity_results))
    audit.add("Statistics", "activity-bh05", int(activity_results["BH05"].sum()) == 173, int(activity_results["BH05"].sum()))
    audit.add("Statistics", "diagnosis-tests", int(diagnosis["PValue"].notna().sum()) == 202, int(diagnosis["PValue"].notna().sum()))
    audit.add("Statistics", "diagnosis-no-bh05", int(diagnosis["BH05"].sum()) == 0, int(diagnosis["BH05"].sum()))
    pseudo = sensitivity.loc[sensitivity["Analysis"].eq("Pseudocount"), "RankSpearman"]
    audit.add("Statistics", "pseudocount-rank-sensitivity", pseudo.between(0.70, 0.80).all(), pseudo.tolist())

    audit_chapter(args.chapter.resolve(), audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    stage_figures(root, frozen, args.stage_dir.resolve(), args.python.resolve(), args.figure_dir.resolve(), audit)

    raise SystemExit(
        finish(
            article=64,
            audit=audit,
            output=args.output_dir.resolve(),
            payload={
                "analysis_samples": metrics["analysis_samples"],
                "independent_subjects": metrics["independent_subjects"],
                "selected_pathways": metrics["selected_pathways"],
                "relative_activity_bh05": metrics["relative_activity_bh05"],
                "diagnosis_bh05": metrics["diagnosis_bh05"],
            },
        )
    )


if __name__ == "__main__":
    main()
