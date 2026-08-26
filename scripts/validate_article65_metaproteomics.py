#!/usr/bin/env python3
"""Fail-closed validation for Article 65 metaproteomic evidence."""

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
    "65-evidence-threshold-audit",
    "65-protein-namespace-audit",
    "65-multiomic-sample-alignment",
    "65-three-layer-concordance",
    "65-ec-cross-subject-correlation",
    "65-threshold-richness-sensitivity",
)

EXPECTED_RESOURCES = {
    "hmp2-metadata.csv": (9_074_342, "656b7bd97660ddb875548805e30bede31f2d1208293f7170d2d5755e33862ec9"),
    "mgx-ecs-rela.tsv.gz": (88_776_917, "23f81ec1b6a995cfed83224816dc89002f7f7a2880afd796233b1f5386fab220"),
    "mtx-ecs-rela.tsv.gz": (21_323_580, "0a5da078c521bd1f7c62d9583b9923cd05f686f344f2fdad75e234134657b3a0"),
    "mpx-1pep-1pct.tsv.gz": (1_097_492, "b4ceee58817059b77a2b010061420f39ce30e1598bd831cdef5cdd771d32dc02"),
    "mpx-1pep-5pct.tsv.gz": (1_710_498, "254cde477a5f95bbd5cf62de0364b0442959fc4e343a8903585f728729cbd111"),
    "mpx-2pep-1pct.tsv.gz": (424_364, "a001f9dfe9417f1d85004a43ec97199a16ad7e55862459f36da769da8632b12b"),
    "mpx-2pep-5pct.tsv.gz": (480_606, "b064d1b6ef3ee49914a5375a656f1d239634254b3694f02da87ead76af925c26"),
    "mpx-ecs.tsv.gz": (91_574, "6328db9bcb40dbe0e927747e7df3dcd0d9974d1fe5e35a9631f3620ae3565a7a"),
    "mpx-kos.tsv.gz": (146_771, "e3af45100537405a54b574ec5c8994d837d2f353f629970365d0cd3fe1d42e39"),
    "lloyd-price-supp-fig1.pdf": (4_278_512, "78d9845d62be38019e3f58daf2a9b085eda5dbb5a2a1c40491130bcd33e80760"),
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
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "bibliography": "bibliography: ../references.bib" in text,
        "analysis-seed": "65001" in text,
        "plot-seed": "20260765" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "hardware-contract": all(token in text for token in ("RAM", "CPU", "磁盘", "耗时")),
        "real-data": all(token in text for token in ("IBDMDB", "Lloyd-Price", "10.1038/s41586-019-1237-9")),
        "sample-contract": all(token in text for token in ("450", "281", "187", "186", "76", "154")),
        "feature-contract": all(token in text for token in ("864", "263", "10%", "789")),
        "threshold-contract": all(token in text for token in ("10,469", "66,606", "373", "1,413.5")),
        "results-contract": all(token in text for token in ("0.731", "0.109", "0.238", "0.476", "0.152", "0.097")),
        "fdr-contract": all(token in text for token in ("203", "66", "45", "5.4\\times10^{-6}")),
        "evidence-chain": all(token in text for token in ("PSM", "shared peptide", "protein group", "target-decoy")),
        "software-contract": all(token in text for token in ("MSGF+ v10072", "DIAMOND v0.8.22.84", "FragPipe 24.0", "MSFragger 4.4")),
        "interpretation-boundary": all(token in text for token in ("absolute protein", "translation", "enzyme activity", "flux", "causality")),
        "methods-template": "**Metaproteomic processing and multi-omic alignment.**" in text,
        "results-template": "**Results template.**" in text,
        "frozen-input": "data/small/65-metaproteomics-frozen" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "anchor": "../figures/65-lloyd-price-supp-fig1-page1-original.png" in text,
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
    audit.add("Matrix", f"{path.name}-finite", np.isfinite(values).all(), "all finite")
    audit.add("Matrix", f"{path.name}-nonnegative", (values >= 0).all(), float(values.min()))
    audit.add("Matrix", f"{path.name}-positive-row-sums", (values.sum(axis=1) > 0).all(), float(values.sum(axis=1).min()))
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
            str(frozen / "scripts/plot_article65_metaproteomics.py"),
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
        final = final_figures / f"{stem}.png"
        status = staged.is_file() and final.is_file() and image_pixel_sha(staged) == image_pixel_sha(final)
        audit.add(
            "Reanalysis",
            f"{stem}-pixel-identical",
            status,
            image_pixel_sha(staged) if staged.is_file() else "MISSING",
        )


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_checksums(frozen, audit)

    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    expected_metrics = {
        "article": 65,
        "source_mgx_profiles": 1638,
        "source_mtx_profiles": 735,
        "source_mpx_profiles": 450,
        "source_mgx_unstratified_ecs": 2173,
        "source_mtx_unstratified_ecs": 2082,
        "source_mpx_ec_rows": 910,
        "shared_ecs": 864,
        "selected_ecs": 263,
        "raw_triple_samples": 187,
        "technical_replicates_excluded": 0,
        "zero_layer_samples_excluded": 1,
        "analysis_samples": 186,
        "independent_subjects": 76,
        "raw_four_layer_samples": 155,
        "complete_four_layer_samples": 154,
        "ec_correlation_tests": 789,
        "ec_correlation_bh05_total": 314,
    }
    for key, expected in expected_metrics.items():
        audit.add("Metric", key, metrics.get(key) == expected, metrics.get(key))
    expected_pairs = {
        "DNA–RNA": (0.7307673229121173, 0.4756254272043746, 203),
        "DNA–Protein": (0.10927148902866204, 0.15245264118302965, 66),
        "RNA–Protein": (0.23774152898764994, 0.09655651472072058, 45),
    }
    for pair, (sample_median, ec_median, hits) in expected_pairs.items():
        record = metrics.get("pair_metrics", {}).get(pair, {})
        audit.add("Metric", f"{pair}-subject-sample-median", near(record.get("median_subject_sample_spearman"), sample_median), record)
        audit.add("Metric", f"{pair}-ec-median", near(record.get("median_ec_cross_subject_spearman"), ec_median), record)
        audit.add("Metric", f"{pair}-bh05", record.get("ec_bh05") == hits, record)

    manifest = json.loads((frozen / "source-manifest.json").read_text(encoding="utf-8"))
    resources = manifest.get("resources", {})
    for name, (size, digest) in EXPECTED_RESOURCES.items():
        record = resources.get(name, {})
        audit.add("Source", f"{name}-bytes", record.get("bytes") == size, record.get("bytes"))
        audit.add("Source", f"{name}-sha256", record.get("sha256") == digest, record.get("sha256"))
        url = str(record.get("url", ""))
        official = (
            url.startswith("https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/")
            if name != "lloyd-price-supp-fig1.pdf"
            else url.startswith("https://static-content.springer.com/")
        )
        audit.add("Source", f"{name}-official-url", official, url)

    threshold = pd.read_csv(frozen / "threshold-audit.tsv", sep="\t").set_index("Threshold")
    expected_thresholds = {
        "1 peptide · 1% FDR": (35045, 449, 31876, 1076.0, 1310576.0),
        "1 peptide · 5% FDR": (74721, 450, 66606, 1413.5, 1589978.0),
        "2 peptides · 1% FDR": (11284, 450, 10469, 373.0, 926768.0),
        "2 peptides · 5% FDR": (12582, 450, 11666, 447.0, 1083447.0),
    }
    audit.add("Threshold", "four-gates", set(threshold.index) == set(expected_thresholds), list(threshold.index))
    for label, expected in expected_thresholds.items():
        observed = threshold.loc[label]
        values = (
            int(observed["ProteinTableRows"]),
            int(observed["Profiles"]),
            int(observed["DetectedProteinIDs"]),
            float(observed["MedianDetectedPerSample"]),
            float(observed["TotalReportedCounts"]),
        )
        audit.add("Threshold", label, values == expected, values)
    richness = pd.read_csv(frozen / "richness-correlations.tsv", sep="\t")
    audit.add("Threshold", "six-richness-correlations", len(richness) == 6, len(richness))
    audit.add("Threshold", "rank-stability-range", richness["SpearmanRho"].between(0.975, 0.999).all(), richness["SpearmanRho"].tolist())
    overlap = pd.read_csv(frozen / "protein-id-overlap.tsv", sep="\t")
    audit.add("Threshold", "nonnested-id-universes", (overlap["OnlyA"].gt(0) & overlap["OnlyB"].gt(0)).any(), overlap[["OnlyA", "OnlyB"]].to_dict("records"))

    namespace = pd.read_csv(frozen / "protein-namespace-audit.tsv", sep="\t").set_index("Namespace")
    expected_namespace = {
        "Generic accession": (10515, 10461, 884797.0),
        "Taxon-prefixed reference": (523, 0, 0.0),
        "Explicit host namespace": (236, 0, 0.0),
        "Contaminant namespace": (10, 8, 41971.0),
    }
    for label, expected in expected_namespace.items():
        observed = namespace.loc[label]
        values = (int(observed["ProteinIDs"]), int(observed["DetectedProteinIDs"]), float(observed["ReportedCounts"]))
        audit.add("Namespace", label, values == expected, values)
    audit.add("Namespace", "count-shares-sum-one", near(namespace["ReportedCountShare"].sum(), 1.0), namespace["ReportedCountShare"].sum())

    metadata = pd.read_csv(frozen / "sample-metadata.tsv", sep="\t")
    audit.add("Sample", "rows-unique-samples", len(metadata) == metadata["SampleID"].nunique() == 186, len(metadata))
    audit.add("Sample", "subjects", metadata["SubjectID"].nunique() == 76, metadata["SubjectID"].nunique())
    audit.add("Sample", "positive-layer-sums", metadata[["DNAECSum", "RNAECSum", "ProteinECSum"]].gt(0).all().all(), metadata[["DNAECSum", "RNAECSum", "ProteinECSum"]].min().to_dict())
    mbx = metadata["MetabolomicsAvailable"].astype(str).str.lower().isin(["true", "1", "t", "yes"])
    audit.add("Sample", "mbx-exact-id-availability", int(mbx.sum()) == 154, int(mbx.sum()))
    subject_counts = metadata.drop_duplicates("SubjectID")["Diagnosis"].value_counts().to_dict()
    audit.add("Sample", "diagnosis-subject-counts", subject_counts == {"CD": 33, "Control": 22, "UC": 21}, subject_counts)
    attrition = pd.read_csv(frozen / "sample-attrition.tsv", sep="\t").set_index("Stage")["Count"].to_dict()
    expected_attrition = {
        "MPX profiles": 450,
        "Exact MPX + MGX": 281,
        "Exact MPX + MGX + MTX": 187,
        "Three-layer EC-complete": 186,
        "Triple profiles with MBX product": 155,
        "Four-layer & EC-complete": 154,
    }
    audit.add("Sample", "attrition-contract", attrition == expected_attrition, attrition)

    feature = pd.read_csv(frozen / "ec-feature-audit.tsv", sep="\t")
    selected_flag = feature["Selected"].astype(str).str.lower().isin(["true", "1", "t", "yes"])
    audit.add("Feature", "shared-ecs", len(feature) == 864, len(feature))
    audit.add("Feature", "selected-ecs", int(selected_flag.sum()) == 263, int(selected_flag.sum()))
    audit.add("Feature", "selected-gate", feature.loc[selected_flag, "MinimumPrevalence"].ge(0.10).all(), feature.loc[selected_flag, "MinimumPrevalence"].min())
    audit.add("Feature", "unselected-below-gate", feature.loc[~selected_flag, "MinimumPrevalence"].lt(0.10).all(), feature.loc[~selected_flag, "MinimumPrevalence"].max())

    dna = audit_matrix(frozen / "dna-ec-relative.tsv.gz", 186, 264, "SampleID", audit)
    rna = audit_matrix(frozen / "rna-ec-relative.tsv.gz", 186, 264, "SampleID", audit)
    protein = audit_matrix(frozen / "protein-ec-relative.tsv.gz", 186, 264, "SampleID", audit)
    subject_dna = audit_matrix(frozen / "subject-dna-ec-relative.tsv.gz", 76, 264, "SubjectID", audit)
    subject_rna = audit_matrix(frozen / "subject-rna-ec-relative.tsv.gz", 76, 264, "SubjectID", audit)
    subject_protein = audit_matrix(frozen / "subject-protein-ec-relative.tsv.gz", 76, 264, "SubjectID", audit)
    audit.add("Matrix", "sample-order-three-layers", dna["SampleID"].tolist() == rna["SampleID"].tolist() == protein["SampleID"].tolist(), "aligned")
    audit.add("Matrix", "subject-order-three-layers", subject_dna["SubjectID"].tolist() == subject_rna["SubjectID"].tolist() == subject_protein["SubjectID"].tolist(), "aligned")
    audit.add("Matrix", "selected-feature-order", dna.columns[1:].tolist() == rna.columns[1:].tolist() == protein.columns[1:].tolist(), "aligned")

    summary = pd.read_csv(frozen / "concordance-summary.tsv", sep="\t").set_index("LayerPair")
    for pair, (expected, _, _) in expected_pairs.items():
        audit.add("Statistics", f"{pair}-summary", near(summary.loc[pair, "MedianSpearmanRho"], expected), summary.loc[pair].to_dict())
        audit.add("Statistics", f"{pair}-ci-brackets-estimate", summary.loc[pair, "CILow"] <= expected <= summary.loc[pair, "CIHigh"], summary.loc[pair].to_dict())
    correlations = pd.read_csv(frozen / "ec-correlations.tsv", sep="\t")
    bh05 = correlations["BH05"].astype(str).str.lower().isin(["true", "1", "t", "yes"])
    audit.add("Statistics", "789-tests", len(correlations) == 789, len(correlations))
    audit.add("Statistics", "314-bh05", int(bh05.sum()) == 314, int(bh05.sum()))
    pair_hits = correlations.assign(BH05Parsed=bh05).groupby("LayerPair")["BH05Parsed"].sum().to_dict()
    audit.add("Statistics", "pair-bh05", pair_hits == {"DNA–Protein": 66, "DNA–RNA": 203, "RNA–Protein": 45}, pair_hits)
    top = correlations.loc[correlations["LayerPair"].eq("RNA–Protein")].nlargest(1, "SpearmanRho").iloc[0]
    audit.add("Statistics", "top-rna-protein-ec", top["EC"] == "2.7.1.45" and near(top["SpearmanRho"], 0.5320194562591866), top.to_dict())

    audit_chapter(args.chapter.resolve(), audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    anchor = args.figure_dir.resolve() / "65-lloyd-price-supp-fig1-page1-original.png"
    audit.add("Figure", "anchor-exists", anchor.is_file() and anchor.stat().st_size > 0, str(anchor))
    audit.add("Figure", "anchor-sha256", anchor.is_file() and hashlib.sha256(anchor.read_bytes()).hexdigest() == "d35f590c1f22d734e8c828fbcf1d6a98f1f0714465171f2bab89db5343b37ca7", str(anchor))
    stage_figures(root, frozen, args.stage_dir.resolve(), args.python.resolve(), args.figure_dir.resolve(), audit)

    raise SystemExit(
        finish(
            article=65,
            audit=audit,
            output=args.output_dir.resolve(),
            payload={
                "analysis_samples": metrics["analysis_samples"],
                "independent_subjects": metrics["independent_subjects"],
                "selected_ecs": metrics["selected_ecs"],
                "complete_four_layer_samples": metrics["complete_four_layer_samples"],
                "ec_correlation_bh05_total": metrics["ec_correlation_bh05_total"],
            },
        )
    )


if __name__ == "__main__":
    main()
