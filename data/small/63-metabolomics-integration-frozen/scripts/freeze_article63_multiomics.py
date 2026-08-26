#!/usr/bin/env python3
"""Create a checksum-covered Article 63 paired-multiomics evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


ROOT_FILES = (
    "analysis-contract.json",
    "covariate-design.tsv",
    "feature-attrition.tsv",
    "feature-audit.tsv",
    "metabolome-intensity.tsv.gz",
    "metabolome-log1p.tsv.gz",
    "microbiome-clr.tsv.gz",
    "microbiome-relative.tsv.gz",
    "prepare-software-versions.json",
    "sample-metadata.tsv",
)

GLOBAL_FILES = (
    "global-concordance.tsv",
    "ordination-variance.tsv",
    "prism-ordination-models.rds",
    "procrustes-scores.tsv",
    "session-info.txt",
)

HALLA_INPUT_FILES = (
    "prism-metabolome-adjusted.tsv",
    "prism-metabolome-raw.tsv",
    "prism-microbiome-adjusted.tsv",
    "prism-microbiome-raw.tsv",
    "validation-metabolome-adjusted.tsv",
    "validation-metabolome-raw.tsv",
    "validation-microbiome-adjusted.tsv",
    "validation-microbiome-raw.tsv",
)

HALLA_ROOT_FILES = (
    "adjusted-halla.log",
    "raw-halla.log",
    "run-manifest.json",
)

HALLA_BRANCH_FILES = (
    "X_linkage.npy",
    "Y_linkage.npy",
    "all_associations.txt",
    "performance.txt",
    "sig_clusters.txt",
)

DIABLO_FILES = (
    "bootstrap-checkpoint.rds",
    "bootstrap-feature-stability.tsv",
    "external-class-metrics.tsv",
    "external-confusion.tsv",
    "external-metrics.tsv",
    "external-predictions.tsv",
    "final-model.rds",
    "final-selected-features.tsv",
    "label-permutation-null.tsv",
    "label-permutation-summary.tsv",
    "latent-correlations.tsv",
    "latent-scores.tsv",
    "permutation-checkpoint.rds",
    "session-info.txt",
    "tuning-object.rds",
    "tuning-summary.tsv",
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

SCRIPT_FILES = (
    "download_article63_multiomics_data.py",
    "prepare_article63_multiomics.py",
    "run_article63_global.R",
    "run_article63_halla.py",
    "run_article63_diablo.R",
    "summarize_article63_multiomics.py",
    "plot_article63_multiomics.py",
    "freeze_article63_multiomics.py",
    "validate_article63_multiomics.py",
    "article42_44_validation_utils.py",
)

ENV_FILES = (
    "multiomics-python.yml",
    "multiomics-r-packages.tsv",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy(source: Path, target: Path) -> None:
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def copy_group(source: Path, target: Path, names: tuple[str, ...]) -> None:
    for name in names:
        copy(source / name, target / name)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    cache = args.cache_dir.resolve()
    target = args.output_dir.resolve()

    metrics = json.loads(
        (work / "summary/analysis-metrics.json").read_text(encoding="utf-8")
    )
    source_manifest = json.loads(
        (cache / "download-manifest.json").read_text(encoding="utf-8")
    )
    if metrics.get("article") != 63 or source_manifest.get("article") != 63:
        raise RuntimeError("Article identity mismatch")

    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)

    copy_group(work, target, ROOT_FILES)
    copy_group(work / "global", target / "global", GLOBAL_FILES)
    copy_group(work / "halla", target / "halla", HALLA_INPUT_FILES)
    copy_group(work / "halla-results", target / "halla-results", HALLA_ROOT_FILES)
    for branch in ("adjusted", "raw"):
        copy_group(
            work / "halla-results" / branch,
            target / "halla-results" / branch,
            HALLA_BRANCH_FILES,
        )
    copy_group(work / "diablo", target / "diablo", DIABLO_FILES)
    copy_group(work / "summary", target / "summary", SUMMARY_FILES)
    copy(cache / "download-manifest.json", target / "source/download-manifest.json")
    for name in SCRIPT_FILES:
        copy(root / "scripts" / name, target / "scripts" / name)
    for name in ENV_FILES:
        copy(root / "env" / name, target / "env" / name)

    contract = {
        "article": 63,
        "created_from": str(work.relative_to(root)),
        "dataset": "FRANZOSA_IBD_2019",
        "source_commit": source_manifest["commit"],
        "source_tables_included": False,
        "source_manifest_and_checksums_included": True,
        "transformed_analysis_matrices_included": True,
        "samples": metrics["samples"],
        "independent_subjects": metrics["independent_subjects"],
        "discovery_samples": metrics["discovery_samples"],
        "external_validation_samples": metrics["external_validation_samples"],
        "selected_microbes": metrics["selected_microbes"],
        "selected_metabolites": metrics["selected_metabolites"],
        "seed": 63001,
        "global_permutations": 9999,
        "global_bootstraps": 2000,
        "diablo_training_bootstraps": 100,
        "diablo_label_permutations": 200,
        "diablo_metric_bootstraps": 2000,
        "interpretation": "paired association and external prediction; no production, consumption, or causal claim",
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    notice = """Article 63 frozen evidence bundle

The paired inputs derive from the Franzosa et al. inflammatory-bowel-disease
microbiome/metabolome study as standardized in the Muller et al. curated
collection. The source commit, URLs, byte counts, and SHA-256 checksums are
retained, while the four approximately 31-MB source tables are fetched by the
included downloader rather than duplicated. The exact transformed matrices,
feature and sample ledgers, HAllA inputs/results, global tests, DIABLO model
objects/results, summaries, environments, and regeneration scripts are stored.

All 155 PRISM subjects form the discovery/training cohort. The 65-subject
Validation cohort is projected through PRISM-fitted ordination axes, re-tests
prespecified HAllA pairs, and receives one final DIABLO evaluation; it is never
used for feature filtering or tuning. Microbe-metabolite association and
cross-block prediction do not establish microbial production, consumption,
metabolic handoff, directionality, or causality.
"""
    (target / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    payloads = sorted(
        path
        for path in target.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    lines = [
        f"{sha256(path)}  {path.relative_to(target).as_posix()}" for path in payloads
    ]
    (target / "file-checksums.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(f"Frozen Article 63 bundle: {len(payloads)} payload files in {target}")


if __name__ == "__main__":
    main()
