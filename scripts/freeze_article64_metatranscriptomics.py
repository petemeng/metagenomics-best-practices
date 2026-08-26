#!/usr/bin/env python3
"""Create a checksum-covered Article 64 metatranscriptomics evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "activity-log2-rna-dna.tsv.gz",
    "activity-results.tsv",
    "analysis-metrics.json",
    "concordance-summary.tsv",
    "diagnosis-results.tsv",
    "dna-relative.tsv.gz",
    "exclusion-ledger.json",
    "feature-audit.tsv",
    "rna-relative.tsv.gz",
    "sample-attrition.tsv",
    "sample-concordance.tsv",
    "sample-metadata.tsv",
    "sensitivity-summary.tsv",
    "software-versions.json",
    "source-manifest.json",
    "subject-activity.tsv.gz",
    "subject-concordance.tsv",
)

SCRIPT_FILES = (
    "download_article64_metatranscriptomics_data.py",
    "prepare_article64_metatranscriptomics.py",
    "plot_article64_metatranscriptomics.py",
    "freeze_article64_metatranscriptomics.py",
    "validate_article64_metatranscriptomics.py",
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


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    target = args.output_dir.resolve()
    metrics = json.loads((work / "analysis-metrics.json").read_text(encoding="utf-8"))
    if metrics.get("article") != 64:
        raise RuntimeError("Article identity mismatch")

    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    for name in WORK_FILES:
        copy(work / name, target / name)
    for name in SCRIPT_FILES:
        copy(root / "scripts" / name, target / "scripts" / name)
    for name in ENV_FILES:
        copy(root / "env" / name, target / "env" / name)

    contract = {
        "article": 64,
        "dataset": "IBDMDB HMP2 MGX/MTX",
        "doi": "10.1038/s41586-019-1237-9",
        "source_tables_included": False,
        "source_manifest_and_checksums_included": True,
        "transformed_matrices_included": True,
        "analysis_samples": metrics["analysis_samples"],
        "independent_subjects": metrics["independent_subjects"],
        "selected_pathways": metrics["selected_pathways"],
        "technical_replicates_excluded": metrics["technical_replicates_excluded"],
        "zero_layer_samples_excluded": metrics["zero_layer_samples_excluded"],
        "seed": metrics["seed"],
        "plot_seed": metrics["plot_seed"],
        "interpretation": "relative DNA/RNA allocation; not absolute transcription, per-cell regulation, flux, or causality",
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    notice = """Article 64 frozen evidence bundle

The source data are the official IBDMDB HMP2 merged HUMAnN pathway-relative-
abundance tables for metagenomes and metatranscriptomes plus the 2018-08-20
metadata release. Their exact URLs, byte counts, and SHA-256 checksums are
retained in source-manifest.json; the approximately 20 MB of original source
tables are fetched by the included downloader rather than duplicated here.

The bundle includes the exact 711-sample paired transformed matrices, the
104-subject ledgers, all pathway-level tests and sensitivity results, scripts,
and environment records. Paired technical-replicate columns and four samples
with a zero pathway layer were excluded before analysis. DNA/RNA ratios are
relative allocation indices and do not estimate absolute transcription rates,
per-cell regulation, pathway flux, or causality.
"""
    (target / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    payloads = sorted(
        path for path in target.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    (target / "file-checksums.sha256").write_text(
        "\n".join(f"{sha256(path)}  {path.relative_to(target).as_posix()}" for path in payloads) + "\n",
        encoding="utf-8",
    )
    print(f"Frozen Article 64 bundle: {len(payloads)} payload files in {target}")


if __name__ == "__main__":
    main()
