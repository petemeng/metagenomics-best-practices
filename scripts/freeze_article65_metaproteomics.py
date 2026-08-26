#!/usr/bin/env python3
"""Create a checksum-covered Article 65 metaproteomics evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "analysis-metrics.json",
    "concordance-summary.tsv",
    "dna-ec-relative.tsv.gz",
    "ec-correlations.tsv",
    "ec-feature-audit.tsv",
    "exclusion-ledger.json",
    "mgx-ec-unstratified.tsv.gz",
    "mtx-ec-unstratified.tsv.gz",
    "protein-ec-relative.tsv.gz",
    "protein-id-overlap.tsv",
    "protein-namespace-audit.tsv",
    "richness-correlations.tsv",
    "rna-ec-relative.tsv.gz",
    "sample-attrition.tsv",
    "sample-concordance.tsv",
    "sample-metadata.tsv",
    "sample-protein-richness.tsv",
    "software-versions.json",
    "source-manifest.json",
    "subject-concordance.tsv",
    "subject-dna-ec-relative.tsv.gz",
    "subject-protein-ec-relative.tsv.gz",
    "subject-rna-ec-relative.tsv.gz",
    "threshold-audit.tsv",
)

SCRIPT_FILES = (
    "download_article65_metaproteomics_data.py",
    "prepare_article65_metaproteomics.py",
    "plot_article65_metaproteomics.py",
    "freeze_article65_metaproteomics.py",
    "validate_article65_metaproteomics.py",
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
    if metrics.get("article") != 65:
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
        "article": 65,
        "dataset": "IBDMDB HMP2 MGX/MTX/MPX plus exact-ID MBX product availability",
        "doi": "10.1038/s41586-019-1237-9",
        "source_tables_included": False,
        "source_manifest_and_checksums_included": True,
        "streamed_unstratified_ec_tables_included": True,
        "transformed_matrices_included": True,
        "analysis_samples": metrics["analysis_samples"],
        "independent_subjects": metrics["independent_subjects"],
        "selected_ecs": metrics["selected_ecs"],
        "complete_four_layer_samples": metrics["complete_four_layer_samples"],
        "seed": metrics["seed"],
        "plot_seed": metrics["plot_seed"],
        "interpretation": "relative functional concordance; not absolute protein expression, translation rate, enzyme activity, flux, or causality",
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    notice = """Article 65 frozen evidence bundle

The source data are the official IBDMDB HMP2 merged metaproteomic protein-
filter and EC products, merged metagenomic and metatranscriptomic HUMAnN EC
relative-abundance tables, the 2018-08-20 metadata release, and the published
Supplementary Figure 1 PDF. Exact URLs, byte counts, and SHA-256 checksums are
retained in source-manifest.json. The approximately 122 MB of original source
files are fetched by the included downloader rather than duplicated here.

The bundle includes streamed unstratified MGX/MTX EC tables, the exact 186-
sample/76-subject transformed matrices, all four protein-threshold audits,
namespace and identifier-overlap ledgers, all 789 EC correlation tests, scripts,
and environment records. One triple-assay profile with zero MPX EC evidence was
excluded. Metabolomics is represented only by exact biospecimen-level product
availability in this chapter; no metabolite values are imputed. Associations
among independently normalized profiles do not estimate absolute protein
expression, translation rate, enzyme activity, metabolic flux, or causality.
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
    print(f"Frozen Article 65 bundle: {len(payloads)} payload files in {target}")


if __name__ == "__main__":
    main()
