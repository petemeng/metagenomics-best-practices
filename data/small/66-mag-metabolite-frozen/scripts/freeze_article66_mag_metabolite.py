#!/usr/bin/env python3
"""Create a checksum-covered Article 66 MAG-metabolite evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "analysis-metrics.json",
    "combined-donor-pathways.tsv",
    "community-gapfill-audit.tsv",
    "community-model-audit.tsv",
    "individual-model-audit.tsv",
    "individual-model-summary.tsv",
    "mag-ledger.tsv",
    "mag-quality-summary.tsv",
    "majzoub-fig2-original.jpg",
    "metabolite-overlap-audit.tsv",
    "metabolomics-coverage.tsv",
    "methods-contract.json",
    "model-comparison-summary.tsv",
    "pathway-counts.tsv",
    "pathways-between-donors.tsv",
    "pathways-combined.tsv",
    "pathways-within-donor.tsv",
    "phenotype-evidence.tsv",
    "quality-function-correlations.tsv",
    "software-versions.json",
    "source-manifest.json",
    "source-anomaly-ledger.tsv",
)

SCRIPT_FILES = (
    "download_article66_mag_metabolite_data.py",
    "prepare_article66_mag_metabolite.py",
    "plot_article66_mag_metabolite.py",
    "freeze_article66_mag_metabolite.py",
    "validate_article66_mag_metabolite.py",
    "article42_44_validation_utils.py",
)

ENV_FILES = ("multiomics-python.yml", "multiomics-r-packages.tsv")


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
    if metrics.get("article") != 66:
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
        "article": 66,
        "study": "Majzoub et al. mSystems 2024",
        "doi": "10.1128/msystems.00746-24",
        "pmcid": "PMC11406951",
        "human_metagenomics_accession": "PRJEB50699",
        "source_workbook_included": False,
        "source_figure_included": True,
        "source_manifest_and_checksums_included": True,
        "derived_mag_rows": metrics["mag_total"],
        "normalized_mag_rows": metrics["normalized_mag_total"],
        "independent_donor_units": 2,
        "longitudinal_metabolomics_samples": metrics["metabolomics_longitudinal_samples"],
        "seed": metrics["seed"],
        "plot_seed": metrics["plot_seed"],
        "interpretation": "genome-content-derived metabolic predictions and donor-level presence/absence confirmation; not measured microbial flux, metabolite source attribution, or causal phenotype mediation",
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    notice = """Article 66 frozen evidence bundle

The source study is Majzoub et al., mSystems 2024 (PMCID PMC11406951,
DOI 10.1128/msystems.00746-24). The official full-text XML and supplementary
archive are fetched from Europe PMC with exact URL, byte-count and SHA-256
contracts in source-manifest.json. The original workbook is not duplicated in
this bundle; the included downloader retrieves it. The open article Figure 2
JPEG is retained with its checksum for the visual anchor.

Derived tables retain all 170 human MAG quality rows, 142 normalized individual
GEM rows, all 14 Table 1 community-model configurations, Figure 2 overlap
arithmetic, the longitudinal metabolomics coverage contract, all published
pathway lists from Tables S5-S7, and the two-donor phenotype evidence boundary.
The study's published high-quality MAG label uses CheckM completeness and
contamination only. It cannot be upgraded to full MIMAG high-quality status
without rRNA and tRNA evidence. Metabolite confirmation is donor-level
presence/absence and does not estimate sensitivity, specificity, accuracy,
microbial source attribution, metabolic flux, phenotype mediation, or causality.
"""
    (target / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    payloads = sorted(
        path
        for path in target.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    (target / "file-checksums.sha256").write_text(
        "\n".join(
            f"{sha256(path)}  {path.relative_to(target).as_posix()}" for path in payloads
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Frozen Article 66 bundle: {len(payloads)} payload files in {target}")


if __name__ == "__main__":
    main()
