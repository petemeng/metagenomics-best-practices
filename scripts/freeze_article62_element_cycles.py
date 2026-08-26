#!/usr/bin/env python3
"""Create a checksum-covered, self-contained Article 62 evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


ANALYSIS_FILES = (
    "analysis-contract.json",
    "analysis-metrics.json",
    "carrier-phylum-summary.tsv",
    "community-mag-concordance.tsv",
    "environment-associations.tsv",
    "mag-carrier-summary.tsv",
    "mag-process-evidence.tsv.gz",
    "mag-recovery-ceiling.tsv",
    "missing-ko-columns.tsv",
    "nitrogen-step-summary.tsv",
    "process-rules.tsv",
    "run-ledger.tsv",
    "sample-carrier-fraction.tsv.gz",
    "sample-process-index.tsv.gz",
    "selected-sample-metadata.tsv",
    "software-versions.json",
    "spring-carrier-fraction.tsv",
    "spring-process-index.tsv",
    "temperature-regime-summary.tsv",
)

INPUT_FILES = (
    "download-manifest.json",
    "sample-metadata.tsv",
    "ko-proportions-in-metagenomes.tsv.gz",
    "mag-metadata.tsv",
    "mag-abundances-per-sample.biom",
    "kos-in-mags.tsv.gz",
    "diting-pathway-formulas-v0.3.txt",
)

SCRIPT_FILES = (
    "download_article62_element_data.py",
    "run_article62_element_cycles.py",
    "plot_article62_element_cycles.py",
    "freeze_article62_element_cycles.py",
    "validate_article62_element_cycles.py",
    "article42_44_validation_utils.py",
)

ENV_FILES = (
    "drep.yml",
    "drep-linux-64.lock",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--analysis-dir", type=Path, required=True)
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


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    analysis = args.analysis_dir.resolve()
    cache = args.cache_dir.resolve()
    target = args.output_dir.resolve()

    for name in ANALYSIS_FILES:
        if not (analysis / name).is_file():
            raise FileNotFoundError(f"Article 62 analysis output missing: {name}")
    for name in INPUT_FILES:
        if not (cache / name).is_file():
            raise FileNotFoundError(f"Article 62 source input missing: {name}")

    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)

    for name in ANALYSIS_FILES:
        copy(analysis / name, target / name)
    for name in INPUT_FILES:
        copy(cache / name, target / "inputs" / name)
    copy(
        root / "data/small/62-element-cycle-marker-rules.tsv",
        target / "inputs/62-element-cycle-marker-rules.tsv",
    )
    for name in SCRIPT_FILES:
        copy(root / "scripts" / name, target / "scripts" / name)
    for name in ENV_FILES:
        copy(root / "env" / name, target / "env" / name)

    metrics = json.loads((analysis / "analysis-metrics.json").read_text(encoding="utf-8"))
    contract = {
        "article": 62,
        "created_from": str(analysis.relative_to(root)),
        "dataset_doi": "10.6084/m9.figshare.30284068.v2",
        "paper_doi": "10.1038/s41597-026-07139-w",
        "samples": metrics["samples"],
        "independent_hot_springs": metrics["hot_springs"],
        "mags": metrics["mags"],
        "processes": metrics["processes"],
        "source_tables_included": True,
        "seed": metrics["seed"],
        "permutations": 9999,
        "bootstraps": 2000,
        "primary_inference_unit": "hot-spring median",
        "interpretation": "genetic potential; no activity, direction, or rate claim",
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    notice = """Article 62 frozen evidence bundle

The biological inputs are the exact Figshare v2 tables accompanying the 2026
Scientific Data survey of 500 shotgun metagenomes from 56 western-US hot
springs. They comprise published community KO proportions, metadata for 780
quality-filtered MAGs, MAG KO counts, the 780 x 500 MAG abundance BIOM table,
and sample metadata. Exact byte counts and SHA-256 checksums are retained.

The 20 process rules are explicit, checksum-covered marker screens informed by
DiTing v0.3 and strengthened with carrier-completeness and homology caveats.
Community associations use one median per hot spring. MAG abundance fractions
are conditional on the recovered-MAG pool because the median whole-metagenome
read recruitment is low. Marker presence indicates genetic potential; it does
not establish expression, reaction direction, process rate, or metabolite
handoff.
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
    print(f"Frozen Article 62 bundle: {len(payloads)} payload files in {target}")


if __name__ == "__main__":
    main()
