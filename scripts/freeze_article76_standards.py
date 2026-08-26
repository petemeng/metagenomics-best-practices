#!/usr/bin/env python3
"""Create a checksum-covered Article 76 reporting-standards bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "analysis-metrics.json",
    "article44-mimag-compliance.tsv",
    "article54-miuvig-compliance.tsv",
    "checklist-items.tsv",
    "checklist-section-counts.tsv",
    "data-NOTICE.txt",
    "field-responsibility-matrix.tsv",
    "figure-manifest.json",
    "methods-contract.json",
    "mimag-quality-criteria.tsv",
    "miuvig-mandatory-metadata.tsv",
    "miuvig-quality-categories.tsv",
    "not-applicable-ledger.tsv",
    "reporting-layer-map.tsv",
    "source-artifact-manifest.tsv",
    "source-manifest.json",
    "standard-selection-matrix.tsv",
    "standards-crosswalk.tsv",
    "streams-figure1-original.png",
    "submission-readiness.tsv",
)

SCRIPT_FILES = (
    "download_article76_standards.py",
    "prepare_article76_standards.py",
    "plot_article76_standards.py",
    "freeze_article76_standards.py",
    "validate_article76_standards.py",
)

ENV_FILES = ("multiomics-python.yml", "renv.lock")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    output = args.output_dir.resolve()
    staging = output.with_name(output.name + ".staging")
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    (staging / "scripts").mkdir()
    (staging / "env").mkdir()

    missing = [name for name in WORK_FILES if not (work / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing Article 76 work files: {missing}")
    for name in WORK_FILES:
        shutil.copy2(work / name, staging / name)
    for name in SCRIPT_FILES:
        source = root / "scripts" / name
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, staging / "scripts" / name)
    for name in ENV_FILES:
        source = root / "env" / name
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, staging / "env" / name)

    manifest = {
        "article": 76,
        "payload_files": len(WORK_FILES),
        "script_files": len(SCRIPT_FILES),
        "environment_files": len(ENV_FILES),
        "contract": (
            "STORMS v1.03, STREAMS v1.0, MIMAG and MIUViG source identities, "
            "checklist rows, crosswalks, real MAG/UViG audits, ownership ledger, "
            "N/A decisions and readiness records are frozen."
        ),
        "source_boundary": (
            "Article 44 MAGs and Article 54 virus fixtures are independent public "
            "examples and cannot be combined into one biological claim."
        ),
    }
    (staging / "bundle-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    paths = sorted(path for path in staging.rglob("*") if path.is_file())
    lines = [f"{sha256(path)}  {path.relative_to(staging).as_posix()}" for path in paths]
    (staging / "file-checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if output.exists():
        shutil.rmtree(output)
    staging.replace(output)
    print(f"frozen\t{output}\t{len(lines)} checksum-covered files")


if __name__ == "__main__":
    main()
