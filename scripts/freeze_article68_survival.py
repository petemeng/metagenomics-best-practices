#!/usr/bin/env python3
"""Create a checksum-covered Article 68 survival-analysis evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "adjusted-cox-model.rds",
    "analysis-metrics.json",
    "calibration-365d.tsv",
    "cohort-attrition.tsv",
    "composition-audit.tsv",
    "cox-model-audit.tsv",
    "cox-model-estimates.tsv",
    "cutoff-evaluation.tsv",
    "cutoff-km-curves.tsv",
    "cutoff-risk-table.tsv",
    "cutoff-search-full-leaky.tsv",
    "cutoff-search-training.tsv",
    "cutoff-split-ledger.tsv",
    "cv-performance-summary.tsv",
    "faecalibacterium-feature-audit.tsv",
    "faecalibacterium-schoenfeld.tsv",
    "incremental-model-test.tsv",
    "leave-one-out-influence.tsv",
    "lkt-feature-map.tsv",
    "lkt-survival-ppm.tsv.gz",
    "metadata-variable-audit.tsv",
    "methods-contract.json",
    "model-metrics.json",
    "nonlinearity-test.tsv",
    "performance-bootstrap.tsv",
    "proportional-hazards-tests.tsv",
    "r-session-info.txt",
    "repeated-cv-fold-audit.tsv",
    "repeated-cv-predictions-long.tsv",
    "repeated-cv-predictions.tsv",
    "software-versions-python.json",
    "software-versions-r.tsv",
    "source-manifest.json",
    "spencer-data-dictionary.xlsx",
    "spencer-fig3a-original.png",
    "spline-cox-model.rds",
    "spline-effect-curve.tsv",
    "survival-cohort.tsv",
    "time-dependent-roc-365d.tsv",
)

SCRIPT_FILES = (
    "download_article68_survival_data.py",
    "prepare_article68_survival.py",
    "run_article68_survival_models.R",
    "plot_article68_survival.py",
    "freeze_article68_survival.py",
    "validate_article68_survival.py",
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
        raise FileNotFoundError(f"Missing Article 68 work files: {missing}")
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

    bundle = {
        "article": 68,
        "payload_files": len(WORK_FILES),
        "script_files": len(SCRIPT_FILES),
        "environment_files": len(ENV_FILES),
        "source_work_dir": str(work),
        "contract": "All derived statistics are frozen; raw 22.7 MB public workbook is checksum-pinned and downloaded separately.",
    }
    (staging / "bundle-manifest.json").write_text(
        json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    checksum_paths = sorted(path for path in staging.rglob("*") if path.is_file())
    lines = [f"{sha256(path)}  {path.relative_to(staging).as_posix()}" for path in checksum_paths]
    (staging / "file-checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")

    if output.exists():
        shutil.rmtree(output)
    staging.replace(output)
    print(f"frozen\t{output}\t{len(lines)} checksum-covered files")


if __name__ == "__main__":
    main()
