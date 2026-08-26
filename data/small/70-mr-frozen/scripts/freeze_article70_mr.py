#!/usr/bin/env python3
"""Create a checksum-covered Article 70 Mendelian-randomization bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "analysis-metrics.json",
    "bmi-instruments-raw.tsv.gz",
    "chd-associations-raw.tsv.gz",
    "design-audit.tsv",
    "egger-intercept.tsv",
    "harmonisation-audit.tsv",
    "harmonised-instruments.rds",
    "harmonised-instruments.tsv.gz",
    "hemani-figure1-original.png",
    "input-attrition.tsv",
    "instrument-strength-summary.tsv",
    "leave-one-out.tsv",
    "methods-contract.json",
    "model-metrics.json",
    "mr-analysis-object.rds",
    "mr-estimates.tsv",
    "mr-heterogeneity.tsv",
    "mr-presso-estimates.tsv",
    "mr-presso-object.rds",
    "mr-presso-outliers.tsv",
    "mr-presso-tests.tsv",
    "presso-outlier-exclusion-estimates.tsv",
    "r-session-info.txt",
    "radial-ivw-estimates.tsv",
    "radial-ivw-outliers.tsv",
    "single-snp-estimates.tsv.gz",
    "software-versions-analysis.tsv",
    "software-versions-preparation.tsv",
    "source-manifest.json",
    "steiger-directionality.tsv",
    "steiger-per-snp.tsv",
    "twosamplemr-vig-perform-mr.RData",
)

SCRIPT_FILES = (
    "download_article70_mr_data.py",
    "prepare_article70_mr_inputs.R",
    "run_article70_mr_analysis.R",
    "plot_article70_mr.py",
    "freeze_article70_mr.py",
    "validate_article70_mr.py",
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
        raise FileNotFoundError(f"Missing Article 70 work files: {missing}")
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
        "article": 70,
        "payload_files": len(WORK_FILES),
        "script_files": len(SCRIPT_FILES),
        "environment_files": len(ENV_FILES),
        "source_work_dir": str(work),
        "contract": (
            "The exact TwoSampleMR vignette associations and all derived MR "
            "diagnostics are frozen; OpenGWAS API access is not required."
        ),
    }
    (staging / "bundle-manifest.json").write_text(
        json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    checksum_paths = sorted(path for path in staging.rglob("*") if path.is_file())
    lines = [
        f"{sha256(path)}  {path.relative_to(staging).as_posix()}"
        for path in checksum_paths
    ]
    (staging / "file-checksums.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    if output.exists():
        shutil.rmtree(output)
    staging.replace(output)
    print(f"frozen\t{output}\t{len(lines)} checksum-covered files")


if __name__ == "__main__":
    main()
