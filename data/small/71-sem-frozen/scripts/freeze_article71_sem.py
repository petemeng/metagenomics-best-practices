#!/usr/bin/env python3
"""Create a checksum-covered Article 71 piecewise-SEM evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "all-sample-metrics.tsv",
    "analysis-metrics.json",
    "antibiotic-overlap-by-diagnosis.tsv",
    "directed-separation-claims.tsv",
    "exposure-by-cohort-diagnosis.tsv",
    "faecalibacterium-path-bootstrap.tsv.gz",
    "franzosa-fig1-original.png",
    "leave-one-out-paths.tsv",
    "local-model-diagnostics.tsv",
    "local-path-coefficients-hc3.tsv",
    "methods-contract.json",
    "model-metrics.json",
    "node-contract.tsv",
    "outcome-path-transport.tsv",
    "path-effect-summary.tsv",
    "propensity-positivity-audit.tsv",
    "r-session-info.txt",
    "sample-attrition.tsv",
    "sem-fit-comparison.tsv",
    "sem-model-objects.rds",
    "sem-path-bootstrap.tsv.gz",
    "sem-primary-cohort.tsv",
    "sem-validation-cohort.tsv",
    "software-versions-python.json",
    "software-versions-r.tsv",
    "source-manifest.json",
    "standardization-parameters.tsv",
    "variable-summary.tsv",
    "variance-inflation.tsv",
)

SCRIPT_FILES = (
    "download_article71_sem_data.py",
    "prepare_article71_sem.py",
    "run_article71_sem_models.R",
    "plot_article71_sem.py",
    "freeze_article71_sem.py",
    "validate_article71_sem.py",
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
        raise FileNotFoundError(f"Missing Article 71 work files: {missing}")
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
        "article": 71,
        "payload_files": len(WORK_FILES),
        "script_files": len(SCRIPT_FILES),
        "environment_files": len(ENV_FILES),
        "source_work_dir": str(work),
        "contract": (
            "The pinned Franzosa genus profiles, complete-case cohorts, local "
            "path models, 5,000 bootstrap refits, positivity diagnostics, and "
            "transport checks are frozen for offline rendering and audit."
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
