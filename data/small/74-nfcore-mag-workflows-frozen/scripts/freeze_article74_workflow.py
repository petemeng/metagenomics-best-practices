#!/usr/bin/env python3
"""Create a checksum-covered Article 74 workflow engineering bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "analysis-metrics.json",
    "database-lock.tsv",
    "execution-unit-matrix.tsv",
    "figure-manifest.json",
    "funcscan-branch-contract.tsv",
    "funcscan-metromap-original.png",
    "hardware-envelope.tsv",
    "hpc.slurm.config",
    "input-manifest.tsv",
    "mag-metromap-original.png",
    "methods-contract.json",
    "parameter-precedence.tsv",
    "params.publication.yml",
    "params.stub.yml",
    "pipeline-defaults-audit.tsv",
    "production-command.sh",
    "profile-contract.tsv",
    "profile-parse-evidence.tsv",
    "provenance-bundle.tsv",
    "release-lock.tsv",
    "runtime-environment.tsv",
    "runtime-summary.tsv",
    "samplesheet.csv",
    "source-manifest.json",
    "stub-runtime-trace.tsv",
    "stub-scope.tsv",
)

SCRIPT_FILES = (
    "download_article74_workflow.py",
    "prepare_article74_workflow.py",
    "plot_article74_workflow.py",
    "freeze_article74_workflow.py",
    "validate_article74_workflow.py",
)

ENV_FILES = ("assembly.yml", "multiomics-python.yml")


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
    (staging / "source").mkdir()
    (staging / "scripts").mkdir()
    (staging / "env").mkdir()

    missing = [name for name in WORK_FILES if not (work / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing Article 74 work files: {missing}")
    source_files = sorted(path for path in (work / "source").iterdir() if path.is_file())
    if len(source_files) != 19:
        raise ValueError(f"Expected 19 selected source snapshots, observed {len(source_files)}")

    for name in WORK_FILES:
        shutil.copy2(work / name, staging / name)
    for path in source_files:
        shutil.copy2(path, staging / "source" / path.name)
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
        "article": 74,
        "payload_files": len(WORK_FILES),
        "selected_source_files": len(source_files),
        "script_files": len(SCRIPT_FILES),
        "environment_files": len(ENV_FILES),
        "contract": (
            "nf-core/mag 5.5.0, nf-core/funcscan 4.0.0, Nextflow 26.04.0, "
            "the scientific parameter file, SLURM/Apptainer infrastructure, database "
            "checksums, input checksums, stub-run boundary, and resume trace are frozen."
        ),
    }
    (staging / "bundle-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    paths = sorted(path for path in staging.rglob("*") if path.is_file())
    lines = [f"{sha256(path)}  {path.relative_to(staging).as_posix()}" for path in paths]
    (staging / "file-checksums.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    if output.exists():
        shutil.rmtree(output)
    staging.replace(output)
    print(f"frozen\t{output}\t{len(lines)} checksum-covered files")


if __name__ == "__main__":
    main()
