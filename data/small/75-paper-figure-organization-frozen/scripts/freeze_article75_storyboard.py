#!/usr/bin/env python3
"""Create a checksum-covered Article 75 publication-storyboard bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "analysis-metrics.json",
    "claim-evidence-matrix.tsv",
    "data-NOTICE.txt",
    "evidence-language-ladder.tsv",
    "figure-manifest.json",
    "figure-style-contract.tsv",
    "main-figure-storyboard.tsv",
    "main-supplement-map.tsv",
    "methods-contract.json",
    "panel-register.tsv",
    "paper-source-manifest.json",
    "result-traceability-ledger.tsv",
    "reviewer-attack-map.tsv",
    "sensitivity-matrix.tsv",
    "series-evidence-metrics.tsv",
    "source-artifact-manifest.tsv",
    "version-ledger-example.tsv",
    "wirbel-figure1-original.jpg",
    "wirbel-main-figure-ledger.tsv",
)

SCRIPT_FILES = (
    "download_article75_paper.py",
    "prepare_article75_storyboard.py",
    "plot_article75_storyboard.py",
    "freeze_article75_storyboard.py",
    "validate_article75_storyboard.py",
)

ENV_FILES = ("multiomics-python.yml", "renv.lock")


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
        raise FileNotFoundError(f"Missing Article 75 work files: {missing}")
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
        "article": 75,
        "payload_files": len(WORK_FILES),
        "script_files": len(SCRIPT_FILES),
        "environment_files": len(ENV_FILES),
        "contract": (
            "The Wirbel 2019 figure arc, five-figure claim storyboard, panel units, "
            "main/supplement placement, claim-evidence matrix, sensitivity plan, "
            "version ledger, result traceability, and public-use boundary are frozen."
        ),
        "source_boundary": (
            "Numerical examples come from different public tutorial datasets and "
            "cannot be combined into one biological claim."
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
