#!/usr/bin/env python3
"""Create a checksum-covered Article 73 public-resource bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "accession-crosswalk.tsv",
    "analysis-metrics.json",
    "download-decisions.tsv",
    "gem-biome-summary.tsv",
    "gem-figure1-original.png",
    "gem-map-visualization-sample.tsv",
    "gem-metadata-preview.tsv",
    "gem-quality-visualization-sample.tsv",
    "gtdb-release-history.tsv",
    "gtdb-selected-files.tsv",
    "methods-contract.json",
    "mgnify-catalogues.tsv",
    "resource-registry.tsv",
    "source-manifest.json",
    "tara-analyses.tsv",
    "tara-analysis-lineage.tsv",
    "tara-metadata-completeness.tsv",
    "tara-samples.tsv",
    "tara-study-downloads.tsv",
)

SOURCE_FILES = (
    "gem-figure1-original.png",
    "gem-genome-metadata.tsv",
    "gem-readme.md",
    "gtdb-r226-release-notes.txt",
    "gtdb-r226-version.txt",
    "gtdb-r232-md5.txt",
    "gtdb-r232-release-notes.txt",
    "gtdb-r232-version.txt",
    "mgnify-analyses-page1.json",
    "mgnify-analyses-page2.json",
    "mgnify-analyses-page3.json",
    "mgnify-catalogues.json",
    "mgnify-marine-v2.json",
    "mgnify-samples-page1.json",
    "mgnify-samples-page2.json",
    "mgnify-study.json",
)

SCRIPT_FILES = (
    "download_article73_resources.py",
    "prepare_article73_resources.py",
    "plot_article73_resources.py",
    "freeze_article73_resources.py",
    "validate_article73_resources.py",
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
    (staging / "source").mkdir()
    (staging / "scripts").mkdir()
    (staging / "env").mkdir()

    missing = [name for name in WORK_FILES if not (work / name).is_file()]
    missing += [f"source/{name}" for name in SOURCE_FILES if not (work / "source" / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing Article 73 work files: {missing}")

    for name in WORK_FILES:
        shutil.copy2(work / name, staging / name)
    for name in SOURCE_FILES:
        shutil.copy2(work / "source" / name, staging / "source" / name)
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
        "article": 73,
        "payload_files": len(WORK_FILES) + len(SOURCE_FILES),
        "script_files": len(SCRIPT_FILES),
        "environment_files": len(ENV_FILES),
        "source_work_dir": str(work),
        "contract": (
            "The MGnify v2 study/sample/run/analysis lineage, GEM metadata, "
            "catalogue boundaries, GTDB R232 release identity, official source "
            "responses, and paper anchor are frozen for offline audit."
        ),
    }
    (staging / "bundle-manifest.json").write_text(
        json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    paths = sorted(path for path in staging.rglob("*") if path.is_file())
    lines = [
        f"{sha256(path)}  {path.relative_to(staging).as_posix()}" for path in paths
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
