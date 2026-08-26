#!/usr/bin/env python3
"""Create a checksum-covered Article 77 local release-readiness bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "analysis-metrics.json",
    "availability-statements.md",
    "code-availability-components.tsv",
    "container-ledger.tsv",
    "database-manifest-public.tsv",
    "data-availability-components.tsv",
    "data-NOTICE.txt",
    "figure-manifest.json",
    "identifier-registry.tsv",
    "methods-contract.json",
    "object-relationship.tsv",
    "package-index.tsv",
    "policy-assertions.tsv",
    "policy-source-registry.tsv",
    "release-artifact-manifest.tsv",
    "release-gate-ledger.tsv",
    "release-readiness.tsv",
    "repository-routing-matrix.tsv",
    "source-artifact-manifest.tsv",
    "source-manifest.json",
    "tenhoopen-figure2-original.jpg",
)

SCRIPT_FILES = (
    "download_article77_release_sources.py",
    "prepare_article77_release.py",
    "plot_article77_release.py",
    "freeze_article77_release.py",
    "validate_article77_release.py",
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
        raise FileNotFoundError(f"Missing Article 77 work files: {missing}")
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
        "article": 77,
        "payload_files": len(WORK_FILES),
        "script_files": len(SCRIPT_FILES),
        "environment_files": len(ENV_FILES),
        "contract": (
            "Repository routing, object/accession lineage, local artifact checksums, "
            "database and container ledgers, release gates, and Data/Code Availability "
            "drafts are frozen without performing an external submission."
        ),
        "publication_boundary": (
            "This is a local readiness packet, not an SRA/ENA, Zenodo, GitHub or OCI "
            "release. No pending token is a public identifier."
        ),
        "ownership_boundary": (
            "PRJEB52977 tutorial inputs are third-party public records and are cited, "
            "not resubmitted as investigator-owned data."
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
