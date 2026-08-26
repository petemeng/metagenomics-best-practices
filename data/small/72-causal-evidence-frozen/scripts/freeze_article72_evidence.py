#!/usr/bin/env python3
"""Create a checksum-covered Article 72 causal-evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


WORK_FILES = (
    "analysis-metrics.json",
    "buffie-figure4-original.jpg",
    "claim-contracts.tsv",
    "claim-downgrade-examples.tsv",
    "evidence-ledger.tsv",
    "evidence-packet.tsv",
    "human-intervention-outcomes.tsv",
    "methods-contract.json",
    "publication-metadata.tsv",
    "rung-definitions.tsv",
    "source-manifest.json",
    "study-domain-coverage.tsv",
)

SCRIPT_FILES = (
    "download_article72_evidence.py",
    "prepare_article72_evidence.py",
    "plot_article72_evidence.py",
    "freeze_article72_evidence.py",
    "validate_article72_evidence.py",
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
        raise FileNotFoundError(f"Missing Article 72 work files: {missing}")
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
        "article": 72,
        "payload_files": len(WORK_FILES),
        "script_files": len(SCRIPT_FILES),
        "environment_files": len(ENV_FILES),
        "source_work_dir": str(work),
        "contract": (
            "Eleven DOI identities, ten primary-study evidence rows, three "
            "claim contracts, seven non-additive domains, randomized-trial "
            "numerators, and the official anchor are frozen for offline audit."
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
