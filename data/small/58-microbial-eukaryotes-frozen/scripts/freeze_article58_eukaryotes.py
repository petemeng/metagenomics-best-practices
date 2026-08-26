#!/usr/bin/env python3
"""Create the compact, checksum-covered Article 58 frozen bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


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


def copy(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    target = args.output_dir.resolve()
    if not (work / ".article58-summary-complete").is_file():
        raise FileNotFoundError("Article 58 summary sentinel is missing")
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)

    for path in sorted((work / "summary").glob("*")):
        if path.is_file():
            copy(path, target / path.name)

    for name in (
        "asset-check-audit.tsv",
        "prepared-fastq-audit.tsv",
        "reference-sequence-ledger.tsv",
        "eukrep-fragment-ledger.tsv",
        "run-contract.json",
        "command-log.tsv",
        "tool-versions.tsv",
        "database-audit.tsv",
    ):
        copy(work / name, target / name)

    eukdetect_root = work / "results/eukdetect"
    for relative in (
        "Zymo_D6300_filtered_hits_table.txt",
        "Zymo_D6300_filtered_hits_eukfrac.txt",
        "Zymo_D6300.normalized.tsv",
        "filtering/Zymo_D6300_all_hits_table.txt",
    ):
        copy(eukdetect_root / relative, target / "eukdetect" / relative)

    for branch in ("eukrep-reference", "eukrep-assembly"):
        source = work / "results" / branch
        for path in sorted(source.glob("*.names")):
            copy(path, target / branch / path.name)
    for name in ("fastp-eukdetect.json", "fastp-assembly.json"):
        copy(work / "results/qc" / name, target / "qc" / name)

    for path in sorted((work / "logs").glob("*")):
        if path.is_file():
            copy(path, target / "logs" / path.name)

    script_names = (
        "download_article58_eukdetect2_db.sh",
        "download_article58_eukaryote_sources.sh",
        "prepare_article58_eukaryotes.py",
        "run_article58_eukaryotes.py",
        "summarize_article58_eukaryotes.py",
        "plot_article58_eukaryotes.R",
        "freeze_article58_eukaryotes.py",
        "validate_article58_eukaryotes.py",
        "article42_44_validation_utils.py",
    )
    for name in script_names:
        copy(root / "scripts" / name, target / "scripts" / name)

    for name in (
        "microbial-eukaryotes.yml",
        "microbial-eukaryotes-linux-64.lock",
        "microbial-eukaryotes-pip-lock.txt",
        "eukrep-legacy.yml",
        "eukrep-legacy-linux-64.lock",
        "assembly.yml",
        "assembly-linux-64.lock",
        "read-qc.yml",
        "read-qc-linux-64.lock",
    ):
        copy(root / "env" / name, target / "env" / name)

    copy(
        root / "data/small/58-eukaryote-source-manifest.tsv",
        target / "source-manifest.tsv",
    )
    copy(
        root / "data/small/58-eukdetect2-database-manifest.tsv",
        target / "eukdetect2-database-manifest.tsv",
    )

    notice = """Article 58 frozen evidence bundle

The complete SRR12324253 paired FASTQ archives, the official Zymo v2 reference
archive, the 20-million-pair prepared FASTQ branch, MEGAHIT assembly, and the
7.12-GB EukDetect2 database are not duplicated here. Their public identities,
byte counts, immutable checksums, commands, environment locks, resource logs,
and compact derived evidence tables are retained.

EukDetect2 is represented by its 2026 preprint and the database release dated
16 March 2026. The installed conda record reports 2.0.2, Python distribution
metadata reports 2.0.0, and the CLI banner reports v2.0.1; all three values are
preserved rather than silently unified.

The EukRep reference benchmark is deterministic: 80 evenly spaced,
non-overlapping fragments per species at each of 3, 5, 10, and 20 kb from ten
official Zymo genomes. It is an exact mock-specific audit, not a universal
performance estimate.
"""
    (target / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    contract = {
        "article": 58,
        "created_from": str(work.relative_to(root)),
        "large_inputs_included": False,
        "random_output_requested": False,
        "seed": 20260758,
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    files = sorted(
        path
        for path in target.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    lines = [f"{sha256(path)}  {path.relative_to(target).as_posix()}" for path in files]
    (target / "file-checksums.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(f"Frozen Article 58 bundle: {len(files)} payload files in {target}")


if __name__ == "__main__":
    main()
