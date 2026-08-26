#!/usr/bin/env python3
"""Freeze compact, checksum-covered Article 36 annotation evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
from pathlib import Path


SEED = 20260736
SCRIPT_NAMES = (
    "bootstrap_article36_eggnog_database.sh",
    "prepare_article36_eggnog_inputs.py",
    "run_article36_eggnog_annotation.py",
    "summarize_article36_eggnog_annotation.py",
    "freeze_article36_eggnog_annotation.py",
    "validate_article36_eggnog_annotation.py",
)
ROOT_FILES = (
    "input-audit.tsv",
    "input-lineage.tsv",
    "database-audit.tsv",
    "tool-versions.tsv",
    "command-log.tsv",
    "run-contract.json",
    "run-metadata.json",
)
SUMMARY_FILES = (
    "gene-functional-annotation.tsv.gz",
    "annotation-state-summary.tsv",
    "field-coverage-summary.tsv",
    "annotation-evidence-overlap.tsv",
    "completeness-annotation-summary.tsv",
    "length-annotation-summary.tsv",
    "cog-category-summary.tsv",
    "ko-term-summary.tsv",
    "go-term-summary.tsv",
    "fractional-allocation-audit.tsv",
    "go-evidence-audit.tsv",
    "delgado-table3-eggnog-anchor.tsv",
    "resource-usage.tsv",
    "run-summary.json",
)
LOG_FILES = (
    "emapper-main.stdout.log",
    "emapper-main.stderr.log",
    "emapper-main.time.txt",
    "emapper-go-all.stdout.log",
    "emapper-go-all.stderr.log",
    "emapper-go-all.time.txt",
)


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


def copy_required(source: Path, destination: Path) -> None:
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def copy_log(source: Path, destination: Path) -> None:
    """Copy a run log even when the audited command emitted no text."""
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def deterministic_gzip(source: Path, destination: Path) -> None:
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as src, destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as dst:
            shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    output = args.output_dir.resolve()
    if not (work / ".article36-summary-complete").is_file():
        raise FileNotFoundError("Article 36 summary is incomplete")
    stage = output.parent / f".{output.name}.staging"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)

    for name in ROOT_FILES:
        copy_required(work / name, stage / name)
    for name in SUMMARY_FILES:
        copy_required(work / "summary" / name, stage / name)
    for name in LOG_FILES:
        copy_log(work / "logs" / name, stage / "logs" / name)
    for name in SCRIPT_NAMES:
        copy_required(root / "scripts" / name, stage / "scripts" / name)
    copy_required(root / "env/eggnog-annotation.yml", stage / "env/eggnog-annotation.yml")
    copy_required(root / "env/eggnog-annotation-linux-64.lock", stage / "env/eggnog-annotation-linux-64.lock")
    copy_required(
        root / "data/small/36-eggnog-database-manifest.tsv",
        stage / "database-manifest.tsv",
    )
    copy_required(
        root / "data/small/34-nonredundant-gene-catalog-frozen/sources/PMC9074274.fullTextXML",
        stage / "sources/PMC9074274.fullTextXML",
    )
    copy_required(
        root / "data/small/34-nonredundant-gene-catalog-frozen/data-NOTICE.txt",
        stage / "sources/article34-data-NOTICE.txt",
    )
    copy_required(
        root / "data/small/35-gene-abundance-frozen/data-NOTICE.txt",
        stage / "sources/article35-data-NOTICE.txt",
    )

    deterministic_gzip(
        work / "annotation/main/catalog-main.emapper.annotations",
        stage / "raw/catalog-main.emapper.annotations.gz",
    )
    deterministic_gzip(
        work / "annotation/main/catalog-main.emapper.seed_orthologs",
        stage / "raw/catalog-main.emapper.seed_orthologs.gz",
    )
    deterministic_gzip(
        work / "annotation/go-all/catalog-go-all.emapper.annotations",
        stage / "raw/catalog-go-all.emapper.annotations.gz",
    )

    run_summary = json.loads((work / "summary/run-summary.json").read_text(encoding="utf-8"))
    run_contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))
    database_rows = list(
        csv.DictReader(
            (root / "data/small/36-eggnog-database-manifest.tsv").open(encoding="utf-8"),
            delimiter="\t",
        )
    )
    source_rows = [
        {
            "Source": "Delgado and Andersson paper XML",
            "Identifier": "PMC9074274 / DOI 10.1186/s40168-022-01259-2",
            "Release": "published 2022-05-06",
            "SHA256": sha256(root / "data/small/34-nonredundant-gene-catalog-frozen/sources/PMC9074274.fullTextXML"),
            "Role": "historical method and Table 3 anchor",
        },
        {
            "Source": "Article 34 protein catalog",
            "Identifier": "MEGAHIT mix primary representatives",
            "Release": "Article 34 frozen bundle",
            "SHA256": run_contract["input_faa_sha256"],
            "Role": "93,782 query proteins",
        },
        {
            "Source": "Article 35 gene counts",
            "Identifier": "MOCK1 and MOCK2 audited raw counts",
            "Release": "Article 35 frozen bundle",
            "SHA256": run_contract["abundance_sha256"],
            "Role": "read-weighted missingness",
        },
    ]
    for row in database_rows:
        source_rows.append(
            {
                "Source": f"eggNOG {row['Asset']}",
                "Identifier": row["SourceURL"],
                "Release": row["Release"],
                "SHA256": row["InstalledSHA256"],
                "Role": row["Notes"],
            }
        )
    write_tsv(stage / "source-audit.tsv", source_rows)

    frozen_contract = {
        "article": 36,
        "seed": SEED,
        "frozen_at": "2026-07-28",
        "paper_doi": "10.1186/s40168-022-01259-2",
        "catalog_genes": run_summary["catalog_genes"],
        "seed_ortholog_genes": run_summary["seed_ortholog_genes"],
        "annotation_states": run_summary["annotation_states"],
        "eggnog_mapper": run_contract["eggnog_mapper"],
        "eggnog_database": run_contract["eggnog_database"],
        "diamond": run_contract["search"]["diamond"],
        "go_evidence_primary": run_contract["annotation"]["go_evidence_primary"],
        "go_evidence_sensitivity": run_contract["annotation"]["go_evidence_sensitivity"],
        "large_assets_excluded": [
            "eggNOG annotation SQLite database",
            "eggNOG taxonomy database",
            "eggNOG DIAMOND protein index",
            "uncompressed input FAA",
            "raw DIAMOND hits table",
        ],
        "functional_annotations_are_predictions": True,
        "absence_is_not_gene_absence": True,
    }
    (stage / "frozen-contract.json").write_text(json.dumps(frozen_contract, indent=2) + "\n", encoding="utf-8")
    notice = """Article 36 frozen-data notice

This bundle contains complete eggNOG-mapper annotation rows for the 93,782
real catalog proteins, a GO all-evidence sensitivity output, compact gene-level
annotation states, KO/GO/COG summaries, exact commands, resource logs, software
locks, database identities, and source provenance. The roughly 52-GB installed
eggNOG 5.0.2 database, uncompressed query FAA, and raw DIAMOND hits are excluded.

All function calls are orthology-based predictions. A missing KO, GO, COG, EC,
or seed ortholog is database- and policy-specific missingness, not evidence that
the biological gene or activity is absent. Tutorial code is MIT; tutorial prose
is CC BY 4.0 under the repository LICENSE. The Delgado and Andersson article is
CC BY 4.0.
"""
    (stage / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    checksum_rows = []
    for path in sorted(item for item in stage.rglob("*") if item.is_file() and item.name != "file-checksums.sha256"):
        checksum_rows.append(f"{sha256(path)}  {path.relative_to(stage).as_posix()}")
    (stage / "file-checksums.sha256").write_text("\n".join(checksum_rows) + "\n", encoding="utf-8")
    if output.exists():
        shutil.rmtree(output)
    stage.replace(output)
    print(json.dumps({"output": str(output), "files": len(checksum_rows) + 1, "contract": frozen_contract}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
