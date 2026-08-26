#!/usr/bin/env python3
"""Freeze compact, checksum-locked Article 35 evidence without raw reads or databases."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


SEED = 20260735
SCRIPT_NAMES = (
    "download_article35_gene_abundance_sources.sh",
    "prepare_article35_gene_abundance_inputs.py",
    "parse_article35_sam.py",
    "run_article35_gene_abundance.py",
    "summarize_article35_gene_abundance.py",
    "freeze_article35_gene_abundance.py",
    "validate_article35_gene_abundance.py",
)
ROOT_FILES = ("input-lineage.tsv", "run-contract.json", "tool-versions.tsv", "command-log.tsv")
SUMMARY_FILES = (
    "legacy-mapping-audit.tsv",
    "mapping-policy-summary.tsv",
    "alignment-quality-histogram.tsv",
    "gene-functional-annotation.tsv.gz",
    "gene-abundance-long.tsv.gz",
    "unit-rank-audit.tsv",
    "normalization-audit.tsv",
    "gene-length-bin-summary.tsv",
    "gene-completeness-summary.tsv",
    "annotation-search-audit.tsv",
    "reaction-abundance-long.tsv",
    "functional-aggregation-audit.tsv",
    "resource-usage.tsv",
    "run-summary.json",
)
SELECTED_LOGS = (
    "build-bowtie2-index.time.txt",
    "legacy-map-MOCK1.stderr.log",
    "legacy-map-MOCK2.stderr.log",
    "legacy-htseq-MOCK1.time.txt",
    "legacy-htseq-MOCK2.time.txt",
    "audited-map-MOCK1.bowtie2.stderr.log",
    "audited-map-MOCK1.time.txt",
    "audited-map-MOCK2.bowtie2.stderr.log",
    "audited-map-MOCK2.time.txt",
    "annotate-uniref90.stderr.log",
    "annotate-uniref90.time.txt",
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


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    output = args.output_dir.resolve()
    if not (work / ".article35-summary-complete").is_file():
        raise FileNotFoundError("Article 35 summary is incomplete")

    stage = output.parent / f".{output.name}.staging"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)

    for name in ROOT_FILES:
        copy_required(work / name, stage / name)
    for name in SUMMARY_FILES:
        copy_required(work / "summary" / name, stage / name)
    for name in SELECTED_LOGS:
        copy_required(work / "logs" / name, stage / "logs" / name)
    for name in SCRIPT_NAMES:
        copy_required(root / "scripts" / name, stage / "scripts" / name)
    copy_required(root / "env/gene-abundance.yml", stage / "env/gene-abundance.yml")
    copy_required(root / "env/gene-abundance-linux-64.lock", stage / "env/gene-abundance-linux-64.lock")
    copy_required(
        root / "data/small/34-nonredundant-gene-catalog-frozen/data-NOTICE.txt",
        stage / "sources/article34-data-NOTICE.txt",
    )
    copy_required(
        root / "data/small/30-short-read-assembly-frozen/data-NOTICE.txt",
        stage / "sources/article30-data-NOTICE.txt",
    )
    copy_required(
        root / "data/small/19-humann3-frozen/database-audit.tsv",
        stage / "sources/article19-database-audit.tsv",
    )
    copy_required(
        root / "data/small/34-nonredundant-gene-catalog-frozen/sources/PMC9074274.fullTextXML",
        stage / "sources/PMC9074274.fullTextXML",
    )

    run_summary = json.loads((work / "summary/run-summary.json").read_text(encoding="utf-8"))
    run_contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))
    lineage = list(csv.DictReader((work / "input-lineage.tsv").open(encoding="utf-8"), delimiter="\t"))
    source_rows = [
        {
            "Source": row["Asset"],
            "Identifier": row["Source"],
            "VersionOrRelease": (
                "Article 34 frozen catalog" if row["Asset"].startswith("primary-gene-catalog")
                else "PRJEB52977 clean subset" if row["Asset"].startswith(("MOCK1", "MOCK2"))
                else "HUMAnN/UniRef90 v201901b" if "uniref90" in row["Asset"]
                else "published 2022-05-06"
            ),
            "Bytes": row["Bytes"],
            "SHA256": row["SHA256"],
            "IdentityGate": row["IdentityGate"],
            "Role": row["Role"],
        }
        for row in lineage
    ]
    write_tsv(stage / "source-audit.tsv", source_rows)

    notice = """Article 35 frozen-data notice

The compact bundle contains complete gene-level counts and normalized units,
best-hit functional annotations, reaction-level summaries, mapping ledgers,
selected logs, exact software locks, and source identities. It intentionally
does not contain raw/clean FASTQ, Bowtie2 indexes, SAM/BAM files, temporary
DIAMOND files, or the 36.3-GB UniRef90 database.

The two mock communities are technical calibration inputs. Their reads also
contributed to the Article 34 catalog, so mapping and detection are not an
independent predictive evaluation. The Delgado and Andersson XML is CC BY 4.0.
Original tutorial code is MIT and tutorial prose is CC BY 4.0 under the root
LICENSE.
"""
    (stage / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    contract = {
        "article": 35,
        "seed": SEED,
        "frozen_at": "2026-07-28",
        "paper_doi": "10.1186/s40168-022-01259-2",
        "catalog_genes": run_summary["catalog_genes"],
        "samples": run_summary["samples"],
        "input_reads": run_summary["input_reads"],
        "main_assigned_reads": run_summary["main_assigned_reads"],
        "uniref90_hit_genes": run_summary["uniref90_hit_genes"],
        "uniref90_release": run_contract["annotation_assets"]["uniref90_release"],
        "uniref90_db_sha256": run_contract["annotation_assets"]["uniref90_sha256"],
        "reaction_map_sha256": run_contract["annotation_assets"]["reaction_map_sha256"],
        "read_unit": "R1 and R2 mapped independently as reads",
        "primary_mapping_policy": run_contract["primary_branch"],
        "functional_policy": run_contract["functional_branch"],
        "large_assets_excluded": ["FASTQ", "Bowtie2 index", "SAM/BAM", "UniRef90 DIAMOND database"],
    }
    (stage / "frozen-contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")

    checksum_rows: list[str] = []
    for path in sorted(item for item in stage.rglob("*") if item.is_file() and item.name != "file-checksums.sha256"):
        checksum_rows.append(f"{sha256(path)}  {path.relative_to(stage)}")
    (stage / "file-checksums.sha256").write_text("\n".join(checksum_rows) + "\n", encoding="utf-8")

    if output.exists():
        shutil.rmtree(output)
    stage.replace(output)
    print(json.dumps({"output": str(output), "files": len(checksum_rows) + 1, "contract": contract}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
