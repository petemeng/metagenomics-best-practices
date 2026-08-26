#!/usr/bin/env python3
"""Create a compact, checksum-covered Article 59 evidence bundle."""

from __future__ import annotations

import argparse
import gzip
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


def gzip_copy(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as source_handle, target.open("wb") as raw_target:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_target, mtime=0) as gzip_handle:
            shutil.copyfileobj(source_handle, gzip_handle, length=8 * 1024 * 1024)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    target = args.output_dir.resolve()
    if not (work / ".article59-summary-complete").is_file():
        raise FileNotFoundError("Article 59 summary sentinel is missing")
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)

    for path in sorted((work / "summary").glob("*")):
        if path.is_file():
            copy(path, target / path.name)

    for name in (
        "database-audit.tsv",
        "dram-database-build.tsv",
        "dram-database-file-audit.tsv",
        "metabolic-database-build.tsv",
        "metabolic-database-file-audit.tsv",
        "metabolic-protein-id-audit.tsv",
        "kofam-compatibility-audit.tsv",
        "merops-normalization-audit.tsv",
        "input-mag-ledger.tsv",
        "truncation-ledger.tsv",
        "preparation-contract.json",
        "run-contract.json",
        "command-log.tsv",
        "tool-smoke.tsv",
        "dram-config.json",
        "dram-annotate-commands.txt",
        "metabolic-prodigal-commands.txt",
    ):
        copy(work / name, target / name)

    gzip_copy(work / "dram-annotation/annotations.tsv", target / "dram/annotations.tsv.gz")
    for name in ("genome_stats.tsv", "product.tsv", "metabolism_summary.xlsx", "distill.log"):
        copy(work / "dram-distill" / name, target / "dram" / name)
    copy(
        work / "database/dram-kofam-2026-06-01/module_step_form.tsv",
        target / "dram/module_step_form.tsv",
    )

    metabolic = work / "metabolic-output"
    for index in range(1, 7):
        name = f"METABOLIC_result_worksheet{index}.tsv"
        copy(
            metabolic / "METABOLIC_result_each_spreadsheet" / name,
            target / "metabolic/worksheets" / name,
        )
    for path in sorted((metabolic / "KEGG_identifier_result").glob("*.txt")):
        copy(path, target / "metabolic/KEGG_identifier_result" / path.name)
    for name in ("METABOLIC_result.xlsx", "METABOLIC_log.log", "METABOLIC_run.log"):
        copy(metabolic / name, target / "metabolic" / name)

    for path in sorted((work / "logs").glob("*")):
        if path.is_file():
            copy(path, target / "logs" / path.name)

    script_names = (
        "download_article59_metabolism_db.py",
        "prepare_article59_metabolism.py",
        "setup_article59_dram_database.py",
        "setup_article59_metabolic_database.py",
        "run_article59_metabolism.py",
        "summarize_article59_metabolism.py",
        "plot_article59_metabolism.py",
        "freeze_article59_metabolism.py",
        "validate_article59_metabolism.py",
        "article42_44_validation_utils.py",
    )
    for name in script_names:
        copy(root / "scripts" / name, target / "scripts" / name)

    for name in (
        "dram.yml",
        "dram-linux-64.lock",
        "metabolic.yml",
        "metabolic-linux-64.lock",
        "drep.yml",
        "drep-linux-64.lock",
    ):
        copy(root / "env" / name, target / "env" / name)
    copy(
        root / "data/small/59-metabolism-database-manifest.tsv",
        target / "database-manifest.tsv",
    )

    notice = """Article 59 frozen evidence bundle

The 24 primary inputs are real dRep representatives reconstructed from the
PRJEB52977 uneven mock study and classified with GTDB-Tk R232 in Articles 45
and 46. Four additional FASTAs are deterministic SGB_002 sequence-retention
sensitivities generated with seed 59002. Full MAG nucleotide/protein files are
not duplicated here; their source representative, decompressed SHA-256,
quality, taxonomy, retained bases and contig counts are in input-mag-ledger.tsv.

The 1.56-GB KOfam archive (7.2 GB extracted profiles), 21-GB DRAM merged/indexed
database, 458-MB MEROPS sequence library, dbCAN HMM database and transient
DRAM/METABOLIC work directories are external-only. Their immutable sources,
byte counts, checksums, releases, compatibility decisions and measured build
resources are retained. K18513 was requested by METABOLIC v4.0 but absent from
KOfam 2026-06-01; the 2,643 available profiles were used and the omission is
fail-closed in kofam-compatibility-audit.tsv.

The bundle contains compressed DRAM annotations, DRAM distillation tables,
all six METABOLIC worksheets, per-genome KO result/hit tables, commands,
environment locks, the globally unique genome-prefixed protein-ID audit,
resource logs and derived evidence tables. Gene hits,
module coverage and curated trait rules describe genomic potential; they are
not measurements of expression, metabolite turnover or flux.
"""
    (target / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    contract = {
        "article": 59,
        "created_from": str(work.relative_to(root)),
        "primary_real_mags": 24,
        "deterministic_sensitivity_genomes": 4,
        "large_databases_included": False,
        "full_mag_fastas_included": False,
        "compressed_dram_annotations_included": True,
        "metabolic_worksheets_included": 6,
        "seed": 59002,
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
    print(f"Frozen Article 59 bundle: {len(files)} payload files in {target}")


if __name__ == "__main__":
    main()
