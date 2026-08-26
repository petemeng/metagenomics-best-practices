#!/usr/bin/env python3
"""Freeze compact, checksum-locked Article 34 evidence and the primary catalog."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
from pathlib import Path
from typing import Iterator, TextIO


SEED = 20260734
SCRIPT_NAMES = (
    "download_article34_gene_catalog_sources.sh",
    "prepare_article34_gene_catalog_inputs.py",
    "run_article34_gene_catalog.py",
    "summarize_article34_gene_catalog.py",
    "freeze_article34_gene_catalog.py",
    "validate_article34_gene_catalog.py",
)
SUMMARY_FILES = (
    "input-lineage.tsv",
    "truth-genomes.tsv",
    "prepare-audit.json",
    "run-contract.json",
    "tool-versions.tsv",
    "command-log.tsv",
    "gene-prediction-summary.tsv",
    "catalog-summary.tsv",
    "truth-audit-summary.tsv",
    "per-genome-recovery.tsv",
    "mix-origin-summary.tsv",
    "cluster-size-bins.tsv",
    "gene-length-histogram.tsv",
    "method-agreement.tsv",
    "resource-usage.tsv",
    "run-summary.json",
    "primary-catalog-membership.tsv.gz",
    "primary-catalog-representatives.tsv.gz",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--benchmark-repo", type=Path, required=True)
    parser.add_argument("--paper-xml", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_record(handle: TextIO, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def copy_required(source: Path, destination: Path) -> None:
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def gzip_deterministic(source: Path, destination: Path) -> None:
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as src, destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as dst:
            shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)


def write_head(source: Path, destination: Path, records: int = 12) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        for index, (header, sequence) in enumerate(fasta_records(source)):
            if index >= records:
                break
            write_record(handle, header, sequence)


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    repo = args.benchmark_repo.resolve()
    paper = args.paper_xml.resolve()
    output = args.output_dir.resolve()
    if not (work / ".article34-summary-complete").is_file():
        raise FileNotFoundError("Article 34 summary is incomplete")
    if "10.1186/s40168-022-01259-2" not in paper.read_text(encoding="utf-8"):
        raise ValueError("The supplied paper XML does not contain the locked DOI")

    stage = output.parent / f".{output.name}.staging"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)

    for name in SUMMARY_FILES:
        copy_required(work / "summary" / name, stage / name)

    primary_dir = work / "catalogs/megahit/mix-primary"
    gzip_deterministic(primary_dir / "representatives.faa", stage / "catalog/megahit-mix-primary.faa.gz")
    gzip_deterministic(primary_dir / "representatives.fna", stage / "catalog/megahit-mix-primary.fna.gz")
    write_head(primary_dir / "representatives.faa", stage / "catalog/megahit-mix-primary.head.faa")
    write_head(primary_dir / "representatives.fna", stage / "catalog/megahit-mix-primary.head.fna")

    selected_logs = (
        "prodigal-megahit-m1.time.txt",
        "megahit-individual-primary-cluster.time.txt",
        "megahit-mix-primary-cluster.time.txt",
        "megahit-mix-cdhit.time.txt",
        "truth-nr-primary-cluster.time.txt",
        "truth-megahit-mix-primary-catalog-to-truth-recovery-search.time.txt",
        "truth-megahit-mix-primary-catalog-to-truth-support-search.time.txt",
    )
    for name in selected_logs:
        copy_required(work / "logs" / name, stage / "logs" / name)

    for name in SCRIPT_NAMES:
        copy_required(root / "scripts" / name, stage / "scripts" / name)
    copy_required(root / "env/gene-catalog.yml", stage / "env/gene-catalog.yml")
    copy_required(root / "env/gene-catalog-linux-64.lock", stage / "env/gene-catalog-linux-64.lock")
    copy_required(repo / "LICENSE", stage / "sources/benchmark-repository-LICENSE.txt")
    copy_required(repo / "README.md", stage / "sources/benchmark-repository-README.md")
    copy_required(paper, stage / "sources/PMC9074274.fullTextXML")
    copy_required(
        root / "data/small/30-short-read-assembly-frozen/data-NOTICE.txt",
        stage / "sources/article30-data-NOTICE.txt",
    )

    run_summary = json.loads((work / "summary/run-summary.json").read_text(encoding="utf-8"))
    source_rows = [
        {
            "Source": "Delgado and Andersson 2022 full text",
            "Identifier": "doi:10.1186/s40168-022-01259-2; PMCID:PMC9074274",
            "URL": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9074274/fullTextXML",
            "VersionOrCommit": "published 2022-05-06",
            "LocalFile": "sources/PMC9074274.fullTextXML",
            "SHA256": sha256(paper),
            "Role": "anchor design and paper-compatible Prodigal/MMseqs2 parameters",
        },
        {
            "Source": "Meslier et al. benchmark repository",
            "Identifier": "PRJEB52977 support code and exact references",
            "URL": "https://forgemia.inra.fr/metagenopolis/benchmark_mock",
            "VersionOrCommit": "a429a3724d4593f35b8d7323b20252a6be90e1cd",
            "LocalFile": "sources/benchmark-repository-README.md",
            "SHA256": sha256(repo / "README.md"),
            "Role": "87 exact MOCK2 genomes and Supplementary Table S1 lineage",
        },
        {
            "Source": "Article 30 frozen assemblies",
            "Identifier": "ERR9765746 + ERR9765747; six >=1 kb FASTAs",
            "URL": "https://www.ebi.ac.uk/ena/browser/view/PRJEB52977",
            "VersionOrCommit": "MEGAHIT 1.2.9; metaSPAdes 4.3.0",
            "LocalFile": "input-lineage.tsv",
            "SHA256": sha256(work / "summary/input-lineage.tsv"),
            "Role": "individual/co-assembly inputs with upstream checksum verification",
        },
    ]
    write_tsv(stage / "source-audit.tsv", source_rows)

    notice = """Article 34 frozen-data notice

The six assembly FASTAs are not duplicated here; their identities and upstream
checksums are recorded in input-lineage.tsv. The compact bundle includes the
primary MEGAHIT mix protein/nucleotide catalog, expanded cluster membership,
summary tables, selected execution logs, exact environment locks, and source
licenses. Raw FASTQ and the 87 full reference genomes remain outside Git.

The benchmark repository is redistributed under its included license. The
Delgado and Andersson article XML is CC BY 4.0. Original tutorial code is MIT;
tutorial prose is CC BY 4.0 under the repository LICENSE.
"""
    (stage / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    contract = {
        "article": 34,
        "seed": SEED,
        "frozen_at": "2026-07-28",
        "paper_doi": "10.1186/s40168-022-01259-2",
        "benchmark_commit": "a429a3724d4593f35b8d7323b20252a6be90e1cd",
        "assembly_branches": 6,
        "truth_genomes": 87,
        "primary_catalog": run_summary["primary_catalog"],
        "catalog_payloads": [
            "catalog/megahit-mix-primary.faa.gz",
            "catalog/megahit-mix-primary.fna.gz",
            "primary-catalog-membership.tsv.gz",
            "primary-catalog-representatives.tsv.gz",
        ],
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
