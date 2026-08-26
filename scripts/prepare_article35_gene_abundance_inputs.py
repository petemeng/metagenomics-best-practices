#!/usr/bin/env python3
"""Prepare checksum-audited reads, gene reference, and annotation assets for Article 35."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
from pathlib import Path
from typing import Iterator, TextIO


SEED = 20260735
CATALOG_GENES = 93_782
UNIREF90_BYTES = 36_332_007_356
UNIREF90_SHA256 = "67a00a99ead2a00c737b4b9cb7e64ecc9085c2539bbd21a2d0c92913936995a8"
REACTION_MAP_SHA256 = "8419ce78a62ca9130914f2c347a9708111cedc7de52ba274659ce51ec7de7752"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--article30-raw-dir", type=Path, required=True)
    parser.add_argument("--uniref90-db", type=Path, required=True)
    parser.add_argument("--reaction-map", type=Path, required=True)
    parser.add_argument("--verify-large-db-sha256", action="store_true")
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
            elif header is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_checksums(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.strip():
            digest, relative = raw.split(maxsplit=1)
            values[relative.strip().lstrip("*")] = digest
    return values


def replace_symlink(link: Path, target: Path) -> None:
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.is_symlink() or link.exists():
        if link.is_symlink() and link.resolve() == target.resolve():
            return
        if link.is_dir() and not link.is_symlink():
            raise IsADirectoryError(link)
        link.unlink()
    link.symlink_to(target.resolve())


def prepare_catalog(root: Path, work: Path) -> tuple[list[dict], dict]:
    bundle = root / "data/small/34-nonredundant-gene-catalog-frozen"
    checksums = read_checksums(bundle / "file-checksums.sha256")
    fna_gz = bundle / "catalog/megahit-mix-primary.fna.gz"
    faa_gz = bundle / "catalog/megahit-mix-primary.faa.gz"
    metadata_gz = bundle / "primary-catalog-representatives.tsv.gz"
    for path in (fna_gz, faa_gz, metadata_gz):
        relative = str(path.relative_to(bundle))
        if sha256(path) != checksums.get(relative):
            raise ValueError(f"Article 34 checksum mismatch: {relative}")

    reference = work / "reference"
    reference.mkdir(parents=True, exist_ok=True)
    fna = reference / "megahit-mix-primary.fna"
    faa = reference / "megahit-mix-primary.faa"
    gff = reference / "megahit-mix-primary.gff"
    metadata = reference / "catalog-metadata.tsv"
    with gzip.open(fna_gz, "rb") as source, fna.open("wb") as target:
        shutil.copyfileobj(source, target, 8 * 1024 * 1024)
    with gzip.open(faa_gz, "rb") as source, faa.open("wb") as target:
        shutil.copyfileobj(source, target, 8 * 1024 * 1024)
    with gzip.open(metadata_gz, "rt", encoding="utf-8") as source, metadata.open("w", encoding="utf-8", newline="\n") as target:
        shutil.copyfileobj(source, target)

    proteins = {header.split()[0]: sequence.rstrip("*") for header, sequence in fasta_records(faa)}
    metadata_rows: dict[str, dict] = {}
    with metadata.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metadata_rows[row["RepresentativeID"]] = row
    catalog_rows: list[dict] = []
    seen: set[str] = set()
    total_bases = 0
    with gff.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##gff-version 3\n")
        for header, sequence in fasta_records(fna):
            gene = header.split()[0]
            if gene in seen:
                raise ValueError(f"Duplicate catalog identifier: {gene}")
            seen.add(gene)
            protein = proteins.get(gene)
            row = metadata_rows.get(gene)
            if protein is None or row is None:
                raise ValueError(f"Unpaired Article 34 catalog identifier: {gene}")
            if len(sequence) != int(row["NtLength"]) or len(protein) != int(row["AaLength"]):
                raise ValueError(f"Catalog length mismatch: {gene}")
            total_bases += len(sequence)
            attributes = f"ID={gene};Name={gene}"
            handle.write(f"{gene}\tArticle35\tCDS\t1\t{len(sequence)}\t.\t+\t0\t{attributes}\n")
            catalog_rows.append(
                {
                    "GeneID": gene,
                    "NtLength": len(sequence),
                    "AaLength": len(protein),
                    "Completeness": row["Completeness"],
                    "PartialCode": row["PartialCode"],
                    "RepresentativeOrigin": row["RepresentativeOrigin"],
                    "ClusterSize": row["ClusterSize"],
                }
            )
    if len(seen) != CATALOG_GENES or set(proteins) != seen or set(metadata_rows) != seen:
        raise ValueError(
            f"Expected {CATALOG_GENES} paired catalog genes; observed nt={len(seen)}, aa={len(proteins)}, metadata={len(metadata_rows)}"
        )
    write_tsv(reference / "catalog-audit.tsv", catalog_rows)
    summary = {
        "catalog_genes": len(seen),
        "catalog_bases": total_bases,
        "fna_gz_sha256": sha256(fna_gz),
        "faa_gz_sha256": sha256(faa_gz),
        "metadata_gz_sha256": sha256(metadata_gz),
        "fna_sha256": sha256(fna),
        "faa_sha256": sha256(faa),
        "gff_sha256": sha256(gff),
    }
    lineage = [
        {
            "Asset": "primary-gene-catalog-fna",
            "Role": "Bowtie2 nucleotide reference",
            "Source": str(fna_gz.relative_to(root)),
            "Bytes": fna_gz.stat().st_size,
            "SHA256": summary["fna_gz_sha256"],
            "IdentityGate": "Article34 frozen manifest",
        },
        {
            "Asset": "primary-gene-catalog-faa",
            "Role": "DIAMOND protein query",
            "Source": str(faa_gz.relative_to(root)),
            "Bytes": faa_gz.stat().st_size,
            "SHA256": summary["faa_gz_sha256"],
            "IdentityGate": "Article34 frozen manifest",
        },
        {
            "Asset": "primary-gene-catalog-metadata",
            "Role": "length, completeness, origin, and cluster-size metadata",
            "Source": str(metadata_gz.relative_to(root)),
            "Bytes": metadata_gz.stat().st_size,
            "SHA256": summary["metadata_gz_sha256"],
            "IdentityGate": "Article34 frozen manifest",
        },
    ]
    return lineage, summary


def prepare_reads(root: Path, raw: Path, work: Path) -> tuple[list[dict], dict]:
    run_summary = json.loads((root / "data/small/30-short-read-assembly-frozen/run-summary.json").read_text())
    specs = (
        ("MOCK1", "ERR9765746"),
        ("MOCK2", "ERR9765747"),
    )
    lineage: list[dict] = []
    summary: dict[str, dict] = {}
    for sample, run in specs:
        summary[sample] = {}
        for mate in ("R1", "R2"):
            source = raw / "clean" / f"{run}_clean_{mate}.fastq.gz"
            if not source.is_file() or source.stat().st_size == 0:
                raise FileNotFoundError(source)
            upstream = run_summary["clean_fastq_audit"][sample][mate]
            observed = sha256(source)
            if observed != upstream["CompressedSHA256"] or source.stat().st_size != upstream["CompressedBytes"]:
                raise ValueError(f"Article 30 clean FASTQ mismatch: {source}")
            link = work / "inputs" / f"{sample}_{mate}.fastq.gz"
            replace_symlink(link, source)
            summary[sample][mate] = {
                "records": int(upstream["Records"]),
                "bases": int(upstream["Bases"]),
                "bytes": source.stat().st_size,
                "sha256": observed,
                "path": str(source.resolve()),
            }
            lineage.append(
                {
                    "Asset": f"{sample}-{mate}",
                    "Role": "full clean read mapping input",
                    "Source": str(source.resolve()),
                    "Bytes": source.stat().st_size,
                    "SHA256": observed,
                    "IdentityGate": "Article30 frozen run-summary",
                }
            )
        if summary[sample]["R1"]["records"] != summary[sample]["R2"]["records"]:
            raise ValueError(f"Mate count mismatch for {sample}")
    return lineage, summary


def prepare_annotation_assets(root: Path, args: argparse.Namespace) -> tuple[list[dict], dict]:
    db = args.uniref90_db.resolve()
    reaction = args.reaction_map.resolve()
    if not db.is_file() or db.stat().st_size != UNIREF90_BYTES:
        raise ValueError(f"UniRef90 DIAMOND database byte mismatch: {db}")
    db_observed = sha256(db) if args.verify_large_db_sha256 else UNIREF90_SHA256
    db_gate = "observed full SHA-256" if args.verify_large_db_sha256 else "Article19 frozen upstream SHA-256"
    if db_observed != UNIREF90_SHA256:
        raise ValueError("UniRef90 DIAMOND database SHA-256 mismatch")
    if not reaction.is_file() or sha256(reaction) != REACTION_MAP_SHA256:
        raise ValueError("HUMAnN UniRef90-to-MetaCyc reaction mapping checksum mismatch")
    paper = root / "data/small/34-nonredundant-gene-catalog-frozen/sources/PMC9074274.fullTextXML"
    if not paper.is_file():
        raise FileNotFoundError(paper)
    rows = [
        {
            "Asset": "uniref90-diamond-db",
            "Role": "best-hit functional label",
            "Source": str(db),
            "Bytes": db.stat().st_size,
            "SHA256": UNIREF90_SHA256,
            "IdentityGate": db_gate,
        },
        {
            "Asset": "uniref90-metacyc-reaction-map",
            "Role": "one-to-many functional crosswalk",
            "Source": str(reaction),
            "Bytes": reaction.stat().st_size,
            "SHA256": REACTION_MAP_SHA256,
            "IdentityGate": "Article19 frozen regroup audit",
        },
        {
            "Asset": "delgado-paper-xml",
            "Role": "Fig. 3 anchor and historical mapping method",
            "Source": str(paper.relative_to(root)),
            "Bytes": paper.stat().st_size,
            "SHA256": sha256(paper),
            "IdentityGate": "Article34 frozen source",
        },
    ]
    return rows, {
        "uniref90_db": str(db),
        "uniref90_release": "v201901b",
        "uniref90_sha256": UNIREF90_SHA256,
        "reaction_map": str(reaction),
        "reaction_map_sha256": REACTION_MAP_SHA256,
        "paper_xml": str(paper),
        "paper_xml_sha256": sha256(paper),
    }


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    raw = args.article30_raw_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)

    catalog_lineage, catalog = prepare_catalog(root, work)
    read_lineage, reads = prepare_reads(root, raw, work)
    annotation_lineage, annotation = prepare_annotation_assets(root, args)
    lineage = catalog_lineage + read_lineage + annotation_lineage
    write_tsv(work / "input-lineage.tsv", lineage)

    contract = {
        "article": 35,
        "seed_namespace": SEED,
        "historical_branch": {
            "selection": "seqtk sample -s100, 10000 forward reads per sample",
            "mapping": "Bowtie2 --local, one reported primary alignment",
            "counting": "HTSeq -f bam -r pos -t CDS -i ID -s no -a 0",
        },
        "primary_branch": {
            "read_unit": "R1 and R2 mapped independently as reads",
            "mapping": "Bowtie2 --very-sensitive-local -k 2 --seed 20260735",
            "identity_formula": "1 - NM / sum(CIGAR M,I,D,=,X alignment columns)",
            "query_coverage_formula": "sum(CIGAR M,I,=,X query-aligned bases) / full query length",
            "filters": {"minimum_mapq": 10, "minimum_identity": 0.95, "minimum_query_coverage": 0.80},
            "normalizations": ["raw count", "CPM", "RPKM", "TPM"],
        },
        "functional_branch": {
            "annotation": "DIAMOND blastp best UniRef90 v201901b hit",
            "iterate": ["faster", "sensitive"],
            "max_target_seqs": 5,
            "masking": 1,
            "best_hit_order": "bitscore desc, evalue asc, identity desc, UniRef90 ID asc",
            "thresholds": {"minimum_identity_pct": 50, "minimum_query_coverage_pct": 80, "maximum_evalue": 1e-5},
            "crosswalk": "HUMAnN 3.9 MetaCyc reaction level4ec UniRef mapping",
            "primary_one_to_many_policy": "equal split with mass conservation",
            "sensitivity_policy": "copy full abundance to every mapped reaction",
        },
        "catalog": catalog,
        "reads": reads,
        "annotation_assets": annotation,
    }
    (work / "run-contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    (work / ".article35-inputs-complete").write_text("prepared\n", encoding="utf-8")
    print(json.dumps({"status": "prepared", "catalog_genes": catalog["catalog_genes"], "samples": list(reads), "lineage_rows": len(lineage)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
