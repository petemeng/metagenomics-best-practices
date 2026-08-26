#!/usr/bin/env python3
"""Prepare checksum-locked Zymo inputs for Article 58."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import BinaryIO, Iterator


EXPECTED_SPECIES = (
    "Bacillus_subtilis",
    "Cryptococcus_neoformans",
    "Enterococcus_faecalis",
    "Escherichia_coli",
    "Lactobacillus_fermentum",
    "Listeria_monocytogenes",
    "Pseudomonas_aeruginosa",
    "Saccharomyces_cerevisiae",
    "Salmonella_enterica",
    "Staphylococcus_aureus",
)
EUKARYOTES = {"Cryptococcus_neoformans", "Saccharomyces_cerevisiae"}
FRAGMENT_LENGTHS = (3000, 5000, 10000, 20000)
FRAGMENTS_PER_SPECIES = 80


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--eukdetect-pairs", type=int, default=1_000_000)
    parser.add_argument("--assembly-pairs", type=int, default=20_000_000)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header = line[1:]
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                parts.append(line)
    if header is not None:
        yield header, "".join(parts).upper()


def safe_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "sequence"


def read_record(handle: BinaryIO, mate: int, index: int) -> tuple[list[bytes], int]:
    record = [handle.readline() for _ in range(4)]
    if any(line == b"" for line in record):
        raise RuntimeError(
            f"Compressed byte prefix ended before pair {index:,}, mate R{mate}; "
            "increase the locked byte range"
        )
    if not record[0].startswith(b"@") or not record[2].startswith(b"+"):
        raise RuntimeError(f"Malformed FASTQ record at pair {index:,}, mate R{mate}")
    sequence = record[1].rstrip(b"\r\n")
    quality = record[3].rstrip(b"\r\n")
    if len(sequence) != len(quality):
        raise RuntimeError(f"Sequence/quality length mismatch at pair {index:,}, R{mate}")
    return record, len(sequence)


def read_name(header: bytes) -> bytes:
    token = header[1:].split(None, 1)[0]
    if token.endswith(b"/1") or token.endswith(b"/2"):
        token = token[:-2]
    return token


def materialize_fastq_prefixes(
    r1_path: Path,
    r2_path: Path,
    output_dir: Path,
    eukdetect_pairs: int,
    assembly_pairs: int,
) -> dict[str, object]:
    if eukdetect_pairs > assembly_pairs:
        raise ValueError("eukdetect_pairs cannot exceed assembly_pairs")
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "assembly_r1": output_dir / f"SRR12324253.first-{assembly_pairs}-pairs.R1.fastq",
        "assembly_r2": output_dir / f"SRR12324253.first-{assembly_pairs}-pairs.R2.fastq",
        "eukdetect_r1": output_dir / f"SRR12324253.first-{eukdetect_pairs}-pairs.R1.fastq",
        "eukdetect_r2": output_dir / f"SRR12324253.first-{eukdetect_pairs}-pairs.R2.fastq",
    }
    base_counts = {"assembly_r1": 0, "assembly_r2": 0, "eukdetect_r1": 0, "eukdetect_r2": 0}
    with r1_path.open("rb") as raw1, r2_path.open("rb") as raw2:
        with gzip.GzipFile(fileobj=raw1, mode="rb") as in1, gzip.GzipFile(
            fileobj=raw2, mode="rb"
        ) as in2:
            with paths["assembly_r1"].open("wb") as asm1, paths["assembly_r2"].open(
                "wb"
            ) as asm2, paths["eukdetect_r1"].open("wb") as det1, paths[
                "eukdetect_r2"
            ].open("wb") as det2:
                for index in range(1, assembly_pairs + 1):
                    rec1, len1 = read_record(in1, 1, index)
                    rec2, len2 = read_record(in2, 2, index)
                    if read_name(rec1[0]) != read_name(rec2[0]):
                        raise RuntimeError(f"Mate identifiers diverge at pair {index:,}")
                    asm1.writelines(rec1)
                    asm2.writelines(rec2)
                    base_counts["assembly_r1"] += len1
                    base_counts["assembly_r2"] += len2
                    if index <= eukdetect_pairs:
                        det1.writelines(rec1)
                        det2.writelines(rec2)
                        base_counts["eukdetect_r1"] += len1
                        base_counts["eukdetect_r2"] += len2
    return {
        "paths": {key: str(path) for key, path in paths.items()},
        "base_counts": base_counts,
        "eukdetect_total_bases": base_counts["eukdetect_r1"] + base_counts["eukdetect_r2"],
        "assembly_total_bases": base_counts["assembly_r1"] + base_counts["assembly_r2"],
    }


def evenly_spaced_indices(total: int, requested: int) -> list[int]:
    if total < requested:
        raise RuntimeError(f"Only {total} eligible fragments; need {requested}")
    if requested == 1:
        return [0]
    indices = [(i * (total - 1)) // (requested - 1) for i in range(requested)]
    if len(set(indices)) != requested:
        raise RuntimeError("Deterministic fragment index selection produced duplicates")
    return indices


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    raw = args.raw_dir.resolve()
    source_manifest = args.source_manifest.resolve()
    marker = work / ".article58-inputs-complete"
    if marker.is_file() and not args.force:
        print(f"Article 58 inputs already prepared: {work}")
        return
    work.mkdir(parents=True, exist_ok=True)

    rows = read_tsv(source_manifest)
    required_roles = {"zymo-reference-v2", "srr12324253-r1-fastq", "srr12324253-r2-fastq"}
    observed_roles = {row["Role"] for row in rows}
    if observed_roles != required_roles:
        raise RuntimeError(f"Unexpected source roles: {sorted(observed_roles)}")
    asset_audit: list[dict[str, object]] = []
    by_role: dict[str, Path] = {}
    for row in rows:
        path = root / row["LocalPath"]
        if not path.is_file():
            raise FileNotFoundError(path)
        observed_bytes = path.stat().st_size
        algorithm = row["ChecksumAlgorithm"]
        observed_checksum = digest(path, algorithm)
        passed = observed_bytes == int(row["LocalBytes"]) and observed_checksum == row["LocalChecksum"]
        asset_audit.append(
            {
                "Role": row["Role"],
                "LocalPath": row["LocalPath"],
                "ObservedBytes": observed_bytes,
                "ExpectedBytes": row["LocalBytes"],
                "ChecksumAlgorithm": algorithm,
                "ObservedChecksum": observed_checksum,
                "ExpectedChecksum": row["LocalChecksum"],
                "ChecksumPass": str(passed).lower(),
                "RemoteBytes": row["RemoteBytes"],
                "RemoteMD5": row["RemoteMD5"],
                "ByteRange": row["ByteRange"],
                "URL": row["URL"],
            }
        )
        if not passed:
            raise RuntimeError(f"Source checksum gate failed for {path}")
        by_role[row["Role"]] = path
    write_tsv(work / "asset-check-audit.tsv", asset_audit)

    reference_dir = work / "references" / "ZymoBIOMICS.STD.refseq.v2"
    reference_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(by_role["zymo-reference-v2"]) as archive:
        genome_members = sorted(
            name for name in archive.namelist() if "/Genomes/" in name and name.endswith(".fasta")
        )
        if len(genome_members) != 10:
            raise RuntimeError(f"Expected 10 Zymo genomes, observed {len(genome_members)}")
        for member in genome_members:
            target = reference_dir / Path(member).name
            with archive.open(member) as source, target.open("wb") as destination:
                for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                    destination.write(block)

    genomes: dict[str, list[tuple[str, str]]] = {}
    reference_ledger: list[dict[str, object]] = []
    for path in sorted(reference_dir.glob("*.fasta")):
        species = next((name for name in EXPECTED_SPECIES if path.name.startswith(name)), None)
        if species is None:
            raise RuntimeError(f"Cannot infer species from {path.name}")
        records = list(parse_fasta(path))
        if not records:
            raise RuntimeError(f"No reference sequences in {path}")
        genomes[species] = records
        domain = "Eukaryote" if species in EUKARYOTES else "Prokaryote"
        for index, (header, sequence) in enumerate(records, start=1):
            reference_ledger.append(
                {
                    "ReferenceID": f"ZYMOV2|{species}|{index:05d}|{safe_token(header.split()[0])}",
                    "Species": species.replace("_", " "),
                    "DomainTruth": domain,
                    "SourceFile": path.name,
                    "SourceHeader": header,
                    "LengthBp": len(sequence),
                }
            )
    if set(genomes) != set(EXPECTED_SPECIES):
        raise RuntimeError("Zymo reference species inventory is incomplete")
    write_tsv(work / "reference-sequence-ledger.tsv", reference_ledger)

    combined = work / "references" / "zymo-v2-all-references.fna"
    with combined.open("w", encoding="utf-8") as handle:
        for species in EXPECTED_SPECIES:
            for index, (header, sequence) in enumerate(genomes[species], start=1):
                ref_id = f"ZYMOV2|{species}|{index:05d}|{safe_token(header.split()[0])}"
                handle.write(f">{ref_id}\n{sequence}\n")

    benchmark_dir = work / "benchmarks"
    benchmark_dir.mkdir(parents=True, exist_ok=True)
    fragment_ledger: list[dict[str, object]] = []
    for fragment_length in FRAGMENT_LENGTHS:
        fasta_path = benchmark_dir / f"eukrep-fragments-{fragment_length}.fna"
        with fasta_path.open("w", encoding="utf-8") as fasta:
            for species in EXPECTED_SPECIES:
                candidates: list[tuple[str, int, str]] = []
                for contig_index, (header, sequence) in enumerate(genomes[species], start=1):
                    for start in range(0, len(sequence) - fragment_length + 1, fragment_length):
                        candidates.append((f"{contig_index:05d}|{safe_token(header.split()[0])}", start, sequence[start : start + fragment_length]))
                selected = evenly_spaced_indices(len(candidates), FRAGMENTS_PER_SPECIES)
                domain = "Eukaryote" if species in EUKARYOTES else "Prokaryote"
                for rank, candidate_index in enumerate(selected, start=1):
                    source_id, start, sequence = candidates[candidate_index]
                    fragment_id = f"FRAG|{fragment_length}|{species}|{domain}|{rank:03d}|{source_id}|{start + 1}"
                    fasta.write(f">{fragment_id}\n{sequence}\n")
                    fragment_ledger.append(
                        {
                            "FragmentID": fragment_id,
                            "FragmentLength": fragment_length,
                            "Species": species.replace("_", " "),
                            "DomainTruth": domain,
                            "WithinSpeciesRank": rank,
                            "SourceSequence": source_id,
                            "Start1": start + 1,
                            "End1": start + fragment_length,
                            "Selection": "80 evenly spaced non-overlapping exact-length fragments per species",
                        }
                    )
    write_tsv(work / "eukrep-fragment-ledger.tsv", fragment_ledger)

    fastq = materialize_fastq_prefixes(
        by_role["srr12324253-r1-fastq"],
        by_role["srr12324253-r2-fastq"],
        work / "inputs",
        args.eukdetect_pairs,
        args.assembly_pairs,
    )
    fastq_audit: list[dict[str, object]] = []
    for role, value in fastq["paths"].items():
        path = Path(str(value))
        fastq_audit.append(
            {
                "Role": role,
                "Path": str(path.relative_to(root)),
                "Bytes": path.stat().st_size,
                "SHA256": digest(path, "sha256"),
                "Bases": fastq["base_counts"][role],
                "Pairs": args.eukdetect_pairs if role.startswith("eukdetect") else args.assembly_pairs,
            }
        )
    write_tsv(work / "prepared-fastq-audit.tsv", fastq_audit)
    write_tsv(
        work / "library-sizes-pre-qc.tsv",
        [{"Sample": "Zymo_D6300", "TotalBases": fastq["eukdetect_total_bases"]}],
    )
    run_contract = {
        "article": 58,
        "seed": 20260758,
        "random_output_requested": False,
        "run_accession": "SRR12324253",
        "bioproject": "PRJNA648136",
        "reference_doi": "10.5281/zenodo.3935737",
        "reference_genomes": 10,
        "expected_eukaryotes": 2,
        "fragment_lengths": list(FRAGMENT_LENGTHS),
        "fragments_per_species_length": FRAGMENTS_PER_SPECIES,
        "eukdetect_pairs": args.eukdetect_pairs,
        "assembly_pairs": args.assembly_pairs,
        "eukdetect_total_bases_pre_qc": fastq["eukdetect_total_bases"],
        "assembly_total_bases_pre_qc": fastq["assembly_total_bases"],
    }
    (work / "run-contract.json").write_text(
        json.dumps(run_contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    marker.write_text("verified\n", encoding="utf-8")
    print(f"Article 58 inputs prepared: {work}")


if __name__ == "__main__":
    main()
