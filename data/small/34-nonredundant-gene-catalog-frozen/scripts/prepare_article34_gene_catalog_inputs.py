#!/usr/bin/env python3
"""Prepare checksum-audited assemblies and exact mock references for Article 34."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import subprocess
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator, TextIO

import openpyxl


SEED = 20260734

ALIASES = {
    "Chlorobaculum_tepidum_TLS": "Chlorobium_tepidum_TLS",
    "Trichormus_sp._PCC_7120": "Nostoc_sp._PCC_7120",
    "Sulfurisphaera_tokodaii_str._7": "Sulfolobus_tokodaii_str._7",
    "Phocaeicola_vulgatus_ATCC_8482": "Bacteroides_vulgatus_ATCC_8482",
    "Micromonospora_tropica_CNB-440": "Salinispora_tropica_CNB-440",
    "Micromonospora_arenicola_CNS-205": "Salinispora_arenicola_CNS-205",
    "Desulfovibrio_mercurii_ND_132": "Desulfovibrio_desulfuricans_ND132",
    "Bordetella_pertussis_RB50": "Bordetella_bronchiseptica_RB50",
    "Nitratidesulfovibrio_vulgaris_Hildenborough": "Desulfovibrio_vulgaris_Hildenborough",
    "Cereibacter_A_sphaeroides_ATCC_17029": "Rhodobacter_sphaeroides_ATCC_17029",
    "Pauljensenia_odontolyticus_ATCC_17982": "Actinomyces_odontolyticus_ATCC_17982",
    "Cutibacterium_acnes_ATCC_11828": "Propionibacterium_acnes_ATCC_11828",
}

BRANCHES = (
    {
        "Branch": "megahit-m1",
        "Assembler": "MEGAHIT",
        "Strategy": "Individual",
        "Mock": "MOCK1",
        "InputPairs": 1_999_853,
        "RelativePath": "contigs/megahit-single-MOCK1.ge1000.fna.gz",
    },
    {
        "Branch": "megahit-m2",
        "Assembler": "MEGAHIT",
        "Strategy": "Individual",
        "Mock": "MOCK2",
        "InputPairs": 1_999_888,
        "RelativePath": "contigs/megahit-single-MOCK2.ge1000.fna.gz",
    },
    {
        "Branch": "megahit-co",
        "Assembler": "MEGAHIT",
        "Strategy": "Co-assembly",
        "Mock": "MOCK1+MOCK2",
        "InputPairs": 3_999_741,
        "RelativePath": "contigs/megahit-coassembly.ge1000.fna.gz",
    },
    {
        "Branch": "metaspades-m1",
        "Assembler": "metaSPAdes",
        "Strategy": "Individual",
        "Mock": "MOCK1",
        "InputPairs": 1_999_853,
        "RelativePath": "contigs/metaspades-single-MOCK1.ge1000.fna.gz",
    },
    {
        "Branch": "metaspades-m2",
        "Assembler": "metaSPAdes",
        "Strategy": "Individual",
        "Mock": "MOCK2",
        "InputPairs": 1_999_888,
        "RelativePath": "contigs/metaspades-single-MOCK2.ge1000.fna.gz",
    },
    {
        "Branch": "metaspades-co",
        "Assembler": "metaSPAdes",
        "Strategy": "Co-assembly",
        "Mock": "MOCK1+MOCK2",
        "InputPairs": 3_999_741,
        "RelativePath": "contigs/metaspades-coassembly.ge1000.fna.gz",
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--benchmark-repo", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:]
                chunks = []
            elif name is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def write_record(handle: TextIO, name: str, sequence: str) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def deterministic_gzip_text(path: Path) -> TextIO:
    raw = path.open("wb")
    gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(gz, encoding="utf-8", newline="")


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = deterministic_gzip_text(path) if path.suffix == ".gz" else path.open("w", encoding="utf-8", newline="")
    with handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_checksums(path: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        digest, relative = raw.split(maxsplit=1)
        checksums[relative.strip().lstrip("*")] = digest
    return checksums


def canonical_digest(sequence_hashes: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in sorted(sequence_hashes):
        digest.update(value.encode())
        digest.update(b"\n")
    return digest.hexdigest()


def normalize_assemblies(root: Path, work: Path) -> tuple[list[dict], list[dict]]:
    bundle = root / "data/small/30-short-read-assembly-frozen"
    manifest_path = bundle / "file-checksums.sha256"
    expected = read_checksums(manifest_path)
    lineage: list[dict] = []
    contig_rows: list[dict] = []
    out_dir = work / "assemblies"
    out_dir.mkdir(parents=True, exist_ok=True)

    for spec in BRANCHES:
        source = bundle / spec["RelativePath"]
        observed_sha = sha256(source)
        expected_sha = expected.get(spec["RelativePath"])
        if expected_sha is None or observed_sha != expected_sha:
            raise ValueError(f"Upstream checksum mismatch: {source}")
        output = out_dir / f"{spec['Branch']}.fna"
        sequence_hashes: list[str] = []
        total_bases = 0
        with output.open("w", encoding="utf-8", newline="\n") as handle:
            for index, (original_header, sequence) in enumerate(fasta_records(source), start=1):
                if len(sequence) < 1_000:
                    raise ValueError(f"{source} contains a contig below the locked 1 kb threshold")
                normalized_id = f"{spec['Branch']}__c{index:06d}"
                sequence_sha = hashlib.sha256(sequence.encode()).hexdigest()
                sequence_hashes.append(sequence_sha)
                total_bases += len(sequence)
                write_record(handle, normalized_id, sequence)
                contig_rows.append(
                    {
                        "Branch": spec["Branch"],
                        "NormalizedContigID": normalized_id,
                        "OriginalHeader": original_header,
                        "LengthBp": len(sequence),
                        "SequenceSHA256": sequence_sha,
                    }
                )
        lineage.append(
            {
                **spec,
                "SourcePath": str(source.relative_to(root)),
                "SourceCompressedBytes": source.stat().st_size,
                "SourceCompressedSHA256": observed_sha,
                "UpstreamManifest": str(manifest_path.relative_to(root)),
                "UpstreamChecksumVerified": "yes",
                "NormalizedPath": str(output.relative_to(work)),
                "Contigs": len(sequence_hashes),
                "AssemblyBases": total_bases,
                "MinimumContigBp": 1_000,
                "CanonicalSequenceSHA256": canonical_digest(sequence_hashes),
                "NormalizedFASTA_SHA256": sha256(output),
            }
        )
    return lineage, contig_rows


def workbook_rows(repo: Path) -> tuple[list, dict[str, tuple]]:
    supplement = repo / "script_r/Supplementary_Table_S1.xlsx"
    workbook = openpyxl.load_workbook(supplement, read_only=True, data_only=True)
    values = list(workbook.active.iter_rows(values_only=True))
    return list(values[0]), {str(row[3]).strip(): row for row in values[1:] if row[3]}


def abundance(value: object) -> float:
    return 0.0 if value in (None, "NP") else float(value)


def normalize_truth(repo: Path, work: Path) -> tuple[list[dict], list[dict], dict]:
    headers, accession_rows = workbook_rows(repo)
    m1 = [line.strip() for line in (repo / "profiling/MOCK_001.list").read_text().splitlines() if line.strip()]
    m2 = [line.strip() for line in (repo / "profiling/MOCK_002.list").read_text().splitlines() if line.strip()]
    if len(m1) != 71 or len(m2) != 87 or not set(m1) < set(m2):
        raise ValueError("Expected MOCK1=71, MOCK2=87, with MOCK1 a strict subset of MOCK2")

    genome_dir = repo / "reference/all_genomes_listed"
    output = work / "truth/mock2-exact-genomes.fna"
    output.parent.mkdir(parents=True, exist_ok=True)
    genome_rows: list[dict] = []
    contig_rows: list[dict] = []
    aggregate_hashes: Counter[str] = Counter()

    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in m2:
            row = accession_rows.get(accession)
            if row is None:
                raise ValueError(f"No Supplementary Table S1 row for {accession}")
            current_label = str(row[2]).strip()
            repository_label = ALIASES.get(current_label, current_label)
            source = genome_dir / f"{repository_label}.fna.gz"
            if not source.is_file():
                raise FileNotFoundError(source)
            a1 = abundance(row[14])
            a2 = abundance(row[15])
            genome_bases = 0
            genome_contigs = 0
            source_hashes: list[str] = []
            for index, (original_header, sequence) in enumerate(fasta_records(source), start=1):
                normalized_id = f"truth__{accession}__c{index:03d}"
                sequence_sha = hashlib.sha256(sequence.encode()).hexdigest()
                aggregate_hashes[sequence_sha] += 1
                source_hashes.append(sequence_sha)
                genome_bases += len(sequence)
                genome_contigs += 1
                write_record(handle, normalized_id, sequence)
                contig_rows.append(
                    {
                        "GenBankAssembly": accession,
                        "Reference": current_label,
                        "RepositoryLabel": repository_label,
                        "NormalizedContigID": normalized_id,
                        "OriginalHeader": original_header,
                        "LengthBp": len(sequence),
                        "SequenceSHA256": sequence_sha,
                    }
                )
            genome_rows.append(
                {
                    "GenBankAssembly": accession,
                    "Reference": current_label,
                    "RepositoryLabel": repository_label,
                    "InMOCK1": "yes" if accession in set(m1) else "no",
                    "InMOCK2": "yes",
                    "ExpectedAbundanceMOCK1Pct": f"{a1:.12g}",
                    "ExpectedAbundanceMOCK2Pct": f"{a2:.12g}",
                    "ExpectedAbundanceMeanPct": f"{(a1 + a2) / 2:.12g}",
                    "ReferenceContigs": genome_contigs,
                    "ReferenceBases": genome_bases,
                    "CompressedBytes": source.stat().st_size,
                    "CompressedSHA256": sha256(source),
                    "CanonicalSequenceSHA256": canonical_digest(source_hashes),
                    "TaxonomyRenameAlias": "yes" if current_label != repository_label else "no",
                }
            )

    combined = repo / "reference/MOCK_002.fasta.gz"
    combined_hashes = Counter(hashlib.sha256(sequence.encode()).hexdigest() for _, sequence in fasta_records(combined))
    if aggregate_hashes != combined_hashes:
        raise ValueError("The 87 normalized genomes do not equal repository MOCK_002.fasta.gz")

    audit = {
        "seed": SEED,
        "mock1_genomes": len(m1),
        "mock2_genomes": len(m2),
        "mock1_is_strict_subset_of_mock2": True,
        "truth_fasta": str(output.relative_to(work)),
        "truth_fasta_sha256": sha256(output),
        "truth_contigs": len(contig_rows),
        "truth_bases": sum(int(row["ReferenceBases"]) for row in genome_rows),
        "combined_mock2_sha256": sha256(combined),
        "supplement_s1_sha256": sha256(repo / "script_r/Supplementary_Table_S1.xlsx"),
        "supplement_headers": [str(item) for item in headers[:17]],
        "benchmark_commit": subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            check=True,
            text=True,
            capture_output=True,
        ).stdout.strip(),
    }
    return genome_rows, contig_rows, audit


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    repo = args.benchmark_repo.resolve()
    work = args.work_dir.resolve()
    summary = work / "summary"
    summary.mkdir(parents=True, exist_ok=True)

    lineage, contig_rows = normalize_assemblies(root, work)
    truth_genomes, truth_contigs, truth_audit = normalize_truth(repo, work)

    write_tsv(summary / "input-lineage.tsv", lineage)
    write_tsv(summary / "contig-map.tsv.gz", contig_rows)
    write_tsv(summary / "truth-genomes.tsv", truth_genomes)
    write_tsv(summary / "truth-contig-map.tsv.gz", truth_contigs)
    audit = {
        "seed": SEED,
        "assembly_branches": len(lineage),
        "assembly_contigs": len(contig_rows),
        "assembly_bases": sum(int(row["AssemblyBases"]) for row in lineage),
        "comparison_design": "two individual 2M-pair assemblies versus one 4M-pair co-assembly per assembler; mix is cascaded clustering",
        "truth": truth_audit,
    }
    (summary / "prepare-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    (work / ".article34-inputs-complete").write_text(json.dumps(audit, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
