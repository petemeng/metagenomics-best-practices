#!/usr/bin/env python3
"""Prepare checksum-audited assemblies, truth sets, and diagnostic controls for Article 33."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator, TextIO

import openpyxl


SEED = 20260733
FRAGMENT_BP = 50_000
CHIMERA_BLOCK_BP = 100_000
CHIMERA_CONTIGS = 20

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
        "Branch": "sr-megahit-m1",
        "Display": "MEGAHIT · M1 single",
        "EvaluationSet": "MOCK1",
        "Family": "Short read",
        "Assembler": "MEGAHIT 1.2.9",
        "Strategy": "single · 2M pairs",
        "SourceArticle": "30",
        "SourceRelativePath": "data/small/30-short-read-assembly-frozen/contigs/megahit-single-MOCK1.ge1000.fna.gz",
        "BundleRelativePath": "contigs/megahit-single-MOCK1.ge1000.fna.gz",
    },
    {
        "Branch": "sr-metaspades-m1",
        "Display": "metaSPAdes · M1 single",
        "EvaluationSet": "MOCK1",
        "Family": "Short read",
        "Assembler": "metaSPAdes 4.3.0",
        "Strategy": "single · 2M pairs",
        "SourceArticle": "30",
        "SourceRelativePath": "data/small/30-short-read-assembly-frozen/contigs/metaspades-single-MOCK1.ge1000.fna.gz",
        "BundleRelativePath": "contigs/metaspades-single-MOCK1.ge1000.fna.gz",
    },
    {
        "Branch": "sr-megahit-m2",
        "Display": "MEGAHIT · M2 single",
        "EvaluationSet": "MOCK2",
        "Family": "Short read",
        "Assembler": "MEGAHIT 1.2.9",
        "Strategy": "single · 2M pairs",
        "SourceArticle": "30",
        "SourceRelativePath": "data/small/30-short-read-assembly-frozen/contigs/megahit-single-MOCK2.ge1000.fna.gz",
        "BundleRelativePath": "contigs/megahit-single-MOCK2.ge1000.fna.gz",
    },
    {
        "Branch": "sr-metaspades-m2",
        "Display": "metaSPAdes · M2 single",
        "EvaluationSet": "MOCK2",
        "Family": "Short read",
        "Assembler": "metaSPAdes 4.3.0",
        "Strategy": "single · 2M pairs",
        "SourceArticle": "30",
        "SourceRelativePath": "data/small/30-short-read-assembly-frozen/contigs/metaspades-single-MOCK2.ge1000.fna.gz",
        "BundleRelativePath": "contigs/metaspades-single-MOCK2.ge1000.fna.gz",
    },
    {
        "Branch": "sr-megahit-co",
        "Display": "MEGAHIT · co-assembly",
        "EvaluationSet": "MOCK1+MOCK2",
        "Family": "Short read",
        "Assembler": "MEGAHIT 1.2.9",
        "Strategy": "co · 4M pairs",
        "SourceArticle": "30",
        "SourceRelativePath": "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz",
        "BundleRelativePath": "contigs/megahit-coassembly.ge1000.fna.gz",
    },
    {
        "Branch": "sr-metaspades-co",
        "Display": "metaSPAdes · co-assembly",
        "EvaluationSet": "MOCK1+MOCK2",
        "Family": "Short read",
        "Assembler": "metaSPAdes 4.3.0",
        "Strategy": "co · 4M pairs",
        "SourceArticle": "30",
        "SourceRelativePath": "data/small/30-short-read-assembly-frozen/contigs/metaspades-coassembly.ge1000.fna.gz",
        "BundleRelativePath": "contigs/metaspades-coassembly.ge1000.fna.gz",
    },
    {
        "Branch": "lr-flye-ont",
        "Display": "Flye · ONT R9",
        "EvaluationSet": "MOCK1",
        "Family": "Long read",
        "Assembler": "Flye 2.9.6",
        "Strategy": "ONT · full run",
        "SourceArticle": "31",
        "SourceRelativePath": "data/small/31-long-read-assembly-frozen/assemblies/flye-ont-r9.ge1000.fna.gz",
        "BundleRelativePath": "assemblies/flye-ont-r9.ge1000.fna.gz",
    },
    {
        "Branch": "lr-flye-hifi",
        "Display": "Flye · HiFi",
        "EvaluationSet": "MOCK1",
        "Family": "Long read",
        "Assembler": "Flye 2.9.6",
        "Strategy": "HiFi · full run",
        "SourceArticle": "31",
        "SourceRelativePath": "data/small/31-long-read-assembly-frozen/assemblies/flye-hifi.ge1000.fna.gz",
        "BundleRelativePath": "assemblies/flye-hifi.ge1000.fna.gz",
    },
    {
        "Branch": "lr-hifiasm-hifi",
        "Display": "hifiasm-meta · HiFi",
        "EvaluationSet": "MOCK1",
        "Family": "Long read",
        "Assembler": "hifiasm-meta 0.3.2-r686",
        "Strategy": "HiFi · full run",
        "SourceArticle": "31",
        "SourceRelativePath": "data/small/31-long-read-assembly-frozen/assemblies/hifiasm-meta-hifi.ge1000.fna.gz",
        "BundleRelativePath": "assemblies/hifiasm-meta-hifi.ge1000.fna.gz",
    },
    {
        "Branch": "lr-metamdbg-hifi",
        "Display": "metaMDBG · HiFi",
        "EvaluationSet": "MOCK1",
        "Family": "Long read",
        "Assembler": "metaMDBG 1.0",
        "Strategy": "HiFi · full run",
        "SourceArticle": "31",
        "SourceRelativePath": "data/small/31-long-read-assembly-frozen/assemblies/metamdbg-hifi.ge1000.fna.gz",
        "BundleRelativePath": "assemblies/metamdbg-hifi.ge1000.fna.gz",
    },
    {
        "Branch": "sr-spades-10m",
        "Display": "SPAdes · short-only",
        "EvaluationSet": "MOCK1",
        "Family": "Short read",
        "Assembler": "SPAdes 4.3.0",
        "Strategy": "single · 10M pairs",
        "SourceArticle": "32",
        "SourceRelativePath": "data/small/32-hybrid-assembly-polishing-frozen/assemblies/spades-short-only.ge1000.fna.gz",
        "BundleRelativePath": "assemblies/spades-short-only.ge1000.fna.gz",
    },
    {
        "Branch": "hybrid-illumina-ont",
        "Display": "SPAdes · Illumina+ONT",
        "EvaluationSet": "MOCK1",
        "Family": "Hybrid",
        "Assembler": "SPAdes 4.3.0",
        "Strategy": "10M pairs + ONT",
        "SourceArticle": "32",
        "SourceRelativePath": "data/small/32-hybrid-assembly-polishing-frozen/assemblies/spades-illumina-ont.ge1000.fna.gz",
        "BundleRelativePath": "assemblies/spades-illumina-ont.ge1000.fna.gz",
    },
    {
        "Branch": "hybrid-illumina-hifi",
        "Display": "SPAdes · Illumina+HiFi",
        "EvaluationSet": "MOCK1",
        "Family": "Hybrid",
        "Assembler": "SPAdes 4.3.0",
        "Strategy": "10M pairs + HiFi",
        "SourceArticle": "32",
        "SourceRelativePath": "data/small/32-hybrid-assembly-polishing-frozen/assemblies/spades-illumina-hifi.ge1000.fna.gz",
        "BundleRelativePath": "assemblies/spades-illumina-hifi.ge1000.fna.gz",
    },
    {
        "Branch": "lr-ont-polypolish",
        "Display": "Flye ONT · Polypolish",
        "EvaluationSet": "MOCK1",
        "Family": "Polished long read",
        "Assembler": "Flye 2.9.6 + Polypolish 0.6.1",
        "Strategy": "ONT + 10M pairs",
        "SourceArticle": "32",
        "SourceRelativePath": "data/small/32-hybrid-assembly-polishing-frozen/assemblies/flye-ont-polypolish-default.ge1000.fna.gz",
        "BundleRelativePath": "assemblies/flye-ont-polypolish-default.ge1000.fna.gz",
    },
    {
        "Branch": "lr-ont-careful",
        "Display": "Flye ONT · careful",
        "EvaluationSet": "MOCK1",
        "Family": "Polished long read",
        "Assembler": "Flye 2.9.6 + Polypolish 0.6.1",
        "Strategy": "ONT + 10M pairs",
        "SourceArticle": "32",
        "SourceRelativePath": "data/small/32-hybrid-assembly-polishing-frozen/assemblies/flye-ont-polypolish-careful.ge1000.fna.gz",
        "BundleRelativePath": "assemblies/flye-ont-polypolish-careful.ge1000.fna.gz",
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
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


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
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def write_tsv(path: Path, rows: Iterable[dict], fields: list[str] | None = None) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    fields = fields or list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_checksum_manifest(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        digest, relative = line.split("  ", 1)
        result[relative] = digest
    return result


def canonical_hash(records: Iterable[tuple[str, str]]) -> str:
    digest = hashlib.sha256()
    for name, sequence in records:
        digest.update(name.split(None, 1)[0].encode())
        digest.update(b"\0")
        digest.update(sequence.encode())
        digest.update(b"\n")
    return digest.hexdigest()


def nx(lengths: list[int], fraction: float) -> tuple[int, int]:
    target = sum(lengths) * fraction
    cumulative = 0
    for index, length in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += length
        if cumulative >= target:
            return length, index
    return 0, 0


def metrics(records: list[tuple[str, str]]) -> dict[str, object]:
    identifiers = [name.split(None, 1)[0] for name, _ in records]
    if not records or len(identifiers) != len(set(identifiers)):
        raise ValueError("Empty FASTA or duplicate first-field identifiers")
    lengths = [len(sequence) for _, sequence in records]
    n50, l50 = nx(lengths, 0.5)
    n90, l90 = nx(lengths, 0.9)
    alphabet = "ACGTNRYKMSWBDHV"
    counts: Counter[str] = Counter()
    for _, sequence in records:
        for base in alphabet:
            counts[base] += sequence.count(base)
    invalid = sum(lengths) - sum(counts.values())
    informative = sum(counts[base] for base in "ACGT")
    return {
        "Contigs": len(records),
        "TotalBp": sum(lengths),
        "LargestBp": max(lengths),
        "N50Bp": n50,
        "L50": l50,
        "N90Bp": n90,
        "L90": l90,
        "GCPercent": f"{100 * (counts['G'] + counts['C']) / informative:.9f}" if informative else "NA",
        "NBases": counts["N"],
        "InvalidIUPACBases": invalid,
        "ContigsGe10kb": sum(length >= 10_000 for length in lengths),
        "BasesGe10kb": sum(length for length in lengths if length >= 10_000),
        "ContigsGe100kb": sum(length >= 100_000 for length in lengths),
        "BasesGe100kb": sum(length for length in lengths if length >= 100_000),
        "CanonicalSHA256": canonical_hash(records),
    }


def fragment_records(records: list[tuple[str, str]]) -> list[tuple[str, str]]:
    fragmented: list[tuple[str, str]] = []
    for name, sequence in records:
        identifier = name.split(None, 1)[0]
        chunks = [sequence[start : start + FRAGMENT_BP] for start in range(0, len(sequence), FRAGMENT_BP)]
        if len(chunks) > 1 and len(chunks[-1]) < 1000:
            chunks[-2] += chunks[-1]
            chunks.pop()
        for index, chunk in enumerate(chunks, 1):
            if len(chunk) < 1000:
                raise ValueError(f"Diagnostic fragment below 1 kb: {identifier}")
            fragmented.append((f"frag_{identifier}_{index:04d}", chunk))
    return fragmented


def chimeric_records(records: list[tuple[str, str]]) -> tuple[list[tuple[str, str]], list[dict]]:
    eligible = sorted(
        ((index, name, sequence) for index, (name, sequence) in enumerate(records) if len(sequence) >= 3 * CHIMERA_BLOCK_BP),
        key=lambda item: (-len(item[2]), item[1]),
    )[:CHIMERA_CONTIGS]
    if len(eligible) != CHIMERA_CONTIGS:
        raise ValueError(f"Expected {CHIMERA_CONTIGS} long contigs for the diagnostic control, found {len(eligible)}")
    output = list(records)
    blocks = []
    for index, name, sequence in eligible:
        start = (len(sequence) - CHIMERA_BLOCK_BP) // 2
        blocks.append((index, name, start, sequence[start : start + CHIMERA_BLOCK_BP]))
    audit = []
    for position, (index, name, start, old_block) in enumerate(blocks):
        donor_index, donor_name, _, new_block = blocks[(position - 1) % len(blocks)]
        sequence = output[index][1]
        output[index] = (output[index][0], sequence[:start] + new_block + sequence[start + CHIMERA_BLOCK_BP :])
        audit.append(
            {
                "TargetContig": name.split(None, 1)[0],
                "TargetStart0": start,
                "BlockBp": CHIMERA_BLOCK_BP,
                "DonorContig": donor_name.split(None, 1)[0],
                "OriginalBlockSHA256": hashlib.sha256(old_block.encode()).hexdigest(),
                "ReplacementBlockSHA256": hashlib.sha256(new_block.encode()).hexdigest(),
                "TargetRecordIndex": index,
                "DonorRecordIndex": donor_index,
            }
        )
    return output, audit


def sequence_counter(records: Iterable[tuple[str, str]]) -> Counter[str]:
    """Count IUPAC symbols using C-level str.count, not per-base Python updates."""
    counter: Counter[str] = Counter()
    for _, sequence in records:
        for base in "ACGTNRYKMSWBDHV":
            counter[base] += sequence.count(base)
    return counter


def fasta_inventory(path: Path) -> tuple[int, int, Counter[str]]:
    records = list(fasta_records(path))
    return len(records), sum(len(sequence) for _, sequence in records), Counter(
        hashlib.sha256(sequence.encode()).hexdigest() for _, sequence in records
    )


def build_truth(repo: Path, work: Path) -> tuple[list[dict], dict]:
    s1 = repo / "script_r/Supplementary_Table_S1.xlsx"
    genome_dir = repo / "reference/all_genomes_listed"
    workbook = openpyxl.load_workbook(s1, read_only=True, data_only=True)
    values = list(workbook.active.iter_rows(values_only=True))
    headers = values[0]
    accession_to_row = {str(row[3]).strip(): row for row in values[1:] if row[3]}
    set_accessions = {
        "MOCK1": [line.strip() for line in (repo / "profiling/MOCK_001.list").read_text().splitlines() if line.strip()],
        "MOCK2": [line.strip() for line in (repo / "profiling/MOCK_002.list").read_text().splitlines() if line.strip()],
    }
    if not set(set_accessions["MOCK1"]) < set(set_accessions["MOCK2"]):
        raise ValueError("Expected MOCK1 truth accessions to be a strict subset of MOCK2")
    truth_rows: list[dict] = []
    set_spec = (("MOCK1", "MOCK1", 14), ("MOCK2", "MOCK2", 15), ("MOCK1+MOCK2", "MOCK2", None))
    for evaluation_set, member_set, abundance_column in set_spec:
        refs_dir = work / "truth/references" / evaluation_set
        refs_dir.mkdir(parents=True, exist_ok=True)
        for accession in set_accessions[member_set]:
            row = accession_to_row.get(accession)
            if row is None:
                raise ValueError(f"No Supplementary Table S1 row for {accession}")
            current_label = str(row[2]).strip()
            repository_label = ALIASES.get(current_label, current_label)
            source = genome_dir / f"{repository_label}.fna.gz"
            if not source.is_file():
                raise FileNotFoundError(source)
            link = refs_dir / f"{current_label}.fna.gz"
            if link.is_symlink():
                link.unlink()
            elif link.exists():
                raise ValueError(f"Truth path is not a symlink: {link}")
            os.symlink(source.resolve(), link)
            abundance_m1 = 0.0 if row[14] in (None, "NP") else float(row[14])
            abundance_m2 = 0.0 if row[15] in (None, "NP") else float(row[15])
            abundance = (
                abundance_m1 if abundance_column == 14 else abundance_m2 if abundance_column == 15 else (abundance_m1 + abundance_m2) / 2
            )
            record_count, bases, _ = fasta_inventory(source)
            truth_rows.append(
                {
                    "EvaluationSet": evaluation_set,
                    "Reference": current_label,
                    "RepositoryLabel": repository_label,
                    "GenBankAssembly": accession,
                    "ExpectedAbundancePct": f"{abundance:.12g}",
                    "ExpectedAbundanceMOCK1Pct": f"{abundance_m1:.12g}",
                    "ExpectedAbundanceMOCK2Pct": f"{abundance_m2:.12g}",
                    "ReferenceRecords": record_count,
                    "ReferenceBases": bases,
                    "CompressedBytes": source.stat().st_size,
                    "CompressedSHA256": sha256(source),
                    "TaxonomyRenameAlias": "yes" if current_label != repository_label else "no",
                }
            )
    for set_name, expected in (("MOCK1", 71), ("MOCK2", 87), ("MOCK1+MOCK2", 87)):
        observed = sum(row["EvaluationSet"] == set_name for row in truth_rows)
        if observed != expected:
            raise ValueError(f"Expected {expected} references for {set_name}, found {observed}")
    for set_name, combined_name in (("MOCK1", "MOCK_001.fasta.gz"), ("MOCK2", "MOCK_002.fasta.gz")):
        combined = repo / "reference" / combined_name
        aggregate_records = aggregate_bases = 0
        aggregate_hashes: Counter[str] = Counter()
        for path in sorted((work / "truth/references" / set_name).glob("*.fna.gz")):
            records, bases, hashes = fasta_inventory(path)
            aggregate_records += records
            aggregate_bases += bases
            aggregate_hashes.update(hashes)
        observed = fasta_inventory(combined)
        if (aggregate_records, aggregate_bases, aggregate_hashes) != observed:
            raise ValueError(f"Per-genome references do not equal {combined_name}")
    audit = {
        "seed": SEED,
        "supplement_s1_sha256": sha256(s1),
        "mock1_truth_genomes": 71,
        "mock2_truth_genomes": 87,
        "mock1_is_strict_subset_of_mock2": True,
        "coassembly_truth_policy": "union of MOCK1 and MOCK2; equal to MOCK2 membership",
        "coassembly_abundance_policy": "arithmetic mean of MOCK1 and MOCK2 expected DNA percentages; NP treated as zero",
        "supplement_headers": [str(value) for value in headers[:17]],
    }
    return truth_rows, audit


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    repo = args.benchmark_repo.resolve()
    work = args.work_dir.resolve()
    assemblies_dir = work / "assemblies"
    summary_dir = work / "summary"
    assemblies_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    bundle_manifests = {
        "30": read_checksum_manifest(root / "data/small/30-short-read-assembly-frozen/file-checksums.sha256"),
        "31": read_checksum_manifest(root / "data/small/31-long-read-assembly-frozen/file-checksums.sha256"),
        "32": read_checksum_manifest(root / "data/small/32-hybrid-assembly-polishing-frozen/file-checksums.sha256"),
    }
    lineage_rows = []
    source_audit = []
    flye_hifi_records: list[tuple[str, str]] | None = None
    for spec in BRANCHES:
        source = root / spec["SourceRelativePath"]
        if not source.is_file():
            raise FileNotFoundError(source)
        observed = sha256(source)
        expected = bundle_manifests[spec["SourceArticle"]].get(spec["BundleRelativePath"])
        if expected is None or observed != expected:
            raise ValueError(f"Upstream checksum mismatch for {source}: expected {expected}, observed {observed}")
        target = assemblies_dir / f"{spec['Branch']}.fasta"
        if target.is_file() and target.stat().st_size > 0:
            records = list(fasta_records(target))
            print(f"Reusing prepared {spec['Branch']}", flush=True)
        else:
            records = list(fasta_records(source))
            write_fasta(target, records)
            print(f"Prepared {spec['Branch']}", flush=True)
        if spec["Branch"] == "lr-flye-hifi":
            flye_hifi_records = records
        row = {
            **{key: value for key, value in spec.items() if key != "BundleRelativePath"},
            "EvidenceClass": "biological",
            "PreparedRelativePath": target.relative_to(work).as_posix(),
            "SourceCompressedBytes": source.stat().st_size,
            "SourceCompressedSHA256": observed,
            **metrics(records),
        }
        lineage_rows.append(row)
        source_audit.append(
            {
                "Branch": spec["Branch"],
                "SourceArticle": spec["SourceArticle"],
                "SourceRelativePath": spec["SourceRelativePath"],
                "BundleRelativePath": spec["BundleRelativePath"],
                "ExpectedSHA256": expected,
                "ObservedSHA256": observed,
                "Status": "PASS",
            }
        )

    if flye_hifi_records is None:
        raise ValueError("The Flye HiFi source branch was not prepared")
    base = flye_hifi_records
    fragments = fragment_records(base)
    chimera, block_audit = chimeric_records(base)
    print("Built diagnostic sequences", flush=True)
    if sequence_counter(base) != sequence_counter(fragments):
        raise ValueError("Fragmented diagnostic control does not preserve the assembly base multiset")
    if sequence_counter(base) != sequence_counter(chimera):
        raise ValueError("Chimeric diagnostic control does not preserve the assembly base multiset")
    if sorted(len(sequence) for _, sequence in base) != sorted(len(sequence) for _, sequence in chimera):
        raise ValueError("Chimeric diagnostic control changed the contig length multiset")
    control_specs = (
        ("diagnostic-fragmented-50kb", "Diagnostic · fragmented", fragments, f"deterministic split at {FRAGMENT_BP} bp"),
        ("diagnostic-chimeric-rotation", "Diagnostic · block rotation", chimera, f"rotate {CHIMERA_BLOCK_BP}-bp blocks across {CHIMERA_CONTIGS} contigs"),
    )
    for branch, display, records, strategy in control_specs:
        target = assemblies_dir / f"{branch}.fasta"
        write_fasta(target, records)
        print(f"Prepared {branch}", flush=True)
        lineage_rows.append(
            {
                "Branch": branch,
                "Display": display,
                "EvaluationSet": "MOCK1",
                "Family": "Diagnostic control",
                "Assembler": "derived from Flye HiFi",
                "Strategy": strategy,
                "SourceArticle": "33",
                "SourceRelativePath": "derived from lr-flye-hifi",
                "EvidenceClass": "diagnostic",
                "PreparedRelativePath": target.relative_to(work).as_posix(),
                "SourceCompressedBytes": "NA",
                "SourceCompressedSHA256": "NA",
                **metrics(records),
            }
        )

    truth_rows, truth_audit = build_truth(repo, work)
    print("Verified MOCK1/MOCK2 truth sets", flush=True)
    write_tsv(summary_dir / "input-lineage.tsv", lineage_rows)
    write_tsv(summary_dir / "source-bundle-audit.tsv", source_audit)
    write_tsv(summary_dir / "control-block-audit.tsv", block_audit)
    write_tsv(summary_dir / "truth-manifest.tsv", truth_rows)
    (summary_dir / "truth-audit.json").write_text(json.dumps(truth_audit, indent=2) + "\n", encoding="utf-8")
    run_summary = {
        "seed": SEED,
        "biological_assemblies": len(BRANCHES),
        "diagnostic_controls": 2,
        "evaluation_branches": len(lineage_rows),
        "evaluation_sets": {name: sum(row["EvaluationSet"] == name for row in lineage_rows) for name in ("MOCK1", "MOCK2", "MOCK1+MOCK2")},
        "truth_genomes": {"MOCK1": 71, "MOCK2": 87, "MOCK1+MOCK2": 87},
        "fragment_bp": FRAGMENT_BP,
        "chimera_block_bp": CHIMERA_BLOCK_BP,
        "chimera_contigs": CHIMERA_CONTIGS,
    }
    (summary_dir / "prepare-summary.json").write_text(json.dumps(run_summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(run_summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
