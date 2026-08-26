#!/usr/bin/env python3
"""Prepare checksum-gated real MAGs and shared proteins for Article 60."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import random
import shlex
import subprocess
import time
from pathlib import Path
from typing import Iterable, Iterator


PRIMARY_PANEL = (
    "SGB_002",  # complete nitrifier
    "SGB_006",  # contamination stress case
    "SGB_008",  # incomplete archaeal MAG
    "SGB_010",  # complete methanogen
    "SGB_016",  # complete anaerobic heterotroph
    "SGB_018",  # reduced obligate symbiont
    "SGB_021",  # large, mildly contaminated genome
    "SGB_024",  # complete thermophilic heterotroph
)
TRUNCATION_PARENT = "SGB_002"
TRUNCATION_SEED = 59002
TRUNCATION_FRACTIONS = (0.50, 0.70, 0.90, 1.00)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--prodigal-env", type=Path, required=True)
    parser.add_argument("--gapseq-env", type=Path, required=True)
    parser.add_argument("--carveme-env", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path, algorithm: str = "sha256") -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def digest_gzip_payload(path: Path) -> str:
    hasher = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    header: str | None = None
    sequence: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def record_set_digest(records: Iterable[tuple[str, str]]) -> str:
    hasher = hashlib.sha256()
    for header, sequence in sorted(records):
        hasher.update(header.encode())
        hasher.update(b"\0")
        hasher.update(sequence.encode())
        hasher.update(b"\0")
    return hasher.hexdigest()


def sequence_multiset_digest(records: Iterable[tuple[str, str]]) -> str:
    hasher = hashlib.sha256()
    for sequence in sorted(sequence for _, sequence in records):
        hasher.update(sequence.encode())
        hasher.update(b"\0")
    return hasher.hexdigest()


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> tuple[int, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    residues = 0
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            count += 1
            residues += len(sequence)
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    if count == 0 or residues == 0:
        raise RuntimeError(f"Empty FASTA output: {path}")
    return count, residues


def capture(command: list[str], env: dict[str, str] | None = None) -> str:
    completed = subprocess.run(
        command,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if completed.returncode != 0 or not lines:
        raise RuntimeError(f"Version command failed: {shlex.join(command)}\n{completed.stdout}")
    return lines[0]


def run_logged(
    command: list[str], log_prefix: Path, env: dict[str, str] | None = None
) -> tuple[float, int]:
    log_prefix.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    elapsed = time.perf_counter() - started
    log_prefix.with_suffix(".stdout.log").write_text(completed.stdout, encoding="utf-8")
    log_prefix.with_suffix(".stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed ({completed.returncode}): {shlex.join(command)}\n{completed.stderr}"
        )
    return elapsed, completed.returncode


def prepare_genomes(root: Path, work: Path) -> list[dict[str, object]]:
    taxonomy_path = root / "data/small/46-gtdbtk-taxonomy-frozen/taxonomy-summary.tsv"
    genome_dir = root / "data/small/45-drep-dereplication-frozen/representative-genomes"
    taxonomy = {row["SGB"]: row for row in read_tsv(taxonomy_path)}
    if set(PRIMARY_PANEL) - set(taxonomy):
        raise RuntimeError("The locked taxonomy table is missing one or more panel genomes")

    rows: list[dict[str, object]] = []
    source_records: dict[str, list[tuple[str, str]]] = {}
    for genome in PRIMARY_PANEL:
        row = taxonomy[genome]
        source = genome_dir / f"{row['Representative']}.gz"
        if not source.is_file():
            raise FileNotFoundError(source)
        payload_sha = digest_gzip_payload(source)
        if payload_sha != row["RepresentativeSHA256"]:
            raise RuntimeError(f"Representative payload checksum mismatch for {genome}")
        records = list(fasta_records(source))
        source_records[genome] = records
        output = work / "inputs/genomes" / f"{genome}.fna"
        contigs, bases = write_fasta(output, records)
        rows.append(
            {
                "Genome": genome,
                "AnalysisSet": "Primary real MAG",
                "Representative": row["Representative"],
                "RepresentativeSHA256": payload_sha,
                "SequenceSetSHA256": record_set_digest(records),
                "InputFileSHA256": digest(output),
                "GTDBRelease": row["GTDBRelease"],
                "Species": row["Species"].removeprefix("s__") or "Unresolved",
                "Domain": row["Domain"].removeprefix("d__"),
                "Phylum": row["Phylum"].removeprefix("p__"),
                "CompletenessPct": row["Completeness"],
                "ContaminationPct": row["Contamination"],
                "MIMAGQuality": row["MIMAGQuality"],
                "Contigs": contigs,
                "GenomeBp": bases,
                "RetentionTargetPct": "100",
                "RetentionObservedPct": "100",
                "ParentGenome": genome,
                "DeterministicSeed": "NA",
            }
        )

    parent_records = source_records[TRUNCATION_PARENT]
    parent_bp = sum(len(sequence) for _, sequence in parent_records)
    order = list(range(len(parent_records)))
    random.Random(TRUNCATION_SEED).shuffle(order)
    parent_row = next(row for row in rows if row["Genome"] == TRUNCATION_PARENT)
    for fraction in TRUNCATION_FRACTIONS:
        retained: list[tuple[str, str]] = []
        retained_bp = 0
        for index in order:
            retained.append(parent_records[index])
            retained_bp += len(parent_records[index][1])
            if retained_bp >= parent_bp * fraction:
                break
        genome = f"TRUNC_{int(fraction * 100):03d}"
        output = work / "inputs/genomes" / f"{genome}.fna"
        contigs, bases = write_fasta(output, retained)
        rows.append(
            {
                "Genome": genome,
                "AnalysisSet": "Deterministic truncation sensitivity",
                "Representative": parent_row["Representative"],
                "RepresentativeSHA256": parent_row["RepresentativeSHA256"],
                "SequenceSetSHA256": record_set_digest(retained),
                "InputFileSHA256": digest(output),
                "GTDBRelease": parent_row["GTDBRelease"],
                "Species": parent_row["Species"],
                "Domain": parent_row["Domain"],
                "Phylum": parent_row["Phylum"],
                "CompletenessPct": "NA",
                "ContaminationPct": "NA",
                "MIMAGQuality": "In-silico truncation; not a MIMAG category",
                "Contigs": contigs,
                "GenomeBp": bases,
                "RetentionTargetPct": f"{fraction * 100:.0f}",
                "RetentionObservedPct": f"{bases / parent_bp * 100:.6f}",
                "ParentGenome": TRUNCATION_PARENT,
                "DeterministicSeed": TRUNCATION_SEED,
            }
        )
    write_tsv(work / "input-mag-ledger.tsv", rows)
    write_tsv(
        work / "truncation-ledger.tsv",
        [row for row in rows if row["AnalysisSet"].startswith("Deterministic")],
    )
    return rows


def call_genes(
    rows: list[dict[str, object]], work: Path, prodigal_env: Path
) -> list[dict[str, object]]:
    prodigal = prodigal_env / "bin/prodigal"
    if not prodigal.is_file():
        raise FileNotFoundError(prodigal)
    audits: list[dict[str, object]] = []
    commands: list[dict[str, object]] = []
    for row in rows:
        genome = str(row["Genome"])
        nucleotide = work / "inputs/genomes" / f"{genome}.fna"
        raw_proteins = work / "inputs/proteins-raw" / f"{genome}.faa"
        genes = work / "inputs/genes" / f"{genome}.fna"
        gff = work / "inputs/gff" / f"{genome}.gff"
        raw_proteins.parent.mkdir(parents=True, exist_ok=True)
        genes.parent.mkdir(parents=True, exist_ok=True)
        gff.parent.mkdir(parents=True, exist_ok=True)
        command = [
            str(prodigal), "-q", "-p", "meta", "-i", str(nucleotide),
            "-a", str(raw_proteins), "-d", str(genes), "-f", "gff", "-o", str(gff),
        ]
        elapsed, return_code = run_logged(
            command, work / "logs" / f"prodigal-{genome}"
        )
        raw_records = list(fasta_records(raw_proteins))
        prefixed = [(f"{genome}__{header}", sequence) for header, sequence in raw_records]
        proteins = work / "inputs/proteins" / f"{genome}.faa"
        genes_count, amino_acids = write_fasta(proteins, prefixed)
        identifiers = [header.split()[0] for header, _ in prefixed]
        unique_ids = len(identifiers) == len(set(identifiers))
        prefix_pass = all(identifier.startswith(f"{genome}__") for identifier in identifiers)
        if not unique_ids or not prefix_pass:
            raise RuntimeError(f"Protein identifier audit failed for {genome}")
        audits.append(
            {
                "Genome": genome,
                "Genes": genes_count,
                "AminoAcids": amino_acids,
                "ProteinFileSHA256": digest(proteins),
                "ProteinSequenceMultisetSHA256": sequence_multiset_digest(prefixed),
                "UniqueProteinIDs": str(unique_ids).lower(),
                "GenomePrefixPass": str(prefix_pass).lower(),
                "ProdigalMode": "meta",
            }
        )
        commands.append(
            {
                "Genome": genome,
                "Command": shlex.join(command),
                "ElapsedSeconds": f"{elapsed:.3f}",
                "ReturnCode": return_code,
            }
        )
    write_tsv(work / "protein-id-audit.tsv", audits)
    write_tsv(work / "prodigal-command-log.tsv", commands)
    parent_hash = next(
        row["ProteinSequenceMultisetSHA256"] for row in audits if row["Genome"] == "SGB_002"
    )
    control_hash = next(
        row["ProteinSequenceMultisetSHA256"] for row in audits if row["Genome"] == "TRUNC_100"
    )
    if parent_hash != control_hash:
        raise RuntimeError("TRUNC_100 protein sequences do not reproduce SGB_002")
    return audits


def tool_smoke(args: argparse.Namespace, work: Path) -> list[dict[str, object]]:
    gapseq_env = os.environ.copy()
    gapseq_env["PATH"] = f"{args.gapseq_env / 'bin'}:/usr/bin:/bin"
    gapseq_env["HOME"] = str(work / "home/gapseq")
    carveme_env = os.environ.copy()
    carveme_env["PATH"] = f"{args.carveme_env / 'bin'}:/usr/bin:/bin"
    rows = [
        {
            "Tool": "gapseq",
            "VersionEvidence": capture([str(args.gapseq_env / "bin/gapseq"), "-v"], gapseq_env),
        },
        {
            "Tool": "gapseq DIAMOND",
            "VersionEvidence": capture([str(args.gapseq_env / "bin/diamond"), "version"], gapseq_env),
        },
        {
            "Tool": "GLPK",
            "VersionEvidence": capture([str(args.gapseq_env / "bin/glpsol"), "--version"], gapseq_env),
        },
        {
            "Tool": "CarveMe",
            "VersionEvidence": capture(
                [str(args.carveme_env / "bin/python"), "-c", "import carveme; print(carveme.__version__)"],
                carveme_env,
            ),
        },
        {
            "Tool": "CarveMe DIAMOND",
            "VersionEvidence": capture([str(args.carveme_env / "bin/diamond"), "version"], carveme_env),
        },
        {
            "Tool": "SCIP",
            "VersionEvidence": capture([str(args.carveme_env / "bin/scip"), "--version"], carveme_env),
        },
        {
            "Tool": "Prodigal",
            "VersionEvidence": capture([str(args.prodigal_env / "bin/prodigal"), "-v"]),
        },
    ]
    write_tsv(work / "tool-smoke.tsv", rows)
    return rows


def main() -> None:
    args = parse_args()
    for attribute in (
        "project_root", "work_dir", "prodigal_env", "gapseq_env", "carveme_env"
    ):
        setattr(args, attribute, getattr(args, attribute).resolve())
    args.work_dir.mkdir(parents=True, exist_ok=True)
    rows = prepare_genomes(args.project_root, args.work_dir)
    proteins = call_genes(rows, args.work_dir, args.prodigal_env)
    versions = tool_smoke(args, args.work_dir)
    contract = {
        "article": 60,
        "primary_real_mags": len(PRIMARY_PANEL),
        "sensitivity_genomes": len(TRUNCATION_FRACTIONS),
        "truncation_parent": TRUNCATION_PARENT,
        "truncation_seed": TRUNCATION_SEED,
        "protein_caller": "Prodigal 2.6.3",
        "protein_caller_mode": "meta",
        "shared_protein_inputs_across_tools": True,
        "protein_files": len(proteins),
        "tool_smoke_rows": len(versions),
        "gtdb_release": "R232",
    }
    (args.work_dir / "preparation-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (args.work_dir / ".article60-inputs-complete").write_text(
        "checksum-gated Article 60 inputs prepared\n", encoding="utf-8"
    )
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

