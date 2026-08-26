#!/usr/bin/env python3
"""Normalize Article 31 assembler outputs and build circular-junction targets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import re
from pathlib import Path
from typing import Iterable, Iterator, TextIO


BRANCHES = (
    {
        "branch": "flye-ont-r9",
        "assembler": "Flye",
        "platform": "ONT R9",
        "relative_fasta": "assemblies/flye-ont-r9/assembly.fasta",
        "circular_source": "flye",
    },
    {
        "branch": "flye-hifi",
        "assembler": "Flye",
        "platform": "PacBio HiFi",
        "relative_fasta": "assemblies/flye-hifi/assembly.fasta",
        "circular_source": "flye",
    },
    {
        "branch": "hifiasm-meta-hifi",
        "assembler": "hifiasm-meta",
        "platform": "PacBio HiFi",
        "relative_fasta": "normalized/hifiasm-meta-hifi.primary.fasta",
        "circular_source": "header",
    },
    {
        "branch": "metamdbg-hifi",
        "assembler": "metaMDBG",
        "platform": "PacBio HiFi",
        "relative_fasta": "assemblies/metamdbg-hifi/contigs.fasta.gz",
        "circular_source": "header",
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


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
            else:
                if name is None:
                    raise ValueError(f"FASTA sequence before header: {path}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def write_fasta_record(handle: TextIO, name: str, sequence: str) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def gfa_to_fasta(gfa: Path, fasta: Path) -> None:
    if not gfa.is_file():
        raise FileNotFoundError(gfa)
    fasta.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gfa.open(encoding="utf-8") as source, fasta.open(
        "w", encoding="utf-8", newline="\n"
    ) as target:
        for raw in source:
            if not raw.startswith("S\t"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[2] == "*":
                raise ValueError(f"Sequence-free segment in primary GFA: {gfa}")
            circular = "yes" if fields[1].endswith("c") else "no"
            write_fasta_record(
                target, f"{fields[1]} circular={circular} source=primary_gfa", fields[2]
            )
            count += 1
    if count == 0:
        raise ValueError(f"No sequence segments found in {gfa}")


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def nx(lengths: list[int], fraction: float) -> tuple[int, int]:
    if not lengths:
        return 0, 0
    target = sum(lengths) * fraction
    cumulative = 0
    for rank, length in enumerate(sorted(lengths, reverse=True), start=1):
        cumulative += length
        if cumulative >= target:
            return length, rank
    raise AssertionError("Nx invariant failed")


def flye_circular_map(path: Path) -> dict[str, bool]:
    mapping: dict[str, bool] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 4:
                raise ValueError(f"Malformed Flye assembly_info row: {raw.rstrip()}")
            mapping[fields[0]] = fields[3].upper() == "Y"
    return mapping


def header_circular(header: str) -> bool:
    token = header.split(None, 1)[0]
    explicit = re.search(r"(?:^|\s)circular=(yes|no)(?:\s|$)", header, re.I)
    if explicit:
        return explicit.group(1).lower() == "yes"
    return token.endswith("c")


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    normalized = work / "normalized"
    normalized.mkdir(parents=True, exist_ok=True)

    hifiasm_gfa = work / "assemblies/hifiasm-meta-hifi/asm.p_ctg.gfa"
    hifiasm_fasta = normalized / "hifiasm-meta-hifi.primary.fasta"
    gfa_to_fasta(hifiasm_gfa, hifiasm_fasta)

    inventory_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    path_rows: list[dict[str, object]] = []

    for meta in BRANCHES:
        branch = str(meta["branch"])
        assembler = str(meta["assembler"])
        platform = str(meta["platform"])
        fasta = work / str(meta["relative_fasta"])
        if not fasta.is_file() or fasta.stat().st_size == 0:
            raise FileNotFoundError(fasta)
        circular_by_id: dict[str, bool] = {}
        if meta["circular_source"] == "flye":
            circular_by_id = flye_circular_map(fasta.parent / "assembly_info.txt")

        records = list(fasta_records(fasta))
        if not records:
            raise ValueError(f"No assembly records: {fasta}")
        readback_fasta = normalized / f"{branch}.ge1000.fasta"
        with readback_fasta.open("w", encoding="utf-8", newline="\n") as handle:
            for header, sequence in records:
                if len(sequence) >= 1_000:
                    write_fasta_record(handle, header, sequence)
        if readback_fasta.stat().st_size == 0:
            raise ValueError(f"No >=1 kb contigs for read-back: {branch}")
        seen: set[str] = set()
        branch_candidates: list[tuple[str, str, str]] = []
        for index, (header, sequence) in enumerate(records, start=1):
            contig_id = header.split(None, 1)[0]
            if contig_id in seen:
                raise ValueError(f"Duplicate contig ID in {fasta}: {contig_id}")
            seen.add(contig_id)
            if not sequence or set(sequence) - set("ACGTNRYKMSWBDHV"):
                raise ValueError(f"Unexpected assembly alphabet: {branch}/{contig_id}")
            if meta["circular_source"] == "flye":
                if contig_id not in circular_by_id:
                    raise ValueError(f"Flye metadata missing contig: {branch}/{contig_id}")
                circular = circular_by_id[contig_id]
            else:
                circular = header_circular(header)
            acgt = sum(sequence.count(base) for base in "ACGT")
            gc = sequence.count("G") + sequence.count("C")
            inventory_rows.append(
                {
                    "Branch": branch,
                    "Assembler": assembler,
                    "Platform": platform,
                    "ContigIndex": index,
                    "ContigID": contig_id,
                    "LengthBp": len(sequence),
                    "GCPct": f"{100 * gc / acgt:.6f}" if acgt else "NA",
                    "CircularCandidate": "TRUE" if circular else "FALSE",
                    "SequenceSHA256": hashlib.sha256(sequence.encode()).hexdigest(),
                    "OriginalHeader": header,
                }
            )
            if circular:
                branch_candidates.append((contig_id, header, sequence))

        for threshold in (0, 1_000, 10_000, 100_000, 1_000_000):
            selected = [
                (header, sequence)
                for header, sequence in records
                if len(sequence) >= threshold
            ]
            lengths = [len(sequence) for _, sequence in selected]
            n50, l50 = nx(lengths, 0.5)
            n90, l90 = nx(lengths, 0.9)
            selected_ids = {header.split(None, 1)[0] for header, _ in selected}
            metric_rows.append(
                {
                    "Branch": branch,
                    "Assembler": assembler,
                    "Platform": platform,
                    "ThresholdBp": threshold,
                    "ContigCount": len(lengths),
                    "TotalBp": sum(lengths),
                    "N50Bp": n50,
                    "L50": l50,
                    "N90Bp": n90,
                    "L90": l90,
                    "LargestBp": max(lengths, default=0),
                    "CircularCandidates": sum(
                        contig_id in selected_ids
                        for contig_id, _, _ in branch_candidates
                    ),
                }
            )

        junction_path = normalized / f"{branch}.junctions.fasta"
        with junction_path.open("w", encoding="utf-8", newline="\n") as handle:
            for candidate_index, (contig_id, header, sequence) in enumerate(
                branch_candidates, start=1
            ):
                junction_id = f"{branch}__junction_{candidate_index:06d}"
                eligible = len(sequence) >= 10_000
                flank = min(5_000, len(sequence) // 2) if eligible else 0
                if eligible:
                    junction = sequence[-flank:] + sequence[:flank]
                    write_fasta_record(handle, junction_id, junction)
                candidate_rows.append(
                    {
                        "Branch": branch,
                        "Assembler": assembler,
                        "Platform": platform,
                        "JunctionID": junction_id,
                        "ContigID": contig_id,
                        "LengthBp": len(sequence),
                        "JunctionEligible": "TRUE" if eligible else "FALSE",
                        "JunctionFlankBp": flank,
                        "SequenceSHA256": hashlib.sha256(sequence.encode()).hexdigest(),
                        "OriginalHeader": header,
                    }
                )

        path_rows.append(
            {
                "Branch": branch,
                "Assembler": assembler,
                "Platform": platform,
                "AssemblyFASTA": str(fasta),
                "ReadbackFASTA": str(readback_fasta),
                "JunctionFASTA": str(junction_path),
                "ContigCount": len(records),
                "CircularCandidates": len(branch_candidates),
            }
        )

    write_tsv(
        normalized / "contig-inventory.tsv",
        inventory_rows,
        [
            "Branch",
            "Assembler",
            "Platform",
            "ContigIndex",
            "ContigID",
            "LengthBp",
            "GCPct",
            "CircularCandidate",
            "SequenceSHA256",
            "OriginalHeader",
        ],
    )
    write_tsv(
        normalized / "assembly-metrics.tsv",
        metric_rows,
        [
            "Branch",
            "Assembler",
            "Platform",
            "ThresholdBp",
            "ContigCount",
            "TotalBp",
            "N50Bp",
            "L50",
            "N90Bp",
            "L90",
            "LargestBp",
            "CircularCandidates",
        ],
    )
    write_tsv(
        normalized / "circular-candidates.tsv",
        candidate_rows,
        [
            "Branch",
            "Assembler",
            "Platform",
            "JunctionID",
            "ContigID",
            "LengthBp",
            "JunctionEligible",
            "JunctionFlankBp",
            "SequenceSHA256",
            "OriginalHeader",
        ],
    )
    write_tsv(
        normalized / "assembly-paths.tsv",
        path_rows,
        [
            "Branch",
            "Assembler",
            "Platform",
            "AssemblyFASTA",
            "ReadbackFASTA",
            "JunctionFASTA",
            "ContigCount",
            "CircularCandidates",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
