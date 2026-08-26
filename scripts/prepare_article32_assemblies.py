#!/usr/bin/env python3
"""Normalize Article 32 assemblies at one explicit length threshold."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", action="append", required=True, help="BRANCH=FASTA")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-length", type=int, default=1000)
    return parser.parse_args()


def parse_fasta(path: Path):
    header = None
    chunks: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header, chunks = line[1:], []
            else:
                if header is None:
                    raise ValueError(f"Sequence before header at {path}:{line_number}")
                if re.search(r"[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]", line):
                    raise ValueError(f"Invalid nucleotide at {path}:{line_number}")
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks).upper()


def nx(lengths: list[int], fraction: float) -> tuple[int, int]:
    target = sum(lengths) * fraction
    cumulative = 0
    for index, length in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += length
        if cumulative >= target:
            return length, index
    return 0, 0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    inputs: list[tuple[str, Path]] = []
    for value in args.assembly:
        if "=" not in value:
            raise SystemExit(f"Assembly must be BRANCH=FASTA: {value}")
        branch, raw_path = value.split("=", 1)
        path = Path(raw_path).resolve()
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing assembly: {path}")
        inputs.append((branch, path))
    if len({branch for branch, _ in inputs}) != len(inputs):
        raise SystemExit("Duplicate branch label")

    rows = []
    for branch, source in inputs:
        destination = args.output_dir / f"{branch}.ge{args.min_length}.fasta"
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        lengths: list[int] = []
        source_records = retained_records = source_bases = retained_bases = 0
        gc = ns = 0
        identifiers: set[str] = set()
        with temporary.open("w", encoding="utf-8", newline="\n") as out:
            for header, sequence in parse_fasta(source):
                source_records += 1
                source_bases += len(sequence)
                if len(sequence) < args.min_length:
                    continue
                identifier = header.split()[0]
                if identifier in identifiers:
                    raise ValueError(f"Duplicate FASTA identifier in {source}: {identifier}")
                identifiers.add(identifier)
                retained_records += 1
                retained_bases += len(sequence)
                lengths.append(len(sequence))
                gc += sequence.count("G") + sequence.count("C")
                ns += sequence.count("N")
                out.write(f">{header}\n")
                for offset in range(0, len(sequence), 80):
                    out.write(sequence[offset : offset + 80] + "\n")
        if not lengths:
            temporary.unlink(missing_ok=True)
            raise ValueError(f"No sequence >= {args.min_length} bp in {source}")
        temporary.replace(destination)
        n50, l50 = nx(lengths, 0.5)
        n90, l90 = nx(lengths, 0.9)
        rows.append(
            {
                "Branch": branch,
                "SourceRecords": source_records,
                "SourceBases": source_bases,
                "MinimumLengthBp": args.min_length,
                "Sequences": retained_records,
                "TotalLengthBp": retained_bases,
                "LargestBp": max(lengths),
                "N50Bp": n50,
                "L50": l50,
                "N90Bp": n90,
                "L90": l90,
                "SequencesGe10kb": sum(length >= 10_000 for length in lengths),
                "SequencesGe100kb": sum(length >= 100_000 for length in lengths),
                "SequencesGe1Mb": sum(length >= 1_000_000 for length in lengths),
                "NBases": ns,
                "GCPctExcludingN": f"{100 * gc / (retained_bases - ns):.6f}" if retained_bases > ns else "NA",
                "SHA256": sha256(destination),
                "NormalizedFile": destination.name,
            }
        )
    fields = list(rows[0])
    with (args.output_dir / "assembly-structure.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
