#!/usr/bin/env python3
"""Prepare checksum-gated public inputs for Article 57 plasmid/MGE analysis."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterator, TextIO


SEED = 20260757
BENCHMARK_COMMIT = "a429a3724d4593f35b8d7323b20252a6be90e1cd"
REFERENCE_SHA256 = "5617cc377fc503141d7a27d7c52ce874e3393e3939a9fdcbbd43fe0268c6092c"
COASSEMBLY_SHA256 = "904f92521ff0ce9f12bd52d153bb249ec816fc900051e06b4b12bc5da74a270a"
RGI_COASSEMBLY_SHA256 = "ce859071af7596e9c80ec28b37601edcb448441a1549ab41fe0949f99db76fd2"
RGI_STAPH_SHA256 = "b67d3033a6fdd3789978dba97be8ea8a020a5f30ad18a8f8aa6642ebdf7c6175"
STAPH_ASSEMBLY = "GCA_000013465.1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--benchmark-repo", type=Path, required=True)
    parser.add_argument("--coassembly", type=Path, required=True)
    parser.add_argument("--rgi-coassembly", type=Path, required=True)
    parser.add_argument("--rgi-staphylococcus", type=Path, required=True)
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
    return path.open("r", encoding="utf-8")


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
                raise ValueError(f"Sequence before header in {path}")
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


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


def dump_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def safe_link(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_symlink() or target.exists():
        target.unlink()
    target.symlink_to(source.resolve())


def compact_rgi(source: Path, target: Path, source_label: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tier = row["Cut_Off"].strip()
            if tier not in {"Perfect", "Strict"}:
                continue
            rows.append(
                {
                    "Source": source_label,
                    "ORF_ID": row["ORF_ID"],
                    "Contig": row["Contig"],
                    "Start": int(row["Start"]),
                    "Stop": int(row["Stop"]),
                    "Orientation": row["Orientation"],
                    "EvidenceTier": tier,
                    "BestHitARO": row["Best_Hit_ARO"],
                    "ARO": row["ARO"],
                    "ModelType": row["Model_type"],
                    "IdentityPercent": row["Best_Identities"],
                    "ReferenceLengthPercent": row["Percentage Length of Reference Sequence"],
                    "DrugClasses": row["Drug Class"],
                    "ResistanceMechanisms": row["Resistance Mechanism"],
                    "AMRGeneFamilies": row["AMR Gene Family"],
                    "Nudged": row["Nudged"],
                }
            )
    write_tsv(target, rows)
    return rows


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    repo = args.benchmark_repo.resolve()
    coassembly = args.coassembly.resolve()
    rgi_coassembly = args.rgi_coassembly.resolve()
    rgi_staph = args.rgi_staphylococcus.resolve()
    inputs = work / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)

    reference = repo / "reference/MOCK_002.fasta.gz"
    required = (reference, coassembly, rgi_coassembly, rgi_staph)
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    observed_commit = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()
    if observed_commit != BENCHMARK_COMMIT:
        raise RuntimeError(f"Unexpected benchmark commit: {observed_commit}")

    expected = {
        reference: REFERENCE_SHA256,
        coassembly: COASSEMBLY_SHA256,
        rgi_coassembly: RGI_COASSEMBLY_SHA256,
        rgi_staph: RGI_STAPH_SHA256,
    }
    audit_rows: list[dict[str, object]] = []
    for path, digest in expected.items():
        observed = sha256(path)
        passed = observed == digest
        audit_rows.append(
            {
                "Asset": path.name,
                "SourcePath": str(path.relative_to(root)) if path.is_relative_to(root) else str(path),
                "Bytes": path.stat().st_size,
                "ExpectedSHA256": digest,
                "ObservedSHA256": observed,
                "ChecksumPass": str(passed).lower(),
            }
        )
        if not passed:
            raise RuntimeError(f"Checksum gate failed: {path}")

    safe_link(reference, inputs / "mock2-all-references.fna.gz")
    safe_link(coassembly, inputs / "coassembly.fna.gz")

    reference_rows: list[dict[str, object]] = []
    staph_path = inputs / "staphylococcus-USA300.fna"
    reference_bases = 0
    with staph_path.open("w", encoding="utf-8", newline="\n") as staph_handle:
        for header, sequence in fasta_records(reference):
            seq_name = header.split()[0]
            original_lower = header.lower()
            is_plasmid = "plasmid" in original_lower
            if is_plasmid:
                sequence_class = "GenBank-labelled plasmid"
            elif any(
                marker in original_lower
                for marker in ("complete genome", "chromosome", "complete sequence")
            ):
                sequence_class = "Complete cellular replicon"
            else:
                sequence_class = "Draft scaffold/other sequence"
            reference_rows.append(
                {
                    "SeqName": seq_name,
                    "Assembly": seq_name.split("|")[0],
                    "ReferenceLabel": "Plasmid" if is_plasmid else "Other replicon",
                    "ReferenceSequenceClass": sequence_class,
                    "OriginalHeader": header,
                    "LengthBp": len(sequence),
                    "SequenceSHA256": hashlib.sha256(sequence.encode()).hexdigest(),
                }
            )
            reference_bases += len(sequence)
            if seq_name.startswith(STAPH_ASSEMBLY + "|"):
                staph_handle.write(f">{header}\n")
                for start in range(0, len(sequence), 80):
                    staph_handle.write(sequence[start : start + 80] + "\n")

    if len(reference_rows) != 399:
        raise RuntimeError(f"Expected 399 MOCK2 reference replicons, observed {len(reference_rows)}")
    if sum(row["ReferenceLabel"] == "Plasmid" for row in reference_rows) != 43:
        raise RuntimeError("Expected 43 GenBank-labelled plasmid replicons")
    if sum(row["Assembly"] == STAPH_ASSEMBLY for row in reference_rows) != 4:
        raise RuntimeError("Expected four USA300 replicons")
    write_tsv(work / "reference-replicon-labels.tsv", reference_rows)

    plasmid_rows = sorted(
        (row for row in reference_rows if row["ReferenceLabel"] == "Plasmid"),
        key=lambda row: (int(row["LengthBp"]), str(row["SeqName"])),
    )
    unmatched = {
        str(row["SeqName"]): row
        for row in reference_rows
        if row["ReferenceLabel"] == "Other replicon"
    }
    matched_names: set[str] = {str(row["SeqName"]) for row in plasmid_rows}
    match_rows: list[dict[str, object]] = []
    for plasmid in plasmid_rows:
        p_length = int(plasmid["LengthBp"])
        match_name, match = min(
            unmatched.items(),
            key=lambda item: (
                abs(math.log(int(item[1]["LengthBp"]) / p_length)),
                abs(int(item[1]["LengthBp"]) - p_length),
                item[0],
            ),
        )
        unmatched.pop(match_name)
        matched_names.add(match_name)
        match_rows.extend(
            [
                {
                    **plasmid,
                    "PairID": f"pair-{len(match_rows) // 2 + 1:02d}",
                    "PairRole": "Plasmid",
                    "MatchedTo": match_name,
                },
                {
                    **match,
                    "PairID": f"pair-{len(match_rows) // 2 + 1:02d}",
                    "PairRole": "Length-matched other replicon",
                    "MatchedTo": plasmid["SeqName"],
                },
            ]
        )
    write_tsv(work / "reference-benchmark-labels.tsv", match_rows)
    matched_fasta = inputs / "matched-reference-replicons.fna.gz"
    with matched_fasta.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gz_handle:
            with io.TextIOWrapper(gz_handle, encoding="utf-8", newline="\n") as out_handle:
                written = 0
                for header, sequence in fasta_records(reference):
                    seq_name = header.split()[0]
                    if seq_name not in matched_names:
                        continue
                    out_handle.write(f">{header}\n")
                    for start in range(0, len(sequence), 80):
                        out_handle.write(sequence[start : start + 80] + "\n")
                    written += 1
    if written != 86:
        raise RuntimeError(f"Expected 86 matched benchmark replicons, observed {written}")

    coassembly_records = 0
    coassembly_bases = 0
    for _header, sequence in fasta_records(coassembly):
        coassembly_records += 1
        coassembly_bases += len(sequence)
        if len(sequence) < 1_000:
            raise RuntimeError("Co-assembly contains a contig below the locked 1 kb gate")
    if coassembly_records != 18_354 or coassembly_bases != 84_811_518:
        raise RuntimeError(
            f"Unexpected co-assembly dimensions: {coassembly_records}, {coassembly_bases}"
        )

    co_rows = compact_rgi(
        rgi_coassembly, work / "rgi-primary-coassembly.tsv", "MEGAHIT co-assembly"
    )
    staph_rows = compact_rgi(
        rgi_staph, work / "rgi-primary-staphylococcus.tsv", "USA300 reference"
    )
    if len(co_rows) != 34:
        raise RuntimeError(f"Expected 34 primary co-assembly RGI calls, observed {len(co_rows)}")
    if len(staph_rows) != 21:
        raise RuntimeError(f"Expected 21 primary USA300 RGI calls, observed {len(staph_rows)}")

    write_tsv(work / "asset-check-audit.tsv", audit_rows)
    write_tsv(
        work / "input-lineage.tsv",
        [
            {
                "Input": "mock2-all-references.fna.gz",
                "PublicIdentity": "PRJEB52977 MOCK2 exact references",
                "Origin": "Meslier benchmark_mock fixed commit",
                "Records": len(reference_rows),
                "Bases": reference_bases,
                "Role": "full identity inventory and alignment reference",
            },
            {
                "Input": "matched-reference-replicons.fna.gz",
                "PublicIdentity": "43 plasmids plus 43 length-matched other replicons",
                "Origin": "deterministic no-replacement matching within MOCK2 references",
                "Records": len(match_rows),
                "Bases": sum(int(row["LengthBp"]) for row in match_rows),
                "Role": "balanced geNomad reference benchmark",
            },
            {
                "Input": "coassembly.fna.gz",
                "PublicIdentity": "PRJEB52977 MOCK1+MOCK2 4M-pair subset",
                "Origin": "MEGAHIT 1.2.9 meta-sensitive; contigs >=1 kb",
                "Records": coassembly_records,
                "Bases": coassembly_bases,
                "Role": "metagenomic plasmid discovery and ARG context",
            },
            {
                "Input": "staphylococcus-USA300.fna",
                "PublicIdentity": "GCA_000013465.1",
                "Origin": "exact MOCK2 reference set",
                "Records": 4,
                "Bases": sum(
                    int(row["LengthBp"])
                    for row in reference_rows
                    if row["Assembly"] == STAPH_ASSEMBLY
                ),
                "Role": "real plasmid-borne ARG positive control",
            },
        ],
    )
    dump_json(
        work / "run-contract.json",
        {
            "article": 57,
            "seed": SEED,
            "benchmark_commit": BENCHMARK_COMMIT,
            "reference_replicons": len(reference_rows),
            "reference_plasmid_labels": 43,
            "reference_benchmark_replicons": len(match_rows),
            "reference_benchmark_matching": "minimum absolute log-length ratio without replacement; deterministic tie-break by absolute bp difference and SeqName",
            "coassembly_contigs": coassembly_records,
            "coassembly_bases": coassembly_bases,
            "rgi_primary_coassembly_calls": len(co_rows),
            "rgi_primary_staphylococcus_calls": len(staph_rows),
            "genomad_min_score": 0.7,
            "genomad_score_calibration": False,
            "random_output_requested": False,
        },
    )
    (work / ".article57-inputs-complete").write_text("verified\n", encoding="utf-8")
    print(
        "Article 57 inputs verified: "
        f"{len(reference_rows)} references, {coassembly_records} co-assembly contigs, "
        f"{len(co_rows)} co-assembly primary ARG calls"
    )


if __name__ == "__main__":
    main()
