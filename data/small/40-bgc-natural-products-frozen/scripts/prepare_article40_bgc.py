#!/usr/bin/env python3
"""Prepare checksum-gated real genomes and a deterministic metagenome subset for Article 40."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


SEED = 20260740
INPUT_HASHES = {
    "salinispora": "f6b2f3c9170344c992c13ef1c915780a1295908a8f2f0dd6661aa12fba205c00",
    "nostoc": "6a30f675c3d94cd6aba7de06a1c27ef2fbd5a38e5424603d790669cf13894bb9",
    "coassembly": "904f92521ff0ce9f12bd52d153bb249ec816fc900051e06b4b12bc5da74a270a",
    "membership": "c39a16f375e25fc04aa38e11e9d64cedc12d4351cf6a84e0672bd34fe25e6ad9",
    "abundance": "8dc7ea7c1f0a6e61625fd1aada492d252b60eed050a91af6560986987e43236b",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--antismash-database", type=Path, required=True)
    parser.add_argument("--antismash-env", type=Path, required=True)
    parser.add_argument("--gecco-env", type=Path, required=True)
    parser.add_argument("--minimum-contig", type=int, default=20_000)
    parser.add_argument("--fragment-size", type=int, default=20_000)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence).upper()))
                name = line[1:].split(None, 1)[0]
                sequence = []
            else:
                sequence.append(line.strip())
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def version(command: list[str], runtime: dict[str, str]) -> str:
    result = subprocess.run(command, env=runtime, text=True, capture_output=True, check=True)
    return " | ".join(x.strip() for x in (result.stdout, result.stderr) if x.strip()).replace("\t", " ")


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    database = args.antismash_database.resolve()
    antismash_env = args.antismash_env.resolve()
    gecco_env = args.gecco_env.resolve()
    sources = {
        "salinispora": root / "data/raw/article32/benchmark_mock/reference/all_genomes_listed/Salinispora_tropica_CNB-440.fna.gz",
        "nostoc": root / "data/raw/article32/benchmark_mock/reference/all_genomes_listed/Nostoc_sp._PCC_7120.fna.gz",
        "coassembly": root / "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz",
        "membership": root / "data/small/34-nonredundant-gene-catalog-frozen/primary-catalog-membership.tsv.gz",
        "abundance": root / "data/small/35-gene-abundance-frozen/gene-abundance-long.tsv.gz",
    }
    for label, source in sources.items():
        if not source.is_file() or sha256(source) != INPUT_HASHES[label]:
            raise ValueError(f"Missing or checksum-mismatched {label}: {source}")

    manifest = root / "data/small/40-bgc-database-manifest.tsv"
    with manifest.open(encoding="utf-8") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    database_audit: list[dict[str, object]] = []
    for row in manifest_rows:
        path = database / row["File"]
        observed_bytes = path.stat().st_size if path.is_file() else 0
        observed_hash = sha256(path) if path.is_file() else "MISSING"
        status = observed_bytes == int(row["Bytes"]) and observed_hash == row["SHA256"]
        database_audit.append({
            "Database": row["Database"], "Release": row["Release"], "File": row["File"],
            "ExpectedBytes": row["Bytes"], "ObservedBytes": observed_bytes,
            "ExpectedSHA256": row["SHA256"], "ObservedSHA256": observed_hash,
            "Status": "PASS" if status else "FAIL",
        })
        if not status:
            raise ValueError(f"antiSMASH database identity failed: {path}")

    if work.exists():
        shutil.rmtree(work)
    for child in ("inputs", "antismash", "gecco", "summary", "logs", "tmp"):
        (work / child).mkdir(parents=True, exist_ok=True)

    salinispora = read_fasta(sources["salinispora"])
    nostoc = read_fasta(sources["nostoc"])
    coassembly = read_fasta(sources["coassembly"])
    if len(salinispora) != 1 or sum(len(s) for _, s in salinispora) != 5_183_331:
        raise ValueError("Unexpected Salinispora tropica assembly")
    if len(nostoc) != 7 or sum(len(s) for _, s in nostoc) != 7_211_789:
        raise ValueError("Unexpected Nostoc PCC 7120 assembly")
    if len(coassembly) != 18_354 or sum(len(s) for _, s in coassembly) != 84_811_518:
        raise ValueError("Unexpected Article 30 co-assembly")

    sal_records = [("salinispora-chromosome", salinispora[0][1])]
    nostoc_records = [(f"nostoc-record-{index:02d}", sequence) for index, (_, sequence) in enumerate(nostoc, 1)]
    co_records: list[tuple[str, str]] = []
    co_map: list[dict[str, object]] = []
    for index, (original, sequence) in enumerate(coassembly, 1):
        normalized = f"megahit-co__c{index:06d}"
        if len(sequence) >= args.minimum_contig:
            co_records.append((normalized, sequence))
            co_map.append({"OriginalContigID": original, "NormalizedContigID": normalized, "LengthBp": len(sequence), "Selected": "yes"})
    if len(co_records) != 675 or sum(len(s) for _, s in co_records) != 39_665_540:
        raise ValueError("Deterministic >=20 kb co-assembly subset changed")

    fragments: list[tuple[str, str]] = []
    fragment_map: list[dict[str, object]] = []
    sequence = salinispora[0][1]
    for index, start in enumerate(range(0, len(sequence), args.fragment_size), 1):
        stop = min(start + args.fragment_size, len(sequence))
        name = f"salinispora-fragment-{index:04d}"
        fragments.append((name, sequence[start:stop]))
        fragment_map.append({"FragmentID": name, "GenomeStart1": start + 1, "GenomeEnd1": stop, "LengthBp": stop - start})

    write_fasta(work / "inputs/salinispora-full.fna", sal_records)
    write_fasta(work / "inputs/salinispora-fragmented.fna", fragments)
    write_fasta(work / "inputs/nostoc.fna", nostoc_records)
    write_fasta(work / "inputs/coassembly-ge20kb.fna", co_records)
    write_tsv(work / "inputs/coassembly-id-map.tsv", co_map)
    write_tsv(work / "inputs/salinispora-fragment-map.tsv", fragment_map)

    antismash_runtime = os.environ.copy()
    antismash_runtime.update({"PATH": f"{antismash_env / 'bin'}:/usr/bin:/bin", "MPLCONFIGDIR": str(work / "tmp/antismash-mpl"), "PYTHONHASHSEED": "0", "PYTHONNOUSERSITE": "1", "PYTHONPATH": ""})
    gecco_runtime = os.environ.copy()
    gecco_runtime.update({"PATH": f"{gecco_env / 'bin'}:/usr/bin:/bin", "PYTHONHASHSEED": "0", "PYTHONNOUSERSITE": "1", "PYTHONPATH": ""})
    antismash = antismash_env / "bin/antismash"
    gecco = gecco_env / "bin/gecco"
    prodigal = antismash_env / "bin/prodigal"
    for binary in (antismash, gecco, prodigal):
        if not binary.is_file():
            raise FileNotFoundError(binary)
    prereq = subprocess.run([str(antismash), "--databases", str(database), "--check-prereqs", "--taxon", "bacteria"], env=antismash_runtime, text=True, capture_output=True, check=True)
    if "All prerequisites satisfied" not in prereq.stdout + prereq.stderr:
        raise ValueError("antiSMASH prerequisite audit did not pass")
    tool_rows = [
        {"Tool": "antiSMASH", "VersionEvidence": version([str(antismash), "--version"], antismash_runtime)},
        {"Tool": "GECCO", "VersionEvidence": version([str(gecco), "--version"], gecco_runtime)},
        {"Tool": "Prodigal", "VersionEvidence": version([str(prodigal), "-v"], antismash_runtime)},
    ]
    write_tsv(work / "database-audit.tsv", database_audit)
    write_tsv(work / "tool-versions.tsv", tool_rows)
    write_tsv(work / "input-audit.tsv", [
        {"Asset": "Salinispora tropica CNB-440", "Records": 1, "Bases": 5_183_331, "SelectedRecords": 1, "SelectedBases": 5_183_331, "CompressedSHA256": INPUT_HASHES["salinispora"]},
        {"Asset": "Salinispora deterministic 20 kb fragments", "Records": len(fragments), "Bases": 5_183_331, "SelectedRecords": len(fragments), "SelectedBases": 5_183_331, "CompressedSHA256": "derived"},
        {"Asset": "Nostoc sp. PCC 7120", "Records": 7, "Bases": 7_211_789, "SelectedRecords": 7, "SelectedBases": 7_211_789, "CompressedSHA256": INPUT_HASHES["nostoc"]},
        {"Asset": "MEGAHIT co-assembly >=1 kb", "Records": 18_354, "Bases": 84_811_518, "SelectedRecords": len(co_records), "SelectedBases": sum(len(s) for _, s in co_records), "CompressedSHA256": INPUT_HASHES["coassembly"]},
    ])
    write_tsv(work / "input-lineage.tsv", [
        {"Asset": "Salinispora tropica CNB-440", "Identifier": "GCA_000016425.1 / CP000667.1", "Role": "complete-genome BGC positive control"},
        {"Asset": "Nostoc sp. PCC 7120", "Identifier": "GCA_000009705.1 / BA000019.2 plus plasmids", "Role": "cyanobacterial BGC positive control"},
        {"Asset": "MEGAHIT co-assembly", "Identifier": "Article 30 checksum-locked co-assembly; deterministic contigs >=20 kb", "Role": "real fragmented-metagenome branch"},
        {"Asset": "Gene catalog membership", "Identifier": "Article 34 checksum-locked primary membership", "Role": "map co-assembly ORFs to representative genes"},
        {"Asset": "Gene abundance", "Identifier": "Article 35 checksum-locked assigned-read ledger", "Role": "read-weighted BGC profile"},
    ])
    contract = {
        "article": 40, "seed": SEED, "antismash": "8.0.4", "gecco": "0.10.3",
        "prodigal": "2.6.3", "pfam": "35.0", "mibig": "4.0", "mite": "1.3",
        "minimum_metagenome_contig_bp": args.minimum_contig, "fragment_size_bp": args.fragment_size,
        "antismash_modules": ["core detection", "cc-mibig", "cb-knownclusters"],
        "gecco_threshold": 0.8, "gecco_minimum_cds": 3,
        "catalog_membership_sha256": INPUT_HASHES["membership"], "abundance_sha256": INPUT_HASHES["abundance"],
        "novelty_claim": "similarity-to-MIBiG only; not compound novelty",
    }
    (work / "run-contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    (work / ".article40-inputs-complete").write_text(json.dumps(contract, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
