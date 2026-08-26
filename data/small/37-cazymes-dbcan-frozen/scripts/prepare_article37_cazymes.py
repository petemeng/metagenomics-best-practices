#!/usr/bin/env python3
"""Identity-gate the real catalog, abundance ledger, B. theta control, and dbCAN database."""

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


SEED = 20260737
CATALOG_GENES = 93_782
FAA_SHA256 = "3db88ff78a548dddfc48caa8a17f04bcbb58dcafe2345d9dd31bc4e12f2a3569"
METADATA_SHA256 = "677479da11ef41a6f27f11798c38dfd9b5830b5c564d9a8f436951413acd7c09"
ABUNDANCE_SHA256 = "8dc7ea7c1f0a6e61625fd1aada492d252b60eed050a91af6560986987e43236b"
BTHETA_GZ_SHA256 = "d4e886c57d6148df6d9b3410f15bbbe94979b5adafba7bec4d787e235cd51068"
BTHETA_CANONICAL_SHA256 = "82c3491d50fdebcb40c761e01d49987601df244ee625343a91e4279ef5ca2b9e"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--database-dir", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--database-manifest", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_fasta_sha256(path: Path) -> tuple[str, int, int]:
    digest = hashlib.sha256()
    records = bases = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            else:
                seq = line.strip().upper()
                if seq:
                    bases += len(seq)
                    digest.update(seq.encode())
    return digest.hexdigest(), records, bases


def read_tsv(path: Path, compressed: bool = False) -> list[dict[str, str]]:
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def version(command: list[str], env: dict[str, str]) -> str:
    result = subprocess.run(command, env=env, text=True, capture_output=True, check=True)
    return " | ".join(x.strip() for x in (result.stdout, result.stderr) if x.strip()).replace("\t", " ")


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    db = args.database_dir.resolve()
    env_prefix = args.env_prefix.resolve()
    manifest = (args.database_manifest or root / "data/small/37-dbcan-database-manifest.tsv").resolve()

    faa_gz = root / "data/small/34-nonredundant-gene-catalog-frozen/catalog/megahit-mix-primary.faa.gz"
    metadata_gz = root / "data/small/34-nonredundant-gene-catalog-frozen/primary-catalog-representatives.tsv.gz"
    abundance_gz = root / "data/small/35-gene-abundance-frozen/gene-abundance-long.tsv.gz"
    btheta_gz = root / "data/raw/article32/benchmark_mock/reference/all_genomes_listed/Bacteroides_thetaiotaomicron_VPI-5482.fna.gz"
    expected = {
        faa_gz: FAA_SHA256,
        metadata_gz: METADATA_SHA256,
        abundance_gz: ABUNDANCE_SHA256,
        btheta_gz: BTHETA_GZ_SHA256,
    }
    for path, digest in expected.items():
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Missing or checksum-mismatched input: {path}")
    canonical, btheta_records, btheta_bases = canonical_fasta_sha256(btheta_gz)
    if canonical != BTHETA_CANONICAL_SHA256 or btheta_records != 2 or btheta_bases != 6_293_399:
        raise ValueError("B. thetaiotaomicron canonical sequence identity failed")

    database_rows = read_tsv(manifest)
    database_audit: list[dict[str, object]] = []
    for row in database_rows:
        path = db / row["File"]
        observed_bytes = path.stat().st_size if path.is_file() else 0
        observed_hash = sha256(path) if path.is_file() else "MISSING"
        passed = observed_bytes == int(row["Bytes"]) and observed_hash == row["SHA256"]
        database_audit.append({
            "File": row["File"], "Release": row["Release"],
            "ExpectedBytes": row["Bytes"], "ObservedBytes": observed_bytes,
            "ExpectedSHA256": row["SHA256"], "ObservedSHA256": observed_hash,
            "Status": "PASS" if passed else "FAIL",
        })
        if not passed:
            raise ValueError(f"dbCAN database identity failed: {row['File']}")

    if work.exists():
        shutil.rmtree(work)
    for child in ("inputs", "catalog", "btheta", "summary", "logs", "tmp"):
        (work / child).mkdir(parents=True, exist_ok=True)
    catalog_faa = work / "inputs/catalog.faa"
    with gzip.open(faa_gz, "rt", encoding="utf-8") as source, catalog_faa.open("w", encoding="utf-8") as target:
        shutil.copyfileobj(source, target)
    with gzip.open(btheta_gz, "rt", encoding="utf-8") as source, (work / "inputs/btheta.fna").open("w", encoding="utf-8") as target:
        shutil.copyfileobj(source, target)

    fasta_ids: list[str] = []
    residues = 0
    with catalog_faa.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                fasta_ids.append(line[1:].split(None, 1)[0])
            else:
                residues += len(line.strip().rstrip("*"))
    if len(fasta_ids) != CATALOG_GENES or len(set(fasta_ids)) != CATALOG_GENES:
        raise ValueError("Catalog protein count/uniqueness failed")
    metadata = read_tsv(metadata_gz, compressed=True)
    if {r["RepresentativeID"] for r in metadata} != set(fasta_ids):
        raise ValueError("Catalog metadata IDs do not match FASTA")
    abundance = read_tsv(abundance_gz, compressed=True)
    sample_counts: dict[str, int] = {}
    sample_ids: dict[str, set[str]] = {}
    for row in abundance:
        sample_counts[row["Sample"]] = sample_counts.get(row["Sample"], 0) + int(row["RawCount"])
        sample_ids.setdefault(row["Sample"], set()).add(row["GeneID"])
    if sample_counts != {"MOCK1": 2_784_234, "MOCK2": 2_777_443}:
        raise ValueError(f"Unexpected assigned-read ledger: {sample_counts}")
    if any(ids != set(fasta_ids) for ids in sample_ids.values()):
        raise ValueError("Abundance IDs do not match catalog")

    runtime = os.environ.copy()
    runtime.update({"PATH": f"{env_prefix / 'bin'}:/usr/bin:/bin", "PYTHONHASHSEED": "0", "PYTHONNOUSERSITE": "1", "PYTHONPATH": ""})
    python = env_prefix / "bin/python"
    run_dbcan = env_prefix / "bin/run_dbcan"
    diamond = env_prefix / "bin/diamond"
    for path in (python, run_dbcan, diamond):
        if not path.is_file():
            raise FileNotFoundError(path)
    tool_rows = [
        {"Tool": "Python", "VersionEvidence": version([str(python), "--version"], runtime)},
        {"Tool": "dbCAN", "VersionEvidence": version([str(python), "-c", "import importlib.metadata as m; print(m.version('dbcan'))"], runtime)},
        {"Tool": "DIAMOND", "VersionEvidence": version([str(diamond), "version"], runtime)},
        {"Tool": "pyHMMER", "VersionEvidence": version([str(python), "-c", "import pyhmmer; print(pyhmmer.__version__)"], runtime)},
        {"Tool": "Pyrodigal", "VersionEvidence": version([str(python), "-c", "import pyrodigal; print(pyrodigal.__version__)"], runtime)},
    ]
    write_tsv(work / "database-audit.tsv", database_audit)
    write_tsv(work / "tool-versions.tsv", tool_rows)
    write_tsv(work / "input-audit.tsv", [{
        "CatalogProteins": len(fasta_ids), "UniqueCatalogIDs": len(set(fasta_ids)),
        "CatalogResidues": residues, "MOCK1AssignedReads": sample_counts["MOCK1"],
        "MOCK2AssignedReads": sample_counts["MOCK2"], "BthetaRecords": btheta_records,
        "BthetaBases": btheta_bases,
    }])
    write_tsv(work / "input-lineage.tsv", [
        {"Asset": "Catalog proteins", "Source": str(faa_gz.relative_to(root)), "SHA256": FAA_SHA256, "Role": "93,782 protein queries"},
        {"Asset": "Representative metadata", "Source": str(metadata_gz.relative_to(root)), "SHA256": METADATA_SHA256, "Role": "ORF completeness strata"},
        {"Asset": "Gene abundance", "Source": str(abundance_gz.relative_to(root)), "SHA256": ABUNDANCE_SHA256, "Role": "read-weighted profiles"},
        {"Asset": "B. thetaiotaomicron VPI-5482", "Source": str(btheta_gz.relative_to(root)), "SHA256": BTHETA_GZ_SHA256, "Role": "real CGC/substrate control"},
    ])
    contract = {
        "article": 37, "seed": SEED, "catalog_proteins": CATALOG_GENES,
        "dbcan": "5.2.9", "database_release": "db_v5-2-9_5-5-2026",
        "primary_call": "at least two of diamond,hmm,dbCANsub",
        "single_tool_branch": "sensitivity only", "threads_catalog": 32,
        "hmm_evalue": 1e-15, "hmm_coverage": 0.35, "diamond_evalue": 1e-102,
        "cgc_control": "Bacteroides thetaiotaomicron VPI-5482 / GCA_000011065.1",
        "catalog_faa_sha256": FAA_SHA256, "metadata_sha256": METADATA_SHA256,
        "abundance_sha256": ABUNDANCE_SHA256, "btheta_canonical_sha256": canonical,
    }
    (work / "run-contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    (work / ".article37-inputs-complete").write_text(json.dumps(contract, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
