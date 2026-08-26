#!/usr/bin/env python3
"""Prepare and identity-gate the real Article 36 eggNOG-mapper inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
from pathlib import Path


SEED = 20260736
CATALOG_GENES = 93_782
FAA_SHA256 = "3db88ff78a548dddfc48caa8a17f04bcbb58dcafe2345d9dd31bc4e12f2a3569"
METADATA_SHA256 = "677479da11ef41a6f27f11798c38dfd9b5830b5c564d9a8f436951413acd7c09"
ABUNDANCE_SHA256 = "8dc7ea7c1f0a6e61625fd1aada492d252b60eed050a91af6560986987e43236b"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--database-dir", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--database-manifest", type=Path)
    parser.add_argument(
        "--skip-large-hash",
        action="store_true",
        help="Development-only: check database byte counts and SQLite release without rehashing large files.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def command_version(command: list[str], env: dict[str, str]) -> str:
    completed = subprocess.run(command, env=env, text=True, capture_output=True, check=False)
    text = "\n".join(part.strip() for part in (completed.stdout, completed.stderr) if part.strip())
    if completed.returncode != 0:
        raise RuntimeError(f"Version command failed ({completed.returncode}): {' '.join(command)}\n{text}")
    return text.replace("\t", " ").replace("\n", " | ")


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    db_dir = args.database_dir.resolve()
    env_prefix = args.env_prefix.resolve()
    manifest_path = (
        args.database_manifest.resolve()
        if args.database_manifest
        else root / "data/small/36-eggnog-database-manifest.tsv"
    )
    faa_gz = root / "data/small/34-nonredundant-gene-catalog-frozen/catalog/megahit-mix-primary.faa.gz"
    metadata_gz = root / "data/small/34-nonredundant-gene-catalog-frozen/primary-catalog-representatives.tsv.gz"
    abundance_gz = root / "data/small/35-gene-abundance-frozen/gene-abundance-long.tsv.gz"

    for path in (faa_gz, metadata_gz, abundance_gz, manifest_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    expected_inputs = {
        faa_gz: FAA_SHA256,
        metadata_gz: METADATA_SHA256,
        abundance_gz: ABUNDANCE_SHA256,
    }
    for path, expected in expected_inputs.items():
        observed = sha256(path)
        if observed != expected:
            raise ValueError(f"Input checksum mismatch for {path}: {observed}")

    database_rows = read_tsv(manifest_path)
    if not database_rows:
        raise ValueError("Database manifest is empty")
    required_rows = [row for row in database_rows if row["RequiredForRun"] == "yes"]
    if len(required_rows) != 4:
        raise ValueError(f"Expected four required database assets, observed {len(required_rows)}")
    database_audit: list[dict[str, object]] = []
    for row in database_rows:
        installed = db_dir / row["InstalledFile"]
        observed_bytes = installed.stat().st_size if installed.is_file() else 0
        bytes_ok = observed_bytes == int(row["InstalledBytes"])
        observed_hash = "NOT_REHASHED"
        hash_ok = True
        if row["RequiredForRun"] == "yes" and not args.skip_large_hash:
            observed_hash = sha256(installed) if installed.is_file() else "MISSING"
            hash_ok = observed_hash == row["InstalledSHA256"]
        status = bytes_ok and hash_ok
        database_audit.append(
            {
                "Asset": row["Asset"],
                "Release": row["Release"],
                "InstalledFile": row["InstalledFile"],
                "ExpectedBytes": row["InstalledBytes"],
                "ObservedBytes": observed_bytes,
                "ExpectedSHA256": row["InstalledSHA256"],
                "ObservedSHA256": observed_hash,
                "RequiredForRun": row["RequiredForRun"],
                "Status": "PASS" if status else "FAIL",
            }
        )
        if not status:
            raise ValueError(f"Database identity gate failed: {row['Asset']}")

    with sqlite3.connect(db_dir / "eggnog.db") as connection:
        db_version = connection.execute("SELECT version FROM version").fetchone()[0]
    if db_version != "5.0.2":
        raise ValueError(f"Expected eggNOG 5.0.2, observed {db_version}")

    if work.exists():
        shutil.rmtree(work)
    (work / "inputs").mkdir(parents=True)
    (work / "annotation/main").mkdir(parents=True)
    (work / "annotation/go-all").mkdir(parents=True)
    (work / "summary").mkdir(parents=True)
    (work / "logs").mkdir(parents=True)
    (work / "tmp").mkdir(parents=True)

    catalog_faa = work / "inputs/catalog.faa"
    ids: list[str] = []
    seen: set[str] = set()
    sequence_count = 0
    residue_count = 0
    terminal_stop_count = 0
    current_has_sequence = False
    with gzip.open(faa_gz, "rt", encoding="utf-8") as source, catalog_faa.open("w", encoding="utf-8") as target:
        for line in source:
            target.write(line)
            if line.startswith(">"):
                if current_has_sequence:
                    sequence_count += 1
                identifier = line[1:].split(None, 1)[0]
                if identifier in seen:
                    raise ValueError(f"Duplicate FASTA identifier: {identifier}")
                seen.add(identifier)
                ids.append(identifier)
                current_has_sequence = False
            else:
                sequence = line.strip()
                if sequence:
                    current_has_sequence = True
                    residue_count += len(sequence.rstrip("*"))
                    terminal_stop_count += int(sequence.endswith("*"))
        if current_has_sequence:
            sequence_count += 1
    if sequence_count != CATALOG_GENES or len(ids) != CATALOG_GENES:
        raise ValueError(f"Catalog FASTA count mismatch: {sequence_count}/{len(ids)}")

    metadata = read_tsv_gz(metadata_gz)
    metadata_ids = [row["RepresentativeID"] for row in metadata]
    if len(metadata) != CATALOG_GENES or set(metadata_ids) != seen:
        raise ValueError("Catalog FASTA and representative metadata IDs do not match")
    completeness = {label: 0 for label in ("Complete", "Partial", "Incomplete")}
    for row in metadata:
        completeness[row["Completeness"]] += 1

    abundance_rows = read_tsv_gz(abundance_gz)
    sample_gene_ids: dict[str, set[str]] = {}
    sample_counts: dict[str, int] = {}
    for row in abundance_rows:
        sample = row["Sample"]
        sample_gene_ids.setdefault(sample, set()).add(row["GeneID"])
        sample_counts[sample] = sample_counts.get(sample, 0) + int(row["RawCount"])
    if set(sample_gene_ids) != {"MOCK1", "MOCK2"}:
        raise ValueError(f"Unexpected abundance samples: {sorted(sample_gene_ids)}")
    for sample, gene_ids in sample_gene_ids.items():
        if gene_ids != seen:
            raise ValueError(f"{sample} abundance IDs do not match the catalog")

    runtime_env = os.environ.copy()
    runtime_env.update(
        {
            "PATH": f"{env_prefix / 'bin'}:/usr/bin:/bin",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "PYTHONPATH": "",
        }
    )
    emapper = env_prefix / "bin/emapper.py"
    python = env_prefix / "bin/python"
    diamond = env_prefix / "bin/diamond"
    for path in (emapper, python, diamond):
        if not path.is_file():
            raise FileNotFoundError(path)
    tool_rows = [
        {"Tool": "Python", "VersionEvidence": command_version([str(python), "--version"], runtime_env)},
        {"Tool": "DIAMOND", "VersionEvidence": command_version([str(diamond), "version"], runtime_env)},
        {
            "Tool": "eggNOG-mapper",
            "VersionEvidence": command_version(
                [str(emapper), "--version", "--data_dir", str(db_dir)], runtime_env
            ),
        },
    ]

    write_tsv(work / "database-audit.tsv", database_audit)
    write_tsv(work / "tool-versions.tsv", tool_rows)
    write_tsv(
        work / "input-audit.tsv",
        [
            {
                "CatalogProteins": sequence_count,
                "UniqueIDs": len(seen),
                "ResiduesExcludingTerminalStop": residue_count,
                "ProteinsWithTerminalStop": terminal_stop_count,
                "Complete": completeness["Complete"],
                "Partial": completeness["Partial"],
                "Incomplete": completeness["Incomplete"],
                "MOCK1RawReads": sample_counts["MOCK1"],
                "MOCK2RawReads": sample_counts["MOCK2"],
                "EggNOGRelease": db_version,
            }
        ],
    )
    lineage_rows = [
        {
            "Asset": "catalog-proteins",
            "Source": "Article 34 MEGAHIT mix primary representatives",
            "Bytes": faa_gz.stat().st_size,
            "SHA256": FAA_SHA256,
            "Role": "93,782 protein queries",
        },
        {
            "Asset": "catalog-metadata",
            "Source": "Article 34 primary representative metadata",
            "Bytes": metadata_gz.stat().st_size,
            "SHA256": METADATA_SHA256,
            "Role": "length and ORF completeness strata",
        },
        {
            "Asset": "gene-abundance",
            "Source": "Article 35 audited raw read counts",
            "Bytes": abundance_gz.stat().st_size,
            "SHA256": ABUNDANCE_SHA256,
            "Role": "read-weighted annotation missingness",
        },
    ]
    write_tsv(work / "input-lineage.tsv", lineage_rows)

    contract = {
        "article": 36,
        "seed": SEED,
        "catalog_genes": CATALOG_GENES,
        "eggnog_mapper": "2.1.15",
        "eggnog_database": db_version,
        "search": {
            "mode": "diamond",
            "diamond": "2.0.15",
            "sensitivity": "sensitive",
            "iterate": "yes",
            "evalue": 0.001,
            "seed_ortholog_evalue": 0.001,
            "identity_filter": None,
            "query_coverage_filter": None,
            "subject_coverage_filter": None,
            "outfmt_short": True,
        },
        "annotation": {
            "tax_scope": "auto",
            "tax_scope_mode": "inner_narrowest",
            "target_orthologs": "all",
            "go_evidence_primary": "non-electronic",
            "go_evidence_sensitivity": "all",
            "pfam_realign": "none",
        },
        "dark_matter_states": [
            "No seed ortholog",
            "Orthology only",
            "Broad/family only",
            "Specific identifier",
        ],
        "input_faa_sha256": FAA_SHA256,
        "metadata_sha256": METADATA_SHA256,
        "abundance_sha256": ABUNDANCE_SHA256,
        "database_manifest": manifest_path.name,
        "large_database_rehashed": not args.skip_large_hash,
    }
    (work / "run-contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    (work / ".article36-inputs-ready").write_text("ready\n", encoding="utf-8")
    print(json.dumps({"work_dir": str(work), "contract": contract, "input_audit": str(work / 'input-audit.tsv')}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
