#!/usr/bin/env python3
"""Identity-gate Article 38 catalog, co-assembly, controls, and CARD 4.0.1."""

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


SEED = 20260738
CATALOG_GENES = 93_782
INPUT_HASHES = {
    "catalog_faa": "3db88ff78a548dddfc48caa8a17f04bcbb58dcafe2345d9dd31bc4e12f2a3569",
    "catalog_fna": "56f0be1fa7230517318dd745deba55be204da473ee3b6abbc24bd56ccaf3ceb6",
    "metadata": "677479da11ef41a6f27f11798c38dfd9b5830b5c564d9a8f436951413acd7c09",
    "abundance": "8dc7ea7c1f0a6e61625fd1aada492d252b60eed050a91af6560986987e43236b",
    "coassembly": "904f92521ff0ce9f12bd52d153bb249ec816fc900051e06b4b12bc5da74a270a",
    "pseudomonas": "698d281e6146e58b763e1f3ca8999258f7a2a0661bef69ef6d0959ed8ad9768d",
    "staphylococcus": "dedc519dfa13b276676cece80acae4949ee37a532bfd25a4a6af6988c55188d1",
}
CARD_ARCHIVE_BYTES = 4_620_455
CARD_ARCHIVE_SHA256 = "d80520c51ef3bae1098cd00e0c3d1d013e5785f8d8df787dcc5f435927557df9"
CARD_JSON_SHA256 = "dee4dcdb0d9c7f79107452d64211d816d1eab55289ddf8dc5f1e99ddfdc5e111"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--database-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_stats(path: Path) -> tuple[int, int]:
    opener = gzip.open if path.suffix == ".gz" else open
    records = bases = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            else:
                bases += len(line.strip().rstrip("*"))
    return records, bases


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def tool_version(command: list[str], env: dict[str, str]) -> str:
    result = subprocess.run(command, env=env, text=True, capture_output=True, check=True)
    return " | ".join(part.strip() for part in (result.stdout, result.stderr) if part.strip()).replace("\t", " ")


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    database = args.database_root.resolve()
    env_prefix = args.env_prefix.resolve()
    paths = {
        "catalog_faa": root / "data/small/34-nonredundant-gene-catalog-frozen/catalog/megahit-mix-primary.faa.gz",
        "catalog_fna": root / "data/small/34-nonredundant-gene-catalog-frozen/catalog/megahit-mix-primary.fna.gz",
        "metadata": root / "data/small/34-nonredundant-gene-catalog-frozen/primary-catalog-representatives.tsv.gz",
        "abundance": root / "data/small/35-gene-abundance-frozen/gene-abundance-long.tsv.gz",
        "coassembly": root / "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz",
        "pseudomonas": root / "data/raw/article32/benchmark_mock/reference/all_genomes_listed/Pseudomonas_aeruginosa_ATCC_9027.fna.gz",
        "staphylococcus": root / "data/raw/article32/benchmark_mock/reference/all_genomes_listed/Staphylococcus_aureus_USA300_FPR3757.fna.gz",
    }
    for key, path in paths.items():
        if not path.is_file() or sha256(path) != INPUT_HASHES[key]:
            raise ValueError(f"Missing or checksum-mismatched {key}: {path}")

    archive = database / "archive/broadstreet-v4.0.1.tar.bz2"
    card_json = database / "rgi-db/card.json"
    loaded = database / "rgi-db/loaded_databases.json"
    if archive.stat().st_size != CARD_ARCHIVE_BYTES or sha256(archive) != CARD_ARCHIVE_SHA256:
        raise ValueError("CARD archive identity failed")
    if not card_json.is_file() or sha256(card_json) != CARD_JSON_SHA256:
        raise ValueError("Loaded CARD card.json identity failed")
    card = json.loads(card_json.read_text(encoding="utf-8"))
    if card.get("_version") != "4.0.1":
        raise ValueError(f"Expected CARD 4.0.1, observed {card.get('_version')}")
    if not loaded.is_file() or json.loads(loaded.read_text(encoding="utf-8"))["card_json"]["data_version"] != "4.0.1":
        raise ValueError("RGI loaded-database ledger failed")

    if work.exists():
        shutil.rmtree(work)
    for child in ("inputs", "catalog", "coassembly", "pseudomonas", "staphylococcus", "summary", "logs", "tmp"):
        (work / child).mkdir(parents=True, exist_ok=True)
    for source_key, target_name in (
        ("catalog_faa", "catalog.faa"),
        ("catalog_fna", "catalog.fna"),
        ("coassembly", "coassembly.fna"),
        ("pseudomonas", "pseudomonas.fna"),
        ("staphylococcus", "staphylococcus.fna"),
    ):
        with gzip.open(paths[source_key], "rt", encoding="utf-8") as source, (work / "inputs" / target_name).open("w", encoding="utf-8") as target:
            shutil.copyfileobj(source, target)
    stats = {key: fasta_stats(path) for key, path in paths.items() if key in {"catalog_faa", "catalog_fna", "coassembly", "pseudomonas", "staphylococcus"}}
    if stats["catalog_faa"][0] != CATALOG_GENES or stats["catalog_fna"][0] != CATALOG_GENES:
        raise ValueError(f"Catalog record count failed: {stats}")

    runtime = os.environ.copy()
    runtime.update({
        "PATH": f"{env_prefix / 'bin'}:/usr/bin:/bin",
        "DATA_PATH": str(database / "rgi-db"),
        "MPLCONFIGDIR": str(work / "tmp/matplotlib"),
        "XDG_CACHE_HOME": str(work / "tmp/xdg"),
        "PYTHONHASHSEED": "0", "PYTHONNOUSERSITE": "1", "PYTHONPATH": "",
    })
    (work / "tmp/matplotlib").mkdir(parents=True, exist_ok=True)
    (work / "tmp/xdg/fontconfig").mkdir(parents=True, exist_ok=True)
    rgi = env_prefix / "bin/rgi"
    if not rgi.is_file():
        raise FileNotFoundError(rgi)
    tool_rows = [
        {"Tool": "RGI", "VersionEvidence": tool_version([str(rgi), "main", "-v"], runtime)},
        {"Tool": "DIAMOND", "VersionEvidence": tool_version([str(env_prefix / "bin/diamond"), "version"], runtime)},
        {"Tool": "BLAST+", "VersionEvidence": tool_version([str(env_prefix / "bin/blastp"), "-version"], runtime)},
        {"Tool": "Pyrodigal", "VersionEvidence": tool_version([str(env_prefix / "bin/python"), "-c", "import pyrodigal; print(pyrodigal.__version__)"], runtime)},
    ]
    input_rows = []
    for key in ("catalog_faa", "catalog_fna", "coassembly", "pseudomonas", "staphylococcus"):
        records, bases = stats[key]
        input_rows.append({"Asset": key, "Records": records, "ResiduesOrBases": bases, "CompressedSHA256": INPUT_HASHES[key]})
    write_tsv(work / "input-audit.tsv", input_rows)
    write_tsv(work / "database-audit.tsv", [
        {"Asset": "CARD canonical archive", "Release": "4.0.1", "Bytes": CARD_ARCHIVE_BYTES, "SHA256": CARD_ARCHIVE_SHA256, "Status": "PASS"},
        {"Asset": "CARD card.json", "Release": card["_version"], "Bytes": card_json.stat().st_size, "SHA256": CARD_JSON_SHA256, "Status": "PASS"},
    ])
    write_tsv(work / "tool-versions.tsv", tool_rows)
    write_tsv(work / "input-lineage.tsv", [
        {"Asset": "Catalog proteins", "Source": str(paths["catalog_faa"].relative_to(root)), "Identifier": "Article 34 primary catalog", "Role": "gene-level ARG calls and Article 35 abundance join"},
        {"Asset": "Co-assembly", "Source": str(paths["coassembly"].relative_to(root)), "Identifier": "Article 30 MEGAHIT co-assembly >=1 kb", "Role": "nucleotide/model/context sensitivity"},
        {"Asset": "P. aeruginosa ATCC 9027", "Source": str(paths["pseudomonas"].relative_to(root)), "Identifier": "GCA_002563335.1; benchmark repository commit a429a372", "Role": "Gram-negative positive control"},
        {"Asset": "S. aureus USA300 FPR3757", "Source": str(paths["staphylococcus"].relative_to(root)), "Identifier": "GCA_000013465.1; benchmark repository commit a429a372", "Role": "Gram-positive positive control"},
        {"Asset": "CARD", "Source": "https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2", "Identifier": "CARD 4.0.1", "Role": "canonical resistance ontology and model thresholds"},
    ])
    contract = {
        "article": 38, "seed": SEED, "rgi": "6.0.8", "card": "4.0.1",
        "catalog_proteins": CATALOG_GENES, "alignment": "DIAMOND",
        "primary_tiers": ["Perfect", "Strict"], "sensitivity_tier": "Loose",
        "include_loose": True, "include_nudge": False,
        "coassembly_low_quality": True, "orf_finder": "PYRODIGAL",
        "catalog_faa_sha256": INPUT_HASHES["catalog_faa"],
        "coassembly_sha256": INPUT_HASHES["coassembly"],
        "card_archive_sha256": CARD_ARCHIVE_SHA256,
        "card_json_sha256": CARD_JSON_SHA256,
    }
    (work / "run-contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    (work / ".article38-inputs-complete").write_text(json.dumps(contract, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
