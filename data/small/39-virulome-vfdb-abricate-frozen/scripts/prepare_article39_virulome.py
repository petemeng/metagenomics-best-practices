#!/usr/bin/env python3
"""Identity-gate Article 39 inputs and official VFDB 2026-07-24 snapshots."""

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


SEED = 20260739
CATALOG_GENES = 93_782
HASHES = {
    "catalog_fna": "56f0be1fa7230517318dd745deba55be204da473ee3b6abbc24bd56ccaf3ceb6",
    "metadata": "677479da11ef41a6f27f11798c38dfd9b5830b5c564d9a8f436951413acd7c09",
    "abundance": "8dc7ea7c1f0a6e61625fd1aada492d252b60eed050a91af6560986987e43236b",
    "coassembly": "904f92521ff0ce9f12bd52d153bb249ec816fc900051e06b4b12bc5da74a270a",
    "pseudomonas": "698d281e6146e58b763e1f3ca8999258f7a2a0661bef69ef6d0959ed8ad9768d",
    "staphylococcus": "dedc519dfa13b276676cece80acae4949ee37a532bfd25a4a6af6988c55188d1",
    "vfdb_core_raw": "0bc8a522660bc13b1891f41cccc8eb30d3dfe0f13fb236ef9dc1245365113afb",
    "vfdb_full_raw": "1dd2df13b77d791bd29ec8aefb5eaaabfd3d9007ba751545e91e5516edcdef60",
    "vfdb_core_formatted": "62d721cc53e7694b78f02d8f58287ed1563c2db943b4c08f9e46414516be761b",
    "vfdb_full_formatted": "00bd797524195ec5be0ff44b69af279e8abd21f27e273ee4f02ff5ec1d3e16dc",
}
RAW_BYTES = {"vfdb_core_raw": 2_011_507, "vfdb_full_raw": 10_683_762}
FORMATTED_BYTES = {"vfdb_core_formatted": 7_094_381, "vfdb_full_formatted": 55_969_046}


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
    with opener(path, "rt", encoding="latin-1") as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            else:
                bases += len(line.strip())
    return records, bases


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def version(command: list[str], env: dict[str, str]) -> str:
    result = subprocess.run(command, env=env, text=True, capture_output=True, check=True)
    return " | ".join(part.strip() for part in (result.stdout, result.stderr) if part.strip()).replace("\t", " ")


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    database = args.database_root.resolve()
    env_prefix = args.env_prefix.resolve()
    paths = {
        "catalog_fna": root / "data/small/34-nonredundant-gene-catalog-frozen/catalog/megahit-mix-primary.fna.gz",
        "metadata": root / "data/small/34-nonredundant-gene-catalog-frozen/primary-catalog-representatives.tsv.gz",
        "abundance": root / "data/small/35-gene-abundance-frozen/gene-abundance-long.tsv.gz",
        "coassembly": root / "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz",
        "pseudomonas": root / "data/raw/article32/benchmark_mock/reference/all_genomes_listed/Pseudomonas_aeruginosa_ATCC_9027.fna.gz",
        "staphylococcus": root / "data/raw/article32/benchmark_mock/reference/all_genomes_listed/Staphylococcus_aureus_USA300_FPR3757.fna.gz",
        "vfdb_core_raw": database / "raw/VFDB_setA_nt.fas.gz",
        "vfdb_full_raw": database / "raw/VFDB_setB_nt.fas.gz",
        "vfdb_core_formatted": database / "abricate/vfdb-core/sequences",
        "vfdb_full_formatted": database / "abricate/vfdb-full/sequences",
    }
    for key, path in paths.items():
        if not path.is_file() or sha256(path) != HASHES[key]:
            raise ValueError(f"Missing or checksum-mismatched {key}: {path}")
    for key, expected in {**RAW_BYTES, **FORMATTED_BYTES}.items():
        if paths[key].stat().st_size != expected:
            raise ValueError(f"Byte-count mismatch for {key}")
    core_records, core_bases = fasta_stats(paths["vfdb_core_raw"])
    full_records, full_bases = fasta_stats(paths["vfdb_full_raw"])
    if (core_records, full_records) != (4_895, 35_188):
        raise ValueError("Unexpected VFDB sequence counts")
    for label, expected_count in (("vfdb-core", 4_895), ("vfdb-full", 35_188)):
        manifest = json.loads((database / f"abricate/{label}/build-manifest.json").read_text(encoding="utf-8"))
        if manifest["database"] != label or manifest["sequences"] != expected_count:
            raise ValueError(f"VFDB build manifest failed: {label}")

    if work.exists():
        shutil.rmtree(work)
    for child in ("inputs", "catalog-core-primary", "catalog-core-sensitive", "catalog-full-primary", "coassembly", "pseudomonas", "staphylococcus", "summary", "logs", "tmp"):
        (work / child).mkdir(parents=True, exist_ok=True)
    for source_key, target_name in (
        ("catalog_fna", "catalog.fna"), ("coassembly", "coassembly.fna"),
        ("pseudomonas", "pseudomonas.fna"), ("staphylococcus", "staphylococcus.fna"),
    ):
        with gzip.open(paths[source_key], "rt", encoding="utf-8") as source, (work / "inputs" / target_name).open("w", encoding="utf-8") as target:
            shutil.copyfileobj(source, target)
    stats = {key: fasta_stats(paths[key]) for key in ("catalog_fna", "coassembly", "pseudomonas", "staphylococcus")}
    if stats["catalog_fna"][0] != CATALOG_GENES:
        raise ValueError("Catalog nucleotide record count failed")

    runtime = os.environ.copy()
    runtime.update({"PATH": f"{env_prefix / 'bin'}:/usr/bin:/bin", "LC_ALL": "C", "PYTHONHASHSEED": "0"})
    abricate = env_prefix / "bin/abricate"
    tool_rows = [
        {"Tool": "ABRicate", "VersionEvidence": version([str(abricate), "--version"], runtime)},
        {"Tool": "BLAST+", "VersionEvidence": version([str(env_prefix / "bin/blastn"), "-version"], runtime)},
    ]
    write_tsv(work / "tool-versions.tsv", tool_rows)
    write_tsv(work / "input-audit.tsv", [
        {"Asset": key, "Records": stats[key][0], "Bases": stats[key][1], "CompressedSHA256": HASHES[key]}
        for key in ("catalog_fna", "coassembly", "pseudomonas", "staphylococcus")
    ])
    write_tsv(work / "database-audit.tsv", [
        {"Database": "VFDB core set A", "Release": "2026-07-24", "Sequences": core_records, "Bases": core_bases, "RawBytes": RAW_BYTES["vfdb_core_raw"], "RawSHA256": HASHES["vfdb_core_raw"], "FormattedSHA256": HASHES["vfdb_core_formatted"], "Status": "PASS"},
        {"Database": "VFDB full set B", "Release": "2026-07-24", "Sequences": full_records, "Bases": full_bases, "RawBytes": RAW_BYTES["vfdb_full_raw"], "RawSHA256": HASHES["vfdb_full_raw"], "FormattedSHA256": HASHES["vfdb_full_formatted"], "Status": "PASS"},
    ])
    write_tsv(work / "input-lineage.tsv", [
        {"Asset": "Gene catalog nucleotide FASTA", "Identifier": "Article 34 / 93,782 representatives", "Role": "gene-level virulence calls and abundance join"},
        {"Asset": "MEGAHIT co-assembly", "Identifier": "Article 30 / contigs >=1 kb", "Role": "contig-context branch"},
        {"Asset": "P. aeruginosa ATCC 9027", "Identifier": "GCA_002563335.1; repository commit a429a372", "Role": "Gram-negative positive control"},
        {"Asset": "S. aureus USA300 FPR3757", "Identifier": "GCA_000013465.1; repository commit a429a372", "Role": "Gram-positive positive control"},
        {"Asset": "VFDB set A", "Identifier": "https://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz", "Role": "primary curated core reference"},
        {"Asset": "VFDB set B", "Identifier": "https://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz", "Role": "pre-specified database sensitivity branch"},
    ])
    contract = {
        "article": 39, "seed": SEED, "abricate": "1.4.0", "blast": "2.17.0",
        "vfdb_release": "2026-07-24", "catalog_genes": CATALOG_GENES,
        "primary": {"database": "vfdb-core", "minimum_identity_percent": 90, "minimum_reference_coverage_percent": 80},
        "threshold_sensitivity": {"database": "vfdb-core", "minimum_identity_percent": 80, "minimum_reference_coverage_percent": 80},
        "database_sensitivity": {"database": "vfdb-full", "minimum_identity_percent": 90, "minimum_reference_coverage_percent": 80},
        "catalog_fna_sha256": HASHES["catalog_fna"], "coassembly_sha256": HASHES["coassembly"],
        "vfdb_core_sha256": HASHES["vfdb_core_raw"], "vfdb_full_sha256": HASHES["vfdb_full_raw"],
    }
    (work / "run-contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    (work / ".article39-inputs-complete").write_text(json.dumps(contract, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

