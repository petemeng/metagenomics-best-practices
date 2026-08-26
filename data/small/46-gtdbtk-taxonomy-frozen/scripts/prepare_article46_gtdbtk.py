#!/usr/bin/env python3
"""Prepare checksum-identified Article 45 representatives for GTDB-Tk R232."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
from pathlib import Path

from article41_44_utils import dump_json, fasta_summary, read_tsv, sha256, write_tsv


def cluster_key(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in value.split("_"))


def verify_manifest(root: Path) -> None:
    manifest = root / "file-checksums.sha256"
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        path = root / relative
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Frozen checksum mismatch: {relative}")


def installed_size(path: Path) -> tuple[int, int]:
    files, total = 0, 0
    for directory, _, names in os.walk(path):
        for name in names:
            candidate = Path(directory) / name
            if candidate.is_file():
                files += 1
                total += candidate.stat().st_size
    return files, total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--database-dir", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    database, archive = args.database_dir.resolve(), args.archive.resolve()
    frozen45 = root / "data/small/45-drep-dereplication-frozen"
    verify_manifest(frozen45)

    manifest_path = root / "data/small/46-gtdbtk-database-manifest.tsv"
    manifest_rows = read_tsv(manifest_path)
    if len(manifest_rows) != 1:
        raise ValueError("Article 46 database manifest must contain one row")
    manifest = manifest_rows[0]
    if manifest["release_id"] != "R232" or manifest["tool_version"] != "2.7.2":
        raise ValueError("Article 46 requires GTDB-Tk 2.7.2 with R232")
    sentinel = Path(f"{archive}.{manifest['checksum_algorithm']}.verified")
    if not archive.is_file() or archive.stat().st_size != int(manifest["expected_compressed_bytes"]):
        raise ValueError("GTDB R232 archive byte count is not verified")
    if not sentinel.is_file():
        raise FileNotFoundError(f"Run db/download_db.sh verify gtdbtk-r232 first: {sentinel}")
    sentinel_values = {}
    with sentinel.open(encoding="utf-8") as handle:
        for line in handle:
            key, value = line.rstrip("\n").split("\t", 1)
            sentinel_values[key] = value
    if sentinel_values.get("checksum") != manifest["expected_checksum"] or sentinel_values.get("bytes") != manifest["expected_compressed_bytes"]:
        raise ValueError("GTDB R232 archive verification sentinel disagrees with manifest")
    expected_dirs = {
        "markers", "masks", "metadata", "mrca_red", "msa",
        "pplacer", "radii", "skani", "split", "taxonomy",
    }
    observed_dirs = {path.name for path in database.iterdir() if path.is_dir()} if database.is_dir() else set()
    if not expected_dirs.issubset(observed_dirs):
        raise ValueError(f"GTDB R232 extraction incomplete: missing {sorted(expected_dirs - observed_dirs)}")
    database_files, database_bytes = installed_size(database)
    if database_bytes != int(manifest["expected_installed_bytes"]) or database_files != 263:
        raise ValueError(
            f"GTDB R232 installed-size audit failed: {database_files} files / {database_bytes} bytes"
        )
    if manifest["validation_status"] != "VERIFIED_LOCAL_CHECK_INSTALL_PASS":
        raise ValueError("GTDB R232 manifest does not record a passing gtdbtk check_install")

    genomes_dir = work / "inputs/genomes"
    if genomes_dir.exists():
        shutil.rmtree(genomes_dir)
    genomes_dir.mkdir(parents=True)
    membership = [
        row for row in read_tsv(frozen45 / "cluster-membership.tsv.gz")
        if row["Branch"] == "Species 95% ANI" and row["IsRepresentative"].lower() == "true"
    ]
    if len(membership) != 24:
        raise ValueError("Expected 24 Article 45 representatives")
    genomes: list[dict[str, object]] = []
    for index, row in enumerate(sorted(membership, key=lambda item: cluster_key(item["Cluster"])), start=1):
        sgb = f"SGB_{index:03d}"
        source = frozen45 / "representative-genomes" / f"{row['Genome']}.gz"
        target = genomes_dir / f"{sgb}.fna"
        with gzip.open(source, "rb") as input_handle, target.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle)
        if sha256(target) != row["SHA256"]:
            raise ValueError(f"Representative checksum mismatch: {row['Genome']}")
        fasta, _ = fasta_summary(target)
        genomes.append({
            "SGB": sgb,
            "dRepCluster": row["Cluster"],
            "Representative": row["Genome"],
            "RepresentativeSHA256": row["SHA256"],
            "GenomeBp": fasta["TotalBp"],
            "Contigs": fasta["Contigs"],
            "N50Bp": fasta["N50Bp"],
            "Completeness": float(row["Completeness"]),
            "Contamination": float(row["Contamination"]),
            "MIMAGQuality": row["MIMAGQuality"],
            "Path": str(target),
        })
    write_tsv(work / "genome-ledger.tsv", genomes)
    write_tsv(work / "database-audit.tsv", [{
        "DatabaseID": manifest["database_id"],
        "Tool": manifest["tool"],
        "ToolVersion": manifest["tool_version"],
        "Release": manifest["release_id"],
        "ReleaseDate": manifest["release_date"],
        "Archive": str(archive),
        "ArchiveBytes": archive.stat().st_size,
        "ChecksumAlgorithm": manifest["checksum_algorithm"],
        "ArchiveChecksum": sentinel_values["checksum"],
        "InstalledPath": str(database),
        "InstalledFiles": database_files,
        "InstalledBytes": database_bytes,
        "ManifestStatus": manifest["validation_status"],
        "LocalStatus": "ARCHIVE_MD5_BYTES_EXTRACTION_AND_GTDBTK_CHECK_INSTALL_PASS",
    }])
    write_tsv(work / "input-lineage.tsv", [{
        "Output": "Article 46 GTDB-Tk input catalog",
        "ImmediateInput": "24 Article 45 95%-ANI representatives",
        "Transformation": "checksum-exact decompression and stable SGB ID assignment",
        "TruthUsed": "No",
        "Evidence": "genome-ledger.tsv",
    }])
    dump_json(work / "run-contract.json", {
        "article": 46,
        "seed": 20260746,
        "random_process": False,
        "input_genomes": 24,
        "input_order": [row["SGB"] for row in genomes],
        "gtdbtk_version": "2.7.2",
        "gtdb_release": "R232",
        "full_tree": False,
        "minimum_percent_aa": 10,
        "minimum_alignment_fraction": 0.5,
        "pplacer_cpus": 1,
        "write_single_copy_genes": True,
        "keep_intermediates": True,
        "taxonomy_source": "independent GTDB-Tk classification",
        "truth_used_for_taxonomy": False,
        "database_path": str(database),
    })
    (work / ".article46-inputs-complete").write_text("complete\n", encoding="utf-8")
    print(f"Prepared {len(genomes)} SGBs against GTDB {manifest['release_id']}")


if __name__ == "__main__":
    main()
