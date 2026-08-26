#!/usr/bin/env python3
"""Download, extract, and verify Article 60 gapseq/CarveMe model assets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import tarfile
from pathlib import Path


GAPSEQ_ASSETS = ("Bacteria.tar.gz", "Archaea.tar.gz", "md5sums.txt")
CARVEME_PATHS = {
    "universe_bacteria.xml.gz": "data/generated/universe_bacteria.xml.gz",
    "universe_archaea.xml.gz": "data/generated/universe_archaea.xml.gz",
    "bigg_proteins.dmnd": "data/generated/bigg_proteins.dmnd",
    "media_db.tsv": "data/input/media_db.tsv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--carveme-env", type=Path, required=True)
    parser.add_argument("mode", choices=("download", "verify"))
    return parser.parse_args()


def file_digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 7 or len({row["Asset"] for row in rows}) != 7:
        raise RuntimeError("Article 60 database manifest must contain seven unique assets")
    return {row["Asset"]: row for row in rows}


def run(command: list[str]) -> None:
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(command)}")


def verify_file(path: Path, row: dict[str, str]) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    observed_bytes = path.stat().st_size
    observed_sha = file_digest(path, "sha256")
    if observed_bytes != int(row["ExpectedBytes"]) or observed_sha != row["SHA256"]:
        raise RuntimeError(
            f"Identity gate failed for {path}: bytes={observed_bytes}, sha256={observed_sha}"
        )
    expected_md5 = row["PublisherMD5"]
    if expected_md5 != "NA":
        observed_md5 = file_digest(path, "md5")
        if observed_md5 != expected_md5:
            raise RuntimeError(
                f"Publisher MD5 gate failed for {path}: expected={expected_md5}, observed={observed_md5}"
            )


def download_asset(target: Path, row: dict[str, str]) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_file():
        verify_file(target, row)
        return
    part = target.with_name(target.name + ".part")
    run(
        [
            "curl", "--fail", "--location", "--retry", "50", "--retry-all-errors",
            "--retry-delay", "5", "--connect-timeout", "30", "--speed-time", "180",
            "--speed-limit", "1024", "--continue-at", "-", "--output", str(part),
            row["Source"],
        ]
    )
    verify_file(part, row)
    part.replace(target)


def safe_extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with tarfile.open(archive, "r:gz") as handle:
        for member in handle.getmembers():
            target = (destination / member.name).resolve()
            if root != target and root not in target.parents:
                raise RuntimeError(f"Unsafe archive member: {member.name}")
            if member.issym() or member.islnk():
                raise RuntimeError(f"Links are not accepted in database archives: {member.name}")
        handle.extractall(destination, filter="data")


def nested_checksums(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        fields = raw.split(maxsplit=1)
        if len(fields) != 2:
            raise RuntimeError(f"Malformed publisher checksum row: {raw}")
        rows.append((fields[0], fields[1].removeprefix("./")))
    if len(rows) != 6:
        raise RuntimeError(f"Expected six nested database checksums, observed {len(rows)}")
    return rows


def extract_gapseq(cache: Path, rows: dict[str, dict[str, str]]) -> None:
    base = cache / "gapseq-seqdb-v1.5"
    archives = base / "archives"
    seqdb = base / "seqdb"
    for asset in ("Bacteria.tar.gz", "Archaea.tar.gz"):
        safe_extract(archives / asset, seqdb)
    checksums = nested_checksums(archives / "md5sums.txt")
    for expected, relative in checksums:
        nested = seqdb / relative
        if not nested.is_file() or file_digest(nested, "md5") != expected:
            raise RuntimeError(f"Nested gapseq archive identity gate failed: {relative}")
        safe_extract(nested, nested.parent)
    stamp = {
        "zenodoID": 20446806,
        "version": "1.5",
        "date": "2026-05-29",
        "archive_sha256": {
            domain: rows[f"{domain}.tar.gz"]["SHA256"]
            for domain in ("Bacteria", "Archaea")
        },
    }
    for domain in ("Bacteria", "Archaea"):
        (seqdb / domain / "version_seqDB.json").write_text(
            json.dumps(stamp, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )


def carveme_root(environment: Path) -> Path:
    candidates = sorted(
        {path.resolve() for path in (environment / "lib").glob("python*/site-packages/carveme")}
    )
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one CarveMe package root, observed {candidates}")
    return candidates[0]


def verify_all(
    cache: Path, environment: Path, rows: dict[str, dict[str, str]]
) -> dict[str, object]:
    base = cache / "gapseq-seqdb-v1.5"
    archives = base / "archives"
    for asset in GAPSEQ_ASSETS:
        verify_file(archives / asset, rows[asset])
        print(f"PASS\tgapseq\t{asset}")

    seqdb = base / "seqdb"
    nested = nested_checksums(archives / "md5sums.txt")
    nested_passed = 0
    fasta_counts: dict[str, int] = {}
    for expected, relative in nested:
        path = seqdb / relative
        if not path.is_file() or file_digest(path, "md5") != expected:
            raise RuntimeError(f"Nested archive checksum mismatch: {path}")
        nested_passed += 1
    for domain in ("Bacteria", "Archaea"):
        version = seqdb / domain / "version_seqDB.json"
        metadata = json.loads(version.read_text(encoding="utf-8"))
        if metadata.get("zenodoID") != 20446806 or metadata.get("version") != "1.5":
            raise RuntimeError(f"Wrong sequence database stamp: {version}")
        fasta_counts[domain] = sum(1 for _ in (seqdb / domain).rglob("*.fasta"))
        if fasta_counts[domain] == 0:
            raise RuntimeError(f"No extracted FASTA files for {domain}")

    package_root = carveme_root(environment)
    for asset, relative in CARVEME_PATHS.items():
        verify_file(package_root / relative, rows[asset])
        print(f"PASS\tCarveMe\t{asset}")
    result = {
        "article": 60,
        "gapseq_release": "2.1.0",
        "gapseq_sequence_db_version": "1.5",
        "gapseq_zenodo_record": 20446806,
        "gapseq_archive_assets": len(GAPSEQ_ASSETS),
        "gapseq_nested_archives": nested_passed,
        "gapseq_extracted_fasta": fasta_counts,
        "carveme_release": "1.6.6",
        "carveme_assets": len(CARVEME_PATHS),
    }
    (base / "database-audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    cache = args.cache_dir.resolve()
    environment = args.carveme_env.resolve()
    rows = read_manifest(root / "data/small/60-gem-database-manifest.tsv")
    if args.mode == "download":
        archive_dir = cache / "gapseq-seqdb-v1.5/archives"
        for asset in GAPSEQ_ASSETS:
            download_asset(archive_dir / asset, rows[asset])
        extract_gapseq(cache, rows)
    verify_all(cache, environment, rows)


if __name__ == "__main__":
    main()
