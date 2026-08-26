#!/usr/bin/env python3
"""Download and verify the locked AGORA/MICOM resources for Article 61."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import zipfile
from pathlib import Path

import pandas as pd


FILES = {
    "AGORA201RefSeq216Species": "agora201_refseq216_species_1.qza",
    "AGORA2RefSeqSpeciesManifest": "agora2_refseq_species.tsv",
    "WesternDietGutAGORA": "western_diet_gut_agora.qza",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("mode", choices=("download", "verify"))
    return parser.parse_args()


def digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if set(FILES) != {row["Asset"] for row in rows}:
        raise RuntimeError("Article 61 resource manifest identity mismatch")
    return {row["Asset"]: row for row in rows}


def verify_file(path: Path, row: dict[str, str]) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    result = {
        "asset": row["Asset"],
        "file": path.name,
        "bytes": path.stat().st_size,
        "md5": digest(path, "md5"),
        "sha256": digest(path, "sha256"),
    }
    if result["bytes"] != int(row["ExpectedBytes"]):
        raise RuntimeError(f"Byte-count mismatch: {path}")
    if result["sha256"] != row["SHA256"]:
        raise RuntimeError(f"SHA-256 mismatch: {path}")
    if row["PublisherMD5"] != "NA" and result["md5"] != row["PublisherMD5"]:
        raise RuntimeError(f"Publisher MD5 mismatch: {path}")
    return result


def download(path: Path, row: dict[str, str]) -> None:
    if path.is_file():
        verify_file(path, row)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    command = [
        "curl", "--fail", "--location", "--retry", "50", "--retry-all-errors",
        "--retry-delay", "5", "--connect-timeout", "30", "--speed-time", "180",
        "--speed-limit", "1024", "--continue-at", "-", "--output", str(part),
        row["Source"],
    ]
    subprocess.run(command, check=True)
    verify_file(part, row)
    part.replace(path)


def zip_integrity(path: Path) -> None:
    with zipfile.ZipFile(path) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise RuntimeError(f"Corrupt archive member in {path}: {bad}")


def qza_manifest(path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as archive:
        members = [name for name in archive.namelist() if name.endswith("/data/manifest.csv")]
        if len(members) != 1:
            raise RuntimeError(f"Expected one model manifest in {path}")
        return pd.read_csv(archive.open(members[0]))


def qza_medium(path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as archive:
        members = [name for name in archive.namelist() if name.endswith("/data/medium.csv")]
        if len(members) != 1:
            raise RuntimeError(f"Expected one medium table in {path}")
        return pd.read_csv(archive.open(members[0]))


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    cache = args.cache_dir.resolve()
    rows = read_manifest(root / "data/small/61-community-database-manifest.tsv")
    if args.mode == "download":
        for asset, name in FILES.items():
            download(cache / name, rows[asset])

    audits = [verify_file(cache / name, rows[asset]) for asset, name in FILES.items()]
    model_path = cache / FILES["AGORA201RefSeq216Species"]
    medium_path = cache / FILES["WesternDietGutAGORA"]
    zip_integrity(model_path)
    zip_integrity(medium_path)
    models = qza_manifest(model_path)
    source = pd.read_csv(cache / FILES["AGORA2RefSeqSpeciesManifest"], sep="\t", low_memory=False)
    medium = qza_medium(medium_path)
    if len(models) != 1746 or models.species.nunique() != 1746:
        raise RuntimeError("AGORA QIIME artifact does not contain 1746 unique species models")
    if len(source) != 7302 or source.species.nunique() != 1746:
        raise RuntimeError("Locked AGORA source manifest count mismatch")
    if len(medium) != 171 or medium.reaction.duplicated().any():
        raise RuntimeError("Locked Western-diet medium contract mismatch")
    result = {
        "article": 61,
        "assets": audits,
        "agora_models": len(models),
        "agora_species": models.species.nunique(),
        "agora_source_strains": len(source),
        "medium_reactions": len(medium),
        "positive_medium_reactions": int(medium.flux.gt(0).sum()),
    }
    (cache / "database-audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
