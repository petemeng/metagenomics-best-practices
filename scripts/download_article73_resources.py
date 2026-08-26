#!/usr/bin/env python3
"""Download and lock the public-resource snapshot used by Article 73."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import time
import urllib.request
from pathlib import Path

from PIL import Image


SNAPSHOT_DATE = "2026-08-23"
URLS = {
    "mgnify-study.json": (
        "https://www.ebi.ac.uk/metagenomics/api/v2/studies/MGYS00000410"
    ),
    "mgnify-samples-page1.json": (
        "https://www.ebi.ac.uk/metagenomics/api/v2/studies/"
        "MGYS00000410/samples/?page=1&page_size=100"
    ),
    "mgnify-samples-page2.json": (
        "https://www.ebi.ac.uk/metagenomics/api/v2/studies/"
        "MGYS00000410/samples/?page=2&page_size=100"
    ),
    "mgnify-analyses-page1.json": (
        "https://www.ebi.ac.uk/metagenomics/api/v2/studies/"
        "MGYS00000410/analyses/?page=1&page_size=100"
    ),
    "mgnify-analyses-page2.json": (
        "https://www.ebi.ac.uk/metagenomics/api/v2/studies/"
        "MGYS00000410/analyses/?page=2&page_size=100"
    ),
    "mgnify-analyses-page3.json": (
        "https://www.ebi.ac.uk/metagenomics/api/v2/studies/"
        "MGYS00000410/analyses/?page=3&page_size=100"
    ),
    "mgnify-catalogues.json": (
        "https://www.ebi.ac.uk/metagenomics/api/v2/genomes/catalogues/"
        "?page=1&page_size=100"
    ),
    "mgnify-marine-v2.json": (
        "https://www.ebi.ac.uk/metagenomics/api/v2/genomes/catalogues/"
        "marine-v2-0"
    ),
    "gtdb-r232-version.txt": (
        "https://data.gtdb.ecogenomic.org/releases/release232/232.0/VERSION.txt"
    ),
    "gtdb-r232-release-notes.txt": (
        "https://data.gtdb.ecogenomic.org/releases/release232/232.0/"
        "RELEASE_NOTES.txt"
    ),
    "gtdb-r232-md5.txt": (
        "https://data.gtdb.ecogenomic.org/releases/release232/232.0/MD5SUM.txt"
    ),
    "gtdb-r226-version.txt": (
        "https://data.gtdb.ecogenomic.org/releases/release226/226.0/VERSION.txt"
    ),
    "gtdb-r226-release-notes.txt": (
        "https://data.gtdb.ecogenomic.org/releases/release226/226.0/"
        "RELEASE_NOTES.txt"
    ),
    "gem-readme.md": "https://portal.nersc.gov/GEM/README.md",
    "gem-genome-metadata.tsv": (
        "https://portal.nersc.gov/GEM/genomes/genome_metadata.tsv"
    ),
    "gem-figure1-original.png": (
        "https://media.springernature.com/full/springer-static/image/"
        "art%3A10.1038%2Fs41587-020-0718-6/MediaObjects/"
        "41587_2020_718_Fig1_HTML.png"
    ),
}

# Snapshot acquired on 2026-08-23. Every resource is checked before preparation.
LOCKS: dict[str, tuple[int, str]] = {
    "mgnify-study.json": (3_088, "610ae775d4a7958f44285df1a6d4fcbac29b47023376d4533ec1d8e3b092e385"),
    "mgnify-samples-page1.json": (474_961, "6d7dec16fde6bc83e9eccf3ea55d795d17b15fb9fa79439fe9e018d951b6e1e4"),
    "mgnify-samples-page2.json": (169_616, "b09c27ae61c1c19273f201889cb9b19be158d5e4a8ac240f6555472d759e7030"),
    "mgnify-analyses-page1.json": (60_273, "cc2525849f194c1eecbfec3e654fc7c378a93d84f06e0f1dc1311b67f13c85c7"),
    "mgnify-analyses-page2.json": (60_174, "5e2261673ef284c24bc7d8f315bbaeaad32d01f30d0cd3bfb73b73b6d461cdad"),
    "mgnify-analyses-page3.json": (29_460, "8ec9c246e980253b201619e6538c87fa974a72f3aa1f0c91a8406895bc82fee8"),
    "mgnify-catalogues.json": (29_351, "e2571384392ee5686ca807b9d02d1f9f9e8187e8b2319a1f599891ffd83bef40"),
    "mgnify-marine-v2.json": (1_463, "2d190ed7eab56b4c85fd17f3e37266d0966696295b482543fe8a3e5cdffee5ec"),
    "gtdb-r232-version.txt": (29, "7aaa7dca8b101daaab5635cad0895051ff427223292b04913f8791ad8c53d591"),
    "gtdb-r232-release-notes.txt": (1_363, "c6fb891abcbbec1ac753d1f9d8bc920a8adb060ad0f51b41b99094daee62f2f7"),
    "gtdb-r232-md5.txt": (5_482, "a5c2cc52b7d319e70bb678bbe52acc6d4a697dde8f6e130b2ce39e706ba8939d"),
    "gtdb-r226-version.txt": (29, "f1a34e5c882e437196f2823c2c34281100596f273660c8682da63b529f7d62de"),
    "gtdb-r226-release-notes.txt": (1_809, "c93851a64cc75c425e73fb9b0a472b3cd9f2546c034fb1e4c4ff1a2e7e4ae1f4"),
    "gem-readme.md": (2_850, "29fcef8f24cdf6ae2d3ee0905d3a5af18cfd3abcb478042985625bba75a1ad2a"),
    "gem-genome-metadata.tsv": (12_350_919, "fd0ad382e4ec9dbc07915333b6c2e4b53257f6d3a9f47aad7da1d2cad6e83e37"),
    "gem-figure1-original.png": (428_457, "7e36a5130f753401362c08d86892c7ce53977642318619f0cba27249dd87ec11"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, default=Path("db/evidence/article73"))
    parser.add_argument("--seed-dir", type=Path)
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, target: Path, attempts: int = 8) -> None:
    partial = target.with_suffix(target.suffix + ".part")
    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "User-Agent": (
                        "metagenomics-best-practices/1.0 "
                        "(public-resource snapshot; contact repository maintainer)"
                    )
                },
            )
            with urllib.request.urlopen(request, timeout=240) as response, partial.open("wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
            partial.replace(target)
            return
        except Exception:
            partial.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            time.sleep(min(2 * attempt, 16))


def seed(target: Path, seed_dir: Path | None) -> None:
    if target.is_file() or seed_dir is None:
        return
    candidate = seed_dir / target.name
    if candidate.is_file():
        shutil.copy2(candidate, target)


def validate_content(cache: Path) -> None:
    study = json.loads((cache / "mgnify-study.json").read_text())
    if study.get("accession") != "MGYS00000410":
        raise ValueError("Unexpected MGnify study")
    if set(study.get("ena_accessions", [])) != {"ERP001736", "PRJEB1787"}:
        raise ValueError("Unexpected MGnify/ENA study crosswalk")

    sample_pages = [
        json.loads((cache / f"mgnify-samples-page{page}.json").read_text())
        for page in (1, 2)
    ]
    if {payload.get("count") for payload in sample_pages} != {136}:
        raise ValueError("MGnify sample count changed")
    if [len(payload.get("items", [])) for payload in sample_pages] != [100, 36]:
        raise ValueError("Unexpected MGnify sample pagination")

    analysis_pages = [
        json.loads((cache / f"mgnify-analyses-page{page}.json").read_text())
        for page in (1, 2, 3)
    ]
    if {payload.get("count") for payload in analysis_pages} != {249}:
        raise ValueError("MGnify analysis count changed")
    if [len(payload.get("items", [])) for payload in analysis_pages] != [100, 100, 49]:
        raise ValueError("Unexpected MGnify analysis pagination")

    catalogues = json.loads((cache / "mgnify-catalogues.json").read_text())
    if catalogues.get("count", 0) < 15 or len(catalogues.get("items", [])) < 15:
        raise ValueError("MGnify catalogue response is incomplete")
    marine = json.loads((cache / "mgnify-marine-v2.json").read_text())
    if marine.get("catalogue_id") != "marine-v2-0":
        raise ValueError("Unexpected MGnify Marine catalogue")

    version = (cache / "gtdb-r232-version.txt").read_text()
    release_notes = (cache / "gtdb-r232-release-notes.txt").read_text()
    if "v232" not in version or "Released Apr 15, 2026" not in version:
        raise ValueError("Unexpected GTDB R232 version")
    if "901,341 genomes" not in release_notes or "199,923 species clusters" not in release_notes:
        raise ValueError("Unexpected GTDB R232 release notes")

    with (cache / "gem-genome-metadata.tsv").open(encoding="utf-8") as handle:
        rows = sum(1 for _ in handle) - 1
    if rows != 52_515:
        raise ValueError(f"Unexpected GEM metadata rows: {rows}")
    with Image.open(cache / "gem-figure1-original.png") as image:
        if image.width < 1400 or image.height < 800:
            raise ValueError(f"GEM Figure 1 is too small: {image.size}")


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    seed_dir = args.seed_dir.resolve() if args.seed_dir else None

    for name, url in URLS.items():
        target = cache / name
        seed(target, seed_dir)
        if not target.is_file():
            if args.verify_only:
                raise FileNotFoundError(target)
            download(url, target)

    records: dict[str, dict[str, object]] = {}
    for name, url in URLS.items():
        path = cache / name
        observed = (path.stat().st_size, sha256(path))
        if name in LOCKS and observed != LOCKS[name]:
            raise RuntimeError(
                f"Locked resource drift for {name}: {observed} != {LOCKS[name]}"
            )
        records[name] = {
            "url": url,
            "bytes": observed[0],
            "sha256": observed[1],
        }

    validate_content(cache)
    with Image.open(cache / "gem-figure1-original.png") as image:
        records["gem-figure1-original.png"]["width"] = image.width
        records["gem-figure1-original.png"]["height"] = image.height

    manifest = {
        "article": 73,
        "snapshot_date": SNAPSHOT_DATE,
        "mgnify_api": "v2",
        "mgnify_study": "MGYS00000410",
        "ena_project": "PRJEB1787",
        "gtdb_release": "R11-RS232",
        "gtdb_release_date": "2026-04-15",
        "gem_paper_doi": "10.1038/s41587-020-0718-6",
        "resource_count": len(records),
        "resources": records,
    }
    (cache / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
