#!/usr/bin/env python3
"""Download and checksum the public HMP2 inputs used by Article 67."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path


RESOURCES = {
    "hmp2-metadata.csv": {
        "url": "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/metadata/hmp2_metadata_2018-08-20.csv",
        "bytes": 9_074_342,
        "sha256": "656b7bd97660ddb875548805e30bede31f2d1208293f7170d2d5755e33862ec9",
    },
    "taxonomic_profiles.tsv.gz": {
        "url": "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles.tsv.gz",
        "bytes": 666_634,
        "sha256": "5728531a6a1236371be6795b1da84ff5f6dd029035179d24d3c3be72d814e72c",
    },
    "lloyd-price-fig3-original.png": {
        "url": "https://media.springernature.com/lw1200/springer-static/image/art%3A10.1038%2Fs41586-019-1237-9/MediaObjects/41586_2019_1237_Fig3_HTML.png",
        "bytes": 349_963,
        "sha256": "b5ea6e10c6036945c97377705a3128c2d64164de0e84dc9d3874c2d242e5eed7",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--seed-metadata", type=Path)
    parser.add_argument("--seed-profiles", type=Path)
    parser.add_argument("--seed-figure", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid(path: Path, expected: dict[str, object]) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == expected["bytes"]
        and sha256(path) == expected["sha256"]
    )


def download(url: str, target: Path) -> None:
    partial = target.with_suffix(target.suffix + ".part")
    partial.unlink(missing_ok=True)
    request = urllib.request.Request(
        url, headers={"User-Agent": "metagenomics-best-practices/1.0"}
    )
    with urllib.request.urlopen(request, timeout=300) as response, partial.open("wb") as output:
        shutil.copyfileobj(response, output, length=8 * 1024 * 1024)
    partial.replace(target)


def materialize(target: Path, expected: dict[str, object], seed: Path | None) -> None:
    if valid(target, expected):
        return
    if seed is not None and valid(seed.resolve(), expected):
        shutil.copy2(seed.resolve(), target)
    else:
        download(str(expected["url"]), target)
    if not valid(target, expected):
        raise RuntimeError(f"Checksum or byte-count mismatch: {target}")


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    seeds = {
        "hmp2-metadata.csv": args.seed_metadata,
        "taxonomic_profiles.tsv.gz": args.seed_profiles,
        "lloyd-price-fig3-original.png": args.seed_figure,
    }
    manifest: dict[str, object] = {
        "article": 67,
        "study": "Lloyd-Price et al. Nature 2019",
        "title": "Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases",
        "doi": "10.1038/s41586-019-1237-9",
        "pmcid": "PMC6650278",
        "cohort": "Integrative Human Microbiome Project / IBDMDB (HMP2)",
        "profile_release": "2018-05-04",
        "metadata_release": "2018-08-20",
        "license_note": "HMP2 public data terms apply; the article figure is used for scholarly quotation and visual anchoring.",
        "resources": {},
    }
    for name, expected in RESOURCES.items():
        target = cache / name
        materialize(target, expected, seeds[name])
        manifest["resources"][name] = {**expected, "downloaded_path": name}
        print(f"verified\t{name}\t{target.stat().st_size}\t{sha256(target)}")
    (cache / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
