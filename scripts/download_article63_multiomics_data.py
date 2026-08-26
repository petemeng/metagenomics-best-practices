#!/usr/bin/env python3
"""Download and verify the paired Franzosa 2019 data used in Article 63."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.request
from pathlib import Path


COMMIT = "89a519d8c832008fbc6e650453e83e2f04858d02"
BASE_URL = (
    "https://raw.githubusercontent.com/borenstein-lab/"
    f"microbiome-metabolome-curated-data/{COMMIT}/"
    "data/processed_data/FRANZOSA_IBD_2019"
)
RESOURCES = (
    {
        "name": "genera.tsv",
        "bytes": 18101016,
        "sha256": "c4a541fe198a147beccd72d52fb2ebbf75a8cdf75cb3df75f823290971409d3f",
    },
    {
        "name": "mtb.tsv",
        "bytes": 12485385,
        "sha256": "528b5e5953bd3697dd1ecf551d810d536c0679bd922e3fa3a6956c1412c6288c",
    },
    {
        "name": "mtb.map.tsv",
        "bytes": 695256,
        "sha256": "0dcdcce04a4e9b2b9b1632a410959baa4802ea9e14fc7c44f63bc17f699e5c65",
    },
    {
        "name": "metadata.tsv",
        "bytes": 39838,
        "sha256": "f7396e3d6838b3b30f78b02bd568753757f84c956cd351966dbe654d50285376",
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("db/multiomics-cache/franzosa-89a519d8"),
    )
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify(path: Path, resource: dict[str, object]) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    observed_bytes = path.stat().st_size
    observed_sha256 = sha256(path)
    if observed_bytes != resource["bytes"]:
        raise RuntimeError(
            f"Size mismatch for {path.name}: {observed_bytes} != {resource['bytes']}"
        )
    if observed_sha256 != resource["sha256"]:
        raise RuntimeError(
            f"SHA-256 mismatch for {path.name}: {observed_sha256} != {resource['sha256']}"
        )
    return {
        "name": path.name,
        "bytes": observed_bytes,
        "sha256": observed_sha256,
        "url": f"{BASE_URL}/{path.name}",
    }


def download(url: str, target: Path, attempts: int = 10) -> None:
    partial = target.with_suffix(target.suffix + ".part")
    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "metagenomics-best-practices/1.0"},
            )
            with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
            partial.replace(target)
            return
        except Exception:
            if attempt == attempts:
                raise
            time.sleep(min(2 * attempt, 20))


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []
    for resource in RESOURCES:
        target = cache / str(resource["name"])
        if not target.is_file():
            if args.verify_only:
                raise FileNotFoundError(target)
            download(f"{BASE_URL}/{target.name}", target)
        manifest.append(verify(target, resource))

    payload = {
        "article": 63,
        "dataset": "FRANZOSA_IBD_2019",
        "repository": "borenstein-lab/microbiome-metabolome-curated-data",
        "commit": COMMIT,
        "paper_doi": "10.1038/s41564-018-0306-4",
        "resource_doi": "10.1038/s41522-022-00345-5",
        "license": "CC BY 4.0 for the curated resource article; source data terms retained",
        "resources": manifest,
    }
    (cache / "download-manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
