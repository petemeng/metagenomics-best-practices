#!/usr/bin/env python3
"""Download and checksum-lock the public Article 71 SEM inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import time
import urllib.request
from pathlib import Path


COMMIT = "89a519d8c832008fbc6e650453e83e2f04858d02"
BASE_URL = (
    "https://raw.githubusercontent.com/borenstein-lab/"
    f"microbiome-metabolome-curated-data/{COMMIT}/"
    "data/processed_data/FRANZOSA_IBD_2019"
)
ANCHOR_URL = (
    "https://media.springernature.com/full/springer-static/image/"
    "art%3A10.1038%2Fs41564-018-0306-4/MediaObjects/"
    "41564_2018_306_Fig1_HTML.png"
)
RESOURCES = (
    {
        "name": "genera.tsv",
        "bytes": 18_101_016,
        "sha256": "c4a541fe198a147beccd72d52fb2ebbf75a8cdf75cb3df75f823290971409d3f",
        "url": f"{BASE_URL}/genera.tsv",
    },
    {
        "name": "metadata.tsv",
        "bytes": 39_838,
        "sha256": "f7396e3d6838b3b30f78b02bd568753757f84c956cd351966dbe654d50285376",
        "url": f"{BASE_URL}/metadata.tsv",
    },
    {
        "name": "franzosa-fig1-original.png",
        "bytes": 170_306,
        "sha256": "7b81b865ae65659ad476d6f5210a3bc383b4eadad50d2ad7793a0b99df2450eb",
        "url": ANCHOR_URL,
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cache-dir", type=Path, default=Path("db/sem/article71")
    )
    parser.add_argument(
        "--seed-resource-dir",
        type=Path,
        help="Optional directory containing verified genera.tsv and metadata.tsv",
    )
    parser.add_argument(
        "--seed-anchor",
        type=Path,
        help="Optional verified copy of the official Figure 1 PNG",
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
    observed_size = path.stat().st_size
    observed_sha = sha256(path)
    if observed_size != resource["bytes"]:
        raise RuntimeError(
            f"Size mismatch for {path.name}: {observed_size} != {resource['bytes']}"
        )
    if observed_sha != resource["sha256"]:
        raise RuntimeError(
            f"SHA-256 mismatch for {path.name}: {observed_sha} != {resource['sha256']}"
        )
    return {
        "name": path.name,
        "bytes": observed_size,
        "sha256": observed_sha,
        "url": resource["url"],
    }


def download(url: str, target: Path, attempts: int = 10) -> None:
    partial = target.with_suffix(target.suffix + ".part")
    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "metagenomics-best-practices/1.0"},
            )
            with urllib.request.urlopen(
                request, timeout=180
            ) as response, partial.open("wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
            partial.replace(target)
            return
        except Exception:
            partial.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            time.sleep(min(2 * attempt, 20))


def seed_if_available(
    target: Path,
    resource: dict[str, object],
    resource_dir: Path | None,
    anchor: Path | None,
) -> None:
    if target.is_file():
        return
    if target.name == "franzosa-fig1-original.png" and anchor is not None:
        candidate = anchor
    elif resource_dir is not None:
        candidate = resource_dir / target.name
    else:
        return
    if candidate.is_file():
        verify(candidate, resource)
        shutil.copy2(candidate, target)


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    seed_dir = (
        args.seed_resource_dir.resolve()
        if args.seed_resource_dir is not None
        else None
    )
    seed_anchor = (
        args.seed_anchor.resolve() if args.seed_anchor is not None else None
    )
    records: list[dict[str, object]] = []
    for resource in RESOURCES:
        target = cache / str(resource["name"])
        seed_if_available(target, resource, seed_dir, seed_anchor)
        if not target.is_file():
            if args.verify_only:
                raise FileNotFoundError(target)
            download(str(resource["url"]), target)
        records.append(verify(target, resource))

    manifest = {
        "article": 71,
        "dataset": "FRANZOSA_IBD_2019",
        "repository": "borenstein-lab/microbiome-metabolome-curated-data",
        "repository_commit": COMMIT,
        "paper_doi": "10.1038/s41564-018-0306-4",
        "resource_doi": "10.1038/s41522-022-00345-5",
        "anchor_figure": "Franzosa et al. 2019 Figure 1",
        "license": (
            "CC BY 4.0 for the curated resource article; source data and "
            "Nature figure terms remain applicable"
        ),
        "resources": {record["name"]: record for record in records},
    }
    (cache / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
