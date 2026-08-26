#!/usr/bin/env python3
"""Download and verify the public hot-spring tables used in Article 62."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path


RESOURCES = (
    {
        "id": 61153429,
        "name": "Data_01 - Sample_metadata.tsv",
        "local": "sample-metadata.tsv",
        "bytes": 133160,
        "md5": "402095b04e2c5518cbec462f41528d4f",
    },
    {
        "id": 61153438,
        "name": "Data_05 - KO_proportions_in_metagenomes.tsv.gz",
        "local": "ko-proportions-in-metagenomes.tsv.gz",
        "bytes": 15712433,
        "md5": "7eeca888e63b137b8dbf69e3e1bf4c1e",
    },
    {
        "id": 62623690,
        "name": "Data_07 - MAG_metadata.tsv",
        "local": "mag-metadata.tsv",
        "bytes": 3750601,
        "md5": "80dc62e4b3b0bf8c70b42b446d2979aa",
    },
    {
        "id": 61153444,
        "name": "Data_08 - MAG_abundances_per_sample.biom",
        "local": "mag-abundances-per-sample.biom",
        "bytes": 2017689,
        "md5": "0c1762f31a5c4473dec78571a7b74287",
    },
    {
        "id": 61153468,
        "name": "Data_15 - KOs_in_MAGs.tsv.gz",
        "local": "kos-in-mags.tsv.gz",
        "bytes": 684507,
        "md5": "35c3732f5b5e3c8c406fe9402f0a789e",
    },
)

GITHUB_RESOURCES = (
    {
        "name": "DiTing Pathway_formulas.txt v0.3",
        "local": "diting-pathway-formulas-v0.3.txt",
        "bytes": 70799,
        "git_blob_sha1": "68fa7a0e1f1462e18a6d0d4978643729c2ffc3e3",
        "url": (
            "https://raw.githubusercontent.com/xuechunxu/DiTing/"
            "53e1d3edb84be08b7aacb79ac588be671250b477/Pathway_formulas.txt"
        ),
        "commit": "53e1d3edb84be08b7aacb79ac588be671250b477",
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, default=Path("db/element-cycle-cache"))
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def digest(path: Path, algorithm: str) -> str:
    checksum = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def verify(path: Path, resource: dict[str, object]) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    observed_bytes = path.stat().st_size
    observed_md5 = digest(path, "md5")
    if observed_bytes != resource["bytes"]:
        raise RuntimeError(
            f"Size mismatch for {path.name}: {observed_bytes} != {resource['bytes']}"
        )
    if observed_md5 != resource["md5"]:
        raise RuntimeError(
            f"MD5 mismatch for {path.name}: {observed_md5} != {resource['md5']}"
        )
    return {
        "figshare_file_id": resource["id"],
        "source_name": resource["name"],
        "local_file": path.name,
        "bytes": observed_bytes,
        "md5": observed_md5,
        "sha256": digest(path, "sha256"),
        "url": f"https://ndownloader.figshare.com/files/{resource['id']}",
    }


def git_blob_sha1(path: Path) -> str:
    payload = path.read_bytes()
    header = f"blob {len(payload)}\0".encode()
    return hashlib.sha1(header + payload).hexdigest()


def verify_github(path: Path, resource: dict[str, object]) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    observed_bytes = path.stat().st_size
    observed_blob = git_blob_sha1(path)
    if observed_bytes != resource["bytes"]:
        raise RuntimeError(
            f"Size mismatch for {path.name}: {observed_bytes} != {resource['bytes']}"
        )
    if observed_blob != resource["git_blob_sha1"]:
        raise RuntimeError(
            f"Git blob mismatch for {path.name}: {observed_blob} != {resource['git_blob_sha1']}"
        )
    return {
        "source_name": resource["name"],
        "local_file": path.name,
        "bytes": observed_bytes,
        "git_blob_sha1": observed_blob,
        "sha256": digest(path, "sha256"),
        "url": resource["url"],
        "commit": resource["commit"],
    }


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []
    for resource in RESOURCES:
        target = cache / str(resource["local"])
        if not target.exists():
            if args.verify_only:
                raise FileNotFoundError(target)
            partial = target.with_suffix(target.suffix + ".part")
            request = urllib.request.Request(
                f"https://ndownloader.figshare.com/files/{resource['id']}",
                headers={"User-Agent": "metagenomics-best-practices/1.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
            partial.replace(target)
        manifest.append(verify(target, resource))

    for resource in GITHUB_RESOURCES:
        target = cache / str(resource["local"])
        if not target.exists():
            if args.verify_only:
                raise FileNotFoundError(target)
            partial = target.with_suffix(target.suffix + ".part")
            request = urllib.request.Request(
                str(resource["url"]),
                headers={"User-Agent": "metagenomics-best-practices/1.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
            partial.replace(target)
        manifest.append(verify_github(target, resource))

    payload = {
        "article": 62,
        "dataset_doi": "10.6084/m9.figshare.30284068.v2",
        "paper_doi": "10.1038/s41597-026-07139-w",
        "figshare_article_id": 30284068,
        "figshare_version": 2,
        "resources": manifest,
    }
    (cache / "download-manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
