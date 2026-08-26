#!/usr/bin/env python3
"""Download and checksum the official IBDMDB inputs used in Article 64."""

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
    "mgx-pathabundance-rela.tsv.gz": {
        "url": "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/pathabundance_relab.tsv.gz",
        "bytes": 9_037_246,
        "sha256": "e840fb86ed8049bc697a20a3904da11d29878ed02f43f50764a28c93d2111216",
    },
    "mtx-pathabundance-rela.tsv.gz": {
        "url": "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MTX/2017-12-14/pathabundance_relab.tsv.gz",
        "bytes": 1_715_915,
        "sha256": "ee2a1afb69b66bbdac014b48db6a692dc09146d539e4a93d3fea0c2d9903ac08",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
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
    temporary = target.with_suffix(target.suffix + ".part")
    temporary.unlink(missing_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "metagenomics-best-practices/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as output:
        shutil.copyfileobj(response, output, length=8 * 1024 * 1024)
    temporary.replace(target)


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, object] = {
        "article": 64,
        "dataset": "IBDMDB_HMP2",
        "study": "Lloyd-Price et al. Nature 2019",
        "doi": "10.1038/s41586-019-1237-9",
        "resources": {},
    }
    for name, expected in RESOURCES.items():
        target = cache / name
        if not valid(target, expected):
            download(str(expected["url"]), target)
        if not valid(target, expected):
            raise RuntimeError(f"Checksum or byte-count mismatch: {target}")
        manifest["resources"][name] = {
            **expected,
            "downloaded_path": name,
        }
        print(f"verified\t{name}\t{target.stat().st_size}\t{sha256(target)}")
    (cache / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
