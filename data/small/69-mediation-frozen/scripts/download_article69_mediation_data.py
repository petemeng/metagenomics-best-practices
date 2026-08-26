#!/usr/bin/env python3
"""Download and checksum the public Spencer et al. inputs used by Article 69."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path


COMMIT = "a95b1a020b890dbe93a960dec946e197b63b7d15"
RAW = f"https://raw.githubusercontent.com/mda-primetr/Spencer_et_al_2021/{COMMIT}"
RESOURCES = {
    "spencer-human-wgs.xlsx": {
        "url": f"{RAW}/data/PD1_Wargo_Human_WGS_Relabund_and_metadata_light_filtering.xlsx",
        "bytes": 22_718_390,
        "sha256": "4fa937513f33a9fe2e1127554caf89c4287698a422123213428c73d0f2bb0968",
    },
    "spencer-data-dictionary.xlsx": {
        "url": f"{RAW}/docs/DataDictionary.xlsx",
        "bytes": 101_198,
        "sha256": "dc200443076f6fc252eed6b9239446712cfbf65aded581896a7498c0827de15c",
    },
    "spencer-paper.pdf": {
        "url": "https://www.statnews.com/wp-content/uploads/2022/01/science.aaz7015.pdf",
        "bytes": 3_071_441,
        "sha256": "414c9350d6fef5c1916b464c43b271992d14d6c814ae9817dfa8dd51b642c707",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--seed-workbook", type=Path)
    parser.add_argument("--seed-dictionary", type=Path)
    parser.add_argument("--seed-paper", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid(path: Path, expected: dict[str, object]) -> bool:
    return path.is_file() and path.stat().st_size == expected["bytes"] and sha256(path) == expected["sha256"]


def download(url: str, target: Path) -> None:
    partial = target.with_suffix(target.suffix + ".part")
    partial.unlink(missing_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "metagenomics-best-practices/1.0"})
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
        "spencer-human-wgs.xlsx": args.seed_workbook,
        "spencer-data-dictionary.xlsx": args.seed_dictionary,
        "spencer-paper.pdf": args.seed_paper,
    }
    manifest: dict[str, object] = {
        "article": 69,
        "study": "Spencer et al. Science 2021",
        "title": "Dietary fiber and probiotics influence the gut microbiome and melanoma immunotherapy response",
        "doi": "10.1126/science.aaz7015",
        "pmcid": "PMC8970537",
        "bioproject": "PRJNA770295",
        "repository": "https://github.com/mda-primetr/Spencer_et_al_2021",
        "repository_commit": COMMIT,
        "profile_type": "JAMS last-known-taxon relative abundance in parts per million",
        "license_note": "Public study data terms apply; the article panels are reproduced as a scholarly quotation for visual anchoring.",
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
