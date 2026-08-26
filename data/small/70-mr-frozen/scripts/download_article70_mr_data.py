#!/usr/bin/env python3
"""Download and checksum the public TwoSampleMR inputs used by Article 70."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path


TWOSAMPLEMR_VERSION = "0.7.9"
TWOSAMPLEMR_COMMIT = "3d119f20d6fc164b0c7f710f5590fee9580f2c7b"
RAW = f"https://raw.githubusercontent.com/MRCIEU/TwoSampleMR/{TWOSAMPLEMR_COMMIT}"
RESOURCES = {
    "twosamplemr-vig-perform-mr.RData": {
        "url": f"{RAW}/inst/extdata/vig_perform_mr.RData",
        "bytes": 900_980,
        "sha256": "7a13b142efeafc0fee9b80888d60d1db4b06f175afe02e3a3c45c0bc11d63502",
    },
    "hemani-mrbase-paper.pdf": {
        "url": "https://elifesciences.org/articles/34408.pdf",
        "bytes": 2_256_456,
        "sha256": "f4db594e534dc54417755c5368a782ff139d429e6b1f17b64af97c77d074876e",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--seed-rdata", type=Path)
    parser.add_argument("--seed-paper", type=Path)
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
        "twosamplemr-vig-perform-mr.RData": args.seed_rdata,
        "hemani-mrbase-paper.pdf": args.seed_paper,
    }
    manifest: dict[str, object] = {
        "article": 70,
        "analysis_example": "Body mass index to coronary heart disease",
        "design": "two-sample Mendelian randomization using summary associations",
        "package": "TwoSampleMR",
        "package_version": TWOSAMPLEMR_VERSION,
        "repository": "https://github.com/MRCIEU/TwoSampleMR",
        "repository_commit": TWOSAMPLEMR_COMMIT,
        "repository_tag": f"v{TWOSAMPLEMR_VERSION}",
        "exposure": {
            "opengwas_id": "ieu-a-2",
            "trait": "Body mass index",
            "study": "Locke et al. Nature 2015",
            "doi": "10.1038/nature14177",
            "maximum_sample_size": 339_224,
        },
        "outcome": {
            "opengwas_id": "ieu-a-7",
            "trait": "Coronary heart disease",
            "study": "Nikpay et al. Nature Genetics 2015",
            "doi": "10.1038/ng.3396",
            "cases": 60_801,
            "controls": 123_504,
            "sample_size": 184_305,
        },
        "anchor": {
            "study": "Hemani et al. eLife 2018",
            "doi": "10.7554/eLife.34408",
            "figure": "Figure 1",
        },
        "license_note": (
            "TwoSampleMR is MIT licensed; the eLife article is open access. "
            "Original GWAS data-use terms remain applicable."
        ),
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
