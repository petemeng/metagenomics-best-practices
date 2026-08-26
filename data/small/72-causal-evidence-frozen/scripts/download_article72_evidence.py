#!/usr/bin/env python3
"""Download canonical DOI metadata and the Article 72 anchor figure."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import time
import urllib.parse
import urllib.request
from pathlib import Path

from PIL import Image


PUBLICATIONS = (
    ("vannood2013fmt", "10.1056/NEJMoa1205037"),
    ("feuerstadt2022ser109", "10.1056/NEJMoa2106516"),
    ("smillie2018engraftment", "10.1016/j.chom.2018.01.003"),
    ("li2016fmtstrains", "10.1126/science.aad8852"),
    ("buffie2015cscindens", "10.1038/nature13828"),
    ("theriot2014metabolome", "10.1038/ncomms4114"),
    ("wirbel2019crc", "10.1038/s41591-019-0406-6"),
    ("kostic2013fusobacterium", "10.1016/j.chom.2013.07.007"),
    ("rubinstein2013fada", "10.1016/j.chom.2013.07.012"),
    ("bullman2017fusobacterium", "10.1126/science.aal5240"),
    ("neville2018commensalkoch", "10.1016/j.mib.2017.10.001"),
)

ANCHOR_URLS = (
    "https://media.springernature.com/full/springer-static/image/"
    "art%3A10.1038%2Fnature13828/MediaObjects/"
    "41586_2015_Article_BFnature13828_Fig4_HTML.jpg",
    "https://media.springernature.com/lw1200/springer-static/image/"
    "art%3A10.1038%2Fnature13828/MediaObjects/"
    "41586_2015_Article_BFnature13828_Fig4_HTML.jpg",
    "https://media.springernature.com/m685/springer-static/image/"
    "art%3A10.1038%2Fnature13828/MediaObjects/"
    "41586_2015_Article_BFnature13828_Fig4_HTML.jpg",
)

# Filled after the first verified acquisition; these locks prevent silent drift.
METADATA_BYTES = 23_575
METADATA_SHA256 = "c48e3ac9ef4ecf0da7c1f92e8dbb59f93b42f3841a79960284de5a9cc2825a53"
ANCHOR_BYTES = 143_655
ANCHOR_SHA256 = "aaa42762d0f3c9287a266a2681774a90b908b64d9614dbf8ec1d1b93588792c9"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, default=Path("db/evidence/article72"))
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
                        "(mailto:example@example.org)"
                    )
                },
            )
            with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
            partial.replace(target)
            return
        except Exception:
            partial.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            time.sleep(min(attempt * 2, 16))


def canonical_crossref(message: dict[str, object], key: str, doi: str) -> dict[str, object]:
    observed_doi = str(message.get("DOI", "")).lower()
    if observed_doi != doi.lower():
        raise ValueError(f"Crossref DOI mismatch: {observed_doi} != {doi}")
    published = message.get("published-print") or message.get("published-online") or message.get("issued")
    date_parts = published["date-parts"][0]
    authors = []
    for author in message.get("author", []):
        authors.append(
            {
                "family": author.get("family", ""),
                "given": author.get("given", ""),
                "orcid": author.get("ORCID", ""),
            }
        )
    return {
        "citation_key": key,
        "doi": observed_doi,
        "title": (message.get("title") or [""])[0],
        "journal": (message.get("container-title") or [""])[0],
        "year": int(date_parts[0]),
        "volume": str(message.get("volume", "")),
        "issue": str(message.get("issue", "")),
        "pages": str(message.get("page", "")),
        "publisher": str(message.get("publisher", "")),
        "type": str(message.get("type", "")),
        "authors": authors,
        "source_url": f"https://doi.org/{observed_doi}",
    }


def acquire_metadata(target: Path) -> None:
    records = []
    for key, doi in PUBLICATIONS:
        encoded = urllib.parse.quote(doi, safe="")
        url = f"https://api.crossref.org/works/{encoded}"
        request = urllib.request.Request(
            url,
            headers={
                "User-Agent": (
                    "metagenomics-best-practices/1.0 "
                    "(mailto:example@example.org)"
                )
            },
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = json.load(response)
        records.append(canonical_crossref(payload["message"], key, doi))
    target.write_text(
        json.dumps(records, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def verify_locked(path: Path, expected_bytes: int | None, expected_sha: str | None) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    observed_bytes = path.stat().st_size
    observed_sha = sha256(path)
    if expected_bytes is not None and observed_bytes != expected_bytes:
        raise RuntimeError(f"Size mismatch for {path.name}: {observed_bytes} != {expected_bytes}")
    if expected_sha is not None and observed_sha != expected_sha:
        raise RuntimeError(f"SHA-256 mismatch for {path.name}: {observed_sha} != {expected_sha}")
    return {"name": path.name, "bytes": observed_bytes, "sha256": observed_sha}


def acquire_anchor(target: Path) -> str:
    errors = []
    for url in ANCHOR_URLS:
        try:
            download(url, target)
            with Image.open(target) as image:
                if image.width < 900 or image.height < 500:
                    raise ValueError(f"Anchor is too small: {image.size}")
                image.verify()
            return url
        except Exception as error:
            target.unlink(missing_ok=True)
            errors.append(f"{url}: {error}")
    raise RuntimeError("Could not acquire anchor figure:\n" + "\n".join(errors))


def seed(target: Path, seed_dir: Path | None) -> None:
    if target.is_file() or seed_dir is None:
        return
    candidate = seed_dir / target.name
    if candidate.is_file():
        shutil.copy2(candidate, target)


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    seed_dir = args.seed_dir.resolve() if args.seed_dir is not None else None

    metadata = cache / "publication-metadata.json"
    anchor = cache / "buffie-figure4-original.jpg"
    seed(metadata, seed_dir)
    seed(anchor, seed_dir)

    if not metadata.is_file():
        if args.verify_only:
            raise FileNotFoundError(metadata)
        acquire_metadata(metadata)
    metadata_record = verify_locked(metadata, METADATA_BYTES, METADATA_SHA256)

    anchor_url = ANCHOR_URLS[0]
    if not anchor.is_file():
        if args.verify_only:
            raise FileNotFoundError(anchor)
        anchor_url = acquire_anchor(anchor)
    anchor_record = verify_locked(anchor, ANCHOR_BYTES, ANCHOR_SHA256)
    with Image.open(anchor) as image:
        anchor_record["width"] = image.width
        anchor_record["height"] = image.height

    manifest = {
        "article": 72,
        "metadata_source": "Crossref REST API",
        "publication_count": len(PUBLICATIONS),
        "doi_set": [doi.lower() for _, doi in PUBLICATIONS],
        "anchor": {
            **anchor_record,
            "url": anchor_url,
            "paper_doi": "10.1038/nature13828",
            "figure": "Figure 4",
        },
        "resources": {
            "publication-metadata.json": {
                **metadata_record,
                "url_template": "https://api.crossref.org/works/{doi}",
            }
        },
    }
    (cache / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
