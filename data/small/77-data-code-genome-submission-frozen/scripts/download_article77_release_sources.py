#!/usr/bin/env python3
"""Download and verify the version-locked primary source for Article 77."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from PIL import Image


ARTICLE = 77
SNAPSHOT_DATE = "2026-08-23"
XML = {
    "filename": "PMC5737865-fulltext.xml",
    "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5737865/fullTextXML",
    "bytes": 118079,
    "sha256": "8952fa2f6bc3c8960e67a71be152a16d6700aa81405e8a44c984c6cb4b7f0d90",
}
FIGURE = {
    "filename": "gix047fig2.jpg",
    "archive_url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5737865/supplementaryFiles",
    "archive_member": "gix047fig2.jpg",
    "bytes": 69518,
    "sha256": "1a7ebe7ca72de90d9de9a6f7d5daec8b641a73fece7b3458d40a46fc45237477",
    "width": 767,
    "height": 508,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--seed-dir",
        type=Path,
        help="Optional verified local cache; network remains the fallback.",
    )
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def valid(path: Path, record: dict[str, object]) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == record["bytes"]
        and sha256(path) == record["sha256"]
    )


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(
        url, headers={"User-Agent": "metagenomics-best-practices/article-77"}
    )
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                while block := response.read(8 * 1024 * 1024):
                    handle.write(block)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(destination)


def seed_candidates(seed_dir: Path | None, filename: str) -> list[Path]:
    if seed_dir is None:
        return []
    aliases = {
        XML["filename"]: ("article77-lifecycle.xml", XML["filename"]),
        FIGURE["filename"]: ("article77-lifecycle-fig2.jpg", FIGURE["filename"]),
    }
    return [seed_dir / name for name in aliases[filename]]


def use_seed(seed_dir: Path | None, destination: Path, record: dict[str, object]) -> bool:
    for candidate in seed_candidates(seed_dir, destination.name):
        if valid(candidate, record):
            shutil.copy2(candidate, destination)
            return True
    return False


def obtain_xml(output: Path, seed_dir: Path | None, verify_only: bool) -> None:
    destination = output / str(XML["filename"])
    if valid(destination, XML):
        return
    if verify_only:
        raise RuntimeError(f"Missing or checksum-mismatched source: {destination}")
    if not use_seed(seed_dir, destination, XML):
        download(str(XML["url"]), destination)
    if not valid(destination, XML):
        raise RuntimeError(f"XML failed byte lock: {destination}")


def obtain_figure(output: Path, seed_dir: Path | None, verify_only: bool) -> None:
    destination = output / str(FIGURE["filename"])
    if valid(destination, FIGURE):
        return
    if verify_only:
        raise RuntimeError(f"Missing or checksum-mismatched source: {destination}")
    if use_seed(seed_dir, destination, FIGURE):
        return
    archive = output / "europepmc-supplementaryFiles.tmp.zip"
    download(str(FIGURE["archive_url"]), archive)
    try:
        with zipfile.ZipFile(archive) as bundle:
            with bundle.open(str(FIGURE["archive_member"])) as source:
                with destination.open("wb") as target:
                    shutil.copyfileobj(source, target)
    finally:
        archive.unlink(missing_ok=True)
    if not valid(destination, FIGURE):
        raise RuntimeError(f"Selected archive member failed byte lock: {destination}")


def semantic_checks(output: Path) -> None:
    root = ET.parse(output / str(XML["filename"])).getroot()
    doi = root.findtext('.//article-id[@pub-id-type="doi"]')
    if doi != "10.1093/gigascience/gix047":
        raise ValueError(f"Unexpected source DOI: {doi!r}")
    license_text = " ".join("".join(root.find(".//permissions").itertext()).split())
    if "Creative Commons Attribution License" not in license_text:
        raise ValueError("The expected CC BY statement is absent")
    hrefs = {
        value
        for node in root.findall(".//graphic")
        for key, value in node.attrib.items()
        if key.endswith("href")
    }
    if FIGURE["archive_member"] not in hrefs:
        raise ValueError("Figure 2 is not linked by the full-text XML")
    with Image.open(output / str(FIGURE["filename"])) as image:
        expected = (FIGURE["width"], FIGURE["height"])
        if image.size != expected:
            raise ValueError(f"Unexpected Figure 2 dimensions: {image.size}")


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    seed_dir = args.seed_dir.resolve() if args.seed_dir else None
    obtain_xml(output, seed_dir, args.verify_only)
    obtain_figure(output, seed_dir, args.verify_only)
    semantic_checks(output)
    manifest = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "source": {
            "identity": "ten Hoopen et al. 2017; DOI 10.1093/gigascience/gix047",
            "pmcid": "PMC5737865",
            "license": "CC BY 4.0",
            "xml": XML,
            "selected_figure_member": FIGURE,
        },
        "archive_policy": (
            "The Europe PMC ZIP wrapper is transport-only and is not locked because "
            "its timestamps may change; the selected Figure 2 member is byte-locked."
        ),
    }
    (output / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"verified\t{XML['filename']}\t{XML['sha256']}")
    print(f"verified\t{FIGURE['filename']}\t{FIGURE['sha256']}")


if __name__ == "__main__":
    main()
