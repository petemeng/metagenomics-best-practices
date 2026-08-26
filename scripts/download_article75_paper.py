#!/usr/bin/env python3
"""Download or verify the public paper evidence used by Article 75.

The Europe PMC supplementary endpoint creates a transport ZIP whose timestamps may
change between requests.  Reproducibility is therefore defined by the SHA-256 of
the selected members, not by the byte identity of the outer ZIP container.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


ARTICLE = 75
SNAPSHOT_DATE = "2026-08-23"
PMCID = "PMC7984229"
DOI = "10.1038/s41591-019-0406-6"
TITLE = (
    "Meta-analysis of fecal metagenomes reveals global microbial signatures "
    "that are specific for colorectal cancer"
)
XML_RECORD = {
    "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7984229/fullTextXML",
    "bytes": 168853,
    "sha256": "b6363f2dcf652a352546e8550b43673ffc6518047e35bd722d1f60653e5793a8",
}
ZIP_URL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7984229/"
    "supplementaryFiles?includeInlineImage=true"
)
FIGURE_MEMBERS = {
    "emss-81948-f001.jpg": {
        "bytes": 430124,
        "sha256": "af9131933d8cdee04b209cae966ff7bda430dd8d7fcc6b7d832cbbae3222c357",
        "md5": "b9488928e0770401e85eb27faf123496",
    },
    "emss-81948-f002.jpg": {
        "bytes": 431025,
        "sha256": "8a36654583df120f956bfbcefb2878f0b8261bf230bbdc192b3845d7ea536802",
        "md5": "55b02c83605054304cbf36faa3c8ce09",
    },
    "emss-81948-f003.jpg": {
        "bytes": 299113,
        "sha256": "8a2ede4a5346955e9de94ab3d483c0b947804d537764687a87e037edc2782a79",
        "md5": "8224d93016240da48859bfbb2e0a7f51",
    },
    "emss-81948-f004.jpg": {
        "bytes": 395299,
        "sha256": "140c46c99bae74d8a658c3e2c07601651d4422d7c80f97fd572e06878ccde085",
        "md5": "78de46e85c0dd5e339bb1508d63295b7",
    },
    "emss-81948-f005.jpg": {
        "bytes": 209051,
        "sha256": "8bed178f11fa58804d2dcbe61767b4fe17755b4c7cc593e9f710ad096c0606be",
        "md5": "6a607e6dd698a13a72bf5a3190e53532",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_xml(path: Path) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == XML_RECORD["bytes"]
        and sha256_file(path) == XML_RECORD["sha256"]
    )


def valid_zip(path: Path) -> bool:
    if not path.is_file():
        return False
    try:
        with zipfile.ZipFile(path) as archive:
            for name, record in FIGURE_MEMBERS.items():
                payload = archive.read(name)
                if len(payload) != record["bytes"]:
                    return False
                if sha256_bytes(payload) != record["sha256"]:
                    return False
                if hashlib.md5(payload).hexdigest() != record["md5"]:
                    return False
    except (KeyError, OSError, zipfile.BadZipFile):
        return False
    return True


def download(url: str, destination: Path) -> None:
    headers = {"User-Agent": "metagenomics-best-practices/article-75"}
    request = urllib.request.Request(url, headers=headers)
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


def validate_semantics(xml_path: Path) -> None:
    root = ET.parse(xml_path).getroot()
    doi = root.findtext('.//article-id[@pub-id-type="doi"]')
    title_node = root.find(".//article-title")
    title = "" if title_node is None else " ".join("".join(title_node.itertext()).split())
    if doi != DOI or title != TITLE:
        raise ValueError(f"Unexpected paper identity: DOI={doi!r}; title={title!r}")
    figures: dict[str, str] = {}
    xlink = "{http://www.w3.org/1999/xlink}href"
    for figure in root.findall(".//fig"):
        label = figure.find("label")
        caption = figure.find("caption")
        graphic = figure.find("graphic")
        if label is None or caption is None or graphic is None:
            continue
        figure_label = " ".join("".join(label.itertext()).split())
        if figure_label in {f"Figure {number}" for number in range(1, 6)}:
            figures[figure_label] = graphic.attrib.get(xlink, "")
    expected = {
        f"Figure {number}": f"emss-81948-f{number:03d}.jpg"
        for number in range(1, 6)
    }
    if figures != expected:
        raise ValueError(f"Unexpected main-figure mapping: {figures}")


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    xml_path = output / "PMC7984229-fulltext.xml"
    zip_path = output / "PMC7984229-supplementary.zip"

    if not valid_xml(xml_path):
        if args.verify_only:
            raise RuntimeError(f"Missing or checksum-mismatched source: {xml_path}")
        download(str(XML_RECORD["url"]), xml_path)
    if not valid_xml(xml_path):
        raise RuntimeError("Downloaded Europe PMC XML failed its checksum lock")

    if not valid_zip(zip_path):
        if args.verify_only:
            raise RuntimeError(f"Missing or member-mismatched source: {zip_path}")
        download(ZIP_URL, zip_path)
    if not valid_zip(zip_path):
        raise RuntimeError("Downloaded Europe PMC ZIP failed selected-member locks")

    validate_semantics(xml_path)
    manifest = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "paper": {"pmcid": PMCID, "doi": DOI, "title": TITLE},
        "full_text_xml": XML_RECORD,
        "supplementary_endpoint": ZIP_URL,
        "transport_note": (
            "The ZIP wrapper is not byte-locked because server-generated member "
            "timestamps may change; the five selected JPEG members are byte-locked."
        ),
        "selected_members": FIGURE_MEMBERS,
    }
    (output / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"verified\t{xml_path}\t{XML_RECORD['sha256']}")
    print(f"verified\t{zip_path}\t{len(FIGURE_MEMBERS)} selected figure members")
    print(f"manifest\t{output / 'download-manifest.json'}")


if __name__ == "__main__":
    main()
