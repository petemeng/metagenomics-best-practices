#!/usr/bin/env python3
"""Download and checksum the open Article 66 paper and supplement."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.request
import zipfile
from pathlib import Path


RESOURCES = {
    "paper.xml": {
        "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11406951/fullTextXML",
        "bytes": 122_202,
        "sha256": "6b3fcd4db9bbb0a3c6dd7409db3a70ee7fa5154d04a06d2112626281a34eccec",
    },
    "supplementary-files.zip": {
        "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11406951/supplementaryFiles",
        "bytes": 955_771,
        "sha256": "49c674c1bad788ea15e3d48ed8b2909b0d67c22e318276103abe4b0d307f28ae",
    },
}

MEMBERS = {
    "supplementary-tables.xlsx": {
        "member": "msystems.00746-24-s0002.xlsx",
        "bytes": 101_735,
        "sha256": "c740876ed3dce897ac9105a5434a88792afa8bc6bc072ae588e5cc1173f6d416",
    },
    "majzoub-fig2-original.jpg": {
        "member": "msystems.00746-24.f002.jpg",
        "bytes": 41_669,
        "sha256": "e96f19623d9337a6cb41c1ea0294b13ea69426df1166b88cf19fa55c661a284b",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--seed-xml", type=Path)
    parser.add_argument("--seed-zip", type=Path)
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


def materialize(
    name: str,
    expected: dict[str, object],
    target: Path,
    seed: Path | None,
) -> None:
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
    seeds = {"paper.xml": args.seed_xml, "supplementary-files.zip": args.seed_zip}

    manifest: dict[str, object] = {
        "article": 66,
        "study": "Majzoub et al. mSystems 2024",
        "title": "Refining microbial community metabolic models derived from metagenomics using reference-based taxonomic profiling",
        "doi": "10.1128/msystems.00746-24",
        "pmcid": "PMC11406951",
        "human_metagenomics_accession": "PRJEB50699",
        "license": "CC BY 4.0 for the article; supplement redistributed by the publisher under the authors' licence",
        "retrieval_service": "Europe PMC",
        "resources": {},
        "extracted_members": {},
    }
    for name, expected in RESOURCES.items():
        target = cache / name
        materialize(name, expected, target, seeds[name])
        manifest["resources"][name] = {**expected, "downloaded_path": name}
        print(f"verified\t{name}\t{target.stat().st_size}\t{sha256(target)}")

    archive = cache / "supplementary-files.zip"
    with zipfile.ZipFile(archive) as handle:
        names = set(handle.namelist())
        for output_name, expected in MEMBERS.items():
            member = str(expected["member"])
            if member not in names:
                raise RuntimeError(f"Missing archive member: {member}")
            target = cache / output_name
            with handle.open(member) as source, target.open("wb") as output:
                shutil.copyfileobj(source, output)
            if not valid(target, expected):
                raise RuntimeError(f"Extracted member mismatch: {target}")
            manifest["extracted_members"][output_name] = {
                **expected,
                "extracted_path": output_name,
            }
            print(f"verified\t{output_name}\t{target.stat().st_size}\t{sha256(target)}")

    (cache / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
