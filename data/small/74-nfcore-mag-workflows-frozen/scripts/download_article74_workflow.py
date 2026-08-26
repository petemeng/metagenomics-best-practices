#!/usr/bin/env python3
"""Download or verify the immutable official inputs used by Article 74."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
import urllib.request
from pathlib import Path


ARTICLE = 74
SNAPSHOT_DATE = "2026-08-24"

SOURCES = {
    "mag-release.json": {
        "url": "https://api.github.com/repos/nf-core/mag/releases/tags/5.5.0",
        "bytes": 12421,
        "sha256": "420560790437a4efbf763a7cc4b462fa9e2ee4fd26ff72259d76c2fa0675eb5e",
    },
    "mag-tag.json": {
        "url": "https://api.github.com/repos/nf-core/mag/git/ref/tags/5.5.0",
        "bytes": 368,
        "sha256": "8b7d951bac555ee2d813c6ac0e660856b9bcce23b1c5e19f39492360f1fee3c2",
    },
    "mag-5.5.0.tar.gz": {
        "url": "https://github.com/nf-core/mag/archive/refs/tags/5.5.0.tar.gz",
        "bytes": 4466883,
        "sha256": "871528217457f88677e1dcd5a3822bed9d1234c3f07366652cb5f49b7a2700b3",
    },
    "funcscan-release.json": {
        "url": "https://api.github.com/repos/nf-core/funcscan/releases/tags/4.0.0",
        "bytes": 5288,
        "sha256": "10f9e3f8ee140af4eb714bb10e413821c730fc8730a1c1edc4c990066e8ff8ce",
    },
    "funcscan-tag.json": {
        "url": "https://api.github.com/repos/nf-core/funcscan/git/ref/tags/4.0.0",
        "bytes": 369,
        "sha256": "f3c0cc06aab4a68fd190691213fff6d55806361831553667f2309ff42d18ed76",
    },
    "funcscan-4.0.0.tar.gz": {
        "url": "https://github.com/nf-core/funcscan/archive/refs/tags/4.0.0.tar.gz",
        "bytes": 2768719,
        "sha256": "1d387cdd7265529ac6cebceda890334578a35d71755378ccdbe32fd9b5c8f224",
    },
    "nextflow-release.json": {
        "url": "https://api.github.com/repos/nextflow-io/nextflow/releases/tags/v26.04.0",
        "bytes": 25059,
        "sha256": "954c127fe696886cedb6c389f50a53c2927b8cf1af00505b9977d3d01f9abf0b",
    },
    "nextflow": {
        "url": "https://github.com/nextflow-io/nextflow/releases/download/v26.04.0/nextflow",
        "bytes": 17246,
        "sha256": "2f0e68fa22df782bbebad4a964138756f2cf19a5544c07af2913a5f730646e44",
    },
    "nextflow-26.04.0-dist": {
        "url": "https://github.com/nextflow-io/nextflow/releases/download/v26.04.0/nextflow-26.04.0-dist",
        "bytes": 41712507,
        "sha256": "5e2b4a354b4d7634d7211b71417d61606878fb49e9b224b50ded6e2c69114870",
    },
    "checkm2-record.json": {
        "url": "https://zenodo.org/api/records/14897628",
        "bytes": 3707,
        "sha256": "3ffa00d3993d656db8ca6ebe51e0638565283e327ab4a129117c954e81546dfd",
    },
    "gunc-progenomes2.1.dmnd.gz.md5": {
        "url": "https://black.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz.md5",
        "bytes": 56,
        "sha256": "82e2160c120dc3015a87cb019d5e2217ae2702396ebd9b2034076d3833dedecc",
    },
    "gunc-progenomes2.1.dmnd.md5": {
        "url": "https://black.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.md5",
        "bytes": 53,
        "sha256": "48919518491557dc1903764f234b3076dd287f59d55216075fb7c5e387d82a41",
    },
    "gunc-database-v1.1.0.py": {
        "url": "https://raw.githubusercontent.com/grp-bork/gunc/v1.1.0/gunc/gunc_database.py",
        "bytes": 6382,
        "sha256": "f3514cfd0ad7662311108b6c1dad0d7c01089e7db6ebc20a9e1fb9fb92f11778",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid(path: Path, record: dict[str, object]) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == record["bytes"]
        and sha256(path) == record["sha256"]
    )


def download(url: str, destination: Path) -> None:
    headers = {"User-Agent": "metagenomics-best-practices/article-74"}
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


def validate_semantics(output: Path) -> None:
    mag_release = json.loads((output / "mag-release.json").read_text())
    mag_tag = json.loads((output / "mag-tag.json").read_text())
    func_release = json.loads((output / "funcscan-release.json").read_text())
    func_tag = json.loads((output / "funcscan-tag.json").read_text())
    nf_release = json.loads((output / "nextflow-release.json").read_text())
    checkm2 = json.loads((output / "checkm2-record.json").read_text())

    assert mag_release["tag_name"] == "5.5.0"
    assert mag_tag["object"]["sha"] == "56abab5b023ce953c9c43fe21090d156ad0e18af"
    assert func_release["tag_name"] == "4.0.0"
    assert func_tag["object"]["sha"] == "aee3dc965eb0c77267435544dda30da858763913"
    assert nf_release["tag_name"] == "v26.04.0"
    assert checkm2["id"] == 14897628
    assert checkm2["files"][0]["size"] == 1735095710
    assert checkm2["files"][0]["checksum"] == "md5:07c10655620843b517d0df0c160d911f"
    assert "bc93a855e0760aad5c4e5f2d0e26da46" in (
        output / "gunc-progenomes2.1.dmnd.gz.md5"
    ).read_text()
    assert "447c9330056b02f29f30fe81fe4af4eb" in (
        output / "gunc-progenomes2.1.dmnd.md5"
    ).read_text()


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    records: dict[str, dict[str, object]] = {}
    for name, record in SOURCES.items():
        path = output / name
        if not valid(path, record):
            if args.verify_only:
                raise RuntimeError(f"Missing or checksum-mismatched source: {path}")
            download(str(record["url"]), path)
        if not valid(path, record):
            raise RuntimeError(f"Downloaded source failed checksum lock: {path}")
        if name in {"nextflow", "nextflow-26.04.0-dist"}:
            os.chmod(path, 0o755)
        records[name] = dict(record)
        print(f"verified\t{name}\t{record['bytes']}\t{record['sha256']}")

    validate_semantics(output)
    manifest = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "resource_count": len(records),
        "nfcore_mag_release": "5.5.0",
        "nfcore_mag_commit": "56abab5b023ce953c9c43fe21090d156ad0e18af",
        "nfcore_funcscan_release": "4.0.0",
        "nfcore_funcscan_commit": "aee3dc965eb0c77267435544dda30da858763913",
        "nextflow_release": "26.04.0",
        "resources": records,
    }
    (output / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"manifest\t{output / 'download-manifest.json'}")


if __name__ == "__main__":
    main()
