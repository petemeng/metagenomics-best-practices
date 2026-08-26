#!/usr/bin/env python3
"""Download and checksum the official IBDMDB inputs used in Article 65."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path


IBDMDB = "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb"
MPX = f"{IBDMDB}/products/HMP2/MPX/2017-03-20"
RESOURCES = {
    "hmp2-metadata.csv": {
        "url": f"{IBDMDB}/metadata/hmp2_metadata_2018-08-20.csv",
        "bytes": 9_074_342,
        "sha256": "656b7bd97660ddb875548805e30bede31f2d1208293f7170d2d5755e33862ec9",
    },
    "mgx-ecs-rela.tsv.gz": {
        "url": f"{IBDMDB}/products/HMP2/MGX/2018-05-04/ecs_relab.tsv.gz",
        "bytes": 88_776_917,
        "sha256": "23f81ec1b6a995cfed83224816dc89002f7f7a2880afd796233b1f5386fab220",
    },
    "mtx-ecs-rela.tsv.gz": {
        "url": f"{IBDMDB}/products/HMP2/MTX/2017-12-14/ecs_relab.tsv.gz",
        "bytes": 21_323_580,
        "sha256": "0a5da078c521bd1f7c62d9583b9923cd05f686f344f2fdad75e234134657b3a0",
    },
    "mpx-1pep-1pct.tsv.gz": {
        "url": f"{MPX}/HMP2_proteomics_1pep-pro_1p_FDR_5ppm.tsv.gz",
        "bytes": 1_097_492,
        "sha256": "b4ceee58817059b77a2b010061420f39ce30e1598bd831cdef5cdd771d32dc02",
    },
    "mpx-1pep-5pct.tsv.gz": {
        "url": f"{MPX}/HMP2_proteomics_1pep-pro_5p_FDR_5ppm.tsv.gz",
        "bytes": 1_710_498,
        "sha256": "254cde477a5f95bbd5cf62de0364b0442959fc4e343a8903585f728729cbd111",
    },
    "mpx-2pep-1pct.tsv.gz": {
        "url": f"{MPX}/HMP2_proteomics_2pep-pro_1p_FDR_5ppm.tsv.gz",
        "bytes": 424_364,
        "sha256": "a001f9dfe9417f1d85004a43ec97199a16ad7e55862459f36da769da8632b12b",
    },
    "mpx-2pep-5pct.tsv.gz": {
        "url": f"{MPX}/HMP2_proteomics_2pep-pro_5p_FDR_5ppm.tsv.gz",
        "bytes": 480_606,
        "sha256": "b064d1b6ef3ee49914a5375a656f1d239634254b3694f02da87ead76af925c26",
    },
    "mpx-ecs.tsv.gz": {
        "url": f"{MPX}/HMP2_proteomics_ecs.tsv.gz",
        "bytes": 91_574,
        "sha256": "6328db9bcb40dbe0e927747e7df3dcd0d9974d1fe5e35a9631f3620ae3565a7a",
    },
    "mpx-kos.tsv.gz": {
        "url": f"{MPX}/HMP2_proteomics_kos.tsv.gz",
        "bytes": 146_771,
        "sha256": "e3af45100537405a54b574ec5c8994d837d2f353f629970365d0cd3fe1d42e39",
    },
    "lloyd-price-supp-fig1.pdf": {
        "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1237-9/MediaObjects/41586_2019_1237_MOESM3_ESM.pdf",
        "bytes": 4_278_512,
        "sha256": "78d9845d62be38019e3f58daf2a9b085eda5dbb5a2a1c40491130bcd33e80760",
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
    request = urllib.request.Request(
        url, headers={"User-Agent": "metagenomics-best-practices/1.0"}
    )
    with urllib.request.urlopen(request, timeout=300) as response, temporary.open("wb") as output:
        shutil.copyfileobj(response, output, length=8 * 1024 * 1024)
    temporary.replace(target)


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, object] = {
        "article": 65,
        "dataset": "IBDMDB_HMP2",
        "study": "Lloyd-Price et al. Nature 2019",
        "doi": "10.1038/s41586-019-1237-9",
        "proteomics_release": "2017-03-20",
        "resources": {},
    }
    for name, expected in RESOURCES.items():
        target = cache / name
        if not valid(target, expected):
            download(str(expected["url"]), target)
        if not valid(target, expected):
            raise RuntimeError(f"Checksum or byte-count mismatch: {target}")
        manifest["resources"][name] = {**expected, "downloaded_path": name}
        print(f"verified\t{name}\t{target.stat().st_size}\t{sha256(target)}")
    (cache / "download-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
