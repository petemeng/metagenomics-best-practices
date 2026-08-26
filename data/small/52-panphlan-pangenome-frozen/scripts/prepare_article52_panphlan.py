#!/usr/bin/env python3
"""Prepare checksum-gated official PanPhlAn 3 E. rectale tutorial assets."""

from __future__ import annotations

import argparse
import bz2
import hashlib
import shutil
import subprocess
import tarfile
from pathlib import Path

from article41_44_utils import dump_json, sha256, write_tsv


PANGENOME_URL = "https://www.dropbox.com/s/oz6uu7tjxuhfkei/Eubacterium_rectale.tar.bz2?dl=1"
MAPS_URL = "https://www.dropbox.com/s/z81e87dvtzp6pu3/panphlan_tutorial_map_results.tar.bz2?dl=1"
PANGENOME_MD5 = "87b111ad1612930d34ab17f9e0f09e34"
PANGENOME_SHA256 = "411d8736790ef591e8131d58e096d68a15c477bc21eb39b41c337af950d4dfc6"
MAPS_SHA256 = "5f9bfecbe3e5459e06f48677076b49a4a37263935ac8cff8422cf1ef6bc9bc50"
PANPHLAN_COMMIT = "4294c3f5c92be9b9ef7d61b69df43e7f27c51601"
WIKI_COMMIT = "328be1e6074015c5183668913ed5e0cf0f879fe5"
SOFTWARE = {
    "panphlan_profiling.py": "901b4dc3710a26145dca064216d62769772b743db55683a956b66253565e7480",
    "misc.py": "d8f7283e847c506178205d4e05da878c5d6a06e1eb767118ca509be4a4af2c15",
}
SAMPLES = (
    ("CCMD34381688ST-21-0", "CCMD34381688ST-21-0", "ZellerG_2014", "DEU"),
    ("G78505", "G78505", "VatanenT_2016", "RUS"),
    ("G88884", "G88884", "SchirmerM_2016", "NLD"),
    ("G88970", "G88970", "SchirmerM_2016", "NLD"),
    ("G89027", "G89027", "SchirmerM_2016", "NLD"),
    ("H2M514903", "H2M514903", "LiJ_2017", "CHN"),
    ("H3M518116", "H3M518116", "LiJ_2017", "CHN"),
    ("HD-1", "HD-1", "QinN_2014", "CHN"),
    ("HD-5", "HD-5", "QinN_2014", "CHN"),
    ("HV-6", "HV-6", "QinN_2014", "CHN"),
    ("LD-48", "LD-48", "QinN_2014", "CHN"),
    ("M1", "M1.48.ST", "BritoIL_2016", "FJI"),
    ("M2_48_ST", "M2.48.ST", "BritoIL_2016", "FJI"),
    ("M2_58_ST", "M2.58.ST", "BritoIL_2016", "FJI"),
    ("N021", "N021", "WenC_2017", "CHN"),
    ("RA023", "RA023", "WenC_2017", "CHN"),
    ("S353", "S353", "KarlssonFH_2013", "SWE"),
    ("SID530054", "SID530054", "FengQ_2015", "AUT"),
    ("SRS011302", "SRS011302", "HMP_2012", "USA"),
    ("SZAXPI003417-4", "SZAXPI003417-4", "YuJ_2015", "CHN"),
    ("T2D-025", "T2D-025", "QinJ_2012", "CHN"),
    ("T2D-063", "T2D-063", "QinJ_2012", "CHN"),
    ("T2D-105", "T2D-105", "QinJ_2012", "CHN"),
    ("W1", "W1.27.ST", "BritoIL_2016", "FJI"),
    ("YSZC12003_36795", "YSZC12003_36795", "XieH_2016", "GBR"),
)


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "curl", "-fL", "--retry", "8", "--retry-all-errors",
            "--connect-timeout", "30", "--max-time", "1800",
            "--output", str(target), url,
        ],
        check=True,
    )


def stage(source: Path | None, url: str, target: Path) -> None:
    if source is None:
        download(url, target)
    else:
        shutil.copy2(source.resolve(), target)


def safe_extract(archive: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    root = output.resolve()
    with tarfile.open(archive, "r:bz2") as handle:
        for member in handle.getmembers():
            candidate = (output / member.name).resolve()
            if root not in candidate.parents and candidate != root:
                raise ValueError(f"Unsafe archive member: {member.name}")
        handle.extractall(output, filter="data")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--pangenome-archive", type=Path)
    parser.add_argument("--map-archive", type=Path)
    parser.add_argument("--panphlan-repo", type=Path)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if work.exists():
        shutil.rmtree(work)
    for folder in ("downloads", "pangenome", "map-archive", "decoded-maps", "software", "logs", "output"):
        (work / folder).mkdir(parents=True)

    pangenome_archive = work / "downloads/Eubacterium_rectale.tar.bz2"
    maps_archive = work / "downloads/panphlan_tutorial_map_results.tar.bz2"
    stage(args.pangenome_archive, PANGENOME_URL, pangenome_archive)
    stage(args.map_archive, MAPS_URL, maps_archive)
    if pangenome_archive.stat().st_size != 83_327_158 or md5(pangenome_archive) != PANGENOME_MD5 or sha256(pangenome_archive) != PANGENOME_SHA256:
        raise ValueError("Official E. rectale pangenome archive checksum/size mismatch")
    if maps_archive.stat().st_size != 3_714_090 or sha256(maps_archive) != MAPS_SHA256:
        raise ValueError("Official PanPhlAn mapping archive checksum/size mismatch")
    safe_extract(pangenome_archive, work / "pangenome")
    safe_extract(maps_archive, work / "map-archive")

    if args.panphlan_repo:
        repository = args.panphlan_repo.resolve()
        for name in SOFTWARE:
            shutil.copy2(repository / name, work / "software" / name)
    else:
        for name in SOFTWARE:
            download(
                f"https://raw.githubusercontent.com/SegataLab/panphlan/{PANPHLAN_COMMIT}/{name}",
                work / "software" / name,
            )
    for name, expected in SOFTWARE.items():
        if sha256(work / "software" / name) != expected:
            raise ValueError(f"Pinned PanPhlAn source checksum mismatch: {name}")

    archived_maps = work / "map-archive/map_results_erectale"
    map_rows: list[dict[str, object]] = []
    for archive_id, sample, study, country in SAMPLES:
        source = archived_maps / f"{archive_id}_erectale.csv.bz2"
        target = work / "decoded-maps" / sample
        if not source.is_file():
            raise FileNotFoundError(f"Missing official map result: {source.name}")
        with bz2.open(source, "rb") as input_handle, target.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle)
        map_rows.append({
            "ArchiveSampleID": archive_id,
            "Sample": sample,
            "Study": study,
            "Country": country,
            "CompressedFile": source.name,
            "CompressedBytes": source.stat().st_size,
            "CompressedSHA256": sha256(source),
            "DecodedBytes": target.stat().st_size,
            "DecodedSHA256": sha256(target),
            "CompatibilityTransformation": "bzip2 decode; canonical tutorial sample ID; no value transformation",
        })
    if len(map_rows) != 25 or len({row["Sample"] for row in map_rows}) != 25:
        raise ValueError("Expected 25 unique tutorial mapping profiles")

    pangenome_dir = work / "pangenome/Eubacterium_rectale"
    required = (
        "Eubacterium_rectale_pangenome.tsv",
        "Eubacterium_rectale_pangenome_contigs.fna",
        "panphlan_Eubacterium_rectale_annot.tsv",
    )
    if any(not (pangenome_dir / name).is_file() for name in required):
        raise FileNotFoundError("Official pangenome archive is incomplete")
    pangenome_rows = []
    for path in sorted(pangenome_dir.iterdir()):
        if path.is_file():
            pangenome_rows.append({
                "File": path.name,
                "Bytes": path.stat().st_size,
                "SHA256": sha256(path),
            })

    write_tsv(work / "sample-metadata.tsv", map_rows)
    write_tsv(work / "pangenome-file-manifest.tsv", pangenome_rows)
    write_tsv(work / "asset-manifest.tsv", [
        {
            "Asset": "Eubacterium rectale pangenome archive",
            "URL": PANGENOME_URL, "Bytes": pangenome_archive.stat().st_size,
            "MD5": md5(pangenome_archive), "SHA256": sha256(pangenome_archive),
        },
        {
            "Asset": "25 official tutorial mapping profiles",
            "URL": MAPS_URL, "Bytes": maps_archive.stat().st_size,
            "MD5": md5(maps_archive), "SHA256": sha256(maps_archive),
        },
    ])
    write_tsv(work / "software-provenance.tsv", [
        {
            "Software": "PanPhlAn", "Version": "3.1",
            "Commit": PANPHLAN_COMMIT, "WikiCommit": WIKI_COMMIT,
            "ProfilingSHA256": SOFTWARE["panphlan_profiling.py"],
        }
    ])
    write_tsv(work / "input-lineage.tsv", [
        {
            "Output": "25 decoded E. rectale gene-coverage profiles",
            "ImmediateInput": "official PanPhlAn 3 tutorial mapping-results archive",
            "Transformation": "checksum-gated bzip2 decode and canonical sample-ID normalization",
            "Evidence": "asset-manifest.tsv; sample-metadata.tsv",
        },
        {
            "Output": "E. rectale 15-reference pangenome",
            "ImmediateInput": "official PanPhlAn pangenome archive",
            "Transformation": "safe tar extraction after official MD5 and pinned SHA-256 validation",
            "Evidence": "asset-manifest.tsv; pangenome-file-manifest.tsv",
        },
    ])
    dump_json(work / "run-contract.json", {
        "article": 52,
        "seed": 20260752,
        "species": "Eubacterium rectale",
        "panphlan_version": "3.1",
        "panphlan_commit": PANPHLAN_COMMIT,
        "wiki_commit": WIKI_COMMIT,
        "samples": 25,
        "studies": 13,
        "countries": 9,
        "reference_genomes": 15,
        "primary_thresholds": {
            "min_coverage": 2.0, "left_max": 1.25, "right_min": 0.75,
            "th_non_present": 0.25, "th_present": 0.5, "th_multicopy": 1.5,
        },
        "sensitivity_thresholds": {
            "min_coverage": 1.0, "left_max": 1.70, "right_min": 0.30,
            "th_non_present": 0.25, "th_present": 0.5, "th_multicopy": 1.5,
        },
        "random_output_requested": False,
        "metadata_used_for_profiling": False,
    })
    (work / ".article52-inputs-complete").write_text("complete\n", encoding="utf-8")
    print("Prepared 25 official PanPhlAn maps and the 15-reference E. rectale pangenome")


if __name__ == "__main__":
    main()
