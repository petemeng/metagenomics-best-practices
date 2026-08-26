#!/usr/bin/env python3
"""Download and checksum the official 25-sample StrainPhlAn 4 tutorial assets."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tarfile
from pathlib import Path

from article41_44_utils import dump_json, sha256, write_tsv


BASE = "https://cmprod1.cibio.unitn.it/biobakery4/biobakery_strainphlan4"
DATABASE_URL = (
    "https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/"
    "mpa_vJan21_CHOCOPhlAnSGB_202103.tar"
)
DATABASE_NAME = "mpa_vJan21_CHOCOPhlAnSGB_202103.pkl"
SAMPLES = (
    ("CMD64776337ST-21-0", "ZellerG_2014", "DEU"),
    ("G78505", "VatanenT_2016", "RUS"),
    ("G88884", "SchirmerM_2016", "NLD"),
    ("G88970", "SchirmerM_2016", "NLD"),
    ("G89027", "SchirmerM_2016", "NLD"),
    ("H2M514903", "LiJ_2017", "CHN"),
    ("H3M518116", "LiJ_2017", "CHN"),
    ("HD-1", "QinN_2014", "CHN"),
    ("HD-5", "QinN_2014", "CHN"),
    ("HV-6", "QinN_2014", "CHN"),
    ("LD-48", "QinN_2014", "CHN"),
    ("M1.42.ST", "BritoIL_2016", "FJI"),
    ("M2.48.ST", "BritoIL_2016", "FJI"),
    ("M2.58.ST", "BritoIL_2016", "FJI"),
    ("N021", "WenC_2017", "CHN"),
    ("RA023", "WenC_2017", "CHN"),
    ("S353", "KarlssonFH_2013", "SWE"),
    ("SID530054", "FengQ_2015", "AUT"),
    ("SRS011302", "HMP_2012", "USA"),
    ("SZAXPI003417-4", "YuJ_2015", "CHN"),
    ("T2D-025", "QinJ_2012", "CHN"),
    ("T2D-063", "QinJ_2012", "CHN"),
    ("T2D-105", "QinJ_2012", "CHN"),
    ("W1.27.ST", "BritoIL_2016", "FJI"),
    ("YSZC12003_36795", "XieH_2016", "GBR"),
)
REFERENCES = (
    "GCA_000020625.fna.bz2",
    "GCA_000209935.fna.bz2",
    "GCA_001404855.fna.bz2",
    "GCA_001405295.fna.bz2",
    "GCA_001406375.fna.bz2",
    "GCA_001406835.fna.bz2",
)
BASELINE_OUTPUTS = (
    "RAxML_bestTree.t__SGB4933_group.StrainPhlAn4.tre",
    "t__SGB4933_group.StrainPhlAn4_concatenated.aln",
    "t__SGB4933_group.info",
    "t__SGB4933_group.polymorphic",
)


def download(url: str, target: Path, *, byte_range: str | None = None) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "curl", "-fsSL", "--show-error", "--retry", "8", "--retry-all-errors",
        "--connect-timeout", "30", "--max-time", "1800",
    ]
    if byte_range:
        command.extend(["--range", byte_range])
    command.extend(["--output", str(target), url])
    subprocess.run(command, check=True)
    if not target.is_file() or target.stat().st_size == 0:
        raise ValueError(f"Empty download: {url}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if work.exists():
        shutil.rmtree(work)
    for name in ("consensus_markers", "clade_markers", "reference_genomes", "official_baseline", "database", "logs"):
        (work / name).mkdir(parents=True)

    assets: list[dict[str, object]] = []
    for sample, study, country in SAMPLES:
        url = f"{BASE}/consensus_markers/{sample}.pkl"
        target = work / "consensus_markers" / f"{sample}.pkl"
        download(url, target)
        assets.append({
            "AssetType": "consensus marker profile",
            "Name": target.name,
            "URL": url,
            "Bytes": target.stat().st_size,
            "SHA256": sha256(target),
            "OfficialLastModified": "2022-08-25",
        })

    marker_name = "t__SGB4933_group.fna"
    marker_url = f"{BASE}/db_markers/{marker_name}"
    marker = work / "clade_markers" / marker_name
    download(marker_url, marker)
    assets.append({
        "AssetType": "clade marker FASTA", "Name": marker.name,
        "URL": marker_url, "Bytes": marker.stat().st_size,
        "SHA256": sha256(marker), "OfficialLastModified": "2022-08-25",
    })

    for name in REFERENCES:
        url = f"{BASE}/reference_genomes/{name}"
        target = work / "reference_genomes" / name
        download(url, target)
        assets.append({
            "AssetType": "reference genome", "Name": name, "URL": url,
            "Bytes": target.stat().st_size, "SHA256": sha256(target),
            "OfficialLastModified": "2022-08-25",
        })

    for name in BASELINE_OUTPUTS:
        url = f"{BASE}/output/{name}"
        target = work / "official_baseline" / name
        download(url, target)
        assets.append({
            "AssetType": "official precomputed output", "Name": name,
            "URL": url, "Bytes": target.stat().st_size,
            "SHA256": sha256(target), "OfficialLastModified": "2022-08-25",
        })

    # The uncompressed 2.6-GB archive stores its 55,986,225-byte metadata
    # pickle first. A verified HTTP range therefore retrieves all bytes needed
    # by StrainPhlAn without downloading Bowtie2 indices or marker sequences.
    head_tar = work / "database/mpa_vJan21.head.tar"
    download(DATABASE_URL, head_tar, byte_range="0-67108863")
    with tarfile.open(head_tar, mode="r:") as archive:
        member = archive.next()
        if member is None or member.name != DATABASE_NAME:
            raise ValueError(f"Unexpected first Jan21 tar member: {member}")
        if member.size != 55_986_225:
            raise ValueError(f"Unexpected Jan21 database metadata size: {member.size}")
        source = archive.extractfile(member)
        if source is None:
            raise ValueError("Could not extract Jan21 database metadata")
        database = work / "database" / DATABASE_NAME
        with database.open("wb") as output:
            shutil.copyfileobj(source, output)
    head_tar.unlink()
    assets.append({
        "AssetType": "MetaPhlAn metadata extracted by HTTP range",
        "Name": database.name, "URL": DATABASE_URL,
        "Bytes": database.stat().st_size, "SHA256": sha256(database),
        "OfficialLastModified": "2022-04-01",
    })

    write_tsv(work / "asset-manifest.tsv", assets)
    write_tsv(
        work / "sample-metadata.tsv",
        [
            {"Sample": sample, "Study": study, "Country": country,
             "MarkerFile": f"consensus_markers/{sample}.pkl"}
            for sample, study, country in SAMPLES
        ],
    )
    (work / "sample-list.txt").write_text(
        "".join(str(work / "consensus_markers" / f"{sample}.pkl") + "\n" for sample, _, _ in SAMPLES),
        encoding="utf-8",
    )
    (work / "reference-list.txt").write_text(
        "".join(str(work / "reference_genomes" / name) + "\n" for name in REFERENCES),
        encoding="utf-8",
    )
    write_tsv(
        work / "input-lineage.tsv",
        [
            {
                "Output": "25 consensus-marker profiles",
                "ImmediateInput": "bioBakery StrainPhlAn 4 official tutorial assets",
                "Transformation": "HTTPS download only; files are not deserialized during preparation",
                "Evidence": "asset-manifest.tsv; sample-metadata.tsv",
            },
            {
                "Output": "Jan21 MetaPhlAn metadata pickle",
                "ImmediateInput": "official uncompressed Jan21 database tar",
                "Transformation": "HTTP bytes 0-67108863; extract first complete tar member only",
                "Evidence": "asset-manifest.tsv; fixed member size 55,986,225 bytes",
            },
        ],
    )
    dump_json(
        work / "run-contract.json",
        {
            "article": 51,
            "seed": 20260751,
            "tutorial_wiki_commit": "db372805bafe634d5c5f104f33cc9bec11d8bcdf",
            "tutorial_assets_last_modified": "2022-08-25",
            "database_release": "mpa_vJan21_CHOCOPhlAnSGB_202103",
            "target_clade": "t__SGB4933_group",
            "samples": 25,
            "studies": 13,
            "countries": 9,
            "references": 6,
            "trim_sequences": 50,
            "sample_with_n_markers": 20,
            "sample_with_n_markers_perc": 25,
            "marker_in_n_samples_perc": 50,
            "sample_with_n_markers_after_filt": 20,
            "sample_with_n_markers_after_filt_perc": 25,
            "breadth_thres": 80,
            "phylophlan_mode": "fast",
            "official_baseline_threshold_sensitivity": {
                "sample_with_n_markers": 20,
                "sample_with_n_markers_perc": 80,
                "marker_in_n_samples_perc": 80,
                "sample_with_n_markers_after_filt": 20,
                "sample_with_n_markers_after_filt_perc": 80,
                "breadth_thres": 80,
            },
            "truth_used_for_tree": False,
        },
    )
    (work / ".article51-inputs-complete").write_text("complete\n", encoding="utf-8")
    print(f"Prepared {len(SAMPLES)} marker profiles, {len(REFERENCES)} references, and Jan21 metadata")


if __name__ == "__main__":
    main()
