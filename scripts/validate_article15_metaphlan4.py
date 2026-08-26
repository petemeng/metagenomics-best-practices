#!/usr/bin/env python3
"""Validate Article 15 MetaPhlAn 4 profiling and frozen evidence.

Initialization mode is used once with the ignored FASTQ, checksum-locked full
database, and mapout files. Routine mode is network-free: it audits only the
small frozen bundle, regenerates publication tables/figures, and is the QA
entry point declared in tutorial.yaml.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import gc
import gzip
import hashlib
import json
import math
import os
import pickle
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any, BinaryIO, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article15-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


INDEX_NAME = "mpa_vJan26_CHOCOPhlAnSGB_202605"
SUBSAMPLING_SEED = 20260721
EXPECTED_INPUT_PAIRS = 99_991
EXPECTED_INPUT_READS = 199_982
EXPECTED_PROFILING_READS = 199_929
EXPECTED_SHORT_READS = 53
EXPECTED_TAXONOMY_ENTRIES = 72_000
EXPECTED_MARKER_ENTRIES = 13_907_686
EXPECTED_KNOWN_SGBS = 27_891
EXPECTED_UNKNOWN_SGBS = 44_109
EXPECTED_VERSIONS = {
    "MetaPhlAn": "4.2.5",
    "Bowtie2": "2.5.5",
    "Python": "3.12.13",
}
EXPECTED_ENV_YAML_SHA256 = (
    "a314c7d33025e2057ab609f4f1e101b6ab3921338e55c587d318cf2d6d65f874"
)
EXPECTED_ENV_LOCK_SHA256 = (
    "d2884ede2be3ae1c27300505d5a1b6c784a4b30a9f259a3efd9b628a66d302d2"
)
EXPECTED_ENV_PACKAGES = 554
EXPECTED_SOURCE_MANIFEST_SHA256 = (
    "9c0928e60b2431ca549786b5f6302844507862c0f9768ce42e9188ce72a06501"
)
EXPECTED_DATABASE_MANIFEST_SHA256 = (
    "de817d33059b10110b63a16515d630d616eb3fae8790f682eb50d0ceb8fe5713"
)
EXPECTED_NOTICE_SHA256 = (
    "1ec10bdbd9bf314f10a9a099a670814a7d39a931fd9b9b20b324637a858bc8d3"
)
EXPECTED_ARTICLE13_SUMMARY_SHA256 = (
    "f6dcc51b6535247de7f370dc2334994dd85dca61d3f86d252294127faa3460fe"
)
EXPECTED_ARCHIVES = {
    "marker_metadata": {
        "bytes": 6_014_095_360,
        "md5": "7162b0c3493663dce9abef08ccc06aea",
        "sha256": (
            "f945834fb9b8204a0b2235c51a41cd6084cc1965939a9768dd045e5c61c8c946"
        ),
    },
    "bowtie2_index": {
        "bytes": 41_742_510_080,
        "md5": "ac93e1e9c0829629266f5b6ab19c318d",
        "sha256": "0df03a3d877be434c2e89b6cc3201d92c387215700ad66b44585e8883e87bb26",
    },
}
EXPECTED_CLEAN = {
    "R1": {
        "records": 99_991,
        "bases": 14_974_589,
        "shorter_than_70": 22,
        "bytes": 8_661_319,
        "compressed_sha256": (
            "ce438de2814a6e605cc24b00e53b697d018f325bc2a9694cece5be158ff32101"
        ),
        "uncompressed_sha256": (
            "0705a8d7a2c9a4382a49f0b234cb35252ffd539aa53b7608d47b73295f7ba9ca"
        ),
    },
    "R2": {
        "records": 99_991,
        "bases": 14_835_184,
        "shorter_than_70": 31,
        "bytes": 10_045_722,
        "compressed_sha256": (
            "207d080cd5dc98bd10fd4af8c492fb3f2000f73e1458a9d1d77e57347f4fe459"
        ),
        "uncompressed_sha256": (
            "78ed92a8be5eb6e06055e490d179e7e303e1ef481b0c699c6334859ee7d4e8cb"
        ),
    },
}
EXPECTED_PAIR_ID_SHA256 = (
    "457cef6e9d603790dfbc26b716b0498169b54c31bc903d067d449d8dcc86770d"
)
EXPECTED_LOCK_PREFIXES = {
    "metaphlan": "metaphlan-4.2.5-",
    "bowtie2": "bowtie2-2.5.5-",
    "python": "python-3.12.13-",
    "matplotlib-base": "matplotlib-base-3.11.0-",
    "pillow": "pillow-12.3.0-",
    "numpy": "numpy-2.5.1-",
}
FIGURE_STEMS = (
    "15-metaphlan-composition",
    "15-sgb-marker-support",
    "15-detection-quantification-sensitivity",
)
TEXT_SUFFIXES = {".csv", ".json", ".log", ".sh", ".tsv", ".txt"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--environment-prefix", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path)
    parser.add_argument("--raw-dir", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--metadata-archive", type=Path)
    parser.add_argument("--bowtie2-archive", type=Path)
    parser.add_argument("--database-dir", type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--initialize-frozen", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    return parser.parse_args()


def hash_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_file_multi(path: Path, algorithms: tuple[str, ...]) -> dict[str, str]:
    digests = {algorithm: hashlib.new(algorithm) for algorithm in algorithms}
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            for digest in digests.values():
                digest.update(block)
    return {algorithm: digest.hexdigest() for algorithm, digest in digests.items()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path,
    rows: Iterable[dict[str, Any]],
    fields: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
        encoding="utf-8",
    )


def command_output(command: list[str]) -> tuple[int, str]:
    process = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env={**os.environ, "LC_ALL": "C", "TZ": "UTC"},
    )
    return process.returncode, process.stdout.strip()


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict[str, str]] = []

    def add(
        self,
        category: str,
        check_id: str,
        passed: bool,
        detail: str,
    ) -> None:
        self.rows.append(
            {
                "Category": category,
                "CheckID": check_id,
                "Status": "PASS" if passed else "FAIL",
                "Detail": detail,
            }
        )

    @property
    def passed(self) -> int:
        return sum(row["Status"] == "PASS" for row in self.rows)

    @property
    def failed(self) -> int:
        return sum(row["Status"] == "FAIL" for row in self.rows)

    def category(self, name: str) -> list[dict[str, str]]:
        return [row for row in self.rows if row["Category"] == name]


def parse_explicit_lock(path: Path) -> tuple[int, dict[str, str]]:
    packages = [
        line.strip().rsplit("/", 1)[-1].split("#", 1)[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith("https://")
    ]
    targets: dict[str, str] = {}
    for name, prefix in EXPECTED_LOCK_PREFIXES.items():
        matches = [package for package in packages if package.startswith(prefix)]
        targets[name] = matches[0] if len(matches) == 1 else ""
    return len(packages), targets


def observed_tool_versions(
    environment_prefix: Path,
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    commands = {
        "MetaPhlAn": [str(environment_prefix / "bin/metaphlan"), "--version"],
        "Bowtie2": [str(environment_prefix / "bin/bowtie2"), "--version"],
        "Python": [str(environment_prefix / "bin/python"), "--version"],
    }
    patterns = {
        "MetaPhlAn": r"MetaPhlAn version ([0-9.]+)",
        "Bowtie2": r" version ([0-9.]+)",
        "Python": r"Python ([0-9.]+)",
    }
    versions: dict[str, str] = {}
    rows: list[dict[str, Any]] = []
    for tool, command in commands.items():
        return_code, output = command_output(command)
        match = re.search(patterns[tool], output)
        observed = match.group(1) if match else ""
        expected = EXPECTED_VERSIONS[tool]
        versions[tool] = observed
        rows.append(
            {
                "Tool": tool,
                "ExpectedVersion": expected,
                "ObservedVersion": observed,
                "ReturnCode": return_code,
                "Status": (
                    "PASS"
                    if return_code == 0 and observed == expected
                    else "FAIL"
                ),
            }
        )
    return versions, rows


def normalized_read_id(header: bytes) -> bytes:
    token = header[1:].split(None, 1)[0]
    if token.endswith(b"/1") or token.endswith(b"/2"):
        token = token[:-2]
    return token


def read_fastq_record(
    handle: BinaryIO,
    digest: Any,
    path: Path,
    record_number: int,
) -> tuple[bytes, bytes] | None:
    lines = [handle.readline() for _ in range(4)]
    if lines[0] == b"":
        if any(lines[1:]):
            raise ValueError(f"Malformed terminal FASTQ record in {path}")
        return None
    if any(line == b"" for line in lines[1:]):
        raise ValueError(f"Truncated FASTQ record {record_number} in {path}")
    for line in lines:
        digest.update(line)
    header, sequence, plus, quality = [line.rstrip(b"\r\n") for line in lines]
    if not header.startswith(b"@") or not plus.startswith(b"+"):
        raise ValueError(f"Malformed FASTQ record {record_number} in {path}")
    if len(sequence) != len(quality):
        raise ValueError(f"Sequence/quality length mismatch in {path}")
    return header, sequence


def audit_fastq_pair(read1_path: Path, read2_path: Path) -> dict[str, Any]:
    pair_digest = hashlib.sha256()
    mate_digests = {"R1": hashlib.sha256(), "R2": hashlib.sha256()}
    mate_stats = {
        "R1": {"records": 0, "bases": 0, "shorter_than_70": 0},
        "R2": {"records": 0, "bases": 0, "shorter_than_70": 0},
    }
    with gzip.open(read1_path, "rb") as read1, gzip.open(read2_path, "rb") as read2:
        record_number = 0
        while True:
            record1 = read_fastq_record(
                read1, mate_digests["R1"], read1_path, record_number + 1
            )
            record2 = read_fastq_record(
                read2, mate_digests["R2"], read2_path, record_number + 1
            )
            if record1 is None and record2 is None:
                break
            if record1 is None or record2 is None:
                raise ValueError("R1/R2 FASTQ record counts differ")
            read_id1 = normalized_read_id(record1[0])
            read_id2 = normalized_read_id(record2[0])
            if read_id1 != read_id2:
                raise ValueError(
                    f"Unsynchronized paired IDs at record {record_number + 1}"
                )
            pair_digest.update(read_id1)
            pair_digest.update(b"\n")
            record_number += 1
            for mate, record in (("R1", record1), ("R2", record2)):
                sequence = record[1]
                mate_stats[mate]["records"] += 1
                mate_stats[mate]["bases"] += len(sequence)
                mate_stats[mate]["shorter_than_70"] += len(sequence) < 70
    for mate, path in (("R1", read1_path), ("R2", read2_path)):
        mate_stats[mate].update(
            {
                "bytes": path.stat().st_size,
                "compressed_sha256": hash_file(path),
                "uncompressed_sha256": mate_digests[mate].hexdigest(),
            }
        )
    return {
        "pairs": mate_stats["R1"]["records"],
        "mates_synchronized": True,
        "normalized_pair_id_sha256": pair_digest.hexdigest(),
        "R1": mate_stats["R1"],
        "R2": mate_stats["R2"],
    }


def parse_profile(path: Path) -> dict[str, Any]:
    headers: list[str] = []
    fields: list[str] | None = None
    rows: list[dict[str, str]] = []
    processed_reads: int | None = None
    estimated_mapped_reads: int | None = None
    sample_id = ""
    database = ""
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                headers.append(line)
                if line == f"#{INDEX_NAME}":
                    database = INDEX_NAME
                match = re.fullmatch(r"#([0-9]+) reads processed", line)
                if match:
                    processed_reads = int(match.group(1))
                match = re.fullmatch(
                    r"#Estimated reads mapped to known clades: ([0-9]+)", line
                )
                if match:
                    estimated_mapped_reads = int(match.group(1))
                if line.startswith("#SampleID\t"):
                    sample_id = line.split("\t", 1)[1]
                if line.startswith("#clade_name\t"):
                    fields = line[1:].split("\t")
                continue
            if fields is None:
                raise ValueError(f"Profile rows precede a column header in {path}")
            values = line.split("\t")
            if len(values) != len(fields):
                raise ValueError(f"Unexpected profile column count in {path}: {line}")
            rows.append(dict(zip(fields, values)))
    if processed_reads is None or estimated_mapped_reads is None or fields is None:
        raise ValueError(f"Incomplete MetaPhlAn profile headers in {path}")
    return {
        "path": path,
        "headers": headers,
        "rows": rows,
        "processed_reads": processed_reads,
        "estimated_mapped_reads": estimated_mapped_reads,
        "sample_id": sample_id,
        "database": database,
    }


def lineage_component(lineage: str, prefix: str) -> str:
    for component in lineage.split("|"):
        if component.startswith(prefix):
            return component
    return ""


def exact_rank_rows(profile: dict[str, Any], rank: str) -> list[dict[str, str]]:
    prefix = f"{rank}__"
    return [
        row
        for row in profile["rows"]
        if row["clade_name"].split("|")[-1].startswith(prefix)
    ]


def unclassified_row(profile: dict[str, Any]) -> dict[str, str] | None:
    return next(
        (
            row
            for row in profile["rows"]
            if row["clade_name"] == "UNCLASSIFIED"
        ),
        None,
    )


def sgb_class(lineage: str) -> str:
    species = lineage_component(lineage, "s__")
    return "Unknown SGB" if "_SGB" in species else "Known SGB"


def display_taxon(label: str) -> str:
    if "__" in label:
        label = label.split("__", 1)[1]
    return label.replace("_", " ")


def parse_mapout(path: Path) -> dict[str, Any]:
    opener = bz2.open if path.suffix == ".bz2" else open
    marker_hits: Counter[str] = Counter()
    mapped_records = 0
    nreads: int | None = None
    average_read_length: float | None = None
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            read_id, value = line.split("\t", 1)
            if read_id == "#nreads":
                nreads = int(value)
            elif read_id == "#avg_read_length":
                average_read_length = float(value)
            else:
                marker_hits[value] += 1
                mapped_records += 1
    if nreads is None or average_read_length is None:
        raise ValueError(f"Mapout footer is incomplete: {path}")
    return {
        "path": path,
        "nreads": nreads,
        "average_read_length": average_read_length,
        "mapped_records": mapped_records,
        "distinct_markers": len(marker_hits),
        "marker_hits": marker_hits,
        "bytes": path.stat().st_size,
        "sha256": hash_file(path),
    }


def profile_table_rows(
    profile: dict[str, Any],
    rank: str,
    pair_depth: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in exact_rank_rows(profile, rank):
        lineage = row["clade_name"]
        rows.append(
            {
                "SampleID": profile["sample_id"],
                "PairDepth": pair_depth,
                "ProfilingReads": profile["processed_reads"],
                "Rank": "Species" if rank == "s" else "SGB",
                "CladeName": lineage,
                "TaxID": row["clade_taxid"],
                "SpeciesLabel": lineage_component(lineage, "s__"),
                "SGBLabel": lineage_component(lineage, "t__"),
                "SGBClass": sgb_class(lineage),
                "RelativeAbundancePct": float(row["relative_abundance"]),
                "Coverage": float(row["coverage"]),
                "EstimatedReads": int(row["estimated_number_of_reads_from_the_clade"]),
            }
        )
    return rows


def build_species_and_sgb_tables(
    profile: dict[str, Any], pair_depth: int
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    species_rows = profile_table_rows(profile, "s", pair_depth)
    sgb_rows = profile_table_rows(profile, "t", pair_depth)
    unclassified = unclassified_row(profile)
    if unclassified is not None:
        species_rows.append(
            {
                "SampleID": profile["sample_id"],
                "PairDepth": pair_depth,
                "ProfilingReads": profile["processed_reads"],
                "Rank": "Unclassified",
                "CladeName": "UNCLASSIFIED",
                "TaxID": "-1",
                "SpeciesLabel": "UNCLASSIFIED",
                "SGBLabel": "",
                "SGBClass": "Unclassified",
                "RelativeAbundancePct": float(unclassified["relative_abundance"]),
                "Coverage": "",
                "EstimatedReads": int(
                    unclassified["estimated_number_of_reads_from_the_clade"]
                ),
            }
        )
    return species_rows, sgb_rows


PROFILE_FIELDS = [
    "SampleID",
    "PairDepth",
    "ProfilingReads",
    "Rank",
    "CladeName",
    "TaxID",
    "SpeciesLabel",
    "SGBLabel",
    "SGBClass",
    "RelativeAbundancePct",
    "Coverage",
    "EstimatedReads",
]


def build_threshold_sensitivity(
    species_rows: list[dict[str, Any]],
    sgb_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for rank, rows in (("Species", species_rows), ("SGB", sgb_rows)):
        classified = [row for row in rows if row["Rank"] != "Unclassified"]
        classified_total = sum(
            float(row["RelativeAbundancePct"]) for row in classified
        )
        for threshold in (0.0, 0.01, 0.1, 1.0):
            retained = [
                row
                for row in classified
                if float(row["RelativeAbundancePct"]) >= threshold
            ]
            retained_abundance = sum(
                float(row["RelativeAbundancePct"]) for row in retained
            )
            output.append(
                {
                    "Rank": rank,
                    "ThresholdPct": threshold,
                    "DetectedClades": len(retained),
                    "RetainedAbundancePct": retained_abundance,
                    "ExcludedAbundancePct": classified_total - retained_abundance,
                    "RetainedFractionOfClassified": (
                        retained_abundance / classified_total
                        if classified_total
                        else 0.0
                    ),
                    "Rule": "post_profile_only",
                }
            )
    return output


THRESHOLD_FIELDS = [
    "Rank",
    "ThresholdPct",
    "DetectedClades",
    "RetainedAbundancePct",
    "ExcludedAbundancePct",
    "RetainedFractionOfClassified",
    "Rule",
]


def build_depth_sensitivity(
    profiles: dict[int, dict[str, Any]],
    subsample_audits: dict[int, dict[str, Any]],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for pair_depth in sorted(profiles):
        profile = profiles[pair_depth]
        species_rows, sgb_rows = build_species_and_sgb_tables(
            profile, pair_depth
        )
        classified_species = [
            row for row in species_rows if row["Rank"] == "Species"
        ]
        unclassified = next(
            row for row in species_rows if row["Rank"] == "Unclassified"
        )
        top_species = max(
            classified_species,
            key=lambda row: float(row["RelativeAbundancePct"]),
        )
        audit = subsample_audits[pair_depth]
        output.append(
            {
                "PairDepth": pair_depth,
                "InputReads": pair_depth * 2,
                "ProfilingReads": profile["processed_reads"],
                "ShortReadsExcluded": pair_depth * 2 - profile["processed_reads"],
                "DetectedSpecies": len(classified_species),
                "DetectedSGBs": len(sgb_rows),
                "KnownSGBs": sum(row["SGBClass"] == "Known SGB" for row in sgb_rows),
                "UnknownSGBs": sum(row["SGBClass"] == "Unknown SGB" for row in sgb_rows),
                "UnclassifiedPct": float(unclassified["RelativeAbundancePct"]),
                "ClassifiedPct": 100.0 - float(unclassified["RelativeAbundancePct"]),
                "TopSpecies": top_species["SpeciesLabel"],
                "TopSpeciesAbundancePct": float(
                    top_species["RelativeAbundancePct"]
                ),
                "SubsamplingSeed": (
                    SUBSAMPLING_SEED
                    if pair_depth < EXPECTED_INPUT_PAIRS
                    else "not_applicable_all_reads"
                ),
                "PairedSubsamplingArgument": (
                    pair_depth * 2
                    if pair_depth < EXPECTED_INPUT_PAIRS
                    else "not_applicable_all_reads"
                ),
                "NormalizedPairIDSHA256": audit["normalized_pair_id_sha256"],
            }
        )
    return output


DEPTH_FIELDS = [
    "PairDepth",
    "InputReads",
    "ProfilingReads",
    "ShortReadsExcluded",
    "DetectedSpecies",
    "DetectedSGBs",
    "KnownSGBs",
    "UnknownSGBs",
    "UnclassifiedPct",
    "ClassifiedPct",
    "TopSpecies",
    "TopSpeciesAbundancePct",
    "SubsamplingSeed",
    "PairedSubsamplingArgument",
    "NormalizedPairIDSHA256",
]


def load_database_evidence(
    pkl_path: Path,
    profiles: dict[int, dict[str, Any]],
    mapout_audits: dict[int, dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    with bz2.open(pkl_path, "rb") as handle:
        database = pickle.load(handle)

    taxonomy = database["taxonomy"]
    markers = database["markers"]
    taxonomy_sgb_labels = {
        lineage_component(lineage, "t__") for lineage in taxonomy
    }
    known_sgbs = 0
    unknown_sgbs = 0
    non_sgb_leaves = 0
    for lineage in taxonomy:
        if not lineage_component(lineage, "t__"):
            non_sgb_leaves += 1
        elif "_SGB" in lineage_component(lineage, "s__"):
            unknown_sgbs += 1
        else:
            known_sgbs += 1

    marker_totals: Counter[str] = Counter()
    for marker_info in markers.values():
        marker_totals[marker_info["clade"]] += 1

    support_rows: list[dict[str, Any]] = []
    for pair_depth, profile in sorted(profiles.items()):
        observed_markers: Counter[str] = Counter()
        observed_hits: Counter[str] = Counter()
        viral_marker_names = 0
        viral_marker_hits = 0
        metadata_excluded_sgb_marker_names = 0
        metadata_excluded_sgb_marker_hits = 0
        metadata_excluded_sgb_labels: set[str] = set()
        unexpected_marker_names = 0
        for marker_name, hit_count in mapout_audits[pair_depth][
            "marker_hits"
        ].items():
            marker_info = markers.get(marker_name)
            if marker_info is None:
                if marker_name.startswith("VDB|"):
                    viral_marker_names += 1
                    viral_marker_hits += hit_count
                else:
                    sgb_match = re.search(r"\|(SGB[0-9]+(?:_group)?)$", marker_name)
                    sgb_label = sgb_match.group(1) if sgb_match else ""
                    if f"t__{sgb_label}" in taxonomy_sgb_labels:
                        metadata_excluded_sgb_marker_names += 1
                        metadata_excluded_sgb_marker_hits += hit_count
                        metadata_excluded_sgb_labels.add(sgb_label)
                    else:
                        unexpected_marker_names += 1
                continue
            clade = marker_info["clade"]
            observed_markers[clade] += 1
            observed_hits[clade] += hit_count

        _, sgb_rows = build_species_and_sgb_tables(profile, pair_depth)
        for row in sgb_rows:
            sgb_label = row["SGBLabel"]
            total_markers = marker_totals[sgb_label]
            nonzero_markers = observed_markers[sgb_label]
            support_rows.append(
                {
                    "PairDepth": pair_depth,
                    "SGBLineage": row["CladeName"],
                    "SGBLabel": sgb_label,
                    "SpeciesLabel": row["SpeciesLabel"],
                    "SGBClass": row["SGBClass"],
                    "TotalDatabaseMarkers": total_markers,
                    "NonzeroMarkers": nonzero_markers,
                    "MarkerSupportFraction": (
                        nonzero_markers / total_markers if total_markers else 0.0
                    ),
                    "MappedMarkerHits": observed_hits[sgb_label],
                    "Coverage": row["Coverage"],
                    "RelativeAbundancePct": row["RelativeAbundancePct"],
                    "EstimatedReads": row["EstimatedReads"],
                    "PercNonzero": 0.33,
                    "PercNonzeroRole": (
                        "quasi_marker_disambiguation_not_detection_limit"
                    ),
                    "ViralMarkerNamesInMapout": viral_marker_names,
                    "ViralMarkerHitsInMapout": viral_marker_hits,
                    "MetadataExcludedSGBMarkerNamesInMapout": (
                        metadata_excluded_sgb_marker_names
                    ),
                    "MetadataExcludedSGBMarkerHitsInMapout": (
                        metadata_excluded_sgb_marker_hits
                    ),
                    "MetadataExcludedSGBLabelsInMapout": ";".join(
                        sorted(metadata_excluded_sgb_labels)
                    ),
                    "UnexpectedMarkerNamesInMapout": unexpected_marker_names,
                }
            )

    summary = {
        "taxonomy_entries": len(taxonomy),
        "marker_entries": len(markers),
        "known_sgbs": known_sgbs,
        "unknown_sgbs": unknown_sgbs,
        "non_sgb_leaves": non_sgb_leaves,
    }
    del database, taxonomy, markers, marker_totals
    gc.collect()
    return summary, support_rows


MARKER_FIELDS = [
    "PairDepth",
    "SGBLineage",
    "SGBLabel",
    "SpeciesLabel",
    "SGBClass",
    "TotalDatabaseMarkers",
    "NonzeroMarkers",
    "MarkerSupportFraction",
    "MappedMarkerHits",
    "Coverage",
    "RelativeAbundancePct",
    "EstimatedReads",
    "PercNonzero",
    "PercNonzeroRole",
    "ViralMarkerNamesInMapout",
    "ViralMarkerHitsInMapout",
    "MetadataExcludedSGBMarkerNamesInMapout",
    "MetadataExcludedSGBMarkerHitsInMapout",
    "MetadataExcludedSGBLabelsInMapout",
    "UnexpectedMarkerNamesInMapout",
]


def resource_summary(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    values: dict[str, Any] = {}
    patterns = {
        "percent_cpu": r"Percent of CPU this job got:\s*(.+)",
        "elapsed": r"Elapsed \(wall clock\) time \([^)]*\):\s*(.+)",
        "maximum_rss_kb": r"Maximum resident set size \(kbytes\):\s*([0-9]+)",
    }
    text = path.read_text(encoding="utf-8", errors="replace")
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            values[key] = (
                int(match.group(1))
                if key == "maximum_rss_kb"
                else match.group(1).strip()
            )
    return values


def write_database_file_manifest(database_dir: Path, output_path: Path) -> dict[str, Any]:
    rows: list[tuple[str, int, str]] = []
    for path in sorted(item for item in database_dir.iterdir() if item.is_file()):
        rows.append((hash_file(path), path.stat().st_size, path.name))
    output_path.write_text(
        "".join(
            f"{digest}  {size}  ${{METAPHLAN_DB_DIR}}/{name}\n"
            for digest, size, name in rows
        ),
        encoding="utf-8",
    )
    return {
        "file_count": len(rows),
        "total_bytes": sum(size for _, size, _ in rows),
        "files": [name for _, _, name in rows],
    }


def parse_database_file_manifest(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        digest, size, name = line.split("  ", 2)
        rows.append({"SHA256": digest, "Bytes": int(size), "Path": name})
    return rows


def normalize_frozen_paths(
    frozen_dir: Path,
    replacements: list[tuple[Path | str, str]],
) -> dict[str, Any]:
    string_replacements = sorted(
        ((str(source), target) for source, target in replacements),
        key=lambda item: len(item[0]),
        reverse=True,
    )
    changed_files = 0
    replacement_count = 0
    for path in sorted(item for item in frozen_dir.rglob("*") if item.is_file()):
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="strict")
        updated = text
        for source, target in string_replacements:
            count = updated.count(source)
            if count:
                replacement_count += count
                updated = updated.replace(source, target)
        if updated != text:
            path.write_text(updated, encoding="utf-8")
            changed_files += 1
    return {
        "status": "passed",
        "changed_files": changed_files,
        "replacement_count": replacement_count,
        "placeholders": sorted({target for _, target in string_replacements}),
    }


def host_path_hits(frozen_dir: Path, patterns: list[str]) -> list[str]:
    hits: list[str] = []
    for path in sorted(item for item in frozen_dir.rglob("*") if item.is_file()):
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for pattern in patterns:
            if pattern and pattern in text:
                hits.append(f"{path.relative_to(frozen_dir)}:{pattern}")
    return hits


def write_frozen_checksums(frozen_dir: Path) -> None:
    checksum_path = frozen_dir / "file-checksums.sha256"
    rows: list[str] = []
    for path in sorted(item for item in frozen_dir.rglob("*") if item.is_file()):
        if path == checksum_path:
            continue
        rows.append(f"{hash_file(path)}  {path.relative_to(frozen_dir).as_posix()}\n")
    checksum_path.write_text("".join(rows), encoding="utf-8")


def verify_frozen_checksums(frozen_dir: Path) -> tuple[int, list[str]]:
    checksum_path = frozen_dir / "file-checksums.sha256"
    failures: list[str] = []
    entries = 0
    listed_paths: set[str] = set()
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        entries += 1
        listed_paths.add(relative)
        path = frozen_dir / relative
        if not path.is_file():
            failures.append(f"missing:{relative}")
        elif hash_file(path) != digest:
            failures.append(f"sha256:{relative}")
    actual_paths = {
        path.relative_to(frozen_dir).as_posix()
        for path in frozen_dir.rglob("*")
        if path.is_file() and path != checksum_path
    }
    failures.extend(
        f"unexpected:{relative}" for relative in sorted(actual_paths - listed_paths)
    )
    return entries, failures


def archive_audit(path: Path) -> dict[str, Any]:
    digests = hash_file_multi(path, ("md5", "sha256"))
    return {
        "bytes": path.stat().st_size,
        "md5": digests["md5"],
        "sha256": digests["sha256"],
        "name": path.name,
    }


def source_manifest_checks(project_root: Path, checks: Checks) -> None:
    source_path = project_root / "data/small/15-source-manifest.tsv"
    rows = read_tsv(source_path)
    checks.add(
        "source",
        "source-manifest-sha256",
        hash_file(source_path) == EXPECTED_SOURCE_MANIFEST_SHA256,
        hash_file(source_path),
    )
    checks.add("source", "source-manifest-two-mates", len(rows) == 2, str(len(rows)))
    for field, expected in (
        ("ProjectAccession", "PRJEB52977"),
        ("SampleAccession", "SAMEA14435832"),
        ("RunAccession", "ERR9765746"),
        ("Layout", "PAIRED"),
        ("CleanPairs", str(EXPECTED_INPUT_PAIRS)),
        ("NormalizedPairIDSHA256", EXPECTED_PAIR_ID_SHA256),
    ):
        values = {row[field] for row in rows}
        checks.add(
            "source",
            f"source-{field}",
            values == {expected},
            ",".join(sorted(values)),
        )
    article13_path = project_root / "data/small/13-qc-frozen/run-summary.json"
    article13 = json.loads(article13_path.read_text(encoding="utf-8"))
    checks.add(
        "source",
        "article13-summary-sha256",
        hash_file(article13_path) == EXPECTED_ARTICLE13_SUMMARY_SHA256,
        hash_file(article13_path),
    )
    checks.add(
        "source",
        "article13-retained-pairs",
        article13.get("retained_pairs") == EXPECTED_INPUT_PAIRS,
        str(article13.get("retained_pairs")),
    )


def database_manifest_checks(project_root: Path, checks: Checks) -> list[dict[str, str]]:
    path = project_root / "data/small/15-database-manifest.tsv"
    rows = read_tsv(path)
    observed_hash = hash_file(path)
    checks.add(
        "database",
        "database-manifest-sha256",
        observed_hash == EXPECTED_DATABASE_MANIFEST_SHA256,
        observed_hash,
    )
    checks.add("database", "database-manifest-two-assets", len(rows) == 2, str(len(rows)))
    by_role = {row["AssetRole"]: row for row in rows}
    for role, expected in EXPECTED_ARCHIVES.items():
        row = by_role.get(role, {})
        checks.add(
            "database",
            f"{role}-index",
            row.get("IndexName") == INDEX_NAME,
            row.get("IndexName", "missing"),
        )
        checks.add(
            "database",
            f"{role}-bytes",
            row.get("ArchiveBytes") == str(expected["bytes"]),
            row.get("ArchiveBytes", "missing"),
        )
        checks.add(
            "database",
            f"{role}-md5",
            row.get("OfficialMD5") == expected["md5"],
            row.get("OfficialMD5", "missing"),
        )
        checks.add(
            "database",
            f"{role}-sha256",
            row.get("ArchiveSHA256") == expected["sha256"],
            row.get("ArchiveSHA256", "missing"),
        )
        checks.add(
            "database",
            f"{role}-status",
            row.get("Status") == "VERIFIED_AND_EXTRACTED",
            row.get("Status", "missing"),
        )
    return rows


def environment_checks(
    project_root: Path,
    environment_prefix: Path,
    checks: Checks,
) -> tuple[dict[str, str], list[dict[str, Any]], dict[str, str]]:
    versions, tool_rows = observed_tool_versions(environment_prefix)
    for row in tool_rows:
        checks.add(
            "tool",
            f"tool-{row['Tool'].lower()}",
            row["Status"] == "PASS",
            f"expected={row['ExpectedVersion']};observed={row['ObservedVersion']}",
        )
    env_yaml = project_root / "env/biobakery.yml"
    env_lock = project_root / "env/biobakery-linux-64.lock"
    checks.add(
        "tool",
        "environment-yaml-sha256",
        hash_file(env_yaml) == EXPECTED_ENV_YAML_SHA256,
        hash_file(env_yaml),
    )
    checks.add(
        "tool",
        "environment-lock-sha256",
        hash_file(env_lock) == EXPECTED_ENV_LOCK_SHA256,
        hash_file(env_lock),
    )
    package_count, lock_targets = parse_explicit_lock(env_lock)
    checks.add(
        "tool",
        "environment-lock-package-count",
        package_count == EXPECTED_ENV_PACKAGES,
        str(package_count),
    )
    for package, observed in lock_targets.items():
        checks.add(
            "tool",
            f"lock-{package}",
            bool(observed),
            observed or "missing_or_duplicated",
        )
    return versions, tool_rows, lock_targets


def fastq_checks(audit: dict[str, Any], checks: Checks, prefix: str = "clean") -> None:
    checks.add(
        "source",
        f"{prefix}-pair-count",
        audit["pairs"] == EXPECTED_INPUT_PAIRS,
        str(audit["pairs"]),
    )
    checks.add(
        "source",
        f"{prefix}-paired-id-sha256",
        audit["normalized_pair_id_sha256"] == EXPECTED_PAIR_ID_SHA256,
        audit["normalized_pair_id_sha256"],
    )
    for mate in ("R1", "R2"):
        for key, expected in EXPECTED_CLEAN[mate].items():
            observed = audit[mate][key]
            checks.add(
                "source",
                f"{prefix}-{mate.lower()}-{key.replace('_', '-')}",
                observed == expected,
                str(observed),
            )


def required_initialization_args(args: argparse.Namespace) -> list[Path]:
    values = [
        args.frozen_dir,
        args.raw_dir,
        args.work_dir,
        args.metadata_archive,
        args.bowtie2_archive,
        args.database_dir,
    ]
    if any(value is None for value in values):
        raise ValueError(
            "Initialization/preflight requires frozen, raw, work, archive, and database paths"
        )
    return [Path(value).resolve() for value in values]


def run_preflight(args: argparse.Namespace) -> None:
    (
        frozen_dir,
        raw_dir,
        work_dir,
        metadata_archive,
        bowtie2_archive,
        database_dir,
    ) = required_initialization_args(args)
    del frozen_dir, raw_dir, work_dir
    project_root = args.project_root.resolve()
    environment_prefix = args.environment_prefix.resolve()
    checks = Checks()
    source_manifest_checks(project_root, checks)
    database_manifest_checks(project_root, checks)
    environment_checks(project_root, environment_prefix, checks)

    clean_audit = audit_fastq_pair(
        project_root / "data/raw/article13/ERR9765746_clean_R1.fastq.gz",
        project_root / "data/raw/article13/ERR9765746_clean_R2.fastq.gz",
    )
    fastq_checks(clean_audit, checks)

    for role, archive in (
        ("marker_metadata", metadata_archive),
        ("bowtie2_index", bowtie2_archive),
    ):
        checks.add(
            "database",
            f"preflight-{role}-bytes",
            archive.is_file()
            and archive.stat().st_size == EXPECTED_ARCHIVES[role]["bytes"],
            str(archive.stat().st_size if archive.exists() else "missing"),
        )
    pkl_path = database_dir / f"{INDEX_NAME}.pkl"
    index_files = sorted(database_dir.glob(f"{INDEX_NAME}*.bt2*"))
    checks.add("database", "preflight-pkl", pkl_path.is_file(), pkl_path.name)
    checks.add(
        "database",
        "preflight-six-bowtie2-files",
        len(index_files) == 6,
        str(len(index_files)),
    )
    payload = {
        "status": "passed" if checks.failed == 0 else "failed",
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "checks": checks.rows,
    }
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if checks.failed:
        raise SystemExit(1)


def profile_checks(
    profiles: dict[int, dict[str, Any]],
    classified_only: dict[str, Any],
    checks: Checks,
) -> None:
    full = profiles[EXPECTED_INPUT_PAIRS]
    checks.add(
        "profile",
        "full-profile-database",
        full["database"] == INDEX_NAME,
        full["database"],
    )
    checks.add(
        "profile",
        "full-profile-processed-reads",
        full["processed_reads"] == EXPECTED_PROFILING_READS,
        str(full["processed_reads"]),
    )
    command_header = next(
        (header for header in full["headers"] if " --input_type " in header),
        "",
    )
    for required in (
        "--input_type fastq",
        "--db_dir",
        f"--index {INDEX_NAME}",
        "--offline",
        "--mapout",
        "--nproc 8",
        "--read_min_len 70",
        "--perc_nonzero 0.33",
        "-t rel_ab_w_read_stats",
    ):
        checks.add(
            "profile",
            f"main-command-{re.sub('[^a-z0-9]+', '-', required.lower()).strip('-')}",
            required in command_header,
            required,
        )
    unclassified = unclassified_row(full)
    checks.add(
        "profile",
        "default-unclassified-present",
        unclassified is not None,
        "present" if unclassified else "missing",
    )
    if unclassified is not None:
        for rank in ("s", "t"):
            total = sum(
                float(row["relative_abundance"])
                for row in exact_rank_rows(full, rank)
            ) + float(unclassified["relative_abundance"])
            checks.add(
                "profile",
                f"default-{rank}-sum-100",
                math.isclose(total, 100.0, abs_tol=1e-3),
                f"{total:.8f}",
            )

    checks.add(
        "profile",
        "classified-only-same-processed-reads",
        classified_only["processed_reads"] == full["processed_reads"],
        str(classified_only["processed_reads"]),
    )
    checks.add(
        "profile",
        "classified-only-no-unclassified-row",
        unclassified_row(classified_only) is None,
        "absent" if unclassified_row(classified_only) is None else "present",
    )
    for rank in ("s", "t"):
        default_rows = {row["clade_name"]: row for row in exact_rank_rows(full, rank)}
        known_rows = {
            row["clade_name"]: row
            for row in exact_rank_rows(classified_only, rank)
        }
        checks.add(
            "profile",
            f"classified-only-{rank}-same-clades",
            set(default_rows) == set(known_rows),
            f"default={len(default_rows)};classified_only={len(known_rows)}",
        )
        total = sum(float(row["relative_abundance"]) for row in known_rows.values())
        checks.add(
            "profile",
            f"classified-only-{rank}-sum-100",
            math.isclose(total, 100.0, abs_tol=1e-3),
            f"{total:.8f}",
        )
        evidence_equal = all(
            math.isclose(
                float(default_rows[name]["coverage"]),
                float(known_rows[name]["coverage"]),
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
            and default_rows[name]["estimated_number_of_reads_from_the_clade"]
            == known_rows[name]["estimated_number_of_reads_from_the_clade"]
            for name in set(default_rows) & set(known_rows)
        )
        checks.add(
            "profile",
            f"classified-only-{rank}-same-evidence",
            evidence_equal,
            "coverage_and_estimated_reads",
        )

    for pair_depth, profile in profiles.items():
        expected_input_reads = pair_depth * 2
        checks.add(
            "sensitivity",
            f"depth-{pair_depth}-processed-not-over-input",
            0 < profile["processed_reads"] <= expected_input_reads,
            f"processed={profile['processed_reads']};input={expected_input_reads}",
        )
        checks.add(
            "sensitivity",
            f"depth-{pair_depth}-species-positive",
            len(exact_rank_rows(profile, "s")) > 0,
            str(len(exact_rank_rows(profile, "s"))),
        )
        checks.add(
            "sensitivity",
            f"depth-{pair_depth}-sgb-positive",
            len(exact_rank_rows(profile, "t")) > 0,
            str(len(exact_rank_rows(profile, "t"))),
        )


def initialize_frozen(args: argparse.Namespace) -> None:
    (
        frozen_dir,
        raw_dir,
        work_dir,
        metadata_archive,
        bowtie2_archive,
        database_dir,
    ) = required_initialization_args(args)
    project_root = args.project_root.resolve()
    environment_prefix = args.environment_prefix.resolve()
    checks = Checks()

    source_manifest_checks(project_root, checks)
    database_manifest_checks(project_root, checks)
    versions, _, lock_targets = environment_checks(
        project_root, environment_prefix, checks
    )
    notice_path = project_root / "data/small/15-data-NOTICE.txt"
    checks.add(
        "source",
        "notice-sha256",
        hash_file(notice_path) == EXPECTED_NOTICE_SHA256,
        hash_file(notice_path),
    )

    clean_r1 = project_root / "data/raw/article13/ERR9765746_clean_R1.fastq.gz"
    clean_r2 = project_root / "data/raw/article13/ERR9765746_clean_R2.fastq.gz"
    clean_audit = audit_fastq_pair(clean_r1, clean_r2)
    fastq_checks(clean_audit, checks)

    archive_audits = {
        "marker_metadata": archive_audit(metadata_archive),
        "bowtie2_index": archive_audit(bowtie2_archive),
    }
    for role, observed in archive_audits.items():
        expected = EXPECTED_ARCHIVES[role]
        for key in ("bytes", "md5", "sha256"):
            checks.add(
                "database",
                f"archive-{role}-{key}",
                observed[key] == expected[key],
                str(observed[key]),
            )

    database_file_summary = write_database_file_manifest(
        database_dir, frozen_dir / "database-files.sha256"
    )
    database_files = parse_database_file_manifest(
        frozen_dir / "database-files.sha256"
    )
    checks.add(
        "database",
        "database-pkl-present",
        any(row["Path"].endswith(f"/{INDEX_NAME}.pkl") for row in database_files),
        str(len(database_files)),
    )
    index_files = [row for row in database_files if ".bt2" in row["Path"]]
    checks.add(
        "database",
        "database-six-bowtie2-files",
        len(index_files) == 6,
        str(len(index_files)),
    )

    full_profile = parse_profile(frozen_dir / "profile-all.tsv")
    classified_only = parse_profile(
        frozen_dir / "profile-classified-only.tsv"
    )
    profiles = {
        20_000: parse_profile(work_dir / "profile-depth-20000.tsv"),
        50_000: parse_profile(work_dir / "profile-depth-50000.tsv"),
        EXPECTED_INPUT_PAIRS: full_profile,
    }
    profile_checks(profiles, classified_only, checks)

    mapout_paths = {
        20_000: work_dir / "ERR9765746-depth-20000.mapout.bz2",
        50_000: work_dir / "ERR9765746-depth-50000.mapout.bz2",
        EXPECTED_INPUT_PAIRS: work_dir / "ERR9765746-full.mapout.bz2",
    }
    mapout_audits = {
        pair_depth: parse_mapout(path)
        for pair_depth, path in mapout_paths.items()
    }
    for pair_depth, audit in mapout_audits.items():
        checks.add(
            "profile",
            f"mapout-{pair_depth}-profile-read-denominator",
            audit["nreads"] == profiles[pair_depth]["processed_reads"],
            f"mapout={audit['nreads']};profile={profiles[pair_depth]['processed_reads']}",
        )
        checks.add(
            "profile",
            f"mapout-{pair_depth}-marker-hits-positive",
            audit["mapped_records"] > 0,
            str(audit["mapped_records"]),
        )

    subsample_audits = {EXPECTED_INPUT_PAIRS: clean_audit}
    for pair_depth in (20_000, 50_000):
        audit = audit_fastq_pair(
            work_dir / f"subsample-{pair_depth}.fastq.R1.gz",
            work_dir / f"subsample-{pair_depth}.fastq.R2.gz",
        )
        subsample_audits[pair_depth] = audit
        checks.add(
            "sensitivity",
            f"subsample-{pair_depth}-pair-count",
            audit["pairs"] == pair_depth,
            str(audit["pairs"]),
        )

    database_summary, marker_rows = load_database_evidence(
        database_dir / f"{INDEX_NAME}.pkl",
        profiles,
        mapout_audits,
    )
    checks.add(
        "database",
        "database-taxonomy-entries",
        database_summary["taxonomy_entries"] == EXPECTED_TAXONOMY_ENTRIES,
        str(database_summary["taxonomy_entries"]),
    )
    checks.add(
        "database",
        "database-marker-entries",
        database_summary["marker_entries"] == EXPECTED_MARKER_ENTRIES,
        str(database_summary["marker_entries"]),
    )
    checks.add(
        "database",
        "database-known-and-unknown-sgbs",
        database_summary["known_sgbs"] == EXPECTED_KNOWN_SGBS
        and database_summary["unknown_sgbs"] == EXPECTED_UNKNOWN_SGBS,
        f"known={database_summary['known_sgbs']};unknown={database_summary['unknown_sgbs']}",
    )

    species_rows, sgb_rows = build_species_and_sgb_tables(
        full_profile, EXPECTED_INPUT_PAIRS
    )
    checks.add(
        "sgb",
        "detected-species-positive",
        sum(row["Rank"] == "Species" for row in species_rows) > 0,
        str(sum(row["Rank"] == "Species" for row in species_rows)),
    )
    checks.add(
        "sgb",
        "detected-sgbs-positive",
        len(sgb_rows) > 0,
        str(len(sgb_rows)),
    )
    write_tsv(frozen_dir / "species-profile.tsv", species_rows, PROFILE_FIELDS)
    write_tsv(frozen_dir / "sgb-profile.tsv", sgb_rows, PROFILE_FIELDS)
    write_tsv(frozen_dir / "marker-support.tsv", marker_rows, MARKER_FIELDS)

    marker_valid = all(
        int(row["TotalDatabaseMarkers"]) > 0
        and 0 <= int(row["NonzeroMarkers"]) <= int(row["TotalDatabaseMarkers"])
        and 0.0 <= float(row["MarkerSupportFraction"]) <= 1.0
        and row["PercNonzeroRole"]
        == "quasi_marker_disambiguation_not_detection_limit"
        and int(row["ViralMarkerNamesInMapout"]) >= 0
        and int(row["ViralMarkerHitsInMapout"]) >= 0
        and int(row["MetadataExcludedSGBMarkerNamesInMapout"]) >= 0
        and int(row["MetadataExcludedSGBMarkerHitsInMapout"]) >= 0
        and int(row["UnexpectedMarkerNamesInMapout"]) == 0
        for row in marker_rows
    )
    checks.add(
        "marker",
        "marker-support-valid",
        bool(marker_rows) and marker_valid,
        f"rows={len(marker_rows)}",
    )
    non_profile_marker_audit: dict[str, dict[str, Any]] = {}
    for pair_depth in sorted(profiles):
        row = next(
            item
            for item in marker_rows
            if int(item["PairDepth"]) == pair_depth
        )
        non_profile_marker_audit[str(pair_depth)] = {
            "distinct_vdb_marker_names": int(row["ViralMarkerNamesInMapout"]),
            "vdb_marker_hits": int(row["ViralMarkerHitsInMapout"]),
            "distinct_metadata_excluded_sgb_marker_names": int(
                row["MetadataExcludedSGBMarkerNamesInMapout"]
            ),
            "metadata_excluded_sgb_marker_hits": int(
                row["MetadataExcludedSGBMarkerHitsInMapout"]
            ),
            "metadata_excluded_sgb_labels": [
                label
                for label in row["MetadataExcludedSGBLabelsInMapout"].split(";")
                if label
            ],
            "unexpected_marker_names": int(
                row["UnexpectedMarkerNamesInMapout"]
            ),
        }
    non_profile_marker_valid = all(
        evidence["distinct_vdb_marker_names"] > 0
        and evidence["vdb_marker_hits"] >= evidence["distinct_vdb_marker_names"]
        and evidence["distinct_metadata_excluded_sgb_marker_names"] > 0
        and evidence["metadata_excluded_sgb_marker_hits"]
        >= evidence["distinct_metadata_excluded_sgb_marker_names"]
        and evidence["metadata_excluded_sgb_labels"]
        and evidence["unexpected_marker_names"] == 0
        for evidence in non_profile_marker_audit.values()
    )
    checks.add(
        "marker",
        "non-pkl-sgb-markers-accounted-separately",
        len(non_profile_marker_audit) == len(profiles)
        and non_profile_marker_valid,
        json.dumps(non_profile_marker_audit, sort_keys=True),
    )

    threshold_rows = build_threshold_sensitivity(species_rows, sgb_rows)
    depth_rows = build_depth_sensitivity(profiles, subsample_audits)
    write_tsv(
        frozen_dir / "threshold-sensitivity.tsv",
        threshold_rows,
        THRESHOLD_FIELDS,
    )
    write_tsv(
        frozen_dir / "depth-sensitivity.tsv",
        depth_rows,
        DEPTH_FIELDS,
    )
    checks.add(
        "sensitivity",
        "threshold-grid",
        {(row["Rank"], row["ThresholdPct"]) for row in threshold_rows}
        == {
            (rank, threshold)
            for rank in ("Species", "SGB")
            for threshold in (0.0, 0.01, 0.1, 1.0)
        },
        str(len(threshold_rows)),
    )
    checks.add(
        "sensitivity",
        "depth-grid",
        {int(row["PairDepth"]) for row in depth_rows}
        == {20_000, 50_000, EXPECTED_INPUT_PAIRS},
        ",".join(str(row["PairDepth"]) for row in depth_rows),
    )

    resource_files = {
        "full": frozen_dir / "logs/metaphlan-full.resources.txt",
        "classified_only": (
            frozen_dir / "logs/metaphlan-classified-only.resources.txt"
        ),
        "depth_20000": frozen_dir / "logs/metaphlan-depth-20000.resources.txt",
        "depth_50000": frozen_dir / "logs/metaphlan-depth-50000.resources.txt",
        "extract_metadata": frozen_dir / "logs/extract-metadata.resources.txt",
        "extract_bowtie2": frozen_dir / "logs/extract-bowtie2.resources.txt",
    }
    resources = {
        label: resource_summary(path)
        for label, path in resource_files.items()
        if path.exists()
    }

    unclassified = unclassified_row(full_profile)
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "evidence_date": "2026-07-21",
        "source_project": "PRJEB52977",
        "source_sample": "SAMEA14435832",
        "run_accession": "ERR9765746",
        "input_pairs": EXPECTED_INPUT_PAIRS,
        "input_reads": EXPECTED_INPUT_READS,
        "profiling_reads": full_profile["processed_reads"],
        "short_reads_excluded": EXPECTED_INPUT_READS
        - full_profile["processed_reads"],
        "read_min_len": 70,
        "metaphlan_version": versions["MetaPhlAn"],
        "bowtie2_version": versions["Bowtie2"],
        "python_version": versions["Python"],
        "database_release": INDEX_NAME,
        "database_metadata_md5": EXPECTED_ARCHIVES["marker_metadata"]["md5"],
        "database_bowtie2_md5": EXPECTED_ARCHIVES["bowtie2_index"]["md5"],
        "archive_audits": archive_audits,
        "database_files": database_file_summary,
        "database_catalog": database_summary,
        "estimated_reads_mapped_to_known_clades": full_profile[
            "estimated_mapped_reads"
        ],
        "unclassified_pct": float(unclassified["relative_abundance"])
        if unclassified
        else None,
        "detected_species": sum(
            row["Rank"] == "Species" for row in species_rows
        ),
        "detected_sgbs": len(sgb_rows),
        "known_sgbs_detected": sum(
            row["SGBClass"] == "Known SGB" for row in sgb_rows
        ),
        "unknown_sgbs_detected": sum(
            row["SGBClass"] == "Unknown SGB" for row in sgb_rows
        ),
        "default_unclassified_estimation": True,
        "classified_only_profile_is_separate": True,
        "perc_nonzero": 0.33,
        "perc_nonzero_role": "quasi_marker_disambiguation_not_detection_limit",
        "subsampling_seed": SUBSAMPLING_SEED,
        "subsampling_pair_depths": [20_000, 50_000],
        "full_depth_uses_all_reads_without_subsampling": True,
        "clean_fastq_audit": clean_audit,
        "subsample_audits": {
            str(depth): audit
            for depth, audit in subsample_audits.items()
            if depth != EXPECTED_INPUT_PAIRS
        },
        "mapout_audits": {
            str(depth): {
                key: value
                for key, value in audit.items()
                if key not in {"path", "marker_hits"}
            }
            for depth, audit in mapout_audits.items()
        },
        "non_profile_marker_audit": non_profile_marker_audit,
        "resource_summary": resources,
        "environment_yaml_sha256": EXPECTED_ENV_YAML_SHA256,
        "environment_lock_sha256": EXPECTED_ENV_LOCK_SHA256,
        "environment_lock_packages": EXPECTED_ENV_PACKAGES,
        "environment_lock_targets": lock_targets,
        "raw_fastq_committed": False,
        "database_archives_committed": False,
        "database_files_committed": False,
        "mapout_committed": False,
        "one_time_network_access": True,
        "qa_network_access": False,
        "initialization_checks_passed": checks.passed,
        "initialization_checks_failed": checks.failed,
        "checksum_failures": 0,
    }
    write_json(frozen_dir / "run-summary.json", summary)
    write_tsv(
        frozen_dir / "initialization-audit.tsv",
        checks.rows,
        ["Category", "CheckID", "Status", "Detail"],
    )

    normalization = normalize_frozen_paths(
        frozen_dir,
        [
            (metadata_archive, "${METAPHLAN_METADATA_ARCHIVE}"),
            (bowtie2_archive, "${METAPHLAN_BOWTIE2_ARCHIVE}"),
            (database_dir, "${METAPHLAN_DB_DIR}"),
            (work_dir, "${WORK_DIR}"),
            (raw_dir, "${RAW_DIR}"),
            (environment_prefix, "${BIOBAKERY_ENV_PREFIX}"),
            (project_root, "${PROJECT_ROOT}"),
            (Path.home(), "${HOME}"),
        ],
    )
    hits = host_path_hits(
        frozen_dir,
        [str(project_root), str(environment_prefix), str(raw_dir), str(database_dir)],
    )
    if hits:
        raise RuntimeError(f"Host-specific paths remain in frozen evidence: {hits[:5]}")

    summary_path = frozen_dir / "run-summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["path_normalization"] = {
        **normalization,
        "host_specific_paths_retained": False,
        "remaining_hits": 0,
    }
    write_json(summary_path, summary)
    write_frozen_checksums(frozen_dir)
    _, checksum_failures = verify_frozen_checksums(frozen_dir)
    if checksum_failures:
        raise RuntimeError(f"Frozen checksum failures: {checksum_failures}")
    if checks.failed:
        raise RuntimeError(
            f"Article 15 initialization failed {checks.failed} checks; see initialization-audit.tsv"
        )


def publication_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = figure_dir / f"{stem}.pdf"
    png_path = figure_dir / f"{stem}.png"
    tiff_path = figure_dir / f"{stem}.tiff"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=350, bbox_inches="tight", facecolor="white")
    with Image.open(png_path) as image:
        image.save(
            tiff_path,
            format="TIFF",
            compression="tiff_lzw",
            dpi=(350, 350),
        )
    plt.close(fig)


def plot_composition(
    species_rows: list[dict[str, str]], figure_dir: Path
) -> None:
    classified = [row for row in species_rows if row["Rank"] == "Species"]
    unclassified = next(
        row for row in species_rows if row["Rank"] == "Unclassified"
    )
    classified.sort(
        key=lambda row: float(row["RelativeAbundancePct"]), reverse=True
    )
    top = classified[:10]
    other = sum(float(row["RelativeAbundancePct"]) for row in classified[10:])
    labels = [display_taxon(row["SpeciesLabel"]) for row in top]
    values = [float(row["RelativeAbundancePct"]) for row in top]
    colors = ["#2A6F97"] * len(top)
    if other > 0:
        labels.append("Other classified")
        values.append(other)
        colors.append("#8D99AE")
    labels.append("Unclassified")
    values.append(float(unclassified["RelativeAbundancePct"]))
    colors.append("#D9D9D9")

    order = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(8.2, max(4.6, 0.39 * len(labels))))
    ax.barh(order, values, color=colors, edgecolor="white", linewidth=0.6)
    ax.set_yticks(order, labels)
    ax.set_xlabel("Relative abundance (%)")
    ax.set_title("Marker-based species-level composition of MOCK1")
    ax.set_xlim(0, max(values) * 1.18 if values else 1)
    for y, value in zip(order, values):
        ax.text(value + max(values) * 0.012, y, f"{value:.3f}", va="center", fontsize=8)
    ax.text(
        1.0,
        -0.12,
        "Default MetaPhlAn denominator includes estimated unclassified reads",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color="#555555",
        fontsize=8,
    )
    save_figure(fig, figure_dir, "15-metaphlan-composition")


def plot_marker_support(
    marker_rows: list[dict[str, str]],
    database_catalog: dict[str, Any],
    figure_dir: Path,
) -> None:
    full_rows = [
        row
        for row in marker_rows
        if int(row["PairDepth"]) == EXPECTED_INPUT_PAIRS
    ]
    palette = {"Known SGB": "#2A6F97", "Unknown SGB": "#E76F51"}
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.35), gridspec_kw={"width_ratios": [0.8, 1.6]})

    catalog_labels = ["Known SGB", "Unknown SGB"]
    catalog_values = [
        int(database_catalog["known_sgbs"]),
        int(database_catalog["unknown_sgbs"]),
    ]
    axes[0].bar(
        catalog_labels,
        catalog_values,
        color=[palette[label] for label in catalog_labels],
        width=0.62,
    )
    axes[0].set_ylabel("SGBs in database")
    axes[0].set_title("Database catalog")
    axes[0].tick_params(axis="x", rotation=20)
    for index, value in enumerate(catalog_values):
        axes[0].text(index, value, f"{value:,}", ha="center", va="bottom", fontsize=8)

    for group in catalog_labels:
        rows = [row for row in full_rows if row["SGBClass"] == group]
        if not rows:
            continue
        x = [max(float(row["RelativeAbundancePct"]), 1e-7) for row in rows]
        y = [float(row["MarkerSupportFraction"]) for row in rows]
        sizes = [28 + 16 * math.log10(1 + int(row["MappedMarkerHits"])) for row in rows]
        axes[1].scatter(
            x,
            y,
            s=sizes,
            alpha=0.78,
            color=palette[group],
            edgecolor="white",
            linewidth=0.45,
            label=group,
        )
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Relative abundance (%)")
    axes[1].set_ylabel("Raw nonzero marker fraction")
    axes[1].set_ylim(-0.02, 1.04)
    axes[1].set_title("Detected SGB evidence")
    axes[1].legend(loc="lower right")
    label_rows = sorted(
        full_rows,
        key=lambda row: float(row["RelativeAbundancePct"]),
        reverse=True,
    )[:5]
    label_rows.sort(
        key=lambda row: float(row["MarkerSupportFraction"]), reverse=True
    )
    label_positions = np.linspace(0.985, 0.745, len(label_rows))
    for row, label_y in zip(label_rows, label_positions):
        axes[1].annotate(
            display_taxon(row["SpeciesLabel"]),
            (
                max(float(row["RelativeAbundancePct"]), 1e-7),
                float(row["MarkerSupportFraction"]),
            ),
            xycoords="data",
            xytext=(0.985, label_y),
            textcoords="axes fraction",
            ha="right",
            va="center",
            fontsize=6.4,
            arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.55},
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.4},
        )
    fig.text(
        0.99,
        0.01,
        "Marker fraction is descriptive; --perc_nonzero 0.33 controls quasi-marker disambiguation",
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#555555",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save_figure(fig, figure_dir, "15-sgb-marker-support")


def plot_sensitivity(
    depth_rows: list[dict[str, str]],
    threshold_rows: list[dict[str, str]],
    figure_dir: Path,
) -> None:
    depth_rows = sorted(depth_rows, key=lambda row: int(row["PairDepth"]))
    colors = {"Species": "#2A6F97", "SGB": "#E76F51"}
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.25))
    depths = [int(row["PairDepth"]) for row in depth_rows]
    for rank, field in (("Species", "DetectedSpecies"), ("SGB", "DetectedSGBs")):
        axes[0].plot(
            depths,
            [int(row[field]) for row in depth_rows],
            marker="o" if rank == "Species" else "s",
            linestyle="-" if rank == "Species" else "--",
            linewidth=4.2 if rank == "Species" else 2.0,
            markersize=6.8 if rank == "Species" else 5.2,
            color=colors[rank],
            label=rank,
        )
    axes[0].set_xlabel("Input read pairs")
    axes[0].set_ylabel("Detected clades")
    axes[0].set_title("Paired-depth sensitivity")
    axes[0].legend()
    axes[0].ticklabel_format(style="plain", axis="x")
    if all(
        int(row["DetectedSpecies"]) == int(row["DetectedSGBs"])
        for row in depth_rows
    ):
        axes[0].text(
            0.98,
            0.04,
            "Species and SGB counts coincide",
            transform=axes[0].transAxes,
            ha="right",
            va="bottom",
            fontsize=7.5,
            color="#555555",
        )
    for row in depth_rows:
        axes[0].text(
            int(row["PairDepth"]),
            int(row["DetectedSGBs"]),
            f"  {float(row['UnclassifiedPct']):.1f}% unclassified",
            fontsize=7,
            va="bottom",
            color="#555555",
        )

    threshold_labels = ["0", "0.01", "0.1", "1"]
    x = np.arange(len(threshold_labels))
    threshold_by_rank: dict[str, list[dict[str, str]]] = {}
    for rank in ("Species", "SGB"):
        rows = sorted(
            (row for row in threshold_rows if row["Rank"] == rank),
            key=lambda row: float(row["ThresholdPct"]),
        )
        threshold_by_rank[rank] = rows
        axes[1].plot(
            x,
            [100 * float(row["RetainedFractionOfClassified"]) for row in rows],
            marker="o" if rank == "Species" else "s",
            linestyle="-" if rank == "Species" else "--",
            linewidth=4.2 if rank == "Species" else 2.0,
            markersize=6.8 if rank == "Species" else 5.2,
            color=colors[rank],
            label=rank,
        )
    curves_coincide = all(
        species["DetectedClades"] == sgb["DetectedClades"]
        and math.isclose(
            float(species["RetainedFractionOfClassified"]),
            float(sgb["RetainedFractionOfClassified"]),
            abs_tol=1e-12,
        )
        for species, sgb in zip(
            threshold_by_rank["Species"], threshold_by_rank["SGB"]
        )
    )
    if curves_coincide:
        for xpos, row in zip(x, threshold_by_rank["SGB"]):
            axes[1].text(
                xpos,
                100 * float(row["RetainedFractionOfClassified"]) - 2.3,
                str(row["DetectedClades"]),
                ha="center",
                va="top",
                fontsize=7,
                color="#555555",
            )
    else:
        for rank, vertical_offset in (("Species", 1.2), ("SGB", -2.3)):
            for xpos, row in zip(x, threshold_by_rank[rank]):
                axes[1].text(
                    xpos,
                    100 * float(row["RetainedFractionOfClassified"])
                    + vertical_offset,
                    str(row["DetectedClades"]),
                    ha="center",
                    va="bottom" if vertical_offset > 0 else "top",
                    fontsize=7,
                    color=colors[rank],
                )
    axes[1].set_xticks(x, threshold_labels)
    axes[1].set_xlabel("Post-profile abundance threshold (%)")
    axes[1].set_ylabel("Classified abundance retained (%)")
    axes[1].set_ylim(0, 103)
    axes[1].set_title("Detection versus quantification")
    axes[1].legend(loc="lower left")
    axes[1].text(
        0.98,
        0.52,
        (
            "Species and SGB curves coincide; numbers show retained clades"
            if curves_coincide
            else "Numbers show retained clades"
        ),
        transform=axes[1].transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#555555",
    )
    fig.tight_layout()
    save_figure(fig, figure_dir, "15-detection-quantification-sensitivity")


def routine_validate(args: argparse.Namespace) -> None:
    if args.output_dir is None or args.figure_dir is None:
        raise ValueError("Routine validation requires --output-dir and --figure-dir")
    project_root = args.project_root.resolve()
    environment_prefix = args.environment_prefix.resolve()
    frozen_dir = (
        args.frozen_dir.resolve()
        if args.frozen_dir
        else project_root / "data/small/15-metaphlan-frozen"
    )
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    source_manifest_checks(project_root, checks)
    database_manifest_checks(project_root, checks)
    versions, _, _ = environment_checks(project_root, environment_prefix, checks)
    notice_path = project_root / "data/small/15-data-NOTICE.txt"
    checks.add(
        "source",
        "notice-sha256",
        hash_file(notice_path) == EXPECTED_NOTICE_SHA256,
        hash_file(notice_path),
    )

    required_frozen = [
        "run-summary.json",
        "tool-versions.tsv",
        "commands.sh",
        "profile-all.tsv",
        "profile-classified-only.tsv",
        "species-profile.tsv",
        "sgb-profile.tsv",
        "marker-support.tsv",
        "depth-sensitivity.tsv",
        "threshold-sensitivity.tsv",
        "database-files.sha256",
        "file-checksums.sha256",
    ]
    for relative in required_frozen:
        checks.add(
            "profile",
            f"frozen-{relative.replace('/', '-').replace('.', '-')}",
            (frozen_dir / relative).is_file(),
            relative,
        )

    checksum_entries, checksum_failures = verify_frozen_checksums(frozen_dir)
    checks.add(
        "profile",
        "frozen-checksums",
        not checksum_failures,
        f"entries={checksum_entries};failures={len(checksum_failures)}",
    )
    run_summary = json.loads(
        (frozen_dir / "run-summary.json").read_text(encoding="utf-8")
    )
    checks.add(
        "profile",
        "frozen-run-status",
        run_summary.get("status") == "passed",
        str(run_summary.get("status")),
    )
    checks.add(
        "source",
        "frozen-input-pairs",
        run_summary.get("input_pairs") == EXPECTED_INPUT_PAIRS,
        str(run_summary.get("input_pairs")),
    )
    checks.add(
        "source",
        "frozen-input-reads",
        run_summary.get("input_reads") == EXPECTED_INPUT_READS,
        str(run_summary.get("input_reads")),
    )
    checks.add(
        "profile",
        "frozen-profiling-reads",
        run_summary.get("profiling_reads") == EXPECTED_PROFILING_READS,
        str(run_summary.get("profiling_reads")),
    )
    checks.add(
        "profile",
        "frozen-short-read-ledger",
        run_summary.get("short_reads_excluded") == EXPECTED_SHORT_READS,
        str(run_summary.get("short_reads_excluded")),
    )

    full_profile = parse_profile(frozen_dir / "profile-all.tsv")
    classified_only = parse_profile(
        frozen_dir / "profile-classified-only.tsv"
    )
    profile_checks(
        {EXPECTED_INPUT_PAIRS: full_profile}, classified_only, checks
    )
    species_rows = read_tsv(frozen_dir / "species-profile.tsv")
    sgb_rows = read_tsv(frozen_dir / "sgb-profile.tsv")
    marker_rows = read_tsv(frozen_dir / "marker-support.tsv")
    depth_rows = read_tsv(frozen_dir / "depth-sensitivity.tsv")
    threshold_rows = read_tsv(frozen_dir / "threshold-sensitivity.tsv")

    detected_species = sum(row["Rank"] == "Species" for row in species_rows)
    detected_sgbs = len(sgb_rows)
    checks.add("sgb", "species-table-positive", detected_species > 0, str(detected_species))
    checks.add("sgb", "sgb-table-positive", detected_sgbs > 0, str(detected_sgbs))
    checks.add(
        "sgb",
        "sgb-class-values",
        {row["SGBClass"] for row in sgb_rows}
        <= {"Known SGB", "Unknown SGB"},
        ",".join(sorted({row["SGBClass"] for row in sgb_rows})),
    )
    checks.add(
        "marker",
        "marker-rows-cover-all-depths",
        {int(row["PairDepth"]) for row in marker_rows}
        == {20_000, 50_000, EXPECTED_INPUT_PAIRS},
        str(len(marker_rows)),
    )
    marker_valid = all(
        int(row["TotalDatabaseMarkers"]) > 0
        and 0 <= int(row["NonzeroMarkers"]) <= int(row["TotalDatabaseMarkers"])
        and 0 <= float(row["MarkerSupportFraction"]) <= 1
        and row["PercNonzeroRole"]
        == "quasi_marker_disambiguation_not_detection_limit"
        and int(row["ViralMarkerNamesInMapout"]) >= 0
        and int(row["ViralMarkerHitsInMapout"]) >= 0
        and int(row["MetadataExcludedSGBMarkerNamesInMapout"]) >= 0
        and int(row["MetadataExcludedSGBMarkerHitsInMapout"]) >= 0
        and int(row["UnexpectedMarkerNamesInMapout"]) == 0
        for row in marker_rows
    )
    checks.add("marker", "marker-table-valid", marker_valid, str(len(marker_rows)))
    non_profile_marker_audit = run_summary.get("non_profile_marker_audit", {})
    checks.add(
        "marker",
        "non-pkl-sgb-markers-accounted-separately",
        set(non_profile_marker_audit) == {"20000", "50000", "99991"}
        and all(
            int(evidence["distinct_vdb_marker_names"]) > 0
            and int(evidence["vdb_marker_hits"])
            >= int(evidence["distinct_vdb_marker_names"])
            and int(evidence["distinct_metadata_excluded_sgb_marker_names"])
            > 0
            and int(evidence["metadata_excluded_sgb_marker_hits"])
            >= int(evidence["distinct_metadata_excluded_sgb_marker_names"])
            and bool(evidence["metadata_excluded_sgb_labels"])
            and int(evidence["unexpected_marker_names"]) == 0
            for evidence in non_profile_marker_audit.values()
        ),
        json.dumps(non_profile_marker_audit, sort_keys=True),
    )
    checks.add(
        "sensitivity",
        "depth-grid",
        {int(row["PairDepth"]) for row in depth_rows}
        == {20_000, 50_000, EXPECTED_INPUT_PAIRS},
        ",".join(row["PairDepth"] for row in depth_rows),
    )
    checks.add(
        "sensitivity",
        "subsampling-seed",
        all(
            row["SubsamplingSeed"] == str(SUBSAMPLING_SEED)
            for row in depth_rows
            if int(row["PairDepth"]) < EXPECTED_INPUT_PAIRS
        ),
        str(SUBSAMPLING_SEED),
    )
    checks.add(
        "sensitivity",
        "threshold-grid",
        {(row["Rank"], float(row["ThresholdPct"])) for row in threshold_rows}
        == {
            (rank, threshold)
            for rank in ("Species", "SGB")
            for threshold in (0.0, 0.01, 0.1, 1.0)
        },
        str(len(threshold_rows)),
    )

    database_files = parse_database_file_manifest(
        frozen_dir / "database-files.sha256"
    )
    checks.add(
        "database",
        "frozen-database-pkl",
        any(row["Path"].endswith(f"/{INDEX_NAME}.pkl") for row in database_files),
        str(len(database_files)),
    )
    checks.add(
        "database",
        "frozen-database-six-index-files",
        sum(".bt2" in row["Path"] for row in database_files) == 6,
        str(sum(".bt2" in row["Path"] for row in database_files)),
    )
    database_catalog = run_summary["database_catalog"]
    checks.add(
        "database",
        "frozen-database-catalog",
        database_catalog["taxonomy_entries"] == EXPECTED_TAXONOMY_ENTRIES
        and database_catalog["marker_entries"] == EXPECTED_MARKER_ENTRIES
        and database_catalog["known_sgbs"] == EXPECTED_KNOWN_SGBS
        and database_catalog["unknown_sgbs"] == EXPECTED_UNKNOWN_SGBS,
        (
            f"taxonomy={database_catalog['taxonomy_entries']};"
            f"markers={database_catalog['marker_entries']};"
            f"known={database_catalog['known_sgbs']};"
            f"unknown={database_catalog['unknown_sgbs']}"
        ),
    )
    path_hits = host_path_hits(
        frozen_dir,
        [str(project_root.parent), str(environment_prefix), str(Path.home())],
    )
    checks.add(
        "profile",
        "frozen-path-normalization",
        not path_hits,
        f"hits={len(path_hits)}",
    )

    publication_style()
    plot_composition(species_rows, figure_dir)
    plot_marker_support(marker_rows, database_catalog, figure_dir)
    plot_sensitivity(depth_rows, threshold_rows, figure_dir)
    for stem in FIGURE_STEMS:
        for suffix in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            checks.add(
                "figure",
                f"figure-{stem}-{suffix}",
                path.is_file() and path.stat().st_size > 0,
                str(path.stat().st_size if path.exists() else "missing"),
            )

    audit_files = {
        "tool-audit.tsv": "tool",
        "source-audit.tsv": "source",
        "database-audit.tsv": "database",
        "profile-audit.tsv": "profile",
        "sgb-audit.tsv": "sgb",
        "marker-audit.tsv": "marker",
        "sensitivity-audit.tsv": "sensitivity",
    }
    for filename, category in audit_files.items():
        write_tsv(
            output_dir / filename,
            checks.category(category),
            ["Category", "CheckID", "Status", "Detail"],
        )
    log_lines = [
        f"{row['Status']}\t{row['Category']}\t{row['CheckID']}\t{row['Detail']}\n"
        for row in checks.rows
    ]
    (output_dir / "validation.log").write_text(
        "".join(log_lines), encoding="utf-8"
    )

    unclassified = next(
        row for row in species_rows if row["Rank"] == "Unclassified"
    )
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "evidence_date": run_summary["evidence_date"],
        "input_pairs": EXPECTED_INPUT_PAIRS,
        "input_reads": EXPECTED_INPUT_READS,
        "profiling_reads": EXPECTED_PROFILING_READS,
        "short_reads_excluded": EXPECTED_SHORT_READS,
        "metaphlan_version": versions["MetaPhlAn"],
        "bowtie2_version": versions["Bowtie2"],
        "database_release": INDEX_NAME,
        "database_metadata_md5": EXPECTED_ARCHIVES["marker_metadata"]["md5"],
        "database_bowtie2_md5": EXPECTED_ARCHIVES["bowtie2_index"]["md5"],
        "detected_species": detected_species,
        "detected_sgbs": detected_sgbs,
        "known_sgbs_detected": sum(
            row["SGBClass"] == "Known SGB" for row in sgb_rows
        ),
        "unknown_sgbs_detected": sum(
            row["SGBClass"] == "Unknown SGB" for row in sgb_rows
        ),
        "default_unclassified_estimation": True,
        "unclassified_pct": float(unclassified["RelativeAbundancePct"]),
        "perc_nonzero": 0.33,
        "perc_nonzero_role": "quasi_marker_disambiguation_not_detection_limit",
        "subsampling_seed": SUBSAMPLING_SEED,
        "full_depth_uses_all_reads_without_subsampling": True,
        "checksum_entries": checksum_entries,
        "checksum_failures": len(checksum_failures),
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
    }
    write_json(output_dir / "validation-summary.json", summary)
    if checks.failed:
        raise SystemExit(1)


def main() -> None:
    args = parse_args()
    selected_modes = sum(
        (args.preflight_only, args.initialize_frozen)
    )
    if selected_modes > 1:
        raise ValueError("Choose only one of --preflight-only and --initialize-frozen")
    if args.preflight_only:
        run_preflight(args)
    elif args.initialize_frozen:
        initialize_frozen(args)
    else:
        routine_validate(args)


if __name__ == "__main__":
    main()
