#!/usr/bin/env python3
"""Validate Article 16 Kraken2/Bracken profiling and frozen evidence.

Initialization mode is used once with ignored FASTQ, the checksum-locked
Standard-8 database, and per-fragment Kraken output. Routine mode is fully
offline: it audits the frozen aggregate bundle and regenerates publication
tables and figures for the tutorial QA workflow.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article16-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


EXPECTED_INPUT_PAIRS = 99_991
EXPECTED_INPUT_READS = 199_982
EXPECTED_PAIR_ID_SHA256 = (
    "457cef6e9d603790dfbc26b716b0498169b54c31bc903d067d449d8dcc86770d"
)
EXPECTED_CLEAN = {
    "R1": {
        "records": 99_991,
        "bases": 14_974_589,
        "bytes": 8_661_319,
        "sha256": "ce438de2814a6e605cc24b00e53b697d018f325bc2a9694cece5be158ff32101",
    },
    "R2": {
        "records": 99_991,
        "bases": 14_835_184,
        "bytes": 10_045_722,
        "sha256": "207d080cd5dc98bd10fd4af8c492fb3f2000f73e1458a9d1d77e57347f4fe459",
    },
}
EXPECTED_ARCHIVE_BYTES = 5_946_578_575
EXPECTED_ARCHIVE_MD5 = "7685f43cce057c2ca18511c925399b72"
EXPECTED_ARCHIVE_SHA256 = (
    "b17d05ca1459564b49b63d014c4b2ee6ebe1ca0143e6c184e0dfd6d940a55981"
)
EXPECTED_VERSIONS = {
    "Kraken2": "2.17.1",
    "BrackenCLI": "3.0.1",
    "BrackenPackage": "3.1p1",
    "Python": "3.12.13",
}
EXPECTED_ENV_YAML_SHA256 = (
    "51e81b16c47386ac58982445f686466b193999269ab9aef42342e15758a0c9af"
)
EXPECTED_ENV_LOCK_SHA256 = (
    "7ea37cd9e01dee9d6091ca42c0da7e514b637090be9e3f3362e98a0ef8b0789e"
)
EXPECTED_ENV_PACKAGES = 187
EXPECTED_SOURCE_MANIFEST_SHA256 = (
    "8bd3fddad0754056104d91885b16d67dc51c9bf6c46b9837a71df5f70653766e"
)
EXPECTED_DATABASE_MANIFEST_SHA256 = (
    "b880bf02c890dbb0391035403a8b658d052025e723b457feafde5d7192d55dc6"
)
EXPECTED_DATABASE_FILES_MD5_SHA256 = (
    "8691836fe757e975828ce709d6ff0cee668102f76421d27235a762adf3844ba7"
)
EXPECTED_NOTICE_SHA256 = (
    "8813ed349a70b4045971f33264d96b35be420a409118446a9ab43b5cb45b14bb"
)
EXPECTED_ARTICLE13_SUMMARY_SHA256 = (
    "f6dcc51b6535247de7f370dc2334994dd85dca61d3f86d252294127faa3460fe"
)
EXPECTED_LOCK_PREFIXES = {
    "kraken2": "kraken2-2.17.1-",
    "bracken": "bracken-3.1p1-",
    "python": "python-3.12.13-",
    "matplotlib-base": "matplotlib-base-3.10.5-",
    "pillow": "pillow-12.3.0-",
    "pyyaml": "pyyaml-6.0.3-",
}
FIGURE_STEMS = (
    "16-kraken-classification-ledger",
    "16-bracken-redistribution",
    "16-bracken-parameter-sensitivity",
)
TEXT_SUFFIXES = {".csv", ".json", ".log", ".md5", ".sh", ".tsv", ".txt"}
RANK_LABELS = {
    "U": "Unclassified",
    "R": "Root",
    "D": "Domain",
    "K": "Kingdom",
    "P": "Phylum",
    "C": "Class",
    "O": "Order",
    "F": "Family",
    "G": "Genus",
    "S": "Species",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--environment-prefix", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path)
    parser.add_argument("--raw-dir", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--database-archive", type=Path)
    parser.add_argument("--database-dir", type=Path)
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


def hash_file_multi(path: Path) -> dict[str, str | int]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return {"bytes": path.stat().st_size, "md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict[str, str]] = []

    def add(self, category: str, check_id: str, passed: bool, detail: Any) -> None:
        self.rows.append(
            {
                "Category": category,
                "CheckID": check_id,
                "Status": "PASS" if passed else "FAIL",
                "Detail": str(detail),
            }
        )

    @property
    def passed(self) -> int:
        return sum(row["Status"] == "PASS" for row in self.rows)

    @property
    def failed(self) -> int:
        return sum(row["Status"] == "FAIL" for row in self.rows)


def static_contract_checks(project_root: Path, checks: Checks) -> None:
    expected = {
        "env/kraken.yml": EXPECTED_ENV_YAML_SHA256,
        "env/kraken-linux-64.lock": EXPECTED_ENV_LOCK_SHA256,
        "data/small/16-source-manifest.tsv": EXPECTED_SOURCE_MANIFEST_SHA256,
        "data/small/16-database-manifest.tsv": EXPECTED_DATABASE_MANIFEST_SHA256,
        "data/small/16-standard8-files.md5": EXPECTED_DATABASE_FILES_MD5_SHA256,
        "data/small/16-data-NOTICE.txt": EXPECTED_NOTICE_SHA256,
        "data/small/13-qc-frozen/run-summary.json": EXPECTED_ARTICLE13_SUMMARY_SHA256,
    }
    for relative, digest in expected.items():
        path = project_root / relative
        observed = hash_file(path) if path.is_file() else "missing"
        checks.add("contract", f"sha256-{Path(relative).name}", observed == digest, observed)

    source_rows = read_tsv(project_root / "data/small/16-source-manifest.tsv")
    checks.add("source", "source-two-mates", len(source_rows) == 2, len(source_rows))
    checks.add(
        "source",
        "source-identifiers",
        {row["RunAccession"] for row in source_rows} == {"ERR9765746"}
        and {row["ProjectAccession"] for row in source_rows} == {"PRJEB52977"}
        and {row["SampleAccession"] for row in source_rows} == {"SAMEA14435832"},
        "ERR9765746/PRJEB52977/SAMEA14435832",
    )
    checks.add(
        "source",
        "source-paired-ledger",
        all(int(row["Records"]) == EXPECTED_INPUT_PAIRS for row in source_rows)
        and sum(int(row["Records"]) for row in source_rows) == EXPECTED_INPUT_READS,
        sum(int(row["Records"]) for row in source_rows),
    )

    database_rows = read_tsv(project_root / "data/small/16-database-manifest.tsv")
    checks.add("database", "one-database-row", len(database_rows) == 1, len(database_rows))
    if database_rows:
        row = database_rows[0]
        checks.add(
            "database",
            "database-release-fixed",
            row["release_id"] == "Standard-8-20260626"
            and "/latest/" not in row["archive_url"],
            row["release_id"],
        )
        checks.add(
            "database",
            "database-archive-contract",
            row["expected_checksum"] == EXPECTED_ARCHIVE_MD5
            and int(row["expected_compressed_bytes"]) == EXPECTED_ARCHIVE_BYTES,
            f"{row['expected_checksum']};{row['expected_compressed_bytes']}",
        )
        checks.add(
            "database",
            "database-cap-disclosed",
            "capped at 8 GB" in row["notes"]
            and "sensitivity/accuracy" in row["notes"],
            row["notes"],
        )


def environment_checks(project_root: Path, prefix: Path, checks: Checks) -> dict[str, str]:
    lock_lines = (project_root / "env/kraken-linux-64.lock").read_text(encoding="utf-8").splitlines()
    package_lines = [line for line in lock_lines if line.startswith("http")]
    checks.add("environment", "lock-package-count", len(package_lines) == EXPECTED_ENV_PACKAGES, len(package_lines))
    for name, marker in EXPECTED_LOCK_PREFIXES.items():
        matches = [line for line in package_lines if marker in line]
        checks.add("environment", f"lock-{name}", len(matches) == 1, len(matches))

    env = os.environ.copy()
    env["PATH"] = f"{prefix / 'bin'}:{env.get('PATH', '')}"
    kraken_text = subprocess.run(
        [str(prefix / "bin/kraken2"), "--version"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout
    bracken_text = subprocess.run(
        [str(prefix / "bin/bracken"), "-v"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout
    python_text = subprocess.run(
        [str(prefix / "bin/python"), "--version"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout
    versions = {
        "Kraken2": re.search(r"Kraken version ([0-9.]+)", kraken_text).group(1),
        "BrackenCLI": re.search(r"Bracken v([^\s]+)", bracken_text).group(1),
        "Python": re.search(r"Python ([0-9.]+)", python_text).group(1),
    }
    meta_files = list((prefix / "conda-meta").glob("bracken-*.json"))
    versions["BrackenPackage"] = (
        json.loads(meta_files[0].read_text(encoding="utf-8"))["version"]
        if len(meta_files) == 1
        else "missing"
    )
    for key, expected in EXPECTED_VERSIONS.items():
        checks.add("environment", f"version-{key.lower()}", versions[key] == expected, versions[key])
    checks.add(
        "environment",
        "bracken-version-string-discrepancy-recorded",
        versions["BrackenCLI"] == "3.0.1" and versions["BrackenPackage"] == "3.1p1",
        f"cli={versions['BrackenCLI']};package={versions['BrackenPackage']}",
    )
    return versions


def fastq_records(path: Path):
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not plus or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield header.rstrip("\r\n"), sequence


def normalized_id(header: str) -> str:
    value = header[1:].split()[0]
    if value.endswith("/1") or value.endswith("/2"):
        value = value[:-2]
    return value


def audit_fastq_pair(r1: Path, r2: Path) -> dict[str, Any]:
    stats = {
        "R1": {"records": 0, "bases": 0, "bytes": r1.stat().st_size, "sha256": hash_file(r1)},
        "R2": {"records": 0, "bases": 0, "bytes": r2.stat().st_size, "sha256": hash_file(r2)},
    }
    pair_digest = hashlib.sha256()
    iterator1 = fastq_records(r1)
    iterator2 = fastq_records(r2)
    while True:
        record1 = next(iterator1, None)
        record2 = next(iterator2, None)
        if record1 is None or record2 is None:
            if record1 is not None or record2 is not None:
                raise ValueError("Paired FASTQ files have unequal record counts")
            break
        id1, id2 = normalized_id(record1[0]), normalized_id(record2[0])
        if id1 != id2:
            raise ValueError(f"Unsynchronized pair: {id1} != {id2}")
        pair_digest.update(id1.encode("ascii"))
        pair_digest.update(b"\n")
        for mate, record in (("R1", record1), ("R2", record2)):
            stats[mate]["records"] += 1
            stats[mate]["bases"] += len(record[1])
    return {
        "pairs": stats["R1"]["records"],
        "mates_synchronized": True,
        "normalized_pair_id_sha256": pair_digest.hexdigest(),
        "R1": stats["R1"],
        "R2": stats["R2"],
    }


def fastq_checks(audit: dict[str, Any], checks: Checks) -> None:
    checks.add("source", "clean-pair-count", audit["pairs"] == EXPECTED_INPUT_PAIRS, audit["pairs"])
    checks.add(
        "source",
        "clean-pair-id-hash",
        audit["normalized_pair_id_sha256"] == EXPECTED_PAIR_ID_SHA256,
        audit["normalized_pair_id_sha256"],
    )
    for mate in ("R1", "R2"):
        for key, expected in EXPECTED_CLEAN[mate].items():
            checks.add("source", f"clean-{mate.lower()}-{key}", audit[mate][key] == expected, audit[mate][key])


def parse_kraken_report(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            values = line.rstrip("\n").split("\t")
            if len(values) != 6:
                raise ValueError(f"Unexpected Kraken report row: {line[:120]}")
            raw_name = values[5]
            rows.append(
                {
                    "percentage": float(values[0]),
                    "clade_fragments": int(values[1]),
                    "direct_fragments": int(values[2]),
                    "rank": values[3],
                    "taxid": values[4],
                    "name": raw_name.strip(),
                    "indent_spaces": len(raw_name) - len(raw_name.lstrip(" ")),
                }
            )
    return rows


def audit_kraken_output(path: Path) -> dict[str, Any]:
    status = Counter()
    taxids = Counter()
    pair_digest = hashlib.sha256()
    length_pipe_rows = 0
    paired_hit_separator_rows = 0
    rows = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            values = line.rstrip("\n").split("\t")
            if len(values) != 5:
                raise ValueError(f"Unexpected Kraken output row: {line[:120]}")
            call, read_id, taxid, lengths, hit_list = values
            status[call] += 1
            taxids[taxid] += 1
            pair_digest.update(normalized_id("@" + read_id).encode("ascii"))
            pair_digest.update(b"\n")
            length_pipe_rows += "|" in lengths
            paired_hit_separator_rows += "|:|" in hit_list
            rows += 1
    return {
        "rows": rows,
        "classified_fragments": status["C"],
        "unclassified_fragments": status["U"],
        "other_status_rows": rows - status["C"] - status["U"],
        "distinct_assigned_taxids": len([taxid for taxid in taxids if taxid != "0"]),
        "normalized_pair_id_sha256": pair_digest.hexdigest(),
        "paired_length_rows": length_pipe_rows,
        "paired_hit_separator_rows": paired_hit_separator_rows,
    }


def parse_bracken(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in read_tsv(path):
        rows.append(
            {
                "name": row["name"],
                "taxonomy_id": row["taxonomy_id"],
                "taxonomy_lvl": row["taxonomy_lvl"],
                "kraken_assigned_reads": int(row["kraken_assigned_reads"]),
                "added_reads": int(row["added_reads"]),
                "new_est_reads": int(row["new_est_reads"]),
                "fraction_total_reads": float(row["fraction_total_reads"]),
            }
        )
    return rows


def build_classification_tables(report: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    unclassified = next(row for row in report if row["rank"] == "U")
    root = next(row for row in report if row["rank"] == "R" and row["taxid"] == "1")
    total = unclassified["direct_fragments"] + root["clade_fragments"]
    ledger = [
        {"Status": "Classified", "Fragments": root["clade_fragments"], "Percent": 100 * root["clade_fragments"] / total},
        {"Status": "Unclassified", "Fragments": unclassified["direct_fragments"], "Percent": 100 * unclassified["direct_fragments"] / total},
    ]
    counts = Counter()
    for row in report:
        rank = row["rank"]
        label = RANK_LABELS.get(rank, "Other rank")
        if rank == "U":
            continue
        counts[label] += row["direct_fragments"]
    order = ["Root", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Other rank"]
    rank_rows = [
        {
            "AssignmentRank": label,
            "DirectFragments": counts[label],
            "PercentOfAllFragments": 100 * counts[label] / total,
            "PercentOfClassifiedFragments": 100 * counts[label] / root["clade_fragments"] if root["clade_fragments"] else 0.0,
        }
        for label in order
        if counts[label] > 0
    ]
    return ledger, rank_rows


def bracken_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "taxa": len(rows),
        "sum_kraken_assigned": sum(row["kraken_assigned_reads"] for row in rows),
        "sum_added": sum(row["added_reads"] for row in rows),
        "sum_estimated": sum(row["new_est_reads"] for row in rows),
        "fraction_sum": sum(row["fraction_total_reads"] for row in rows),
    }


def total_variation(first: list[dict[str, Any]], second: list[dict[str, Any]]) -> float:
    left = {row["taxonomy_id"]: row["fraction_total_reads"] for row in first}
    right = {row["taxonomy_id"]: row["fraction_total_reads"] for row in second}
    return 0.5 * sum(abs(left.get(key, 0.0) - right.get(key, 0.0)) for key in set(left) | set(right))


def build_redistribution(report: list[dict[str, Any]], primary: list[dict[str, Any]]) -> list[dict[str, Any]]:
    species = {row["taxid"]: row for row in report if row["rank"] == "S"}
    rows: list[dict[str, Any]] = []
    for row in primary:
        source = species[row["taxonomy_id"]]
        rows.append(
            {
                "Name": row["name"],
                "TaxonomyID": row["taxonomy_id"],
                "KrakenSpeciesCladeFragments": source["clade_fragments"],
                "KrakenSpeciesDirectFragments": source["direct_fragments"],
                "BrackenKrakenAssignedFragments": row["kraken_assigned_reads"],
                "BrackenAddedFragments": row["added_reads"],
                "BrackenEstimatedFragments": row["new_est_reads"],
                "BrackenFraction": row["fraction_total_reads"],
            }
        )
    return rows


def build_sensitivity(branches: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    primary = branches["species-r150-t10"]
    config = {
        "species-r150-t10": ("Species", 150, 10),
        "species-r150-t0": ("Species", 150, 0),
        "species-r100-t10": ("Species", 100, 10),
        "genus-r150-t10": ("Genus", 150, 10),
    }
    rows: list[dict[str, Any]] = []
    for label, values in branches.items():
        metrics = bracken_metrics(values)
        rank, read_length, threshold = config[label]
        rows.append(
            {
                "Configuration": label,
                "Rank": rank,
                "ReadLengthBp": read_length,
                "ThresholdFragments": threshold,
                "Taxa": metrics["taxa"],
                "KrakenAssignedFragments": metrics["sum_kraken_assigned"],
                "AddedFragments": metrics["sum_added"],
                "EstimatedFragments": metrics["sum_estimated"],
                "FractionSum": metrics["fraction_sum"],
                "TotalVariationFromPrimary": 0.0 if label == "species-r150-t10" else (
                    total_variation(primary, values) if rank == "Species" else "NA"
                ),
            }
        )
    return rows


def resource_summary(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    result: dict[str, Any] = {}
    for key, pattern in {
        "percent_cpu": r"Percent of CPU this job got:\s*(.+)",
        "elapsed": r"Elapsed \(wall clock\) time \([^)]*\):\s*(.+)",
        "maximum_rss_kb": r"Maximum resident set size \(kbytes\):\s*([0-9]+)",
    }.items():
        match = re.search(pattern, text)
        if match:
            result[key] = int(match.group(1)) if key == "maximum_rss_kb" else match.group(1).strip()
    return result


def expected_database_md5(project_root: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in (project_root / "data/small/16-standard8-files.md5").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        values[name] = digest
    return values


def audit_database_files(project_root: Path, database_dir: Path, checks: Checks) -> dict[str, Any]:
    expected = expected_database_md5(project_root)
    rows: list[dict[str, Any]] = []
    for name, md5 in expected.items():
        path = database_dir / name
        observed = hash_file(path, "md5") if path.is_file() else "missing"
        checks.add("database", f"file-md5-{name}", observed == md5, observed)
        if path.is_file():
            rows.append({"name": name, "bytes": path.stat().st_size, "md5": observed, "sha256": hash_file(path)})
    bracken_lengths = sorted(
        int(match.group(1))
        for row in rows
        if (match := re.fullmatch(r"database([0-9]+)mers\.kmer_distrib", row["name"]))
    )
    checks.add("database", "seventeen-published-files", len(rows) == 17, len(rows))
    checks.add("database", "bracken-length-grid", bracken_lengths == [50, 75, 100, 150, 200, 250, 300], bracken_lengths)
    return {
        "file_count": len(rows),
        "total_bytes": sum(row["bytes"] for row in rows),
        "bracken_read_lengths": bracken_lengths,
        "files": rows,
    }


def profile_checks(
    report: list[dict[str, Any]],
    output_audit: dict[str, Any],
    branches: dict[str, list[dict[str, Any]]],
    checks: Checks,
) -> None:
    unclassified = next(row for row in report if row["rank"] == "U")
    root = next(row for row in report if row["rank"] == "R" and row["taxid"] == "1")
    total = unclassified["direct_fragments"] + root["clade_fragments"]
    checks.add("classification", "report-total-pairs", total == EXPECTED_INPUT_PAIRS, total)
    checks.add("classification", "output-one-row-per-pair", output_audit["rows"] == EXPECTED_INPUT_PAIRS, output_audit["rows"])
    checks.add(
        "classification",
        "output-status-conservation",
        output_audit["classified_fragments"] + output_audit["unclassified_fragments"] == EXPECTED_INPUT_PAIRS
        and output_audit["other_status_rows"] == 0,
        f"C={output_audit['classified_fragments']};U={output_audit['unclassified_fragments']}",
    )
    checks.add("classification", "report-output-classified-agree", root["clade_fragments"] == output_audit["classified_fragments"], root["clade_fragments"])
    checks.add("classification", "report-output-unclassified-agree", unclassified["direct_fragments"] == output_audit["unclassified_fragments"], unclassified["direct_fragments"])
    checks.add("classification", "paired-length-field", output_audit["paired_length_rows"] == EXPECTED_INPUT_PAIRS, output_audit["paired_length_rows"])
    checks.add("classification", "paired-hit-separator", output_audit["paired_hit_separator_rows"] == EXPECTED_INPUT_PAIRS, output_audit["paired_hit_separator_rows"])
    checks.add(
        "classification",
        "paired-id-lineage",
        output_audit["normalized_pair_id_sha256"] == EXPECTED_PAIR_ID_SHA256,
        output_audit["normalized_pair_id_sha256"],
    )
    checks.add(
        "classification",
        "direct-assignment-conservation",
        sum(row["direct_fragments"] for row in report) == EXPECTED_INPUT_PAIRS,
        sum(row["direct_fragments"] for row in report),
    )
    checks.add("classification", "species-positive", sum(row["rank"] == "S" for row in report) > 0, sum(row["rank"] == "S" for row in report))

    primary = branches["species-r150-t10"]
    primary_ids = {row["taxonomy_id"] for row in primary}
    species_report = {row["taxid"]: row for row in report if row["rank"] == "S"}
    checks.add("bracken", "primary-species-positive", len(primary) > 0, len(primary))
    checks.add("bracken", "primary-taxids-in-report", primary_ids <= set(species_report), len(primary_ids - set(species_report)))
    checks.add(
        "bracken",
        "primary-kraken-clade-parity",
        all(row["kraken_assigned_reads"] == species_report[row["taxonomy_id"]]["clade_fragments"] for row in primary),
        len(primary),
    )
    for label, rows in branches.items():
        metrics = bracken_metrics(rows)
        checks.add("bracken", f"{label}-rows-positive", metrics["taxa"] > 0, metrics["taxa"])
        # Bracken prints each fraction to five decimals. The row-wise rounding
        # error therefore grows with the number of reported taxa.
        checks.add("bracken", f"{label}-fraction-sum", math.isclose(metrics["fraction_sum"], 1.0, abs_tol=2e-3), metrics["fraction_sum"])
        checks.add(
            "bracken",
            f"{label}-row-arithmetic",
            all(row["kraken_assigned_reads"] + row["added_reads"] == row["new_est_reads"] for row in rows),
            metrics["sum_estimated"],
        )
    t0_ids = {row["taxonomy_id"] for row in branches["species-r150-t0"]}
    checks.add("sensitivity", "threshold-zero-superset", primary_ids <= t0_ids, len(t0_ids - primary_ids))
    checks.add("sensitivity", "read-length-same-threshold", len(branches["species-r100-t10"]) == len(primary), f"100={len(branches['species-r100-t10'])};150={len(primary)}")
    checks.add(
        "sensitivity",
        "read-length-changes-estimates",
        total_variation(primary, branches["species-r100-t10"]) > 0,
        total_variation(primary, branches["species-r100-t10"]),
    )


def normalize_frozen_paths(frozen_dir: Path, replacements: list[tuple[Path | str, str]]) -> dict[str, Any]:
    pairs = sorted(((str(source), target) for source, target in replacements), key=lambda item: len(item[0]), reverse=True)
    changed_files = 0
    replacement_count = 0
    for path in sorted(item for item in frozen_dir.rglob("*") if item.is_file()):
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="strict")
        updated = text
        for source, target in pairs:
            count = updated.count(source)
            replacement_count += count
            updated = updated.replace(source, target)
        if updated != text:
            path.write_text(updated, encoding="utf-8")
            changed_files += 1
    return {"status": "passed", "changed_files": changed_files, "replacement_count": replacement_count, "placeholders": sorted({target for _, target in pairs})}


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
    target = frozen_dir / "file-checksums.sha256"
    rows = []
    for path in sorted(item for item in frozen_dir.rglob("*") if item.is_file() and item != target):
        rows.append(f"{hash_file(path)}  {path.relative_to(frozen_dir).as_posix()}\n")
    target.write_text("".join(rows), encoding="utf-8")


def verify_frozen_checksums(frozen_dir: Path) -> tuple[int, list[str]]:
    target = frozen_dir / "file-checksums.sha256"
    failures: list[str] = []
    listed: set[str] = set()
    for line in target.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        listed.add(relative)
        path = frozen_dir / relative
        if not path.is_file():
            failures.append(f"missing:{relative}")
        elif hash_file(path) != digest:
            failures.append(f"sha256:{relative}")
    actual = {
        path.relative_to(frozen_dir).as_posix()
        for path in frozen_dir.rglob("*")
        if path.is_file() and path != target
    }
    failures.extend(f"unexpected:{relative}" for relative in sorted(actual - listed))
    failures.extend(f"listed-only:{relative}" for relative in sorted(listed - actual))
    return len(listed), failures


def publication_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / f"{stem}.png", dpi=350, bbox_inches="tight")
    fig.savefig(figure_dir / f"{stem}.tiff", dpi=350, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def plot_classification(ledger: list[dict[str, Any]], ranks: list[dict[str, Any]], figure_dir: Path) -> None:
    publication_style()
    colors = {"Classified": "#2F6B7C", "Unclassified": "#C7CED1"}
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3), gridspec_kw={"width_ratios": [0.78, 1.35]})
    ax = axes[0]
    bottom = 0.0
    for row in ledger:
        ax.bar([0], [row["Percent"]], bottom=bottom, width=0.56, color=colors[row["Status"]], label=row["Status"])
        if row["Percent"] >= 5:
            ax.text(0, bottom + row["Percent"] / 2, f"{row['Percent']:.1f}%\n{row['Fragments']:,}", ha="center", va="center", color="white" if row["Status"] == "Classified" else "#27323A", fontsize=9)
        bottom += row["Percent"]
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.6, 0.6)
    ax.set_xticks([0], ["Paired fragments"])
    ax.set_ylabel("Percent of 99,991 paired fragments")
    ax.set_title("A  Classification ledger", loc="left", fontweight="bold")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)

    ax = axes[1]
    shown = list(reversed(ranks))
    labels = [row["AssignmentRank"] for row in shown]
    values = [row["PercentOfClassifiedFragments"] for row in shown]
    bars = ax.barh(labels, values, color="#D7835F")
    ax.set_xlabel("Direct assignments (% of classified fragments)")
    ax.set_title("B  Where Kraken placed the fragment", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#E5E7E9", linewidth=0.7)
    ax.set_axisbelow(True)
    for bar, value in zip(bars, values):
        ax.text(value + max(values) * 0.012, bar.get_y() + bar.get_height() / 2, f"{value:.1f}%", va="center", fontsize=8)
    ax.set_xlim(0, max(values) * 1.18)
    fig.suptitle("Kraken2 reports one taxonomic call per synchronized read pair", y=1.03, fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[0])


def plot_redistribution(rows: list[dict[str, Any]], figure_dir: Path) -> None:
    publication_style()
    top = sorted(rows, key=lambda row: row["BrackenFraction"], reverse=True)[:12]
    raw_total = sum(row["KrakenSpeciesCladeFragments"] for row in rows)
    top = list(reversed(top))
    raw = [100 * row["KrakenSpeciesCladeFragments"] / raw_total for row in top]
    estimated = [100 * row["BrackenFraction"] for row in top]
    labels = [re.sub(r"^(Candidatus )", "Ca. ", row["Name"]) for row in top]
    y = np.arange(len(top))
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    for idx, (left, right) in enumerate(zip(raw, estimated)):
        ax.plot([left, right], [idx, idx], color="#B7C0C5", linewidth=2, zorder=1)
    ax.scatter(raw, y, s=45, color="#70858F", label="Kraken species-clade share", zorder=2)
    ax.scatter(estimated, y, s=48, color="#D7835F", label="Bracken estimate", zorder=3)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Within species-resolved mass (%)")
    ax.set_title("Bracken redistributes compatible higher-rank fragments", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#E5E7E9", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right")
    ax.text(0.0, -0.15, "Raw and re-estimated values use separate within-rank denominators; lines show direction, not uncertainty.", transform=ax.transAxes, fontsize=8, color="#4D5A61")
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[1])


def plot_sensitivity(branches: dict[str, list[dict[str, Any]]], sensitivity: list[dict[str, Any]], figure_dir: Path) -> None:
    publication_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.5), gridspec_kw={"width_ratios": [0.85, 1.25]})
    labels = ["S:150/t10", "S:150/t0", "S:100/t10", "G:150/t10"]
    taxa = [int(row["Taxa"]) for row in sensitivity]
    colors = ["#2F6B7C", "#7AA6B2", "#D7835F", "#9B7EBD"]
    bars = axes[0].bar(labels, taxa, color=colors)
    axes[0].set_ylabel("Reported taxa")
    axes[0].set_title("A  Rank and threshold alter feature count", loc="left", fontweight="bold")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].grid(axis="y", color="#E5E7E9", linewidth=0.7)
    axes[0].set_axisbelow(True)
    for bar, value in zip(bars, taxa):
        axes[0].text(bar.get_x() + bar.get_width() / 2, value, str(value), ha="center", va="bottom", fontsize=8)

    primary = {row["taxonomy_id"]: row["fraction_total_reads"] for row in branches["species-r150-t10"]}
    comparisons = [
        ("species-r100-t10", "100 bp distribution", "#D7835F", "o"),
        ("species-r150-t0", "Threshold 0", "#7AA6B2", "s"),
    ]
    maxima = [0.0]
    for label, legend, color, marker in comparisons:
        other = {row["taxonomy_id"]: row["fraction_total_reads"] for row in branches[label]}
        ids = sorted(set(primary) | set(other))
        x = np.array([100 * primary.get(taxid, 0.0) for taxid in ids])
        y = np.array([100 * other.get(taxid, 0.0) for taxid in ids])
        axes[1].scatter(x, y, s=24, alpha=0.75, color=color, marker=marker, label=legend)
        maxima.extend(x.tolist())
        maxima.extend(y.tolist())
    limit = max(maxima) * 1.06
    axes[1].plot([0, limit], [0, limit], linestyle="--", color="#6B7377", linewidth=1)
    axes[1].set_xlim(0, limit)
    axes[1].set_ylim(0, limit)
    axes[1].set_xlabel("Primary estimate: 150 bp, threshold 10 (%)")
    axes[1].set_ylabel("Alternative estimate (%)")
    axes[1].set_title("B  Read-length model changes abundance", loc="left", fontweight="bold")
    axes[1].grid(color="#E5E7E9", linewidth=0.7)
    axes[1].set_axisbelow(True)
    axes[1].legend(loc="lower right")
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[2])


def validate_figures(figure_dir: Path, checks: Checks) -> None:
    for stem in FIGURE_STEMS:
        for suffix in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            checks.add("figure", f"{stem}-{suffix}", path.is_file() and path.stat().st_size > 1000, path.stat().st_size if path.is_file() else 0)
        png = figure_dir / f"{stem}.png"
        if png.is_file():
            with Image.open(png) as image:
                checks.add("figure", f"{stem}-dimensions", image.width >= 2000 and image.height >= 1000, f"{image.width}x{image.height}")


def frozen_inputs(frozen_dir: Path) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    report = parse_kraken_report(frozen_dir / "kraken-report.tsv")
    branches = {
        label: parse_bracken(frozen_dir / f"bracken-{label}.tsv")
        for label in ("species-r150-t10", "species-r150-t0", "species-r100-t10", "genus-r150-t10")
    }
    return report, branches


def initialize_frozen(args: argparse.Namespace) -> None:
    required = [args.frozen_dir, args.raw_dir, args.work_dir, args.database_archive, args.database_dir]
    if any(value is None for value in required):
        raise ValueError("Initialization requires frozen/raw/work/archive/database paths")
    project_root = args.project_root.resolve()
    prefix = args.environment_prefix.resolve()
    frozen_dir = args.frozen_dir.resolve()
    raw_dir = args.raw_dir.resolve()
    work_dir = args.work_dir.resolve()
    archive = args.database_archive.resolve()
    database_dir = args.database_dir.resolve()
    checks = Checks()

    static_contract_checks(project_root, checks)
    versions = environment_checks(project_root, prefix, checks)
    clean_r1 = project_root / "data/raw/article13/ERR9765746_clean_R1.fastq.gz"
    clean_r2 = project_root / "data/raw/article13/ERR9765746_clean_R2.fastq.gz"
    clean_audit = audit_fastq_pair(clean_r1, clean_r2)
    fastq_checks(clean_audit, checks)

    archive_audit = hash_file_multi(archive)
    checks.add("database", "archive-bytes", archive_audit["bytes"] == EXPECTED_ARCHIVE_BYTES, archive_audit["bytes"])
    checks.add("database", "archive-md5", archive_audit["md5"] == EXPECTED_ARCHIVE_MD5, archive_audit["md5"])
    checks.add("database", "archive-sha256", archive_audit["sha256"] == EXPECTED_ARCHIVE_SHA256, archive_audit["sha256"])
    database_audit = audit_database_files(project_root, database_dir, checks)

    report, branches = frozen_inputs(frozen_dir)
    kraken_output_audit = audit_kraken_output(work_dir / "ERR9765746.kraken.output")
    profile_checks(report, kraken_output_audit, branches, checks)
    ledger, ranks = build_classification_tables(report)
    redistribution = build_redistribution(report, branches["species-r150-t10"])
    sensitivity = build_sensitivity(branches)
    write_tsv(frozen_dir / "classification-ledger.tsv", ledger, ["Status", "Fragments", "Percent"])
    write_tsv(frozen_dir / "rank-assignment.tsv", ranks, ["AssignmentRank", "DirectFragments", "PercentOfAllFragments", "PercentOfClassifiedFragments"])
    write_tsv(
        frozen_dir / "redistribution-audit.tsv",
        redistribution,
        ["Name", "TaxonomyID", "KrakenSpeciesCladeFragments", "KrakenSpeciesDirectFragments", "BrackenKrakenAssignedFragments", "BrackenAddedFragments", "BrackenEstimatedFragments", "BrackenFraction"],
    )
    write_tsv(
        frozen_dir / "parameter-sensitivity.tsv",
        sensitivity,
        ["Configuration", "Rank", "ReadLengthBp", "ThresholdFragments", "Taxa", "KrakenAssignedFragments", "AddedFragments", "EstimatedFragments", "FractionSum", "TotalVariationFromPrimary"],
    )

    primary_metrics = bracken_metrics(branches["species-r150-t10"])
    unclassified = next(row for row in report if row["rank"] == "U")
    root = next(row for row in report if row["rank"] == "R" and row["taxid"] == "1")
    resources = {
        "kraken2_full": resource_summary(frozen_dir / "logs/kraken2-full.resources.txt"),
        "bracken_species_r150_t10": resource_summary(frozen_dir / "logs/bracken-species-r150-t10.resources.txt"),
        "extract_database": resource_summary(frozen_dir / "logs/extract-database.resources.txt"),
    }
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "evidence_date": "2026-07-21",
        "source_project": "PRJEB52977",
        "source_sample": "SAMEA14435832",
        "run_accession": "ERR9765746",
        "input_pairs": EXPECTED_INPUT_PAIRS,
        "input_reads": EXPECTED_INPUT_READS,
        "classification_unit": "paired_fragment",
        "classified_fragments": root["clade_fragments"],
        "unclassified_fragments": unclassified["direct_fragments"],
        "classified_fraction": root["clade_fragments"] / EXPECTED_INPUT_PAIRS,
        "kraken2_version": versions["Kraken2"],
        "bracken_cli_version": versions["BrackenCLI"],
        "bracken_package_version": versions["BrackenPackage"],
        "python_version": versions["Python"],
        "database_release": "Standard-8-20260626",
        "database_archive_bytes": archive_audit["bytes"],
        "database_archive_md5": archive_audit["md5"],
        "database_archive_sha256": archive_audit["sha256"],
        "database_files": database_audit,
        "confidence": 0.0,
        "minimum_hit_groups": 2,
        "primary_bracken_rank": "S",
        "primary_bracken_read_length": 150,
        "primary_bracken_threshold": 10,
        "kraken_species_rows": sum(row["rank"] == "S" for row in report),
        "bracken_species": primary_metrics["taxa"],
        "bracken_kraken_assigned_fragments": primary_metrics["sum_kraken_assigned"],
        "bracken_added_fragments": primary_metrics["sum_added"],
        "bracken_estimated_fragments": primary_metrics["sum_estimated"],
        "bracken_fraction_sum": primary_metrics["fraction_sum"],
        "threshold_zero_species": len(branches["species-r150-t0"]),
        "read_length_100_species": len(branches["species-r100-t10"]),
        "genus_estimates": len(branches["genus-r150-t10"]),
        "read_length_total_variation": total_variation(branches["species-r150-t10"], branches["species-r100-t10"]),
        "threshold_total_variation": total_variation(branches["species-r150-t10"], branches["species-r150-t0"]),
        "clean_fastq_audit": clean_audit,
        "kraken_output_audit": kraken_output_audit,
        "resource_summary": resources,
        "environment_yaml_sha256": EXPECTED_ENV_YAML_SHA256,
        "environment_lock_sha256": EXPECTED_ENV_LOCK_SHA256,
        "environment_lock_packages": EXPECTED_ENV_PACKAGES,
        "raw_fastq_committed": False,
        "database_archive_committed": False,
        "database_index_committed": False,
        "per_fragment_output_committed": False,
        "one_time_network_access": True,
        "qa_network_access": False,
        "initialization_checks_passed": checks.passed,
        "initialization_checks_failed": checks.failed,
        "checksum_failures": 0,
    }
    write_json(frozen_dir / "run-summary.json", summary)
    write_tsv(frozen_dir / "initialization-audit.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])

    normalization = normalize_frozen_paths(
        frozen_dir,
        [
            (archive, "${KRAKEN_DB_ARCHIVE}"),
            (database_dir, "${KRAKEN_DB_DIR}"),
            (work_dir, "${WORK_DIR}"),
            (raw_dir, "${RAW_DIR}"),
            (prefix, "${KRAKEN_ENV_PREFIX}"),
            (project_root, "${PROJECT_ROOT}"),
            (Path.home(), "${HOME}"),
        ],
    )
    hits = host_path_hits(frozen_dir, [str(project_root), str(prefix), str(raw_dir), str(database_dir), str(archive)])
    if hits:
        raise RuntimeError(f"Host-specific paths remain in frozen evidence: {hits[:5]}")
    summary_path = frozen_dir / "run-summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["path_normalization"] = {**normalization, "host_specific_paths_retained": False, "remaining_hits": 0}
    write_json(summary_path, summary)
    write_frozen_checksums(frozen_dir)
    _, checksum_failures = verify_frozen_checksums(frozen_dir)
    if checksum_failures:
        raise RuntimeError(f"Frozen checksum failures: {checksum_failures}")
    if checks.failed:
        raise RuntimeError(f"Article 16 initialization failed {checks.failed} checks")


def routine_validation(args: argparse.Namespace) -> None:
    project_root = args.project_root.resolve()
    prefix = args.environment_prefix.resolve()
    frozen_dir = (args.frozen_dir or project_root / "data/small/16-kraken-bracken-frozen").resolve()
    output_dir = (args.output_dir or project_root / "results/16-kraken2-bracken").resolve()
    figure_dir = (args.figure_dir or project_root / "figures").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    static_contract_checks(project_root, checks)
    versions = environment_checks(project_root, prefix, checks)
    checksum_entries, checksum_failures = verify_frozen_checksums(frozen_dir)
    checks.add("frozen", "checksums", not checksum_failures, ";".join(checksum_failures) if checksum_failures else checksum_entries)
    run_summary = json.loads((frozen_dir / "run-summary.json").read_text(encoding="utf-8"))
    checks.add("frozen", "initialization-status", run_summary.get("status") == "passed", run_summary.get("status"))
    checks.add("frozen", "initialization-checks", run_summary.get("initialization_checks_failed") == 0, run_summary.get("initialization_checks_failed"))
    checks.add("frozen", "archive-contract", run_summary.get("database_archive_bytes") == EXPECTED_ARCHIVE_BYTES and run_summary.get("database_archive_md5") == EXPECTED_ARCHIVE_MD5 and run_summary.get("database_archive_sha256") == EXPECTED_ARCHIVE_SHA256, run_summary.get("database_archive_md5"))
    checks.add("frozen", "network-free-qa", run_summary.get("qa_network_access") is False, run_summary.get("qa_network_access"))
    checks.add("frozen", "no-large-inputs-committed", all(run_summary.get(key) is False for key in ("raw_fastq_committed", "database_archive_committed", "database_index_committed", "per_fragment_output_committed")), "FASTQ/database/per-fragment excluded")

    report, branches = frozen_inputs(frozen_dir)
    output_audit = run_summary["kraken_output_audit"]
    profile_checks(report, output_audit, branches, checks)
    checks.add("environment", "routine-version-parity", versions["Kraken2"] == run_summary["kraken2_version"] and versions["BrackenCLI"] == run_summary["bracken_cli_version"] and versions["BrackenPackage"] == run_summary["bracken_package_version"], versions)

    ledger, ranks = build_classification_tables(report)
    redistribution = build_redistribution(report, branches["species-r150-t10"])
    sensitivity = build_sensitivity(branches)
    plot_classification(ledger, ranks, figure_dir)
    plot_redistribution(redistribution, figure_dir)
    plot_sensitivity(branches, sensitivity, figure_dir)
    validate_figures(figure_dir, checks)

    write_tsv(output_dir / "source-audit.tsv", [
        {"Metric": "Input paired fragments", "Value": EXPECTED_INPUT_PAIRS, "Status": "PASS"},
        {"Metric": "Input reads", "Value": EXPECTED_INPUT_READS, "Status": "PASS"},
        {"Metric": "Normalized pair ID SHA-256", "Value": EXPECTED_PAIR_ID_SHA256, "Status": "PASS"},
    ], ["Metric", "Value", "Status"])
    write_tsv(output_dir / "database-audit.tsv", [
        {"Metric": "Release", "Value": run_summary["database_release"], "Status": "PASS"},
        {"Metric": "Archive bytes", "Value": run_summary["database_archive_bytes"], "Status": "PASS"},
        {"Metric": "Archive MD5", "Value": run_summary["database_archive_md5"], "Status": "PASS"},
        {"Metric": "Internal files", "Value": run_summary["database_files"]["file_count"], "Status": "PASS"},
    ], ["Metric", "Value", "Status"])
    write_tsv(output_dir / "classification-audit.tsv", ranks, ["AssignmentRank", "DirectFragments", "PercentOfAllFragments", "PercentOfClassifiedFragments"])
    write_tsv(output_dir / "redistribution-audit.tsv", redistribution, ["Name", "TaxonomyID", "KrakenSpeciesCladeFragments", "KrakenSpeciesDirectFragments", "BrackenKrakenAssignedFragments", "BrackenAddedFragments", "BrackenEstimatedFragments", "BrackenFraction"])
    write_tsv(output_dir / "sensitivity-audit.tsv", sensitivity, ["Configuration", "Rank", "ReadLengthBp", "ThresholdFragments", "Taxa", "KrakenAssignedFragments", "AddedFragments", "EstimatedFragments", "FractionSum", "TotalVariationFromPrimary"])
    write_tsv(output_dir / "tool-audit.tsv", [
        {"Tool": "Kraken2", "CLI": versions["Kraken2"], "Package": "2.17.1", "Status": "PASS"},
        {"Tool": "Bracken", "CLI": versions["BrackenCLI"], "Package": versions["BrackenPackage"], "Status": "PASS"},
        {"Tool": "Python", "CLI": versions["Python"], "Package": "3.12.13", "Status": "PASS"},
    ], ["Tool", "CLI", "Package", "Status"])

    payload = {
        "status": "passed" if checks.failed == 0 else "failed",
        "input_pairs": EXPECTED_INPUT_PAIRS,
        "input_reads": EXPECTED_INPUT_READS,
        "classification_unit": "paired_fragment",
        "classified_fragments": run_summary["classified_fragments"],
        "unclassified_fragments": run_summary["unclassified_fragments"],
        "classified_fraction": run_summary["classified_fraction"],
        "kraken2_version": versions["Kraken2"],
        "bracken_cli_version": versions["BrackenCLI"],
        "bracken_package_version": versions["BrackenPackage"],
        "database_release": run_summary["database_release"],
        "database_archive_md5": EXPECTED_ARCHIVE_MD5,
        "bracken_species": run_summary["bracken_species"],
        "bracken_added_fragments": run_summary["bracken_added_fragments"],
        "bracken_estimated_fragments": run_summary["bracken_estimated_fragments"],
        "threshold_zero_species": run_summary["threshold_zero_species"],
        "read_length_100_species": run_summary["read_length_100_species"],
        "genus_estimates": run_summary["genus_estimates"],
        "read_length_total_variation": run_summary["read_length_total_variation"],
        "frozen_checksum_entries": checksum_entries,
        "checksum_failures": len(checksum_failures),
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
    }
    write_json(output_dir / "validation-summary.json", payload)
    write_tsv(output_dir / "validation-audit.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])
    (output_dir / "validation.log").write_text(
        f"Article 16 validation {payload['status']}: {checks.passed} passed, {checks.failed} failed\n",
        encoding="utf-8",
    )
    if checks.failed:
        failed = [row for row in checks.rows if row["Status"] == "FAIL"]
        raise RuntimeError(f"Article 16 routine validation failed: {failed[:8]}")


def main() -> None:
    args = parse_args()
    if args.initialize_frozen:
        initialize_frozen(args)
    else:
        routine_validation(args)


if __name__ == "__main__":
    main()
