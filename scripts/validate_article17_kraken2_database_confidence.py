#!/usr/bin/env python3
"""Validate Article 17 database/confidence evidence and render audit figures.

Initialization mode is a one-time operation over ignored FASTQ, extracted
Kraken databases and per-fragment output. Routine mode is network-free and
reads only the checksum-locked frozen reports and aggregate audit tables.
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
import shutil
import subprocess
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article17-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


EXPECTED_VERSIONS = {
    "Kraken2": "2.17.1",
    "BrackenCLI": "3.0.1",
    "BrackenPackage": "3.1p1",
    "Python": "3.12.13",
}
EXPECTED_ENV_PACKAGES = 187
EXPECTED_LOCK_PREFIXES = {
    "kraken2": "kraken2-2.17.1-",
    "bracken": "bracken-3.1p1-",
    "python": "python-3.12.13-",
    "matplotlib-base": "matplotlib-base-3.10.5-",
    "pillow": "pillow-12.3.0-",
    "pyyaml": "pyyaml-6.0.3-",
}
EXPECTED_STATIC_SHA256 = {
    "env/kraken.yml": "51e81b16c47386ac58982445f686466b193999269ab9aef42342e15758a0c9af",
    "env/kraken-linux-64.lock": "7ea37cd9e01dee9d6091ca42c0da7e514b637090be9e3f3362e98a0ef8b0789e",
    "data/small/13-qc-frozen/run-summary.json": "f6dcc51b6535247de7f370dc2334994dd85dca61d3f86d252294127faa3460fe",
    "data/small/14-host-removal-frozen/run-summary.json": "0a9899f953e131fce652df5e6975c4ffcaa15c5d15b772d27f2f30b2a5d04c5a",
    "data/small/16-standard8-files.md5": "8691836fe757e975828ce709d6ff0cee668102f76421d27235a762adf3844ba7",
    "data/small/17-source-manifest.tsv": "8c21bc44c8521f2eaf2a70109b779f989a0534e7a53459d46d4750444e0269f4",
    "data/small/17-standard16-files.md5": "0cf5fadc553124c0a893901a783ec6ea541b141b950144c0f956417cb2bfa531",
    "data/small/17-pluspf8-files.md5": "e53903f879dcc0faac6b46516578e014e9d954f0b740be1e938b86ed34216ee9",
    "data/small/17-database-manifest.tsv": "ac55c8c7fa9cd1b0adfc971eaecf83b4fded9827bd3c65cdf75636921de057e1",
    "data/small/17-mock1-truth.tsv": "35354d36ebc966915ad8866eb803a6f066a2777abfb548ae8855c274f837341a",
    "data/small/17-ncbi-genome-snapshot.jsonl": "67cca077df50b81c2435a416f99eec4829f1a05e2cee8ec49cf85fec94bbc364",
    "data/small/17-ncbi-sequence-snapshot.jsonl": "f6cc020a858f2def09df94ea386b99e4ad5a9ec4bd1356169c571e32b34f1f47",
    "data/small/17-truth-provenance.tsv": "ba20c4f09cba0054802f6416d98e0c4c16e6afc4a41a4efe4a8bdcfa0320c691",
    "data/small/17-data-NOTICE.txt": "68226212f22a1004150d23c7b0eb4ad693fe8f8641f450c56604016180a0670f",
}
DATABASES = {
    "standard8": {
        "release": "Standard-8-20260626",
        "archive_rel": "archives/kraken2-standard8-20260626/k2_standard_08_GB_20260626.tar.gz",
        "dir_rel": "standard-8-20260626",
        "bytes": 5_946_578_575,
        "md5": "7685f43cce057c2ca18511c925399b72",
        "files_md5": "data/small/16-standard8-files.md5",
        "content": "Standard",
        "cap_gb": 8,
    },
    "standard16": {
        "release": "Standard-16-20260626",
        "archive_rel": "archives/kraken2-standard16-20260626/k2_standard_16_GB_20260626.tar.gz",
        "dir_rel": "standard-16-20260626",
        "bytes": 11_995_707_291,
        "md5": "f130daa49fd0befa688330b288a623de",
        "files_md5": "data/small/17-standard16-files.md5",
        "content": "Standard",
        "cap_gb": 16,
    },
    "pluspf8": {
        "release": "PlusPF-8-20260626",
        "archive_rel": "archives/kraken2-pluspf8-20260626/k2_pluspf_08_GB_20260626.tar.gz",
        "dir_rel": "pluspf-8-20260626",
        "bytes": 5_933_654_083,
        "md5": "79a153b99f045bc2ae95e6d57c17a02d",
        "files_md5": "data/small/17-pluspf8-files.md5",
        "content": "Standard+RefSeq fungi+protozoa",
        "cap_gb": 8,
    },
}
CONTROL_EXPECTED = {
    "mock": {
        "pairs": 99_991,
        "pair_sha256": "457cef6e9d603790dfbc26b716b0498169b54c31bc903d067d449d8dcc86770d",
        "R1": {
            "bytes": 8_661_319,
            "bases": 14_974_589,
            "sha256": "ce438de2814a6e605cc24b00e53b697d018f325bc2a9694cece5be158ff32101",
        },
        "R2": {
            "bytes": 10_045_722,
            "bases": 14_835_184,
            "sha256": "207d080cd5dc98bd10fd4af8c492fb3f2000f73e1458a9d1d77e57347f4fe459",
        },
    },
    "human": {
        "pairs": 20_000,
        "pair_sha256": "f6006da01154f735c5d80dce94be41711c58e33fd1c0b4290c7fbe734a8b9a18",
        "R1": {
            "bytes": 1_202_366,
            "bases": 2_020_000,
            "sha256": "906a7ed253f496e29847a826e426201e6310b0980070d8664d1b58f57ddf8c8e",
        },
        "R2": {
            "bytes": 1_240_150,
            "bases": 2_020_000,
            "sha256": "5686d2a6080b8b9103854885d2f06507db05ae332faab1fe2d9cd18691f9b60a",
        },
    },
}
EXPECTED_TRUTH_ASSEMBLIES = 71
EXPECTED_TRUTH_SUM = 100.01215654459679
CONFIDENCES = (0.0, 0.05, 0.10, 0.20, 0.50)
FIGURE_STEMS = (
    "17-reference-coverage",
    "17-confidence-tradeoff",
    "17-database-stability",
    "17-hit-group-control",
)
TEXT_SUFFIXES = {".csv", ".json", ".log", ".md5", ".sh", ".tsv", ".txt"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--environment-prefix", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path)
    parser.add_argument("--raw-dir", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--database-root", type=Path)
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


def hash_file_multi(path: Path) -> dict[str, Any]:
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


def read_key_value_tsv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        key, value = line.split("\t", 1)
        values[key] = value
    return values


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
    for relative, expected in EXPECTED_STATIC_SHA256.items():
        path = project_root / relative
        observed = hash_file(path) if path.is_file() else "missing"
        checks.add("contract", f"sha256-{Path(relative).name}", observed == expected, observed)

    source = read_tsv(project_root / "data/small/17-source-manifest.tsv")
    checks.add("source", "four-mate-rows", len(source) == 4, len(source))
    checks.add(
        "source",
        "two-control-identities",
        {row["ControlID"] for row in source} == {"mock_truth", "human_method_control"}
        and {row["RunAccession"] for row in source} == {"ERR9765746", "ERR194147"},
        sorted({row["RunAccession"] for row in source}),
    )
    checks.add(
        "source",
        "bracken-read-length-contract",
        {(row["ControlID"], row["BrackenReadLength"]) for row in source}
        == {("mock_truth", "150"), ("human_method_control", "100")},
        sorted({(row["ControlID"], row["BrackenReadLength"]) for row in source}),
    )

    databases = read_tsv(project_root / "data/small/17-database-manifest.tsv")
    checks.add("database", "three-database-rows", len(databases) == 3, len(databases))
    checks.add(
        "database",
        "same-release",
        {row["release_date"] for row in databases} == {"2026-06-26"}
        and all("/latest/" not in row["archive_url"] for row in databases),
        sorted({row["release_id"] for row in databases}),
    )
    manifest_ids = {
        "kraken2-standard8-20260626": "standard8",
        "kraken2-standard16-20260626": "standard16",
        "kraken2-pluspf8-20260626": "pluspf8",
    }
    expected_installed_bytes = {
        "standard8": 8_638_898_739,
        "standard16": 16_636_247_365,
        "pluspf8": 8_643_693_403,
    }
    checks.add("database", "database-identities", set(manifest_ids) == {row["database_id"] for row in databases}, sorted(row["database_id"] for row in databases))
    for row in databases:
        expected = DATABASES[manifest_ids[row["database_id"]]]
        short_id = manifest_ids[row["database_id"]]
        checks.add(
            "database",
            f"archive-{row['release_id']}",
            int(row["expected_compressed_bytes"]) == expected["bytes"]
            and row["expected_checksum"] == expected["md5"],
            f"{row['expected_compressed_bytes']};{row['expected_checksum']}",
        )
        checks.add(
            "database",
            f"installed-{row['release_id']}",
            int(row["expected_installed_bytes"]) == expected_installed_bytes[short_id]
            and row["validation_status"] == "VERIFIED_UPSTREAM_AND_LOCAL",
            f"{row['expected_installed_bytes']};{row['validation_status']}",
        )

    truth = read_tsv(project_root / "data/small/17-mock1-truth.tsv")
    truth_sum = math.fsum(float(row["ExpectedGenomePercent"]) for row in truth)
    checks.add("truth", "truth-assembly-count", len(truth) == EXPECTED_TRUTH_ASSEMBLIES, len(truth))
    checks.add("truth", "truth-accessions-unique", len({row["AssemblyAccession"] for row in truth}) == 71, len({row["AssemblyAccession"] for row in truth}))
    checks.add("truth", "truth-source-sum-preserved", abs(truth_sum - EXPECTED_TRUTH_SUM) < 1e-12, repr(truth_sum))
    checks.add("truth", "genome-snapshot-rows", sum(1 for _ in (project_root / "data/small/17-ncbi-genome-snapshot.jsonl").open()) == 71, 71)
    checks.add("truth", "sequence-snapshot-nonempty", (project_root / "data/small/17-ncbi-sequence-snapshot.jsonl").stat().st_size > 50_000, (project_root / "data/small/17-ncbi-sequence-snapshot.jsonl").stat().st_size)


def environment_checks(project_root: Path, prefix: Path, checks: Checks) -> dict[str, str]:
    lock_lines = (project_root / "env/kraken-linux-64.lock").read_text(encoding="utf-8").splitlines()
    package_lines = [line for line in lock_lines if line.startswith("http")]
    checks.add("environment", "lock-package-count", len(package_lines) == EXPECTED_ENV_PACKAGES, len(package_lines))
    for name, marker in EXPECTED_LOCK_PREFIXES.items():
        matches = [line for line in package_lines if marker in line]
        checks.add("environment", f"lock-{name}", len(matches) == 1, len(matches))

    env = os.environ.copy()
    env["PATH"] = f"{prefix / 'bin'}:{env.get('PATH', '')}"
    kraken = subprocess.run([str(prefix / "bin/kraken2"), "--version"], check=True, capture_output=True, text=True, env=env).stdout
    bracken = subprocess.run([str(prefix / "bin/bracken"), "-v"], check=True, capture_output=True, text=True, env=env).stdout
    python = subprocess.run([str(prefix / "bin/python"), "--version"], check=True, capture_output=True, text=True, env=env).stdout
    versions = {
        "Kraken2": re.search(r"Kraken version ([0-9.]+)", kraken).group(1),
        "BrackenCLI": re.search(r"Bracken v([^\s]+)", bracken).group(1),
        "Python": re.search(r"Python ([0-9.]+)", python).group(1),
    }
    meta = list((prefix / "conda-meta").glob("bracken-*.json"))
    versions["BrackenPackage"] = json.loads(meta[0].read_text(encoding="utf-8"))["version"] if len(meta) == 1 else "missing"
    for key, expected in EXPECTED_VERSIONS.items():
        checks.add("environment", f"version-{key.lower()}", versions[key] == expected, versions[key])
    checks.add(
        "environment",
        "bracken-cli-package-discrepancy-recorded",
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
    iterator1, iterator2 = fastq_records(r1), fastq_records(r2)
    while True:
        record1, record2 = next(iterator1, None), next(iterator2, None)
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


def audit_fastq_controls(project_root: Path, checks: Checks) -> dict[str, Any]:
    paths = {
        "mock": (
            project_root / "data/raw/article13/ERR9765746_clean_R1.fastq.gz",
            project_root / "data/raw/article13/ERR9765746_clean_R2.fastq.gz",
        ),
        "human": (
            project_root / "data/raw/article14/ERR194147_prefix20k_R1.fastq.gz",
            project_root / "data/raw/article14/ERR194147_prefix20k_R2.fastq.gz",
        ),
    }
    audits: dict[str, Any] = {}
    for control, (r1, r2) in paths.items():
        audit = audit_fastq_pair(r1, r2)
        expected = CONTROL_EXPECTED[control]
        checks.add("source", f"{control}-pairs", audit["pairs"] == expected["pairs"], audit["pairs"])
        checks.add("source", f"{control}-pair-id-sha256", audit["normalized_pair_id_sha256"] == expected["pair_sha256"], audit["normalized_pair_id_sha256"])
        for mate in ("R1", "R2"):
            for key in ("bytes", "bases", "sha256"):
                checks.add("source", f"{control}-{mate.lower()}-{key}", audit[mate][key] == expected[mate][key], audit[mate][key])
        audits[control] = audit
    return audits


def parse_nodes(path: Path) -> tuple[dict[int, int], dict[int, str]]:
    parents: dict[int, int] = {}
    ranks: dict[int, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t|\t")
            taxid = int(fields[0])
            parents[taxid] = int(fields[1])
            ranks[taxid] = fields[2]
    return parents, ranks


def species_ancestor(taxid: int, parents: dict[int, int], ranks: dict[int, str], cache: dict[int, int | None]) -> int | None:
    if taxid in cache:
        return cache[taxid]
    trail: list[int] = []
    seen: set[int] = set()
    current = taxid
    answer: int | None = None
    while current in parents and current not in seen:
        if current in cache:
            answer = cache[current]
            break
        trail.append(current)
        seen.add(current)
        if ranks.get(current) == "species":
            answer = current
            break
        parent = parents[current]
        if parent == current:
            break
        current = parent
    for item in trail:
        cache[item] = answer
    return answer


def ancestor_taxids(taxid: int, parents: dict[int, int]) -> set[int]:
    ancestors: set[int] = set()
    current = parents.get(taxid)
    while current is not None and current not in ancestors:
        ancestors.add(current)
        parent = parents.get(current)
        if parent is None or parent == current:
            break
        current = parent
    return ancestors


def selected_scientific_names(path: Path, wanted: set[int]) -> dict[int, str]:
    names: dict[int, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t|\t")
            taxid = int(fields[0])
            if taxid in wanted and len(fields) > 3 and fields[3].startswith("scientific name"):
                names[taxid] = fields[1]
    return names


def parse_expected_md5(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        checksum, filename = line.split(maxsplit=1)
        rows[filename.strip()] = checksum
    return rows


def build_reference_crosswalk(
    project_root: Path,
    database_root: Path,
    checks: Checks,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    truth = read_tsv(project_root / "data/small/17-mock1-truth.tsv")
    first_db = database_root / DATABASES["standard8"]["dir_rel"]
    parents, ranks = parse_nodes(first_db / "nodes.dmp")
    cache: dict[int, int | None] = {}
    truth_species: dict[str, int] = {}
    for row in truth:
        species = species_ancestor(int(row["NCBITaxID"]), parents, ranks, cache)
        if species is None:
            raise ValueError(f"No species ancestor for {row['AssemblyAccession']}")
        truth_species[row["AssemblyAccession"]] = species
    species_names = selected_scientific_names(first_db / "names.dmp", set(truth_species.values()) | {9606})

    crosswalk: list[dict[str, Any]] = []
    database_audit: list[dict[str, Any]] = []
    archive_summaries: dict[str, Any] = {}
    for database_id, contract in DATABASES.items():
        archive = database_root / contract["archive_rel"]
        database_dir = database_root / contract["dir_rel"]
        archive_summary = hash_file_multi(archive)
        archive_summaries[database_id] = archive_summary
        checks.add("database", f"{database_id}-archive-bytes", archive_summary["bytes"] == contract["bytes"], archive_summary["bytes"])
        checks.add("database", f"{database_id}-archive-md5", archive_summary["md5"] == contract["md5"], archive_summary["md5"])

        expected_files = parse_expected_md5(project_root / contract["files_md5"])
        file_failures: list[str] = []
        installed_bytes = 0
        for filename, expected_md5 in expected_files.items():
            path = database_dir / filename
            if not path.is_file():
                file_failures.append(f"missing:{filename}")
                continue
            installed_bytes += path.stat().st_size
            observed = hash_file(path, "md5")
            if observed != expected_md5:
                file_failures.append(f"md5:{filename}")
        checks.add("database", f"{database_id}-internal-md5", not file_failures, ";".join(file_failures) if file_failures else len(expected_files))

        seqids: set[str] = set()
        reference_taxids: set[int] = set()
        with (database_dir / "seqid2taxid.map").open(encoding="utf-8") as handle:
            for line in handle:
                sequence_id, taxid = line.rstrip("\n").split("\t")
                seqids.add(sequence_id.split("|")[-1])
                reference_taxids.add(int(taxid))
        represented_species = {
            species
            for taxid in reference_taxids
            if (species := species_ancestor(taxid, parents, ranks, cache)) is not None
        }

        status_counts: Counter[str] = Counter()
        status_abundance: defaultdict[str, float] = defaultdict(float)
        exact_matching_sequences = 0
        for row in truth:
            accession = row["AssemblyAccession"]
            species = truth_species[accession]
            candidates = {
                item
                for item in (row["RefSeqSequenceAccessions"] + ";" + row["GenBankSequenceAccessions"]).split(";")
                if item
            }
            matched = sorted(candidates & seqids)
            if matched:
                status = "Exact expected assembly"
            elif species in represented_species:
                status = "Alternate same-species reference"
            else:
                status = "No same-species reference"
            expected_abundance = float(row["ExpectedGenomePercent"])
            status_counts[status] += 1
            status_abundance[status] += expected_abundance
            exact_matching_sequences += len(matched)
            crosswalk.append(
                {
                    "DatabaseID": database_id,
                    "ReleaseID": contract["release"],
                    "TruthIndex": int(row["TruthIndex"]),
                    "ExpectedOrganismGTDBRS207": row["ExpectedOrganismGTDBRS207"],
                    "StrainName": row["StrainName"],
                    "AssemblyAccession": accession,
                    "NCBITaxID": int(row["NCBITaxID"]),
                    "NCBIOrganismName": row["NCBIOrganismName"],
                    "SpeciesTaxID": species,
                    "DatabaseSpeciesName": species_names.get(species, ""),
                    "ExpectedGenomePercent": expected_abundance,
                    "CandidateSequenceAccessions": len(candidates),
                    "ExactMatchingSequences": len(matched),
                    "ReferenceStatus": status,
                }
            )

        present_species = {
            int(row["SpeciesTaxID"])
            for row in crosswalk
            if row["DatabaseID"] == database_id and row["ReferenceStatus"] != "No same-species reference"
        }
        expected_species = {
            int(row["SpeciesTaxID"])
            for row in crosswalk
            if row["DatabaseID"] == database_id
        }
        database_audit.append(
            {
                "DatabaseID": database_id,
                "ReleaseID": contract["release"],
                "Content": contract["content"],
                "CapGB": contract["cap_gb"],
                "ArchiveBytes": archive_summary["bytes"],
                "ArchiveMD5": archive_summary["md5"],
                "ArchiveSHA256": archive_summary["sha256"],
                "InstalledFiles": len(expected_files),
                "InstalledBytes": installed_bytes,
                "SequenceMapRows": len(seqids),
                "ReferenceTaxIDs": len(reference_taxids),
                "TruthAssemblies": len(truth),
                "TruthSpecies": len(expected_species),
                "ReferencePresentTruthSpecies": len(present_species),
                "ExactAssemblies": status_counts["Exact expected assembly"],
                "AlternateReferenceAssemblies": status_counts["Alternate same-species reference"],
                "NoSameSpeciesReferenceAssemblies": status_counts["No same-species reference"],
                "ExactExpectedPercent": status_abundance["Exact expected assembly"],
                "AlternateExpectedPercent": status_abundance["Alternate same-species reference"],
                "NoReferenceExpectedPercent": status_abundance["No same-species reference"],
                "ExactMatchingSequences": exact_matching_sequences,
                "InternalMD5Failures": len(file_failures),
            }
        )
        checks.add("reference", f"{database_id}-crosswalk-rows", sum(row["DatabaseID"] == database_id for row in crosswalk) == 71, status_counts)
        checks.add("reference", f"{database_id}-three-status-total", sum(status_counts.values()) == 71, dict(status_counts))
        checks.add("reference", f"{database_id}-abundance-total", abs(sum(status_abundance.values()) - EXPECTED_TRUTH_SUM) < 1e-10, sum(status_abundance.values()))
        checks.add(
            "reference",
            f"{database_id}-reference-mapping-present",
            status_counts["Exact expected assembly"]
            + status_counts["Alternate same-species reference"]
            > 0,
            status_counts["Exact expected assembly"]
            + status_counts["Alternate same-species reference"],
        )
    return crosswalk, database_audit, archive_summaries


def parse_kraken_report(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"Unexpected Kraken report row in {path}: {line[:120]}")
            rows.append(
                {
                    "percentage": float(fields[0]),
                    "clade": int(fields[1]),
                    "direct": int(fields[2]),
                    "rank": fields[3],
                    "taxid": int(fields[4]),
                    "name": fields[5].strip(),
                }
            )
    return rows


def parse_bracken(path: Path) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for row in read_tsv(path):
        output.append(
            {
                "name": row["name"],
                "taxid": int(row["taxonomy_id"]),
                "rank": row["taxonomy_lvl"],
                "kraken_assigned": int(row["kraken_assigned_reads"]),
                "added": int(row["added_reads"]),
                "estimated": int(row["new_est_reads"]),
                "fraction": float(row["fraction_total_reads"]),
            }
        )
    return output


def audit_kraken_output(path: Path, with_pair_hash: bool = False) -> dict[str, Any]:
    statuses: Counter[str] = Counter()
    pair_digest = hashlib.sha256()
    rows = 0
    paired_length_rows = 0
    paired_hit_rows = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ValueError(f"Unexpected Kraken per-fragment row in {path}")
            statuses[fields[0]] += 1
            rows += 1
            paired_length_rows += "|" in fields[3]
            paired_hit_rows += "|:|" in fields[4]
            if with_pair_hash:
                pair_digest.update(fields[1].encode("ascii"))
                pair_digest.update(b"\n")
    return {
        "rows": rows,
        "classified": statuses["C"],
        "unclassified": statuses["U"],
        "paired_length_rows": paired_length_rows,
        "paired_hit_rows": paired_hit_rows,
        "pair_id_sha256": pair_digest.hexdigest() if with_pair_hash else "not_recomputed",
    }


def report_totals(report: list[dict[str, Any]]) -> dict[str, int]:
    unclassified_rows = [row for row in report if row["rank"] == "U" and row["taxid"] == 0]
    root_rows = [row for row in report if row["rank"] == "R" and row["taxid"] == 1]
    if len(unclassified_rows) > 1 or len(root_rows) > 1 or not (unclassified_rows or root_rows):
        raise ValueError("Kraken report must contain at most one unclassified row and one root row")
    return {
        "unclassified": unclassified_rows[0]["clade"] if unclassified_rows else 0,
        "classified": root_rows[0]["clade"] if root_rows else 0,
    }


def truth_sets(crosswalk: list[dict[str, Any]], database_id: str) -> tuple[set[int], set[int], dict[int, float]]:
    rows = [row for row in crosswalk if row["DatabaseID"] == database_id]
    all_species = {int(row["SpeciesTaxID"]) for row in rows}
    present_species = {
        int(row["SpeciesTaxID"])
        for row in rows
        if row["ReferenceStatus"] != "No same-species reference"
    }
    abundance: defaultdict[int, float] = defaultdict(float)
    for row in rows:
        abundance[int(row["SpeciesTaxID"])] += float(row["ExpectedGenomePercent"])
    return all_species, present_species, dict(abundance)


def total_variation_dict(first: dict[int, float], second: dict[int, float]) -> float:
    first_total, second_total = sum(first.values()), sum(second.values())
    if first_total <= 0 or second_total <= 0:
        return math.nan
    keys = set(first) | set(second)
    return 0.5 * sum(abs(first.get(key, 0.0) / first_total - second.get(key, 0.0) / second_total) for key in keys)


def branch_metrics(
    design: list[dict[str, str]],
    frozen_dir: Path,
    work_dir: Path,
    crosswalk: list[dict[str, Any]],
    human_ancestor_taxids: set[int],
    checks: Checks,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[int, float]]]:
    classifications: list[dict[str, Any]] = []
    truth_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    distributions: dict[str, dict[int, float]] = {}
    representative_hashed: set[str] = set()
    for branch in design:
        branch_id = branch["BranchID"]
        database_id = branch["DatabaseID"]
        control = branch["ControlID"]
        input_pairs = CONTROL_EXPECTED[control]["pairs"]
        report_path = frozen_dir / "reports" / f"{branch_id}.kreport.tsv"
        bracken_path = frozen_dir / "bracken" / f"{branch_id}.tsv"
        output_path = work_dir / "per-fragment" / f"{branch_id}.output"
        report = parse_kraken_report(report_path)
        bracken = parse_bracken(bracken_path)
        bracken_sentinel = read_key_value_tsv(
            work_dir / "sentinels" / f"{branch_id}.bracken.done.tsv"
        )
        bracken_outcome = bracken_sentinel.get("outcome", "estimated")
        totals = report_totals(report)
        with_hash = control not in representative_hashed
        output = audit_kraken_output(output_path, with_pair_hash=with_hash)
        representative_hashed.add(control)
        checks.add("branch", f"{branch_id}-input-conservation", totals["classified"] + totals["unclassified"] == input_pairs, f"{totals['classified']}+{totals['unclassified']}")
        checks.add("branch", f"{branch_id}-output-rows", output["rows"] == input_pairs, output["rows"])
        checks.add("branch", f"{branch_id}-output-report-parity", output["classified"] == totals["classified"] and output["unclassified"] == totals["unclassified"], f"{output['classified']};{output['unclassified']}")
        checks.add("branch", f"{branch_id}-paired-structure", output["paired_length_rows"] == input_pairs and output["paired_hit_rows"] == input_pairs, f"{output['paired_length_rows']};{output['paired_hit_rows']}")
        if with_hash:
            checks.add("branch", f"{control}-output-pair-hash", output["pair_id_sha256"] == CONTROL_EXPECTED[control]["pair_sha256"], output["pair_id_sha256"])
        checks.add("branch", f"{branch_id}-bracken-rank", all(row["rank"] == "S" for row in bracken), len(bracken))
        expected_outcome = (
            "estimated" if bracken else "no_eligible_species_above_threshold"
        )
        checks.add(
            "branch",
            f"{branch_id}-bracken-outcome",
            bracken_outcome == expected_outcome,
            bracken_outcome,
        )

        species_rows = [row for row in report if row["rank"] == "S"]
        species_counts = {row["taxid"]: row["clade"] for row in species_rows}
        species_total = sum(species_counts.values())
        bracken_distribution = {row["taxid"]: row["estimated"] for row in bracken if row["estimated"] > 0}
        distributions[branch_id] = bracken_distribution
        bracken_total = sum(bracken_distribution.values())
        human_clade = next((row["clade"] for row in species_rows if row["taxid"] == 9606), 0)
        human_ancestor_fragments = sum(
            row["direct"] for row in report if row["taxid"] in human_ancestor_taxids
        )
        classifications.append(
            {
                **branch,
                "InputFragments": input_pairs,
                "ClassifiedFragments": totals["classified"],
                "UnclassifiedFragments": totals["unclassified"],
                "ClassifiedFraction": totals["classified"] / input_pairs,
                "SpeciesResolvedFragments": species_total,
                "SpeciesResolvedFraction": species_total / input_pairs,
                "DetectedSpecies": sum(count > 0 for count in species_counts.values()),
                "BrackenSpecies": len(bracken_distribution),
                "BrackenEstimatedFragments": bracken_total,
                "BrackenOutcome": bracken_outcome,
                "HumanCladeFragments": human_clade,
                "HumanCompatibleAncestorFragments": human_ancestor_fragments,
            }
        )

        all_truth, present_truth, expected_abundance = truth_sets(crosswalk, database_id)
        if control == "mock":
            detected_present = {taxid for taxid in present_truth if species_counts.get(taxid, 0) > 0}
            detected_truth = {taxid for taxid in all_truth if species_counts.get(taxid, 0) > 0}
            nontruth = {taxid: count for taxid, count in species_counts.items() if taxid not in all_truth}
            truth_fragments = sum(species_counts.get(taxid, 0) for taxid in all_truth)
            nontruth_fragments = sum(nontruth.values())
            bracken_truth = sum(bracken_distribution.get(taxid, 0) for taxid in all_truth)
            bracken_nontruth = sum(value for taxid, value in bracken_distribution.items() if taxid not in all_truth)
            expected_present_distribution = {taxid: expected_abundance[taxid] for taxid in present_truth}
            observed_present_distribution = {taxid: bracken_distribution.get(taxid, 0) for taxid in present_truth}
            truth_rows.append(
                {
                    **branch,
                    "TruthAssemblies": EXPECTED_TRUTH_ASSEMBLIES,
                    "TruthSpecies": len(all_truth),
                    "ReferencePresentTruthSpecies": len(present_truth),
                    "DetectedReferencePresentSpecies": len(detected_present),
                    "DetectedTruthSpecies": len(detected_truth),
                    "ReferenceAwareSpeciesRecall": len(detected_present) / len(present_truth),
                    "DetectedSpeciesPrecision": len(detected_truth) / len(species_counts) if species_counts else math.nan,
                    "TruthSpeciesFragments": truth_fragments,
                    "NonTruthSpeciesFragments": nontruth_fragments,
                    "TruthSpeciesFragmentPrecision": truth_fragments / species_total if species_total else math.nan,
                    "NonTruthSpeciesPer10k": nontruth_fragments / input_pairs * 10_000,
                    "BrackenTruthFragments": bracken_truth,
                    "BrackenNonTruthFragments": bracken_nontruth,
                    "BrackenTruthFraction": bracken_truth / bracken_total if bracken_total else math.nan,
                    "BrackenReferenceAwareTV": total_variation_dict(expected_present_distribution, observed_present_distribution),
                    "ReferenceCoveredExpectedPercent": sum(expected_abundance[taxid] for taxid in present_truth),
                    "NoReferenceExpectedPercent": sum(expected_abundance[taxid] for taxid in all_truth - present_truth),
                }
            )
        else:
            unsupported_classified = (
                totals["classified"] - human_clade - human_ancestor_fragments
            )
            checks.add(
                "control",
                f"{branch_id}-human-compatible-partition",
                unsupported_classified >= 0
                and human_clade
                + human_ancestor_fragments
                + unsupported_classified
                == totals["classified"],
                f"{human_clade}+{human_ancestor_fragments}+{unsupported_classified}",
            )
            nonhuman_species = species_total - human_clade
            bracken_human = bracken_distribution.get(9606, 0)
            bracken_nonhuman = bracken_total - bracken_human
            control_rows.append(
                {
                    **branch,
                    "InputFragments": input_pairs,
                    "ClassifiedFragments": totals["classified"],
                    "HumanCladeFragments": human_clade,
                    "HumanCompatibleAncestorFragments": human_ancestor_fragments,
                    "UnsupportedClassifiedFragments": unsupported_classified,
                    "UnsupportedClassifiedPer10k": unsupported_classified / input_pairs * 10_000,
                    "NonHumanSpeciesFragments": nonhuman_species,
                    "NonHumanSpeciesPer10k": nonhuman_species / input_pairs * 10_000,
                    "BrackenHumanFragments": bracken_human,
                    "BrackenNonHumanFragments": bracken_nonhuman,
                    "BrackenNonHumanPer10k": bracken_nonhuman / input_pairs * 10_000,
                    "Boundary": "Human-rich method control; not an extraction blank",
                }
            )
    return classifications, truth_rows, control_rows, distributions


def abundance_stability(
    design: list[dict[str, str]],
    distributions: dict[str, dict[int, float]],
) -> list[dict[str, Any]]:
    primary = [row for row in design if row["PrimaryMatrix"] == "Yes"]
    lookup = {(row["ControlID"], float(row["Confidence"]), row["DatabaseID"]): row["BranchID"] for row in primary}
    output: list[dict[str, Any]] = []
    for control in ("mock", "human"):
        for confidence in CONFIDENCES:
            for first, second in combinations(DATABASES, 2):
                first_id = lookup[(control, confidence, first)]
                second_id = lookup[(control, confidence, second)]
                first_dist, second_dist = distributions[first_id], distributions[second_id]
                first_taxa, second_taxa = set(first_dist), set(second_dist)
                union = first_taxa | second_taxa
                output.append(
                    {
                        "ControlID": control,
                        "Confidence": confidence,
                        "MinimumHitGroups": 2,
                        "DatabaseA": first,
                        "DatabaseB": second,
                        "BranchA": first_id,
                        "BranchB": second_id,
                        "TaxaA": len(first_taxa),
                        "TaxaB": len(second_taxa),
                        "TaxonJaccard": len(first_taxa & second_taxa) / len(union) if union else 1.0,
                        "BrackenTotalVariation": total_variation_dict(first_dist, second_dist),
                    }
                )
    return output


def build_hit_group_table(
    design: list[dict[str, str]],
    classifications: list[dict[str, Any]],
    truth_rows: list[dict[str, Any]],
    control_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    classification_by_branch = {row["BranchID"]: row for row in classifications}
    truth_by_branch = {row["BranchID"]: row for row in truth_rows}
    control_by_branch = {row["BranchID"]: row for row in control_rows}
    output: list[dict[str, Any]] = []
    for row in design:
        if row["HitGroupMatrix"] != "Yes":
            continue
        branch_id = row["BranchID"]
        base = classification_by_branch[branch_id]
        result = {
            "BranchID": branch_id,
            "ControlID": row["ControlID"],
            "Confidence": float(row["Confidence"]),
            "MinimumHitGroups": int(row["MinimumHitGroups"]),
            "ClassifiedFraction": base["ClassifiedFraction"],
            "SpeciesResolvedFraction": base["SpeciesResolvedFraction"],
        }
        if row["ControlID"] == "mock":
            metrics = truth_by_branch[branch_id]
            result.update(
                {
                    "ReferenceAwareSpeciesRecall": metrics["ReferenceAwareSpeciesRecall"],
                    "NonTruthSpeciesPer10k": metrics["NonTruthSpeciesPer10k"],
                    "UnsupportedClassifiedPer10k": "",
                    "NonHumanSpeciesPer10k": "",
                }
            )
        else:
            metrics = control_by_branch[branch_id]
            result.update(
                {
                    "ReferenceAwareSpeciesRecall": "",
                    "NonTruthSpeciesPer10k": "",
                    "UnsupportedClassifiedPer10k": metrics["UnsupportedClassifiedPer10k"],
                    "NonHumanSpeciesPer10k": metrics["NonHumanSpeciesPer10k"],
                }
            )
        output.append(result)
    return sorted(output, key=lambda item: (item["ControlID"], item["MinimumHitGroups"]))


def elapsed_seconds(value: str) -> float:
    pieces = value.strip().split(":")
    if len(pieces) == 2:
        minutes, seconds = pieces
        return int(minutes) * 60 + float(seconds)
    if len(pieces) == 3:
        hours, minutes, seconds = pieces
        return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    raise ValueError(f"Unexpected elapsed time: {value}")


def parse_resource_file(path: Path) -> dict[str, Any]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("Elapsed (wall clock) time"):
            values["elapsed"] = stripped.rsplit(": ", 1)[1]
        elif stripped.startswith("Maximum resident set size (kbytes):"):
            values["rss_kb"] = stripped.rsplit(":", 1)[1].strip()
        elif stripped.startswith("User time (seconds):"):
            values["user_seconds"] = stripped.rsplit(":", 1)[1].strip()
        elif stripped.startswith("System time (seconds):"):
            values["system_seconds"] = stripped.rsplit(":", 1)[1].strip()
        elif stripped.startswith("Exit status:"):
            values["exit_status"] = stripped.rsplit(":", 1)[1].strip()
    required = {"elapsed", "rss_kb", "user_seconds", "system_seconds", "exit_status"}
    if not required <= set(values):
        raise ValueError(f"Incomplete resource file {path}: {values}")
    return {
        "WallSeconds": elapsed_seconds(values["elapsed"]),
        "MaxRSSMiB": int(values["rss_kb"]) / 1024,
        "UserSeconds": float(values["user_seconds"]),
        "SystemSeconds": float(values["system_seconds"]),
        "ExitStatus": int(values["exit_status"]),
    }


def resource_table(design: list[dict[str, str]], work_dir: Path, checks: Checks) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for row in design:
        for tool in ("kraken", "bracken"):
            resource_path = work_dir / "resources" / f"{row['BranchID']}.{tool}.txt"
            metrics = parse_resource_file(resource_path)
            outcome = "classified"
            accepted = metrics["ExitStatus"] == 0
            if tool == "bracken":
                sentinel = read_key_value_tsv(
                    work_dir / "sentinels" / f"{row['BranchID']}.bracken.done.tsv"
                )
                outcome = sentinel.get("outcome", "estimated")
                accepted = (
                    (outcome == "estimated" and metrics["ExitStatus"] == 0)
                    or (
                        outcome == "no_eligible_species_above_threshold"
                        and metrics["ExitStatus"] == 1
                    )
                )
            checks.add("resource", f"{row['BranchID']}-{tool}-exit", accepted, f"{metrics['ExitStatus']};{outcome}")
            output.append({**row, "Tool": "Kraken2" if tool == "kraken" else "Bracken", "Outcome": outcome, **metrics})
    return output


def write_frozen_checksums(frozen_dir: Path) -> int:
    entries: list[tuple[str, str]] = []
    for path in sorted(frozen_dir.rglob("*")):
        if not path.is_file() or path.name == "file-checksums.sha256" or path.name.endswith(".tmp"):
            continue
        entries.append((hash_file(path), path.relative_to(frozen_dir).as_posix()))
    (frozen_dir / "file-checksums.sha256").write_text(
        "".join(f"{digest}  {relative}\n" for digest, relative in entries),
        encoding="utf-8",
    )
    return len(entries)


def verify_frozen_checksums(frozen_dir: Path) -> tuple[int, list[str]]:
    manifest = frozen_dir / "file-checksums.sha256"
    failures: list[str] = []
    entries = 0
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        path = frozen_dir / relative
        entries += 1
        if not path.is_file():
            failures.append(f"missing:{relative}")
        elif hash_file(path) != expected:
            failures.append(f"sha256:{relative}")
    return entries, failures


def set_publication_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(
        figure_dir / f"{stem}.tiff",
        dpi=350,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


DB_ORDER = ("standard8", "standard16", "pluspf8")
DB_LABELS = {"standard8": "Standard-8", "standard16": "Standard-16", "pluspf8": "PlusPF-8"}
DB_COLORS = {"standard8": "#3B6FB6", "standard16": "#D35F5F", "pluspf8": "#2B8C7E"}


def plot_reference_coverage(database_audit: list[dict[str, Any]], figure_dir: Path) -> None:
    set_publication_style()
    categories = (
        ("ExactAssemblies", "ExactExpectedPercent", "Exact assembly", "#3B6FB6"),
        ("AlternateReferenceAssemblies", "AlternateExpectedPercent", "Alternate same species", "#E69F00"),
        ("NoSameSpeciesReferenceAssemblies", "NoReferenceExpectedPercent", "No same-species reference", "#8C8C8C"),
    )
    by_db = {row["DatabaseID"]: row for row in database_audit}
    x = np.arange(len(DB_ORDER))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.7), constrained_layout=True)
    bottoms = np.zeros(len(DB_ORDER))
    for count_key, _, label, color in categories:
        values = np.array([float(by_db[db][count_key]) for db in DB_ORDER])
        axes[0].bar(x, values, bottom=bottoms, color=color, label=label, width=0.68)
        for i, (bottom, value) in enumerate(zip(bottoms, values)):
            if value >= 3:
                axes[0].text(i, bottom + value / 2, f"{int(value)}", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
        bottoms += values
    axes[0].set_xticks(x, [DB_LABELS[db] for db in DB_ORDER])
    axes[0].set_ylabel("Expected assemblies")
    axes[0].set_title("A  Reference representation by assembly")
    axes[0].set_ylim(0, 75)

    bottoms = np.zeros(len(DB_ORDER))
    for _, percent_key, label, color in categories:
        values = np.array([float(by_db[db][percent_key]) for db in DB_ORDER])
        axes[1].bar(x, values, bottom=bottoms, color=color, label=label, width=0.68)
        bottoms += values
    axes[1].set_xticks(x, [DB_LABELS[db] for db in DB_ORDER])
    axes[1].set_ylabel("Expected genome abundance (%)")
    axes[1].set_title("B  Expected mass behind each reference class")
    axes[1].axhline(EXPECTED_TRUTH_SUM, color="#333333", lw=0.8, ls="--")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.05))
    save_figure(fig, figure_dir, "17-reference-coverage")


def plot_confidence_tradeoff(
    truth_rows: list[dict[str, Any]],
    control_rows: list[dict[str, Any]],
    figure_dir: Path,
) -> None:
    set_publication_style()
    truth = [row for row in truth_rows if row["PrimaryMatrix"] == "Yes"]
    controls = [row for row in control_rows if row["PrimaryMatrix"] == "Yes"]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5), constrained_layout=True)
    x_positions = np.arange(len(CONFIDENCES))
    x_labels = ("0", "0.05", "0.10", "0.20", "0.50")
    for db in DB_ORDER:
        subset = sorted((row for row in truth if row["DatabaseID"] == db), key=lambda row: float(row["Confidence"]))
        axes[0].plot(x_positions, [100 * float(row["ReferenceAwareSpeciesRecall"]) for row in subset], marker="o", color=DB_COLORS[db], label=DB_LABELS[db])
        axes[1].plot(x_positions, [float(row["NonTruthSpeciesPer10k"]) for row in subset], marker="o", color=DB_COLORS[db])
        control_subset = sorted((row for row in controls if row["DatabaseID"] == db), key=lambda row: float(row["Confidence"]))
        axes[2].plot(x_positions, [float(row["UnsupportedClassifiedPer10k"]) for row in control_subset], marker="o", color=DB_COLORS[db])
    axes[0].set_title("A  Reference-aware mock recovery")
    axes[0].set_ylabel("Expected species recovered (%)")
    axes[1].set_title("B  Non-truth mock species burden")
    axes[1].set_ylabel("Species-level fragments per 10k")
    axes[2].set_title("C  Human-control unsupported burden")
    axes[2].set_ylabel("Classified fragments per 10k (symlog)")
    axes[2].set_yscale("symlog", linthresh=1)
    for axis in axes:
        axis.set_xlabel("Kraken2 confidence")
        axis.set_xticks(x_positions, x_labels)
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    axes[0].legend(loc="lower left")
    save_figure(fig, figure_dir, "17-confidence-tradeoff")


def plot_database_stability(
    classifications: list[dict[str, Any]],
    stability: list[dict[str, Any]],
    figure_dir: Path,
) -> None:
    set_publication_style()
    primary_mock = [row for row in classifications if row["PrimaryMatrix"] == "Yes" and row["ControlID"] == "mock"]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5), constrained_layout=True)
    x_positions = np.arange(len(CONFIDENCES))
    x_labels = ("0", "0.05", "0.10", "0.20", "0.50")
    for db in DB_ORDER:
        subset = sorted((row for row in primary_mock if row["DatabaseID"] == db), key=lambda row: float(row["Confidence"]))
        axes[0].plot(
            x_positions,
            [100 * float(row["ClassifiedFraction"]) for row in subset],
            marker="o",
            color=DB_COLORS[db],
            label=DB_LABELS[db],
        )
    axes[0].set_title("A  Mock classification yield")
    axes[0].set_ylabel("Classified paired fragments (%)")
    axes[0].set_xlabel("Kraken2 confidence")
    axes[0].set_xticks(x_positions, x_labels)
    axes[0].legend()
    axes[0].grid(axis="y", color="#DDDDDD", linewidth=0.6)

    pair_colors = {("standard8", "standard16"): "#7B3294", ("standard8", "pluspf8"): "#008837", ("standard16", "pluspf8"): "#C51B7D"}
    for pair, color in pair_colors.items():
        subset = sorted(
            (row for row in stability if row["ControlID"] == "mock" and (row["DatabaseA"], row["DatabaseB"]) == pair),
            key=lambda row: float(row["Confidence"]),
        )
        axes[1].plot(
            x_positions,
            [float(row["BrackenTotalVariation"]) for row in subset],
            marker="o",
            color=color,
            label=f"{DB_LABELS[pair[0]]} vs {DB_LABELS[pair[1]]}",
        )
    axes[1].set_title("B  Species-abundance disagreement")
    axes[1].set_ylabel("Bracken total variation")
    axes[1].set_xlabel("Kraken2 confidence")
    axes[1].set_xticks(x_positions, x_labels)
    axes[1].set_ylim(bottom=0)
    axes[1].text(
        3.5,
        0.045,
        "Not estimable\nwhen a branch is empty",
        ha="center",
        va="center",
        color="#666666",
        fontsize=8,
    )
    axes[1].legend(fontsize=7)
    axes[1].grid(axis="y", color="#DDDDDD", linewidth=0.6)
    save_figure(fig, figure_dir, "17-database-stability")


def plot_hit_group_control(hit_groups: list[dict[str, Any]], figure_dir: Path) -> None:
    set_publication_style()
    mock = sorted((row for row in hit_groups if row["ControlID"] == "mock"), key=lambda row: int(row["MinimumHitGroups"]))
    human = sorted((row for row in hit_groups if row["ControlID"] == "human"), key=lambda row: int(row["MinimumHitGroups"]))
    x = [int(row["MinimumHitGroups"]) for row in mock]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    axes[0].plot(x, [100 * float(row["ReferenceAwareSpeciesRecall"]) for row in mock], marker="o", color="#3B6FB6")
    axes[1].plot(x, [float(row["NonTruthSpeciesPer10k"]) for row in mock], marker="o", color="#D35F5F")
    axes[2].plot(x, [float(row["UnsupportedClassifiedPer10k"]) for row in human], marker="o", color="#2B8C7E")
    axes[0].set_ylim(0, 100)
    axes[1].set_ylim(0, 200)
    axes[2].set_ylim(0, 1)
    titles = ("A  Mock recovery", "B  Non-truth mock burden", "C  Human-control burden")
    ylabels = ("Expected species recovered (%)", "Species-level fragments per 10k", "Classified fragments per 10k")
    for axis, title, ylabel in zip(axes, titles, ylabels):
        axis.set_title(title)
        axis.set_xlabel("Minimum hit groups")
        axis.set_ylabel(ylabel)
        axis.set_xticks((1, 2, 3, 4))
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    fig.suptitle("Standard-16 at confidence 0.10", y=1.04, fontsize=11)
    save_figure(fig, figure_dir, "17-hit-group-control")


def validate_figures(figure_dir: Path, checks: Checks) -> None:
    for stem in FIGURE_STEMS:
        for suffix in (".pdf", ".png", ".tiff"):
            path = figure_dir / f"{stem}{suffix}"
            checks.add("figure", f"{stem}{suffix}-exists", path.is_file() and path.stat().st_size > 1_000, path.stat().st_size if path.exists() else "missing")
        png = Image.open(figure_dir / f"{stem}.png")
        checks.add("figure", f"{stem}-png-rgb", png.mode in {"RGB", "RGBA"}, png.mode)
        tiff = Image.open(figure_dir / f"{stem}.tiff")
        dpi = tiff.info.get("dpi", (0, 0))[0]
        checks.add("figure", f"{stem}-tiff-dpi", float(dpi) >= 349, dpi)


def host_path_hits(frozen_dir: Path, needles: list[str]) -> list[str]:
    hits: list[str] = []
    for path in frozen_dir.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES or path.name == "file-checksums.sha256":
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        for needle in needles:
            if needle and needle in text:
                hits.append(f"{path.relative_to(frozen_dir)}:{needle}")
    return hits


def write_aggregate_tables(
    frozen_dir: Path,
    crosswalk: list[dict[str, Any]],
    database_audit: list[dict[str, Any]],
    classifications: list[dict[str, Any]],
    truth_rows: list[dict[str, Any]],
    control_rows: list[dict[str, Any]],
    stability: list[dict[str, Any]],
    hit_groups: list[dict[str, Any]],
    resources: list[dict[str, Any]],
    source_audits: dict[str, Any],
) -> None:
    write_tsv(
        frozen_dir / "reference-crosswalk.tsv",
        crosswalk,
        [
            "DatabaseID", "ReleaseID", "TruthIndex", "ExpectedOrganismGTDBRS207",
            "StrainName", "AssemblyAccession", "NCBITaxID", "NCBIOrganismName",
            "SpeciesTaxID", "DatabaseSpeciesName", "ExpectedGenomePercent",
            "CandidateSequenceAccessions", "ExactMatchingSequences", "ReferenceStatus",
        ],
    )
    write_tsv(frozen_dir / "database-audit.tsv", database_audit, list(database_audit[0]))
    write_tsv(frozen_dir / "classification-summary.tsv", classifications, list(classifications[0]))
    write_tsv(frozen_dir / "truth-performance.tsv", truth_rows, list(truth_rows[0]))
    write_tsv(frozen_dir / "negative-control.tsv", control_rows, list(control_rows[0]))
    write_tsv(frozen_dir / "abundance-stability.tsv", stability, list(stability[0]))
    write_tsv(frozen_dir / "hit-group-sensitivity.tsv", hit_groups, list(hit_groups[0]))
    write_tsv(frozen_dir / "resource-usage.tsv", resources, list(resources[0]))
    source_rows: list[dict[str, Any]] = []
    for control, audit in source_audits.items():
        for mate in ("R1", "R2"):
            source_rows.append(
                {
                    "ControlID": control,
                    "Mate": mate,
                    "Pairs": audit["pairs"],
                    "Records": audit[mate]["records"],
                    "Bases": audit[mate]["bases"],
                    "CompressedBytes": audit[mate]["bytes"],
                    "CompressedSHA256": audit[mate]["sha256"],
                    "PairIDSHA256": audit["normalized_pair_id_sha256"],
                }
            )
    write_tsv(frozen_dir / "source-audit.tsv", source_rows, list(source_rows[0]))


def initialize_frozen(args: argparse.Namespace) -> None:
    project_root = args.project_root.resolve()
    prefix = args.environment_prefix.resolve()
    if not args.frozen_dir or not args.raw_dir or not args.work_dir or not args.database_root:
        raise SystemExit("Initialization requires frozen/raw/work/database paths")
    frozen_dir = args.frozen_dir.resolve()
    raw_dir = args.raw_dir.resolve()
    work_dir = args.work_dir.resolve()
    database_root = args.database_root.resolve()
    checks = Checks()

    static_contract_checks(project_root, checks)
    versions = environment_checks(project_root, prefix, checks)
    source_audits = audit_fastq_controls(project_root, checks)
    design = read_tsv(frozen_dir / "branch-design.tsv")
    primary = [row for row in design if row["PrimaryMatrix"] == "Yes"]
    hit_matrix = [row for row in design if row["HitGroupMatrix"] == "Yes"]
    checks.add("design", "unique-branches", len(design) == len({row["BranchID"] for row in design}) == 36, len(design))
    checks.add("design", "primary-factorial", len(primary) == 30, len(primary))
    checks.add("design", "primary-databases", {row["DatabaseID"] for row in primary} == set(DATABASES), sorted({row["DatabaseID"] for row in primary}))
    checks.add("design", "primary-controls", {row["ControlID"] for row in primary} == {"mock", "human"}, sorted({row["ControlID"] for row in primary}))
    checks.add("design", "primary-confidence-grid", {float(row["Confidence"]) for row in primary} == set(CONFIDENCES), sorted({float(row["Confidence"]) for row in primary}))
    checks.add("design", "primary-hit-groups-fixed", {int(row["MinimumHitGroups"]) for row in primary} == {2}, sorted({row["MinimumHitGroups"] for row in primary}))
    checks.add("design", "hit-group-grid", len(hit_matrix) == 8 and {int(row["MinimumHitGroups"]) for row in hit_matrix} == {1, 2, 3, 4}, sorted({row["MinimumHitGroups"] for row in hit_matrix}))
    checks.add("design", "hit-group-axis-fixed", {row["DatabaseID"] for row in hit_matrix} == {"standard16"} and {float(row["Confidence"]) for row in hit_matrix} == {0.10}, "Standard-16;confidence=0.10")
    checks.add("design", "read-length-by-control", all(int(row["BrackenReadLength"]) == (150 if row["ControlID"] == "mock" else 100) for row in design), "mock=150;human=100")

    crosswalk, database_audit, archive_summaries = build_reference_crosswalk(project_root, database_root, checks)
    taxonomy_parents, _ = parse_nodes(
        database_root / DATABASES["standard8"]["dir_rel"] / "nodes.dmp"
    )
    human_ancestors = ancestor_taxids(9606, taxonomy_parents)
    checks.add(
        "control",
        "human-ancestor-lineage",
        1 in human_ancestors and 9606 not in human_ancestors,
        len(human_ancestors),
    )
    classifications, truth_rows, control_rows, distributions = branch_metrics(
        design,
        frozen_dir,
        work_dir,
        crosswalk,
        human_ancestors,
        checks,
    )
    stability = abundance_stability(design, distributions)
    hit_groups = build_hit_group_table(design, classifications, truth_rows, control_rows)
    resources = resource_table(design, work_dir, checks)
    checks.add("analysis", "classification-rows", len(classifications) == 36, len(classifications))
    checks.add("analysis", "mock-truth-rows", len(truth_rows) == 18, len(truth_rows))
    checks.add("analysis", "human-control-rows", len(control_rows) == 18, len(control_rows))
    checks.add("analysis", "stability-rows", len(stability) == 30, len(stability))
    checks.add("analysis", "hit-group-rows", len(hit_groups) == 8, len(hit_groups))
    checks.add("analysis", "resource-rows", len(resources) == 72, len(resources))

    write_aggregate_tables(
        frozen_dir,
        crosswalk,
        database_audit,
        classifications,
        truth_rows,
        control_rows,
        stability,
        hit_groups,
        resources,
        source_audits,
    )

    baseline_branch = "standard16-mock-c010-h2"
    control_branch = "standard16-human-c010-h2"
    baseline_truth = next(row for row in truth_rows if row["BranchID"] == baseline_branch)
    baseline_classification = next(row for row in classifications if row["BranchID"] == baseline_branch)
    baseline_control = next(row for row in control_rows if row["BranchID"] == control_branch)
    truth_species = len({int(row["SpeciesTaxID"]) for row in crosswalk if row["DatabaseID"] == "standard16"})
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "evidence_date": "2026-07-21",
        "databases": 3,
        "controls": 2,
        "confidence_levels": 5,
        "primary_branches": 30,
        "unique_kraken_runs": 36,
        "hit_group_matrix_rows": 8,
        "truth_assemblies": EXPECTED_TRUTH_ASSEMBLIES,
        "truth_species": truth_species,
        "truth_abundance_percent": round(EXPECTED_TRUTH_SUM, 6),
        "kraken2_version": versions["Kraken2"],
        "bracken_cli_version": versions["BrackenCLI"],
        "bracken_package_version": versions["BrackenPackage"],
        "python_version": versions["Python"],
        "database_release_date": "2026-06-26",
        "database_archives": archive_summaries,
        "database_audit": database_audit,
        "baseline_branch": baseline_branch,
        "baseline_confidence": 0.10,
        "baseline_minimum_hit_groups": 2,
        "baseline_classified_fragments": baseline_classification["ClassifiedFragments"],
        "baseline_classified_fraction": baseline_classification["ClassifiedFraction"],
        "baseline_reference_aware_species_recall": baseline_truth["ReferenceAwareSpeciesRecall"],
        "baseline_nontruth_species_per10k": baseline_truth["NonTruthSpeciesPer10k"],
        "baseline_bracken_truth_fraction": baseline_truth["BrackenTruthFraction"],
        "baseline_human_unsupported_per10k": baseline_control["UnsupportedClassifiedPer10k"],
        "mock_input_pairs": CONTROL_EXPECTED["mock"]["pairs"],
        "human_control_pairs": CONTROL_EXPECTED["human"]["pairs"],
        "classification_unit": "paired_fragment",
        "bracken_rank": "S",
        "bracken_threshold": 10,
        "mock_bracken_read_length": 150,
        "human_bracken_read_length": 100,
        "raw_fastq_committed": False,
        "database_archive_committed": False,
        "database_index_committed": False,
        "per_fragment_output_committed": False,
        "human_control_is_extraction_blank": False,
        "human_control_boundary": "human-rich method control; unsupported-control calls are not an experimental false-positive rate",
        "universal_confidence_claimed": False,
        "one_time_network_access": True,
        "qa_network_access": False,
        "initialization_checks_passed": checks.passed,
        "initialization_checks_failed": checks.failed,
    }
    write_json(frozen_dir / "run-summary.json", summary)
    write_tsv(frozen_dir / "initialization-audit.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])

    hits = host_path_hits(
        frozen_dir,
        [str(project_root), str(prefix), str(raw_dir), str(work_dir), str(database_root), str(Path.home())],
    )
    if hits:
        raise RuntimeError(f"Host-specific paths remain in frozen Article 17 evidence: {hits[:8]}")
    checksum_entries = write_frozen_checksums(frozen_dir)
    _, checksum_failures = verify_frozen_checksums(frozen_dir)
    if checksum_failures:
        raise RuntimeError(f"Frozen checksum failures: {checksum_failures[:8]}")
    summary = json.loads((frozen_dir / "run-summary.json").read_text(encoding="utf-8"))
    summary["frozen_checksum_entries"] = checksum_entries
    summary["host_specific_paths_retained"] = False
    write_json(frozen_dir / "run-summary.json", summary)
    write_frozen_checksums(frozen_dir)
    if checks.failed:
        failed = [row for row in checks.rows if row["Status"] == "FAIL"]
        raise RuntimeError(f"Article 17 initialization failed: {failed[:10]}")


def copy_audit_tables(frozen_dir: Path, output_dir: Path) -> None:
    mapping = {
        "source-audit.tsv": "source-audit.tsv",
        "database-audit.tsv": "database-audit.tsv",
        "reference-crosswalk.tsv": "reference-audit.tsv",
        "classification-summary.tsv": "branch-audit.tsv",
        "truth-performance.tsv": "truth-audit.tsv",
        "negative-control.tsv": "control-audit.tsv",
        "abundance-stability.tsv": "stability-audit.tsv",
        "hit-group-sensitivity.tsv": "hit-group-audit.tsv",
        "resource-usage.tsv": "resource-audit.tsv",
    }
    for source, target in mapping.items():
        shutil.copyfile(frozen_dir / source, output_dir / target)


def routine_validation(args: argparse.Namespace) -> None:
    project_root = args.project_root.resolve()
    prefix = args.environment_prefix.resolve()
    frozen_dir = (args.frozen_dir or project_root / "data/small/17-kraken-database-confidence-frozen").resolve()
    output_dir = (args.output_dir or project_root / "results/17-kraken2-database-confidence").resolve()
    figure_dir = (args.figure_dir or project_root / "figures").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    static_contract_checks(project_root, checks)
    versions = environment_checks(project_root, prefix, checks)
    checksum_entries, checksum_failures = verify_frozen_checksums(frozen_dir)
    checks.add("frozen", "checksums", not checksum_failures, ";".join(checksum_failures) if checksum_failures else checksum_entries)
    summary = json.loads((frozen_dir / "run-summary.json").read_text(encoding="utf-8"))
    checks.add("frozen", "initialization-status", summary.get("status") == "passed", summary.get("status"))
    checks.add("frozen", "initialization-checks", summary.get("initialization_checks_failed") == 0, summary.get("initialization_checks_failed"))
    checks.add("frozen", "design-contract", summary.get("databases") == 3 and summary.get("controls") == 2 and summary.get("primary_branches") == 30 and summary.get("unique_kraken_runs") == 36, f"{summary.get('databases')};{summary.get('controls')};{summary.get('primary_branches')};{summary.get('unique_kraken_runs')}")
    checks.add("frozen", "truth-contract", summary.get("truth_assemblies") == 71 and abs(float(summary.get("truth_abundance_percent")) - 100.012157) < 1e-9, f"{summary.get('truth_assemblies')};{summary.get('truth_abundance_percent')}")
    checks.add("frozen", "network-free-qa", summary.get("qa_network_access") is False, summary.get("qa_network_access"))
    checks.add("frozen", "negative-control-boundary", summary.get("human_control_is_extraction_blank") is False and "not an experimental false-positive rate" in summary.get("human_control_boundary", ""), summary.get("human_control_boundary"))
    checks.add("frozen", "no-universal-threshold", summary.get("universal_confidence_claimed") is False, summary.get("universal_confidence_claimed"))
    checks.add("frozen", "large-inputs-excluded", all(summary.get(key) is False for key in ("raw_fastq_committed", "database_archive_committed", "database_index_committed", "per_fragment_output_committed")), "FASTQ/database/per-fragment excluded")
    checks.add("environment", "version-parity", versions["Kraken2"] == summary["kraken2_version"] and versions["BrackenCLI"] == summary["bracken_cli_version"] and versions["BrackenPackage"] == summary["bracken_package_version"], versions)

    database_audit = read_tsv(frozen_dir / "database-audit.tsv")
    crosswalk = read_tsv(frozen_dir / "reference-crosswalk.tsv")
    classifications = read_tsv(frozen_dir / "classification-summary.tsv")
    truth_rows = read_tsv(frozen_dir / "truth-performance.tsv")
    control_rows = read_tsv(frozen_dir / "negative-control.tsv")
    stability = read_tsv(frozen_dir / "abundance-stability.tsv")
    hit_groups = read_tsv(frozen_dir / "hit-group-sensitivity.tsv")
    checks.add("frozen", "aggregate-row-counts", len(database_audit) == 3 and len(crosswalk) == 213 and len(classifications) == 36 and len(truth_rows) == 18 and len(control_rows) == 18 and len(stability) == 30 and len(hit_groups) == 8, f"{len(database_audit)};{len(crosswalk)};{len(classifications)};{len(truth_rows)};{len(control_rows)};{len(stability)};{len(hit_groups)}")

    plot_reference_coverage(database_audit, figure_dir)
    plot_confidence_tradeoff(truth_rows, control_rows, figure_dir)
    plot_database_stability(classifications, stability, figure_dir)
    plot_hit_group_control(hit_groups, figure_dir)
    validate_figures(figure_dir, checks)
    copy_audit_tables(frozen_dir, output_dir)

    payload = {
        "status": "passed" if checks.failed == 0 else "failed",
        "databases": 3,
        "controls": 2,
        "confidence_levels": 5,
        "primary_branches": 30,
        "unique_kraken_runs": 36,
        "truth_assemblies": 71,
        "truth_species": summary["truth_species"],
        "truth_abundance_percent": 100.012157,
        "kraken2_version": versions["Kraken2"],
        "bracken_cli_version": versions["BrackenCLI"],
        "bracken_package_version": versions["BrackenPackage"],
        "baseline_branch": summary["baseline_branch"],
        "baseline_classified_fragments": summary["baseline_classified_fragments"],
        "baseline_classified_fraction": summary["baseline_classified_fraction"],
        "baseline_reference_aware_species_recall": summary["baseline_reference_aware_species_recall"],
        "baseline_nontruth_species_per10k": summary["baseline_nontruth_species_per10k"],
        "baseline_bracken_truth_fraction": summary["baseline_bracken_truth_fraction"],
        "baseline_human_unsupported_per10k": summary["baseline_human_unsupported_per10k"],
        "frozen_checksum_entries": checksum_entries,
        "checksum_failures": len(checksum_failures),
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
    }
    write_json(output_dir / "validation-summary.json", payload)
    write_tsv(output_dir / "validation-audit.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])
    (output_dir / "validation.log").write_text(
        f"Article 17 validation {payload['status']}: {checks.passed} passed, {checks.failed} failed\n",
        encoding="utf-8",
    )
    if checks.failed:
        failed = [row for row in checks.rows if row["Status"] == "FAIL"]
        raise RuntimeError(f"Article 17 routine validation failed: {failed[:10]}")


def main() -> None:
    args = parse_args()
    if args.initialize_frozen:
        initialize_frozen(args)
    else:
        routine_validation(args)


if __name__ == "__main__":
    main()
