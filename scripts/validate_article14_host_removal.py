#!/usr/bin/env python3
"""Validate Article 14 host removal, complexity, and duplicate evidence.

Initialization mode is used once after running the public FASTQ controls and
the complete Hostile Bowtie2 index. It audits Git-ignored FASTQs, freezes only
aggregate evidence, normalizes local paths, and writes checksums. Routine mode
is network-free and recreates the article tables and figures from frozen files.
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
from typing import Any, Iterator, TextIO

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article14-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


EXPECTED_VERSIONS = {
    "Hostile": "2.0.2",
    "Bowtie2": "2.5.4",
    "Samtools": "1.21",
    "fastp": "1.3.6",
    "SeqKit": "2.10.0",
    "Python": "3.11.15",
}
EXPECTED_ENV_YAML_SHA256 = (
    "e4cba4b08e702e94020f854287d023f47edc313b393631725f8829faed4bf3f3"
)
EXPECTED_ENV_LOCK_SHA256 = (
    "295a9351bba9e9f627bad21696ab2d0c028a41a3836d047a1a3dff29b2cc613a"
)
EXPECTED_SOURCE_MANIFEST_SHA256 = (
    "00c73ea625ded3e4079611c1d7862fb90afca7bea3b7a482e4350c53a38a6afa"
)
EXPECTED_NOTICE_SHA256 = (
    "86e59f2927f6d6d3e2703832eb4911b32a84ab9a8b37815c0035a1f2bf2affc7"
)
EXPECTED_INDEX_SHA256 = (
    "5b584f5c28abeec5dba78bd37b53fa476dd42af57051d2fb7d2f2098e3a2df13"
)
EXPECTED_INDEX_BYTES = 3_934_284_979
EXPECTED_INDEX_COMMIT = "e87e82caf5b34062db48154409658951fbd75c34"
EXPECTED_INDEX_FILES = (
    "human-t2t-hla.1.bt2",
    "human-t2t-hla.2.bt2",
    "human-t2t-hla.3.bt2",
    "human-t2t-hla.4.bt2",
    "human-t2t-hla.rev.1.bt2",
    "human-t2t-hla.rev.2.bt2",
)
EXPECTED_CONTROLS = {
    "human_positive": {
        "label": "Human control",
        "class": "Human",
        "project": "PRJEB3381",
        "sample": "SAMEA1573618",
        "run": "ERR194147",
        "pair_id_sha256": (
            "f6006da01154f735c5d80dce94be41711c58e33fd1c0b4290c7fbe734a8b9a18"
        ),
        "r1": "ERR194147_prefix20k_R1.fastq.gz",
        "r2": "ERR194147_prefix20k_R2.fastq.gz",
    },
    "mock_retention": {
        "label": "Microbial control",
        "class": "Microbial",
        "project": "PRJEB52977",
        "sample": "SAMEA14435832",
        "run": "ERR9765746",
        "pair_id_sha256": (
            "d8ccb5bd01cd5ab8205f08d56fda62cc9ab4d2e779faa6b250bceaa67c242033"
        ),
        "r1": "ERR9765746_prefix20k_R1.fastq.gz",
        "r2": "ERR9765746_prefix20k_R2.fastq.gz",
    },
}
CONTROL_FILE_STEMS = {
    "human_positive": "human",
    "mock_retention": "mock",
}
EXPECTED_ARCHIVE_MD5 = {
    ("human_positive", "R1"): "63634f9d0736e09debdfe5827a5e82a0",
    ("human_positive", "R2"): "81c3ba5eb6da709324cd124a847b4528",
    ("mock_retention", "R1"): "ed0e6e0ee846542531c742a45181cd6f",
    ("mock_retention", "R2"): "5b60ac93cb69dff77ae38cfa501afd06",
}
EXPECTED_ARCHIVE_BYTES = {
    ("human_positive", "R1"): 51_402_085_017,
    ("human_positive", "R2"): 52_551_746_373,
    ("mock_retention", "R1"): 1_740_647_656,
    ("mock_retention", "R2"): 2_104_551_765,
}
EXPECTED_LOCK_PREFIXES = {
    "hostile": "hostile-2.0.2-",
    "bowtie2": "bowtie2-2.5.4-",
    "samtools": "samtools-1.21-",
    "fastp": "fastp-1.3.6-",
    "seqkit": "seqkit-2.10.0-",
    "python": "python-3.11.15-",
    "matplotlib-base": "matplotlib-base-3.11.0-",
    "pillow": "pillow-12.3.0-",
}
COMPLEXITY_THRESHOLDS = (20, 30, 40)
FIGURE_STEMS = (
    "14-host-removal-recall-retention",
    "14-complexity-sensitivity",
    "14-duplicate-decision",
)
TEXT_SUFFIXES = {
    ".csv",
    ".json",
    ".log",
    ".sh",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--environment-prefix", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path)
    parser.add_argument("--raw-dir", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--index-archive", type=Path)
    parser.add_argument("--index-dir", type=Path)
    parser.add_argument("--initialize-frozen", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path,
    rows: list[dict[str, Any]],
    fields: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


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

    def add(self, check_id: str, passed: bool, detail: str) -> None:
        self.rows.append(
            {
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


def observed_tool_versions(
    environment_prefix: Path,
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    commands = {
        "Hostile": [str(environment_prefix / "bin/hostile"), "--version"],
        "Bowtie2": [str(environment_prefix / "bin/bowtie2"), "--version"],
        "Samtools": [str(environment_prefix / "bin/samtools"), "--version"],
        "fastp": [str(environment_prefix / "bin/fastp"), "--version"],
        "SeqKit": [str(environment_prefix / "bin/seqkit"), "version"],
        "Python": [str(environment_prefix / "bin/python"), "--version"],
    }
    patterns = {
        "Hostile": r"^(?:hostile )?([0-9.]+)$",
        "Bowtie2": r"version ([0-9.]+)",
        "Samtools": r"samtools ([0-9.]+)",
        "fastp": r"fastp ([0-9.]+)",
        "SeqKit": r"seqkit v([0-9.]+)",
        "Python": r"Python ([0-9.]+)",
    }
    versions: dict[str, str] = {}
    rows: list[dict[str, Any]] = []
    for tool, command in commands.items():
        return_code, output = command_output(command)
        match = re.search(patterns[tool], output)
        observed = match.group(1) if match else ""
        expected = EXPECTED_VERSIONS[tool]
        status = (
            "PASS"
            if return_code == 0 and observed == expected
            else "FAIL"
        )
        versions[tool] = observed
        rows.append(
            {
                "Tool": tool,
                "ExpectedVersion": expected,
                "ObservedVersion": observed,
                "ReturnCode": return_code,
                "Status": status,
            }
        )
    return versions, rows


def parse_explicit_lock(path: Path) -> tuple[int, dict[str, str]]:
    packages = [
        line.strip().rsplit("/", 1)[-1]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith("https://")
    ]
    targets: dict[str, str] = {}
    for name, prefix in EXPECTED_LOCK_PREFIXES.items():
        matches = [package for package in packages if package.startswith(prefix)]
        targets[name] = matches[0] if len(matches) == 1 else ""
    return len(packages), targets


def normalized_read_id(header: str) -> str:
    token = header[1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def next_fastq_record(
    handle: TextIO,
    path: Path,
    record_index: int,
) -> tuple[str, str, str, str] | None:
    first = handle.readline()
    if first == "":
        return None
    lines = [first, handle.readline(), handle.readline(), handle.readline()]
    if "" in lines[1:]:
        raise ValueError(f"Truncated FASTQ record {record_index} in {path}")
    header, sequence, plus, quality = [
        line.rstrip("\r\n") for line in lines
    ]
    if not header.startswith("@") or not plus.startswith("+"):
        raise ValueError(f"Malformed FASTQ record {record_index} in {path}")
    if len(sequence) != len(quality):
        raise ValueError(f"Sequence/quality length mismatch in {path}")
    return header, sequence, plus, quality


def iter_fastq_pairs(
    read1_path: Path,
    read2_path: Path,
) -> Iterator[tuple[str, str, str]]:
    with gzip.open(read1_path, "rt", encoding="ascii", newline="") as read1:
        with gzip.open(read2_path, "rt", encoding="ascii", newline="") as read2:
            count = 0
            while True:
                record1 = next_fastq_record(read1, read1_path, count + 1)
                record2 = next_fastq_record(read2, read2_path, count + 1)
                if record1 is None and record2 is None:
                    break
                if record1 is None or record2 is None:
                    raise ValueError("R1/R2 FASTQ record counts differ")
                id1 = normalized_read_id(record1[0])
                id2 = normalized_read_id(record2[0])
                if id1 != id2:
                    raise ValueError("R1/R2 FASTQ identifiers are not synchronized")
                count += 1
                yield id1, record1[1].upper(), record2[1].upper()


def transition_complexity(sequence: str) -> float:
    if len(sequence) < 2:
        return 0.0
    transitions = sum(
        left != right for left, right in zip(sequence, sequence[1:])
    )
    return 100.0 * transitions / (len(sequence) - 1)


def normalized_kmer_entropy(sequence: str, k: int = 2) -> float:
    if len(sequence) < k:
        return 0.0
    kmers = [
        sequence[index : index + k]
        for index in range(len(sequence) - k + 1)
        if set(sequence[index : index + k]) <= {"A", "C", "G", "T"}
    ]
    if not kmers:
        return 0.0
    counts = Counter(kmers)
    total = sum(counts.values())
    entropy = -sum(
        (count / total) * math.log2(count / total)
        for count in counts.values()
    )
    return entropy / math.log2(4**k)


def audit_fastq_pair(
    read1_path: Path,
    read2_path: Path,
    *,
    characterize: bool,
) -> dict[str, Any]:
    pair_id_digest = hashlib.sha256()
    pair_sequence_counter: Counter[str] = Counter()
    complexity_values = {"R1": [], "R2": []}
    entropy_values = {"R1": [], "R2": []}
    records = 0
    for read_id, sequence1, sequence2 in iter_fastq_pairs(
        read1_path, read2_path
    ):
        records += 1
        pair_id_digest.update(read_id.encode("ascii"))
        pair_id_digest.update(b"\n")
        if characterize:
            pair_digest = hashlib.sha256()
            pair_digest.update(sequence1.encode("ascii"))
            pair_digest.update(b"\0")
            pair_digest.update(sequence2.encode("ascii"))
            pair_sequence_counter[pair_digest.hexdigest()] += 1
            for mate, sequence in (("R1", sequence1), ("R2", sequence2)):
                complexity_values[mate].append(
                    transition_complexity(sequence)
                )
                entropy_values[mate].append(normalized_kmer_entropy(sequence))
    audit: dict[str, Any] = {
        "pairs": records,
        "mates_synchronized": True,
        "normalized_pair_id_sha256": pair_id_digest.hexdigest(),
        "read1_compressed_sha256": sha256(read1_path),
        "read2_compressed_sha256": sha256(read2_path),
    }
    if characterize:
        unique_pairs = len(pair_sequence_counter)
        duplicate_pairs = records - unique_pairs
        audit.update(
            {
                "unique_sequence_exact_pairs": unique_pairs,
                "sequence_exact_duplicate_pairs": duplicate_pairs,
                "sequence_exact_duplicate_fraction": (
                    duplicate_pairs / records if records else 0.0
                ),
                "maximum_pair_multiplicity": (
                    max(pair_sequence_counter.values())
                    if pair_sequence_counter
                    else 0
                ),
                "transition_complexity": {
                    mate: quantiles(values)
                    for mate, values in complexity_values.items()
                },
                "normalized_2mer_entropy": {
                    mate: quantiles(values)
                    for mate, values in entropy_values.items()
                },
                "complexity_pair_pass": {
                    str(threshold): sum(
                        complexity_values["R1"][index] >= threshold
                        and complexity_values["R2"][index] >= threshold
                        for index in range(records)
                    )
                    for threshold in COMPLEXITY_THRESHOLDS
                },
            }
        )
    return audit


def quantiles(values: list[float]) -> dict[str, float]:
    if not values:
        return {"minimum": 0.0, "q25": 0.0, "median": 0.0, "q75": 0.0, "maximum": 0.0}
    array = np.asarray(values, dtype=float)
    return {
        "minimum": float(np.min(array)),
        "q25": float(np.quantile(array, 0.25)),
        "median": float(np.median(array)),
        "q75": float(np.quantile(array, 0.75)),
        "maximum": float(np.max(array)),
    }


def hostile_record(frozen_dir: Path, stem: str) -> dict[str, Any]:
    payload = json.loads(
        (frozen_dir / f"hostile/{stem}-hostile.json").read_text(
            encoding="utf-8"
        )
    )
    if not isinstance(payload, list) or len(payload) != 1:
        raise ValueError(f"Expected one Hostile result for {stem}")
    return payload[0]


def resolve_initialization_path(path_value: str, work_dir: Path) -> Path:
    return Path(path_value.replace("${WORK_DIR}", str(work_dir)))


def fastp_payload(
    frozen_dir: Path,
    stem: str,
    branch: str,
) -> dict[str, Any]:
    return json.loads(
        (frozen_dir / f"fastp/{stem}-{branch}.json").read_text(
            encoding="utf-8"
        )
    )


def parse_resource_log(path: Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("Elapsed (wall clock) time"):
            values["elapsed"] = stripped.rsplit(": ", 1)[-1]
        elif stripped.startswith("Maximum resident set size (kbytes)"):
            values["maximum_rss_kb"] = int(stripped.rsplit(": ", 1)[-1])
        elif stripped.startswith("Percent of CPU this job got"):
            values["percent_cpu"] = stripped.rsplit(": ", 1)[-1]
    return values


def write_index_file_manifest(index_dir: Path, frozen_dir: Path) -> None:
    lines = []
    for name in EXPECTED_INDEX_FILES:
        path = index_dir / name
        lines.append(f"{sha256(path)}  {path.stat().st_size}  {name}")
    (frozen_dir / "index-files.sha256").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def parse_index_file_manifest(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        digest, size, name = line.split("  ", 2)
        rows.append(
            {
                "File": name,
                "SHA256": digest,
                "Bytes": int(size),
            }
        )
    return rows


def frozen_text_files(frozen_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in frozen_dir.rglob("*")
        if path.is_file()
        and path.name != "file-checksums.sha256"
        and path.suffix.lower() in TEXT_SUFFIXES
    )


def normalize_frozen_paths(
    project_root: Path,
    environment_prefix: Path,
    raw_dir: Path,
    work_dir: Path,
    index_archive: Path,
    index_dir: Path,
    frozen_dir: Path,
) -> None:
    replacements = (
        (str(work_dir), "${WORK_DIR}"),
        (str(raw_dir), "${RAW_DIR}"),
        (str(index_archive), "${HOSTILE_INDEX_ARCHIVE}"),
        (str(index_dir), "${HOSTILE_INDEX_DIR}"),
        (str(frozen_dir), "${FROZEN_DIR}"),
        (str(environment_prefix), "${HOST_REMOVAL_ENV_PREFIX}"),
        (str(project_root), "${PROJECT_ROOT}"),
    )
    for path in frozen_text_files(frozen_dir):
        text = path.read_text(encoding="utf-8")
        normalized = text
        for original, placeholder in replacements:
            normalized = normalized.replace(original, placeholder)
        normalized = re.sub(
            r"/tmp/tmp[A-Za-z0-9._-]+",
            "${TMPDIR}",
            normalized,
        )
        if normalized != text:
            path.write_text(normalized, encoding="utf-8")


def privacy_audit(
    project_root: Path,
    environment_prefix: Path,
    frozen_dir: Path,
    raw_dir: Path | None = None,
    work_dir: Path | None = None,
    index_dir: Path | None = None,
) -> tuple[list[dict[str, Any]], int]:
    text_files = frozen_text_files(frozen_dir)
    forbidden_literals = [str(project_root), str(environment_prefix)]
    for optional in (raw_dir, work_dir, index_dir):
        if optional is not None:
            forbidden_literals.append(str(optional))
    path_hits = []
    temporary_hits = []
    for path in text_files:
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(frozen_dir).as_posix()
        if any(literal in text for literal in forbidden_literals):
            path_hits.append(relative)
        if re.search(r"/tmp/tmp[A-Za-z0-9._-]+", text):
            temporary_hits.append(relative)
    frozen_fastq = [
        path.relative_to(frozen_dir).as_posix()
        for path in frozen_dir.rglob("*")
        if path.is_file()
        and (
            path.name.endswith(".fastq")
            or path.name.endswith(".fastq.gz")
            or path.name.endswith(".fq")
            or path.name.endswith(".fq.gz")
        )
    ]
    summary = json.loads(
        (frozen_dir / "controls-summary.json").read_text(encoding="utf-8")
    )
    unhashed_identifier_keys = []

    def walk(value: Any, key: str = "") -> None:
        if isinstance(value, dict):
            for child_key, child in value.items():
                if (
                    "read_id" in child_key.lower()
                    and not child_key.lower().endswith("sha256")
                ):
                    unhashed_identifier_keys.append(child_key)
                walk(child, child_key)
        elif isinstance(value, list):
            for child in value:
                walk(child, key)

    walk(summary)
    rows = [
        {
            "Boundary": "Frozen FASTQ files",
            "Observed": len(frozen_fastq),
            "Expected": 0,
            "Status": "PASS" if not frozen_fastq else "FAIL",
            "Interpretation": "Sequence-bearing FASTQs remain in controlled scratch",
        },
        {
            "Boundary": "Host-specific absolute paths",
            "Observed": len(set(path_hits + temporary_hits)),
            "Expected": 0,
            "Status": (
                "PASS" if not path_hits and not temporary_hits else "FAIL"
            ),
            "Interpretation": "Frozen text uses portable placeholders",
        },
        {
            "Boundary": "Unhashed raw read identifiers",
            "Observed": len(unhashed_identifier_keys),
            "Expected": 0,
            "Status": "PASS" if not unhashed_identifier_keys else "FAIL",
            "Interpretation": "Only identifier digests are retained",
        },
        {
            "Boundary": "Human benchmark interpretation",
            "Observed": "public_method_control_only",
            "Expected": "public_method_control_only",
            "Status": (
                "PASS"
                if summary["privacy_boundary"]
                == "aggregate_only_no_human_sequences_or_read_ids_frozen"
                else "FAIL"
            ),
            "Interpretation": "No phenotype or individual-genotype inference",
        },
    ]
    failures = sum(row["Status"] != "PASS" for row in rows)
    return rows, failures


def write_checksum_manifest(frozen_dir: Path) -> None:
    files = sorted(
        path
        for path in frozen_dir.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    lines = [
        f"{sha256(path)}  {path.relative_to(frozen_dir).as_posix()}"
        for path in files
    ]
    (frozen_dir / "file-checksums.sha256").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def verify_checksum_manifest(
    frozen_dir: Path,
) -> tuple[int, int, list[dict[str, Any]]]:
    expected = {}
    for line in (frozen_dir / "file-checksums.sha256").read_text(
        encoding="utf-8"
    ).splitlines():
        if line.strip():
            digest, relative = line.split("  ", 1)
            expected[relative] = digest
    actual = {
        path.relative_to(frozen_dir).as_posix()
        for path in frozen_dir.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    }
    rows = []
    failures = 0
    for relative in sorted(set(expected) | actual):
        path = frozen_dir / relative
        observed = sha256(path) if path.is_file() else ""
        status = (
            "PASS"
            if relative in expected
            and relative in actual
            and observed == expected[relative]
            else "FAIL"
        )
        failures += status != "PASS"
        rows.append(
            {
                "File": relative,
                "ExpectedSHA256": expected.get(relative, ""),
                "ObservedSHA256": observed,
                "Status": status,
            }
        )
    return len(rows), failures, rows


def check_project_contract(
    project_root: Path,
    environment_prefix: Path,
    frozen_dir: Path,
    checks: Checks,
) -> dict[str, Any]:
    versions, tool_rows = observed_tool_versions(environment_prefix)
    for row in tool_rows:
        checks.add(
            f"tool-{row['Tool'].lower()}",
            row["Status"] == "PASS",
            f"expected={row['ExpectedVersion']}; observed={row['ObservedVersion']}",
        )

    env_yaml = project_root / "env/host-removal.yml"
    env_lock = project_root / "env/host-removal-linux-64.lock"
    checks.add(
        "environment-yaml-sha256",
        sha256(env_yaml) == EXPECTED_ENV_YAML_SHA256,
        sha256(env_yaml),
    )
    checks.add(
        "environment-lock-sha256",
        sha256(env_lock) == EXPECTED_ENV_LOCK_SHA256,
        sha256(env_lock),
    )
    lock_count, lock_packages = parse_explicit_lock(env_lock)
    checks.add(
        "environment-lock-package-count",
        lock_count == 142,
        f"explicit package rows={lock_count}",
    )
    for package, prefix in EXPECTED_LOCK_PREFIXES.items():
        checks.add(
            f"lock-{package}",
            lock_packages[package].startswith(prefix),
            lock_packages[package],
        )

    source_manifest = project_root / "data/small/14-source-manifest.tsv"
    source_rows = read_tsv(source_manifest)
    checks.add(
        "source-manifest-sha256",
        sha256(source_manifest) == EXPECTED_SOURCE_MANIFEST_SHA256,
        sha256(source_manifest),
    )
    checks.add(
        "source-row-count",
        len(source_rows) == 4,
        f"rows={len(source_rows)}",
    )
    source_audit_rows = []
    for row in source_rows:
        control = row["ControlID"]
        mate = row["Mate"]
        expected = EXPECTED_CONTROLS.get(control, {})
        expected_md5 = EXPECTED_ARCHIVE_MD5.get((control, mate), "")
        expected_bytes = EXPECTED_ARCHIVE_BYTES.get((control, mate), 0)
        status = (
            expected
            and row["ExpectedClass"] == expected["class"]
            and row["ProjectAccession"] == expected["project"]
            and row["SampleAccession"] == expected["sample"]
            and row["RunAccession"] == expected["run"]
            and row["PrefixRecords"] == "20000"
            and row["Layout"] == "PAIRED"
            and row["ENAReportedMD5"] == expected_md5
            and int(row["ENABytes"]) == expected_bytes
            and row["HTTPSURL"].startswith("https://ftp.sra.ebi.ac.uk/")
        )
        checks.add(
            f"source-{control}-{mate.lower()}",
            bool(status),
            (
                f"{row['RunAccession']} {mate}; md5={row['ENAReportedMD5']}; "
                f"bytes={row['ENABytes']}"
            ),
        )
        source_audit_rows.append(
            {
                "Control": expected.get("label", control),
                "ExpectedClass": row["ExpectedClass"],
                "Project": row["ProjectAccession"],
                "Sample": row["SampleAccession"],
                "Run": row["RunAccession"],
                "Mate": mate,
                "PrefixPairs": row["PrefixRecords"],
                "ENAReportedMD5": row["ENAReportedMD5"],
                "ENABytes": row["ENABytes"],
                "Status": "PASS" if status else "FAIL",
            }
        )

    notice = project_root / "data/small/14-data-NOTICE.txt"
    checks.add(
        "data-notice-sha256",
        sha256(notice) == EXPECTED_NOTICE_SHA256,
        sha256(notice),
    )

    index_manifest_path = project_root / "data/small/14-index-manifest.tsv"
    index_rows = read_tsv(index_manifest_path)
    index_row = index_rows[0] if len(index_rows) == 1 else {}
    index_status = (
        len(index_rows) == 1
        and index_row.get("IndexName") == "human-t2t-hla"
        and index_row.get("Aligner") == "Bowtie2"
        and index_row.get("Release") == "2023-07"
        and index_row.get("ReferenceComponents")
        == "T2T-CHM13v2.0 + IPD-IMGT/HLA v3.51"
        and index_row.get("AssetSHA256") == EXPECTED_INDEX_SHA256
        and index_row.get("ArchiveBytes") == str(EXPECTED_INDEX_BYTES)
        and index_row.get("HostileTag") == "2.0.2"
        and index_row.get("HostileCommit") == EXPECTED_INDEX_COMMIT
        and index_row.get("DownloadStatus") == "VERIFIED_AND_EXTRACTED"
    )
    checks.add(
        "index-source-contract",
        index_status,
        (
            f"release={index_row.get('Release', '')}; "
            f"status={index_row.get('DownloadStatus', '')}"
        ),
    )
    index_file_rows = parse_index_file_manifest(
        frozen_dir / "index-files.sha256"
    )
    checks.add(
        "index-file-count",
        [row["File"] for row in index_file_rows]
        == list(EXPECTED_INDEX_FILES),
        f"files={len(index_file_rows)}",
    )
    checks.add(
        "index-files-nonempty",
        all(
            row["Bytes"] > 0
            and re.fullmatch(r"[0-9a-f]{64}", row["SHA256"]) is not None
            for row in index_file_rows
        ),
        f"total bytes={sum(row['Bytes'] for row in index_file_rows)}",
    )
    archive_rows = read_tsv(frozen_dir / "index-archive-audit.tsv")
    archive_row = archive_rows[0] if len(archive_rows) == 1 else {}
    archive_status = (
        len(archive_rows) == 1
        and archive_row.get("ExpectedBytes") == str(EXPECTED_INDEX_BYTES)
        and archive_row.get("ObservedBytes") == str(EXPECTED_INDEX_BYTES)
        and archive_row.get("ExpectedSHA256") == EXPECTED_INDEX_SHA256
        and archive_row.get("ObservedSHA256") == EXPECTED_INDEX_SHA256
        and archive_row.get("Status") == "VERIFIED"
    )
    checks.add(
        "index-archive-verification",
        archive_status,
        archive_row.get("ObservedSHA256", ""),
    )
    return {
        "versions": versions,
        "tool_rows": tool_rows,
        "lock_count": lock_count,
        "lock_packages": lock_packages,
        "source_rows": source_rows,
        "source_audit_rows": source_audit_rows,
        "index_row": index_row,
        "index_file_rows": index_file_rows,
        "archive_row": archive_row,
        "index_manifest_sha256": sha256(index_manifest_path),
    }


def validate_fastp_command(
    payload: dict[str, Any],
    branch: str,
) -> tuple[bool, str]:
    command = payload.get("command", "")
    required = (
        "--disable_adapter_trimming",
        "--disable_quality_filtering",
        "--disable_length_filtering",
        "--disable_trim_poly_g",
    )
    status = all(flag in command for flag in required)
    if branch == "baseline":
        status = status and "--low_complexity_filter" not in command
        status = status and "--dedup" not in command
    elif branch.startswith("complexity-"):
        threshold = branch.rsplit("-", 1)[-1]
        status = (
            status
            and "--low_complexity_filter" in command
            and (
                f"--complexity_threshold={threshold}" in command
                or f"--complexity_threshold {threshold}" in command
            )
            and "--dedup" not in command
        )
    elif branch == "dedup":
        status = (
            status
            and "--dedup" in command
            and (
                "--dup_calc_accuracy=1" in command
                or "--dup_calc_accuracy 1" in command
            )
            and "--low_complexity_filter" not in command
        )
    return status, command


def initialize_frozen(
    project_root: Path,
    environment_prefix: Path,
    frozen_dir: Path,
    raw_dir: Path,
    work_dir: Path,
    index_archive: Path,
    index_dir: Path,
) -> int:
    checks = Checks()
    write_index_file_manifest(index_dir, frozen_dir)
    contract = check_project_contract(
        project_root, environment_prefix, frozen_dir, checks
    )
    controls_summary = json.loads(
        (frozen_dir / "controls-summary.json").read_text(encoding="utf-8")
    )
    checks.add(
        "controls-summary-status",
        controls_summary["status"] == "passed",
        controls_summary["status"],
    )
    checks.add(
        "controls-total-pairs",
        controls_summary["total_control_pairs"] == 40000,
        str(controls_summary["total_control_pairs"]),
    )
    checks.add(
        "controls-archive-md5-boundary",
        controls_summary["complete_archive_md5_verified"] is False,
        "streamed prefixes do not claim complete-file MD5 verification",
    )

    control_evidence: dict[str, dict[str, Any]] = {}
    complexity_rows: list[dict[str, Any]] = []
    duplicate_rows: list[dict[str, Any]] = []
    host_rows: list[dict[str, Any]] = []

    for control_id, expected in EXPECTED_CONTROLS.items():
        stem = CONTROL_FILE_STEMS[control_id]
        input_r1 = raw_dir / expected["r1"]
        input_r2 = raw_dir / expected["r2"]
        input_audit = audit_fastq_pair(
            input_r1, input_r2, characterize=True
        )
        checks.add(
            f"{stem}-input-pairs",
            input_audit["pairs"] == 20000,
            f"pairs={input_audit['pairs']}",
        )
        checks.add(
            f"{stem}-input-sync",
            input_audit["normalized_pair_id_sha256"]
            == expected["pair_id_sha256"],
            input_audit["normalized_pair_id_sha256"],
        )
        summary_control = controls_summary["controls"][control_id]
        checks.add(
            f"{stem}-controls-summary-pairs",
            summary_control["input_pairs"] == input_audit["pairs"],
            str(summary_control["input_pairs"]),
        )
        checks.add(
            f"{stem}-controls-summary-id-digest",
            summary_control["normalized_pair_id_sha256"]
            == input_audit["normalized_pair_id_sha256"],
            summary_control["normalized_pair_id_sha256"],
        )

        hostile = hostile_record(frozen_dir, stem)
        output_r1 = resolve_initialization_path(
            hostile["fastq1_out_path"], work_dir
        )
        output_r2 = resolve_initialization_path(
            hostile["fastq2_out_path"], work_dir
        )
        output_audit = audit_fastq_pair(
            output_r1, output_r2, characterize=False
        )
        input_reads = int(hostile["reads_in"])
        output_reads = int(hostile["reads_out"])
        removed_reads = int(hostile["reads_removed"])
        input_pairs = input_reads // 2
        output_pairs = output_reads // 2
        removed_pairs = removed_reads // 2
        checks.add(
            f"{stem}-hostile-version",
            hostile["version"] == EXPECTED_VERSIONS["Hostile"],
            hostile["version"],
        )
        checks.add(
            f"{stem}-hostile-aligner",
            hostile["aligner"] == "bowtie2",
            hostile["aligner"],
        )
        checks.add(
            f"{stem}-hostile-index",
            hostile["index"] == "human-t2t-hla",
            hostile["index"],
        )
        checks.add(
            f"{stem}-hostile-options",
            set(hostile["options"]) == {"rename", "reorder"},
            ",".join(hostile["options"]),
        )
        checks.add(
            f"{stem}-hostile-read-ledger",
            input_reads == output_reads + removed_reads
            and input_reads == 40000
            and input_reads % 2 == 0
            and output_reads % 2 == 0,
            (
                f"in={input_reads}; out={output_reads}; "
                f"removed={removed_reads}"
            ),
        )
        checks.add(
            f"{stem}-hostile-output-pairs",
            output_audit["pairs"] == output_pairs,
            f"json={output_pairs}; FASTQ={output_audit['pairs']}",
        )
        if control_id == "human_positive":
            success_fraction = removed_pairs / input_pairs
            checks.add(
                "human-removal-threshold",
                success_fraction >= 0.99,
                f"removed fraction={success_fraction:.8f}",
            )
        else:
            success_fraction = output_pairs / input_pairs
            checks.add(
                "mock-retention-threshold",
                success_fraction >= 0.99,
                f"retained fraction={success_fraction:.8f}",
            )
        host_rows.append(
            {
                "Control": expected["label"],
                "Run": expected["run"],
                "InputPairs": input_pairs,
                "RetainedPairs": output_pairs,
                "RemovedPairs": removed_pairs,
                "RetainedFraction": output_pairs / input_pairs,
                "RemovedFraction": removed_pairs / input_pairs,
                "SuccessMetric": (
                    "Host pairs removed"
                    if control_id == "human_positive"
                    else "Microbial pairs retained"
                ),
                "SuccessFraction": success_fraction,
                "LedgerStatus": "PASS",
            }
        )

        baseline = fastp_payload(frozen_dir, stem, "baseline")
        baseline_status, _ = validate_fastp_command(baseline, "baseline")
        baseline_input_reads = int(
            baseline["summary"]["before_filtering"]["total_reads"]
        )
        baseline_output_reads = int(
            baseline["summary"]["after_filtering"]["total_reads"]
        )
        checks.add(
            f"{stem}-fastp-baseline-command",
            baseline_status,
            "diagnostic-only flags verified",
        )
        checks.add(
            f"{stem}-fastp-baseline-ledger",
            baseline_input_reads == 40000
            and baseline_output_reads == baseline_input_reads,
            f"in={baseline_input_reads}; out={baseline_output_reads}",
        )

        for threshold in COMPLEXITY_THRESHOLDS:
            branch = f"complexity-{threshold}"
            payload = fastp_payload(frozen_dir, stem, branch)
            branch_status, _ = validate_fastp_command(payload, branch)
            filtering = payload["filtering_result"]
            passed_reads = int(filtering["passed_filter_reads"])
            low_complexity_reads = int(filtering["low_complexity_reads"])
            output_path1 = (
                work_dir / f"fastp/{stem}-{branch}_R1.fastq.gz"
            )
            output_path2 = (
                work_dir / f"fastp/{stem}-{branch}_R2.fastq.gz"
            )
            branch_audit = audit_fastq_pair(
                output_path1, output_path2, characterize=False
            )
            retained_pairs = branch_audit["pairs"]
            independent_pairs = input_audit["complexity_pair_pass"][
                str(threshold)
            ]
            checks.add(
                f"{stem}-{branch}-command",
                branch_status,
                f"threshold={threshold}",
            )
            checks.add(
                f"{stem}-{branch}-output-ledger",
                retained_pairs * 2 == passed_reads
                and passed_reads + low_complexity_reads == 40000,
                (
                    f"pairs={retained_pairs}; passed reads={passed_reads}; "
                    f"low-complexity reads={low_complexity_reads}"
                ),
            )
            checks.add(
                f"{stem}-{branch}-independent-count",
                retained_pairs == independent_pairs,
                (
                    f"fastp={retained_pairs}; "
                    f"independent={independent_pairs}"
                ),
            )
            complexity_rows.append(
                {
                    "Control": expected["label"],
                    "Run": expected["run"],
                    "ThresholdPercent": threshold,
                    "InputPairs": input_audit["pairs"],
                    "RetainedPairs": retained_pairs,
                    "RemovedPairs": input_audit["pairs"] - retained_pairs,
                    "RetainedFraction": retained_pairs / input_audit["pairs"],
                    "IndependentRetainedPairs": independent_pairs,
                    "CountAgreement": (
                        "PASS"
                        if retained_pairs == independent_pairs
                        else "FAIL"
                    ),
                    "Interpretation": "Sensitivity branch; not the primary pipeline",
                }
            )

        dedup = fastp_payload(frozen_dir, stem, "dedup")
        dedup_status, _ = validate_fastp_command(dedup, "dedup")
        dedup_filtering = dedup["filtering_result"]
        dedup_output_r1 = work_dir / f"fastp/{stem}-dedup_R1.fastq.gz"
        dedup_output_r2 = work_dir / f"fastp/{stem}-dedup_R2.fastq.gz"
        dedup_audit = audit_fastq_pair(
            dedup_output_r1, dedup_output_r2, characterize=False
        )
        dedup_pairs = dedup_audit["pairs"]
        checks.add(
            f"{stem}-dedup-command",
            dedup_status,
            "fastp --dedup sensitivity branch only",
        )
        checks.add(
            f"{stem}-dedup-output-count",
            dedup_pairs * 2
            == int(dedup["summary"]["after_filtering"]["total_reads"]),
            (
                f"pairs={dedup_pairs}; output reads="
                f"{dedup['summary']['after_filtering']['total_reads']}; "
                f"filter-passed reads={dedup_filtering['passed_filter_reads']}"
            ),
        )
        checks.add(
            f"{stem}-dedup-bounded",
            0 <= dedup_pairs <= input_audit["pairs"],
            f"retained pairs={dedup_pairs}",
        )
        duplicate_rows.append(
            {
                "Control": expected["label"],
                "Run": expected["run"],
                "InputPairs": input_audit["pairs"],
                "SequenceExactUniquePairs": input_audit[
                    "unique_sequence_exact_pairs"
                ],
                "SequenceExactDuplicatePairs": input_audit[
                    "sequence_exact_duplicate_pairs"
                ],
                "SequenceExactDuplicateFraction": input_audit[
                    "sequence_exact_duplicate_fraction"
                ],
                "MaximumPairMultiplicity": input_audit[
                    "maximum_pair_multiplicity"
                ],
                "FastpEstimatedDuplicationRate": float(
                    baseline["duplication"]["rate"]
                ),
                "FastpDedupRetainedPairs": dedup_pairs,
                "FastpDedupRemovedPairs": input_audit["pairs"] - dedup_pairs,
                "PrimaryDeduplication": "false",
                "InferenceBoundary": (
                    "Sequence equality alone does not identify PCR or optical origin"
                ),
            }
        )
        control_evidence[control_id] = {
            "input_audit": input_audit,
            "output_audit": output_audit,
            "hostile": hostile,
            "host_row": host_rows[-1],
            "baseline_fastp": baseline,
            "dedup_fastp": dedup,
            "resource_summary": {
                "hostile": parse_resource_log(
                    frozen_dir / f"logs/hostile-{stem}.resources.txt"
                ),
                "fastp_baseline": parse_resource_log(
                    frozen_dir
                    / f"logs/fastp-{stem}-baseline.resources.txt"
                ),
            },
        }

    write_tsv(
        frozen_dir / "host-removal-ledger.tsv",
        host_rows,
        [
            "Control",
            "Run",
            "InputPairs",
            "RetainedPairs",
            "RemovedPairs",
            "RetainedFraction",
            "RemovedFraction",
            "SuccessMetric",
            "SuccessFraction",
            "LedgerStatus",
        ],
    )
    write_tsv(
        frozen_dir / "complexity-sensitivity.tsv",
        complexity_rows,
        [
            "Control",
            "Run",
            "ThresholdPercent",
            "InputPairs",
            "RetainedPairs",
            "RemovedPairs",
            "RetainedFraction",
            "IndependentRetainedPairs",
            "CountAgreement",
            "Interpretation",
        ],
    )
    write_tsv(
        frozen_dir / "duplicate-sensitivity.tsv",
        duplicate_rows,
        [
            "Control",
            "Run",
            "InputPairs",
            "SequenceExactUniquePairs",
            "SequenceExactDuplicatePairs",
            "SequenceExactDuplicateFraction",
            "MaximumPairMultiplicity",
            "FastpEstimatedDuplicationRate",
            "FastpDedupRetainedPairs",
            "FastpDedupRemovedPairs",
            "PrimaryDeduplication",
            "InferenceBoundary",
        ],
    )
    normalize_frozen_paths(
        project_root,
        environment_prefix,
        raw_dir,
        work_dir,
        index_archive,
        index_dir,
        frozen_dir,
    )
    privacy_rows, privacy_failures = privacy_audit(
        project_root,
        environment_prefix,
        frozen_dir,
        raw_dir,
        work_dir,
        index_dir,
    )
    checks.add(
        "privacy-boundaries",
        privacy_failures == 0,
        f"failures={privacy_failures}",
    )
    write_tsv(
        frozen_dir / "privacy-audit.tsv",
        privacy_rows,
        ["Boundary", "Observed", "Expected", "Status", "Interpretation"],
    )

    human = control_evidence["human_positive"]["host_row"]
    mock = control_evidence["mock_retention"]["host_row"]
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "evidence_date": "2026-07-20",
        "paired_policy": "retain_pair_only_if_both_mates_unmapped",
        "hostile_aligner": "Bowtie2",
        "hostile_index": "human-t2t-hla",
        "hostile_index_release": "2023-07",
        "hostile_index_components": (
            "T2T-CHM13v2.0 + IPD-IMGT/HLA v3.51"
        ),
        "hostile_index_tar_sha256": EXPECTED_INDEX_SHA256,
        "hostile_index_tar_bytes": EXPECTED_INDEX_BYTES,
        "hostile_index_files": len(EXPECTED_INDEX_FILES),
        "human_run_accession": EXPECTED_CONTROLS["human_positive"]["run"],
        "human_input_pairs": human["InputPairs"],
        "human_retained_pairs": human["RetainedPairs"],
        "human_removed_pairs": human["RemovedPairs"],
        "human_removed_fraction": human["RemovedFraction"],
        "mock_run_accession": EXPECTED_CONTROLS["mock_retention"]["run"],
        "mock_input_pairs": mock["InputPairs"],
        "mock_retained_pairs": mock["RetainedPairs"],
        "mock_removed_pairs": mock["RemovedPairs"],
        "mock_retained_fraction": mock["RetainedFraction"],
        "primary_low_complexity_filter": False,
        "complexity_sensitivity_thresholds_percent": list(
            COMPLEXITY_THRESHOLDS
        ),
        "primary_deduplication": False,
        "duplicate_evidence_level": "sequence_exact_pair_hash",
        "pcr_duplicate_origin_inferred": False,
        "optical_duplicate_origin_inferred": False,
        "umi_evidence_available": False,
        "raw_fastq_committed": False,
        "filtered_fastq_committed": False,
        "human_sequences_or_read_ids_frozen": False,
        "complete_source_archive_md5_verified": False,
        "one_time_network_access": True,
        "qa_network_access": False,
        "environment_yaml_sha256": sha256(
            project_root / "env/host-removal.yml"
        ),
        "environment_lock_sha256": sha256(
            project_root / "env/host-removal-linux-64.lock"
        ),
        "source_manifest_sha256": sha256(
            project_root / "data/small/14-source-manifest.tsv"
        ),
        "index_manifest_sha256": contract["index_manifest_sha256"],
        "data_notice_sha256": sha256(
            project_root / "data/small/14-data-NOTICE.txt"
        ),
        "tool_versions": contract["versions"],
        "environment_lock_packages": contract["lock_count"],
        "controls": {
            key: {
                "input_audit": value["input_audit"],
                "hostile_output_audit": value["output_audit"],
                "resource_summary": value["resource_summary"],
            }
            for key, value in control_evidence.items()
        },
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "checksum_failures": 0,
    }
    (frozen_dir / "run-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    write_tsv(
        frozen_dir / "initialization-audit.tsv",
        checks.rows,
        ["CheckID", "Status", "Detail"],
    )
    if checks.failed:
        print(json.dumps(summary, indent=2))
        return 1
    write_checksum_manifest(frozen_dir)
    print(json.dumps(summary, indent=2))
    return 0


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "legend.fontsize": 8.5,
            "figure.dpi": 120,
            "savefig.dpi": 350,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = figure_dir / f"{stem}.pdf"
    png_path = figure_dir / f"{stem}.png"
    tiff_path = figure_dir / f"{stem}.tiff"
    fig.savefig(
        pdf_path,
        bbox_inches="tight",
        metadata={
            "Creator": "metagenomics-best-practices",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    fig.savefig(
        png_path,
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": "metagenomics-best-practices"},
    )
    with Image.open(png_path) as image:
        image.convert("RGB").save(
            tiff_path,
            compression="tiff_lzw",
            dpi=(350, 350),
        )
    plt.close(fig)


def plot_host_removal(rows: list[dict[str, str]], figure_dir: Path) -> None:
    labels = [
        f"{row['Control']}\n{row['Run']} · {int(row['InputPairs']):,} pairs"
        for row in rows
    ]
    success = [100 * float(row["SuccessFraction"]) for row in rows]
    errors = [100 - value for value in success]
    colors = ["#D55E00", "#009E73"]
    fig, axis = plt.subplots(figsize=(7.4, 4.4))
    bars = axis.bar(labels, success, color=colors, width=0.58)
    axis.set_ylim(0, 103)
    axis.set_ylabel("Control-specific success (%)")
    axis.set_title(
        "Host removal requires both recall and microbial retention",
        pad=12,
        fontweight="bold",
    )
    axis.grid(axis="y", color="#E6E6E6", linewidth=0.7)
    for index, (bar, value, error, row) in enumerate(
        zip(bars, success, errors, rows)
    ):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.8,
            f"{value:.3f}%",
            ha="center",
            va="bottom",
            fontweight="bold",
        )
        failure_label = (
            "Host pairs retained"
            if index == 0
            else "Microbial pairs removed"
        )
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            3.5,
            f"{failure_label}\n{error:.3f}%",
            ha="center",
            va="bottom",
            color="white",
            fontsize=8.2,
            fontweight="bold",
        )
    axis.tick_params(axis="x", labelsize=8.8, pad=8)
    fig.subplots_adjust(bottom=0.18)
    save_figure(fig, figure_dir, FIGURE_STEMS[0])


def plot_complexity(rows: list[dict[str, str]], figure_dir: Path) -> None:
    fig, axis = plt.subplots(figsize=(7.6, 4.5))
    style = {
        "Human control": ("#D55E00", "o"),
        "Microbial control": ("#0072B2", "s"),
    }
    for control in ("Human control", "Microbial control"):
        selected = [row for row in rows if row["Control"] == control]
        thresholds = [int(row["ThresholdPercent"]) for row in selected]
        retained = [100 * float(row["RetainedFraction"]) for row in selected]
        color, marker = style[control]
        axis.plot(
            thresholds,
            retained,
            color=color,
            marker=marker,
            linewidth=2,
            markersize=6,
            label=control,
        )
        for x, y, row in zip(thresholds, retained, selected):
            removed = int(row["RemovedPairs"])
            offset = 9 if control == "Microbial control" else -15
            vertical_alignment = "bottom" if offset > 0 else "top"
            axis.annotate(
                f"−{removed}",
                (x, y),
                xytext=(0, offset),
                textcoords="offset points",
                ha="center",
                va=vertical_alignment,
                color=color,
                fontsize=8,
            )
    minimum = min(
        100 * float(row["RetainedFraction"]) for row in rows
    )
    axis.set_ylim(max(0, minimum - 2.5), 100.8)
    axis.set_xticks(list(COMPLEXITY_THRESHOLDS))
    axis.set_xlabel("Transition-complexity threshold (%)")
    axis.set_ylabel("Read pairs retained (%)")
    axis.set_title(
        "Low-complexity filtering is a threshold sensitivity analysis",
        pad=12,
        fontweight="bold",
    )
    axis.grid(color="#E6E6E6", linewidth=0.7)
    axis.legend(frameon=False, loc="lower left")
    axis.text(
        0.99,
        0.03,
        "Labels show removed pairs · Primary pipeline: filter disabled",
        transform=axis.transAxes,
        ha="right",
        color="#555555",
        fontsize=8.2,
    )
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[1])


def plot_duplicate_decision(
    rows: list[dict[str, str]],
    figure_dir: Path,
) -> None:
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(10.0, 4.6),
        gridspec_kw={"width_ratios": [1.1, 1.35]},
    )
    labels = [row["Control"] for row in rows]
    exact = [
        100 * float(row["SequenceExactDuplicateFraction"]) for row in rows
    ]
    fastp = [
        100 * float(row["FastpEstimatedDuplicationRate"]) for row in rows
    ]
    x = np.arange(len(labels))
    width = 0.34
    axes[0].bar(
        x - width / 2,
        exact,
        width,
        color="#0072B2",
        label="Exact pair hash",
    )
    axes[0].bar(
        x + width / 2,
        fastp,
        width,
        color="#E69F00",
        label="fastp estimate",
    )
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels)
    axes[0].set_ylabel("Duplicate signal (%)")
    axes[0].set_title("Observed sequence redundancy")
    axes[0].grid(axis="y", color="#E6E6E6", linewidth=0.7)
    axes[0].legend(frameon=False, fontsize=8)
    for position, value in zip(x - width / 2, exact):
        axes[0].text(
            position,
            value + max(exact + fastp + [1]) * 0.025,
            f"{value:.2f}%",
            ha="center",
            fontsize=8,
        )
    for position, value in zip(x + width / 2, fastp):
        axes[0].text(
            position,
            value + max(exact + fastp + [1]) * 0.025,
            f"{value:.2f}%",
            ha="center",
            fontsize=8,
        )

    axes[1].axis("off")
    axes[1].set_title(
        "What the evidence can support",
        loc="left",
        pad=8,
        fontweight="bold",
    )
    decision_rows = [
        ("Identical paired sequence", "Measured", "#0072B2"),
        ("PCR duplicate origin", "Not identifiable", "#D55E00"),
        ("Optical duplicate origin", "No coordinate evidence", "#D55E00"),
        ("UMI molecule identity", "No UMI evidence", "#D55E00"),
        ("Primary abundance table", "Keep reads; audit sensitivity", "#009E73"),
    ]
    y = 0.88
    for evidence, decision, color in decision_rows:
        axes[1].add_patch(
            matplotlib.patches.FancyBboxPatch(
                (0.02, y - 0.08),
                0.96,
                0.12,
                boxstyle="round,pad=0.01",
                facecolor="#F5F5F5",
                edgecolor="#DDDDDD",
                transform=axes[1].transAxes,
            )
        )
        axes[1].text(
            0.05,
            y,
            evidence,
            transform=axes[1].transAxes,
            va="center",
            fontsize=8.8,
        )
        axes[1].text(
            0.95,
            y,
            decision,
            transform=axes[1].transAxes,
            ha="right",
            va="center",
            color=color,
            fontweight="bold",
            fontsize=8.8,
        )
        y -= 0.17
    fig.suptitle(
        "Duplicate signals do not establish duplicate origin",
        y=1.01,
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[2])


def routine_validation(
    project_root: Path,
    environment_prefix: Path,
    frozen_dir: Path,
    output_dir: Path,
    figure_dir: Path,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    checks = Checks()
    checksum_count, checksum_failures, checksum_rows = (
        verify_checksum_manifest(frozen_dir)
    )
    checks.add(
        "frozen-checksums",
        checksum_failures == 0,
        f"files={checksum_count}; failures={checksum_failures}",
    )
    contract = check_project_contract(
        project_root, environment_prefix, frozen_dir, checks
    )
    privacy_rows, privacy_failures = privacy_audit(
        project_root,
        environment_prefix,
        frozen_dir,
    )
    checks.add(
        "privacy-boundaries",
        privacy_failures == 0,
        f"failures={privacy_failures}",
    )
    frozen_summary = json.loads(
        (frozen_dir / "run-summary.json").read_text(encoding="utf-8")
    )
    checks.add(
        "frozen-run-status",
        frozen_summary["status"] == "passed",
        frozen_summary["status"],
    )
    checks.add(
        "human-input-pairs",
        frozen_summary["human_input_pairs"] == 20000,
        str(frozen_summary["human_input_pairs"]),
    )
    checks.add(
        "human-removal-threshold",
        frozen_summary["human_removed_fraction"] >= 0.99,
        f"{frozen_summary['human_removed_fraction']:.8f}",
    )
    checks.add(
        "mock-input-pairs",
        frozen_summary["mock_input_pairs"] == 20000,
        str(frozen_summary["mock_input_pairs"]),
    )
    checks.add(
        "mock-retention-threshold",
        frozen_summary["mock_retained_fraction"] >= 0.99,
        f"{frozen_summary['mock_retained_fraction']:.8f}",
    )
    checks.add(
        "paired-removal-policy",
        frozen_summary["paired_policy"]
        == "retain_pair_only_if_both_mates_unmapped",
        frozen_summary["paired_policy"],
    )
    checks.add(
        "primary-complexity-disabled",
        frozen_summary["primary_low_complexity_filter"] is False,
        "primary_low_complexity_filter=false",
    )
    checks.add(
        "primary-dedup-disabled",
        frozen_summary["primary_deduplication"] is False,
        "primary_deduplication=false",
    )
    checks.add(
        "routine-no-network",
        frozen_summary["qa_network_access"] is False,
        "qa_network_access=false",
    )
    checks.add(
        "index-tar-sha256",
        frozen_summary["hostile_index_tar_sha256"]
        == EXPECTED_INDEX_SHA256,
        frozen_summary["hostile_index_tar_sha256"],
    )

    host_rows = read_tsv(frozen_dir / "host-removal-ledger.tsv")
    complexity_rows = read_tsv(
        frozen_dir / "complexity-sensitivity.tsv"
    )
    duplicate_rows = read_tsv(
        frozen_dir / "duplicate-sensitivity.tsv"
    )
    checks.add(
        "host-ledger-row-count",
        len(host_rows) == 2
        and all(row["LedgerStatus"] == "PASS" for row in host_rows),
        f"rows={len(host_rows)}",
    )
    checks.add(
        "complexity-row-count",
        len(complexity_rows) == 6
        and all(row["CountAgreement"] == "PASS" for row in complexity_rows),
        f"rows={len(complexity_rows)}",
    )
    checks.add(
        "duplicate-row-count",
        len(duplicate_rows) == 2
        and all(row["PrimaryDeduplication"] == "false" for row in duplicate_rows),
        f"rows={len(duplicate_rows)}",
    )

    configure_plotting()
    plot_host_removal(host_rows, figure_dir)
    plot_complexity(complexity_rows, figure_dir)
    plot_duplicate_decision(duplicate_rows, figure_dir)
    for stem in FIGURE_STEMS:
        for suffix in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            checks.add(
                f"figure-{stem}-{suffix}",
                path.is_file() and path.stat().st_size > 0,
                f"{path.name} bytes={path.stat().st_size if path.exists() else 0}",
            )

    index_audit_rows = [
        {
            "Item": "Official archive",
            "Release": contract["index_row"].get("Release", ""),
            "Reference": contract["index_row"].get(
                "ReferenceComponents", ""
            ),
            "Bytes": EXPECTED_INDEX_BYTES,
            "SHA256": EXPECTED_INDEX_SHA256,
            "Status": "PASS",
        }
    ]
    index_audit_rows.extend(
        {
            "Item": row["File"],
            "Release": "2023-07",
            "Reference": "human-t2t-hla Bowtie2 component",
            "Bytes": row["Bytes"],
            "SHA256": row["SHA256"],
            "Status": "PASS",
        }
        for row in contract["index_file_rows"]
    )
    write_tsv(
        output_dir / "tool-audit.tsv",
        contract["tool_rows"],
        [
            "Tool",
            "ExpectedVersion",
            "ObservedVersion",
            "ReturnCode",
            "Status",
        ],
    )
    write_tsv(
        output_dir / "source-audit.tsv",
        contract["source_audit_rows"],
        [
            "Control",
            "ExpectedClass",
            "Project",
            "Sample",
            "Run",
            "Mate",
            "PrefixPairs",
            "ENAReportedMD5",
            "ENABytes",
            "Status",
        ],
    )
    write_tsv(
        output_dir / "index-audit.tsv",
        index_audit_rows,
        ["Item", "Release", "Reference", "Bytes", "SHA256", "Status"],
    )
    write_tsv(
        output_dir / "host-removal-audit.tsv",
        host_rows,
        list(host_rows[0].keys()),
    )
    write_tsv(
        output_dir / "complexity-audit.tsv",
        complexity_rows,
        list(complexity_rows[0].keys()),
    )
    write_tsv(
        output_dir / "duplicate-audit.tsv",
        duplicate_rows,
        list(duplicate_rows[0].keys()),
    )
    write_tsv(
        output_dir / "privacy-audit.tsv",
        privacy_rows,
        ["Boundary", "Observed", "Expected", "Status", "Interpretation"],
    )
    write_tsv(
        output_dir / "checksum-audit.tsv",
        checksum_rows,
        ["File", "ExpectedSHA256", "ObservedSHA256", "Status"],
    )

    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "paired_policy": frozen_summary["paired_policy"],
        "human_run_accession": frozen_summary["human_run_accession"],
        "human_input_pairs": frozen_summary["human_input_pairs"],
        "human_retained_pairs": frozen_summary["human_retained_pairs"],
        "human_removed_pairs": frozen_summary["human_removed_pairs"],
        "human_removed_fraction": frozen_summary[
            "human_removed_fraction"
        ],
        "mock_run_accession": frozen_summary["mock_run_accession"],
        "mock_input_pairs": frozen_summary["mock_input_pairs"],
        "mock_retained_pairs": frozen_summary["mock_retained_pairs"],
        "mock_removed_pairs": frozen_summary["mock_removed_pairs"],
        "mock_retained_fraction": frozen_summary[
            "mock_retained_fraction"
        ],
        "hostile_index": frozen_summary["hostile_index"],
        "hostile_index_release": frozen_summary[
            "hostile_index_release"
        ],
        "hostile_index_tar_sha256": EXPECTED_INDEX_SHA256,
        "primary_low_complexity_filter": False,
        "complexity_sensitivity_thresholds_percent": list(
            COMPLEXITY_THRESHOLDS
        ),
        "primary_deduplication": False,
        "duplicate_evidence_level": "sequence_exact_pair_hash",
        "pcr_duplicate_origin_inferred": False,
        "optical_duplicate_origin_inferred": False,
        "umi_evidence_available": False,
        "raw_fastq_committed": False,
        "filtered_fastq_committed": False,
        "qa_network_access": False,
        "environment_lock_packages": contract["lock_count"],
        "checksum_files": checksum_count,
        "checksum_failures": checksum_failures,
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
    }
    (output_dir / "validation-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    write_tsv(
        output_dir / "validation.log",
        checks.rows,
        ["CheckID", "Status", "Detail"],
    )
    print(json.dumps(summary, indent=2))
    return 0 if checks.failed == 0 else 1


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    environment_prefix = args.environment_prefix.resolve()
    frozen_dir = (
        args.frozen_dir.resolve()
        if args.frozen_dir
        else project_root / "data/small/14-host-removal-frozen"
    )
    if args.initialize_frozen:
        required = {
            "--raw-dir": args.raw_dir,
            "--work-dir": args.work_dir,
            "--index-archive": args.index_archive,
            "--index-dir": args.index_dir,
        }
        missing = [name for name, value in required.items() if value is None]
        if missing:
            raise SystemExit(
                "Initialization requires " + ", ".join(missing)
            )
        return initialize_frozen(
            project_root,
            environment_prefix,
            frozen_dir,
            args.raw_dir.resolve(),
            args.work_dir.resolve(),
            args.index_archive.resolve(),
            args.index_dir.resolve(),
        )
    if args.output_dir is None or args.figure_dir is None:
        raise SystemExit(
            "--output-dir and --figure-dir are required for routine validation"
        )
    return routine_validation(
        project_root,
        environment_prefix,
        frozen_dir,
        args.output_dir.resolve(),
        args.figure_dir.resolve(),
    )


if __name__ == "__main__":
    raise SystemExit(main())
