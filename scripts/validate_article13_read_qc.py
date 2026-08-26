#!/usr/bin/env python3
"""Validate frozen Article 13 FastQC, fastp, and MultiQC evidence.

Initialization mode is called once after the real FASTQ run. It audits the
Git-ignored raw and cleaned FASTQ files, writes a compact run summary, and
freezes checksums. Routine QA is network-free and reads only committed reports.
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
import sys
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterator, TextIO

os.environ.setdefault(
    "MPLCONFIGDIR",
    "/tmp/metagenome-article13-matplotlib",
)
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


EXPECTED_VERSIONS = {
    "FastQC": "0.12.1",
    "fastp": "1.3.6",
    "MultiQC": "1.35",
    "Python": "3.14.6",
}
EXPECTED_SOURCE = {
    "ProjectAccession": "PRJEB52977",
    "SampleAccession": "SAMEA14435832",
    "RunAccession": "ERR9765746",
    "PrefixRecords": "100000",
}
EXPECTED_FASTP_FLAGS = (
    "--thread 4",
    "--compression 6",
    "--detect_adapter_for_pe",
    "--qualified_quality_phred 20",
    "--unqualified_percent_limit 40",
    "--n_base_limit 5",
    "--length_required 50",
    "--disable_trim_poly_g",
    "--overrepresentation_analysis",
    "--overrepresentation_sampling 20",
)
PROHIBITED_FASTP_FLAGS = (
    "--cut_front",
    "--cut_tail",
    "--cut_right",
    "--correction",
    "--dedup",
    "--low_complexity_filter",
)
FASTP_FILTER_KEYS = (
    "passed_filter_reads",
    "low_quality_reads",
    "too_many_N_reads",
    "adapter_dimer_reads",
    "too_short_reads",
    "too_long_reads",
)
FIGURE_STEM_QUALITY = "13-per-cycle-quality"
FIGURE_STEM_FATE = "13-read-pair-fate"
FIGURE_STEM_MODULES = "13-fastqc-module-states"
FROZEN_TEXT_SUFFIXES = {
    ".csv",
    ".html",
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
    parser.add_argument("--initialize-frozen", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(
    path: Path,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def command_output(command: list[str]) -> tuple[int, str]:
    process = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env={
            **os.environ,
            "LC_ALL": "C",
            "TZ": "UTC",
            "MULTIQC_NO_VERSION_CHECK": "1",
        },
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
        "FastQC": [str(environment_prefix / "bin/fastqc"), "--version"],
        "fastp": [str(environment_prefix / "bin/fastp"), "--version"],
        "MultiQC": [str(environment_prefix / "bin/multiqc"), "--version"],
        "Python": [str(environment_prefix / "bin/python"), "--version"],
    }
    patterns = {
        "FastQC": r"FastQC v([0-9.]+)",
        "fastp": r"fastp ([0-9.]+)",
        "MultiQC": r"multiqc, version ([0-9.]+)",
        "Python": r"Python ([0-9.]+)",
    }
    versions: dict[str, str] = {}
    rows: list[dict[str, Any]] = []
    for tool, command in commands.items():
        return_code, output = command_output(command)
        match = re.search(patterns[tool], output)
        version = match.group(1) if match else ""
        versions[tool] = version
        rows.append(
            {
                "Tool": tool,
                "ExpectedVersion": EXPECTED_VERSIONS[tool],
                "ObservedVersion": version,
                "ReturnCode": return_code,
                "Status": (
                    "PASS"
                    if return_code == 0 and version == EXPECTED_VERSIONS[tool]
                    else "FAIL"
                ),
            }
        )
    return versions, rows


def normalized_read_id(header: str) -> str:
    token = header[1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def next_fastq_record(
    handle: TextIO,
    source: Path,
    record_index: int,
) -> tuple[str, str, str, str] | None:
    header = handle.readline()
    if header == "":
        return None
    lines = [header, handle.readline(), handle.readline(), handle.readline()]
    if "" in lines[1:]:
        raise ValueError(f"Truncated record {record_index} in {source}")
    values = [value.rstrip("\r\n") for value in lines]
    header_text, sequence, plus, quality = values
    if not header_text.startswith("@"):
        raise ValueError(f"Bad header at record {record_index} in {source}")
    if not plus.startswith("+"):
        raise ValueError(f"Bad separator at record {record_index} in {source}")
    if len(sequence) != len(quality):
        raise ValueError(
            f"Sequence/quality length mismatch at record {record_index} in {source}"
        )
    return header_text, sequence, plus, quality


def audit_fastq_pair(
    read1_path: Path,
    read2_path: Path,
) -> dict[str, Any]:
    digests = {"R1": hashlib.sha256(), "R2": hashlib.sha256()}
    lengths: dict[str, list[int]] = {"R1": [], "R2": []}
    total_bases = {"R1": 0, "R2": 0}
    pair_id_digest = hashlib.sha256()
    records = 0
    first_id = ""
    last_id = ""
    with gzip.open(read1_path, "rt", encoding="ascii", newline="") as read1:
        with gzip.open(read2_path, "rt", encoding="ascii", newline="") as read2:
            while True:
                record_index = records + 1
                record1 = next_fastq_record(read1, read1_path, record_index)
                record2 = next_fastq_record(read2, read2_path, record_index)
                if record1 is None and record2 is None:
                    break
                if record1 is None or record2 is None:
                    raise ValueError("R1/R2 record counts differ")
                id1 = normalized_read_id(record1[0])
                id2 = normalized_read_id(record2[0])
                if id1 != id2:
                    raise ValueError(
                        f"R1/R2 IDs differ at record {record_index}: {id1} != {id2}"
                    )
                records += 1
                if records == 1:
                    first_id = id1
                last_id = id1
                pair_id_digest.update(id1.encode("ascii"))
                pair_id_digest.update(b"\n")
                for mate, record in (("R1", record1), ("R2", record2)):
                    serialized = (
                        f"{record[0]}\n{record[1]}\n{record[2]}\n{record[3]}\n"
                    ).encode("ascii")
                    digests[mate].update(serialized)
                    lengths[mate].append(len(record[1]))
                    total_bases[mate] += len(record[1])
    return {
        "pairs": records,
        "mates_synchronized": True,
        "first_read_id": first_id,
        "last_read_id": last_id,
        "normalized_pair_id_sha256": pair_id_digest.hexdigest(),
        "R1": {
            "records": records,
            "total_bases": total_bases["R1"],
            "minimum_read_length": min(lengths["R1"]),
            "maximum_read_length": max(lengths["R1"]),
            "uncompressed_fastq_sha256": digests["R1"].hexdigest(),
            "compressed_fastq_sha256": file_sha256(read1_path),
            "compressed_fastq_bytes": read1_path.stat().st_size,
        },
        "R2": {
            "records": records,
            "total_bases": total_bases["R2"],
            "minimum_read_length": min(lengths["R2"]),
            "maximum_read_length": max(lengths["R2"]),
            "uncompressed_fastq_sha256": digests["R2"].hexdigest(),
            "compressed_fastq_sha256": file_sha256(read2_path),
            "compressed_fastq_bytes": read2_path.stat().st_size,
        },
    }


def parse_fastqc_zip(path: Path) -> dict[str, Any]:
    with zipfile.ZipFile(path) as archive:
        candidates = [
            name
            for name in archive.namelist()
            if name.endswith("/fastqc_data.txt")
        ]
        if len(candidates) != 1:
            raise ValueError(f"Expected one fastqc_data.txt in {path}")
        text = archive.read(candidates[0]).decode("utf-8")
    version_match = re.search(r"^##FastQC\t([0-9.]+)$", text, re.MULTILINE)
    if version_match is None:
        raise ValueError(f"FastQC version missing in {path}")
    modules: list[dict[str, str]] = []
    basic_statistics: dict[str, str] = {}
    current_module = ""
    for line in text.splitlines():
        if line.startswith(">>") and line != ">>END_MODULE":
            fields = line[2:].split("\t")
            if len(fields) == 2:
                current_module = fields[0]
                modules.append(
                    {"Module": fields[0], "Status": fields[1].lower()}
                )
            continue
        if (
            current_module == "Basic Statistics"
            and line
            and not line.startswith("#")
            and "\t" in line
        ):
            key, value = line.split("\t", 1)
            basic_statistics[key] = value
    return {
        "path": path,
        "version": version_match.group(1),
        "modules": modules,
        "basic_statistics": basic_statistics,
    }


def discover_fastqc_reports(frozen_dir: Path) -> list[dict[str, Any]]:
    paths = sorted(
        list((frozen_dir / "raw_fastqc").glob("*_fastqc.zip"))
        + list((frozen_dir / "clean_fastqc").glob("*_fastqc.zip"))
    )
    reports: list[dict[str, Any]] = []
    for path in paths:
        parsed = parse_fastqc_zip(path)
        name = path.name
        stage = "Raw" if "prefix100k" in name else "Clean"
        mate_match = re.search(r"_(R[12])_fastqc\.zip$", name)
        if mate_match is None:
            raise ValueError(f"Could not infer mate from {name}")
        parsed["stage"] = stage
        parsed["mate"] = mate_match.group(1)
        parsed["report"] = f"{stage} {mate_match.group(1)}"
        reports.append(parsed)
    return reports


def multiqc_payload(frozen_dir: Path) -> tuple[Path, dict[str, Any]]:
    matches = list(
        (frozen_dir / "multiqc").glob("*_data/multiqc_data.json")
    )
    if len(matches) != 1:
        raise ValueError("Expected one MultiQC multiqc_data.json")
    path = matches[0]
    return path, json.loads(path.read_text(encoding="utf-8"))


def parse_explicit_lock(path: Path) -> tuple[int, dict[str, str]]:
    package_urls = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith("https://")
    ]
    observed: dict[str, str] = {}
    targets = {
        "fastqc": r"^fastqc-",
        "fastp": r"^fastp-",
        "multiqc": r"^multiqc-",
        "python": r"^python-3\.14\.6-",
        "openjdk": r"^openjdk-",
        "matplotlib-base": r"^matplotlib-base-",
    }
    for package, pattern in targets.items():
        candidates = [
            url.rsplit("/", 1)[-1]
            for url in package_urls
            if re.match(pattern, url.rsplit("/", 1)[-1])
        ]
        observed[package] = candidates[0] if len(candidates) == 1 else ""
    return len(package_urls), observed


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


def verify_checksum_manifest(
    frozen_dir: Path,
) -> tuple[int, int, list[dict[str, Any]]]:
    checksum_path = frozen_dir / "file-checksums.sha256"
    expected: dict[str, str] = {}
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        checksum, relative = line.split("  ", 1)
        expected[relative] = checksum
    actual_files = {
        path.relative_to(frozen_dir).as_posix()
        for path in frozen_dir.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    }
    rows: list[dict[str, Any]] = []
    failures = 0
    for relative in sorted(set(expected) | actual_files):
        path = frozen_dir / relative
        observed = file_sha256(path) if path.is_file() else ""
        expected_hash = expected.get(relative, "")
        status = (
            "PASS"
            if relative in actual_files
            and relative in expected
            and observed == expected_hash
            else "FAIL"
        )
        failures += status == "FAIL"
        rows.append(
            {
                "File": relative,
                "ExpectedSHA256": expected_hash,
                "ObservedSHA256": observed,
                "Status": status,
            }
        )
    return len(rows), failures, rows


def write_checksum_manifest(frozen_dir: Path) -> None:
    checksum_path = frozen_dir / "file-checksums.sha256"
    files = sorted(
        path
        for path in frozen_dir.rglob("*")
        if path.is_file() and path.name != checksum_path.name
    )
    lines = [
        f"{file_sha256(path)}  {path.relative_to(frozen_dir).as_posix()}"
        for path in files
    ]
    checksum_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def frozen_text_files(frozen_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in frozen_dir.rglob("*")
        if path.is_file()
        and path.name != "file-checksums.sha256"
        and path.suffix.lower() in FROZEN_TEXT_SUFFIXES
    )


def normalize_frozen_paths(
    project_root: Path,
    environment_prefix: Path,
    frozen_dir: Path,
    raw_dir: Path,
) -> dict[str, Any]:
    replacements = (
        (str(raw_dir), "${RAW_DIR}"),
        (str(frozen_dir), "${FROZEN_DIR}"),
        (str(environment_prefix), "${READ_QC_ENV_PREFIX}"),
        (str(project_root), "${PROJECT_ROOT}"),
    )
    for path in frozen_text_files(frozen_dir):
        text = path.read_text(encoding="utf-8")
        normalized = text
        for source, placeholder in replacements:
            normalized = normalized.replace(source, placeholder)
        normalized = re.sub(
            r"/tmp/tmp[A-Za-z0-9._-]+",
            "${TMPDIR}",
            normalized,
        )
        if normalized != text:
            path.write_text(normalized, encoding="utf-8")
    audit = audit_frozen_paths(
        project_root,
        environment_prefix,
        frozen_dir,
        raw_dir,
    )
    return {
        "status": "passed" if audit["remaining_hits"] == 0 else "failed",
        "placeholders": [
            "${RAW_DIR}",
            "${FROZEN_DIR}",
            "${READ_QC_ENV_PREFIX}",
            "${PROJECT_ROOT}",
            "${TMPDIR}",
        ],
        "host_specific_paths_retained": audit["remaining_hits"] != 0,
        "remaining_hits": audit["remaining_hits"],
    }


def audit_frozen_paths(
    project_root: Path,
    environment_prefix: Path,
    frozen_dir: Path,
    raw_dir: Path | None = None,
) -> dict[str, Any]:
    forbidden_literals = {
        str(project_root),
        str(environment_prefix),
    }
    if raw_dir is not None:
        forbidden_literals.add(str(raw_dir))
    hits: list[str] = []
    for path in frozen_text_files(frozen_dir):
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(frozen_dir).as_posix()
        for literal in sorted(forbidden_literals):
            if literal in text:
                hits.append(f"{relative}:literal")
        if re.search(r"/tmp/tmp[A-Za-z0-9._-]+", text):
            hits.append(f"{relative}:temporary")
    return {
        "remaining_hits": len(hits),
        "files": sorted(set(hits)),
    }


def build_evidence(
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

    lock_path = project_root / "env/read-qc-linux-64.lock"
    lock_count, lock_packages = parse_explicit_lock(lock_path)
    checks.add(
        "environment-lock-package-count",
        lock_count == 171,
        f"linux-64 explicit package rows={lock_count}",
    )
    lock_expectations = {
        "fastqc": "fastqc-0.12.1-",
        "fastp": "fastp-1.3.6-",
        "multiqc": "multiqc-1.35-",
        "python": "python-3.14.6-",
        "openjdk": "openjdk-25.0.2-",
        "matplotlib-base": "matplotlib-base-3.10.8-",
    }
    for package, prefix in lock_expectations.items():
        checks.add(
            f"lock-{package}",
            lock_packages[package].startswith(prefix),
            lock_packages[package],
        )

    source_rows = read_tsv(project_root / "data/small/13-source-manifest.tsv")
    checks.add(
        "source-row-count",
        len(source_rows) == 2,
        f"rows={len(source_rows)}",
    )
    checks.add(
        "source-mates",
        [row["Mate"] for row in source_rows] == ["R1", "R2"],
        ",".join(row["Mate"] for row in source_rows),
    )
    for key, expected in EXPECTED_SOURCE.items():
        observed = {row[key] for row in source_rows}
        checks.add(
            f"source-{key.lower()}",
            observed == {expected},
            f"expected={expected}; observed={sorted(observed)}",
        )
    checks.add(
        "source-archive-md5",
        {row["ENAReportedMD5"] for row in source_rows}
        == {
            "ed0e6e0ee846542531c742a45181cd6f",
            "5b60ac93cb69dff77ae38cfa501afd06",
        },
        "two ENA-reported complete-file MD5 values",
    )

    subset_summary = json.loads(
        (frozen_dir / "subset-summary.json").read_text(encoding="utf-8")
    )
    checks.add(
        "subset-status",
        subset_summary["status"] == "passed",
        subset_summary["status"],
    )
    checks.add(
        "subset-input-pairs",
        subset_summary["input_pairs"] == 100000,
        f"pairs={subset_summary['input_pairs']}",
    )
    checks.add(
        "subset-mate-sync",
        subset_summary["mates_synchronized"] is True,
        subset_summary["normalized_pair_id_sha256"],
    )
    checks.add(
        "subset-raw-not-committed",
        subset_summary["raw_fastq_committed"] is False,
        "raw_fastq_committed=false",
    )
    checks.add(
        "subset-full-md5-boundary",
        subset_summary["complete_archive_md5_verified"] is False,
        "streaming prefix does not recompute complete-file MD5",
    )

    fastp_path = frozen_dir / "fastp/ERR9765746_fastp.json"
    fastp_data = json.loads(fastp_path.read_text(encoding="utf-8"))
    before = fastp_data["summary"]["before_filtering"]
    after = fastp_data["summary"]["after_filtering"]
    filtering = fastp_data["filtering_result"]
    adapter = fastp_data["adapter_cutting"]
    input_reads = int(before["total_reads"])
    filter_sum = sum(int(filtering[key]) for key in FASTP_FILTER_KEYS)
    read_ledger_delta = input_reads - filter_sum
    checks.add(
        "fastp-version-json",
        fastp_data["summary"]["fastp_version"] == "1.3.6",
        fastp_data["summary"]["fastp_version"],
    )
    checks.add(
        "fastp-input-reads",
        input_reads == 200000,
        f"reads={input_reads}",
    )
    checks.add(
        "fastp-read-ledger",
        read_ledger_delta == 0,
        f"delta={read_ledger_delta}",
    )
    checks.add(
        "fastp-passed-paired",
        int(filtering["passed_filter_reads"]) % 2 == 0,
        f"passed reads={filtering['passed_filter_reads']}",
    )
    command = fastp_data["command"]
    for flag in EXPECTED_FASTP_FLAGS:
        checks.add(
            f"fastp-flag-{flag.split()[0].lstrip('-').replace('_', '-')}",
            flag in command,
            flag,
        )
    for flag in PROHIBITED_FASTP_FLAGS:
        checks.add(
            f"fastp-absent-{flag.lstrip('-').replace('_', '-')}",
            flag not in command,
            f"{flag} absent",
        )

    fastqc_reports = discover_fastqc_reports(frozen_dir)
    checks.add(
        "fastqc-report-count",
        len(fastqc_reports) == 4,
        f"reports={len(fastqc_reports)}",
    )
    module_counts = {len(report["modules"]) for report in fastqc_reports}
    checks.add(
        "fastqc-module-count",
        module_counts == {11},
        f"module counts={sorted(module_counts)}",
    )
    for report in fastqc_reports:
        expected_sequences = (
            100000 if report["stage"] == "Raw" else int(after["total_reads"]) // 2
        )
        observed_sequences = int(
            report["basic_statistics"]["Total Sequences"]
        )
        checks.add(
            f"fastqc-{report['stage'].lower()}-{report['mate'].lower()}-sequences",
            observed_sequences == expected_sequences,
            f"expected={expected_sequences}; observed={observed_sequences}",
        )
        checks.add(
            f"fastqc-{report['stage'].lower()}-{report['mate'].lower()}-version",
            report["version"] == "0.12.1",
            report["version"],
        )

    multiqc_path, multiqc_data = multiqc_payload(frozen_dir)
    data_sources = multiqc_data["report_data_sources"]
    fastp_sources = data_sources.get("fastp", {}).get("all_sections", {})
    fastqc_sources = data_sources.get("FastQC", {}).get("all_sections", {})
    checks.add(
        "multiqc-fastp-report-count",
        len(fastp_sources) == 1,
        f"fastp reports={len(fastp_sources)}",
    )
    checks.add(
        "multiqc-fastqc-report-count",
        len(fastqc_sources) == 4,
        f"FastQC reports={len(fastqc_sources)}",
    )
    checks.add(
        "multiqc-unique-sample-names",
        len(set(fastqc_sources)) == 4,
        ",".join(sorted(fastqc_sources)),
    )
    checks.add(
        "multiqc-version",
        str(multiqc_data["config_version"]) == "1.35",
        str(multiqc_data["config_version"]),
    )

    frozen_fastq = [
        path
        for path in frozen_dir.rglob("*")
        if path.is_file()
        and (
            path.name.endswith(".fastq")
            or path.name.endswith(".fastq.gz")
            or path.name.endswith(".fq")
            or path.name.endswith(".fq.gz")
        )
    ]
    checks.add(
        "frozen-no-fastq",
        not frozen_fastq,
        f"FASTQ files in frozen evidence={len(frozen_fastq)}",
    )

    fastqc_rows: list[dict[str, Any]] = []
    status_counter: Counter[str] = Counter()
    for report in fastqc_reports:
        for module in report["modules"]:
            status_counter[module["Status"]] += 1
            fastqc_rows.append(
                {
                    "Report": report["report"],
                    "Stage": report["stage"],
                    "Mate": report["mate"],
                    "Module": module["Module"],
                    "Status": module["Status"].upper(),
                }
            )

    source_audit_rows: list[dict[str, Any]] = []
    subset_by_mate = {
        row["mate"]: row for row in subset_summary["mates"]
    }
    for row in source_rows:
        mate_metrics = subset_by_mate[row["Mate"]]
        source_audit_rows.append(
            {
                "Mate": row["Mate"],
                "RunAccession": row["RunAccession"],
                "PrefixRecords": row["PrefixRecords"],
                "ENAReportedMD5": row["ENAReportedMD5"],
                "ENABytes": row["ENABytes"],
                "HTTPContentLengthMatches": str(
                    mate_metrics["http_content_length_matches_archive"]
                ).lower(),
                "PrefixUncompressedSHA256": mate_metrics[
                    "uncompressed_fastq_sha256"
                ],
                "PrefixCompressedSHA256": mate_metrics[
                    "compressed_fastq_sha256"
                ],
                "FullMD5Verified": "false",
            }
        )

    fastp_audit_rows = [
        {
            "Metric": "Input read pairs",
            "Value": input_reads // 2,
            "Unit": "pairs",
            "Interpretation": "Deterministic file prefix",
        },
        {
            "Metric": "Passed read pairs",
            "Value": int(filtering["passed_filter_reads"]) // 2,
            "Unit": "pairs",
            "Interpretation": "Both mates retained",
        },
        {
            "Metric": "Too-short read pairs",
            "Value": int(filtering["too_short_reads"]) // 2,
            "Unit": "pairs",
            "Interpretation": "At least one mate failed length >=50 bp",
        },
        {
            "Metric": "Adapter-trimmed reads",
            "Value": int(adapter["adapter_trimmed_reads"]),
            "Unit": "reads",
            "Interpretation": "Zero is a valid evidence-based outcome",
        },
        {
            "Metric": "Duplication estimate",
            "Value": 100 * float(fastp_data["duplication"]["rate"]),
            "Unit": "percent",
            "Interpretation": "Reported, not deduplicated",
        },
        {
            "Metric": "Insert-size peak",
            "Value": int(fastp_data["insert_size"]["peak"]),
            "Unit": "bp",
            "Interpretation": "Overlap-derived estimate",
        },
        {
            "Metric": "Q30 before",
            "Value": 100 * float(before["q30_rate"]),
            "Unit": "percent bases",
            "Interpretation": "Combined R1 and R2",
        },
        {
            "Metric": "Q30 after",
            "Value": 100 * float(after["q30_rate"]),
            "Unit": "percent bases",
            "Interpretation": "Combined R1 and R2",
        },
    ]
    multiqc_rows = [
        {
            "Module": "FastQC",
            "ReportsFound": len(fastqc_sources),
            "ExpectedReports": 4,
            "Status": "PASS" if len(fastqc_sources) == 4 else "FAIL",
            "SampleNames": " | ".join(sorted(fastqc_sources)),
        },
        {
            "Module": "fastp",
            "ReportsFound": len(fastp_sources),
            "ExpectedReports": 1,
            "Status": "PASS" if len(fastp_sources) == 1 else "FAIL",
            "SampleNames": " | ".join(sorted(fastp_sources)),
        },
    ]
    return {
        "versions": versions,
        "tool_rows": tool_rows,
        "lock_count": lock_count,
        "lock_packages": lock_packages,
        "source_rows": source_rows,
        "source_audit_rows": source_audit_rows,
        "subset_summary": subset_summary,
        "fastp_data": fastp_data,
        "fastp_audit_rows": fastp_audit_rows,
        "fastqc_reports": fastqc_reports,
        "fastqc_rows": fastqc_rows,
        "fastqc_status_counter": status_counter,
        "multiqc_path": multiqc_path,
        "multiqc_data": multiqc_data,
        "multiqc_rows": multiqc_rows,
        "input_reads": input_reads,
        "before": before,
        "after": after,
        "filtering": filtering,
        "adapter": adapter,
        "read_ledger_delta": read_ledger_delta,
    }


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


def plot_quality(evidence: dict[str, Any], figure_dir: Path) -> None:
    fastp_data = evidence["fastp_data"]
    stages = (
        ("Before filtering", "before_filtering", "#7A7A7A"),
        ("After filtering", "after_filtering", "#0072B2"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.1), sharey=True)
    for axis, mate_number in zip(axes, (1, 2)):
        for label, key, color in stages:
            values = fastp_data[f"read{mate_number}_{key}"][
                "quality_curves"
            ]["mean"]
            cycles = np.arange(1, len(values) + 1)
            axis.plot(cycles, values, color=color, linewidth=2, label=label)
        axis.axhline(30, color="#D55E00", linestyle="--", linewidth=1)
        axis.text(
            149,
            30.3,
            "Q30",
            color="#D55E00",
            ha="right",
            va="bottom",
            fontsize=8,
        )
        axis.set_title(f"Read {mate_number}")
        axis.set_xlabel("Cycle")
        axis.set_xlim(1, 150)
        axis.set_ylim(28, 43)
        axis.grid(axis="y", color="#E6E6E6", linewidth=0.7)
    axes[0].set_ylabel("Mean Phred score")
    axes[1].legend(frameon=False, loc="lower left")
    fig.suptitle(
        "Per-cycle quality changed negligibly after conservative filtering",
        y=1.01,
        fontsize=12,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.01,
        "ERR9765746 deterministic first 100,000 read pairs",
        ha="center",
        color="#555555",
        fontsize=8.5,
    )
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEM_QUALITY)


def plot_read_fate(evidence: dict[str, Any], figure_dir: Path) -> None:
    filtering = evidence["filtering"]
    input_pairs = evidence["input_reads"] // 2
    passed_pairs = int(filtering["passed_filter_reads"]) // 2
    removed_pairs = input_pairs - passed_pairs
    retained_percent = 100 * passed_pairs / input_pairs
    removed_percent = 100 - retained_percent
    reasons = [
        ("Too short", int(filtering["too_short_reads"]) // 2),
        ("Low quality", int(filtering["low_quality_reads"]) // 2),
        ("Too many N", int(filtering["too_many_N_reads"]) // 2),
        ("Adapter dimer", int(filtering["adapter_dimer_reads"]) // 2),
        ("Too long", int(filtering["too_long_reads"]) // 2),
    ]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9.4, 4.0),
        gridspec_kw={"width_ratios": [1.5, 1]},
    )
    axes[0].barh(
        ["All input pairs"],
        [retained_percent],
        color="#009E73",
        height=0.42,
        label="Retained",
    )
    axes[0].barh(
        ["All input pairs"],
        [removed_percent],
        left=[retained_percent],
        color="#D55E00",
        height=0.42,
        label="Removed",
    )
    axes[0].set_xlim(0, 100)
    axes[0].set_xlabel("Read pairs (%)")
    axes[0].text(
        50,
        0,
        f"{passed_pairs:,} retained ({retained_percent:.3f}%)",
        ha="center",
        va="center",
        color="white",
        fontweight="bold",
    )
    axes[0].annotate(
        f"{removed_pairs} removed\n({removed_percent:.3f}%)",
        xy=(retained_percent, 0),
        xytext=(76, 0.48),
        textcoords="data",
        arrowprops={"arrowstyle": "->", "color": "#333333"},
        ha="center",
        va="bottom",
        fontsize=8.5,
    )
    axes[0].legend(frameon=False, loc="lower left", ncol=2)
    axes[0].grid(axis="x", color="#E6E6E6", linewidth=0.7)
    axes[0].set_title("Pair retention")

    labels = [row[0] for row in reasons]
    values = [row[1] for row in reasons]
    bars = axes[1].bar(
        labels,
        values,
        color=["#D55E00", "#999999", "#999999", "#999999", "#999999"],
    )
    axes[1].set_ylabel("Removed read pairs")
    axes[1].set_title("Filter-specific losses")
    axes[1].tick_params(axis="x", rotation=35)
    axes[1].set_ylim(0, max(values + [1]) + 2)
    axes[1].grid(axis="y", color="#E6E6E6", linewidth=0.7)
    for bar, value in zip(bars, values):
        axes[1].text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.18,
            str(value),
            ha="center",
            va="bottom",
            fontsize=8.5,
        )
    fig.suptitle(
        "Length filtering removed nine pairs; adapter trimming removed none",
        y=1.02,
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEM_FATE)


def plot_fastqc_modules(evidence: dict[str, Any], figure_dir: Path) -> None:
    reports = evidence["fastqc_reports"]
    module_order = [row["Module"] for row in reports[0]["modules"]]
    status_value = {"pass": 0, "warn": 1, "fail": 2}
    matrix = np.array(
        [
            [
                status_value[
                    next(
                        row["Status"]
                        for row in report["modules"]
                        if row["Module"] == module
                    )
                ]
                for module in module_order
            ]
            for report in reports
        ],
        dtype=int,
    )
    colors = matplotlib.colors.ListedColormap(
        ["#009E73", "#E69F00", "#D55E00"]
    )
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = matplotlib.colors.BoundaryNorm(bounds, colors.N)
    fig, axis = plt.subplots(figsize=(10.2, 4.4))
    axis.imshow(matrix, cmap=colors, norm=norm, aspect="auto")
    axis.set_xticks(np.arange(len(module_order)))
    axis.set_xticklabels(module_order, rotation=42, ha="right")
    axis.set_yticks(np.arange(len(reports)))
    axis.set_yticklabels([report["report"] for report in reports])
    axis.set_xlabel("FastQC module")
    axis.set_title(
        "FastQC states are diagnostic signals, not automatic discard rules",
        pad=12,
        fontsize=12,
        fontweight="bold",
    )
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            text = ("PASS", "WARN", "FAIL")[matrix[row_index, column_index]]
            axis.text(
                column_index,
                row_index,
                text,
                ha="center",
                va="center",
                color="white" if text != "WARN" else "#222222",
                fontsize=6.8,
                fontweight="bold",
                rotation=90,
            )
    legend_handles = [
        matplotlib.patches.Patch(color="#009E73", label="PASS"),
        matplotlib.patches.Patch(color="#E69F00", label="WARN"),
        matplotlib.patches.Patch(color="#D55E00", label="FAIL"),
    ]
    axis.legend(
        handles=legend_handles,
        frameon=False,
        ncol=3,
        bbox_to_anchor=(0.5, -0.43),
        loc="upper center",
    )
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEM_MODULES)


def initialize_frozen(
    project_root: Path,
    environment_prefix: Path,
    frozen_dir: Path,
    raw_dir: Path,
) -> int:
    checks = Checks()
    path_normalization = normalize_frozen_paths(
        project_root,
        environment_prefix,
        frozen_dir,
        raw_dir,
    )
    checks.add(
        "frozen-path-normalization",
        path_normalization["status"] == "passed",
        (
            "host-specific or transient path hits="
            f"{path_normalization['remaining_hits']}"
        ),
    )
    evidence = build_evidence(
        project_root,
        environment_prefix,
        frozen_dir,
        checks,
    )
    raw_audit = audit_fastq_pair(
        raw_dir / "ERR9765746_prefix100k_R1.fastq.gz",
        raw_dir / "ERR9765746_prefix100k_R2.fastq.gz",
    )
    clean_audit = audit_fastq_pair(
        raw_dir / "ERR9765746_clean_R1.fastq.gz",
        raw_dir / "ERR9765746_clean_R2.fastq.gz",
    )
    subset_summary = evidence["subset_summary"]
    subset_by_mate = {
        row["mate"]: row for row in subset_summary["mates"]
    }
    checks.add(
        "raw-pair-count",
        raw_audit["pairs"] == 100000,
        f"pairs={raw_audit['pairs']}",
    )
    checks.add(
        "raw-pair-id-hash",
        raw_audit["normalized_pair_id_sha256"]
        == subset_summary["normalized_pair_id_sha256"],
        raw_audit["normalized_pair_id_sha256"],
    )
    for mate in ("R1", "R2"):
        checks.add(
            f"raw-{mate.lower()}-uncompressed-hash",
            raw_audit[mate]["uncompressed_fastq_sha256"]
            == subset_by_mate[mate]["uncompressed_fastq_sha256"],
            raw_audit[mate]["uncompressed_fastq_sha256"],
        )
    retained_pairs = int(
        evidence["filtering"]["passed_filter_reads"]
    ) // 2
    checks.add(
        "clean-pair-count",
        clean_audit["pairs"] == retained_pairs,
        f"expected={retained_pairs}; observed={clean_audit['pairs']}",
    )
    checks.add(
        "clean-mates-synchronized",
        clean_audit["mates_synchronized"] is True,
        clean_audit["normalized_pair_id_sha256"],
    )

    resource_summary = {
        stem: parse_resource_log(frozen_dir / f"logs/{stem}.resources.txt")
        for stem in (
            "fastqc-before",
            "fastp",
            "fastqc-after",
            "multiqc",
        )
    }
    filtering = evidence["filtering"]
    status_counter = evidence["fastqc_status_counter"]
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "evidence_date": "2026-07-20",
        "source_project": "PRJEB52977",
        "source_sample": "SAMEA14435832",
        "run_accession": "ERR9765746",
        "selection": "first 100000 complete synchronized read pairs",
        "input_pairs": evidence["input_reads"] // 2,
        "retained_pairs": int(filtering["passed_filter_reads"]) // 2,
        "removed_pairs": (
            evidence["input_reads"]
            - int(filtering["passed_filter_reads"])
        )
        // 2,
        "retained_pair_fraction": (
            int(filtering["passed_filter_reads"])
            / evidence["input_reads"]
        ),
        "low_quality_pairs": int(filtering["low_quality_reads"]) // 2,
        "too_many_n_pairs": int(filtering["too_many_N_reads"]) // 2,
        "adapter_dimer_pairs": int(filtering["adapter_dimer_reads"]) // 2,
        "too_short_pairs": int(filtering["too_short_reads"]) // 2,
        "too_long_pairs": int(filtering["too_long_reads"]) // 2,
        "adapter_trimmed_reads": int(
            evidence["adapter"]["adapter_trimmed_reads"]
        ),
        "adapter_trimmed_bases": int(
            evidence["adapter"]["adapter_trimmed_bases"]
        ),
        "duplication_rate": float(evidence["fastp_data"]["duplication"]["rate"]),
        "insert_size_peak_bp": int(
            evidence["fastp_data"]["insert_size"]["peak"]
        ),
        "q30_rate_before": float(evidence["before"]["q30_rate"]),
        "q30_rate_after": float(evidence["after"]["q30_rate"]),
        "fastp_read_ledger_delta": evidence["read_ledger_delta"],
        "fastqc_reports": len(evidence["fastqc_reports"]),
        "fastqc_module_states_per_report": len(
            evidence["fastqc_reports"][0]["modules"]
        ),
        "fastqc_pass_states": status_counter["pass"],
        "fastqc_warn_states": status_counter["warn"],
        "fastqc_fail_states": status_counter["fail"],
        "multiqc_fastqc_reports": int(
            evidence["multiqc_rows"][0]["ReportsFound"]
        ),
        "multiqc_fastp_reports": int(
            evidence["multiqc_rows"][1]["ReportsFound"]
        ),
        "fastqc_version": evidence["versions"]["FastQC"],
        "fastp_version": evidence["versions"]["fastp"],
        "multiqc_version": evidence["versions"]["MultiQC"],
        "python_version": evidence["versions"]["Python"],
        "environment_lock_packages": evidence["lock_count"],
        "environment_yaml_sha256": file_sha256(
            project_root / "env/read-qc.yml"
        ),
        "environment_lock_sha256": file_sha256(
            project_root / "env/read-qc-linux-64.lock"
        ),
        "source_manifest_sha256": file_sha256(
            project_root / "data/small/13-source-manifest.tsv"
        ),
        "raw_fastq_audit": raw_audit,
        "clean_fastq_audit": clean_audit,
        "raw_fastq_committed": False,
        "clean_fastq_committed": False,
        "complete_archive_md5_verified": False,
        "one_time_network_access": True,
        "qa_network_access": False,
        "resource_summary": resource_summary,
        "path_normalization": path_normalization,
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
    path_audit = audit_frozen_paths(
        project_root,
        environment_prefix,
        frozen_dir,
    )
    checks.add(
        "frozen-path-normalization",
        path_audit["remaining_hits"] == 0,
        f"host-specific or transient path hits={path_audit['remaining_hits']}",
    )
    evidence = build_evidence(
        project_root,
        environment_prefix,
        frozen_dir,
        checks,
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
        "frozen-input-pairs",
        frozen_summary["input_pairs"] == 100000,
        str(frozen_summary["input_pairs"]),
    )
    checks.add(
        "frozen-retained-pairs",
        frozen_summary["retained_pairs"]
        == int(evidence["filtering"]["passed_filter_reads"]) // 2,
        str(frozen_summary["retained_pairs"]),
    )
    checks.add(
        "frozen-clean-sync",
        frozen_summary["clean_fastq_audit"]["mates_synchronized"] is True,
        frozen_summary["clean_fastq_audit"]["normalized_pair_id_sha256"],
    )
    checks.add(
        "routine-no-network",
        frozen_summary["qa_network_access"] is False,
        "qa_network_access=false",
    )

    configure_plotting()
    plot_quality(evidence, figure_dir)
    plot_read_fate(evidence, figure_dir)
    plot_fastqc_modules(evidence, figure_dir)
    for stem in (
        FIGURE_STEM_QUALITY,
        FIGURE_STEM_FATE,
        FIGURE_STEM_MODULES,
    ):
        for suffix in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            checks.add(
                f"figure-{stem}-{suffix}",
                path.is_file() and path.stat().st_size > 0,
                f"{path.name} bytes={path.stat().st_size if path.exists() else 0}",
            )

    write_tsv(
        output_dir / "tool-audit.tsv",
        evidence["tool_rows"],
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
        evidence["source_audit_rows"],
        [
            "Mate",
            "RunAccession",
            "PrefixRecords",
            "ENAReportedMD5",
            "ENABytes",
            "HTTPContentLengthMatches",
            "PrefixUncompressedSHA256",
            "PrefixCompressedSHA256",
            "FullMD5Verified",
        ],
    )
    write_tsv(
        output_dir / "fastp-audit.tsv",
        evidence["fastp_audit_rows"],
        ["Metric", "Value", "Unit", "Interpretation"],
    )
    write_tsv(
        output_dir / "fastqc-module-audit.tsv",
        evidence["fastqc_rows"],
        ["Report", "Stage", "Mate", "Module", "Status"],
    )
    write_tsv(
        output_dir / "multiqc-audit.tsv",
        evidence["multiqc_rows"],
        [
            "Module",
            "ReportsFound",
            "ExpectedReports",
            "Status",
            "SampleNames",
        ],
    )
    write_tsv(
        output_dir / "checksum-audit.tsv",
        checksum_rows,
        ["File", "ExpectedSHA256", "ObservedSHA256", "Status"],
    )

    filtering = evidence["filtering"]
    status_counter = evidence["fastqc_status_counter"]
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "run_accession": "ERR9765746",
        "source_project": "PRJEB52977",
        "source_sample": "SAMEA14435832",
        "input_pairs": evidence["input_reads"] // 2,
        "retained_pairs": int(filtering["passed_filter_reads"]) // 2,
        "removed_pairs": (
            evidence["input_reads"]
            - int(filtering["passed_filter_reads"])
        )
        // 2,
        "retained_pair_fraction": (
            int(filtering["passed_filter_reads"])
            / evidence["input_reads"]
        ),
        "too_short_pairs": int(filtering["too_short_reads"]) // 2,
        "adapter_trimmed_reads": int(
            evidence["adapter"]["adapter_trimmed_reads"]
        ),
        "duplication_rate": float(evidence["fastp_data"]["duplication"]["rate"]),
        "insert_size_peak_bp": int(
            evidence["fastp_data"]["insert_size"]["peak"]
        ),
        "q30_rate_before": float(evidence["before"]["q30_rate"]),
        "q30_rate_after": float(evidence["after"]["q30_rate"]),
        "fastp_read_ledger_delta": evidence["read_ledger_delta"],
        "fastqc_reports": len(evidence["fastqc_reports"]),
        "fastqc_module_states_per_report": len(
            evidence["fastqc_reports"][0]["modules"]
        ),
        "fastqc_pass_states": status_counter["pass"],
        "fastqc_warn_states": status_counter["warn"],
        "fastqc_fail_states": status_counter["fail"],
        "multiqc_fastqc_reports": int(
            evidence["multiqc_rows"][0]["ReportsFound"]
        ),
        "multiqc_fastp_reports": int(
            evidence["multiqc_rows"][1]["ReportsFound"]
        ),
        "fastqc_version": evidence["versions"]["FastQC"],
        "fastp_version": evidence["versions"]["fastp"],
        "multiqc_version": evidence["versions"]["MultiQC"],
        "python_version": evidence["versions"]["Python"],
        "environment_lock_packages": evidence["lock_count"],
        "raw_fastq_committed": False,
        "clean_fastq_committed": False,
        "complete_archive_md5_verified": False,
        "qa_network_access": False,
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
        else project_root / "data/small/13-qc-frozen"
    )
    if args.initialize_frozen:
        if args.raw_dir is None:
            raise SystemExit("--raw-dir is required with --initialize-frozen")
        return initialize_frozen(
            project_root,
            environment_prefix,
            frozen_dir,
            args.raw_dir.resolve(),
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
