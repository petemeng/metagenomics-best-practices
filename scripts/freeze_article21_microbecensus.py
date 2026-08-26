#!/usr/bin/env python3
"""Validate and freeze publishable MicrobeCensus evidence for Article 21."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


EXPECTED_COMMIT = "dfc42d356bfd7943633cde6c0fbfc0b116f29ae2"
EXPECTED_R1_SHA256 = "ce438de2814a6e605cc24b00e53b697d018f325bc2a9694cece5be158ff32101"
EXPECTED_R2_SHA256 = "207d080cd5dc98bd10fd4af8c492fb3f2000f73e1458a9d1d77e57347f4fe459"
EXPECTED_READS = 199_982
EXPECTED_PAIRS = 99_991
EXPECTED_BASES = 29_809_773


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(root: Path, *args: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(root), *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def count_fastq(path: Path) -> tuple[int, int]:
    reads = 0
    bases = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Invalid FASTQ structure in {path}")
            if len(sequence) != len(quality):
                raise ValueError(f"Sequence/quality length mismatch in {path}")
            reads += 1
            bases += len(sequence)
    return reads, bases


def parse_output(path: Path) -> dict[str, str]:
    parsed: dict[str, str] = {}
    section = ""
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line in {"Parameters", "Results"}:
            section = line.lower()
            continue
        key, value = line.split(":\t", 1)
        parsed[f"{section}.{key}"] = value
    return parsed


def parse_resource(path: Path) -> dict[str, str]:
    keep = {
        "User time (seconds)": "user_seconds",
        "System time (seconds)": "system_seconds",
        "Elapsed (wall clock) time (h:mm:ss or m:ss)": "elapsed",
        "Maximum resident set size (kbytes)": "maximum_rss_kib",
        "Exit status": "exit_status",
    }
    parsed: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        for source, target in keep.items():
            prefix = source + ":"
            if line.startswith(prefix):
                parsed[target] = line[len(prefix) :].strip()
    return parsed


def normalized_log(text: str, replacements: dict[str, str]) -> str:
    for old, new in sorted(replacements.items(), key=lambda item: -len(item[0])):
        text = text.replace(old, new)
    return text


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    raw = args.raw_dir.resolve()
    source = args.source_root.resolve()
    frozen = args.frozen_dir.resolve()
    frozen.mkdir(parents=True, exist_ok=True)

    article13 = json.loads(
        (root / "data/small/13-qc-frozen/run-summary.json").read_text(encoding="utf-8")
    )
    r1 = root / "data/raw/article13/ERR9765746_clean_R1.fastq.gz"
    r2 = root / "data/raw/article13/ERR9765746_clean_R2.fastq.gz"
    if sha256(r1) != EXPECTED_R1_SHA256 or sha256(r2) != EXPECTED_R2_SHA256:
        raise SystemExit("Article 21 FASTQ checksum mismatch")
    r1_reads, r1_bases = count_fastq(r1)
    r2_reads, r2_bases = count_fastq(r2)
    if (r1_reads, r2_reads) != (EXPECTED_PAIRS, EXPECTED_PAIRS):
        raise SystemExit("Article 21 FASTQ record count mismatch")
    if r1_bases + r2_bases != EXPECTED_BASES:
        raise SystemExit("Article 21 FASTQ base count mismatch")
    if article13["clean_fastq_audit"]["R1"]["compressed_fastq_sha256"] != EXPECTED_R1_SHA256:
        raise SystemExit("Article 13 R1 lineage mismatch")
    if article13["clean_fastq_audit"]["R2"]["compressed_fastq_sha256"] != EXPECTED_R2_SHA256:
        raise SystemExit("Article 13 R2 lineage mismatch")

    commit = git_value(source, "rev-parse", "HEAD")
    tag = git_value(source, "describe", "--tags", "--exact-match")
    if commit != EXPECTED_COMMIT or tag != "v1.1.1":
        raise SystemExit("MicrobeCensus source identity mismatch")

    branch_rows: list[dict[str, object]] = []
    resource_rows: list[dict[str, object]] = []
    replacements = {
        str(root): "${PROJECT_ROOT}",
        str(raw): "${ARTICLE21_RAW_DIR}",
        str(source): "${MICROBECENSUS_ROOT}",
        str(Path.home()): "${HOME}",
    }
    for read_length in (150, 100):
        stem = f"read-length-{read_length}"
        parsed = parse_output(raw / f"{stem}.output.tsv")
        metadata = json.loads((raw / f"{stem}.metadata.json").read_text(encoding="utf-8"))
        resource = parse_resource(raw / f"{stem}.resource.txt")
        average_genome_size = float(parsed["results.average_genome_size"])
        total_bases = int(parsed["results.total_bases"])
        genome_equivalents = float(parsed["results.genome_equivalents"])
        if total_bases != EXPECTED_BASES:
            raise SystemExit(f"Unexpected total bases for {stem}")
        if abs(genome_equivalents - total_bases / average_genome_size) > 1e-12:
            raise SystemExit(f"Genome-equivalent formula mismatch for {stem}")
        if metadata["source_commit"] != EXPECTED_COMMIT:
            raise SystemExit(f"Source commit mismatch for {stem}")
        if int(resource.get("exit_status", "-1")) != 0:
            raise SystemExit(f"Nonzero resource exit for {stem}")

        log_text = (raw / f"{stem}.log").read_text(encoding="utf-8")
        short_match = re.search(r"(\d+) reads shorter than", log_text)
        hits_match = re.search(r"(\d+) reads hit marker proteins", log_text)
        assigned_match = re.search(r"(\d+) reads assigned to a marker protein", log_text)
        if not (short_match and hits_match and assigned_match):
            raise SystemExit(f"Could not parse MicrobeCensus log for {stem}")
        branch_rows.append(
            {
                "Branch": "Primary" if read_length == 150 else "Sensitivity",
                "ReadLengthBp": read_length,
                "ReadsInSequenceUniverse": EXPECTED_READS,
                "ReadsSampledForAGS": int(parsed["parameters.reads_sampled"]),
                "ReadsTooShort": int(short_match.group(1)),
                "MarkerHits": int(hits_match.group(1)),
                "AssignedMarkerReads": int(assigned_match.group(1)),
                "AverageGenomeSizeBp": average_genome_size,
                "TotalBases": total_bases,
                "GenomeEquivalents": genome_equivalents,
                "NReadsCeiling": int(metadata["nreads_ceiling"]),
                "Threads": int(metadata["threads"]),
                "QualityFiltering": "No",
                "DuplicateFiltering": "No",
            }
        )
        resource_rows.append(
            {
                "Branch": "Primary" if read_length == 150 else "Sensitivity",
                "ReadLengthBp": read_length,
                "UserSeconds": float(resource["user_seconds"]),
                "SystemSeconds": float(resource["system_seconds"]),
                "Elapsed": resource["elapsed"],
                "MaximumRSSKiB": int(resource["maximum_rss_kib"]),
                "ExitStatus": int(resource["exit_status"]),
            }
        )
        (frozen / f"microbecensus-{read_length}.log").write_text(
            normalized_log(log_text, replacements), encoding="utf-8"
        )

    primary = next(row for row in branch_rows if row["Branch"] == "Primary")
    with (frozen / "microbecensus-primary.tsv").open("w", encoding="utf-8") as handle:
        handle.write("Metric\tValue\tUnit\n")
        for metric, value, unit in (
            ("source_reads", EXPECTED_READS, "reads"),
            ("source_pairs", EXPECTED_PAIRS, "paired_fragments"),
            ("reads_sampled_for_ags", primary["ReadsSampledForAGS"], "reads"),
            ("read_length", primary["ReadLengthBp"], "bp"),
            ("average_genome_size", primary["AverageGenomeSizeBp"], "bp_per_genome"),
            ("total_bases", primary["TotalBases"], "bp"),
            ("genome_equivalents", primary["GenomeEquivalents"], "genome_coverage"),
        ):
            handle.write(f"{metric}\t{value}\t{unit}\n")

    columns = list(branch_rows[0])
    with (frozen / "microbecensus-read-length-sensitivity.tsv").open(
        "w", encoding="utf-8"
    ) as handle:
        handle.write("\t".join(columns) + "\n")
        for row in branch_rows:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")

    columns = list(resource_rows[0])
    with (frozen / "resource-usage.tsv").open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in resource_rows:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")

    source_rows = [
        {
            "SourceID": "mock1-clean-r1",
            "Type": "public-fastq-derived-clean-read",
            "Identity": "PRJEB52977/SAMEA14435832/ERR9765746 first 100000 pairs after Article 13 fastp",
            "Version": "Article13-2026-07-20",
            "SHA256": EXPECTED_R1_SHA256,
            "License": "ENA source terms; derived FASTQ excluded from Git",
            "PublicReference": "https://www.ebi.ac.uk/ena/browser/view/ERR9765746",
        },
        {
            "SourceID": "mock1-clean-r2",
            "Type": "public-fastq-derived-clean-read",
            "Identity": "PRJEB52977/SAMEA14435832/ERR9765746 first 100000 pairs after Article 13 fastp",
            "Version": "Article13-2026-07-20",
            "SHA256": EXPECTED_R2_SHA256,
            "License": "ENA source terms; derived FASTQ excluded from Git",
            "PublicReference": "https://www.ebi.ac.uk/ena/browser/view/ERR9765746",
        },
        {
            "SourceID": "microbecensus-source",
            "Type": "official-software-source",
            "Identity": commit,
            "Version": tag,
            "SHA256": sha256(source / "microbe_census/microbe_census.py"),
            "License": "GPL-3.0",
            "PublicReference": "https://github.com/snayfach/MicrobeCensus/tree/v1.1.1",
        },
    ]
    columns = list(source_rows[0])
    with (frozen / "source-manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in source_rows:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")

    python_version = sys.version.split()[0]
    import Bio  # type: ignore
    import numpy  # type: ignore

    versions = [
        ("MicrobeCensus source tag", tag),
        ("MicrobeCensus source commit", commit),
        ("MicrobeCensus internal version", "1.1.0"),
        ("RAPsearch2", "2.15"),
        ("Python", python_version),
        ("Biopython", Bio.__version__),
        ("NumPy", numpy.__version__),
        ("Compatibility shim", "RAPsearch2 preflight bytes-to-text only"),
    ]
    with (frozen / "tool-versions.tsv").open("w", encoding="utf-8") as handle:
        handle.write("Component\tVersion\n")
        for component, version in versions:
            handle.write(f"{component}\t{version}\n")

    commands = """#!/usr/bin/env bash
set -euo pipefail
export MICROBECENSUS_ROOT=${MICROBECENSUS_ROOT}
export BIOBAKERY_ENV_PREFIX=${BIOBAKERY_ENV_PREFIX}
PROJECT_ROOT=${PROJECT_ROOT} \\
  ARTICLE21_INPUT_R1=${PROJECT_ROOT}/data/raw/article13/ERR9765746_clean_R1.fastq.gz \\
  ARTICLE21_INPUT_R2=${PROJECT_ROOT}/data/raw/article13/ERR9765746_clean_R2.fastq.gz \\
  bash ${PROJECT_ROOT}/scripts/run_article21_microbecensus.sh
"""
    (frozen / "commands.sh").write_text(commands, encoding="utf-8")

    summary = {
        "status": "passed",
        "evidence_date": "2026-07-22",
        "source_project": "PRJEB52977",
        "source_sample": "SAMEA14435832",
        "run_accession": "ERR9765746",
        "source_pairs": EXPECTED_PAIRS,
        "source_reads": EXPECTED_READS,
        "source_total_bases": EXPECTED_BASES,
        "source_r1_sha256": EXPECTED_R1_SHA256,
        "source_r2_sha256": EXPECTED_R2_SHA256,
        "microbecensus_source_tag": tag,
        "microbecensus_source_commit": commit,
        "microbecensus_internal_version": "1.1.0",
        "compatibility_shim_scope": "RAPsearch2 preflight bytes-to-text only",
        "nreads_ceiling": 100_000_000,
        "threads": 8,
        "quality_filtering": False,
        "duplicate_filtering": False,
        "branches": branch_rows,
        "primary_read_length_bp": int(primary["ReadLengthBp"]),
        "primary_average_genome_size_bp": float(primary["AverageGenomeSizeBp"]),
        "primary_total_bases": int(primary["TotalBases"]),
        "primary_genome_equivalents": float(primary["GenomeEquivalents"]),
        "raw_fastq_committed": False,
        "qa_reads_fastq": False,
        "qa_network_access": False,
    }
    (frozen / "run-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    checksum_path = frozen / "file-checksums.sha256"
    if checksum_path.exists():
        checksum_path.unlink()
    entries = []
    for path in sorted(item for item in frozen.iterdir() if item.is_file()):
        entries.append(f"{sha256(path)}  {path.name}")
    checksum_path.write_text("\n".join(entries) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
