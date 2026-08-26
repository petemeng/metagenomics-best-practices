#!/usr/bin/env python3
"""Freeze compact, checksum-locked evidence for Article 32."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, TextIO


BRANCHES = (
    "spades-short-only",
    "spades-illumina-ont",
    "spades-illumina-hifi",
    "flye-ont",
    "flye-ont-polypolish-default",
    "flye-ont-polypolish-careful",
    "flye-hifi",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    return parser.parse_args()


def hashes(path: Path, algorithms: tuple[str, ...] = ("sha256",)) -> dict[str, str]:
    digests = {name: hashlib.new(name) for name in algorithms}
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            for digest in digests.values():
                digest.update(block)
    return {name: digest.hexdigest() for name, digest in digests.items()}


def sha256(path: Path) -> str:
    return hashes(path)["sha256"]


def count_fastq(path: Path) -> dict[str, object]:
    records = bases = 0
    minimum: int | None = None
    maximum = 0
    pair_ids = hashlib.sha256()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence_raw = handle.readline()
            plus = handle.readline()
            quality_raw = handle.readline()
            if not sequence_raw or not plus or not quality_raw:
                raise ValueError(f"Truncated FASTQ: {path}")
            sequence = sequence_raw.rstrip("\r\n")
            quality = quality_raw.rstrip("\r\n")
            if (
                not header.startswith("@")
                or not plus.startswith("+")
                or len(sequence) != len(quality)
            ):
                raise ValueError(f"Invalid FASTQ: {path}")
            length = len(sequence)
            records += 1
            bases += length
            minimum = length if minimum is None else min(minimum, length)
            maximum = max(maximum, length)
            token = header[1:].split(None, 1)[0]
            if token.endswith("/1") or token.endswith("/2"):
                token = token[:-2]
            pair_ids.update(token.encode())
            pair_ids.update(b"\n")
    if records == 0:
        raise ValueError(f"Empty FASTQ: {path}")
    return {
        "Records": records,
        "Bases": bases,
        "MinimumLength": minimum,
        "MaximumLength": maximum,
        "MeanLength": f"{bases / records:.12f}",
        "PairIDHash": pair_ids.hexdigest(),
        "CompressedBytes": path.stat().st_size,
        "CompressedSHA256": sha256(path),
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def write_deterministic_gzip(records: Iterable[tuple[str, str]], target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    raw = target.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0)
    text = io.TextIOWrapper(zipped, encoding="utf-8", newline="\n")
    try:
        for name, sequence in records:
            text.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                text.write(sequence[start : start + 80] + "\n")
    finally:
        text.close()
        raw.close()


def normalize_copy(source: Path, target: Path, replacements: dict[str, str]) -> None:
    text = source.read_text(encoding="utf-8", errors="replace")
    for old, new in sorted(replacements.items(), key=lambda item: -len(item[0])):
        text = text.replace(old, new)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")


def source_path(root: Path, row: dict[str, str]) -> Path:
    if row["RunAccession"] == "ERR9765746":
        return root / "data/raw/article32/sources" / f"ERR9765746_{row['Mate']}.fastq.gz"
    return root / "data/raw/article32/sources" / f"{row['RunAccession']}.fastq.gz"


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    env_prefix = args.env_prefix.resolve()
    raw = args.raw_dir.resolve()
    work = args.work_dir.resolve()
    frozen = args.frozen_dir.resolve()
    if frozen.exists() and any(frozen.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty frozen directory: {frozen}")
    frozen.mkdir(parents=True, exist_ok=True)
    for directory in (
        "assemblies",
        "env",
        "logs",
        "metaquast/references",
        "references/genomes",
        "resources",
        "scripts",
        "sources",
    ):
        (frozen / directory).mkdir(parents=True, exist_ok=True)

    required_complete = (
        work / "summary/run-summary.json",
        work / "summary/branch-metrics.tsv",
        work / "summary/per-genome-metaquast.tsv",
        work / "tool-versions.tsv",
        work / "metaquast/.article32-complete",
    )
    for path in required_complete:
        if not path.is_file():
            raise FileNotFoundError(path)

    replacements = {
        str(work): "${ARTICLE32_WORK_DIR}",
        str(raw): "${ARTICLE32_RAW_DIR}",
        str(root): "${PROJECT_ROOT}",
        str(env_prefix): "${HYBRID_ENV_PREFIX}",
        str(Path.home()): "${HOME}",
    }

    direct = {
        root / "data/small/32-source-manifest.tsv": frozen / "source-manifest.tsv",
        root / "data/small/32-reference-manifest.tsv": frozen / "reference-manifest.tsv",
        root / "data/small/32-software-releases.tsv": frozen / "software-releases.tsv",
        root / "data/small/32-branch-contract.tsv": frozen / "branch-contract.tsv",
        root / "data/small/32-data-NOTICE.txt": frozen / "data-NOTICE.txt",
        root / "env/hybrid-assembly.yml": frozen / "env/hybrid-assembly.yml",
        root / "env/hybrid-assembly-linux-64.lock": frozen / "env/hybrid-assembly-linux-64.lock",
        root / "env/long-read-assembly-linux-64.lock": frozen / "env/long-read-assembly-linux-64.lock",
        raw / "selected/ERR9765746_selection-summary.json": frozen / "selection-summary.json",
        work / "fastp.json": frozen / "fastp.json",
        work / "tool-versions.tsv": frozen / "tool-versions.tsv",
        work / "truth/truth-manifest.tsv": frozen / "truth-manifest.tsv",
        work / "truth/truth-audit.json": frozen / "truth-audit.json",
        work / "truth/published-anchor.tsv": frozen / "published-anchor.tsv",
        work / "summary/branch-metrics.tsv": frozen / "branch-metrics.tsv",
        work / "summary/per-genome-metaquast.tsv": frozen / "per-genome-metaquast.tsv",
        work / "summary/abundance-bin-recovery.tsv": frozen / "abundance-bin-recovery.tsv",
        work / "summary/polishing-sequence-audit.tsv": frozen / "polishing-sequence-audit.tsv",
        work / "summary/polypolish-log-audit.tsv": frozen / "polypolish-log-audit.tsv",
        work / "summary/split-scaffold-sensitivity.tsv": frozen / "split-scaffold-sensitivity.tsv",
        work / "summary/resource-usage.tsv": frozen / "resource-usage.tsv",
        work / "summary/run-summary.json": frozen / "computation-summary.json",
        raw / "benchmark_mock/script_r/Supplementary_Table_S1.xlsx": frozen / "sources/Supplementary_Table_S1.xlsx",
        raw / "Supplementary_Table_S2.xlsx": frozen / "sources/Supplementary_Table_S2.xlsx",
        raw / "benchmark_mock/LICENSE": frozen / "sources/benchmark-repository-LICENSE.txt",
        raw / "benchmark_mock/reference/MOCK_001.fasta.gz": frozen / "references/MOCK_001.fasta.gz",
        work / "metaquast/combined_reference/transposed_report.tsv": frozen / "metaquast/combined-reference-transposed-report.tsv",
    }
    for source, target in direct.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

    script_names = (
        "download_article32_hybrid_sources.sh",
        "select_article32_read_pairs.py",
        "prepare_article32_truth.py",
        "prepare_article32_assemblies.py",
        "run_article32_hybrid_assembly.sh",
        "summarize_article32_metaquast.py",
        "freeze_article32_hybrid_assembly.py",
        "validate_article32_hybrid_assembly.py",
    )
    for name in script_names:
        source = root / "scripts" / name
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, frozen / "scripts" / name)

    source_rows = read_tsv(root / "data/small/32-source-manifest.tsv")
    source_audit: list[dict[str, object]] = []
    for row in source_rows:
        path = source_path(root, row)
        observed = hashes(path, ("md5", "sha256"))
        passed = (
            path.stat().st_size == int(row["ENABytes"])
            and observed["md5"] == row["ENAReportedMD5"]
            and observed["sha256"] == row["ObservedSHA256"]
        )
        audit = {
            "RunAccession": row["RunAccession"],
            "Mate": row["Mate"],
            "ObservedBytes": path.stat().st_size,
            "ObservedMD5": observed["md5"],
            "ObservedSHA256": observed["sha256"],
            "IdentityStatus": "PASS" if passed else "FAIL",
        }
        if not passed:
            raise ValueError(f"Source identity failed: {audit}")
        source_audit.append(audit)
    write_tsv(frozen / "source-audit.tsv", source_audit, list(source_audit[0]))

    clean_rows: list[dict[str, object]] = []
    for mate in ("R1", "R2"):
        clean_path = raw / "clean" / f"ERR9765746_clean10m_{mate}.fastq.gz"
        if not clean_path.is_file():
            raise FileNotFoundError(clean_path)
        clean_rows.append({"Mate": mate, **count_fastq(clean_path)})
    if clean_rows[0]["Records"] != clean_rows[1]["Records"]:
        raise ValueError("Clean FASTQ mate record counts differ")
    if clean_rows[0]["PairIDHash"] != clean_rows[1]["PairIDHash"]:
        raise ValueError("Clean FASTQ mate IDs differ")
    fastp = json.loads((work / "fastp.json").read_text(encoding="utf-8"))
    fastp_after = fastp["summary"]["after_filtering"]
    if int(clean_rows[0]["Records"]) * 2 != int(fastp_after["total_reads"]):
        raise ValueError("Clean FASTQ record count differs from fastp JSON")
    if sum(int(row["Bases"]) for row in clean_rows) != int(fastp_after["total_bases"]):
        raise ValueError("Clean FASTQ base count differs from fastp JSON")
    write_tsv(frozen / "clean-fastq-audit.tsv", clean_rows, list(clean_rows[0]))

    truth_rows = read_tsv(work / "truth/truth-manifest.tsv")
    if len(truth_rows) != 71:
        raise ValueError(f"Expected 71 truth genomes, found {len(truth_rows)}")
    for row in truth_rows:
        source = work / "truth/references/MOCK1" / f"{row['CurrentGenomeLabel']}.fna.gz"
        target = frozen / "references/genomes" / source.name
        shutil.copy2(source, target)
        if sha256(target) != row["CompressedSHA256"]:
            raise ValueError(f"Truth reference digest failed: {row['CurrentGenomeLabel']}")

    frozen_assemblies: list[dict[str, object]] = []
    metrics = {row["Branch"]: row for row in read_tsv(work / "summary/branch-metrics.tsv")}
    if set(metrics) != set(BRANCHES):
        raise ValueError(f"Unexpected branch set: {sorted(metrics)}")
    for branch in BRANCHES:
        source = work / "normalized/final" / f"{branch}.ge1000.fasta"
        target = frozen / "assemblies" / f"{branch}.ge1000.fna.gz"
        records = list(fasta_records(source))
        write_deterministic_gzip(records, target)
        observed_contigs = len(records)
        observed_bases = sum(len(sequence) for _, sequence in records)
        if observed_contigs != int(metrics[branch]["SequencesGe1kb"]):
            raise ValueError(f"Frozen contig count differs for {branch}")
        if observed_bases != int(metrics[branch]["TotalLengthBp"]):
            raise ValueError(f"Frozen base count differs for {branch}")
        frozen_assemblies.append(
            {
                "Branch": branch,
                "Contigs": observed_contigs,
                "Bases": observed_bases,
                "CompressedBytes": target.stat().st_size,
                "CompressedSHA256": sha256(target),
            }
        )
    write_tsv(
        frozen / "frozen-assemblies.tsv",
        frozen_assemblies,
        list(frozen_assemblies[0]),
    )

    report_count = 0
    for reference_dir in sorted((work / "metaquast/runs_per_reference").iterdir()):
        report = reference_dir / "transposed_report.tsv"
        if report.is_file():
            shutil.copy2(report, frozen / "metaquast/references" / f"{reference_dir.name}.tsv")
            report_count += 1
    if not 1 <= report_count <= 71:
        raise ValueError(f"Expected 1..71 physical per-reference reports, found {report_count}")

    for path in sorted((work / "logs").glob("*.log")):
        normalize_copy(path, frozen / "logs" / path.name, replacements)
    for path in sorted((work / "resources").glob("*.txt")):
        normalize_copy(path, frozen / "resources" / path.name, replacements)

    summary = {
        "article": 32,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "biological_sample": "SAMEA14435832 / MOCK1 / 71 strains",
        "seed": 20260732,
        "selected_illumina_pairs": 10_000_000,
        "source_archives_verified": len(source_audit),
        "clean_fastq_audit": {row["Mate"]: row for row in clean_rows},
        "truth_genomes": len(truth_rows),
        "branches": list(BRANCHES),
        "per_reference_reports": report_count,
        "interpretation": "assembly-strategy audit; hybrid and polishing gains require truth-aware error checks",
    }
    (frozen / "run-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    checksum_lines = []
    for path in sorted(frozen.rglob("*")):
        if path.is_file() and path.name != "file-checksums.sha256":
            checksum_lines.append(f"{sha256(path)}  {path.relative_to(frozen).as_posix()}")
    (frozen / "file-checksums.sha256").write_text(
        "\n".join(checksum_lines) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
