#!/usr/bin/env python3
"""Freeze checksum-locked, compact evidence for Article 31."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import re
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, TextIO


BRANCHES = (
    ("flye-ont-r9", "Flye", "ONT R9", "assemblies/flye-ont-r9/assembly.fasta"),
    ("flye-hifi", "Flye", "PacBio HiFi", "assemblies/flye-hifi/assembly.fasta"),
    (
        "hifiasm-meta-hifi",
        "hifiasm-meta",
        "PacBio HiFi",
        "normalized/hifiasm-meta-hifi.primary.fasta",
    ),
    (
        "metamdbg-hifi",
        "metaMDBG",
        "PacBio HiFi",
        "assemblies/metamdbg-hifi/contigs.fasta.gz",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    return parser.parse_args()


def hash_file(path: Path, algorithms: tuple[str, ...]) -> dict[str, str]:
    digests = {name: hashlib.new(name) for name in algorithms}
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            for digest in digests.values():
                digest.update(chunk)
    return {name: digest.hexdigest() for name, digest in digests.items()}


def sha256(path: Path) -> str:
    return hash_file(path, ("sha256",))["sha256"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


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
                    raise ValueError(f"FASTA sequence before header: {path}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def write_filtered_gzip(records: list[tuple[str, str]], path: Path, threshold: int) -> None:
    raw = path.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0)
    text = io.TextIOWrapper(zipped, encoding="utf-8", newline="\n")
    try:
        for name, sequence in records:
            if len(sequence) < threshold:
                continue
            text.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                text.write(sequence[start : start + 80] + "\n")
    finally:
        text.close()
        raw.close()


def nx(lengths: list[int], fraction: float = 0.5) -> int:
    if not lengths:
        return 0
    target = sum(lengths) * fraction
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            return length
    raise AssertionError("Nx invariant failed")


def count_fastq(path: Path) -> dict[str, object]:
    records = 0
    bases = 0
    minimum: int | None = None
    maximum = 0
    lengths: list[int] = []
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
            if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                raise ValueError(f"Invalid FASTQ record {records + 1}: {path}")
            length = len(sequence)
            records += 1
            bases += length
            minimum = length if minimum is None else min(minimum, length)
            maximum = max(maximum, length)
            lengths.append(length)
    return {
        "ObservedReadCount": records,
        "ObservedBaseCount": bases,
        "ObservedMinimumReadLengthBp": minimum or 0,
        "ObservedMeanReadLengthBp": bases / records if records else 0.0,
        "ObservedReadN50Bp": nx(lengths),
        "ObservedMaximumReadLengthBp": maximum,
    }


def parse_elapsed(value: str) -> float:
    parts = value.split(":")
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    raise ValueError(f"Unrecognized elapsed time: {value}")


def parse_resource(path: Path) -> dict[str, object]:
    wanted = {
        "User time (seconds)": "UserSeconds",
        "System time (seconds)": "SystemSeconds",
        "Elapsed (wall clock) time (h:mm:ss or m:ss)": "WallSeconds",
        "Maximum resident set size (kbytes)": "PeakRSSKiB",
        "Exit status": "ExitStatus",
    }
    found: dict[str, object] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        for source, target in wanted.items():
            prefix = source + ":"
            if line.startswith(prefix):
                value = line[len(prefix) :].strip()
                if target == "WallSeconds":
                    found[target] = parse_elapsed(value)
                elif target in {"PeakRSSKiB", "ExitStatus"}:
                    found[target] = int(value)
                else:
                    found[target] = float(value)
    missing = set(wanted.values()) - set(found)
    if missing:
        raise ValueError(f"Missing resource fields in {path}: {sorted(missing)}")
    return found


def normalize_and_copy(source: Path, target: Path, replacements: dict[str, str]) -> None:
    text = source.read_text(encoding="utf-8", errors="replace")
    for old, new in sorted(replacements.items(), key=lambda item: -len(item[0])):
        text = text.replace(old, new)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")


def graph_summary(branch: str, assembler: str, path: Path) -> dict[str, object]:
    if not path.is_file():
        return {
            "Branch": branch,
            "Assembler": assembler,
            "GraphFilePresent": "FALSE",
            "Segments": 0,
            "Links": 0,
            "Paths": 0,
            "CircularNamedSegments": 0,
        }
    segments = links = paths = circular = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("S\t"):
                segments += 1
                fields = raw.split("\t", 3)
                circular += int(len(fields) > 1 and fields[1].endswith("c"))
            elif raw.startswith("L\t"):
                links += 1
            elif raw.startswith("P\t"):
                paths += 1
    return {
        "Branch": branch,
        "Assembler": assembler,
        "GraphFilePresent": "TRUE",
        "Segments": segments,
        "Links": links,
        "Paths": paths,
        "CircularNamedSegments": circular,
    }


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
    for directory in ("assemblies", "logs", "resources", "scripts", "env"):
        (frozen / directory).mkdir()

    replacements = {
        str(work): "${ARTICLE31_WORK_DIR}",
        str(raw): "${ARTICLE31_RAW_DIR}",
        str(root): "${PROJECT_ROOT}",
        str(env_prefix): "${LONG_READ_ENV_PREFIX}",
        str(Path.home()): "${HOME}",
    }

    direct_copies = {
        root / "data/small/31-source-manifest.tsv": frozen / "source-manifest.tsv",
        root / "data/small/31-data-NOTICE.txt": frozen / "data-NOTICE.txt",
        root / "data/small/31-software-releases.tsv": frozen / "software-releases.tsv",
        root / "env/long-read-assembly.yml": frozen / "env/long-read-assembly.yml",
        root / "env/long-read-assembly-linux-64.lock": frozen
        / "env/long-read-assembly-linux-64.lock",
        raw / "tools/hifiasm-meta-source.tsv": frozen / "hifiasm-meta-source.tsv",
        work / "tool-versions.tsv": frozen / "tool-versions.tsv",
        work / "normalized/assembly-metrics.tsv": frozen / "assembly-metrics.tsv",
        work / "normalized/contig-inventory.tsv": frozen / "contig-inventory.tsv",
        work / "normalized/circular-candidates.tsv": frozen / "circular-candidates.tsv",
        work / "normalized/junction-support.tsv": frozen / "junction-support.tsv",
    }
    for source, target in direct_copies.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, target)

    script_names = (
        "run_article31_long_read_assembly.sh",
        "download_article31_long_reads.sh",
        "bootstrap_article31_hifiasm_meta.sh",
        "prepare_article31_assemblies.py",
        "summarize_article31_paf.py",
        "summarize_article31_junctions.py",
        "freeze_article31_long_read_assembly.py",
        "validate_article31_long_read_assembly.py",
    )
    for name in script_names:
        shutil.copy2(root / "scripts" / name, frozen / "scripts" / name)

    source_rows = read_tsv(root / "data/small/31-source-manifest.tsv")
    source_audit: list[dict[str, object]] = []
    for row in source_rows:
        path = raw / "full" / f"{row['RunAccession']}.fastq.gz"
        if not path.is_file():
            raise FileNotFoundError(path)
        identities = hash_file(path, ("md5", "sha256"))
        metrics = count_fastq(path)
        observed = {
            **row,
            "ObservedBytes": path.stat().st_size,
            "ObservedMD5": identities["md5"],
            "ObservedSHA256": identities["sha256"],
            **metrics,
        }
        observed["IdentityStatus"] = (
            "PASS"
            if int(row["ENABytes"]) == observed["ObservedBytes"]
            and row["ENAReportedMD5"] == observed["ObservedMD5"]
            and int(row["ENAReadCount"]) == observed["ObservedReadCount"]
            and int(row["ENABaseCount"]) == observed["ObservedBaseCount"]
            else "FAIL"
        )
        if observed["IdentityStatus"] != "PASS":
            raise ValueError(f"Source audit failed: {row['RunAccession']}")
        source_audit.append(observed)
    write_tsv(
        frozen / "source-audit.tsv",
        source_audit,
        list(source_audit[0].keys()),
    )

    benchmark = [
        row
        for row in read_tsv(root / "data/small/08-platform-benchmark.tsv")
        if row["PlatformKey"] in {"ONT", "PacBio"}
    ]
    write_tsv(frozen / "published-anchor.tsv", benchmark, list(benchmark[0].keys()))

    frozen_assembly_rows: list[dict[str, object]] = []
    for branch, assembler, platform, relative in BRANCHES:
        source = work / relative
        records = list(fasta_records(source))
        target = frozen / "assemblies" / f"{branch}.ge1000.fna.gz"
        write_filtered_gzip(records, target, 1_000)
        retained = [(name, sequence) for name, sequence in records if len(sequence) >= 1_000]
        frozen_assembly_rows.append(
            {
                "Branch": branch,
                "Assembler": assembler,
                "Platform": platform,
                "Filter": "length>=1000 bp",
                "RetainedContigs": len(retained),
                "RetainedBases": sum(len(sequence) for _, sequence in retained),
                "CompressedBytes": target.stat().st_size,
                "CompressedSHA256": sha256(target),
            }
        )
    write_tsv(
        frozen / "frozen-assemblies.tsv",
        frozen_assembly_rows,
        list(frozen_assembly_rows[0].keys()),
    )

    mapping_rows: list[dict[str, object]] = []
    for branch, _, _, _ in BRANCHES:
        payload = json.loads((work / "mapping" / f"{branch}.json").read_text(encoding="utf-8"))
        mapping_rows.append(payload)
    write_tsv(frozen / "readback-metrics.tsv", mapping_rows, list(mapping_rows[0].keys()))

    attempt_rows: list[dict[str, object]] = []
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    pattern = re.compile(r"^(?P<label>.+)\.attempt(?P<attempt>[0-9]{3})\.txt$")
    for path in sorted((work / "resources").glob("*.attempt*.txt")):
        match = pattern.fullmatch(path.name)
        if not match:
            raise ValueError(f"Unrecognized resource filename: {path.name}")
        metrics = parse_resource(path)
        row = {
            "Step": match.group("label"),
            "Attempt": int(match.group("attempt")),
            **metrics,
        }
        attempt_rows.append(row)
        grouped[str(row["Step"])].append(row)
        normalize_and_copy(path, frozen / "resources" / path.name, replacements)
    if not attempt_rows:
        raise ValueError("No resource records found")
    write_tsv(
        frozen / "resource-attempts.tsv",
        attempt_rows,
        ["Step", "Attempt", "UserSeconds", "SystemSeconds", "WallSeconds", "PeakRSSKiB", "ExitStatus"],
    )
    aggregate_rows: list[dict[str, object]] = []
    for step, rows in sorted(grouped.items()):
        final = rows[-1]
        aggregate_rows.append(
            {
                "Step": step,
                "Attempts": len(rows),
                "UserSeconds": sum(float(row["UserSeconds"]) for row in rows),
                "SystemSeconds": sum(float(row["SystemSeconds"]) for row in rows),
                "WallSeconds": sum(float(row["WallSeconds"]) for row in rows),
                "PeakRSSKiB": max(int(row["PeakRSSKiB"]) for row in rows),
                "FinalUserSeconds": float(final["UserSeconds"]),
                "FinalSystemSeconds": float(final["SystemSeconds"]),
                "FinalWallSeconds": float(final["WallSeconds"]),
                "FinalPeakRSSKiB": int(final["PeakRSSKiB"]),
                "FinalExitStatus": int(final["ExitStatus"]),
                "AllAttemptsPassed": "TRUE"
                if all(int(row["ExitStatus"]) == 0 for row in rows)
                else "FALSE",
            }
        )
        if int(rows[-1]["ExitStatus"]) != 0:
            raise ValueError(f"Final resource attempt did not pass: {step}")
    write_tsv(
        frozen / "resource-usage.tsv",
        aggregate_rows,
        list(aggregate_rows[0].keys()),
    )

    for path in sorted((work / "logs").glob("*.log")):
        normalize_and_copy(path, frozen / "logs" / path.name, replacements)
    native_logs = {
        work / "assemblies/flye-ont-r9/flye.log": frozen / "logs/flye-ont-r9.native.log",
        work / "assemblies/flye-hifi/flye.log": frozen / "logs/flye-hifi.native.log",
        work / "assemblies/metamdbg-hifi/metaMDBG.log": frozen / "logs/metamdbg-hifi.native.log",
    }
    for source, target in native_logs.items():
        if source.is_file():
            normalize_and_copy(source, target, replacements)

    graph_rows = [
        graph_summary(
            "flye-ont-r9",
            "Flye",
            work / "assemblies/flye-ont-r9/assembly_graph.gfa",
        ),
        graph_summary(
            "flye-hifi",
            "Flye",
            work / "assemblies/flye-hifi/assembly_graph.gfa",
        ),
        graph_summary(
            "hifiasm-meta-hifi",
            "hifiasm-meta",
            work / "assemblies/hifiasm-meta-hifi/asm.p_ctg.gfa",
        ),
        graph_summary(
            "metamdbg-hifi",
            "metaMDBG",
            work / "assemblies/metamdbg-hifi/assembly_graph.gfa",
        ),
    ]
    write_tsv(frozen / "assembly-graph-summary.tsv", graph_rows, list(graph_rows[0].keys()))

    summary = {
        "article": 31,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "biological_sample": "SAMEA14435832 / MOCK1 / 71 strains",
        "source_runs": [row["RunAccession"] for row in source_rows],
        "branches": [row[0] for row in BRANCHES],
        "source_archives_verified": len(source_audit),
        "assembly_payloads": len(frozen_assembly_rows),
        "mapping_summaries": len(mapping_rows),
        "resource_steps": len(aggregate_rows),
        "interpretation": "platform-adapted workflow audit; not a platform-only or universal assembler ranking",
    }
    (frozen / "run-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    checksum_rows: list[str] = []
    for path in sorted(frozen.rglob("*")):
        if path.is_file() and path.name != "file-checksums.sha256":
            checksum_rows.append(f"{sha256(path)}  {path.relative_to(frozen).as_posix()}")
    (frozen / "file-checksums.sha256").write_text(
        "\n".join(checksum_rows) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
