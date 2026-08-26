#!/usr/bin/env python3
"""Shared deterministic I/O helpers for the genome-resolved Articles 41–44."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import shutil
from pathlib import Path
from typing import Iterable, Iterator


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path,
    rows: Iterable[dict[str, object]],
    fieldnames: list[str] | tuple[str, ...] | None = None,
) -> None:
    rows = list(rows)
    if not rows and fieldnames is None:
        raise ValueError(f"Refusing schema-free empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        zipped = gzip.GzipFile(filename=str(path), mode="wb", compresslevel=9, mtime=0)
        handle = io.TextIOWrapper(zipped, encoding="utf-8", newline="")
    else:
        handle = path.open("w", encoding="utf-8", newline="")
    with handle:
        fields = list(fieldnames or rows[0].keys())
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def dump_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def parse_elapsed(value: str) -> float:
    parts = value.split(":")
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    raise ValueError(f"Unrecognized elapsed time: {value}")


def parse_time(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Elapsed \(wall clock\) time \([^)]*\):\s*(\S+)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    exit_status = re.search(r"Exit status:\s*(-?\d+)", text)
    command = re.search(r'Command being timed: "(.*)"', text)
    return {
        "Label": path.name.removesuffix(".time.txt"),
        "WallSeconds": round(parse_elapsed(match.group(1)), 3) if match else math.nan,
        "PeakRAMGiB": round(int(rss.group(1)) / 1024 / 1024, 3) if rss else math.nan,
        "ExitStatus": int(exit_status.group(1)) if exit_status else -1,
        "Command": command.group(1) if command else "",
        "MeasurementStatus": "GNU time exact",
    }


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    name: str | None = None
    chunks: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"FASTA sequence before header: {path}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def fasta_summary(path: Path) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    per_contig: dict[str, dict[str, object]] = {}
    lengths: list[int] = []
    gc_bases = 0
    acgt_bases = 0
    for name, sequence in fasta_records(path):
        if name in per_contig:
            raise ValueError(f"Duplicate FASTA identifier: {name}")
        length = len(sequence)
        gc = sequence.count("G") + sequence.count("C")
        acgt = sum(sequence.count(base) for base in "ACGT")
        per_contig[name] = {
            "Contig": name,
            "LengthBp": length,
            "GCPct": 100 * gc / acgt if acgt else math.nan,
        }
        lengths.append(length)
        gc_bases += gc
        acgt_bases += acgt
    if not lengths:
        raise ValueError(f"No FASTA records: {path}")
    ordered = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    cumulative = 0
    n50 = 0
    for length in ordered:
        cumulative += length
        if cumulative >= half:
            n50 = length
            break
    summary = {
        "Contigs": len(lengths),
        "TotalBp": sum(lengths),
        "MinimumBp": min(lengths),
        "MaximumBp": max(lengths),
        "N50Bp": n50,
        "GCPct": 100 * gc_bases / acgt_bases if acgt_bases else math.nan,
    }
    return summary, per_contig


def normalize_copy(source: Path, target: Path, replacements: dict[str, str]) -> None:
    text = source.read_text(encoding="utf-8", errors="replace")
    for old, new in sorted(replacements.items(), key=lambda item: -len(item[0])):
        text = text.replace(old, new)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")


def checksum_tree(root: Path) -> None:
    lines = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "file-checksums.sha256":
            lines.append(f"{sha256(path)}  {path.relative_to(root).as_posix()}")
    (root / "file-checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def freeze_compact_bundle(
    *,
    root: Path,
    work: Path,
    output: Path,
    article: int,
    slug: str,
    source_files: list[Path],
    extra_files: list[tuple[Path, Path]] | None = None,
) -> None:
    """Freeze summaries, normalized logs, contracts and source without large intermediates."""
    sentinel = work / f".article{article}-summary-complete"
    if not sentinel.is_file():
        raise FileNotFoundError(f"Missing summary sentinel: {sentinel}")
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)
    replacements = {
        str(root.resolve()): "${PROJECT_ROOT}",
        str(work.resolve()): f"${{ARTICLE{article}_WORK_DIR}}",
        str(Path.home().resolve()): "${HOME}",
    }
    for path in sorted((work / "summary").iterdir()):
        if path.is_file():
            shutil.copy2(path, output / path.name)
    for name in (
        "command-log.tsv",
        "qc-command-log.tsv",
        "database-audit.tsv",
        "qc-database-audit.tsv",
        "qc-domain-audit.tsv",
        "input-audit.tsv",
        "input-lineage.tsv",
        "input-binsets.tsv",
        "candidate-reconstruction-audit.tsv",
        "refined-bin-source.tsv",
        "selected-mag-reconstruction-audit.tsv",
        "resource-overrides.tsv",
        "run-contract.json",
        "tool-versions.tsv",
    ):
        path = work / name
        if path.is_file():
            normalize_copy(path, output / name, replacements)
    logs_out = output / "logs"
    logs_out.mkdir()
    for path in sorted((work / "logs").iterdir()):
        if path.is_file():
            normalize_copy(path, logs_out / path.name, replacements)
    scripts_out = output / "scripts"
    scripts_out.mkdir()
    for path in source_files:
        relative = path.resolve().relative_to(root.resolve())
        if relative.parts[0] == "scripts":
            target = scripts_out / path.name
        elif relative.parts[0] == "env":
            target = output / "env" / path.name
        else:
            target = output / relative.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    for source, relative in extra_files or []:
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    dump_json(
        output / "frozen-contract.json",
        {
            "article": article,
            "slug": slug,
            "summary_complete": True,
            "checksum_algorithm": "SHA-256",
            "checksum_scope": "all frozen files except file-checksums.sha256",
            "large_intermediates_excluded": True,
        },
    )
    (output / "data-NOTICE.txt").write_text(
        f"Article {article} frozen evidence for {slug}.\n"
        "Large FASTQ, BAM, database, index and binning work files stay outside Git.\n"
        "Input identities, commands, versions, compact outputs and normalized logs are checksum-covered.\n",
        encoding="utf-8",
    )
    checksum_tree(output)
