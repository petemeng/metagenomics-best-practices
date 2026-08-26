#!/usr/bin/env python3
"""Small standard-library helpers shared by Articles 37–40 audit scripts."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Iterable


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    rows = list(rows)
    if not rows and fieldnames is None:
        raise ValueError(f"Refusing schema-free empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        fields = fieldnames or list(rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_time(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    values: dict[str, str] = {}
    for line in text.splitlines():
        if ":" in line:
            key, value = line.strip().split(":", 1)
            values[key.strip()] = value.strip()
    match = re.search(r"Elapsed \(wall clock\) time \([^)]*\):\s*(\S+)", text)
    elapsed = match.group(1) if match else ""
    rss_kb = int(values.get("Maximum resident set size (kbytes)", "0"))
    return {
        "Label": path.name.removesuffix(".time.txt"),
        "Elapsed": elapsed,
        "PeakRAMGiB": round(rss_kb / 1024 / 1024, 3),
        "FileSystemInputs": int(values.get("File system inputs", "0")),
        "FileSystemOutputs": int(values.get("File system outputs", "0")),
        "ExitStatus": int(values.get("Exit status", "-1")),
    }


def read_abundance(root: Path) -> tuple[dict[str, dict[str, int]], dict[str, int]]:
    path = root / "data/small/35-gene-abundance-frozen/gene-abundance-long.tsv.gz"
    counts: dict[str, dict[str, int]] = {}
    totals: dict[str, int] = {}
    for row in read_tsv(path):
        sample = row["Sample"]
        counts.setdefault(row["GeneID"], {})[sample] = int(row["RawCount"])
        totals[sample] = totals.get(sample, 0) + int(row["RawCount"])
    return counts, totals


def read_metadata(root: Path) -> list[dict[str, str]]:
    return read_tsv(root / "data/small/34-nonredundant-gene-catalog-frozen/primary-catalog-representatives.tsv.gz")


def dump_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def freeze_bundle(
    *, root: Path, work: Path, output: Path, article: int, slug: str,
    source_files: list[Path], selected_raw: list[tuple[Path, Path]] | None = None,
) -> None:
    """Copy compact evidence and source contracts, then checksum every frozen file."""
    if not (work / f".article{article}-summary-complete").is_file():
        raise FileNotFoundError(f"Article {article} summary is incomplete")
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)
    for path in sorted((work / "summary").iterdir()):
        if path.is_file():
            shutil.copy2(path, output / path.name)
    for name in ("command-log.tsv", "database-audit.tsv", "input-audit.tsv", "input-lineage.tsv", "run-contract.json", "tool-versions.tsv"):
        source = work / name
        if source.is_file():
            shutil.copy2(source, output / name)
    logs_out = output / "logs"
    logs_out.mkdir()
    for path in sorted((work / "logs").iterdir()):
        if path.is_file():
            shutil.copy2(path, logs_out / path.name)
    scripts_out = output / "scripts"
    scripts_out.mkdir()
    for path in source_files:
        relative = path.resolve().relative_to(root.resolve())
        if relative.parts[0] == "scripts":
            target = scripts_out / path.name
        elif relative.parts[0] == "env":
            target = output / "env" / path.name
        else:
            target = output / path.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    for source, relative in selected_raw or []:
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    contract = {
        "article": article, "slug": slug, "summary_complete": True,
        "checksum_algorithm": "SHA-256", "checksum_scope": "all frozen files except file-checksums.sha256",
    }
    dump_json(output / "frozen-contract.json", contract)
    (output / "data-NOTICE.txt").write_text(
        f"Article {article} frozen evidence for {slug}.\n"
        "Large databases, FASTA inputs, and transient search indexes are excluded.\n"
        "Input identities, commands, versions, compact outputs, resource logs, and every included file are checksum-covered.\n",
        encoding="utf-8",
    )
    rows = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "file-checksums.sha256":
            rows.append(f"{sha256(path)}  {path.relative_to(output).as_posix()}")
    (output / "file-checksums.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")
