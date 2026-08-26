#!/usr/bin/env python3
"""Freeze compact, checksum-locked evidence for Article 30."""

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
from pathlib import Path
from typing import Iterable, Iterator


BRANCHES = (
    ("megahit-single-MOCK1", "MEGAHIT", "Single", "MOCK1", "final.contigs.fa"),
    ("megahit-single-MOCK2", "MEGAHIT", "Single", "MOCK2", "final.contigs.fa"),
    ("megahit-coassembly", "MEGAHIT", "Co-assembly", "MOCK1+MOCK2", "final.contigs.fa"),
    ("metaspades-single-MOCK1", "metaSPAdes", "Single", "MOCK1", "contigs.fasta"),
    ("metaspades-single-MOCK2", "metaSPAdes", "Single", "MOCK2", "contigs.fasta"),
    ("metaspades-coassembly", "metaSPAdes", "Co-assembly", "MOCK1+MOCK2", "contigs.fasta"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_archive_identity(path: Path) -> dict[str, object]:
    """Compute both archive digests in one pass."""
    md5_digest = hashlib.md5()
    sha256_digest = hashlib.sha256()
    observed_bytes = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            observed_bytes += len(chunk)
            md5_digest.update(chunk)
            sha256_digest.update(chunk)
    return {
        "ObservedBytes": observed_bytes,
        "ObservedMD5": md5_digest.hexdigest(),
        "ObservedSHA256": sha256_digest.hexdigest(),
    }


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
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


def nx(lengths: list[int], fraction: float) -> tuple[int, int]:
    if not lengths:
        return 0, 0
    target = sum(lengths) * fraction
    cumulative = 0
    for rank, length in enumerate(sorted(lengths, reverse=True), start=1):
        cumulative += length
        if cumulative >= target:
            return length, rank
    raise AssertionError("Nx invariant failed")


def metrics_for(records: list[tuple[str, str]], threshold: int) -> dict[str, object]:
    sequences = [sequence for _, sequence in records if len(sequence) >= threshold]
    lengths = [len(sequence) for sequence in sequences]
    n50, l50 = nx(lengths, 0.5)
    n90, l90 = nx(lengths, 0.9)
    total = sum(lengths)
    gc = sum(sequence.count("G") + sequence.count("C") for sequence in sequences)
    acgt = sum(sum(sequence.count(base) for base in "ACGT") for sequence in sequences)
    return {
        f"ContigsGE{threshold}": len(lengths),
        f"TotalBpGE{threshold}": total,
        f"LargestBpGE{threshold}": max(lengths, default=0),
        f"N50GE{threshold}": n50,
        f"L50GE{threshold}": l50,
        f"N90GE{threshold}": n90,
        f"L90GE{threshold}": l90,
        f"GCPctGE{threshold}": 100 * gc / acgt if acgt else math.nan,
    }


def write_filtered_gzip(
    records: list[tuple[str, str]], path: Path, threshold: int = 1000
) -> None:
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
    for raw in path.read_text(encoding="utf-8").splitlines():
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


def count_fastq(path: Path) -> dict[str, object]:
    records = 0
    bases = 0
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
            if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
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
    return {
        "Records": records,
        "Bases": bases,
        "MinimumLength": minimum,
        "MaximumLength": maximum,
        "MeanLength": bases / records,
        "PairIDHash": pair_ids.hexdigest(),
        "CompressedBytes": path.stat().st_size,
        "CompressedSHA256": sha256(path),
    }


def normalize_and_copy(source: Path, target: Path, replacements: dict[str, str]) -> None:
    text = source.read_text(encoding="utf-8", errors="replace")
    for old, new in sorted(replacements.items(), key=lambda item: -len(item[0])):
        text = text.replace(old, new)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    raw = args.raw_dir.resolve()
    work = args.work_dir.resolve()
    frozen = args.frozen_dir.resolve()
    replacements = {
        str(work): "${ARTICLE30_WORK_DIR}",
        str(raw): "${ARTICLE30_RAW_DIR}",
        str(root): "${PROJECT_ROOT}",
        str(Path.home()): "${HOME}",
    }
    if frozen.exists() and any(frozen.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty frozen directory: {frozen}")
    frozen.mkdir(parents=True, exist_ok=True)
    (frozen / "contigs").mkdir()
    (frozen / "fastp").mkdir()
    (frozen / "logs").mkdir()

    shutil.copy2(root / "data/small/30-source-manifest.tsv", frozen / "source-manifest.tsv")
    shutil.copy2(root / "data/small/30-data-NOTICE.txt", frozen / "data-NOTICE.txt")
    shutil.copy2(root / "scripts/run_article30_short_read_assembly.sh", frozen / "commands.sh")

    source_audit_rows: list[dict[str, object]] = []
    with (root / "data/small/30-source-manifest.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            archive = raw / "full" / f"{row['RunAccession']}_{row['Mate']}.fastq.gz"
            if not archive.is_file():
                raise SystemExit(f"Missing source archive: {archive}")
            observed = source_archive_identity(archive)
            expected_bytes = int(row["ENABytes"])
            expected_md5 = row["ENAReportedMD5"]
            if (
                observed["ObservedBytes"] != expected_bytes
                or observed["ObservedMD5"] != expected_md5
            ):
                raise SystemExit(f"Source archive identity mismatch: {archive}")
            source_audit_rows.append(
                {
                    "Mock": row["Mock"],
                    "RunAccession": row["RunAccession"],
                    "Mate": row["Mate"],
                    "ExpectedBytes": expected_bytes,
                    "ObservedBytes": observed["ObservedBytes"],
                    "ExpectedMD5": expected_md5,
                    "ObservedMD5": observed["ObservedMD5"],
                    "ObservedSHA256": observed["ObservedSHA256"],
                    "Status": "PASS",
                }
            )
    if len(source_audit_rows) != 4:
        raise SystemExit(f"Expected four source archives, found {len(source_audit_rows)}")
    write_tsv(
        frozen / "source-archive-audit.tsv",
        source_audit_rows,
        list(source_audit_rows[0]),
    )

    selection_rows: list[dict[str, object]] = []
    fastp_rows: list[dict[str, object]] = []
    resource_rows: list[dict[str, object]] = []
    clean_audits: dict[str, dict[str, dict[str, object]]] = {}
    for mock, run in (("MOCK1", "ERR9765746"), ("MOCK2", "ERR9765747")):
        selection = json.loads(
            (raw / "selected" / f"{run}_selection-summary.json").read_text()
        )
        shutil.copy2(
            raw / "selected" / f"{run}_selection-summary.json",
            frozen / f"{run}_selection-summary.json",
        )
        selection_rows.append(
            {
                "Mock": mock,
                "RunAccession": run,
                "SourcePairs": selection["source_pairs"],
                "SelectedPairs": selection["selected_pairs"],
                "SelectedFraction": selection["selected_fraction"],
                "Seed": selection["seed"],
                "PairIDHash": selection["selected_pair_id_sha256"],
                "R1Bases": selection["mates"]["R1"]["bases"],
                "R2Bases": selection["mates"]["R2"]["bases"],
            }
        )
        fastp_path = work / "fastp" / f"{mock}.json"
        fastp = json.loads(fastp_path.read_text())
        normalize_and_copy(
            fastp_path,
            frozen / "fastp" / f"{mock}.json",
            replacements,
        )
        before = fastp["summary"]["before_filtering"]
        after = fastp["summary"]["after_filtering"]
        fastp_rows.append(
            {
                "Mock": mock,
                "InputPairs": before["total_reads"] // 2,
                "OutputPairs": after["total_reads"] // 2,
                "InputBases": before["total_bases"],
                "OutputBases": after["total_bases"],
                "PairRetention": after["total_reads"] / before["total_reads"],
                "BaseRetention": after["total_bases"] / before["total_bases"],
                "Q30Before": before["q30_rate"],
                "Q30After": after["q30_rate"],
            }
        )
        resource = parse_resource(work / "resources" / f"fastp-{mock}.txt")
        if resource["ExitStatus"] != 0:
            raise SystemExit(f"Non-zero fastp resource status: {mock}")
        resource_rows.append(
            {"StepType": "Read QC", "Branch": "fastp", "Sample": mock, **resource}
        )
        clean_audits[mock] = {}
        for mate in ("R1", "R2"):
            audit = count_fastq(raw / "clean" / f"{run}_clean_{mate}.fastq.gz")
            clean_audits[mock][mate] = audit
        if clean_audits[mock]["R1"]["Records"] != clean_audits[mock]["R2"]["Records"]:
            raise SystemExit(f"Clean mate count mismatch: {mock}")
        if clean_audits[mock]["R1"]["PairIDHash"] != clean_audits[mock]["R2"]["PairIDHash"]:
            raise SystemExit(f"Clean mate ID mismatch: {mock}")

    write_tsv(
        frozen / "read-selection-summary.tsv",
        selection_rows,
        list(selection_rows[0]),
    )
    write_tsv(frozen / "fastp-summary.tsv", fastp_rows, list(fastp_rows[0]))

    assembly_rows: list[dict[str, object]] = []
    for branch, assembler, strategy, inputs, filename in BRANCHES:
        assembly = work / "assemblies" / branch / filename
        if not assembly.is_file() or assembly.stat().st_size == 0:
            raise SystemExit(f"Missing assembly: {assembly}")
        records = list(fasta_records(assembly))
        if not records:
            raise SystemExit(f"Empty assembly: {assembly}")
        row: dict[str, object] = {
            "Branch": branch,
            "Assembler": assembler,
            "Strategy": strategy,
            "InputMocks": inputs,
        }
        row.update(metrics_for(records, 500))
        row.update(metrics_for(records, 1000))
        row.update(metrics_for(records, 10000))
        row["AmbiguousBasesAll"] = sum(sequence.count("N") for _, sequence in records)
        assembly_rows.append(row)
        write_filtered_gzip(records, frozen / "contigs" / f"{branch}.ge1000.fna.gz")

        resource = parse_resource(work / "resources" / f"assemble-{branch}.txt")
        if resource["ExitStatus"] != 0:
            raise SystemExit(f"Non-zero assembly resource status: {branch}")
        resource_rows.append(
            {"StepType": "Assembly", "Branch": branch, "Sample": "", **resource}
        )
        index_resource = parse_resource(work / "resources" / f"index-{branch}.txt")
        if index_resource["ExitStatus"] != 0:
            raise SystemExit(f"Non-zero index resource status: {branch}")
        resource_rows.append(
            {"StepType": "Indexing", "Branch": branch, "Sample": "", **index_resource}
        )

    write_tsv(frozen / "assembly-metrics.tsv", assembly_rows, list(assembly_rows[0]))

    recruitment_rows: list[dict[str, object]] = []
    for mapping_json in sorted((work / "mapping").glob("*.json")):
        payload = json.loads(mapping_json.read_text())
        branch = payload["assembly_branch"]
        metadata = next(item for item in BRANCHES if item[0] == branch)
        recruitment_rows.append(
            {
                "Sample": payload["sample"],
                "Branch": branch,
                "Assembler": metadata[1],
                "Strategy": metadata[2],
                "TotalPairs": payload["total_pairs"],
                "MappedReadFraction": payload["mapped_read_fraction"],
                "BothMappedPairFraction": payload["both_mapped_pair_fraction"],
                "ProperPairFraction": payload["proper_pair_fraction"],
                "MappedPrimaryReads": payload["mapped_primary_reads"],
                "ProperPairs": payload["proper_pairs"],
                "DiscordantPairs": payload["discordant_pairs"],
                "SingletonPairs": payload["singleton_pairs"],
            }
        )
        stem = mapping_json.stem
        resource = parse_resource(work / "resources" / f"map-{stem}.txt")
        if resource["ExitStatus"] != 0:
            raise SystemExit(f"Non-zero mapping resource status: {stem}")
        resource_rows.append(
            {
                "StepType": "Read mapping",
                "Branch": branch,
                "Sample": payload["sample"],
                **resource,
            }
        )
    if len(recruitment_rows) != 8:
        raise SystemExit(f"Expected eight sample-to-assembly mappings, found {len(recruitment_rows)}")
    write_tsv(
        frozen / "read-recruitment.tsv",
        recruitment_rows,
        list(recruitment_rows[0]),
    )
    write_tsv(frozen / "resource-usage.tsv", resource_rows, list(resource_rows[0]))

    branch_contract = []
    for branch, assembler, strategy, inputs, _ in BRANCHES:
        branch_contract.append(
            {
                "Branch": branch,
                "Assembler": assembler,
                "Strategy": strategy,
                "InputMocks": inputs,
                "MinimumAssemblerContigBp": 500,
                "ReportingThresholdBp": 1000,
                "Threads": 16,
                "MemorySetting": (
                    "-m 64 (64 GiB declared limit)"
                    if assembler == "metaSPAdes"
                    else "--memory 68719476736 (64 GiB in bytes)"
                ),
                "ReadSelectionSeed": 20260730,
                "AssemblerSeedPolicy": "no random-seed option exposed",
                "Bowtie2Seed": 20260730,
            }
        )
    write_tsv(frozen / "branch-contract.tsv", branch_contract, list(branch_contract[0]))

    shutil.copy2(work / "tool-versions.tsv", frozen / "tool-versions.tsv")
    for source in sorted((work / "logs").glob("*.log")):
        normalize_and_copy(source, frozen / "logs" / source.name, replacements)
    for source in sorted((work / "resources").glob("*.txt")):
        normalize_and_copy(source, frozen / "logs" / source.name, replacements)

    summary = {
        "status": "completed",
        "seed": 20260730,
        "source_runs": 2,
        "selected_pairs_per_run": 2_000_000,
        "clean_fastq_audit": clean_audits,
        "assembly_branches": len(assembly_rows),
        "mapping_branches": len(recruitment_rows),
        "reporting_contig_threshold_bp": 1000,
    }
    (frozen / "run-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    checksum_files = sorted(
        path for path in frozen.rglob("*") if path.is_file() and path.name != "file-checksums.sha256"
    )
    with (frozen / "file-checksums.sha256").open("w", encoding="utf-8") as handle:
        for path in checksum_files:
            handle.write(f"{sha256(path)}  {path.relative_to(frozen).as_posix()}\n")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
