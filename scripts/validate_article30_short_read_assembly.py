#!/usr/bin/env python3
"""Validate Article 30 frozen short-read assembly evidence and draw figures."""

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
from typing import Any, Iterable, Iterator

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article30-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image


SEED = 20260730
EXPECTED_HASHES = {
    "source_manifest": "4b9834e1d815388fde59d0f770cd4e49715b5be0658d714042b32f8f6a960fb6",
    "data_notice": "f65ae8c62ae272f86fcebabb623d4d9b0d6c0bc1d58141febf1d35ce15fa30d7",
    "assembly_yaml": "1f35e6307f791fea16591d36044d4c63f51d6c6e4bf6b0a84bd3f795d490cdb3",
    "assembly_lock": "b6254121dd77c6c8c83a99bf98d49f59962f9dcebb6ace145eda27d5ef45e2e5",
    "read_qc_yaml": "3d94a0e5e8a8b1ba2bf58de1d1c15be06824ae7163965559fce48ce632e94759",
    "read_qc_lock": "19bd7992dd9480fbd0fea34fda676da866c3775648a693fc6822636da64f2368",
}
EXPECTED_TOOLS = {
    "fastp": "1.3.6",
    "MEGAHIT": "1.2.9",
    "metaSPAdes": "4.3.0",
    "Bowtie2": "2.5.5",
    "Python": "3.10.20",
}
EXPECTED_SOURCES = {
    ("MOCK1", "R1"): {
        "run": "ERR9765746",
        "sample": "SAMEA14435832",
        "bytes": 1_740_647_656,
        "md5": "ed0e6e0ee846542531c742a45181cd6f",
        "read_count": 41_195_050,
        "base_count": 6_136_710_329,
    },
    ("MOCK1", "R2"): {
        "run": "ERR9765746",
        "sample": "SAMEA14435832",
        "bytes": 2_104_551_765,
        "md5": "5b60ac93cb69dff77ae38cfa501afd06",
        "read_count": 41_195_050,
        "base_count": 6_136_710_329,
    },
    ("MOCK2", "R1"): {
        "run": "ERR9765747",
        "sample": "SAMEA14435833",
        "bytes": 2_141_177_143,
        "md5": "c0e6bb8f83a818f3feef7334cbf50b28",
        "read_count": 46_347_928,
        "base_count": 6_894_690_816,
    },
    ("MOCK2", "R2"): {
        "run": "ERR9765747",
        "sample": "SAMEA14435833",
        "bytes": 2_417_824_585,
        "md5": "fac9f4841f519c72b31149b492479def",
        "read_count": 46_347_928,
        "base_count": 6_894_690_816,
    },
}
BRANCHES = (
    ("megahit-single-MOCK1", "MEGAHIT", "Single", "MOCK1"),
    ("megahit-single-MOCK2", "MEGAHIT", "Single", "MOCK2"),
    ("megahit-coassembly", "MEGAHIT", "Co-assembly", "MOCK1+MOCK2"),
    ("metaspades-single-MOCK1", "metaSPAdes", "Single", "MOCK1"),
    ("metaspades-single-MOCK2", "metaSPAdes", "Single", "MOCK2"),
    ("metaspades-coassembly", "metaSPAdes", "Co-assembly", "MOCK1+MOCK2"),
)
BRANCH_META = {row[0]: row[1:] for row in BRANCHES}
MAPPING_COMBINATIONS = {
    (sample, branch)
    for sample, own in (("MOCK1", "MOCK1"), ("MOCK2", "MOCK2"))
    for branch, _, strategy, inputs in BRANCHES
    if (strategy == "Single" and inputs == own) or strategy == "Co-assembly"
}
FIGURE_STEMS = (
    "30-assembly-branch-design",
    "30-contiguity-output",
    "30-recruitment-tradeoff",
    "30-resource-footprint",
)
PALETTE = {"MEGAHIT": "#0072B2", "metaSPAdes": "#D55E00"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--assembly-prefix", type=Path, required=True)
    parser.add_argument("--read-qc-prefix", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    parser.add_argument("--chapter", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


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


def verify_checksum_manifest(
    frozen: Path, checks: Checks
) -> list[dict[str, str]]:
    manifest = frozen / "file-checksums.sha256"
    rows: list[dict[str, str]] = []
    expected_names: set[str] = set()
    for number, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if match is None:
            checks.add("Frozen input", f"checksum-line-{number}", False, line)
            continue
        expected, relative = match.groups()
        expected_names.add(relative)
        path = frozen / relative
        observed = sha256(path) if path.is_file() else "MISSING"
        status = observed == expected
        checks.add("Frozen input", f"sha256-{relative}", status, observed)
        rows.append(
            {
                "File": relative,
                "ExpectedSHA256": expected,
                "ObservedSHA256": observed,
                "Status": "PASS" if status else "FAIL",
            }
        )
    observed_names = {
        path.relative_to(frozen).as_posix()
        for path in frozen.rglob("*")
        if path.is_file() and path.name != manifest.name
    }
    checks.add(
        "Frozen input",
        "checksum-manifest-complete",
        observed_names == expected_names,
        f"payloads={len(observed_names)} entries={len(expected_names)}",
    )
    return rows


def command_output(command: list[str], path_prefix: str) -> tuple[int, str]:
    process = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env={
            **os.environ,
            "PATH": path_prefix,
            "LC_ALL": "C",
            "TZ": "UTC",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "PYTHONPATH": "",
        },
    )
    return process.returncode, process.stdout.strip()


def audit_tools(
    assembly_prefix: Path, read_qc_prefix: Path, checks: Checks
) -> list[dict[str, Any]]:
    path_prefix = f"{assembly_prefix / 'bin'}:{read_qc_prefix / 'bin'}:/usr/bin:/bin"
    commands = {
        "fastp": [str(read_qc_prefix / "bin/fastp"), "--version"],
        "MEGAHIT": [str(assembly_prefix / "bin/megahit"), "--version"],
        "metaSPAdes": [str(assembly_prefix / "bin/metaspades.py"), "--version"],
        "Bowtie2": [str(assembly_prefix / "bin/bowtie2"), "--version"],
        "Python": [str(assembly_prefix / "bin/python"), "--version"],
    }
    patterns = {
        "fastp": r"fastp\s+([0-9.]+)",
        "MEGAHIT": r"MEGAHIT\s+v([0-9.]+)",
        "metaSPAdes": r"SPAdes genome assembler v([0-9.]+)",
        "Bowtie2": r"(?:Bowtie 2|bowtie2-align-s) version\s+([0-9.]+)",
        "Python": r"Python\s+([0-9.]+)",
    }
    rows: list[dict[str, Any]] = []
    for tool, command in commands.items():
        return_code, output = command_output(command, path_prefix)
        match = re.search(patterns[tool], output)
        observed = match.group(1) if match else "UNPARSED"
        expected = EXPECTED_TOOLS[tool]
        status = return_code == 0 and observed == expected
        checks.add("Toolchain", f"tool-{tool}", status, f"observed={observed}")
        rows.append(
            {
                "Tool": tool,
                "ExpectedVersion": expected,
                "ObservedVersion": observed,
                "ReturnCode": return_code,
                "Status": "PASS" if status else "FAIL",
            }
        )
    return rows


def audit_project_contract(root: Path, frozen: Path, checks: Checks) -> None:
    paths = {
        "source_manifest": root / "data/small/30-source-manifest.tsv",
        "data_notice": root / "data/small/30-data-NOTICE.txt",
        "assembly_yaml": root / "env/assembly.yml",
        "assembly_lock": root / "env/assembly-linux-64.lock",
        "read_qc_yaml": root / "env/read-qc.yml",
        "read_qc_lock": root / "env/read-qc-linux-64.lock",
    }
    for name, path in paths.items():
        observed = sha256(path) if path.is_file() else "MISSING"
        checks.add(
            "Project contract",
            f"sha256-{name}",
            observed == EXPECTED_HASHES[name],
            observed,
        )
    checks.add(
        "Project contract",
        "frozen-manifest-identical",
        (frozen / "source-manifest.tsv").read_bytes()
        == paths["source_manifest"].read_bytes(),
        "frozen versus repository source manifest",
    )
    checks.add(
        "Project contract",
        "frozen-notice-identical",
        (frozen / "data-NOTICE.txt").read_bytes() == paths["data_notice"].read_bytes(),
        "frozen versus repository data notice",
    )
    assembly_yaml = paths["assembly_yaml"].read_text(encoding="utf-8")
    for token in (
        "python=3.10",
        "megahit=1.2.9",
        "spades=4.3.0",
        "bowtie2=2.5.5",
    ):
        checks.add("Project contract", f"assembly-pin-{token}", token in assembly_yaml, token)
    read_qc_yaml = paths["read_qc_yaml"].read_text(encoding="utf-8")
    checks.add(
        "Project contract", "read-qc-pin-fastp", "fastp=1.3.6" in read_qc_yaml, "fastp=1.3.6"
    )
    lock_expectations = {
        "assembly": (paths["assembly_lock"], ("megahit-1.2.9-", "spades-4.3.0-", "bowtie2-2.5.5-")),
        "read-qc": (paths["read_qc_lock"], ("fastp-1.3.6-",)),
    }
    for label, (path, tokens) in lock_expectations.items():
        text = path.read_text(encoding="utf-8")
        checks.add("Project contract", f"{label}-explicit-lock", "@EXPLICIT" in text, path.name)
        for token in tokens:
            checks.add("Project contract", f"{label}-lock-{token}", token in text, token)


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    sequence: list[str] = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence).upper()
                name = line[1:]
                sequence = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                sequence.append(line)
    if name is not None:
        yield name, "".join(sequence).upper()


def nx(lengths: list[int], fraction: float) -> tuple[int, int]:
    target = sum(lengths) * fraction
    cumulative = 0
    for rank, length in enumerate(sorted(lengths, reverse=True), start=1):
        cumulative += length
        if cumulative >= target:
            return length, rank
    return 0, 0


def contig_metrics(records: list[tuple[str, str]]) -> dict[str, float | int]:
    lengths = [len(sequence) for _, sequence in records]
    n50, l50 = nx(lengths, 0.5)
    n90, l90 = nx(lengths, 0.9)
    gc = sum(sequence.count("G") + sequence.count("C") for _, sequence in records)
    acgt = sum(sum(sequence.count(base) for base in "ACGT") for _, sequence in records)
    return {
        "ContigsGE1000": len(lengths),
        "TotalBpGE1000": sum(lengths),
        "LargestBpGE1000": max(lengths, default=0),
        "N50GE1000": n50,
        "L50GE1000": l50,
        "N90GE1000": n90,
        "L90GE1000": l90,
        "GCPctGE1000": 100 * gc / acgt if acgt else math.nan,
        "AmbiguousBasesGE1000": sum(sequence.count("N") for _, sequence in records),
    }


def close_float(left: float, right: float, tolerance: float = 1e-9) -> bool:
    return math.isfinite(left) and math.isfinite(right) and abs(left - right) <= tolerance


def audit_sources(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    manifest_rows = read_tsv(frozen / "source-manifest.tsv")
    audit_rows = read_tsv(frozen / "source-archive-audit.tsv")
    checks.add("Source data", "source-manifest-four-files", len(manifest_rows) == 4, len(manifest_rows))
    checks.add("Source data", "source-audit-four-files", len(audit_rows) == 4, len(audit_rows))
    manifest_by_key = {(row["Mock"], row["Mate"]): row for row in manifest_rows}
    audit_by_key = {(row["Mock"], row["Mate"]): row for row in audit_rows}
    for key, expected in EXPECTED_SOURCES.items():
        manifest = manifest_by_key.get(key, {})
        audit = audit_by_key.get(key, {})
        identity_ok = (
            manifest.get("RunAccession") == expected["run"]
            and manifest.get("SampleAccession") == expected["sample"]
            and int(manifest.get("ENABytes", -1)) == expected["bytes"]
            and manifest.get("ENAReportedMD5") == expected["md5"]
            and int(manifest.get("ENAReadCount", -1)) == expected["read_count"]
            and int(manifest.get("ENABaseCount", -1)) == expected["base_count"]
            and manifest.get("LibraryLayout") == "PAIRED"
        )
        checks.add("Source data", f"manifest-{'-'.join(key)}", identity_ok, expected["run"])
        audit_ok = (
            audit.get("Status") == "PASS"
            and int(audit.get("ExpectedBytes", -1)) == expected["bytes"]
            and int(audit.get("ObservedBytes", -1)) == expected["bytes"]
            and audit.get("ExpectedMD5") == expected["md5"]
            and audit.get("ObservedMD5") == expected["md5"]
            and re.fullmatch(r"[0-9a-f]{64}", audit.get("ObservedSHA256", "")) is not None
        )
        checks.add("Source data", f"archive-audit-{'-'.join(key)}", audit_ok, audit.get("ObservedMD5", "missing"))
    return audit_rows


def audit_selection_and_fastp(
    frozen: Path, checks: Checks
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, int]]:
    selection_rows = read_tsv(frozen / "read-selection-summary.tsv")
    fastp_rows = read_tsv(frozen / "fastp-summary.tsv")
    checks.add("Read contract", "selection-two-runs", len(selection_rows) == 2, len(selection_rows))
    checks.add("Read contract", "fastp-two-runs", len(fastp_rows) == 2, len(fastp_rows))
    expected_source_pairs = {"MOCK1": 20_597_525, "MOCK2": 23_173_964}
    output_pairs: dict[str, int] = {}
    pair_hashes: list[str] = []
    for row in selection_rows:
        mock = row["Mock"]
        summary_path = frozen / f"{row['RunAccession']}_selection-summary.json"
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        selected = int(row["SelectedPairs"])
        source = int(row["SourcePairs"])
        pair_hash = row["PairIDHash"]
        pair_hashes.append(pair_hash)
        selection_ok = (
            mock in expected_source_pairs
            and source == expected_source_pairs[mock]
            and selected == 2_000_000
            and int(row["Seed"]) == SEED
            and re.fullmatch(r"[0-9a-f]{64}", pair_hash) is not None
            and payload["source_pairs"] == source
            and payload["selected_pairs"] == selected
            and payload["seed"] == SEED
            and payload["selected_pair_id_sha256"] == pair_hash
            and payload["random_generator"]
            == "Python random.Random MT19937; string seed version=2"
            and payload["selection_algorithm"]
            == "one-pass exact sequential sampling without replacement"
            and re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", payload["python_version"])
            is not None
        )
        checks.add("Read contract", f"selection-{mock}", selection_ok, f"{selected}/{source}")
    checks.add("Read contract", "selection-namespaces-distinct", len(set(pair_hashes)) == 2, pair_hashes)

    selection_by_mock = {row["Mock"]: row for row in selection_rows}
    for row in fastp_rows:
        mock = row["Mock"]
        input_pairs = int(row["InputPairs"])
        clean_pairs = int(row["OutputPairs"])
        output_pairs[mock] = clean_pairs
        contract_ok = (
            mock in selection_by_mock
            and input_pairs == int(selection_by_mock[mock]["SelectedPairs"])
            and 0 < clean_pairs <= input_pairs
            and close_float(float(row["PairRetention"]), clean_pairs / input_pairs)
            and 0 < float(row["BaseRetention"]) <= 1
            and 0 <= float(row["Q30Before"]) <= 1
            and 0 <= float(row["Q30After"]) <= 1
        )
        checks.add("Read contract", f"fastp-{mock}", contract_ok, f"{clean_pairs}/{input_pairs}")
    return selection_rows, fastp_rows, output_pairs


def audit_assemblies(
    frozen: Path, checks: Checks
) -> tuple[list[dict[str, str]], dict[str, dict[str, float | int]]]:
    rows = read_tsv(frozen / "assembly-metrics.tsv")
    checks.add("Assembly", "six-assembly-branches", len(rows) == 6, len(rows))
    rows_by_branch = {row["Branch"]: row for row in rows}
    checks.add("Assembly", "expected-branch-identities", set(rows_by_branch) == set(BRANCH_META), sorted(rows_by_branch))
    recomputed: dict[str, dict[str, float | int]] = {}
    allowed = set("ACGTRYSWKMBDHVN")
    for branch, (assembler, strategy, inputs) in BRANCH_META.items():
        row = rows_by_branch.get(branch, {})
        metadata_ok = (
            row.get("Assembler") == assembler
            and row.get("Strategy") == strategy
            and row.get("InputMocks") == inputs
        )
        checks.add("Assembly", f"metadata-{branch}", metadata_ok, f"{assembler}/{strategy}/{inputs}")
        if not row:
            continue
        counts = [int(row[f"ContigsGE{threshold}"]) for threshold in (500, 1000, 10000)]
        totals = [int(row[f"TotalBpGE{threshold}"]) for threshold in (500, 1000, 10000)]
        monotone = counts[0] >= counts[1] >= counts[2] >= 0 and totals[0] >= totals[1] >= totals[2] >= 0
        checks.add("Assembly", f"threshold-monotonicity-{branch}", monotone, f"counts={counts}; bp={totals}")
        metric_bounds = True
        for threshold in (500, 1000, 10000):
            count = int(row[f"ContigsGE{threshold}"])
            largest = int(row[f"LargestBpGE{threshold}"])
            n50 = int(row[f"N50GE{threshold}"])
            n90 = int(row[f"N90GE{threshold}"])
            l50 = int(row[f"L50GE{threshold}"])
            l90 = int(row[f"L90GE{threshold}"])
            if count:
                metric_bounds &= largest >= n50 >= n90 >= threshold
                metric_bounds &= 1 <= l50 <= l90 <= count
            metric_bounds &= 0 <= float(row[f"GCPctGE{threshold}"]) <= 100
        checks.add("Assembly", f"metric-bounds-{branch}", metric_bounds, "Nx/Lx/GC")

        contig_path = frozen / "contigs" / f"{branch}.ge1000.fna.gz"
        with contig_path.open("rb") as handle:
            header = handle.read(8)
        checks.add("Assembly", f"deterministic-gzip-{branch}", len(header) == 8 and header[4:8] == b"\x00\x00\x00\x00", header.hex())
        records = list(fasta_records(contig_path))
        headers_unique = len({name for name, _ in records}) == len(records)
        sequence_ok = all(len(sequence) >= 1000 and set(sequence) <= allowed for _, sequence in records)
        checks.add("Assembly", f"contig-records-{branch}", bool(records) and headers_unique and sequence_ok, len(records))
        observed = contig_metrics(records)
        recomputed[branch] = observed
        exact_fields = (
            "ContigsGE1000",
            "TotalBpGE1000",
            "LargestBpGE1000",
            "N50GE1000",
            "L50GE1000",
            "N90GE1000",
            "L90GE1000",
        )
        exact_ok = all(int(row[field]) == int(observed[field]) for field in exact_fields)
        gc_ok = close_float(float(row["GCPctGE1000"]), float(observed["GCPctGE1000"]))
        checks.add("Assembly", f"frozen-contig-recalculation-{branch}", exact_ok and gc_ok, observed)
        ambiguous_ok = int(row["AmbiguousBasesAll"]) >= int(observed["AmbiguousBasesGE1000"])
        checks.add("Assembly", f"ambiguous-base-ledger-{branch}", ambiguous_ok, row["AmbiguousBasesAll"])
    return rows, recomputed


def audit_recruitment(
    frozen: Path, output_pairs: dict[str, int], checks: Checks
) -> list[dict[str, str]]:
    rows = read_tsv(frozen / "read-recruitment.tsv")
    combinations = {(row["Sample"], row["Branch"]) for row in rows}
    checks.add("Recruitment", "eight-mapping-branches", len(rows) == 8, len(rows))
    checks.add("Recruitment", "mapping-combinations", combinations == MAPPING_COMBINATIONS, sorted(combinations))
    for row in rows:
        sample = row["Sample"]
        branch = row["Branch"]
        total_pairs = int(row["TotalPairs"])
        mapped_reads = int(row["MappedPrimaryReads"])
        proper = int(row["ProperPairs"])
        discordant = int(row["DiscordantPairs"])
        singleton = int(row["SingletonPairs"])
        both = proper + discordant
        arithmetic = (
            total_pairs == output_pairs.get(sample)
            and 0 <= proper <= both <= total_pairs
            and 0 <= singleton <= total_pairs - both
            and mapped_reads == 2 * both + singleton
            and close_float(float(row["MappedReadFraction"]), mapped_reads / (2 * total_pairs))
            and close_float(float(row["BothMappedPairFraction"]), both / total_pairs)
            and close_float(float(row["ProperPairFraction"]), proper / total_pairs)
        )
        fractions = all(0 <= float(row[field]) <= 1 for field in ("MappedReadFraction", "BothMappedPairFraction", "ProperPairFraction"))
        metadata = BRANCH_META.get(branch)
        metadata_ok = metadata is not None and row["Assembler"] == metadata[0] and row["Strategy"] == metadata[1]
        checks.add("Recruitment", f"ledger-{sample}-{branch}", arithmetic and fractions and metadata_ok, f"pairs={total_pairs}; mapped_reads={mapped_reads}")
    return rows


def audit_resources(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    rows = read_tsv(frozen / "resource-usage.tsv")
    counts = Counter(row["StepType"] for row in rows)
    expected = Counter({"Read QC": 2, "Assembly": 6, "Indexing": 6, "Read mapping": 8})
    checks.add("Resources", "resource-row-contract", len(rows) == 22 and counts == expected, dict(counts))
    all_valid = all(
        int(row["ExitStatus"]) == 0
        and float(row["UserSeconds"]) >= 0
        and float(row["SystemSeconds"]) >= 0
        and float(row["WallSeconds"]) > 0
        and int(row["PeakRSSKiB"]) > 0
        for row in rows
    )
    checks.add("Resources", "resource-values-valid", all_valid, f"rows={len(rows)}")
    assembly_branches = {row["Branch"] for row in rows if row["StepType"] == "Assembly"}
    index_branches = {row["Branch"] for row in rows if row["StepType"] == "Indexing"}
    checks.add("Resources", "resource-branch-coverage", assembly_branches == set(BRANCH_META) and index_branches == set(BRANCH_META), "assembly and indexing")
    return rows


def audit_privacy_and_chapter(frozen: Path, chapter: Path, checks: Checks) -> None:
    forbidden_suffixes = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".sam", ".bam", ".cram")
    forbidden_files = [path.relative_to(frozen).as_posix() for path in frozen.rglob("*") if path.is_file() and path.name.endswith(forbidden_suffixes)]
    checks.add("Boundaries", "no-read-or-alignment-payloads", not forbidden_files, forbidden_files)
    leaks: list[str] = []
    for path in frozen.rglob("*"):
        if not path.is_file() or path.suffix.lower() in {".gz", ".png", ".pdf", ".tiff"}:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if any(token in text for token in ("/media/desk16/", "/workspace/", "tly9658")):
            leaks.append(path.relative_to(frozen).as_posix())
    checks.add("Boundaries", "no-local-path-leaks", not leaks, leaks)

    text = chapter.read_text(encoding="utf-8")
    chapter_tokens = (
        "对应论文里的哪张图",
        "理论",
        "准备工作",
        "可复制代码",
        "审计与升级",
        "出版级美化",
        "常见坑",
        "这段 Methods 怎么写",
        "换成你自己的数据怎么做",
        "参考",
    )
    checks.add("Chapter", "draft-false", re.search(r"^draft:\s*false\s*$", text, re.M) is not None, chapter.name)
    checks.add("Chapter", "upstream-eval-false", re.search(r"^\s*eval:\s*false\s*$", text, re.M) is not None, chapter.name)
    checks.add("Chapter", "nine-section-contract", all(token in text for token in chapter_tokens), chapter_tokens)
    checks.add("Chapter", "inline-plot-functions", all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")), "pal_pub/theme_pub/save_pub")
    checks.add("Chapter", "no-source-theme-dependency", 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text, "inline only")
    checks.add("Chapter", "hardware-four-dimensions", all(token in text for token in ("RAM", "磁盘", "核", "耗时")), "RAM/disk/cores/time")
    checks.add("Chapter", "no-server-alternative", any(token in text for token in ("HPC", "云", "子集演示")), "HPC/cloud/subset")
    checks.add("Chapter", "seed-visible", str(SEED) in text, SEED)
    checks.add("Chapter", "tool-versions-visible", all(version in text for version in EXPECTED_TOOLS.values() if version != "3.10.20"), EXPECTED_TOOLS)
    checks.add("Chapter", "figure-references", all(f"figures/{stem}.png" in text for stem in FIGURE_STEMS), FIGURE_STEMS)
    meta_phrases = ("本篇可独立跑通", "全系列", "接口只学一次", "作者代码通常长这样", "（即本文）")
    checks.add("Chapter", "no-author-meta-copy", not any(phrase in text for phrase in meta_phrases), meta_phrases)
    boundary_tokens = ("不是算法单因素", "不代表普适赢家", "N50", "正确性", "第 33 篇")
    checks.add("Chapter", "interpretation-boundaries", all(token in text for token in boundary_tokens), boundary_tokens)


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 11.5,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 120,
            "savefig.dpi": 350,
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
        metadata={"Creator": "metagenomics-best-practices", "CreationDate": None, "ModDate": None},
    )
    fig.savefig(
        png_path,
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": "metagenomics-best-practices"},
    )
    with Image.open(png_path) as image:
        image.convert("RGB").save(tiff_path, compression="tiff_lzw", dpi=(350, 350))
    plt.close(fig)


def add_box(
    axis: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    face: str,
    edge: str = "#334E5C",
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.025",
        linewidth=1.1,
        facecolor=face,
        edgecolor=edge,
    )
    axis.add_patch(patch)
    axis.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=9.2)


def add_arrow(axis: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=11,
            linewidth=1,
            color="#607D8B",
            connectionstyle="arc3,rad=0",
        )
    )


def plot_branch_design(figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8), gridspec_kw={"wspace": 0.28})
    left, right = axes
    for axis in axes:
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.axis("off")
    left.set_title("Single-sample assembly\n2M selected pairs per mock", pad=13)
    add_box(left, 0.04, 0.62, 0.22, 0.14, "MOCK1\nERR9765746", "#EAF3F8")
    add_box(left, 0.04, 0.24, 0.22, 0.14, "MOCK2\nERR9765747", "#EAF3F8")
    positions = ((0.64, 0.70, "MEGAHIT"), (0.64, 0.50, "metaSPAdes"), (0.64, 0.32, "MEGAHIT"), (0.64, 0.12, "metaSPAdes"))
    for x, y, label in positions:
        add_box(left, x, y, 0.27, 0.12, label, PALETTE[label] + "22", PALETTE[label])
    for start, end in (
        ((0.26, 0.69), (0.64, 0.76)),
        ((0.26, 0.69), (0.64, 0.56)),
        ((0.26, 0.31), (0.64, 0.38)),
        ((0.26, 0.31), (0.64, 0.18)),
    ):
        add_arrow(left, start, end)
    left.text(0.49, 0.94, "Four branches", ha="center", va="center", color="#334E5C", fontweight="bold")

    right.set_title("Co-assembly\nCompatible libraries pooled as paired files", pad=13)
    add_box(right, 0.04, 0.58, 0.22, 0.14, "MOCK1\n2M pairs", "#EAF3F8")
    add_box(right, 0.04, 0.28, 0.22, 0.14, "MOCK2\n2M pairs", "#EAF3F8")
    add_box(right, 0.37, 0.43, 0.22, 0.14, "Pooled input\n4M pairs", "#F6F0DA", "#9C7A14")
    add_box(right, 0.69, 0.58, 0.27, 0.12, "MEGAHIT", PALETTE["MEGAHIT"] + "22", PALETTE["MEGAHIT"])
    add_box(right, 0.69, 0.30, 0.27, 0.12, "metaSPAdes", PALETTE["metaSPAdes"] + "22", PALETTE["metaSPAdes"])
    add_arrow(right, (0.26, 0.65), (0.37, 0.52))
    add_arrow(right, (0.26, 0.35), (0.37, 0.48))
    add_arrow(right, (0.59, 0.50), (0.69, 0.64))
    add_arrow(right, (0.59, 0.50), (0.69, 0.36))
    right.text(0.50, 0.94, "Two branches", ha="center", va="center", color="#334E5C", fontweight="bold")
    fig.suptitle("Controlled six-branch short-read assembly design", fontsize=13, fontweight="bold", y=1.02)
    fig.text(0.5, -0.01, "Same QC and reporting thresholds; co-assemblies contain twice the selected read budget.", ha="center", fontsize=9, color="#455A64")
    save_figure(fig, figure_dir, FIGURE_STEMS[0])


def ordered_metric_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_branch = {row["Branch"]: row for row in rows}
    order = (
        "megahit-single-MOCK1",
        "metaspades-single-MOCK1",
        "megahit-single-MOCK2",
        "metaspades-single-MOCK2",
        "megahit-coassembly",
        "metaspades-coassembly",
    )
    return [by_branch[branch] for branch in order]


def plot_contiguity(rows: list[dict[str, str]], figure_dir: Path) -> None:
    ordered = ordered_metric_rows(rows)
    labels = ["MOCK1\nsingle", "MOCK1\nsingle", "MOCK2\nsingle", "MOCK2\nsingle", "MOCK1+2\nco-assembly", "MOCK1+2\nco-assembly"]
    colors = [PALETTE[row["Assembler"]] for row in ordered]
    x = np.arange(len(ordered))
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.2))
    panels = (
        ("TotalBpGE1000", 1e6, "Assembled sequence ≥1 kb (Mbp)"),
        ("N50GE1000", 1e3, "N50 at ≥1 kb (kbp)"),
        ("ContigsGE10000", 1, "Contigs ≥10 kb (count)"),
    )
    for axis, (field, divisor, ylabel) in zip(axes, panels):
        values = [float(row[field]) / divisor for row in ordered]
        bars = axis.bar(x, values, color=colors, width=0.72)
        axis.set_xticks(x, labels, rotation=35, ha="right")
        axis.set_ylabel(ylabel)
        axis.grid(axis="y", color="#E4EAED", linewidth=0.7)
        axis.set_axisbelow(True)
        maximum = max(values) if values else 1
        axis.set_ylim(0, maximum * 1.17 if maximum else 1)
        for bar, value in zip(bars, values):
            label = f"{value:.2f}" if divisor != 1 else f"{int(value)}"
            axis.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + maximum * 0.025, label, ha="center", va="bottom", fontsize=7.6, rotation=90 if len(label) > 5 else 0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=PALETTE[name]) for name in ("MEGAHIT", "metaSPAdes")]
    fig.legend(handles, ("MEGAHIT", "metaSPAdes"), loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle("Contiguity must be read as a metric panel, not an N50 contest", y=1.12, fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[1])


def plot_recruitment(rows: list[dict[str, str]], figure_dir: Path) -> None:
    lookup = {(row["Sample"], row["Assembler"], row["Strategy"]): row for row in rows}
    groups = (("MOCK1", "MEGAHIT"), ("MOCK1", "metaSPAdes"), ("MOCK2", "MEGAHIT"), ("MOCK2", "metaSPAdes"))
    x = np.arange(len(groups))
    width = 0.36
    fig, axis = plt.subplots(figsize=(8.5, 4.7))
    single = [100 * float(lookup[(sample, assembler, "Single")]["MappedReadFraction"]) for sample, assembler in groups]
    pooled = [100 * float(lookup[(sample, assembler, "Co-assembly")]["MappedReadFraction"]) for sample, assembler in groups]
    bars_single = axis.bar(x - width / 2, single, width, label="Single", color="#56B4E9")
    bars_co = axis.bar(x + width / 2, pooled, width, label="Co-assembly", color="#E69F00")
    labels = [f"{sample}\n{assembler}" for sample, assembler in groups]
    axis.set_xticks(x, labels)
    axis.set_ylabel("Primary reads recruited (%)")
    axis.set_ylim(0, 105)
    axis.grid(axis="y", color="#E4EAED", linewidth=0.7)
    axis.set_axisbelow(True)
    axis.legend(frameon=False, ncol=2, loc="lower right")
    for bars in (bars_single, bars_co):
        for bar in bars:
            axis.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.1, f"{bar.get_height():.1f}", ha="center", va="bottom", fontsize=8)
    axis.set_title("Read recruitment adds evidence that contiguity alone cannot provide", pad=12, fontweight="bold")
    fig.text(0.5, 0.01, "Each sample is mapped to its matching single assembly and to the pooled assembly.", ha="center", fontsize=8.8, color="#455A64")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save_figure(fig, figure_dir, FIGURE_STEMS[2])


def plot_resources(rows: list[dict[str, str]], figure_dir: Path) -> None:
    assembly = [row for row in rows if row["StepType"] == "Assembly"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4))
    labels = []
    for row in assembly:
        assembler, strategy, inputs = BRANCH_META[row["Branch"]]
        short = "M1" if inputs == "MOCK1" else "M2" if inputs == "MOCK2" else "M1+2"
        labels.append(f"{short} {'single' if strategy == 'Single' else 'co'}")
    x = np.arange(len(assembly))
    colors = [PALETTE[BRANCH_META[row["Branch"]][0]] for row in assembly]
    wall = [float(row["WallSeconds"]) / 60 for row in assembly]
    rss = [int(row["PeakRSSKiB"]) / 1024**2 for row in assembly]
    for axis, values, ylabel in (
        (axes[0], wall, "Assembly wall time (min)"),
        (axes[1], rss, "Peak resident memory (GiB)"),
    ):
        bars = axis.bar(x, values, color=colors, width=0.72)
        axis.set_xticks(x, labels, rotation=35, ha="right")
        axis.set_ylabel(ylabel)
        axis.grid(axis="y", color="#E4EAED", linewidth=0.7)
        axis.set_axisbelow(True)
        maximum = max(values) if values else 1
        axis.set_ylim(0, maximum * 1.18 if maximum else 1)
        for bar, value in zip(bars, values):
            axis.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + maximum * 0.025, f"{value:.1f}", ha="center", va="bottom", fontsize=7.8)
    handles = [plt.Rectangle((0, 0), 1, 1, color=PALETTE[name]) for name in ("MEGAHIT", "metaSPAdes")]
    fig.legend(handles, ("MEGAHIT", "metaSPAdes"), loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Resource cost is part of the assembly decision", y=1.10, fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[3])


def audit_figures(figure_dir: Path, checks: Checks) -> None:
    for stem in FIGURE_STEMS:
        pdf = figure_dir / f"{stem}.pdf"
        png = figure_dir / f"{stem}.png"
        tiff = figure_dir / f"{stem}.tiff"
        pdf_ok = pdf.is_file() and pdf.stat().st_size > 1000 and pdf.read_bytes()[:4] == b"%PDF"
        checks.add("Figures", f"figure-{stem}-pdf", pdf_ok, pdf.stat().st_size if pdf.exists() else 0)
        try:
            with Image.open(png) as image:
                png_ok = image.format == "PNG" and image.width >= 1800 and image.height >= 1000
                png_detail = f"{image.format} {image.width}x{image.height} {image.mode}"
        except Exception as error:  # pragma: no cover - failure reporting
            png_ok, png_detail = False, repr(error)
        checks.add("Figures", f"figure-{stem}-png", png_ok, png_detail)
        try:
            with Image.open(tiff) as image:
                compression = str(image.info.get("compression", ""))
                dpi = image.info.get("dpi", (0, 0))
                tiff_ok = (
                    image.format == "TIFF"
                    and image.width >= 1800
                    and image.height >= 1000
                    and "lzw" in compression.lower()
                    and min(float(dpi[0]), float(dpi[1])) >= 349
                )
                tiff_detail = f"{image.format} {image.width}x{image.height} {compression} dpi={dpi}"
        except Exception as error:  # pragma: no cover - failure reporting
            tiff_ok, tiff_detail = False, repr(error)
        checks.add("Figures", f"figure-{stem}-tiff", tiff_ok, tiff_detail)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    assembly_prefix = args.assembly_prefix.resolve()
    read_qc_prefix = args.read_qc_prefix.resolve()
    frozen = (args.frozen_dir or root / "data/small/30-short-read-assembly-frozen").resolve()
    output_dir = (args.output_dir or root / "results/30-short-read-assembly").resolve()
    figure_dir = (args.figure_dir or root / "figures").resolve()
    chapter = (args.chapter or root / "chapters/30-short-read-assembly.qmd").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    required = (
        frozen / "file-checksums.sha256",
        frozen / "source-manifest.tsv",
        frozen / "source-archive-audit.tsv",
        frozen / "read-selection-summary.tsv",
        frozen / "fastp-summary.tsv",
        frozen / "assembly-metrics.tsv",
        frozen / "read-recruitment.tsv",
        frozen / "resource-usage.tsv",
        frozen / "branch-contract.tsv",
        frozen / "tool-versions.tsv",
        frozen / "run-summary.json",
        chapter,
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("Missing Article 30 inputs: " + ", ".join(missing))

    checks = Checks()
    checksum_rows = verify_checksum_manifest(frozen, checks)
    audit_project_contract(root, frozen, checks)
    tool_rows = audit_tools(assembly_prefix, read_qc_prefix, checks)
    source_rows = audit_sources(frozen, checks)
    selection_rows, fastp_rows, output_pairs = audit_selection_and_fastp(frozen, checks)
    assembly_rows, _ = audit_assemblies(frozen, checks)
    recruitment_rows = audit_recruitment(frozen, output_pairs, checks)
    resource_rows = audit_resources(frozen, checks)
    audit_privacy_and_chapter(frozen, chapter, checks)

    frozen_versions = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    checks.add("Toolchain", "frozen-tool-versions", frozen_versions == EXPECTED_TOOLS, frozen_versions)
    branch_contract = read_tsv(frozen / "branch-contract.tsv")
    branch_contract_ok = (
        len(branch_contract) == 6
        and {row["Branch"] for row in branch_contract} == set(BRANCH_META)
        and all(int(row["Threads"]) == 16 for row in branch_contract)
        and all(int(row["ReadSelectionSeed"]) == SEED and int(row["Bowtie2Seed"]) == SEED for row in branch_contract)
        and all("64 GiB" in row["MemorySetting"] for row in branch_contract)
    )
    checks.add("Assembly", "branch-contract", branch_contract_ok, f"rows={len(branch_contract)}")
    run_summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    run_summary_ok = (
        run_summary["status"] == "completed"
        and run_summary["seed"] == SEED
        and run_summary["source_runs"] == 2
        and run_summary["selected_pairs_per_run"] == 2_000_000
        and run_summary["assembly_branches"] == 6
        and run_summary["mapping_branches"] == 8
        and run_summary["reporting_contig_threshold_bp"] == 1000
    )
    checks.add("Frozen input", "run-summary-completed", run_summary_ok, run_summary)

    configure_plotting()
    plot_branch_design(figure_dir)
    plot_contiguity(assembly_rows, figure_dir)
    plot_recruitment(recruitment_rows, figure_dir)
    plot_resources(resource_rows, figure_dir)
    audit_figures(figure_dir, checks)

    write_tsv(output_dir / "checksum-audit.tsv", checksum_rows, ["File", "ExpectedSHA256", "ObservedSHA256", "Status"])
    write_tsv(output_dir / "tool-audit.tsv", tool_rows, ["Tool", "ExpectedVersion", "ObservedVersion", "ReturnCode", "Status"])
    write_tsv(output_dir / "source-audit.tsv", source_rows, list(source_rows[0]))
    write_tsv(output_dir / "selection-audit.tsv", selection_rows, list(selection_rows[0]))
    write_tsv(output_dir / "fastp-audit.tsv", fastp_rows, list(fastp_rows[0]))
    write_tsv(output_dir / "assembly-audit.tsv", assembly_rows, list(assembly_rows[0]))
    write_tsv(output_dir / "recruitment-audit.tsv", recruitment_rows, list(recruitment_rows[0]))
    write_tsv(output_dir / "resource-audit.tsv", resource_rows, list(resource_rows[0]))
    write_tsv(output_dir / "validation-checks.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])

    assembly_resources = [row for row in resource_rows if row["StepType"] == "Assembly"]
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "seed": SEED,
        "source_project": "PRJEB52977",
        "source_runs": ["ERR9765746", "ERR9765747"],
        "selected_pairs_per_run": 2_000_000,
        "clean_pairs": output_pairs,
        "assembly_branches": 6,
        "mapping_branches": 8,
        "resource_rows": len(resource_rows),
        "assembly_total_bp_ge1000_range": [
            min(int(row["TotalBpGE1000"]) for row in assembly_rows),
            max(int(row["TotalBpGE1000"]) for row in assembly_rows),
        ],
        "assembly_n50_ge1000_range": [
            min(int(row["N50GE1000"]) for row in assembly_rows),
            max(int(row["N50GE1000"]) for row in assembly_rows),
        ],
        "recruitment_fraction_range": [
            min(float(row["MappedReadFraction"]) for row in recruitment_rows),
            max(float(row["MappedReadFraction"]) for row in recruitment_rows),
        ],
        "assembly_wall_minutes_range": [
            min(float(row["WallSeconds"]) for row in assembly_resources) / 60,
            max(float(row["WallSeconds"]) for row in assembly_resources) / 60,
        ],
        "assembly_peak_rss_gib_range": [
            min(int(row["PeakRSSKiB"]) for row in assembly_resources) / 1024**2,
            max(int(row["PeakRSSKiB"]) for row in assembly_resources) / 1024**2,
        ],
        "algorithm_only_benchmark": False,
        "universal_winner_claimed": False,
        "n50_treated_as_correctness": False,
        "reference_aware_accuracy_evaluated": False,
        "qa_network_access": False,
        "checksum_files": len(checksum_rows),
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "figures": [f"{stem}.{suffix}" for stem in FIGURE_STEMS for suffix in ("pdf", "png", "tiff")],
    }
    (output_dir / "validation-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    with (output_dir / "validation.log").open("w", encoding="utf-8") as handle:
        handle.write("Article 30 short-read assembly validation\n")
        handle.write(f"Status\t{summary['status']}\n")
        handle.write(f"ChecksPassed\t{checks.passed}\n")
        handle.write(f"ChecksFailed\t{checks.failed}\n")
        for row in checks.rows:
            if row["Status"] == "FAIL":
                handle.write(f"FAIL\t{row['Category']}\t{row['CheckID']}\t{row['Detail']}\n")
    print(json.dumps(summary, indent=2))
    return 0 if checks.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
