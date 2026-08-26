#!/usr/bin/env python3
"""Validate Article 32 frozen evidence and draw publication-ready figures."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import site
import statistics
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Iterator, TextIO

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article32-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image


SEED = 20260732
BRANCHES = (
    "spades-short-only",
    "spades-illumina-ont",
    "spades-illumina-hifi",
    "flye-ont",
    "flye-ont-polypolish-default",
    "flye-ont-polypolish-careful",
    "flye-hifi",
)
DISPLAY = {
    "spades-short-only": "Short-only",
    "spades-illumina-ont": "Illumina + ONT",
    "spades-illumina-hifi": "Illumina + HiFi",
    "flye-ont": "ONT draft",
    "flye-ont-polypolish-default": "ONT + Polypolish",
    "flye-ont-polypolish-careful": "ONT + careful",
    "flye-hifi": "HiFi-only",
}
COLORS = {
    "spades-short-only": "#7A8B99",
    "spades-illumina-ont": "#56B4E9",
    "spades-illumina-hifi": "#009E73",
    "flye-ont": "#D55E00",
    "flye-ont-polypolish-default": "#CC79A7",
    "flye-ont-polypolish-careful": "#8C6BB1",
    "flye-hifi": "#E69F00",
}
FIGURE_STEMS = (
    "32-hybrid-branch-design",
    "32-recovery-contiguity",
    "32-consensus-error",
    "32-abundance-resource-audit",
)
EXPECTED_SOURCES = {
    ("ERR9765746", "R1"): (
        1_740_647_656,
        "ed0e6e0ee846542531c742a45181cd6f",
        "8a26ae83abb5bb6ab8c8d30f9fcbfb43c2e04fd7b333eb11c6a40ef6e85c8c2f",
    ),
    ("ERR9765746", "R2"): (
        2_104_551_765,
        "5b60ac93cb69dff77ae38cfa501afd06",
        "7834df2021f96259061bdfcff9a7e7549950ec2509cd2e0730e5c07d17ecfd48",
    ),
    ("ERR9765780", "NA"): (
        3_117_261_341,
        "33eb90ac7437b0039180f03e7a697269",
        "fba6ae446bbd0436ed000d28059f982cde0c508aa6c4dd02f28730e40fcf2916",
    ),
    ("ERR9765783", "NA"): (
        3_982_506_052,
        "02ec4bc541b4e1ec5d0f58e4a519f2cb",
        "0edd53596da282cf947073cf8e4bdd6d980d6b9dd0811af2065d58dfa0624cd6",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    parser.add_argument("--chapter", type=Path)
    return parser.parse_args()


def sha256(path: Path, *, decompress: bool = False) -> str:
    digest = hashlib.sha256()
    opener = gzip.open if decompress else Path.open
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
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


def numeric(value: str | int | float | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text in {"", "-", "NA", "None"}:
        return None
    return float(text)


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


def nx(lengths: list[int], fraction: float = 0.5) -> tuple[int, int]:
    target = sum(lengths) * fraction
    cumulative = 0
    for index, length in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += length
        if cumulative >= target:
            return length, index
    return 0, 0


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


def verify_checksum_manifest(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    manifest = frozen / "file-checksums.sha256"
    rows: list[dict[str, str]] = []
    expected_names: set[str] = set()
    for number, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not line:
            continue
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if match is None:
            checks.add("Frozen input", f"checksum-line-{number}", False, line)
            continue
        expected, relative = match.groups()
        expected_names.add(relative)
        path = frozen / relative
        observed = sha256(path) if path.is_file() else "MISSING"
        checks.add("Frozen input", f"sha256-{relative}", observed == expected, observed)
        rows.append(
            {
                "File": relative,
                "ExpectedSHA256": expected,
                "ObservedSHA256": observed,
                "Status": "PASS" if observed == expected else "FAIL",
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


def command_output(command: list[str], env_prefix: Path) -> tuple[int, str]:
    process = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env={
            **os.environ,
            "PATH": f"{env_prefix}/bin:/usr/bin:/bin",
            "LC_ALL": "C",
            "TZ": "UTC",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "PYTHONPATH": "",
            "MPLCONFIGDIR": os.environ["MPLCONFIGDIR"],
        },
    )
    return process.returncode, process.stdout.strip()


def audit_project(
    root: Path, env_prefix: Path, frozen: Path, checks: Checks
) -> None:
    pairs = {
        "source-manifest": (root / "data/small/32-source-manifest.tsv", frozen / "source-manifest.tsv"),
        "reference-manifest": (root / "data/small/32-reference-manifest.tsv", frozen / "reference-manifest.tsv"),
        "software-releases": (root / "data/small/32-software-releases.tsv", frozen / "software-releases.tsv"),
        "branch-contract": (root / "data/small/32-branch-contract.tsv", frozen / "branch-contract.tsv"),
        "data-notice": (root / "data/small/32-data-NOTICE.txt", frozen / "data-NOTICE.txt"),
        "environment-yaml": (root / "env/hybrid-assembly.yml", frozen / "env/hybrid-assembly.yml"),
        "environment-lock": (root / "env/hybrid-assembly-linux-64.lock", frozen / "env/hybrid-assembly-linux-64.lock"),
        "long-read-environment-lock": (
            root / "env/long-read-assembly-linux-64.lock",
            frozen / "env/long-read-assembly-linux-64.lock",
        ),
    }
    for label, (project, copy) in pairs.items():
        passed = project.is_file() and copy.is_file() and sha256(project) == sha256(copy)
        checks.add("Project contract", f"frozen-{label}-current", passed, sha256(project) if project.is_file() else "MISSING")
    lock = (frozen / "env/hybrid-assembly-linux-64.lock").read_text(encoding="utf-8").splitlines()
    packages = [line for line in lock if line.startswith("https://")]
    required = (
        "fastp-1.3.6",
        "spades-4.3.0",
        "polypolish-0.6.1",
        "bwa-0.7.19",
        "quast-5.3.0",
        "python-3.10.20",
        "matplotlib-base-3.10.9",
        "pillow-12.3.0",
    )
    checks.add("Project contract", "explicit-environment-lock", "@EXPLICIT" in lock and len(packages) >= 190, len(packages))
    checks.add("Project contract", "direct-tools-in-lock", all(any(token in line.lower() for line in packages) for token in required), required)
    long_read_lock = (
        frozen / "env/long-read-assembly-linux-64.lock"
    ).read_text(encoding="utf-8").splitlines()
    checks.add(
        "Project contract",
        "flye-baseline-explicit-lock",
        "@EXPLICIT" in long_read_lock
        and any("/flye-2.9.6-" in line for line in long_read_lock),
        "Flye 2.9.6",
    )
    branch_contract = read_tsv(frozen / "branch-contract.tsv")
    branch_map = {row["Branch"]: row for row in branch_contract}
    branch_contract_ok = (
        set(branch_map) == set(BRANCHES)
        and len(branch_contract) == len(BRANCHES)
        and all(int(row["EvaluationMinimumBp"]) == 1000 for row in branch_contract)
        and all(
            row["SelectionSeed"] == "20260732"
            for row in branch_contract
            if row["ShortReadInput"] != "NA"
        )
        and branch_map["spades-illumina-ont"]["GraphSource"]
        == "Illumina de Bruijn graph"
        and branch_map["flye-ont-polypolish-default"]["GraphSource"]
        == "ONT repeat graph"
    )
    checks.add(
        "Project contract",
        "seven-branch-direction-contract",
        branch_contract_ok,
        sorted(branch_map),
    )
    portable_files = sorted((frozen / "logs").glob("*.log")) + sorted(
        (frozen / "resources").glob("*.txt")
    )
    leaked_paths = []
    for path in portable_files:
        text = path.read_text(encoding="utf-8", errors="replace")
        if str(root) in text or str(env_prefix) in text:
            leaked_paths.append(path.relative_to(frozen).as_posix())
    placeholder_text = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in portable_files
    )
    checks.add(
        "Project contract",
        "portable-logs-no-workspace-path-leak",
        not leaked_paths
        and "${PROJECT_ROOT}" in placeholder_text
        and "${HYBRID_ENV_PREFIX}" in placeholder_text,
        leaked_paths,
    )
    reused_long_read_assemblies = {
        "flye-ont": (
            frozen / "assemblies/flye-ont.ge1000.fna.gz",
            root / "data/small/31-long-read-assembly-frozen/assemblies/flye-ont-r9.ge1000.fna.gz",
        ),
        "flye-hifi": (
            frozen / "assemblies/flye-hifi.ge1000.fna.gz",
            root / "data/small/31-long-read-assembly-frozen/assemblies/flye-hifi.ge1000.fna.gz",
        ),
    }
    for branch, (article32_path, frozen_baseline_path) in reused_long_read_assemblies.items():
        if not article32_path.is_file() or not frozen_baseline_path.is_file():
            checks.add(
                "Project contract",
                f"frozen-baseline-sequence-identity-{branch}",
                False,
                "MISSING",
            )
            continue
        current, current_records, current_bases = sequence_multiset(article32_path)
        prior, prior_records, prior_bases = sequence_multiset(frozen_baseline_path)
        checks.add(
            "Project contract",
            f"frozen-baseline-sequence-identity-{branch}",
            current == prior
            and current_records == prior_records
            and current_bases == prior_bases,
            f"records={current_records}/{prior_records} bases={current_bases}/{prior_bases}",
        )


def audit_tools(env_prefix: Path, frozen: Path, checks: Checks) -> list[dict[str, str]]:
    quast_minimap2_header = (
        env_prefix
        / "lib/python3.10/site-packages/quast_libs/minimap2/minimap.h"
    )
    probes = {
        "fastp": ([str(env_prefix / "bin/fastp"), "--version"], r"fastp 1\.3\.6", True),
        "SPAdes": ([str(env_prefix / "bin/spades.py"), "--version"], r"SPAdes genome assembler v4\.3\.0", True),
        "BWA-MEM": ([str(env_prefix / "bin/bwa")], r"Version: 0\.7\.19", False),
        "Polypolish": ([str(env_prefix / "bin/polypolish"), "--version"], r"Polypolish 0\.6\.1", True),
        "MetaQUAST": ([str(env_prefix / "bin/metaquast.py"), "--version"], r"QUAST v5\.3\.0", True),
        "Python": ([str(env_prefix / "bin/python"), "--version"], r"Python 3\.10\.20", True),
    }
    rows: list[dict[str, str]] = []
    for tool, (command, pattern, require_zero) in probes.items():
        code, output = command_output(command, env_prefix)
        passed = (code == 0 or not require_zero) and re.search(pattern, output) is not None
        checks.add("Toolchain", f"runtime-{tool}", passed, output.replace("\n", ";")[:500])
        rows.append(
            {
                "Tool": tool,
                "ExpectedPattern": pattern,
                "ReturnCode": str(code),
                "Observed": output.replace("\n", ";")[:700],
                "Status": "PASS" if passed else "FAIL",
            }
        )
    bundled_source = (
        quast_minimap2_header.read_text(encoding="utf-8", errors="replace")
        if quast_minimap2_header.is_file()
        else ""
    )
    bundled_ok = '#define MM_VERSION "2.28-r1209"' in bundled_source
    checks.add(
        "Toolchain",
        "metaquast-bundled-minimap2-source-version",
        bundled_ok,
        "2.28-r1209" if bundled_ok else "MISSING",
    )
    rows.append(
        {
            "Tool": "MetaQUAST bundled minimap2 source",
            "ExpectedPattern": 'MM_VERSION "2.28-r1209"',
            "ReturnCode": "NA",
            "Observed": "2.28-r1209" if bundled_ok else "MISSING",
            "Status": "PASS" if bundled_ok else "FAIL",
        }
    )
    isolation_ok = site.ENABLE_USER_SITE is False and os.environ.get("PYTHONNOUSERSITE") == "1"
    checks.add("Toolchain", "python-user-site-disabled", isolation_ok, f"ENABLE_USER_SITE={site.ENABLE_USER_SITE}")
    checks.add("Toolchain", "numpy-version", np.__version__ == "2.2.6", np.__version__)
    checks.add("Toolchain", "matplotlib-version", matplotlib.__version__ == "3.10.9", matplotlib.__version__)
    import PIL

    checks.add("Toolchain", "pillow-version", PIL.__version__ == "12.3.0", PIL.__version__)
    return rows


def audit_sources(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    manifest = read_tsv(frozen / "source-manifest.tsv")
    audit = read_tsv(frozen / "source-audit.tsv")
    manifest_map = {(row["RunAccession"], row["Mate"]): row for row in manifest}
    audit_map = {(row["RunAccession"], row["Mate"]): row for row in audit}
    checks.add("Sources", "source-key-set", set(manifest_map) == set(EXPECTED_SOURCES) == set(audit_map), sorted(manifest_map))
    rows = []
    for key, (size, md5, digest) in EXPECTED_SOURCES.items():
        m = manifest_map.get(key, {})
        a = audit_map.get(key, {})
        passed = (
            int(m.get("ENABytes", -1)) == size
            and m.get("ENAReportedMD5") == md5
            and m.get("ObservedSHA256") == digest
            and int(a.get("ObservedBytes", -1)) == size
            and a.get("ObservedMD5") == md5
            and a.get("ObservedSHA256") == digest
            and a.get("IdentityStatus") == "PASS"
        )
        checks.add("Sources", f"identity-{key[0]}-{key[1]}", passed, a)
        rows.append({"RunAccession": key[0], "Mate": key[1], "Status": "PASS" if passed else "FAIL", "SHA256": a.get("ObservedSHA256", "")})
    return rows


def audit_selection_fastp(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    selection = json.loads((frozen / "selection-summary.json").read_text(encoding="utf-8"))
    expected = {
        "seed": SEED,
        "rng_namespace": "20260732:MOCK1:ERR9765746",
        "source_pairs": 20_597_525,
        "selected_pairs": 10_000_000,
        "selected_pair_id_sha256": "c6578b14f0e643406a5d92e4c6f02a653e758f42ce8394bbccbaec6fc41f2e75",
    }
    for field, value in expected.items():
        checks.add("Selection/QC", f"selection-{field}", selection.get(field) == value, selection.get(field))
    mate_expected = {
        "R1": (1_497_641_905, 850_551_064, "c43ce1ced64c5b8d869cafd858903a6079841627dc8384e3fe5057917ab20ba8"),
        "R2": (1_481_707_773, 1_027_453_379, "7b9b3b6558391fe0beec805fae66cc929bcf6c64463856389fe3af0a9758caa4"),
    }
    rows = []
    for mate, (bases, size, digest) in mate_expected.items():
        observed = selection["mates"][mate]
        passed = observed["bases"] == bases and observed["bytes"] == size and observed["compressed_sha256"] == digest
        checks.add("Selection/QC", f"selection-{mate}-identity", passed, observed)
        rows.append({"Stage": "selection", "Mate": mate, "Reads": selection["selected_pairs"], "Bases": bases, "Status": "PASS" if passed else "FAIL"})
    fastp = json.loads((frozen / "fastp.json").read_text(encoding="utf-8"))
    before = fastp["summary"]["before_filtering"]
    after = fastp["summary"]["after_filtering"]
    expected_fastp = {
        "before": (20_000_000, 2_979_349_678, 0.900057),
        "after": (19_998_640, 2_979_266_463, 0.900059),
    }
    checks.add(
        "Selection/QC",
        "fastp-input-identity",
        (before["total_reads"], before["total_bases"], before["q30_rate"])
        == expected_fastp["before"],
        before,
    )
    checks.add(
        "Selection/QC",
        "fastp-output-identity",
        (after["total_reads"], after["total_bases"], after["q30_rate"])
        == expected_fastp["after"],
        after,
    )
    checks.add("Selection/QC", "fastp-paired-output", after["total_reads"] % 2 == 0 and after["total_reads"] <= before["total_reads"], after["total_reads"])
    checks.add("Selection/QC", "fastp-retains-ge90pct", after["total_reads"] >= 18_000_000, after["total_reads"])
    rows.append({"Stage": "fastp-before", "Mate": "paired", "Reads": before["total_reads"], "Bases": before["total_bases"], "Status": "PASS"})
    rows.append({"Stage": "fastp-after", "Mate": "paired", "Reads": after["total_reads"], "Bases": after["total_bases"], "Status": "PASS"})
    clean = {row["Mate"]: row for row in read_tsv(frozen / "clean-fastq-audit.tsv")}
    checks.add("Selection/QC", "clean-mate-set", set(clean) == {"R1", "R2"}, sorted(clean))
    clean_records = {int(row["Records"]) for row in clean.values()}
    clean_hashes = {row["PairIDHash"] for row in clean.values()}
    clean_bases = sum(int(row["Bases"]) for row in clean.values())
    clean_expected = {
        "R1": (
            843_284_587,
            "691e8190af5206a1b748212f8f0f785b2c9bf516b48920a4265250550baf0a52",
        ),
        "R2": (
            1_013_512_874,
            "f02ad9333d405e0882ecba10dde870facc91fdb43cb1ba778b71244154d4a8ae",
        ),
    }
    clean_ok = (
        clean_records == {after["total_reads"] // 2}
        and len(clean_hashes) == 1
        and clean_bases == after["total_bases"]
        and all(int(row["MinimumLength"]) >= 50 for row in clean.values())
        and all(int(row["CompressedBytes"]) > 0 for row in clean.values())
        and all(re.fullmatch(r"[0-9a-f]{64}", row["CompressedSHA256"]) for row in clean.values())
        and all(
            int(clean[mate]["CompressedBytes"]) == expected_size
            and clean[mate]["CompressedSHA256"] == expected_digest
            for mate, (expected_size, expected_digest) in clean_expected.items()
        )
    )
    checks.add(
        "Selection/QC",
        "clean-fastq-grammar-count-pair-identity",
        clean_ok,
        f"records={sorted(clean_records)} bases={clean_bases} pair_hashes={len(clean_hashes)}",
    )
    for mate in ("R1", "R2"):
        row = clean[mate]
        rows.append(
            {
                "Stage": "fastp-clean",
                "Mate": mate,
                "Reads": row["Records"],
                "Bases": row["Bases"],
                "Status": "PASS" if clean_ok else "FAIL",
            }
        )
    return rows


def sequence_multiset(path: Path) -> tuple[Counter[str], int, int]:
    digests: Counter[str] = Counter()
    records = bases = 0
    for _, sequence in fasta_records(path):
        if re.search(r"[^ACGTRYSWKMBDHVN]", sequence):
            raise ValueError(f"Invalid nucleotide alphabet in {path}")
        digests[hashlib.sha256(sequence.encode()).hexdigest()] += 1
        records += 1
        bases += len(sequence)
    return digests, records, bases


def audit_truth(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    audit = json.loads((frozen / "truth-audit.json").read_text(encoding="utf-8"))
    expected = {
        "benchmark_commit": "a429a3724d4593f35b8d7323b20252a6be90e1cd",
        "supplement_s1_sha256": "937653a56fea7fbfcbe35b3f35c721b4125072ba4ab04c44c9d454697240c6df",
        "supplement_s2_sha256": "74004b9dd5e43b4b7a80e9a92bb17f76412c09d658f8ef77f76486645167f0f3",
        "combined_reference_sha256": "ad2641ac155006387e722ae5ec8592fa077f9ab3cf411cc11f5757430b8e752f",
        "mock1_genomes": 71,
        "mock1_reference_records": 310,
        "mock1_reference_bases": 239_414_479,
        "taxonomy_rename_aliases": 9,
        "sequence_multiset_equal": True,
    }
    for field, value in expected.items():
        checks.add("Truth", f"truth-{field}", audit.get(field) == value, audit.get(field))
    direct_hashes = {
        "supplement-s1": (
            frozen / "sources/Supplementary_Table_S1.xlsx",
            expected["supplement_s1_sha256"],
        ),
        "supplement-s2": (
            frozen / "sources/Supplementary_Table_S2.xlsx",
            expected["supplement_s2_sha256"],
        ),
        "combined-reference": (
            frozen / "references/MOCK_001.fasta.gz",
            expected["combined_reference_sha256"],
        ),
    }
    for label, (path, expected_digest) in direct_hashes.items():
        observed = sha256(path) if path.is_file() else "MISSING"
        checks.add("Truth", f"direct-sha256-{label}", observed == expected_digest, observed)
    truth = read_tsv(frozen / "truth-manifest.tsv")
    checks.add("Truth", "truth-manifest-71", len(truth) == 71, len(truth))
    checks.add("Truth", "truth-abundance-rounding", abs(sum(float(row["ExpectedAbundancePct"]) for row in truth) - 100) < 0.1, sum(float(row["ExpectedAbundancePct"]) for row in truth))
    combined_counter, combined_records, combined_bases = sequence_multiset(frozen / "references/MOCK_001.fasta.gz")
    individual_counter: Counter[str] = Counter()
    rows = []
    individual_records = individual_bases = 0
    for row in truth:
        path = frozen / "references/genomes" / f"{row['CurrentGenomeLabel']}.fna.gz"
        digest_ok = path.is_file() and sha256(path) == row["CompressedSHA256"]
        counter, records, bases = sequence_multiset(path)
        individual_counter.update(counter)
        individual_records += records
        individual_bases += bases
        passed = digest_ok and records == int(row["ReferenceRecords"]) and bases == int(row["ReferenceBases"])
        checks.add("Truth genome", f"reference-{row['CurrentGenomeLabel']}", passed, f"records={records} bases={bases}")
        rows.append({"Reference": row["CurrentGenomeLabel"], "Records": str(records), "Bases": str(bases), "Status": "PASS" if passed else "FAIL"})
    equality = individual_counter == combined_counter and individual_records == combined_records == 310 and individual_bases == combined_bases == 239_414_479
    checks.add("Truth", "separate-equals-combined-sequence-multiset", equality, f"records={combined_records} bases={combined_bases}")
    return rows


def audit_assemblies(frozen: Path, checks: Checks) -> tuple[list[dict[str, Any]], dict[str, dict[str, str]]]:
    metrics = {row["Branch"]: row for row in read_tsv(frozen / "branch-metrics.tsv")}
    frozen_rows = {row["Branch"]: row for row in read_tsv(frozen / "frozen-assemblies.tsv")}
    checks.add("Assemblies", "seven-branch-set", set(metrics) == set(frozen_rows) == set(BRANCHES), sorted(metrics))
    audits: list[dict[str, Any]] = []
    for branch in BRANCHES:
        path = frozen / "assemblies" / f"{branch}.ge1000.fna.gz"
        records = list(fasta_records(path))
        ids = [name.split()[0] for name, _ in records]
        lengths = [len(sequence) for _, sequence in records]
        n50, l50 = nx(lengths)
        n_bases = sum(sequence.count("N") for _, sequence in records)
        decompressed_hash = sha256(path, decompress=True)
        passed = (
            len(ids) == len(set(ids))
            and min(lengths) >= 1000
            and len(records) == int(float(metrics[branch]["SequencesGe1kb"])) == int(frozen_rows[branch]["Contigs"])
            and sum(lengths) == int(float(metrics[branch]["TotalLengthBp"])) == int(frozen_rows[branch]["Bases"])
            and max(lengths) == int(float(metrics[branch]["LargestBp"]))
            and n50 == int(float(metrics[branch]["N50Bp"]))
            and l50 == int(float(metrics[branch]["L50"]))
            and n_bases == int(float(metrics[branch]["NBases"]))
            and sha256(path) == frozen_rows[branch]["CompressedSHA256"]
            and decompressed_hash == metrics[branch]["AssemblySHA256"]
        )
        checks.add("Assemblies", f"assembly-{branch}", passed, f"contigs={len(records)} bases={sum(lengths)} n50={n50}")
        audits.append({"Branch": branch, "Contigs": len(records), "Bases": sum(lengths), "LargestBp": max(lengths), "N50Bp": n50, "L50": l50, "Status": "PASS" if passed else "FAIL"})
    return audits, metrics


def audit_polishing(
    frozen: Path,
    metrics: dict[str, dict[str, str]],
    checks: Checks,
) -> list[dict[str, Any]]:
    rows = read_tsv(frozen / "polishing-sequence-audit.tsv")
    by_mode = {row["Mode"]: row for row in rows}
    checks.add(
        "Polypolish",
        "two-prespecified-modes",
        set(by_mode) == {"default", "careful"}
        and len(rows) == 2
        and int(
            json.loads(
                (frozen / "computation-summary.json").read_text(encoding="utf-8")
            )["polishing_sequence_audit_modes"]
        )
        == 2,
        sorted(by_mode),
    )
    log_path = frozen / "polypolish-log-audit.tsv"
    log_rows = read_tsv(log_path) if log_path.is_file() else []
    logs = {row["Mode"]: row for row in log_rows}
    checks.add(
        "Polypolish",
        "two-nonempty-runtime-logs",
        set(logs) == {"default", "careful"}
        and all(int(row["LogBytes"]) > 0 for row in log_rows),
        sorted(logs),
    )
    bwa_r1 = (frozen / "resources/bwa-align-R1.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    bwa_r2 = (frozen / "resources/bwa-align-R2.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    bwa_contract = (
        " mem -t 32 -a " in bwa_r1
        and " mem -t 32 -a " in bwa_r2
        and "clean10m_R1.fastq.gz" in bwa_r1
        and "clean10m_R2.fastq.gz" not in bwa_r1
        and "clean10m_R2.fastq.gz" in bwa_r2
        and "clean10m_R1.fastq.gz" not in bwa_r2
        and "flye-ont.ge1000.fasta" in bwa_r1
        and "flye-ont.ge1000.fasta" in bwa_r2
    )
    checks.add(
        "Polypolish",
        "bwa-all-alignments-mates-separate",
        bwa_contract,
        "BWA-MEM -a R1/R2",
    )
    filter_resource = (frozen / "resources/polypolish-filter.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    default_resource = (frozen / "resources/polypolish-default.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    careful_resource = (frozen / "resources/polypolish-careful.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    mode_contract = (
        " filter --in1 " in filter_resource
        and " --in2 " in filter_resource
        and " --out1 " in filter_resource
        and " --out2 " in filter_resource
        and " polish " in default_resource
        and "--careful" not in default_resource
        and " polish --careful " in careful_resource
        and all(
            token in default_resource and token in careful_resource
            for token in ("filtered_R1.sam", "filtered_R2.sam")
        )
    )
    checks.add(
        "Polypolish",
        "filter-default-careful-command-contract",
        mode_contract,
        "filter + default + careful",
    )
    audits: list[dict[str, Any]] = []
    draft_metric = metrics["flye-ont"]
    for mode in ("default", "careful"):
        row = by_mode[mode]
        log_row = logs.get(mode, {})
        branch = f"flye-ont-polypolish-{mode}"
        output_metric = metrics[branch]
        shared = int(row["CanonicalSharedSequenceIDs"])
        identical = int(row["IdenticalSequences"])
        changed = int(row["ChangedSequences"])
        passed = (
            row["Branch"] == branch
            and row["CanonicalSequenceIDSetEqual"] == "yes"
            and row["CanonicalSequenceOrderEqual"] == "yes"
            and int(row["LostCanonicalSequenceIDs"]) == 0
            and int(row["NewCanonicalSequenceIDs"]) == 0
            and int(row["ExpectedPolypolishHeaderAnnotations"])
            == int(row["OutputSequences"])
            and int(row["UnexpectedOutputHeaders"]) == 0
            and int(row["DraftSequences"]) == int(float(draft_metric["SequencesGe1kb"]))
            and int(row["OutputSequences"]) == int(float(output_metric["SequencesGe1kb"]))
            and shared == int(row["DraftSequences"]) == int(row["OutputSequences"])
            and identical + changed == shared
            and int(row["DraftTotalBp"]) == int(float(draft_metric["TotalLengthBp"]))
            and int(row["OutputTotalBp"]) == int(float(output_metric["TotalLengthBp"]))
            and int(row["LengthDeltaBp"])
            == int(row["OutputTotalBp"]) - int(row["DraftTotalBp"])
            and row["DraftAssemblySHA256"] == draft_metric["AssemblySHA256"]
            and row["OutputAssemblySHA256"] == output_metric["AssemblySHA256"]
        )
        checks.add(
            "Polypolish",
            f"consensus-only-sequence-contract-{mode}",
            passed,
            f"canonical_shared={shared} changed={changed} length_delta={row['LengthDeltaBp']}",
        )
        log_passed = (
            int(log_row.get("AlignmentFiles", -1)) == 2
            and int(log_row.get("InputAlignments", -1))
            == int(log_row.get("HighQualityAlignmentsKept", -2))
            + int(log_row.get("HighQualityAlignmentsDiscarded", -3))
            and int(log_row.get("PolishedSequences", -1))
            == int(row["DraftSequences"])
            and int(log_row.get("InputBases", -1))
            == int(row["DraftTotalBp"])
            and int(log_row.get("OutputBases", -1))
            == int(row["OutputTotalBp"])
            and int(log_row.get("ChangedPositions", -1)) > 0
            and log_row.get("CarefulFlagVisible")
            == ("yes" if mode == "careful" else "no")
        )
        checks.add(
            "Polypolish",
            f"runtime-log-metrics-{mode}",
            log_passed,
            (
                f"kept={log_row.get('HighQualityAlignmentsKept', 'NA')} "
                f"changed_positions={log_row.get('ChangedPositions', 'NA')}"
            ),
        )
        audits.append(
            {
                **row,
                "LogBytes": log_row.get("LogBytes", "NA"),
                "HighQualityAlignmentsKept": log_row.get(
                    "HighQualityAlignmentsKept", "NA"
                ),
                "HighQualityAlignmentsDiscarded": log_row.get(
                    "HighQualityAlignmentsDiscarded", "NA"
                ),
                "ChangedPositions": log_row.get("ChangedPositions", "NA"),
                "Status": "PASS" if passed and log_passed else "FAIL",
            }
        )
    return audits


def audit_metaquast(frozen: Path, metrics: dict[str, dict[str, str]], checks: Checks) -> list[dict[str, Any]]:
    per_genome = read_tsv(frozen / "per-genome-metaquast.tsv")
    computation = json.loads(
        (frozen / "computation-summary.json").read_text(encoding="utf-8")
    )
    checks.add("MetaQUAST", "per-genome-497-rows", len(per_genome) == 71 * 7, len(per_genome))
    counts = Counter(row["Branch"] for row in per_genome)
    checks.add("MetaQUAST", "each-branch-71-genomes", set(counts) == set(BRANCHES) and set(counts.values()) == {71}, counts)
    reports = list((frozen / "metaquast/references").glob("*.tsv"))
    frozen_summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    checks.add(
        "MetaQUAST",
        "physical-reference-report-inventory",
        1 <= len(reports) <= 71
        and len(reports) == int(frozen_summary["per_reference_reports"]),
        len(reports),
    )
    missing_branch_reports = sum(
        row["BranchReportPresent"] == "no" for row in per_genome
    )
    checks.add(
        "MetaQUAST",
        "explicit-zero-fill-ledger",
        int(computation["physical_reference_reports"]) == len(reports)
        and int(computation["missing_reference_reports"]) == 71 - len(reports)
        and int(computation["missing_branch_reports_zero_recovery"])
        == missing_branch_reports,
        f"reports={len(reports)} missing_branch_rows={missing_branch_reports}",
    )
    missing_rows = [
        row for row in per_genome if row["BranchReportPresent"] == "no"
    ]
    checks.add(
        "MetaQUAST",
        "missing-alignment-zero-recovery-na-error-rates",
        all(
            float(row["GenomeFractionPct"]) == 0
            and numeric(row["MismatchesPer100Kbp"]) is None
            and numeric(row["IndelsPer100Kbp"]) is None
            for row in missing_rows
        ),
        len(missing_rows),
    )
    truth = read_tsv(frozen / "truth-manifest.tsv")
    truth_labels = {row["CurrentGenomeLabel"] for row in truth}
    observed_labels = {row["Reference"] for row in per_genome}
    checks.add("MetaQUAST", "truth-labels-preserved", observed_labels == truth_labels, len(observed_labels))
    audits = []
    for branch in BRANCHES:
        rows = [row for row in per_genome if row["Branch"] == branch]
        recovered90 = sum(row["RecoveredGe90Pct"] == "yes" for row in rows)
        full99 = sum(row["FullGenomeGe99Pct"] == "yes" for row in rows)
        passed = recovered90 == int(float(metrics[branch]["RecoveredGenomesGe90Pct"])) and full99 == int(float(metrics[branch]["FullGenomesGe99Pct"]))
        checks.add("MetaQUAST", f"recovery-counts-{branch}", passed, f">=90={recovered90} >=99={full99}")
        audits.append({"Branch": branch, "RecoveredGe90Pct": recovered90, "FullGenomesGe99Pct": full99, "GenomeFractionPct": metrics[branch]["GenomeFractionPct"], "MismatchesPer100Kbp": metrics[branch]["MismatchesPer100Kbp"], "IndelsPer100Kbp": metrics[branch]["IndelsPer100Kbp"], "Status": "PASS" if passed else "FAIL"})
    bins = read_tsv(frozen / "abundance-bin-recovery.tsv")
    bin_map = {(row["Branch"], row["AbundanceBin"]): row for row in bins}
    for branch in BRANCHES:
        for abundance_bin in ("<0.1%", "0.1-<1%", ">=1%"):
            rows = [row for row in per_genome if row["Branch"] == branch and row["AbundanceBin"] == abundance_bin]
            expected_median = statistics.median(float(row["GenomeFractionPct"]) for row in rows)
            observed = bin_map[(branch, abundance_bin)]
            passed = len(rows) == int(observed["TruthGenomes"]) and math.isclose(expected_median, float(observed["MedianGenomeFractionPct"]), abs_tol=5e-7)
            checks.add("MetaQUAST bins", f"median-{branch}-{abundance_bin}", passed, f"median={expected_median:.6f}")
    return audits


def audit_split_scaffolds(frozen: Path, checks: Checks) -> list[dict[str, Any]]:
    path = frozen / "split-scaffold-sensitivity.tsv"
    rows = read_tsv(path) if path.is_file() else []
    expected_labels = {f"{branch}_broken" for branch in BRANCHES}
    observed_labels = {row["Assembly"] for row in rows}
    parents = [row["ParentBranch"] for row in rows]
    inventory_ok = (
        1 <= len(rows) <= len(BRANCHES)
        and len(observed_labels) == len(rows)
        and observed_labels <= expected_labels
        and len(set(parents)) == len(parents)
        and set(parents) <= set(BRANCHES)
    )
    computation = json.loads(
        (frozen / "computation-summary.json").read_text(encoding="utf-8")
    )
    inventory_ok = inventory_ok and int(
        computation["split_scaffold_sensitivity_branches"]
    ) == len(rows)
    checks.add(
        "MetaQUAST split scaffolds",
        "broken-branch-sensitivity-inventory",
        inventory_ok,
        sorted(observed_labels),
    )
    audited: list[dict[str, Any]] = []
    for row in rows:
        numeric_fields = (
            "SequencesGe1kb",
            "TotalLengthBp",
            "N50Bp",
            "Misassemblies",
            "GenomeFractionPct",
            "MismatchesPer100Kbp",
            "IndelsPer100Kbp",
        )
        passed = all(numeric(row[field]) is not None for field in numeric_fields)
        checks.add(
            "MetaQUAST split scaffolds",
            f"broken-branch-metrics-{row['ParentBranch']}",
            passed,
            row["Assembly"],
        )
        audited.append({**row, "Status": "PASS" if passed else "FAIL"})
    return audited


def audit_resources(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    rows = read_tsv(frozen / "resource-usage.tsv")
    required_prefixes = (
        "fastp-illumina-10m",
        "assemble-spades-short-only.attempt",
        "assemble-spades-illumina-ont.attempt",
        "assemble-spades-illumina-hifi.attempt",
        "bwa-index-flye-ont",
        "bwa-align-R1",
        "bwa-align-R2",
        "polypolish-filter",
        "polypolish-default",
        "polypolish-careful",
        "metaquast-seven-branches",
    )
    successful_steps = {
        prefix
        for prefix in required_prefixes
        if any(
            row["Step"].startswith(prefix)
            and numeric(row["ExitStatus"]) == 0
            and numeric(row["ElapsedSeconds"]) is not None
            and numeric(row["PeakRSSGiB"]) is not None
            for row in rows
        )
    }
    for prefix in required_prefixes:
        checks.add(
            "Resources",
            f"successful-resource-record-{prefix}",
            prefix in successful_steps,
            prefix,
        )
    for row in rows:
        status = "PASS" if numeric(row["ExitStatus"]) == 0 else "HISTORY"
        row["Status"] = status
    return rows


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / f"{stem}.png", dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(figure_dir / f"{stem}.tiff", dpi=350, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def style_axis(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#DCE3E7", linewidth=0.7, zorder=0)
    ax.tick_params(colors="#263238", labelsize=8)


def draw_design(figure_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 6.3))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7.4)
    ax.axis("off")

    def box(x: float, y: float, width: float, height: float, text: str, color: str, fontsize: float = 9) -> None:
        patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.04,rounding_size=0.08", facecolor=color, edgecolor="#37474F", linewidth=1.0)
        ax.add_patch(patch)
        ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=fontsize, color="#172126")

    def arrow(x1: float, y1: float, x2: float, y2: float, rad: float = 0) -> None:
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12, linewidth=1.2, color="#607D8B", connectionstyle=f"arc3,rad={rad}"))

    def panel(x: float, y: float, width: float, height: float, color: str) -> None:
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.05,rounding_size=0.12",
            facecolor=color,
            edgecolor="#90A4AE",
            linewidth=1.0,
            zorder=0,
        )
        ax.add_patch(patch)

    ax.text(0.4, 7.02, "Same MOCK1 DNA; strategy audit, not a causal platform trial", fontsize=13.5, weight="bold", color="#263238")
    ax.text(0.4, 6.62, "Graph source separates hybrid assembly from post-assembly consensus polishing.", fontsize=9.8, color="#455A64")
    ax.text(0.4, 6.22, "Locked inputs", fontsize=8.5, weight="bold", color="#607D8B")
    box(0.4, 5.35, 3.35, 0.68, "10 M Illumina pairs · seed 20260732", "#E9EEF1", fontsize=8.7)
    box(4.32, 5.35, 3.35, 0.68, "Historical ONT R9 · 696,944 reads", "#D9EFF9", fontsize=8.7)
    box(8.25, 5.35, 3.35, 0.68, "PacBio HiFi · 524,805 reads", "#DDF3EA", fontsize=8.7)

    panel(0.4, 1.38, 5.35, 3.62, "#F4F9FC")
    panel(6.15, 1.38, 5.45, 3.62, "#FCF7F2")
    ax.text(0.7, 4.66, "Short-read-first graph", fontsize=11.5, weight="bold", color="#263238")
    ax.text(0.7, 4.34, "SPAdes 4.3.0 · Illumina de Bruijn graph defines topology", fontsize=8.2, color="#546E7A")
    box(0.72, 2.62, 1.42, 1.03, "Short-only\n1 branch", "#E9EEF1", fontsize=8.5)
    box(2.36, 2.62, 1.42, 1.03, "Illumina + ONT\n1 branch", "#D9EFF9", fontsize=8.5)
    box(4.00, 2.62, 1.42, 1.03, "Illumina + HiFi\n1 branch", "#DDF3EA", fontsize=8.5)
    ax.text(0.72, 1.91, "Long reads bridge graph paths; they do not define the initial graph.", fontsize=8.0, color="#455A64")

    ax.text(6.45, 4.66, "Long-read-first graph", fontsize=11.5, weight="bold", color="#263238")
    ax.text(6.45, 4.34, "Flye 2.9.6 repeat graph defines topology", fontsize=8.2, color="#546E7A")
    box(6.48, 3.02, 1.55, 0.88, "ONT draft\n1 branch", "#FBE8D7", fontsize=8.4)
    box(9.66, 3.02, 1.55, 0.88, "HiFi-only\n1 branch", "#FBE8D7", fontsize=8.4)
    box(6.25, 1.78, 1.45, 0.76, "Polypolish\ndefault", "#F1E4F2", fontsize=8.1)
    box(7.93, 1.78, 1.45, 0.76, "Polypolish\ncareful", "#F1E4F2", fontsize=8.1)
    arrow(7.25, 3.02, 6.98, 2.54, rad=0.12)
    arrow(7.45, 3.02, 8.65, 2.54, rad=-0.12)
    ax.text(10.48, 2.15, "Short reads modify\nconsensus only", fontsize=7.6, color="#455A64", ha="center")

    box(3.47, 0.18, 5.06, 0.86, "MetaQUAST 5.3.0 · 7 primary assemblies\n71 references · identity >=97%", "#FFF2CC", fontsize=8.8)
    arrow(3.08, 1.38, 4.55, 1.04)
    arrow(8.88, 1.38, 7.45, 1.04)
    save_figure(fig, figure_dir, "32-hybrid-branch-design")


def draw_recovery(metrics: dict[str, dict[str, str]], figure_dir: Path) -> None:
    order = list(reversed(BRANCHES))
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 5.3), sharey=True)
    specifications = (
        ("N50Bp", "N50 (bp)", True),
        ("GenomeFractionPct", "Genome fraction (%)", False),
        ("FullGenomesGe99Pct", "Genomes with >=99% fraction", False),
    )
    for ax, (field, title, log_scale) in zip(axes, specifications):
        values = [float(metrics[branch][field]) for branch in order]
        y = np.arange(len(order))
        ax.scatter(values, y, s=58, c=[COLORS[branch] for branch in order], edgecolors="white", linewidths=0.7, zorder=3)
        for value, row_y in zip(values, y):
            text = f"{value:,.0f}" if field != "GenomeFractionPct" else f"{value:.1f}"
            ax.annotate(text, (value, row_y), xytext=(5, 0), textcoords="offset points", va="center", fontsize=7)
        if log_scale:
            ax.set_xscale("log")
            ax.xaxis.set_major_locator(mticker.LogLocator(base=10, numticks=6))
            ax.xaxis.set_minor_locator(mticker.NullLocator())
            ax.xaxis.set_major_formatter(
                mticker.FuncFormatter(
                    lambda value, _: (
                        f"{value / 1_000_000:g}M"
                        if value >= 1_000_000
                        else (f"{value / 1_000:g}k" if value >= 1_000 else f"{value:g}")
                    )
                )
            )
        ax.set_title(title, fontsize=10.5, weight="bold")
        ax.set_yticks(y, [DISPLAY[branch] for branch in order])
        style_axis(ax)
    fig.suptitle("Contiguity and recovery must be read together", fontsize=14, weight="bold", x=0.06, ha="left")
    fig.text(0.06, 0.92, "All assemblies were normalized to contigs >=1 kb before the same 71-reference audit.", fontsize=9, color="#455A64")
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    save_figure(fig, figure_dir, "32-recovery-contiguity")


def draw_errors(metrics: dict[str, dict[str, str]], figure_dir: Path) -> None:
    order = list(reversed(BRANCHES))
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 5.3), sharey=True)
    specifications = (
        ("MismatchesPer100Kbp", "Mismatches / 100 kbp"),
        ("IndelsPer100Kbp", "Indels / 100 kbp"),
        ("Misassemblies", "Misassemblies"),
    )
    for ax, (field, title) in zip(axes, specifications):
        values = [float(metrics[branch][field]) for branch in order]
        y = np.arange(len(order))
        ax.barh(y, values, color=[COLORS[branch] for branch in order], alpha=0.88, zorder=2)
        for value, row_y in zip(values, y):
            ax.annotate(f"{value:,.2f}" if field != "Misassemblies" else f"{value:,.0f}", (value, row_y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=7)
        ax.set_xscale("symlog", linthresh=1)
        ax.set_title(title, fontsize=10.5, weight="bold")
        ax.set_yticks(y, [DISPLAY[branch] for branch in order])
        style_axis(ax)
    fig.suptitle("A longer or more recovered assembly can still carry more error", fontsize=14, weight="bold", x=0.06, ha="left")
    fig.text(0.06, 0.92, "Log-like axes retain zero while showing orders-of-magnitude differences.", fontsize=9, color="#455A64")
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    save_figure(fig, figure_dir, "32-consensus-error")


def draw_abundance_resources(frozen: Path, figure_dir: Path) -> None:
    bins = read_tsv(frozen / "abundance-bin-recovery.tsv")
    bin_order = ("<0.1%", "0.1-<1%", ">=1%")
    values = np.array([[float(next(row["MedianGenomeFractionPct"] for row in bins if row["Branch"] == branch and row["AbundanceBin"] == abundance_bin)) for abundance_bin in bin_order] for branch in BRANCHES])
    resources = [
        row
        for row in read_tsv(frozen / "resource-usage.tsv")
        if numeric(row["ElapsedSeconds"]) is not None
        and numeric(row["ExitStatus"]) == 0
    ]
    resources = sorted(resources, key=lambda row: float(row["ElapsedSeconds"]), reverse=True)[:8]
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.8), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    image = ax.imshow(values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            color = "white" if values[row, column] < 55 else "#172126"
            ax.text(column, row, f"{values[row, column]:.1f}", ha="center", va="center", fontsize=8, color=color)
    ax.set_xticks(np.arange(3), bin_order)
    ax.set_yticks(np.arange(7), [DISPLAY[branch] for branch in BRANCHES])
    ax.set_xlabel("Expected abundance bin")
    ax.set_title("Median genome fraction (%)", fontsize=11, weight="bold")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.03)
    colorbar.ax.tick_params(labelsize=7)

    ax = axes[1]
    step_labels = {
        "spades-short-only": "SPAdes short-only",
        "spades-illumina-ont": "SPAdes Illumina + ONT",
        "spades-illumina-hifi": "SPAdes Illumina + HiFi",
        "bwa-index-flye-ont": "BWA index ONT draft",
        "bwa-align-R1": "BWA all-alignments R1",
        "bwa-align-R2": "BWA all-alignments R2",
        "polypolish-filter": "Polypolish pair filter",
        "polypolish-default": "Polypolish default",
        "polypolish-careful": "Polypolish careful",
        "fastp-illumina-10m": "fastp 10 M pairs",
        "select-illumina-10m": "Exact 10 M-pair selection",
        "metaquast-seven-branches": "MetaQUAST 7 branches",
    }
    labels = []
    for row in reversed(resources):
        step = row["Step"].removeprefix("assemble-")
        step = re.sub(r"\.attempt[0-9]+$", "", step)
        labels.append(step_labels.get(step, step))
    hours = [float(row["ElapsedSeconds"]) / 3600 for row in reversed(resources)]
    rss = [float(row["PeakRSSGiB"]) if numeric(row["PeakRSSGiB"]) is not None else 0 for row in reversed(resources)]
    y = np.arange(len(labels))
    ax.barh(y, hours, color="#607D8B", alpha=0.85)
    for value, memory, row_y in zip(hours, rss, y):
        ax.annotate(f"{value:.2f} h · {memory:.1f} GiB", (value, row_y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=7)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Wall time (hours)")
    ax.set_title("Largest recorded workflow steps", fontsize=11, weight="bold")
    style_axis(ax)
    fig.suptitle("Low-abundance recovery and compute cost are part of the decision", fontsize=14, weight="bold", x=0.05, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save_figure(fig, figure_dir, "32-abundance-resource-audit")


def audit_figures(figure_dir: Path, checks: Checks) -> None:
    for stem in FIGURE_STEMS:
        for suffix in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            checks.add("Figures", f"figure-{stem}-{suffix}", path.is_file() and path.stat().st_size > 1000, path.stat().st_size if path.is_file() else "MISSING")
        for suffix in ("png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            try:
                with Image.open(path) as image:
                    width, height = image.size
                    dpi = image.info.get("dpi", (0, 0))
                    compression = image.info.get("compression", "NA")
                    image.verify()
                dpi_x, dpi_y = float(dpi[0]), float(dpi[1])
                passed = (
                    width >= 1600
                    and height >= 900
                    and min(dpi_x, dpi_y) >= 300
                    and (suffix != "tiff" or compression == "tiff_lzw")
                )
                detail = (
                    f"{width}x{height} dpi={dpi_x:.1f}/{dpi_y:.1f} "
                    f"compression={compression}"
                )
            except Exception as error:  # pragma: no cover - diagnostic path
                passed, detail = False, repr(error)
            checks.add("Figures", f"raster-integrity-{stem}-{suffix}", passed, detail)


def audit_chapter(chapter: Path, checks: Checks) -> None:
    text = chapter.read_text(encoding="utf-8")
    checks.add("Chapter", "frontmatter-draft-false", re.search(r"^draft:\s*false\s*$", text, re.MULTILINE) is not None, chapter)
    checks.add("Chapter", "frontmatter-eval-false", re.search(r"^\s*eval:\s*false\s*$", text, re.MULTILINE) is not None, chapter)
    section_tokens = (
        "## 这一步对应论文里的哪张图",
        "## 理论",
        "## 准备工作",
        "## 可复制代码",
        "## 审计与升级",
        "## 出版级美化",
        "## 常见坑",
        "## 这段 Methods 怎么写",
        "## 换成你自己的数据怎么做",
        "## 参考",
    )
    for index, token in enumerate(section_tokens, 1):
        checks.add("Chapter", f"required-section-{index:02d}", token in text, token)
    checks.add("Chapter", "inline-theme-no-source", "theme_pub <- function" in text and 'source("R/theme_pub.R")' not in text, "inline theme")
    checks.add("Chapter", "seed-visible", "20260732" in text, "seed 20260732")
    checks.add(
        "Chapter",
        "actual-results-inserted",
        "ARTICLE32_ACTUAL_TABLE" not in text,
        "actual-results marker must be removed",
    )
    for stem in FIGURE_STEMS:
        checks.add("Chapter", f"figure-reference-{stem}", f"../figures/{stem}.png" in text, stem)
    prohibited = (
        "Planned chapter",
        "Do not publish",
        "本篇可独立跑通",
        "这体现全系列",
        "作者代码通常长这样",
    )
    checks.add("Chapter", "no-reader-facing-meta-writing", all(token not in text for token in prohibited), prohibited)
    method_tokens = ("SPAdes 4.3.0", "Polypolish 0.6.1", "BWA-MEM v0.7.19", "MetaQUAST 5.3.0", "minimap2 2.28-r1209", "a429a3724d4593f35b8d7323b20252a6be90e1cd")
    checks.add("Chapter", "methods-version-contract", all(token in text for token in method_tokens), method_tokens)
    boundary_tokens = (
        "不是平台单因素实验",
        "不能写“Polypolish 保证提升准确性”",
        "操作性回收阈值",
    )
    checks.add(
        "Chapter",
        "interpretation-boundaries",
        all(token in text for token in boundary_tokens),
        boundary_tokens,
    )


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    env_prefix = args.env_prefix.resolve()
    frozen = (args.frozen_dir or root / "data/small/32-hybrid-assembly-polishing-frozen").resolve()
    output = (args.output_dir or root / "results/32-hybrid-assembly-polishing").resolve()
    figure_dir = (args.figure_dir or root / "figures").resolve()
    chapter = (args.chapter or root / "chapters/32-hybrid-assembly-polishing.qmd").resolve()
    output.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    checksum_rows = verify_checksum_manifest(frozen, checks)
    audit_project(root, env_prefix, frozen, checks)
    tool_rows = audit_tools(env_prefix, frozen, checks)
    source_rows = audit_sources(frozen, checks)
    selection_rows = audit_selection_fastp(frozen, checks)
    truth_rows = audit_truth(frozen, checks)
    assembly_rows, metrics = audit_assemblies(frozen, checks)
    polishing_rows = audit_polishing(frozen, metrics, checks)
    metaquast_rows = audit_metaquast(frozen, metrics, checks)
    split_scaffold_rows = audit_split_scaffolds(frozen, checks)
    resource_rows = audit_resources(frozen, checks)

    draw_design(figure_dir)
    draw_recovery(metrics, figure_dir)
    draw_errors(metrics, figure_dir)
    draw_abundance_resources(frozen, figure_dir)
    audit_figures(figure_dir, checks)
    audit_chapter(chapter, checks)

    write_tsv(output / "checksum-audit.tsv", checksum_rows, ["File", "ExpectedSHA256", "ObservedSHA256", "Status"])
    write_tsv(output / "tool-audit.tsv", tool_rows, ["Tool", "ExpectedPattern", "ReturnCode", "Observed", "Status"])
    write_tsv(output / "source-audit.tsv", source_rows, ["RunAccession", "Mate", "Status", "SHA256"])
    write_tsv(output / "selection-fastp-audit.tsv", selection_rows, ["Stage", "Mate", "Reads", "Bases", "Status"])
    write_tsv(output / "truth-reference-audit.tsv", truth_rows, ["Reference", "Records", "Bases", "Status"])
    write_tsv(output / "assembly-audit.tsv", assembly_rows, ["Branch", "Contigs", "Bases", "LargestBp", "N50Bp", "L50", "Status"])
    write_tsv(
        output / "polishing-audit.tsv",
        polishing_rows,
        [
            "Mode",
            "Branch",
            "DraftSequences",
            "OutputSequences",
            "CanonicalSharedSequenceIDs",
            "LostCanonicalSequenceIDs",
            "NewCanonicalSequenceIDs",
            "CanonicalSequenceIDSetEqual",
            "CanonicalSequenceOrderEqual",
            "ExpectedPolypolishHeaderAnnotations",
            "UnexpectedOutputHeaders",
            "IdenticalSequences",
            "ChangedSequences",
            "DraftTotalBp",
            "OutputTotalBp",
            "LengthDeltaBp",
            "DraftAssemblySHA256",
            "OutputAssemblySHA256",
            "LogBytes",
            "HighQualityAlignmentsKept",
            "HighQualityAlignmentsDiscarded",
            "ChangedPositions",
            "Status",
        ],
    )
    write_tsv(output / "metaquast-audit.tsv", metaquast_rows, ["Branch", "RecoveredGe90Pct", "FullGenomesGe99Pct", "GenomeFractionPct", "MismatchesPer100Kbp", "IndelsPer100Kbp", "Status"])
    write_tsv(
        output / "split-scaffold-audit.tsv",
        split_scaffold_rows,
        ["Assembly", "ParentBranch", "SequencesGe1kb", "TotalLengthBp", "N50Bp", "Misassemblies", "GenomeFractionPct", "MismatchesPer100Kbp", "IndelsPer100Kbp", "Status"],
    )
    write_tsv(
        output / "resource-audit.tsv",
        resource_rows,
        ["Step", "ElapsedSeconds", "PeakRSSGiB", "CPUPercent", "FileSystemInputs", "FileSystemOutputs", "ExitStatus", "Status"],
    )
    write_tsv(output / "validation-checks.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "article": 32,
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "checks_total": len(checks.rows),
        "assembly_branches": len(BRANCHES),
        "truth_genomes": 71,
        "seed": SEED,
        "platform_only_experiment": False,
        "universal_hybrid_winner_claimed": False,
        "polishing_guaranteed_to_improve": False,
        "full_genome_threshold_is_proof": False,
        "figures": list(FIGURE_STEMS),
    }
    (output / "validation-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "validation.log").write_text(json.dumps(summary, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    return 0 if checks.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
