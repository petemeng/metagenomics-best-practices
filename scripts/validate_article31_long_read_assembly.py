#!/usr/bin/env python3
"""Validate Article 31 frozen evidence and draw publication-ready figures."""

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
from typing import Any, Iterable, Iterator, TextIO

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article31-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image


SEED = 20260731
EXPECTED_HASHES = {
    "source_manifest": "4124c98a449e3fc7e0e43e7e934fda42912b7960260e8ed2b23e5cdf7a9399f4",
    "data_notice": "eb4240b83227886b013dd3069243e7b47d02a191f381cff9dff57884c83cd026",
    "software_releases": "e5319c2c77c82ddc72353c546efae06d9f3f0b350e1430882a5e9c245d88d19d",
    "environment_yaml": "8ca21481514ea0eaf8e16eead0b45d7f1a9f1c958fa9b1c165efc18b1e4ed563",
    "environment_lock": "3260dea112e0edc8bc56deec85cc2ca2b5e958537f6a4f798362709baa68f7c8",
}
EXPECTED_TOOLS = {
    "Flye": "2.9.6-b1802",
    "metaMDBG": "1.4",
    "minimap2": "2.31-r1302",
    "samtools": "1.23.1",
    "SeqKit": "2.11.0",
    "Python": "3.10.20",
    "hifiasm-meta release": "0.3.5-r81",
    "hifiasm-meta embedded": "ha base version: 0.13-r308;hamt version: 0.3-r079",
}
EXPECTED_SOURCES = {
    "ERR9765780": {
        "platform": "ONT",
        "sample": "SAMEA14435832",
        "bytes": 3_117_261_341,
        "md5": "33eb90ac7437b0039180f03e7a697269",
        "reads": 696_944,
        "bases": 3_125_920_499,
    },
    "ERR9765783": {
        "platform": "HiFi",
        "sample": "SAMEA14435832",
        "bytes": 3_982_506_052,
        "md5": "02ec4bc541b4e1ec5d0f58e4a519f2cb",
        "reads": 524_805,
        "bases": 5_400_038_744,
    },
}
BRANCHES = (
    ("flye-ont-r9", "Flye", "ONT R9"),
    ("flye-hifi", "Flye", "PacBio HiFi"),
    ("hifiasm-meta-hifi", "hifiasm-meta", "PacBio HiFi"),
    ("metamdbg-hifi", "metaMDBG", "PacBio HiFi"),
)
BRANCH_META = {branch: (assembler, platform) for branch, assembler, platform in BRANCHES}
FIGURE_STEMS = (
    "31-long-read-branch-design",
    "31-long-read-contiguity",
    "31-long-read-readback",
    "31-circular-resource-audit",
)
JUNCTION_AUDIT_FIELDS = [
    "Branch",
    "Assembler",
    "Platform",
    "JunctionID",
    "ContigID",
    "LengthBp",
    "JunctionEligible",
    "JunctionFlankBp",
    "SequenceSHA256",
    "OriginalHeader",
    "IdentityFloor",
    "IdentityComputation",
    "AuditedPrimaryPAFRecords",
    "BaseLevelCigarPAFRecords",
    "JunctionSpanningReads",
    "SupportedByGE3Reads",
]
PALETTE = {"Flye": "#0072B2", "hifiasm-meta": "#D55E00", "metaMDBG": "#009E73"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
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


def verify_checksum_manifest(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    manifest = frozen / "file-checksums.sha256"
    rows: list[dict[str, str]] = []
    expected_names: set[str] = set()
    for number, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), start=1):
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
        passed = observed == expected
        checks.add("Frozen input", f"sha256-{relative}", passed, observed)
        rows.append(
            {
                "File": relative,
                "ExpectedSHA256": expected,
                "ObservedSHA256": observed,
                "Status": "PASS" if passed else "FAIL",
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
        },
    )
    return process.returncode, process.stdout.strip()


def audit_project_contract(root: Path, frozen: Path, checks: Checks) -> None:
    paths = {
        "source_manifest": root / "data/small/31-source-manifest.tsv",
        "data_notice": root / "data/small/31-data-NOTICE.txt",
        "software_releases": root / "data/small/31-software-releases.tsv",
        "environment_yaml": root / "env/long-read-assembly.yml",
        "environment_lock": root / "env/long-read-assembly-linux-64.lock",
    }
    frozen_paths = {
        "source_manifest": frozen / "source-manifest.tsv",
        "data_notice": frozen / "data-NOTICE.txt",
        "software_releases": frozen / "software-releases.tsv",
        "environment_yaml": frozen / "env/long-read-assembly.yml",
        "environment_lock": frozen / "env/long-read-assembly-linux-64.lock",
    }
    for label, path in paths.items():
        observed = sha256(path)
        frozen_observed = sha256(frozen_paths[label])
        expected = EXPECTED_HASHES.get(label)
        checks.add(
            "Project contract",
            f"{label}-sha256",
            expected is not None
            and observed == expected
            and frozen_observed == expected,
            f"project={observed}; frozen={frozen_observed}",
        )
    lock_lines = (paths["environment_lock"]).read_text(encoding="utf-8").splitlines()
    package_rows = [line for line in lock_lines if line.startswith("https://")]
    checks.add(
        "Project contract",
        "explicit-lock-104-packages",
        "@EXPLICIT" in lock_lines and len(package_rows) == 104,
        f"packages={len(package_rows)}",
    )
    required_packages = ("flye-2.9.6", "metamdbg-1.4", "minimap2-2.31", "samtools-1.23.1", "seqkit-2.11.0")
    checks.add(
        "Project contract",
        "direct-tools-in-lock",
        all(any(token in line.lower() for line in package_rows) for token in required_packages),
        required_packages,
    )
    release_rows = read_tsv(frozen / "software-releases.tsv")
    release_map = {row["Tool"]: row for row in release_rows}
    release_ok = (
        set(release_map) == {"Flye", "hifiasm-meta", "metaMDBG"}
        and release_map["Flye"]["TagCommit"] == "886b8c17412cdf3a2868a28237bca6c5ad1da156"
        and release_map["hifiasm-meta"]["TagCommit"] == "e4e24f5158091ad901c1ff6f68278559bd41a6b5"
        and release_map["metaMDBG"]["TagCommit"] == "22d7040cdfb384c9898192f0aaa266a05263c5eb"
    )
    checks.add("Project contract", "official-release-commits", release_ok, sorted(release_map))


def audit_tools(env_prefix: Path, frozen: Path, checks: Checks) -> list[dict[str, str]]:
    probes = {
        "Flye": ([str(env_prefix / "bin/flye"), "--version"], r"2\.9\.6-b1802"),
        "metaMDBG": ([str(env_prefix / "bin/metaMDBG")], r"Version:\s+1\.4"),
        "minimap2": ([str(env_prefix / "bin/minimap2"), "--version"], r"2\.31-r1302"),
        "samtools": ([str(env_prefix / "bin/samtools"), "--version"], r"samtools 1\.23\.1"),
        "SeqKit": ([str(env_prefix / "bin/seqkit"), "version"], r"seqkit v2\.11\.0"),
        "Python": ([str(env_prefix / "bin/python"), "--version"], r"Python 3\.10\.20"),
    }
    rows: list[dict[str, str]] = []
    for tool, (command, pattern) in probes.items():
        code, output = command_output(command, env_prefix)
        passed = code == 0 and re.search(pattern, output) is not None
        checks.add("Toolchain", f"runtime-{tool}", passed, output.splitlines()[0] if output else "")
        rows.append(
            {
                "Tool": tool,
                "ExpectedPattern": pattern,
                "ReturnCode": str(code),
                "Observed": output.replace("\n", ";")[:600],
                "Status": "PASS" if passed else "FAIL",
            }
        )
    frozen_versions = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    checks.add("Toolchain", "frozen-version-table", frozen_versions == EXPECTED_TOOLS, frozen_versions)
    source = read_tsv(frozen / "hifiasm-meta-source.tsv")
    source_ok = (
        len(source) == 1
        and source[0]["Tag"] == "hamtv0.3.5"
        and source[0]["Commit"] == "e4e24f5158091ad901c1ff6f68278559bd41a6b5"
        and source[0]["SourceSHA256"] == "8c1c1f394e0d4d3be2c78cb76c4122dd0caf2d088a8b986c161d4d90c194f560"
        and source[0]["VersionOutput"] == EXPECTED_TOOLS["hifiasm-meta embedded"]
        and source[0]["TSNESeedOption"] == "--tsne-seed"
        and source[0]["TSNESeedDefault"] == "42"
        and re.fullmatch(r"[0-9a-f]{64}", source[0]["BinarySHA256"]) is not None
    )
    checks.add("Toolchain", "hifiasm-source-build-ledger", source_ok, source[0] if source else "missing")
    return rows


def audit_sources(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    manifest = read_tsv(frozen / "source-manifest.tsv")
    audit = read_tsv(frozen / "source-audit.tsv")
    manifest_by_run = {row["RunAccession"]: row for row in manifest}
    audit_by_run = {row["RunAccession"]: row for row in audit}
    checks.add(
        "Sources",
        "source-run-set",
        set(manifest_by_run) == set(EXPECTED_SOURCES) == set(audit_by_run),
        sorted(manifest_by_run),
    )
    for run, expected in EXPECTED_SOURCES.items():
        row = audit_by_run.get(run, {})
        passed = bool(row) and (
            row["PlatformKey"] == expected["platform"]
            and row["SampleAccession"] == expected["sample"]
            and int(row["ObservedBytes"]) == expected["bytes"]
            and row["ObservedMD5"] == expected["md5"]
            and int(row["ObservedReadCount"]) == expected["reads"]
            and int(row["ObservedBaseCount"]) == expected["bases"]
            and row["IdentityStatus"] == "PASS"
            and re.fullmatch(r"[0-9a-f]{64}", row["ObservedSHA256"]) is not None
        )
        checks.add("Sources", f"archive-{run}", passed, row.get("ObservedSHA256", "missing"))
    same_sample = len({row["SampleAccession"] for row in audit}) == 1
    unequal_budget = len({int(row["ObservedBaseCount"]) for row in audit}) == 2
    checks.add("Sources", "same-dna-source", same_sample, "SAMEA14435832")
    checks.add("Sources", "unequal-platform-budgets-retained", unequal_budget, [row["ObservedBaseCount"] for row in audit])
    return audit


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
                    raise ValueError(f"Sequence before header: {path}")
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


def metrics_from_inventory(rows: list[dict[str, str]], threshold: int) -> dict[str, int]:
    selected = [row for row in rows if int(row["LengthBp"]) >= threshold]
    lengths = [int(row["LengthBp"]) for row in selected]
    n50, l50 = nx(lengths, 0.5)
    n90, l90 = nx(lengths, 0.9)
    return {
        "ContigCount": len(lengths),
        "TotalBp": sum(lengths),
        "N50Bp": n50,
        "L50": l50,
        "N90Bp": n90,
        "L90": l90,
        "LargestBp": max(lengths, default=0),
        "CircularCandidates": sum(row["CircularCandidate"] == "TRUE" for row in selected),
    }


def audit_assemblies(frozen: Path, checks: Checks) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    inventory = read_tsv(frozen / "contig-inventory.tsv")
    metrics = read_tsv(frozen / "assembly-metrics.tsv")
    expected_branches = set(BRANCH_META)
    checks.add("Assemblies", "inventory-branch-set", {row["Branch"] for row in inventory} == expected_branches, Counter(row["Branch"] for row in inventory))
    checks.add("Assemblies", "metric-grid", len(metrics) == 20 and {int(row["ThresholdBp"]) for row in metrics} == {0, 1000, 10000, 100000, 1000000}, f"rows={len(metrics)}")
    for branch, (assembler, platform) in BRANCH_META.items():
        branch_inventory = [row for row in inventory if row["Branch"] == branch]
        metadata_ok = all(row["Assembler"] == assembler and row["Platform"] == platform for row in branch_inventory)
        checks.add("Assemblies", f"metadata-{branch}", metadata_ok and bool(branch_inventory), f"contigs={len(branch_inventory)}")
        for threshold in (0, 1000, 10000, 100000, 1000000):
            rows = [row for row in metrics if row["Branch"] == branch and int(row["ThresholdBp"]) == threshold]
            expected = metrics_from_inventory(branch_inventory, threshold)
            passed = len(rows) == 1 and all(int(rows[0][field]) == value for field, value in expected.items())
            checks.add("Assemblies", f"metric-ledger-{branch}-{threshold}", passed, expected)
        frozen_path = frozen / "assemblies" / f"{branch}.ge1000.fna.gz"
        observed = list(fasta_records(frozen_path))
        expected_inventory = {row["ContigID"]: row for row in branch_inventory if int(row["LengthBp"]) >= 1000}
        ids = [header.split(None, 1)[0] for header, _ in observed]
        payload_ok = len(ids) == len(set(ids)) == len(expected_inventory) and set(ids) == set(expected_inventory)
        if payload_ok:
            for header, sequence in observed:
                contig_id = header.split(None, 1)[0]
                row = expected_inventory[contig_id]
                payload_ok &= len(sequence) == int(row["LengthBp"])
                payload_ok &= hashlib.sha256(sequence.encode()).hexdigest() == row["SequenceSHA256"]
        checks.add("Assemblies", f"frozen-fasta-{branch}", payload_ok, f"records={len(observed)}")
    return metrics, inventory


def close_float(left: float, right: float, tolerance: float = 1e-8) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def audit_readback(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    rows = read_tsv(frozen / "readback-metrics.tsv")
    checks.add("Read-back", "branch-set", {row["Branch"] for row in rows} == set(BRANCH_META) and len(rows) == 4, [row["Branch"] for row in rows])
    for row in rows:
        branch = row["Branch"]
        assembler, platform = BRANCH_META[branch]
        expected_run = EXPECTED_SOURCES["ERR9765780" if platform == "ONT R9" else "ERR9765783"]
        expected_reads = int(row["ExpectedReads"])
        expected_bases = int(row["ExpectedBases"])
        frozen_fasta = frozen / "assemblies" / f"{branch}.ge1000.fna.gz"
        reference_digest = hashlib.sha256()
        reference_bytes = 0
        with gzip.open(frozen_fasta, "rb") as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                reference_digest.update(block)
                reference_bytes += len(block)
        mapped = int(row["MappedReads"])
        mapq20 = int(row["MapQ20Reads"])
        aligned = int(row["AlignedQueryBasesUnion"])
        exact_matches = int(row["ExactMatchingBases"])
        alignment_blocks = int(row["AlignmentBlockBases"])
        identity = float(row["WeightedAlignmentIdentityPct"])
        identity_floor = 70.0 if platform == "ONT R9" else 95.0
        passed = (
            row["Assembler"] == assembler
            and row["Platform"] == platform
            and expected_reads == expected_run["reads"]
            and expected_bases == expected_run["bases"]
            and int(row["ReferenceThresholdBp"]) == 1_000
            and int(row["ReferenceBytes"]) == reference_bytes
            and row["ReferenceSHA256"] == reference_digest.hexdigest()
            and 0 <= mapq20 <= mapped <= expected_reads
            and mapped > 0
            and 0 < aligned <= expected_bases
            and close_float(float(row["MappedReadPct"]), 100 * mapped / expected_reads)
            and close_float(float(row["MapQ20ReadPct"]), 100 * mapq20 / expected_reads)
            and close_float(float(row["AlignedQueryBasePct"]), 100 * aligned / expected_bases)
            and row["IdentityComputation"] == "minimap2-c-paf-matches-over-block"
            and int(row["BaseLevelCigarPAFRecords"]) == int(row["PrimaryPAFRecords"])
            and 0 < exact_matches <= alignment_blocks
            and close_float(identity, 100 * exact_matches / alignment_blocks)
            and identity_floor <= identity <= 100
            and 0 <= float(row["MultiContigPctOfMapped"]) <= 100
            and 0 <= float(row["EndClippedGE1kbPctOfMapped"]) <= 100
        )
        checks.add("Read-back", f"ledger-{branch}", passed, f"mapped={mapped}/{expected_reads}")
    return rows


def audit_junctions(frozen: Path, inventory: list[dict[str, str]], checks: Checks) -> list[dict[str, str]]:
    candidates = read_tsv(frozen / "circular-candidates.tsv")
    support = read_tsv(frozen / "junction-support.tsv")
    circular_inventory = {(row["Branch"], row["ContigID"], row["SequenceSHA256"]) for row in inventory if row["CircularCandidate"] == "TRUE"}
    candidate_set = {(row["Branch"], row["ContigID"], row["SequenceSHA256"]) for row in candidates}
    checks.add("Circular audit", "candidate-inventory-join", candidate_set == circular_inventory, f"candidates={len(candidate_set)}")
    support_ids = {row["JunctionID"] for row in support}
    candidate_ids = {row["JunctionID"] for row in candidates}
    checks.add("Circular audit", "support-complete-join", support_ids == candidate_ids and len(support) == len(candidates), f"rows={len(support)}")
    for row in support:
        count = int(row["JunctionSpanningReads"])
        audited = int(row["AuditedPrimaryPAFRecords"])
        base_level = int(row["BaseLevelCigarPAFRecords"])
        eligible = row["JunctionEligible"] == "TRUE"
        passed = (
            0 <= count <= audited
            and base_level == audited
            and row["IdentityComputation"] == "minimap2-c-paf-matches-over-block"
            and (row["SupportedByGE3Reads"] == "TRUE") == (count >= 3)
            and (int(row["JunctionFlankBp"]) == 5000 if eligible else int(row["JunctionFlankBp"]) == 0)
        )
        checks.add("Circular audit", f"junction-{row['JunctionID']}", passed, f"support={count}; records={audited}")
    return support


def audit_resources(frozen: Path, support: list[dict[str, str]], checks: Checks) -> list[dict[str, str]]:
    rows = read_tsv(frozen / "resource-usage.tsv")
    labels = {row["Step"] for row in rows}
    required = {f"assemble-{branch}" for branch in BRANCH_META} | {f"map-{branch}" for branch in BRANCH_META}
    eligible_branches = {row["Branch"] for row in support if row["JunctionEligible"] == "TRUE"}
    required |= {f"junction-{branch}" for branch in eligible_branches}
    checks.add("Resources", "required-step-coverage", required <= labels, sorted(labels))
    valid = all(
        int(row["Attempts"]) >= 1
        and float(row["WallSeconds"]) > 0
        and int(row["PeakRSSKiB"]) > 0
        and float(row["FinalWallSeconds"]) > 0
        and int(row["FinalPeakRSSKiB"]) > 0
        and int(row["FinalExitStatus"]) == 0
        for row in rows
    )
    checks.add("Resources", "resource-values", valid, f"rows={len(rows)}")
    return rows


def audit_boundaries_and_chapter(frozen: Path, chapter: Path, checks: Checks) -> None:
    forbidden_suffixes = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".sam", ".bam", ".cram", ".paf")
    forbidden = [path.relative_to(frozen).as_posix() for path in frozen.rglob("*") if path.is_file() and path.name.endswith(forbidden_suffixes)]
    checks.add("Boundaries", "no-read-or-alignment-payload", not forbidden, forbidden)
    leaks: list[str] = []
    for path in frozen.rglob("*"):
        if not path.is_file() or path.suffix.lower() in {".gz", ".png", ".pdf", ".tiff"}:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        local_tokens = (
            "/" + "media/" + "desk16/",
            "/" + "home/" + "tly" + "9658/",
            "tly" + "9658",
        )
        if any(token in text for token in local_tokens):
            leaks.append(path.relative_to(frozen).as_posix())
    checks.add("Boundaries", "no-local-path-leaks", not leaks, leaks)

    text = chapter.read_text(encoding="utf-8")
    tokens = (
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
    checks.add("Chapter", "freeze-auto", re.search(r"^freeze:\s*auto\s*$", text, re.M) is not None, chapter.name)
    checks.add("Chapter", "nine-section-contract", all(token in text for token in tokens), tokens)
    checks.add("Chapter", "inline-plot-functions", all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")), "pal/theme/save")
    checks.add("Chapter", "no-source-theme-dependency", 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text, "inline only")
    checks.add("Chapter", "hardware-four-dimensions", all(token in text for token in ("RAM", "磁盘", "核", "耗时")), "RAM/disk/cores/time")
    checks.add("Chapter", "no-server-alternative", any(token in text for token in ("HPC", "云", "子集演示")), "HPC/cloud/subset")
    checks.add("Chapter", "seed-visible", str(SEED) in text, SEED)
    checks.add("Chapter", "tool-versions-visible", all(token in text for token in ("Flye 2.9.6", "hifiasm-meta 0.3.5", "metaMDBG 1.4", "minimap2 2.31")), "versions")
    checks.add(
        "Chapter",
        "common-readback-reference",
        all(
            token in text
            for token in (
                "--reference-threshold-bp 1000",
                "ReferenceSHA256",
                "minimap2 -c",
                "≥1 kb",
            )
        ),
        "same >=1 kb reference contract",
    )
    checks.add("Chapter", "figure-references", all(f"figures/{stem}.png" in text for stem in FIGURE_STEMS), FIGURE_STEMS)
    checks.add("Chapter", "circular-boundary", all(token in text for token in ("circular candidate", "不等于完整基因组", "junction", "read-back")), "candidate is not proof")
    checks.add("Chapter", "comparison-boundary", all(token in text for token in ("不是平台单因素", "不代表普适赢家", "N50", "正确性", "第 33 篇")), "non-causal boundary")
    meta_phrases = ("本篇可独立跑通", "全系列", "接口只学一次", "作者代码通常长这样", "（即本文）")
    checks.add("Chapter", "no-author-meta-copy", not any(phrase in text for phrase in meta_phrases), meta_phrases)


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
    pdf = figure_dir / f"{stem}.pdf"
    png = figure_dir / f"{stem}.png"
    tiff = figure_dir / f"{stem}.tiff"
    fig.savefig(pdf, bbox_inches="tight", metadata={"Creator": "metagenomics-best-practices", "CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=350, bbox_inches="tight", facecolor="white", metadata={"Software": "metagenomics-best-practices"})
    with Image.open(png) as image:
        image.convert("RGB").save(tiff, compression="tiff_lzw", dpi=(350, 350))
    plt.close(fig)


def add_box(axis: plt.Axes, x: float, y: float, width: float, height: float, text: str, face: str, edge: str = "#334E5C") -> None:
    axis.add_patch(FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.02,rounding_size=0.025", linewidth=1.1, facecolor=face, edgecolor=edge))
    axis.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=9.0)


def add_arrow(axis: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    axis.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=11, linewidth=1, color="#607D8B"))


def plot_design(figure_dir: Path) -> None:
    fig, axis = plt.subplots(figsize=(10.5, 4.5))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    add_box(axis, 0.03, 0.62, 0.20, 0.16, "ONT MinION R9\nERR9765780\n3.13 Gbp", "#E6F2F8")
    add_box(axis, 0.03, 0.22, 0.20, 0.16, "PacBio HiFi\nERR9765783\n5.40 Gbp", "#F0EBF7")
    add_box(axis, 0.36, 0.64, 0.20, 0.12, "Flye 2.9.6\n--nano-raw", "#D9ECF6", PALETTE["Flye"])
    add_box(axis, 0.36, 0.42, 0.20, 0.12, "Flye 2.9.6\n--pacbio-hifi", "#D9ECF6", PALETTE["Flye"])
    add_box(axis, 0.36, 0.23, 0.20, 0.12, "hifiasm-meta 0.3.5\n--force-rs", "#F8E5DA", PALETTE["hifiasm-meta"])
    add_box(axis, 0.36, 0.04, 0.20, 0.12, "metaMDBG 1.4\n--in-hifi", "#DDF2E8", PALETTE["metaMDBG"])
    add_box(axis, 0.72, 0.36, 0.24, 0.25, "Common audit\n≥1 kb contigs\nread-back + split alarms\njunction support\ntime + peak RSS", "#F5F7F8")
    add_arrow(axis, (0.23, 0.70), (0.36, 0.70))
    for end in (0.48, 0.29, 0.10):
        add_arrow(axis, (0.23, 0.30), (0.36, end))
    for start, end in zip((0.70, 0.48, 0.29, 0.10), (0.55, 0.50, 0.45, 0.40)):
        add_arrow(axis, (0.56, start), (0.72, end))
    axis.text(0.5, 0.94, "Platform-adapted long-read assembly audit", ha="center", fontsize=13, fontweight="bold")
    axis.text(0.5, -0.02, "Same DNA source, unequal read budgets and error models; this is not a platform-only experiment.", ha="center", color="#455A64", fontsize=9)
    save_figure(fig, figure_dir, FIGURE_STEMS[0])


def ordered_rows(rows: list[dict[str, str]], threshold: int = 1000) -> list[dict[str, str]]:
    lookup = {(row["Branch"], int(row["ThresholdBp"])): row for row in rows}
    return [lookup[(branch, threshold)] for branch, _, _ in BRANCHES]


def short_label(branch: str) -> str:
    return {
        "flye-ont-r9": "ONT R9\nFlye",
        "flye-hifi": "HiFi\nFlye",
        "hifiasm-meta-hifi": "HiFi\nhifiasm-meta",
        "metamdbg-hifi": "HiFi\nmetaMDBG",
    }[branch]


def plot_contiguity(rows: list[dict[str, str]], figure_dir: Path) -> None:
    ordered = ordered_rows(rows)
    x = np.arange(4)
    colors = [PALETTE[row["Assembler"]] for row in ordered]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.2))
    panels = (("TotalBp", 1e6, "Assembled sequence ≥1 kb (Mbp)"), ("N50Bp", 1e6, "N50 at ≥1 kb (Mbp)"), ("LargestBp", 1e6, "Largest contig (Mbp)"))
    for axis, (field, divisor, ylabel) in zip(axes, panels):
        values = [int(row[field]) / divisor for row in ordered]
        bars = axis.bar(x, values, color=colors, width=0.72)
        axis.set_xticks(x, [short_label(row["Branch"]) for row in ordered], rotation=25, ha="right")
        axis.set_ylabel(ylabel)
        axis.grid(axis="y", color="#E4EAED", linewidth=0.7)
        axis.set_axisbelow(True)
        maximum = max(values) if values else 1
        axis.set_ylim(0, maximum * 1.18 if maximum else 1)
        for bar, value in zip(bars, values):
            axis.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + maximum * 0.025, f"{value:.2f}", ha="center", va="bottom", fontsize=7.8)
    fig.suptitle("Contiguity describes output; it does not establish correctness", y=1.04, fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[1])


def plot_readback(rows: list[dict[str, str]], figure_dir: Path) -> None:
    lookup = {row["Branch"]: row for row in rows}
    label_positions = {
        "flye-ont-r9": (6, 6, "left", "bottom"),
        "flye-hifi": (-8, 13, "right", "bottom"),
        "hifiasm-meta-hifi": (6, 6, "left", "bottom"),
        "metamdbg-hifi": (-8, -13, "right", "top"),
    }
    fig, axis = plt.subplots(figsize=(8.4, 5.0))
    for branch, assembler, platform in BRANCHES:
        row = lookup[branch]
        x = float(row["AlignedQueryBasePct"])
        y = float(row["WeightedAlignmentIdentityPct"])
        split = float(row["MultiContigPctOfMapped"])
        marker = "D" if platform == "ONT R9" else "o"
        axis.scatter(x, y, s=90 + split * 75, marker=marker, color=PALETTE[assembler], edgecolor="white", linewidth=1.1, zorder=3)
        dx, dy, horizontal, vertical = label_positions[branch]
        axis.annotate(
            short_label(branch).replace("\n", " / "),
            (x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=horizontal,
            va=vertical,
            fontsize=8.5,
        )
    axis.set_xlabel("Read bases aligned to ≥1 kb contigs (%)")
    axis.set_ylabel("Weighted alignment identity (%)")
    axis.grid(color="#E4EAED", linewidth=0.7)
    axis.set_axisbelow(True)
    axis.set_xlim(max(0, min(float(row["AlignedQueryBasePct"]) for row in rows) - 5), 101)
    axis.set_ylim(max(0, min(float(row["WeightedAlignmentIdentityPct"]) for row in rows) - 5), 101)
    axis.set_title(
        "Read-back evidence separates consensus agreement from fragmented mappings",
        pad=12,
        fontweight="bold",
    )
    axis.text(0.01, 0.01, "Point area increases with reads mapped across multiple contigs", transform=axis.transAxes, fontsize=8.5, color="#455A64")
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[2])


def plot_circular_resources(support: list[dict[str, str]], resources: list[dict[str, str]], figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5))
    candidate_counts = Counter(row["Branch"] for row in support)
    supported_counts = Counter(row["Branch"] for row in support if row["SupportedByGE3Reads"] == "TRUE")
    x = np.arange(4)
    labels = [short_label(branch) for branch, _, _ in BRANCHES]
    total = [candidate_counts[branch] for branch, _, _ in BRANCHES]
    supported = [supported_counts[branch] for branch, _, _ in BRANCHES]
    width = 0.36
    axes[0].bar(x - width / 2, total, width, label="Software candidates", color="#A9BBC5")
    axes[0].bar(x + width / 2, supported, width, label="≥3 junction reads", color="#CC79A7")
    axes[0].set_xticks(x, labels, rotation=25, ha="right")
    axes[0].set_ylabel("Circular candidates (count)")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(axis="y", color="#E4EAED", linewidth=0.7)
    axes[0].set_axisbelow(True)
    assembly_resources = {row["Step"].removeprefix("assemble-"): row for row in resources if row["Step"].startswith("assemble-")}
    wall = [float(assembly_resources[branch]["WallSeconds"]) / 3600 for branch, _, _ in BRANCHES]
    rss = [int(assembly_resources[branch]["PeakRSSKiB"]) / 1024**2 for branch, _, _ in BRANCHES]
    colors = [PALETTE[assembler] for _, assembler, _ in BRANCHES]
    points = axes[1].scatter(wall, rss, s=95, c=colors, edgecolor="white", linewidth=1.0)
    for index, (branch, _, _) in enumerate(BRANCHES):
        axes[1].annotate(short_label(branch).replace("\n", " / "), (wall[index], rss[index]), xytext=(5, 5), textcoords="offset points", fontsize=8)
    axes[1].set_xlabel("Cumulative wall time across attempts (h)")
    axes[1].set_ylabel("Peak resident memory (GiB)")
    axes[1].grid(color="#E4EAED", linewidth=0.7)
    axes[1].set_axisbelow(True)
    fig.suptitle("Circular flags need junction evidence; resource cost remains explicit", y=1.04, fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[3])


def audit_figures(figure_dir: Path, checks: Checks) -> None:
    for stem in FIGURE_STEMS:
        pdf = figure_dir / f"{stem}.pdf"
        png = figure_dir / f"{stem}.png"
        tiff = figure_dir / f"{stem}.tiff"
        checks.add("Figures", f"{stem}-pdf", pdf.is_file() and pdf.stat().st_size > 1000 and pdf.read_bytes()[:4] == b"%PDF", pdf.stat().st_size if pdf.exists() else 0)
        try:
            with Image.open(png) as image:
                ok = image.format == "PNG" and image.width >= 1800 and image.height >= 1000
                detail = f"{image.format} {image.width}x{image.height}"
        except Exception as error:
            ok, detail = False, repr(error)
        checks.add("Figures", f"{stem}-png", ok, detail)
        try:
            with Image.open(tiff) as image:
                compression = str(image.info.get("compression", ""))
                dpi = image.info.get("dpi", (0, 0))
                ok = image.format == "TIFF" and image.width >= 1800 and image.height >= 1000 and "lzw" in compression.lower() and min(float(dpi[0]), float(dpi[1])) >= 349
                detail = f"{image.format} {image.width}x{image.height} {compression} dpi={dpi}"
        except Exception as error:
            ok, detail = False, repr(error)
        checks.add("Figures", f"{stem}-tiff", ok, detail)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    env_prefix = args.env_prefix.resolve()
    frozen = (args.frozen_dir or root / "data/small/31-long-read-assembly-frozen").resolve()
    output = (args.output_dir or root / "results/31-long-read-assembly").resolve()
    figures = (args.figure_dir or root / "figures").resolve()
    chapter = (args.chapter or root / "chapters/31-long-read-assembly.qmd").resolve()
    output.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    required = (
        frozen / "file-checksums.sha256",
        frozen / "source-manifest.tsv",
        frozen / "source-audit.tsv",
        frozen / "software-releases.tsv",
        frozen / "assembly-metrics.tsv",
        frozen / "contig-inventory.tsv",
        frozen / "readback-metrics.tsv",
        frozen / "junction-support.tsv",
        frozen / "resource-usage.tsv",
        frozen / "tool-versions.tsv",
        frozen / "run-summary.json",
        chapter,
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("Missing Article 31 inputs: " + ", ".join(missing))

    checks = Checks()
    checksum_rows = verify_checksum_manifest(frozen, checks)
    audit_project_contract(root, frozen, checks)
    tool_rows = audit_tools(env_prefix, frozen, checks)
    source_rows = audit_sources(frozen, checks)
    metric_rows, inventory = audit_assemblies(frozen, checks)
    readback_rows = audit_readback(frozen, checks)
    support_rows = audit_junctions(frozen, inventory, checks)
    resource_rows = audit_resources(frozen, support_rows, checks)
    audit_boundaries_and_chapter(frozen, chapter, checks)

    run_summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    run_ok = (
        run_summary["article"] == 31
        and run_summary["biological_sample"] == "SAMEA14435832 / MOCK1 / 71 strains"
        and run_summary["source_runs"] == ["ERR9765780", "ERR9765783"]
        and run_summary["branches"] == [branch for branch, _, _ in BRANCHES]
        and run_summary["assembly_payloads"] == 4
        and run_summary["mapping_summaries"] == 4
    )
    checks.add("Frozen input", "run-summary", run_ok, run_summary)

    configure_plotting()
    plot_design(figures)
    plot_contiguity(metric_rows, figures)
    plot_readback(readback_rows, figures)
    plot_circular_resources(support_rows, resource_rows, figures)
    audit_figures(figures, checks)

    write_tsv(output / "checksum-audit.tsv", checksum_rows, ["File", "ExpectedSHA256", "ObservedSHA256", "Status"])
    write_tsv(output / "tool-audit.tsv", tool_rows, ["Tool", "ExpectedPattern", "ReturnCode", "Observed", "Status"])
    write_tsv(output / "source-audit.tsv", source_rows, list(source_rows[0]))
    write_tsv(output / "assembly-audit.tsv", metric_rows, list(metric_rows[0]))
    write_tsv(output / "readback-audit.tsv", readback_rows, list(readback_rows[0]))
    write_tsv(
        output / "junction-audit.tsv",
        support_rows,
        list(support_rows[0]) if support_rows else JUNCTION_AUDIT_FIELDS,
    )
    write_tsv(output / "resource-audit.tsv", resource_rows, list(resource_rows[0]))
    write_tsv(output / "validation-checks.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])

    assembly_resources = [row for row in resource_rows if row["Step"].startswith("assemble-")]
    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "seed": SEED,
        "source_project": "PRJEB52977",
        "source_sample": "SAMEA14435832",
        "source_runs": list(EXPECTED_SOURCES),
        "assembly_branches": 4,
        "mapping_branches": 4,
        "circular_candidates": len(support_rows),
        "junction_supported_candidates": sum(row["SupportedByGE3Reads"] == "TRUE" for row in support_rows),
        "assembly_wall_hours_range": [min(float(row["WallSeconds"]) for row in assembly_resources) / 3600, max(float(row["WallSeconds"]) for row in assembly_resources) / 3600],
        "assembly_peak_rss_gib_range": [min(int(row["PeakRSSKiB"]) for row in assembly_resources) / 1024**2, max(int(row["PeakRSSKiB"]) for row in assembly_resources) / 1024**2],
        "final_attempt_wall_hours_range": [min(float(row["FinalWallSeconds"]) for row in assembly_resources) / 3600, max(float(row["FinalWallSeconds"]) for row in assembly_resources) / 3600],
        "final_attempt_peak_rss_gib_range": [min(int(row["FinalPeakRSSKiB"]) for row in assembly_resources) / 1024**2, max(int(row["FinalPeakRSSKiB"]) for row in assembly_resources) / 1024**2],
        "platform_only_experiment": False,
        "universal_winner_claimed": False,
        "n50_treated_as_correctness": False,
        "circular_flag_treated_as_complete_genome": False,
        "reference_aware_accuracy_evaluated": False,
        "qa_network_access": False,
        "checksum_files": len(checksum_rows),
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "figures": [f"{stem}.{suffix}" for stem in FIGURE_STEMS for suffix in ("pdf", "png", "tiff")],
    }
    (output / "validation-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    with (output / "validation.log").open("w", encoding="utf-8") as handle:
        handle.write("Article 31 long-read assembly validation\n")
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
