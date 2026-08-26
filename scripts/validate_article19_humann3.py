#!/usr/bin/env python3
"""Initialize and validate Article 19 HUMAnN 3 evidence.

Initialization is the only mode that reads FASTQ-derived work products, large
database inventories, or resource logs. Routine QA is network-free and
database-free: it verifies the frozen checksum manifest, recomputes unit,
regrouping, and stratification invariants from small HUMAnN tables, and renders
four figures.
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
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article19-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


EXPECTED_ENV_SHA256 = {
    "env/biobakery.yml": "a314c7d33025e2057ab609f4f1e101b6ab3921338e55c587d318cf2d6d65f874",
    "env/biobakery-linux-64.lock": "d2884ede2be3ae1c27300505d5a1b6c784a4b30a9f259a3efd9b628a66d302d2",
    "env/relink-biobakery-entrypoints.sh": "3b0b785c61b8618adca4f23acd013a6b3a43063b7ae2849c0332a2c948425f1b",
}
EXPECTED_LOCK_PACKAGES = 554
EXPECTED_COMPAT_SHA256 = "a7fc07ccf0c22f6895f92e891065e54a3fcc4786fb8560c4443d19f7346e1b57"
EXPECTED_CUSTOM_CHOCO_FILES = 37
EXPECTED_HUMANN_SELECTED_SPECIES = 89
EXPECTED_REACTION_FEATURES = 1_402
EXPECTED_REGROUP_UNGROUPED_RPK = 113_812.638612224
EXPECTED_REGROUP_MAPPING = {
    "relative_path": "lib/python3.12/site-packages/humann/data/pathways/metacyc_reactions_level4ec_only.uniref.bz2",
    "bytes": 57_511_575,
    "sha256": "8419ce78a62ca9130914f2c347a9708111cedc7de52ba274659ce51ec7de7752",
    "group": "uniref90_rxn",
}
EXPECTED_INPUT = {
    "pairs": 99_991,
    "reads": 199_982,
    "pair_hash": "457cef6e9d603790dfbc26b716b0498169b54c31bc903d067d449d8dcc86770d",
    "R1": {
        "bytes": 8_661_319,
        "sha256": "ce438de2814a6e605cc24b00e53b697d018f325bc2a9694cece5be158ff32101",
        "bases": 14_974_589,
    },
    "R2": {
        "bytes": 10_045_722,
        "sha256": "207d080cd5dc98bd10fd4af8c492fb3f2000f73e1458a9d1d77e57347f4fe459",
        "bases": 14_835_184,
    },
}
EXPECTED_DATABASES = {
    "metaphlan-vjun23-metadata": {
        "bytes": 3_316_336_640,
        "algorithm": "md5",
        "checksum": "d985de75a217cd319e721863f68e7d33",
        "release": "mpa_vJun23_CHOCOPhlAnSGB_202403-metadata",
        "level": "publisher-checksum",
    },
    "metaphlan-vjun23-bowtie2": {
        "bytes": 22_989_752_320,
        "algorithm": "md5",
        "checksum": "8caae86b4d2931416cbdbb92f5985cef",
        "release": "mpa_vJun23_CHOCOPhlAnSGB_202403-bowtie2",
        "level": "publisher-checksum",
    },
    "humann-chocophlan-full": {
        "bytes": 16_502_062_909,
        "algorithm": "sha256",
        "release": "v201901_v31-full",
        "level": "retrieval-lock",
    },
    "humann-uniref90-full": {
        "bytes": 20_579_913_329,
        "algorithm": "sha256",
        "release": "v201901b-uniref90-annotated-full",
        "level": "retrieval-lock",
    },
}
EXPECTED_TOOL_VERSIONS = {
    "HUMAnN": "3.9",
    "MetaPhlAn": "4.2.5",
    "Bowtie2": "2.5.5",
    "DIAMOND": "2.2.4",
    "Python": "3.12.13",
}
SPECIAL_FEATURES = {"UNMAPPED", "UNINTEGRATED", "UNGROUPED"}
FIGURE_STEMS = (
    "19-read-flow",
    "19-gene-family-stratification",
    "19-pathway-contributions",
    "19-abundance-coverage",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--environment-prefix", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--initialize-frozen", action="store_true")
    parser.add_argument("--cache-root", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--combined-fastq", type=Path)
    parser.add_argument("--vjun-profile", type=Path)
    parser.add_argument("--incompatible-log", type=Path)
    parser.add_argument("--incompatible-status", type=Path)
    parser.add_argument("--gene-rpk", type=Path)
    parser.add_argument("--gene-cpm", type=Path)
    parser.add_argument("--gene-relab", type=Path)
    parser.add_argument("--reaction-rpk", type=Path)
    parser.add_argument("--reaction-cpm", type=Path)
    parser.add_argument("--reaction-relab", type=Path)
    parser.add_argument("--regroup-log", type=Path)
    parser.add_argument("--path-rpk", type=Path)
    parser.add_argument("--path-cpm", type=Path)
    parser.add_argument("--path-relab", type=Path)
    parser.add_argument("--path-coverage", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    return parser.parse_args()


def hash_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
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
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


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


def split_feature(feature: str) -> tuple[str, str]:
    if "|" not in feature:
        return feature, ""
    return tuple(feature.split("|", 1))  # type: ignore[return-value]


def feature_id(feature: str) -> str:
    return feature.split(": ", 1)[0]


def feature_name(feature: str) -> str:
    if ": " in feature:
        return feature.split(": ", 1)[1]
    return feature_id(feature)


def shorten(value: str, width: int = 44) -> str:
    value = value.replace("_", " ")
    return value if len(value) <= width else value[: width - 1] + "…"


def read_humann_table(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("#"):
        raise ValueError(f"Unexpected HUMAnN table header: {path}")
    header = lines[0].split("\t")
    if len(header) != 2:
        raise ValueError(f"Expected one sample column in {path}")
    rows: list[dict[str, Any]] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            raise ValueError(f"Malformed HUMAnN row in {path}: {line[:80]}")
        base, stratum = split_feature(parts[0])
        rows.append(
            {
                "Feature": parts[0],
                "BaseFeature": base,
                "FeatureID": feature_id(base),
                "FeatureName": feature_name(base),
                "Stratum": stratum,
                "Value": float(parts[1]),
                "Special": feature_id(base) in SPECIAL_FEATURES,
            }
        )
    return {"header": header, "rows": rows}


def value_map(table: dict[str, Any]) -> dict[str, float]:
    return {row["Feature"]: row["Value"] for row in table["rows"]}


def table_summary(label: str, table: dict[str, Any]) -> dict[str, Any]:
    rows = table["rows"]
    unstratified = [row for row in rows if not row["Stratum"]]
    stratified = [row for row in rows if row["Stratum"]]
    ordinary = [row for row in unstratified if not row["Special"]]
    special = [row for row in unstratified if row["Special"]]
    return {
        "Table": label,
        "SampleHeader": table["header"][1],
        "Rows": len(rows),
        "UnstratifiedRows": len(unstratified),
        "StratifiedRows": len(stratified),
        "OrdinaryCommunityFeatures": len(ordinary),
        "SpecialCommunityRows": len(special),
        "PositiveOrdinaryFeatures": sum(row["Value"] > 0 for row in ordinary),
        "CommunityValueSum": math.fsum(row["Value"] for row in unstratified),
    }


def gene_stratification_rows(
    rpk: dict[str, Any],
    cpm: dict[str, Any],
    relab: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    cpm_values = value_map(cpm)
    relab_values = value_map(relab)
    strata: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rpk["rows"]:
        if row["Stratum"]:
            strata[row["BaseFeature"]].append(row)
    output: list[dict[str, Any]] = []
    for row in rpk["rows"]:
        if row["Stratum"] or row["Special"]:
            continue
        members = strata.get(row["BaseFeature"], [])
        stratum_rpk = math.fsum(item["Value"] for item in members)
        absolute_difference = row["Value"] - stratum_rpk
        relative_difference = (
            absolute_difference / row["Value"] if row["Value"] else 0.0
        )
        dominant = max(members, key=lambda item: item["Value"], default=None)
        output.append(
            {
                "Feature": row["BaseFeature"],
                "FeatureID": row["FeatureID"],
                "FeatureName": row["FeatureName"],
                "CommunityRPK": row["Value"],
                "StrataRPKSum": stratum_rpk,
                "AbsoluteDifference": absolute_difference,
                "RelativeDifference": relative_difference,
                "CommunityCPM": cpm_values[row["BaseFeature"]],
                "CommunityRelativeAbundance": relab_values[row["BaseFeature"]],
                "Strata": len(members),
                "DominantStratum": dominant["Stratum"] if dominant else "",
                "DominantFraction": (
                    dominant["Value"] / row["Value"]
                    if dominant and row["Value"]
                    else 0.0
                ),
            }
        )
    max_abs = max((abs(row["AbsoluteDifference"]) for row in output), default=0.0)
    max_rel = max((abs(row["RelativeDifference"]) for row in output), default=0.0)
    violations = sum(
        abs(row["AbsoluteDifference"]) > max(1e-6, abs(row["CommunityRPK"]) * 1e-8)
        for row in output
    )
    return output, {
        "features": len(output),
        "max_absolute_difference": max_abs,
        "max_relative_difference": max_rel,
        "violations": violations,
    }


def community_value(table: dict[str, Any], feature: str) -> float:
    matches = [
        row["Value"]
        for row in table["rows"]
        if not row["Stratum"] and row["FeatureID"] == feature
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one unstratified {feature} row; observed {len(matches)}"
        )
    return matches[0]


def regroup_audit_row(
    prefix: Path,
    gene_rpk: dict[str, Any],
    reaction_rpk: dict[str, Any],
) -> dict[str, Any]:
    mapping = prefix / EXPECTED_REGROUP_MAPPING["relative_path"]
    if not mapping.is_file():
        raise FileNotFoundError(mapping)
    source_ordinary = [
        row
        for row in gene_rpk["rows"]
        if not row["Stratum"] and not row["Special"]
    ]
    reactions = [
        row
        for row in reaction_rpk["rows"]
        if not row["Stratum"] and not row["Special"]
    ]
    source_rpk = math.fsum(row["Value"] for row in source_ordinary)
    reaction_rpk_sum = math.fsum(row["Value"] for row in reactions)
    ungrouped_rpk = community_value(reaction_rpk, "UNGROUPED")
    source_unmapped_rpk = community_value(gene_rpk, "UNMAPPED")
    reaction_unmapped_rpk = community_value(reaction_rpk, "UNMAPPED")
    mapped_source_rpk = source_rpk - ungrouped_rpk
    reaction_plus_ungrouped_rpk = reaction_rpk_sum + ungrouped_rpk
    return {
        "MappingGroup": EXPECTED_REGROUP_MAPPING["group"],
        "MappingFile": f"${{BIOBAKERY_ENV_PREFIX}}/{EXPECTED_REGROUP_MAPPING['relative_path']}",
        "MappingBytes": mapping.stat().st_size,
        "MappingSHA256": hash_file(mapping),
        "Function": "sum",
        "UngroupedPolicy": "Y",
        "ProtectedPolicy": "Y",
        "SourceOrdinaryFeatures": len(source_ordinary),
        "SourceOrdinaryRPK": source_rpk,
        "ReactionFeatures": len(reactions),
        "ReactionRPK": reaction_rpk_sum,
        "UngroupedRPK": ungrouped_rpk,
        "MappedSourceRPK": mapped_source_rpk,
        "ReactionPlusUngroupedRPK": reaction_plus_ungrouped_rpk,
        "MappedExpansionFactor": (
            reaction_rpk_sum / mapped_source_rpk if mapped_source_rpk else 0.0
        ),
        "TotalRegroupExpansionFactor": (
            reaction_plus_ungrouped_rpk / source_rpk if source_rpk else 0.0
        ),
        "SourceUnmappedRPK": source_unmapped_rpk,
        "ReactionUnmappedRPK": reaction_unmapped_rpk,
    }


def pathway_contribution_rows(
    rpk: dict[str, Any],
    cpm: dict[str, Any],
    relab: dict[str, Any],
    coverage: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    cpm_values = value_map(cpm)
    relab_values = value_map(relab)
    coverage_values = value_map(coverage)
    strata_rpk: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    strata_cpm: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rpk["rows"]:
        if row["Stratum"]:
            strata_rpk[row["BaseFeature"]].append(row)
    for row in cpm["rows"]:
        if row["Stratum"]:
            strata_cpm[row["BaseFeature"]].append(row)
    output: list[dict[str, Any]] = []
    for row in rpk["rows"]:
        if row["Stratum"] or row["Special"]:
            continue
        members_rpk = strata_rpk.get(row["BaseFeature"], [])
        members_cpm = strata_cpm.get(row["BaseFeature"], [])
        stratum_rpk = math.fsum(item["Value"] for item in members_rpk)
        stratum_cpm = math.fsum(item["Value"] for item in members_cpm)
        dominant = max(members_rpk, key=lambda item: item["Value"], default=None)
        community = row["Value"]
        output.append(
            {
                "Pathway": row["BaseFeature"],
                "PathwayID": row["FeatureID"],
                "PathwayName": row["FeatureName"],
                "CommunityRPK": community,
                "StrataRPKSum": stratum_rpk,
                "RPKDifference": community - stratum_rpk,
                "StrataToCommunityRatio": stratum_rpk / community if community else 0.0,
                "CommunityCPM": cpm_values[row["BaseFeature"]],
                "StrataCPMSum": stratum_cpm,
                "CommunityRelativeAbundance": relab_values[row["BaseFeature"]],
                "Coverage": coverage_values[row["BaseFeature"]],
                "Strata": len(members_rpk),
                "DominantStratum": dominant["Stratum"] if dominant else "",
                "DominantFractionOfStrata": (
                    dominant["Value"] / stratum_rpk
                    if dominant and stratum_rpk
                    else 0.0
                ),
            }
        )
    nonadditive = sum(
        abs(row["RPKDifference"]) > max(1e-6, abs(row["CommunityRPK"]) * 1e-8)
        for row in output
    )
    out_of_range = sum(not 0 <= row["Coverage"] <= 1 for row in output)
    ratios = sorted(row["StrataToCommunityRatio"] for row in output)
    return output, {
        "features": len(output),
        "nonadditive_features": nonadditive,
        "coverage_out_of_range": out_of_range,
        "median_strata_to_community_ratio": (
            ratios[len(ratios) // 2] if ratios else 0.0
        ),
    }


def count_fasta(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(line.startswith(">") for line in handle)


def count_fastq_and_ids(path: Path) -> tuple[int, int]:
    records = 0
    identifiers: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"Truncated FASTQ: {path}")
            records += 1
            identifiers.add(header[1:].strip().replace(" ", ""))
    return records, len(identifiers)


def parse_profile(path: Path) -> tuple[str, dict[str, float]]:
    header_lines: list[str] = []
    species: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            header_lines.append(line)
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        lineage = fields[0]
        taxa = lineage.split("|")
        if not taxa or not taxa[-1].startswith("s__"):
            continue
        species[taxa[-1]] = float(fields[2])
    return "\n".join(header_lines), species


def prescreen_release_rows(
    project_root: Path, vjun_profile: Path
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    header, vjun = parse_profile(vjun_profile)
    vjan_rows = read_tsv(
        project_root / "data/small/15-metaphlan-frozen/species-profile.tsv"
    )
    vjan = {
        row["SpeciesLabel"]: float(row["RelativeAbundancePct"])
        for row in vjan_rows
        if row["SpeciesLabel"] != "UNCLASSIFIED"
    }
    rows: list[dict[str, Any]] = []
    for species in sorted(set(vjun) | set(vjan)):
        vjun_value = vjun.get(species, 0.0)
        vjan_value = vjan.get(species, 0.0)
        selected_vjun = vjun_value >= 0.01
        selected_vjan = vjan_value >= 0.01
        status = (
            "Both"
            if selected_vjun and selected_vjan
            else "vJun23 only"
            if selected_vjun
            else "vJan26 only"
            if selected_vjan
            else "Below threshold"
        )
        rows.append(
            {
                "Species": species,
                "vJun23RelativeAbundancePct": vjun_value,
                "vJan26RelativeAbundancePct": vjan_value,
                "SelectedByvJun23": "Yes" if selected_vjun else "No",
                "SelectedByvJan26": "Yes" if selected_vjan else "No",
                "SelectionStatus": status,
            }
        )
    selected_jun = {key for key, value in vjun.items() if value >= 0.01}
    selected_jan = {key for key, value in vjan.items() if value >= 0.01}
    union = selected_jun | selected_jan
    return rows, {
        "profile_header_has_vjun23": "vJun23" in header,
        "vjun23_species": len(vjun),
        "vjan26_species": len(vjan),
        "vjun23_selected_species": len(selected_jun),
        "vjan26_selected_species": len(selected_jan),
        "selected_intersection": len(selected_jun & selected_jan),
        "selected_union": len(union),
        "selected_jaccard": len(selected_jun & selected_jan) / len(union) if union else 1.0,
    }


def parse_resource(path: Path, stage: str, threads: int) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    elapsed_match = re.search(
        r"^\s*Elapsed \(wall clock\) time.*\):\s*([0-9]+(?::[0-9]+)*(?:\.[0-9]+)?)\s*$",
        text,
        flags=re.MULTILINE,
    )
    rss_match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    elapsed_text = elapsed_match.group(1).strip() if elapsed_match else ""
    seconds = 0.0
    if elapsed_text:
        parts = [float(value) for value in elapsed_text.split(":")]
        for value in parts:
            seconds = seconds * 60 + value
    rss_kib = int(rss_match.group(1)) if rss_match else 0
    return {
        "Stage": stage,
        "Threads": threads,
        "WallSeconds": seconds,
        "MaximumRSSKiB": rss_kib,
        "MaximumRSSGiB": rss_kib / 1024**2,
    }


def tool_versions(prefix: Path) -> list[dict[str, str]]:
    env = os.environ.copy()
    env["PATH"] = f"{prefix / 'bin'}:{env.get('PATH', '')}"
    env["PYTHONNOUSERSITE"] = "1"

    def run(command: list[str]) -> str:
        result = subprocess.run(
            command, capture_output=True, text=True, env=env, check=False
        )
        return result.stdout + result.stderr

    values = {
        "HUMAnN": (
            re.search(r"humann v([0-9.]+)", run([str(prefix / "bin/humann"), "--version"]))
            or [None, "missing"]
        )[1],
        "MetaPhlAn": (
            re.search(
                r"MetaPhlAn version ([0-9.]+)",
                run([str(prefix / "bin/metaphlan"), "--version"]),
            )
            or [None, "missing"]
        )[1],
        "Bowtie2": (
            re.search(
                r"version ([0-9.]+)", run([str(prefix / "bin/bowtie2"), "--version"])
            )
            or [None, "missing"]
        )[1],
        "DIAMOND": (
            re.search(
                r"diamond version ([0-9.]+)",
                run([str(prefix / "bin/diamond"), "version"]),
            )
            or [None, "missing"]
        )[1],
        "Python": (
            re.search(
                r"Python ([0-9.]+)", run([str(prefix / "bin/python"), "--version"])
            )
            or [None, "missing"]
        )[1],
    }
    return [
        {
            "Tool": tool,
            "Version": version,
            "Executable": f"${{BIOBAKERY_ENV_PREFIX}}/bin/{tool.lower()}",
        }
        for tool, version in values.items()
    ]


def normalize_log(text: str, replacements: dict[str, str]) -> str:
    for source, replacement in sorted(
        replacements.items(), key=lambda item: len(item[0]), reverse=True
    ):
        text = text.replace(source, replacement)
    return text


def database_audit_rows(
    project_root: Path, cache_root: Path
) -> list[dict[str, Any]]:
    manifest = {
        row["database_id"]: row
        for row in read_tsv(project_root / "data/small/19-database-manifest.tsv")
    }
    bootstrap = {
        row["database_id"]: row
        for row in read_tsv(cache_root / "manifests/bootstrap-results.tsv")
    }
    inventories = {
        row["database_id"]: row
        for row in read_tsv(cache_root / "manifests/installed-inventory.tsv")
    }
    output: list[dict[str, Any]] = []
    for database_id in EXPECTED_DATABASES:
        row = manifest[database_id]
        acquired = bootstrap[database_id]
        inventory_key = (
            "metaphlan-vjun23"
            if database_id.startswith("metaphlan-")
            else database_id
        )
        inventory = inventories[inventory_key]
        output.append(
            {
                "DatabaseID": database_id,
                "Tool": row["tool"],
                "ToolVersion": row["tool_version"],
                "ReleaseID": row["release_id"],
                "ArchiveURL": row["archive_url"],
                "EffectiveURL": acquired["effective_url"],
                "RemoteContentLength": int(acquired["remote_content_length"]),
                "ArchiveBytes": int(acquired["bytes"]),
                "ChecksumAlgorithm": row["checksum_algorithm"],
                "ExpectedChecksum": row["expected_checksum"],
                "PublisherChecksumAlgorithm": acquired[
                    "publisher_checksum_algorithm"
                ],
                "PublisherExpectedChecksum": acquired["publisher_checksum"],
                "PublisherObservedChecksum": acquired[
                    "observed_publisher_checksum"
                ],
                "ObservedSHA256": acquired["observed_sha256"],
                "VerificationLevel": EXPECTED_DATABASES[database_id]["level"],
                "RetrievedAtUTC": acquired["retrieved_at_utc"],
                "ArchiveIntegrity": acquired["archive_integrity"],
                "InstalledFiles": int(inventory["files"]),
                "InstalledBytes": int(inventory["installed_bytes"]),
                "InstalledPath": inventory["installed_path"].replace(
                    str(cache_root), "${HUMANN_CACHE_ROOT}"
                ),
                "DownloadGate": row["download_gate"],
                "ValidationStatus": row["validation_status"],
            }
        )
    return output


def write_database_file_manifest(cache_root: Path, frozen_dir: Path) -> int:
    source = cache_root / "manifests/installed-files.sha256"
    output = frozen_dir / "database-files.sha256"
    lines: list[str] = []
    for line in source.read_text(encoding="utf-8").splitlines():
        digest, path = line.split(maxsplit=1)
        normalized = path.replace(str(cache_root), "${HUMANN_CACHE_ROOT}")
        lines.append(f"{digest}  {normalized}")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines)


def write_frozen_checksums(frozen_dir: Path) -> int:
    manifest = frozen_dir / "file-checksums.sha256"
    files = sorted(
        path
        for path in frozen_dir.rglob("*")
        if path.is_file() and path != manifest
    )
    lines = [f"{hash_file(path)}  {path.relative_to(frozen_dir).as_posix()}" for path in files]
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines)


def verify_frozen_checksums(frozen_dir: Path, checks: Checks) -> int:
    manifest = frozen_dir / "file-checksums.sha256"
    entries = 0
    seen: set[str] = set()
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        path = frozen_dir / relative
        observed = hash_file(path) if path.is_file() else "missing"
        checks.add("frozen", f"sha256-{relative}", observed == expected, observed)
        entries += 1
        seen.add(relative)
    payloads = {
        path.relative_to(frozen_dir).as_posix()
        for path in frozen_dir.rglob("*")
        if path.is_file() and path != manifest
    }
    checks.add("frozen", "manifest-complete", payloads == seen, f"payloads={len(payloads)};entries={len(seen)}")
    return entries


def initialize_frozen(args: argparse.Namespace) -> None:
    required = {
        "cache_root": args.cache_root,
        "work_dir": args.work_dir,
        "combined_fastq": args.combined_fastq,
        "vjun_profile": args.vjun_profile,
        "incompatible_log": args.incompatible_log,
        "incompatible_status": args.incompatible_status,
        "gene_rpk": args.gene_rpk,
        "gene_cpm": args.gene_cpm,
        "gene_relab": args.gene_relab,
        "reaction_rpk": args.reaction_rpk,
        "reaction_cpm": args.reaction_cpm,
        "reaction_relab": args.reaction_relab,
        "regroup_log": args.regroup_log,
        "path_rpk": args.path_rpk,
        "path_cpm": args.path_cpm,
        "path_relab": args.path_relab,
        "path_coverage": args.path_coverage,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise SystemExit(f"Initialization missing arguments: {', '.join(missing)}")
    root = args.project_root.resolve()
    prefix = args.environment_prefix.resolve()
    frozen = args.frozen_dir.resolve()
    cache = args.cache_root.resolve()
    work = args.work_dir.resolve()
    frozen.mkdir(parents=True, exist_ok=True)

    copies = {
        args.vjun_profile: "prescreen-vJun23.tsv",
        args.gene_rpk: "genefamilies-rpk.tsv",
        args.gene_cpm: "genefamilies-cpm.tsv",
        args.gene_relab: "genefamilies-relab.tsv",
        args.reaction_rpk: "reactions-rpk.tsv",
        args.reaction_cpm: "reactions-cpm.tsv",
        args.reaction_relab: "reactions-relab.tsv",
        args.path_rpk: "pathabundance-rpk.tsv",
        args.path_cpm: "pathabundance-cpm.tsv",
        args.path_relab: "pathabundance-relab.tsv",
        args.path_coverage: "pathcoverage.tsv",
    }
    for source, name in copies.items():
        assert source is not None
        shutil.copy2(source, frozen / name)

    replacements = {
        str(root): "${PROJECT_ROOT}",
        str(prefix): "${BIOBAKERY_ENV_PREFIX}",
        str(cache): "${HUMANN_CACHE_ROOT}",
        str(work): "${ARTICLE19_WORK_DIR}",
    }
    log_sources = {
        work / "logs/metaphlan-vjun23.log": "logs/metaphlan-vjun23.log",
        work / "logs/ERR9765746-humann3.log": "logs/humann3.log",
        args.regroup_log: "logs/humann-regroup-uniref90-rxn.log",
        args.incompatible_log: "logs/humann-vJan26-expected-rejection.log",
        args.incompatible_status: "logs/humann-vJan26-expected-rejection.status",
        work / "logs/metaphlan-vjun23.resources.txt": "logs/metaphlan-vjun23.resources.txt",
        work / "logs/humann3.resources.txt": "logs/humann3.resources.txt",
    }
    for source, relative in log_sources.items():
        assert source is not None
        if not source.is_file():
            raise FileNotFoundError(source)
        target = frozen / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            normalize_log(source.read_text(encoding="utf-8", errors="replace"), replacements),
            encoding="utf-8",
        )

    versions = tool_versions(prefix)
    write_tsv(frozen / "tool-versions.tsv", versions, ["Tool", "Version", "Executable"])

    database = database_audit_rows(root, cache)
    write_tsv(frozen / "database-audit.tsv", database, list(database[0]))
    database_checksum_entries = write_database_file_manifest(cache, frozen)

    profile_rows, profile_summary = prescreen_release_rows(root, args.vjun_profile)
    write_tsv(
        frozen / "prescreen-release-audit.tsv", profile_rows, list(profile_rows[0])
    )

    input_reads, unique_ids = count_fastq_and_ids(args.combined_fastq)
    humann_temp = work / "humann/ERR9765746_MOCK1_humann_temp"
    nucleotide_candidates = sorted(humann_temp.glob("*_bowtie2_unaligned.fa"))
    translated_candidates = sorted(humann_temp.glob("*_diamond_unaligned.fa"))
    if len(nucleotide_candidates) != 1 or len(translated_candidates) != 1:
        raise SystemExit(
            "Expected exactly one nucleotide and one translated unaligned FASTA: "
            f"{nucleotide_candidates}; {translated_candidates}"
        )
    nucleotide_unaligned = count_fasta(nucleotide_candidates[0])
    final_unaligned = count_fasta(translated_candidates[0])
    nucleotide_mapped = input_reads - nucleotide_unaligned
    translated_mapped = nucleotide_unaligned - final_unaligned
    read_flow = [
        {
            "Stage": "Input",
            "Reads": input_reads,
            "FractionOfInput": 1.0,
            "Definition": "Concatenated R1 and R2 individual reads",
        },
        {
            "Stage": "Nucleotide mapped",
            "Reads": nucleotide_mapped,
            "FractionOfInput": nucleotide_mapped / input_reads,
            "Definition": "Accepted ChocoPhlAn nucleotide alignments",
        },
        {
            "Stage": "Translated mapped",
            "Reads": translated_mapped,
            "FractionOfInput": translated_mapped / input_reads,
            "Definition": "Accepted UniRef90 DIAMOND alignments after nucleotide search",
        },
        {
            "Stage": "Unmapped",
            "Reads": final_unaligned,
            "FractionOfInput": final_unaligned / input_reads,
            "Definition": "Reads unaligned after both search tiers",
        },
    ]
    write_tsv(frozen / "read-flow.tsv", read_flow, list(read_flow[0]))

    gene_rpk = read_humann_table(args.gene_rpk)
    gene_cpm = read_humann_table(args.gene_cpm)
    gene_relab = read_humann_table(args.gene_relab)
    reaction_rpk = read_humann_table(args.reaction_rpk)
    reaction_cpm = read_humann_table(args.reaction_cpm)
    reaction_relab = read_humann_table(args.reaction_relab)
    path_rpk = read_humann_table(args.path_rpk)
    path_cpm = read_humann_table(args.path_cpm)
    path_relab = read_humann_table(args.path_relab)
    coverage = read_humann_table(args.path_coverage)

    feature_summary = [
        table_summary("Gene families RPK", gene_rpk),
        table_summary("Gene families CPM", gene_cpm),
        table_summary("Gene families relative abundance", gene_relab),
        table_summary("Regrouped reactions RPK", reaction_rpk),
        table_summary("Regrouped reactions CPM", reaction_cpm),
        table_summary("Regrouped reactions relative abundance", reaction_relab),
        table_summary("Pathway abundance RPK-like", path_rpk),
        table_summary("Pathway abundance CPM", path_cpm),
        table_summary("Pathway abundance relative abundance", path_relab),
        table_summary("Pathway coverage", coverage),
    ]
    write_tsv(frozen / "feature-summary.tsv", feature_summary, list(feature_summary[0]))

    genes, gene_summary = gene_stratification_rows(gene_rpk, gene_cpm, gene_relab)
    reactions, reaction_summary = gene_stratification_rows(
        reaction_rpk, reaction_cpm, reaction_relab
    )
    regroup_audit = regroup_audit_row(prefix, gene_rpk, reaction_rpk)
    write_tsv(frozen / "regroup-audit.tsv", [regroup_audit], list(regroup_audit))
    pathways, pathway_summary = pathway_contribution_rows(
        path_rpk, path_cpm, path_relab, coverage
    )
    top_genes = sorted(genes, key=lambda row: row["CommunityCPM"], reverse=True)[:25]
    top_pathways = sorted(pathways, key=lambda row: row["CommunityCPM"], reverse=True)[:25]
    write_tsv(frozen / "top-genefamilies.tsv", top_genes, list(top_genes[0]))
    write_tsv(frozen / "pathway-contributions.tsv", pathways, list(pathways[0]))
    write_tsv(frozen / "top-pathways.tsv", top_pathways, list(top_pathways[0]))

    stratification = [
        {
            "Output": "Gene families",
            "ExpectedRelationship": "Community equals sum of strata",
            "Features": gene_summary["features"],
            "Violations": gene_summary["violations"],
            "MaximumAbsoluteDifference": gene_summary["max_absolute_difference"],
            "MaximumRelativeDifference": gene_summary["max_relative_difference"],
            "NonadditiveFeatures": 0,
            "MedianStrataToCommunityRatio": 1.0,
        },
        {
            "Output": "Regrouped reactions",
            "ExpectedRelationship": "Community equals sum of strata",
            "Features": reaction_summary["features"],
            "Violations": reaction_summary["violations"],
            "MaximumAbsoluteDifference": reaction_summary[
                "max_absolute_difference"
            ],
            "MaximumRelativeDifference": reaction_summary[
                "max_relative_difference"
            ],
            "NonadditiveFeatures": 0,
            "MedianStrataToCommunityRatio": 1.0,
        },
        {
            "Output": "Pathway abundance",
            "ExpectedRelationship": "Community and strata solved separately",
            "Features": pathway_summary["features"],
            "Violations": 0,
            "MaximumAbsoluteDifference": max(
                (abs(row["RPKDifference"]) for row in pathways), default=0.0
            ),
            "MaximumRelativeDifference": max(
                (
                    abs(row["RPKDifference"]) / row["CommunityRPK"]
                    if row["CommunityRPK"]
                    else 0.0
                    for row in pathways
                ),
                default=0.0,
            ),
            "NonadditiveFeatures": pathway_summary["nonadditive_features"],
            "MedianStrataToCommunityRatio": pathway_summary[
                "median_strata_to_community_ratio"
            ],
        },
    ]
    write_tsv(frozen / "stratification-audit.tsv", stratification, list(stratification[0]))

    units = [
        {
            "Output": "Gene families",
            "NativeUnit": "RPK",
            "DerivedUnits": "CPM; relative abundance",
            "SpecialRows": "UNMAPPED",
            "CanSumStrataToCommunity": "Yes",
            "RenormalizeCoverage": "Not applicable",
            "Interpretation": "Length-normalized gene-family abundance",
        },
        {
            "Output": "Regrouped reactions",
            "NativeUnit": "RPK regrouped by sum",
            "DerivedUnits": "CPM; relative abundance",
            "SpecialRows": "UNMAPPED; UNGROUPED",
            "CanSumStrataToCommunity": "Yes",
            "RenormalizeCoverage": "Not applicable",
            "Interpretation": "UniRef90-to-MetaCyc-reaction alternate feature space; one-to-many mapping can expand totals",
        },
        {
            "Output": "Pathway abundance",
            "NativeUnit": "RPK-derived pathway abundance",
            "DerivedUnits": "CPM; relative abundance",
            "SpecialRows": "UNMAPPED; UNINTEGRATED",
            "CanSumStrataToCommunity": "No",
            "RenormalizeCoverage": "Not applicable",
            "Interpretation": "Complete-pathway copy-like abundance from reaction structure",
        },
        {
            "Output": "Pathway coverage",
            "NativeUnit": "0-1 reaction-evidence score",
            "DerivedUnits": "None",
            "SpecialRows": "UNMAPPED; UNINTEGRATED are sentinels",
            "CanSumStrataToCommunity": "No",
            "RenormalizeCoverage": "No",
            "Interpretation": "Relative pathway evidence, not abundance or read fraction",
        },
    ]
    write_tsv(frozen / "units-contract.tsv", units, list(units[0]))

    resources = [
        parse_resource(work / "logs/metaphlan-vjun23.resources.txt", "MetaPhlAn vJun23 prescreen", 8),
        parse_resource(work / "logs/humann3.resources.txt", "HUMAnN 3 full workflow", 8),
    ]
    write_tsv(frozen / "resource-usage.tsv", resources, list(resources[0]))

    incompatibility_status = int(args.incompatible_status.read_text().strip())
    incompatibility_text = args.incompatible_log.read_text(
        encoding="utf-8", errors="replace"
    )
    humann_log_text = (work / "logs/ERR9765746-humann3.log").read_text(
        encoding="utf-8", errors="replace"
    )
    custom_chocophlan_files = len(
        re.findall(r"Adding file to database:", humann_log_text)
    )
    selected_species_match = re.search(
        r"Total species selected from prescreen:\s*(\d+)", humann_log_text
    )
    humann_selected_species = (
        int(selected_species_match.group(1)) if selected_species_match else 0
    )
    summary = {
        "status": "passed",
        "sample_id": "ERR9765746_MOCK1",
        "input_pairs": EXPECTED_INPUT["pairs"],
        "input_reads": input_reads,
        "combined_unique_read_ids": unique_ids,
        "combined_fastq_sha256": hash_file(args.combined_fastq),
        "humann_version": EXPECTED_TOOL_VERSIONS["HUMAnN"],
        "metaphlan_version": EXPECTED_TOOL_VERSIONS["MetaPhlAn"],
        "metaphlan_database": "mpa_vJun23_CHOCOPhlAnSGB_202403",
        "bowtie2_version": EXPECTED_TOOL_VERSIONS["Bowtie2"],
        "diamond_version": EXPECTED_TOOL_VERSIONS["DIAMOND"],
        "chocophlan_release": "v201901_v31",
        "uniref90_release": "v201901b",
        "metacyc_pathway_database": "metacyc_pathways_structured_filtered_v24_subreactions",
        "threads": 8,
        "pythonhashseed": 0,
        "input_read_flow": {
            "nucleotide_mapped": nucleotide_mapped,
            "translated_mapped": translated_mapped,
            "unmapped": final_unaligned,
        },
        "prescreen": profile_summary,
        "vjan26_direct_input_exit_status": incompatibility_status,
        "vjan26_expected_rejection": (
            incompatibility_status != 0 and "v3 or vJun23" in incompatibility_text
        ),
        "humann_bypass_prescreen": "bypass prescreen = True" in humann_log_text,
        "humann_selected_species_after_alias_mapping": humann_selected_species,
        "custom_chocophlan_files": custom_chocophlan_files,
        "gene_family_features": gene_summary["features"],
        "gene_family_stratification_violations": gene_summary["violations"],
        "regroup_mapping_group": regroup_audit["MappingGroup"],
        "regroup_mapping_file": regroup_audit["MappingFile"],
        "regroup_mapping_bytes": regroup_audit["MappingBytes"],
        "regroup_mapping_sha256": regroup_audit["MappingSHA256"],
        "reaction_features": reaction_summary["features"],
        "reaction_stratification_violations": reaction_summary["violations"],
        "regroup_ungrouped_rpk": regroup_audit["UngroupedRPK"],
        "regroup_mapped_expansion_factor": regroup_audit[
            "MappedExpansionFactor"
        ],
        "regroup_total_expansion_factor": regroup_audit[
            "TotalRegroupExpansionFactor"
        ],
        "pathway_features": pathway_summary["features"],
        "pathway_nonadditive_features": pathway_summary["nonadditive_features"],
        "coverage_out_of_range": pathway_summary["coverage_out_of_range"],
        "database_checksum_entries": database_checksum_entries,
        "resources": resources,
    }
    write_json(frozen / "run-summary.json", summary)
    entries = write_frozen_checksums(frozen)
    print(
        f"PASS Article 19 initialization: {entries} frozen payloads; "
        f"{gene_summary['features']} gene families; "
        f"{reaction_summary['features']} reactions; "
        f"{pathway_summary['features']} pathways"
    )


def environment_checks(project_root: Path, prefix: Path, checks: Checks) -> dict[str, str]:
    for relative, expected in EXPECTED_ENV_SHA256.items():
        path = project_root / relative
        observed = hash_file(path) if path.is_file() else "missing"
        checks.add("environment", f"sha256-{path.name}", observed == expected, observed)
    compat = project_root / "scripts/humann39-compat-bin/metaphlan"
    compat_sha256 = hash_file(compat) if compat.is_file() else "missing"
    checks.add(
        "environment",
        "metaphlan-version-wrapper-sha256",
        compat_sha256 == EXPECTED_COMPAT_SHA256,
        compat_sha256,
    )
    regroup_mapping = prefix / EXPECTED_REGROUP_MAPPING["relative_path"]
    regroup_mapping_bytes = (
        regroup_mapping.stat().st_size if regroup_mapping.is_file() else -1
    )
    regroup_mapping_sha256 = (
        hash_file(regroup_mapping) if regroup_mapping.is_file() else "missing"
    )
    checks.add(
        "environment",
        "regroup-mapping-bytes",
        regroup_mapping_bytes == EXPECTED_REGROUP_MAPPING["bytes"],
        regroup_mapping_bytes,
    )
    checks.add(
        "environment",
        "regroup-mapping-sha256",
        regroup_mapping_sha256 == EXPECTED_REGROUP_MAPPING["sha256"],
        regroup_mapping_sha256,
    )
    lock_lines = (project_root / "env/biobakery-linux-64.lock").read_text(
        encoding="utf-8"
    ).splitlines()
    package_lines = [line for line in lock_lines if line.startswith("http")]
    checks.add(
        "environment",
        "lock-package-count",
        len(package_lines) == EXPECTED_LOCK_PACKAGES,
        len(package_lines),
    )
    markers = {
        "humann": "humann-3.9-",
        "metaphlan": "metaphlan-4.2.5-",
        "bowtie2": "bowtie2-2.5.5-",
        "diamond": "diamond-2.2.4-",
        "python": "/python-3.12.13-",
    }
    for package, marker in markers.items():
        matches = [line for line in package_lines if marker in line]
        checks.add("environment", f"lock-{package}", len(matches) == 1, len(matches))
    versions = {row["Tool"]: row["Version"] for row in tool_versions(prefix)}
    for tool, expected in EXPECTED_TOOL_VERSIONS.items():
        checks.add(
            "environment", f"version-{tool.lower()}", versions.get(tool) == expected, versions.get(tool)
        )
    return versions


def source_checks(project_root: Path, checks: Checks) -> list[dict[str, Any]]:
    source = read_tsv(project_root / "data/small/19-source-manifest.tsv")
    checks.add("source", "two-mate-rows", len(source) == 2, len(source))
    checks.add(
        "source",
        "same-run",
        {row["RunAccession"] for row in source} == {"ERR9765746"},
        sorted({row["RunAccession"] for row in source}),
    )
    checks.add(
        "source",
        "mate-labels",
        {row["Mate"] for row in source} == {"R1", "R2"},
        sorted({row["Mate"] for row in source}),
    )
    for row in source:
        mate = row["Mate"]
        expected = EXPECTED_INPUT[mate]
        checks.add("source", f"{mate}-bytes", int(row["CompressedBytes"]) == expected["bytes"], row["CompressedBytes"])
        checks.add("source", f"{mate}-sha256", row["CompressedSHA256"] == expected["sha256"], row["CompressedSHA256"])
        checks.add("source", f"{mate}-records", int(row["Records"]) == EXPECTED_INPUT["pairs"], row["Records"])
        checks.add("source", f"{mate}-bases", int(row["Bases"]) == expected["bases"], row["Bases"])
        checks.add("source", f"{mate}-pair-hash", row["PairIDHash"] == EXPECTED_INPUT["pair_hash"], row["PairIDHash"])
    qc = json.loads(
        (project_root / "data/small/13-qc-frozen/run-summary.json").read_text(
            encoding="utf-8"
        )
    )
    meta = json.loads(
        (project_root / "data/small/15-metaphlan-frozen/run-summary.json").read_text(
            encoding="utf-8"
        )
    )
    checks.add("source", "article13-passed", qc["status"] == "passed", qc["status"])
    checks.add("source", "article13-pairs", int(qc["retained_pairs"]) == EXPECTED_INPUT["pairs"], qc["retained_pairs"])
    checks.add("source", "article15-passed", meta["status"] == "passed", meta["status"])
    checks.add("source", "article15-pairs", int(meta["input_pairs"]) == EXPECTED_INPUT["pairs"], meta["input_pairs"])
    return [
        {
            "Mate": row["Mate"],
            "RunAccession": row["RunAccession"],
            "Records": int(row["Records"]),
            "Bases": int(row["Bases"]),
            "CompressedBytes": int(row["CompressedBytes"]),
            "CompressedSHA256": row["CompressedSHA256"],
            "PairIDHash": row["PairIDHash"],
        }
        for row in source
    ]


def database_checks(
    project_root: Path, frozen_dir: Path, checks: Checks
) -> list[dict[str, str]]:
    manifest = {
        row["database_id"]: row
        for row in read_tsv(project_root / "data/small/19-database-manifest.tsv")
    }
    audit = read_tsv(frozen_dir / "database-audit.tsv")
    by_id = {row["DatabaseID"]: row for row in audit}
    checks.add("database", "four-assets", set(by_id) == set(EXPECTED_DATABASES), sorted(by_id))
    for database_id, expected in EXPECTED_DATABASES.items():
        row = manifest[database_id]
        frozen = by_id[database_id]
        checks.add("database", f"{database_id}-enabled", row["download_gate"] == "enabled", row["download_gate"])
        checks.add("database", f"{database_id}-release", row["release_id"] == expected["release"], row["release_id"])
        checks.add("database", f"{database_id}-bytes", int(row["expected_compressed_bytes"]) == expected["bytes"] == int(frozen["ArchiveBytes"]), frozen["ArchiveBytes"])
        checks.add(
            "database",
            f"{database_id}-remote-content-length",
            int(frozen["RemoteContentLength"]) == expected["bytes"],
            frozen["RemoteContentLength"],
        )
        checks.add("database", f"{database_id}-algorithm", row["checksum_algorithm"] == expected["algorithm"], row["checksum_algorithm"])
        checksum = row["expected_checksum"]
        valid_checksum = bool(re.fullmatch(r"[0-9a-f]{32}", checksum)) if expected["algorithm"] == "md5" else bool(re.fullmatch(r"[0-9a-f]{64}", checksum))
        checks.add("database", f"{database_id}-checksum-shape", valid_checksum, checksum)
        if "checksum" in expected:
            checks.add("database", f"{database_id}-publisher-checksum", checksum == expected["checksum"], checksum)
            checks.add(
                "database",
                f"{database_id}-observed-publisher-checksum",
                frozen["PublisherObservedChecksum"] == expected["checksum"],
                frozen["PublisherObservedChecksum"],
            )
        else:
            checks.add(
                "database",
                f"{database_id}-retrieval-sha256",
                checksum == frozen["ObservedSHA256"],
                frozen["ObservedSHA256"],
            )
        checks.add("database", f"{database_id}-integrity", frozen["ArchiveIntegrity"] == "tar-list-pass", frozen["ArchiveIntegrity"])
        checks.add("database", f"{database_id}-installed", int(frozen["InstalledFiles"]) > 0 and int(frozen["InstalledBytes"]) > 0, f"{frozen['InstalledFiles']};{frozen['InstalledBytes']}")
        checks.add(
            "database",
            f"{database_id}-installed-bytes",
            int(row["expected_installed_bytes"]) > 0
            and int(row["expected_installed_bytes"]) == int(frozen["InstalledBytes"]),
            frozen["InstalledBytes"],
        )
        checks.add("database", f"{database_id}-verification-level", frozen["VerificationLevel"] == expected["level"], frozen["VerificationLevel"])
        if expected["level"] == "retrieval-lock":
            checks.add("database", f"{database_id}-not-mislabeled-publisher", "RETRIEVAL" in row["validation_status"].upper() and "PUBLISHER" not in row["validation_status"].upper(), row["validation_status"])
    db_hashes = (frozen_dir / "database-files.sha256").read_text(encoding="utf-8").splitlines()
    checks.add("database", "installed-file-hashes", len(db_hashes) > 100, len(db_hashes))
    return audit


def normalization_checks(
    gene_rpk: dict[str, Any],
    gene_cpm: dict[str, Any],
    gene_relab: dict[str, Any],
    reaction_rpk: dict[str, Any],
    reaction_cpm: dict[str, Any],
    reaction_relab: dict[str, Any],
    path_rpk: dict[str, Any],
    path_cpm: dict[str, Any],
    path_relab: dict[str, Any],
    coverage: dict[str, Any],
    checks: Checks,
) -> list[dict[str, Any]]:
    tables = [
        ("Gene families RPK", gene_rpk),
        ("Gene families CPM", gene_cpm),
        ("Gene families relative abundance", gene_relab),
        ("Regrouped reactions RPK", reaction_rpk),
        ("Regrouped reactions CPM", reaction_cpm),
        ("Regrouped reactions relative abundance", reaction_relab),
        ("Pathway abundance RPK-like", path_rpk),
        ("Pathway abundance CPM", path_cpm),
        ("Pathway abundance relative abundance", path_relab),
        ("Pathway coverage", coverage),
    ]
    summaries = [table_summary(label, table) for label, table in tables]
    for summary in summaries:
        checks.add("features", f"{summary['Table']}-nonempty", int(summary["Rows"]) > 0, summary["Rows"])
    checks.add("units", "gene-cpm-sum", abs(summaries[1]["CommunityValueSum"] - 1_000_000) < 1.0, summaries[1]["CommunityValueSum"])
    checks.add("units", "gene-relab-sum", abs(summaries[2]["CommunityValueSum"] - 1.0) < 1e-6, summaries[2]["CommunityValueSum"])
    checks.add("units", "reaction-cpm-sum", abs(summaries[4]["CommunityValueSum"] - 1_000_000) < 1.0, summaries[4]["CommunityValueSum"])
    checks.add("units", "reaction-relab-sum", abs(summaries[5]["CommunityValueSum"] - 1.0) < 1e-6, summaries[5]["CommunityValueSum"])
    checks.add("units", "path-cpm-sum", abs(summaries[7]["CommunityValueSum"] - 1_000_000) < 1.0, summaries[7]["CommunityValueSum"])
    checks.add("units", "path-relab-sum", abs(summaries[8]["CommunityValueSum"] - 1.0) < 1e-6, summaries[8]["CommunityValueSum"])
    coverage_ordinary = [row for row in coverage["rows"] if not row["Special"]]
    out_of_range = sum(not 0 <= row["Value"] <= 1 for row in coverage_ordinary)
    checks.add("units", "coverage-zero-to-one", out_of_range == 0, out_of_range)
    checks.add("units", "coverage-not-normalized", abs(summaries[9]["CommunityValueSum"] - 1.0) > 1e-3, summaries[9]["CommunityValueSum"])
    return summaries


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / f"{stem}.png", dpi=350, bbox_inches="tight")
    fig.savefig(
        figure_dir / f"{stem}.tiff",
        dpi=350,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def render_read_flow(rows: list[dict[str, str]], figure_dir: Path) -> None:
    configure_plot_style()
    mapped = {row["Stage"]: int(row["Reads"]) for row in rows}
    labels = ["Nucleotide mapped", "Translated mapped", "Unmapped"]
    values = [mapped[label] for label in labels]
    colors = ["#2A9D8F", "#457B9D", "#B7B7A4"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.8), gridspec_kw={"width_ratios": [1.1, 1.6]})
    axes[0].barh(labels, values, color=colors, height=0.62)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Individual reads")
    axes[0].set_title("A  Read fate after two search tiers", loc="left", fontweight="bold")
    tick_values = [0, 30_000, 60_000, 90_000, 120_000]
    axes[0].set_xticks(tick_values, [f"{value:,}" for value in tick_values])
    axes[0].set_xlim(0, max(values) * 1.13)
    axes[0].grid(axis="x", color="#DDDDDD", linewidth=0.6)
    axes[0].set_axisbelow(True)
    for index, value in enumerate(values):
        axes[0].text(value + max(values) * 0.02, index, f"{value:,}", va="center")
    left = 0
    input_reads = mapped["Input"]
    for label, value, color in zip(labels, values, colors):
        axes[1].barh([0], [value / input_reads * 100], left=left, color=color, height=0.55, label=label)
        if value / input_reads >= 0.05:
            axes[1].text(left + value / input_reads * 50, 0, f"{value / input_reads:.1%}", ha="center", va="center", color="white", fontweight="bold")
        left += value / input_reads * 100
    axes[1].set_xlim(0, 100)
    axes[1].set_yticks([])
    axes[1].set_xlabel("Percent of input reads")
    axes[1].set_title("B  One denominator: 199,982 reads", loc="left", fontweight="bold")
    axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3)
    fig.suptitle("HUMAnN read flow is distinct from RPK and pathway special rows", y=1.04, fontsize=12, fontweight="bold")
    save_figure(fig, figure_dir, "19-read-flow")


def render_gene_stratification(rows: list[dict[str, Any]], figure_dir: Path) -> None:
    configure_plot_style()
    multistratified = [row for row in rows if int(row["Strata"]) > 1]
    selected = sorted(
        multistratified,
        key=lambda row: float(row["CommunityCPM"]),
        reverse=True,
    )[:12]
    if len(selected) < 12:
        selected_ids = {str(row["FeatureID"]) for row in selected}
        remaining = [row for row in rows if str(row["FeatureID"]) not in selected_ids]
        selected.extend(
            sorted(
                remaining,
                key=lambda row: float(row["CommunityCPM"]),
                reverse=True,
            )[: 12 - len(selected)]
        )
    selected.reverse()
    labels = [shorten(row["FeatureName"], 38) for row in selected]
    community = [float(row["CommunityCPM"]) for row in selected]
    dominant = [float(row["DominantFraction"]) * value for row, value in zip(selected, community)]
    remaining = [value - top for value, top in zip(community, dominant)]
    y = np.arange(len(selected))
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    ax.barh(y, dominant, color="#2A9D8F", label="Dominant stratum")
    ax.barh(y, remaining, left=dominant, color="#A8DADC", label="Other strata")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Copies per million (community denominator)")
    ax.set_title(
        "Multi-stratum gene-family abundance equals the sum of its strata",
        loc="left",
        fontweight="bold",
    )
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right")
    fig.subplots_adjust(left=0.38)
    save_figure(fig, figure_dir, "19-gene-family-stratification")


def render_pathway_contributions(rows: list[dict[str, str]], figure_dir: Path) -> None:
    configure_plot_style()
    selected = sorted(rows, key=lambda row: float(row["CommunityCPM"]), reverse=True)[:12]
    selected.reverse()
    labels = [shorten(row["PathwayName"], 40) for row in selected]
    community = np.array([float(row["CommunityCPM"]) for row in selected])
    strata = np.array([float(row["StrataCPMSum"]) for row in selected])
    y = np.arange(len(selected))
    fig, ax = plt.subplots(figsize=(9.6, 5.8))
    for index, (left, right) in enumerate(zip(community, strata)):
        ax.plot([left, right], [index, index], color="#B0B0B0", linewidth=1.2, zorder=1)
    ax.scatter(community, y, color="#E76F51", s=42, label="Community reconstruction", zorder=3)
    ax.scatter(strata, y, color="#264653", s=42, marker="D", label="Sum of taxon reconstructions", zorder=3)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Copies per million (same normalization denominator)")
    ax.set_title("Pathway abundance is reconstructed separately at community and taxon levels", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2)
    fig.subplots_adjust(left=0.39, bottom=0.18)
    save_figure(fig, figure_dir, "19-pathway-contributions")


def render_abundance_coverage(rows: list[dict[str, str]], figure_dir: Path) -> None:
    configure_plot_style()
    abundance = np.array([float(row["CommunityCPM"]) for row in rows])
    coverage = np.array([float(row["Coverage"]) for row in rows])
    dominance = np.array([float(row["DominantFractionOfStrata"]) for row in rows])
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    scatter = ax.scatter(abundance, coverage, c=dominance, cmap="viridis", s=34, alpha=0.8, edgecolor="white", linewidth=0.3)
    ax.set_xscale("log")
    ax.set_xlabel("Pathway abundance (copies per million, log scale)")
    ax.set_ylabel("Pathway coverage (0–1 evidence score)")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title("Abundance and coverage answer different questions", loc="left", fontweight="bold")
    ax.grid(color="#DDDDDD", linewidth=0.6)
    ax.set_axisbelow(True)
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("Dominant taxon share of stratified abundance")
    candidates = [int(np.argmax(abundance))]
    for lower, upper in ((0.0, 0.20), (0.40, 0.80)):
        eligible = np.where((coverage >= lower) & (coverage < upper))[0]
        if eligible.size:
            candidates.append(int(eligible[np.argmax(abundance[eligible])]))
    offsets = [(-8, -20), (10, -16), (10, 10)]
    alignments = ["right", "left", "left"]
    for index, offset, alignment in zip(candidates, offsets, alignments):
        ax.annotate(
            rows[index]["PathwayID"],
            (abundance[index], coverage[index]),
            xytext=offset,
            textcoords="offset points",
            ha=alignment,
            fontsize=7,
            arrowprops={"arrowstyle": "-", "color": "#666666", "linewidth": 0.6},
        )
    save_figure(fig, figure_dir, "19-abundance-coverage")


def validate_figures(figure_dir: Path, checks: Checks) -> None:
    for stem in FIGURE_STEMS:
        pdf = figure_dir / f"{stem}.pdf"
        png = figure_dir / f"{stem}.png"
        tiff = figure_dir / f"{stem}.tiff"
        checks.add("figure", f"{stem}-pdf", pdf.is_file() and pdf.stat().st_size > 8_000 and pdf.read_bytes()[:4] == b"%PDF", pdf.stat().st_size if pdf.exists() else "missing")
        checks.add("figure", f"{stem}-png", png.is_file() and png.stat().st_size > 15_000, png.stat().st_size if png.exists() else "missing")
        checks.add("figure", f"{stem}-tiff", tiff.is_file() and tiff.stat().st_size > 15_000, tiff.stat().st_size if tiff.exists() else "missing")
        if png.is_file():
            with Image.open(png) as image:
                checks.add("figure", f"{stem}-png-dimensions", image.width >= 1200 and image.height >= 600, f"{image.width}x{image.height}")
        if tiff.is_file():
            with Image.open(tiff) as image:
                dpi = image.info.get("dpi", (0, 0))
                compression = image.info.get("compression", "")
                checks.add("figure", f"{stem}-tiff-dpi", min(dpi) >= 349, dpi)
                checks.add("figure", f"{stem}-tiff-lzw", str(compression).lower() in {"tiff_lzw", "lzw", "5"}, compression)


def routine_validate(args: argparse.Namespace) -> None:
    if args.output_dir is None or args.figure_dir is None:
        raise SystemExit("Routine QA requires --output-dir and --figure-dir")
    root = args.project_root.resolve()
    prefix = args.environment_prefix.resolve()
    frozen = args.frozen_dir.resolve()
    output = args.output_dir.resolve()
    figures = args.figure_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    checksum_entries = verify_frozen_checksums(frozen, checks)
    versions = environment_checks(root, prefix, checks)
    frozen_compat = frozen / "metaphlan-version-compat"
    frozen_compat_sha256 = (
        hash_file(frozen_compat) if frozen_compat.is_file() else "missing"
    )
    checks.add(
        "environment",
        "frozen-metaphlan-version-wrapper",
        frozen_compat_sha256 == EXPECTED_COMPAT_SHA256,
        frozen_compat_sha256,
    )
    source_audit = source_checks(root, checks)
    database_audit = database_checks(root, frozen, checks)

    summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    checks.add("summary", "initialization-passed", summary["status"] == "passed", summary["status"])
    checks.add("summary", "input-pairs", summary["input_pairs"] == EXPECTED_INPUT["pairs"], summary["input_pairs"])
    checks.add("summary", "input-reads", summary["input_reads"] == EXPECTED_INPUT["reads"], summary["input_reads"])
    checks.add("summary", "unique-read-ids", summary["combined_unique_read_ids"] == EXPECTED_INPUT["reads"], summary["combined_unique_read_ids"])
    checks.add("summary", "vjan26-rejected", summary["vjan26_expected_rejection"] is True, summary["vjan26_direct_input_exit_status"])
    checks.add(
        "summary",
        "prescreen-not-bypassed",
        summary["humann_bypass_prescreen"] is False,
        summary["humann_bypass_prescreen"],
    )
    checks.add(
        "summary",
        "humann-selected-species",
        summary["humann_selected_species_after_alias_mapping"]
        == EXPECTED_HUMANN_SELECTED_SPECIES,
        summary["humann_selected_species_after_alias_mapping"],
    )
    checks.add(
        "summary",
        "custom-chocophlan-files",
        summary["custom_chocophlan_files"] == EXPECTED_CUSTOM_CHOCO_FILES,
        summary["custom_chocophlan_files"],
    )
    checks.add("summary", "threads", summary["threads"] == 8, summary["threads"])
    checks.add("summary", "pythonhashseed", summary["pythonhashseed"] == 0, summary["pythonhashseed"])

    prescreen = read_tsv(frozen / "prescreen-release-audit.tsv")
    vjun_selected = sum(row["SelectedByvJun23"] == "Yes" for row in prescreen)
    vjan_selected = sum(row["SelectedByvJan26"] == "Yes" for row in prescreen)
    intersection = sum(row["SelectionStatus"] == "Both" for row in prescreen)
    checks.add("prescreen", "vjun23-header", summary["prescreen"]["profile_header_has_vjun23"] is True, summary["metaphlan_database"])
    checks.add("prescreen", "vjun-selected-count", vjun_selected == summary["prescreen"]["vjun23_selected_species"], vjun_selected)
    checks.add("prescreen", "vjan-selected-count", vjan_selected == summary["prescreen"]["vjan26_selected_species"], vjan_selected)
    checks.add("prescreen", "intersection-count", intersection == summary["prescreen"]["selected_intersection"], intersection)
    checks.add("prescreen", "primary-not-vjan26", summary["metaphlan_database"] == "mpa_vJun23_CHOCOPhlAnSGB_202403", summary["metaphlan_database"])

    read_flow = read_tsv(frozen / "read-flow.tsv")
    flow = {row["Stage"]: int(row["Reads"]) for row in read_flow}
    checks.add("read-flow", "input-denominator", flow["Input"] == EXPECTED_INPUT["reads"], flow["Input"])
    checks.add("read-flow", "fate-conservation", flow["Nucleotide mapped"] + flow["Translated mapped"] + flow["Unmapped"] == flow["Input"], sum(flow[label] for label in ("Nucleotide mapped", "Translated mapped", "Unmapped")))
    checks.add("read-flow", "nonnegative", all(value >= 0 for value in flow.values()), flow)
    checks.add("read-flow", "summary-nucleotide", flow["Nucleotide mapped"] == summary["input_read_flow"]["nucleotide_mapped"], flow["Nucleotide mapped"])
    checks.add("read-flow", "summary-translated", flow["Translated mapped"] == summary["input_read_flow"]["translated_mapped"], flow["Translated mapped"])
    checks.add("read-flow", "summary-unmapped", flow["Unmapped"] == summary["input_read_flow"]["unmapped"], flow["Unmapped"])

    gene_rpk = read_humann_table(frozen / "genefamilies-rpk.tsv")
    gene_cpm = read_humann_table(frozen / "genefamilies-cpm.tsv")
    gene_relab = read_humann_table(frozen / "genefamilies-relab.tsv")
    reaction_rpk = read_humann_table(frozen / "reactions-rpk.tsv")
    reaction_cpm = read_humann_table(frozen / "reactions-cpm.tsv")
    reaction_relab = read_humann_table(frozen / "reactions-relab.tsv")
    path_rpk = read_humann_table(frozen / "pathabundance-rpk.tsv")
    path_cpm = read_humann_table(frozen / "pathabundance-cpm.tsv")
    path_relab = read_humann_table(frozen / "pathabundance-relab.tsv")
    coverage = read_humann_table(frozen / "pathcoverage.tsv")
    feature_audit = normalization_checks(
        gene_rpk,
        gene_cpm,
        gene_relab,
        reaction_rpk,
        reaction_cpm,
        reaction_relab,
        path_rpk,
        path_cpm,
        path_relab,
        coverage,
        checks,
    )
    genes, gene_metrics = gene_stratification_rows(gene_rpk, gene_cpm, gene_relab)
    reactions, reaction_metrics = gene_stratification_rows(
        reaction_rpk, reaction_cpm, reaction_relab
    )
    computed_regroup = regroup_audit_row(prefix, gene_rpk, reaction_rpk)
    frozen_regroup_rows = read_tsv(frozen / "regroup-audit.tsv")
    checks.add(
        "regroup",
        "single-audit-row",
        len(frozen_regroup_rows) == 1,
        len(frozen_regroup_rows),
    )
    if len(frozen_regroup_rows) != 1:
        raise SystemExit("Article 19 regroup audit must contain exactly one row")
    frozen_regroup = frozen_regroup_rows[0]
    pathways, pathway_metrics = pathway_contribution_rows(path_rpk, path_cpm, path_relab, coverage)
    checks.add("stratification", "gene-features", gene_metrics["features"] == summary["gene_family_features"], gene_metrics["features"])
    checks.add("stratification", "gene-additive", gene_metrics["violations"] == 0 == summary["gene_family_stratification_violations"], gene_metrics["violations"])
    checks.add(
        "stratification",
        "reaction-features",
        reaction_metrics["features"]
        == summary["reaction_features"]
        == EXPECTED_REACTION_FEATURES,
        reaction_metrics["features"],
    )
    checks.add(
        "stratification",
        "reaction-additive",
        reaction_metrics["violations"]
        == 0
        == summary["reaction_stratification_violations"],
        reaction_metrics["violations"],
    )
    checks.add(
        "regroup",
        "mapping-group",
        frozen_regroup["MappingGroup"]
        == summary["regroup_mapping_group"]
        == EXPECTED_REGROUP_MAPPING["group"],
        frozen_regroup["MappingGroup"],
    )
    checks.add(
        "regroup",
        "mapping-bytes",
        int(frozen_regroup["MappingBytes"])
        == summary["regroup_mapping_bytes"]
        == EXPECTED_REGROUP_MAPPING["bytes"],
        frozen_regroup["MappingBytes"],
    )
    checks.add(
        "regroup",
        "mapping-sha256",
        frozen_regroup["MappingSHA256"]
        == summary["regroup_mapping_sha256"]
        == EXPECTED_REGROUP_MAPPING["sha256"],
        frozen_regroup["MappingSHA256"],
    )
    checks.add(
        "regroup",
        "ungrouped-present",
        math.isclose(
            computed_regroup["UngroupedRPK"],
            EXPECTED_REGROUP_UNGROUPED_RPK,
            rel_tol=1e-12,
            abs_tol=1e-9,
        )
        and math.isclose(
            float(frozen_regroup["UngroupedRPK"]),
            computed_regroup["UngroupedRPK"],
            rel_tol=1e-12,
            abs_tol=1e-9,
        )
        and math.isclose(
            summary["regroup_ungrouped_rpk"],
            computed_regroup["UngroupedRPK"],
            rel_tol=1e-12,
            abs_tol=1e-9,
        ),
        computed_regroup["UngroupedRPK"],
    )
    checks.add(
        "regroup",
        "protected-unmapped-preserved",
        math.isclose(
            computed_regroup["SourceUnmappedRPK"],
            computed_regroup["ReactionUnmappedRPK"],
            rel_tol=0,
            abs_tol=1e-9,
        ),
        computed_regroup["ReactionUnmappedRPK"],
    )
    checks.add(
        "regroup",
        "mapped-one-to-many-expansion",
        computed_regroup["MappedExpansionFactor"] > 1
        and math.isclose(
            float(frozen_regroup["MappedExpansionFactor"]),
            computed_regroup["MappedExpansionFactor"],
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
        and math.isclose(
            summary["regroup_mapped_expansion_factor"],
            computed_regroup["MappedExpansionFactor"],
            rel_tol=1e-12,
            abs_tol=1e-12,
        ),
        computed_regroup["MappedExpansionFactor"],
    )
    checks.add(
        "regroup",
        "total-feature-space-expansion",
        computed_regroup["TotalRegroupExpansionFactor"] > 1
        and math.isclose(
            float(frozen_regroup["TotalRegroupExpansionFactor"]),
            computed_regroup["TotalRegroupExpansionFactor"],
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
        and math.isclose(
            summary["regroup_total_expansion_factor"],
            computed_regroup["TotalRegroupExpansionFactor"],
            rel_tol=1e-12,
            abs_tol=1e-12,
        ),
        computed_regroup["TotalRegroupExpansionFactor"],
    )
    checks.add("stratification", "pathway-features", pathway_metrics["features"] == summary["pathway_features"], pathway_metrics["features"])
    checks.add("stratification", "pathway-nonadditive", pathway_metrics["nonadditive_features"] == summary["pathway_nonadditive_features"] and pathway_metrics["nonadditive_features"] > 0, pathway_metrics["nonadditive_features"])
    checks.add("stratification", "coverage-range", pathway_metrics["coverage_out_of_range"] == 0 == summary["coverage_out_of_range"], pathway_metrics["coverage_out_of_range"])

    units = read_tsv(frozen / "units-contract.tsv")
    checks.add("units", "four-output-semantics", len(units) == 4, len(units))
    coverage_contract = next(row for row in units if row["Output"] == "Pathway coverage")
    checks.add("units", "coverage-never-renormalized", coverage_contract["RenormalizeCoverage"] == "No", coverage_contract["RenormalizeCoverage"])
    path_contract = next(row for row in units if row["Output"] == "Pathway abundance")
    checks.add("units", "pathway-not-additive", path_contract["CanSumStrataToCommunity"] == "No", path_contract["CanSumStrataToCommunity"])
    reaction_contract = next(
        row for row in units if row["Output"] == "Regrouped reactions"
    )
    checks.add(
        "units",
        "reaction-alternate-feature-space",
        "one-to-many" in reaction_contract["Interpretation"]
        and "UNGROUPED" in reaction_contract["SpecialRows"],
        reaction_contract["Interpretation"],
    )

    render_read_flow(read_flow, figures)
    render_gene_stratification(genes, figures)
    render_pathway_contributions(read_tsv(frozen / "top-pathways.tsv"), figures)
    render_abundance_coverage(read_tsv(frozen / "pathway-contributions.tsv"), figures)
    validate_figures(figures, checks)

    stratification_audit = read_tsv(frozen / "stratification-audit.tsv")
    write_tsv(output / "source-audit.tsv", source_audit, list(source_audit[0]))
    write_tsv(output / "database-audit.tsv", database_audit, list(database_audit[0]))
    write_tsv(output / "prescreen-audit.tsv", prescreen, list(prescreen[0]))
    write_tsv(output / "read-flow-audit.tsv", read_flow, list(read_flow[0]))
    write_tsv(output / "units-audit.tsv", units, list(units[0]))
    write_tsv(output / "stratification-audit.tsv", stratification_audit, list(stratification_audit[0]))
    write_tsv(
        output / "regroup-audit.tsv",
        frozen_regroup_rows,
        list(frozen_regroup_rows[0]),
    )
    write_tsv(output / "feature-audit.tsv", feature_audit, list(feature_audit[0]))
    write_tsv(output / "validation-audit.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])
    validation_summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "input_pairs": summary["input_pairs"],
        "input_reads": summary["input_reads"],
        "humann_version": versions["HUMAnN"],
        "metaphlan_version": versions["MetaPhlAn"],
        "metaphlan_database": summary["metaphlan_database"],
        "bowtie2_version": versions["Bowtie2"],
        "diamond_version": versions["DIAMOND"],
        "chocophlan_release": summary["chocophlan_release"],
        "uniref90_release": summary["uniref90_release"],
        "nucleotide_mapped_reads": flow["Nucleotide mapped"],
        "translated_mapped_reads": flow["Translated mapped"],
        "unmapped_reads": flow["Unmapped"],
        "vjun23_selected_species": vjun_selected,
        "vjan26_selected_species": vjan_selected,
        "humann_selected_species_after_alias_mapping": summary[
            "humann_selected_species_after_alias_mapping"
        ],
        "custom_chocophlan_files": summary["custom_chocophlan_files"],
        "humann_bypass_prescreen": summary["humann_bypass_prescreen"],
        "gene_family_features": gene_metrics["features"],
        "gene_family_stratification_violations": gene_metrics["violations"],
        "regroup_mapping_group": summary["regroup_mapping_group"],
        "regroup_mapping_file": summary["regroup_mapping_file"],
        "regroup_mapping_bytes": summary["regroup_mapping_bytes"],
        "regroup_mapping_sha256": summary["regroup_mapping_sha256"],
        "reaction_features": reaction_metrics["features"],
        "reaction_stratification_violations": reaction_metrics["violations"],
        "regroup_ungrouped_rpk": computed_regroup["UngroupedRPK"],
        "regroup_mapped_expansion_factor": computed_regroup[
            "MappedExpansionFactor"
        ],
        "regroup_total_expansion_factor": computed_regroup[
            "TotalRegroupExpansionFactor"
        ],
        "pathway_features": pathway_metrics["features"],
        "pathway_nonadditive_features": pathway_metrics["nonadditive_features"],
        "coverage_out_of_range": pathway_metrics["coverage_out_of_range"],
        "database_checksum_entries": summary["database_checksum_entries"],
        "frozen_checksum_entries": checksum_entries,
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "qa_network_access": False,
    }
    write_json(output / "validation-summary.json", validation_summary)
    (output / "validation.log").write_text(
        f"Article 19 routine QA: {validation_summary['status']}\n"
        f"Checks: {checks.passed} passed, {checks.failed} failed\n"
        f"Frozen payloads: {checksum_entries}\n"
        f"Read flow: nucleotide={flow['Nucleotide mapped']}; "
        f"translated={flow['Translated mapped']}; unmapped={flow['Unmapped']}\n"
        f"Features: genes={gene_metrics['features']}; "
        f"reactions={reaction_metrics['features']}; "
        f"pathways={pathway_metrics['features']}\n"
        f"Regroup expansion: mapped={computed_regroup['MappedExpansionFactor']:.6f}; "
        f"total={computed_regroup['TotalRegroupExpansionFactor']:.6f}\n",
        encoding="utf-8",
    )
    if checks.failed:
        raise SystemExit(f"Article 19 routine QA failed {checks.failed} checks")
    print(f"PASS Article 19 routine QA: {checks.passed} checks; {checksum_entries} frozen payloads")


def main() -> None:
    args = parse_args()
    if args.initialize_frozen:
        initialize_frozen(args)
    else:
        routine_validate(args)


if __name__ == "__main__":
    main()
