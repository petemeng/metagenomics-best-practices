#!/usr/bin/env python3
"""Build and validate Article 18 profiler-comparison evidence.

Initialization is the only mode that reads the large MetaPhlAn metadata pickle,
the mOTUs database, the mOTUs archive, FASTQ-derived work products, or resource
logs. Routine QA is network-free and database-free: it verifies the frozen
checksum manifest, recomputes all reported metrics from small tables, and
renders the four evidence figures.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import gc
import gzip
import hashlib
import json
import math
import os
import pickle
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article18-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


EXPECTED_ENV_SHA256 = {
    "env/motus.yml": "1c7e446580595e73c58cb19b8d91af07bad025e5cc5df0c0a99b15f510ca2ed3",
    "env/motus-linux-64.lock": "5a7f959c52919d228ee90b6d4c225e18e22c326aaabd1fc23db2338b4d9092d9",
}
EXPECTED_LOCK_PACKAGES = 97
EXPECTED_LOCK_MARKERS = {
    "motus": "motus-4.1.0-",
    "bwa": "bwa-0.7.19-",
    "vsearch": "vsearch-2.31.0-",
    "python": "/python-3.12.13-",
    "matplotlib-base": "matplotlib-base-3.10.5-",
    "numpy": "numpy-2.5.1-",
    "pillow": "pillow-12.3.0-",
}
EXPECTED_INPUT = {
    "pairs": 99_991,
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
EXPECTED_MOTUS_ARCHIVE = {
    "bytes": 5_552_983_255,
    "md5": "471ea128f0c0839f5c4629b949ea5f8a",
    "sha256": "7d2d6382ecf766b23ef362311715cd612243af82c00c936da4597afb4e4df375",
}
EXPECTED_MOTUS_INSTALLED_BYTES = 9_699_381_148
DETECTION_THRESHOLD = 0.0001  # 0.01% of each resolved native profile.
FIGURE_STEMS = (
    "18-feature-space-crosswalk",
    "18-recovery-nontruth",
    "18-composition-agreement",
    "18-resource-footprint",
)
PROFILE_ORDER = (
    "metaphlan-default",
    "kraken-bracken-native",
    "kraken-bracken-control-aware",
    "motus-g1",
    "motus-g3-default",
    "motus-g6",
)
PRIMARY_PROFILES = (
    "metaphlan-default",
    "kraken-bracken-native",
    "motus-g3-default",
)
TEXT_SUFFIXES = {".json", ".log", ".md5", ".sha256", ".sh", ".tsv", ".txt", ".yml"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--environment-prefix", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--motus-database-root", type=Path)
    parser.add_argument("--motus-database-archive", type=Path)
    parser.add_argument("--metaphlan-database-pkl", type=Path)
    parser.add_argument("--initialize-frozen", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    return parser.parse_args()


def hash_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_file_multi(path: Path) -> dict[str, Any]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return {
        "bytes": path.stat().st_size,
        "md5": md5.hexdigest(),
        "sha256": sha256.hexdigest(),
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def as_float(value: Any) -> float:
    if value in (None, "", "nan", "NaN", "NA"):
        return math.nan
    return float(value)


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


def environment_checks(project_root: Path, prefix: Path, checks: Checks) -> dict[str, str]:
    for relative, expected in EXPECTED_ENV_SHA256.items():
        path = project_root / relative
        observed = hash_file(path) if path.is_file() else "missing"
        checks.add("environment", f"sha256-{path.name}", observed == expected, observed)

    lock_lines = (project_root / "env/motus-linux-64.lock").read_text(encoding="utf-8").splitlines()
    package_lines = [line for line in lock_lines if line.startswith("http")]
    checks.add("environment", "lock-package-count", len(package_lines) == EXPECTED_LOCK_PACKAGES, len(package_lines))
    for package, marker in EXPECTED_LOCK_MARKERS.items():
        matches = [line for line in package_lines if marker in line]
        checks.add("environment", f"lock-{package}", len(matches) == 1, len(matches))

    env = os.environ.copy()
    env["PATH"] = f"{prefix / 'bin'}:{env.get('PATH', '')}"
    env["PYTHONNOUSERSITE"] = "1"
    motus_result = subprocess.run(
        [str(prefix / "bin/motus"), "--version"],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    motus_text = motus_result.stdout + motus_result.stderr
    bwa_result = subprocess.run([str(prefix / "bin/bwa")], capture_output=True, text=True, env=env, check=False)
    vsearch_result = subprocess.run([str(prefix / "bin/vsearch"), "--version"], capture_output=True, text=True, env=env, check=False)
    python_result = subprocess.run([str(prefix / "bin/python"), "--version"], capture_output=True, text=True, env=env, check=True)
    versions = {
        "mOTUs": (re.search(r"(?:Version:|mOTUs4:)\s*([0-9.]+)", motus_text) or [None, "missing"])[1],
        "BWA": (re.search(r"Version:\s*([^\s]+)", bwa_result.stderr + bwa_result.stdout) or [None, "missing"])[1],
        "VSEARCH": (re.search(r"vsearch\s+v?([^_\s]+)", vsearch_result.stderr + vsearch_result.stdout, re.I) or [None, "missing"])[1],
        "Python": (re.search(r"Python\s+([0-9.]+)", python_result.stdout + python_result.stderr) or [None, "missing"])[1],
    }
    expected_versions = {
        "mOTUs": "4.1.0",
        "BWA": "0.7.19-r1273",
        "VSEARCH": "2.31.0",
        "Python": "3.12.13",
    }
    for tool, expected in expected_versions.items():
        checks.add("environment", f"version-{tool.lower()}", versions[tool] == expected, versions[tool])
    return versions


def static_source_checks(project_root: Path, checks: Checks) -> tuple[list[dict[str, str]], dict[str, Any]]:
    source = read_tsv(project_root / "data/small/18-source-manifest.tsv")
    checks.add("source", "two-mate-rows", len(source) == 2, len(source))
    checks.add("source", "one-run", {row["RunAccession"] for row in source} == {"ERR9765746"}, sorted({row["RunAccession"] for row in source}))
    checks.add("source", "two-mates", {row["Mate"] for row in source} == {"R1", "R2"}, sorted({row["Mate"] for row in source}))
    for row in source:
        mate = row["Mate"]
        expected = EXPECTED_INPUT[mate]
        checks.add("source", f"{mate.lower()}-bytes", int(row["CompressedBytes"]) == expected["bytes"], row["CompressedBytes"])
        checks.add("source", f"{mate.lower()}-sha256", row["CompressedSHA256"] == expected["sha256"], row["CompressedSHA256"])
        checks.add("source", f"{mate.lower()}-records", int(row["Records"]) == EXPECTED_INPUT["pairs"], row["Records"])
        checks.add("source", f"{mate.lower()}-bases", int(row["Bases"]) == expected["bases"], row["Bases"])
        checks.add("source", f"{mate.lower()}-pair-hash", row["PairIDHash"] == EXPECTED_INPUT["pair_hash"], row["PairIDHash"])

    previous: dict[str, Any] = {}
    for article, relative in {
        "qc": "data/small/13-qc-frozen/run-summary.json",
        "metaphlan": "data/small/15-metaphlan-frozen/run-summary.json",
        "kraken_native": "data/small/16-kraken-bracken-frozen/run-summary.json",
        "kraken_control": "data/small/17-kraken-database-confidence-frozen/run-summary.json",
    }.items():
        summary = json.loads((project_root / relative).read_text(encoding="utf-8"))
        previous[article] = summary
        checks.add("source", f"{article}-status", summary["status"] == "passed", summary["status"])
        pair_value = (
            summary.get("retained_pairs")
            if article == "qc"
            else summary.get("input_pairs", summary.get("mock_input_pairs"))
        )
        checks.add("source", f"{article}-input-pairs", int(pair_value) == EXPECTED_INPUT["pairs"], pair_value)
    checks.add("source", "metaphlan-version", previous["metaphlan"]["metaphlan_version"] == "4.2.5", previous["metaphlan"]["metaphlan_version"])
    checks.add("source", "kraken-version", previous["kraken_native"]["kraken2_version"] == "2.17.1", previous["kraken_native"]["kraken2_version"])
    checks.add("source", "bracken-package", previous["kraken_native"]["bracken_package_version"] == "3.1p1", previous["kraken_native"]["bracken_package_version"])
    return source, previous


def truth_species_table(project_root: Path, checks: Checks) -> tuple[list[dict[str, Any]], dict[int, dict[str, Any]]]:
    rows = [
        row
        for row in read_tsv(project_root / "data/small/17-kraken-database-confidence-frozen/reference-crosswalk.tsv")
        if row["DatabaseID"] == "standard16"
    ]
    checks.add("truth", "assembly-rows", len(rows) == 71, len(rows))
    grouped: dict[int, dict[str, Any]] = {}
    for row in rows:
        taxid = int(row["SpeciesTaxID"])
        if taxid not in grouped:
            grouped[taxid] = {
                "taxid": taxid,
                "name": row["DatabaseSpeciesName"],
                "mass": 0.0,
                "assemblies": [],
                "truth_indices": [],
                "expected_gtdb_names": set(),
                "kraken_represented": False,
                "kraken_reference_statuses": set(),
            }
        item = grouped[taxid]
        item["mass"] += float(row["ExpectedGenomePercent"])
        item["assemblies"].append(row["AssemblyAccession"])
        item["truth_indices"].append(int(row["TruthIndex"]))
        item["expected_gtdb_names"].add(row["ExpectedOrganismGTDBRS207"])
        item["kraken_reference_statuses"].add(row["ReferenceStatus"])
        item["kraken_represented"] = item["kraken_represented"] or row["ReferenceStatus"] != "No same-species reference"
    ordered = [grouped[key] for key in sorted(grouped)]
    checks.add("truth", "current-species", len(ordered) == 69, len(ordered))
    total_mass = math.fsum(item["mass"] for item in ordered)
    checks.add("truth", "printed-mass-preserved", abs(total_mass - 100.01215654459679) < 1e-10, repr(total_mass))
    return ordered, grouped


def species_taxid_from_lineage(lineage: str, taxid_path: str) -> int | None:
    taxa = lineage.split("|")
    taxids = taxid_path.split("|")
    for index, taxon in enumerate(taxa):
        if taxon.startswith("s__") and index < len(taxids) and taxids[index].isdigit():
            return int(taxids[index])
    return None


def sgb_from_lineage(lineage: str) -> str | None:
    for taxon in reversed(lineage.split("|")):
        if taxon.startswith("t__SGB"):
            return taxon
    return None


def build_metaphlan_mapping(
    pkl_path: Path,
    truth_taxids: set[int],
    checks: Checks,
) -> tuple[dict[int, set[str]], dict[str, set[int]], dict[str, Any]]:
    with bz2.open(pkl_path, "rb") as handle:
        database = pickle.load(handle)
    marker_count = len(database["markers"])
    taxonomy = database["taxonomy"]
    merged = database["merged_taxon"]
    taxonomy_count = len(taxonomy)
    merged_count = len(merged)
    database.pop("markers")
    gc.collect()

    species_to_sgbs: defaultdict[int, set[str]] = defaultdict(set)
    sgb_to_species: defaultdict[str, set[int]] = defaultdict(set)
    for lineage, value in taxonomy.items():
        taxid_path = value[0]
        taxid = species_taxid_from_lineage(lineage, taxid_path)
        sgb = sgb_from_lineage(lineage)
        if taxid in truth_taxids and sgb:
            species_to_sgbs[taxid].add(sgb)
            sgb_to_species[sgb].add(taxid)

    merged_truth_links = 0
    for key, alternatives in merged.items():
        target_lineage = key[0] if isinstance(key, tuple) else str(key)
        sgb = sgb_from_lineage(target_lineage)
        if not sgb:
            continue
        for alternative in alternatives:
            if not isinstance(alternative, (tuple, list)) or len(alternative) < 2:
                continue
            taxid = species_taxid_from_lineage(str(alternative[0]), str(alternative[1]))
            if taxid in truth_taxids:
                species_to_sgbs[taxid].add(sgb)
                sgb_to_species[sgb].add(taxid)
                merged_truth_links += 1

    checks.add("metaphlan-database", "taxonomy-entries", taxonomy_count == 72_000, taxonomy_count)
    checks.add("metaphlan-database", "marker-entries", marker_count == 13_907_686, marker_count)
    checks.add("metaphlan-database", "truth-mapping-present", bool(species_to_sgbs), len(species_to_sgbs))
    return dict(species_to_sgbs), dict(sgb_to_species), {
        "taxonomy_entries": taxonomy_count,
        "marker_entries": marker_count,
        "merged_taxon_entries": merged_count,
        "merged_truth_links": merged_truth_links,
    }


def normalize_accession_from_motus(genome_id: str) -> str | None:
    match = re.search(r"(GC[AF])-(\d+)-V(\d+)", genome_id)
    if not match:
        return None
    return f"{match.group(1)}_{match.group(2)}.{int(match.group(3))}"


def normalized_species_name(value: str) -> str:
    value = value.strip()
    if ";s__" in value:
        value = value.rsplit(";s__", 1)[1]
    value = re.sub(r"^[a-z]__", "", value)
    value = value.replace("_", " ")
    value = re.sub(r"\s+", " ", value).strip().lower()
    return value


def read_motus_database_header(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        values[key.strip()] = value.strip()
    return values


def build_motus_mapping(
    project_root: Path,
    database_root: Path,
    truth: dict[int, dict[str, Any]],
    checks: Checks,
) -> tuple[dict[int, set[str]], dict[str, set[int]], dict[int, set[str]], dict[str, str], dict[str, Any]]:
    db_dir = database_root / "db_mOTU"
    taxonomy_path = db_dir / "mOTUsv4.1.gtdb.taxonomy.80mv.tsv.gz"
    taxonomy: dict[str, str] = {}
    taxonomy_versions: set[str] = set()
    with gzip.open(taxonomy_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            taxonomy[row["MOTU"]] = row["GTDB"]
            taxonomy_versions.add(row["GTDB_VERSION"])
    checks.add("motus-database", "taxonomy-r226", taxonomy_versions == {"R226"}, sorted(taxonomy_versions))

    truth_rows = read_tsv(project_root / "data/small/17-mock1-truth.tsv")
    accession_to_species: dict[str, int] = {}
    accession_to_assembly: dict[str, str] = {}
    truth_index_to_species: dict[str, int] = {}
    standard_rows = [
        row
        for row in read_tsv(project_root / "data/small/17-kraken-database-confidence-frozen/reference-crosswalk.tsv")
        if row["DatabaseID"] == "standard16"
    ]
    for row in standard_rows:
        truth_index_to_species[row["TruthIndex"]] = int(row["SpeciesTaxID"])
    for row in truth_rows:
        species = truth_index_to_species[row["TruthIndex"]]
        for field in ("AssemblyAccession", "CurrentAccession", "PairedRefSeqAccession"):
            accession = row.get(field, "")
            if accession:
                accession_to_species[accession] = species
                accession_to_assembly[accession] = row["AssemblyAccession"]

    species_to_motus: defaultdict[int, set[str]] = defaultdict(set)
    motu_to_species: defaultdict[str, set[int]] = defaultdict(set)
    species_exact_assemblies: defaultdict[int, set[str]] = defaultdict(set)
    genome_rows = 0
    exact_rows = 0
    genomes_path = db_dir / "mOTUsv4.1.genomes.tsv.gz"
    with gzip.open(genomes_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genome_rows += 1
            accession = normalize_accession_from_motus(row["GENOME"])
            if accession not in accession_to_species:
                continue
            species = accession_to_species[accession]
            motu = row["MOTU4"]
            species_to_motus[species].add(motu)
            motu_to_species[motu].add(species)
            species_exact_assemblies[species].add(accession_to_assembly[accession])
            exact_rows += 1

    expected_names: defaultdict[str, set[int]] = defaultdict(set)
    for species, item in truth.items():
        for name in item["expected_gtdb_names"]:
            expected_names[normalized_species_name(name)].add(species)
    name_candidates: defaultdict[int, set[str]] = defaultdict(set)
    for motu, lineage in taxonomy.items():
        name = normalized_species_name(lineage)
        if name.startswith("unknown ") or name in {"", "na"}:
            continue
        species_set = expected_names.get(name, set())
        if len(species_set) == 1:
            name_candidates[next(iter(species_set))].add(motu)

    header = read_motus_database_header(db_dir / "mOTUsv4.1.db")
    checks.add("motus-database", "database-version", header.get("database_version") == "4.1", header.get("database_version"))
    checks.add("motus-database", "database-date", header.get("database_date") == "2026-01-26", header.get("database_date"))
    checks.add(
        "motus-database",
        "genome-row-count-discrepancy-recorded",
        genome_rows + 1 == int(header["genomes"]),
        f"data_rows={genome_rows};version_header={header['genomes']}",
    )
    checks.add("motus-database", "exact-truth-assemblies", sum(len(values) for values in species_exact_assemblies.values()) == 66, sum(len(values) for values in species_exact_assemblies.values()))
    checks.add("motus-database", "runtime-count-contract", int(header["assigned_motus"]) + int(header["unassigned_motus"]) == 124_300, int(header["assigned_motus"]) + int(header["unassigned_motus"]))
    return (
        dict(species_to_motus),
        dict(motu_to_species),
        dict(name_candidates),
        taxonomy,
        {
            "header": header,
            "taxonomy_rows": len(taxonomy),
            "taxonomy_versions": sorted(taxonomy_versions),
            "genome_rows": genome_rows,
            "exact_matching_rows": exact_rows,
            "exact_assembly_species": len(species_exact_assemblies),
            "exact_assemblies": sum(len(values) for values in species_exact_assemblies.values()),
            "species_exact_assemblies": {str(key): sorted(value) for key, value in species_exact_assemblies.items()},
        },
    )


def relation_status(features: set[str], reverse: dict[str, set[int]]) -> str:
    if not features:
        return "Absent"
    split = len(features) > 1
    merged = any(len(reverse.get(feature, set())) > 1 for feature in features)
    if split and merged:
        return "Split + merge"
    if split:
        return "Split"
    if merged:
        return "Merge"
    return "One-to-one"


def build_crosswalk(
    truth_ordered: list[dict[str, Any]],
    meta_species_to_features: dict[int, set[str]],
    meta_feature_to_species: dict[str, set[int]],
    motus_species_to_features: dict[int, set[str]],
    motus_feature_to_species: dict[str, set[int]],
    motus_name_candidates: dict[int, set[str]],
    motus_summary: dict[str, Any],
    checks: Checks,
) -> tuple[list[dict[str, Any]], set[int]]:
    rows: list[dict[str, Any]] = []
    exact_assemblies = motus_summary["species_exact_assemblies"]
    for item in truth_ordered:
        taxid = item["taxid"]
        meta_features = set(meta_species_to_features.get(taxid, set()))
        motus_features = set(motus_species_to_features.get(taxid, set()))
        meta_status = relation_status(meta_features, meta_feature_to_species)
        motus_status = relation_status(motus_features, motus_feature_to_species)
        kraken_status = "One-to-one" if item["kraken_represented"] else "Absent"
        common = meta_status == kraken_status == motus_status == "One-to-one"
        rows.append(
            {
                "TruthSpeciesTaxID": taxid,
                "TruthSpeciesName": item["name"],
                "ExpectedGenomePercent": item["mass"],
                "TruthAssemblies": len(item["assemblies"]),
                "AssemblyAccessions": ";".join(sorted(item["assemblies"])),
                "MetaPhlAnStatus": meta_status,
                "MetaPhlAnFeatureCount": len(meta_features),
                "MetaPhlAnFeatures": ";".join(sorted(meta_features)),
                "KrakenBrackenStatus": kraken_status,
                "KrakenReferenceEvidence": ";".join(sorted(item["kraken_reference_statuses"])),
                "mOTUsStatus": motus_status,
                "mOTUsFeatureCount": len(motus_features),
                "mOTUsFeatures": ";".join(sorted(motus_features)),
                "mOTUsExactAssemblies": len(exact_assemblies.get(str(taxid), [])),
                "mOTUsExactAssemblyAccessions": ";".join(exact_assemblies.get(str(taxid), [])),
                "mOTUsExactNameCandidates": ";".join(sorted(motus_name_candidates.get(taxid, set()))),
                "NameOnlyEstablishesReference": "No",
                "CommonOneToOne": "Yes" if common else "No",
            }
        )
    common_domain = {int(row["TruthSpeciesTaxID"]) for row in rows if row["CommonOneToOne"] == "Yes"}
    checks.add("crosswalk", "one-row-per-truth-species", len(rows) == 69, len(rows))
    checks.add("crosswalk", "strict-common-domain-nonempty", len(common_domain) >= 20, len(common_domain))
    checks.add("crosswalk", "name-only-never-reference", all(row["NameOnlyEstablishesReference"] == "No" for row in rows), "all No")
    return rows, common_domain


def parse_motus_profile(path: Path) -> tuple[dict[str, float], str, dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3 or not lines[0].startswith("#tool_version="):
        raise ValueError(f"Unexpected mOTUs profile: {path}")
    metadata: dict[str, str] = {}
    for field in lines[0][1:].split("\t"):
        key, value = field.split("=", 1)
        metadata[key] = value
    reader = csv.DictReader(lines[1:], delimiter="\t")
    sample_field = [field for field in reader.fieldnames or [] if field not in {"mOTU", "Taxonomy"}]
    if len(sample_field) != 1:
        raise ValueError(f"Unexpected mOTUs sample columns: {sample_field}")
    values: dict[str, float] = {}
    unresolved = ""
    for row in reader:
        value = float(row[sample_field[0]])
        values[row["mOTU"]] = value
        if row["mOTU"].endswith("_unassigned"):
            unresolved = row["mOTU"]
    return values, unresolved, metadata


def profile_definitions(
    project_root: Path,
    frozen_dir: Path,
    meta_feature_to_species: dict[str, set[int]],
    motus_feature_to_species: dict[str, set[int]],
    motus_name_candidates: dict[int, set[str]],
    motus_taxonomy: dict[str, str],
    truth: dict[int, dict[str, Any]],
    checks: Checks,
) -> list[dict[str, Any]]:
    profiles: list[dict[str, Any]] = []

    meta_rows = read_tsv(project_root / "data/small/15-metaphlan-frozen/sgb-profile.tsv")
    meta_values = {row["SGBLabel"]: float(row["RelativeAbundancePct"]) for row in meta_rows}
    profiles.append(
        {
            "id": "metaphlan-default",
            "tool": "MetaPhlAn 4",
            "role": "Native default",
            "feature_space": "SGB (species aggregate retained)",
            "database": "mpa_vJan26_CHOCOPhlAnSGB_202605",
            "parameters": "default profiling; --nproc 8",
            "unit": "relative abundance (%)",
            "values": meta_values,
            "unresolved": "",
            "associations": {feature: set(meta_feature_to_species.get(feature, set())) for feature in meta_values},
            "reference_species": set(meta_feature_to_species_species(meta_feature_to_species)),
        }
    )

    kraken_reference = {
        taxid for taxid, item in truth.items() if item["kraken_represented"]
    }
    kraken_specs = (
        (
            "kraken-bracken-native",
            "Native confidence 0",
            "Standard-8-20260626",
            "Kraken confidence 0; hit groups 2; Bracken S/r150/t10",
            project_root / "data/small/16-kraken-bracken-frozen/bracken-species-r150-t10.tsv",
        ),
        (
            "kraken-bracken-control-aware",
            "Control-aware sensitivity",
            "Standard-16-20260626",
            "Kraken confidence 0.10; hit groups 2; Bracken S/r150/t10",
            project_root / "data/small/17-kraken-database-confidence-frozen/bracken/standard16-mock-c010-h2.tsv",
        ),
    )
    for profile_id, role, database, parameters, path in kraken_specs:
        rows = read_tsv(path)
        values = {row["taxonomy_id"]: float(row["fraction_total_reads"]) for row in rows}
        profiles.append(
            {
                "id": profile_id,
                "tool": "Kraken2 + Bracken",
                "role": role,
                "feature_space": "NCBI species taxon",
                "database": database,
                "parameters": parameters,
                "unit": "fraction of Bracken species estimates",
                "values": values,
                "unresolved": "",
                "associations": {
                    feature: ({int(feature)} if int(feature) in truth else set())
                    for feature in values
                },
                "reference_species": kraken_reference,
            }
        )

    exact_name_feature_to_species: defaultdict[str, set[int]] = defaultdict(set)
    for species, candidates in motus_name_candidates.items():
        for feature in candidates:
            exact_name_feature_to_species[feature].add(species)
    strict_motus_reference = set()
    for feature, species_set in motus_feature_to_species.items():
        strict_motus_reference.update(species_set)
    for minimum_markers in (1, 3, 6):
        values, unresolved, metadata = parse_motus_profile(frozen_dir / f"motus-profile-g{minimum_markers}.tsv")
        checks.add("profile", f"motus-g{minimum_markers}-version", metadata["tool_version"] == "4.1.0", metadata["tool_version"])
        checks.add("profile", f"motus-g{minimum_markers}-database", metadata["database_version"] == "4.1", metadata["database_version"])
        checks.add("profile", f"motus-g{minimum_markers}-threshold", int(metadata["min_mgcs"]) == minimum_markers, metadata["min_mgcs"])
        associations: dict[str, set[int]] = {}
        for feature in values:
            strict = set(motus_feature_to_species.get(feature, set()))
            if strict:
                associations[feature] = strict
            else:
                # Exact cross-release GTDB species-name identity may annotate a
                # detected feature, but it never expands the strict denominator.
                associations[feature] = set(exact_name_feature_to_species.get(feature, set()))
        profiles.append(
            {
                "id": f"motus-g{minimum_markers}" if minimum_markers != 3 else "motus-g3-default",
                "tool": "mOTUs 4",
                "role": "Default" if minimum_markers == 3 else ("Recall sensitivity" if minimum_markers == 1 else "Precision sensitivity"),
                "feature_space": "mOTU",
                "database": "mOTUs marker-gene database 4.1 / GTDB R226",
                "parameters": f"75 bp alignment; g={minimum_markers}; INSERT_SCALED",
                "unit": "INSERT_SCALED counts",
                "values": values,
                "unresolved": unresolved,
                "associations": associations,
                "reference_species": strict_motus_reference,
                "taxonomy": motus_taxonomy,
            }
        )

    observed_ids = [profile["id"] for profile in profiles]
    checks.add("profile", "six-profile-branches", tuple(observed_ids) == PROFILE_ORDER, observed_ids)
    return profiles


def meta_feature_to_species_species(mapping: dict[str, set[int]]) -> set[int]:
    output: set[int] = set()
    for species_set in mapping.values():
        output.update(species_set)
    return output


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2.0 + 1.0
        start = end
    return ranks


def spearman_rho(first: np.ndarray, second: np.ndarray) -> float:
    first_ranks = rankdata(first)
    second_ranks = rankdata(second)
    if np.std(first_ranks) == 0 or np.std(second_ranks) == 0:
        return math.nan
    return float(np.corrcoef(first_ranks, second_ranks)[0, 1])


def total_variation(first: np.ndarray, second: np.ndarray) -> float:
    return float(0.5 * np.sum(np.abs(first - second)))


def summarize_profiles(
    profiles: list[dict[str, Any]],
    truth: dict[int, dict[str, Any]],
    common_domain: set[int],
    checks: Checks,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    performance: list[dict[str, Any]] = []
    harmonized: list[dict[str, Any]] = []
    vectors: dict[str, np.ndarray] = {}
    common_order = sorted(common_domain)
    truth_values = np.array([truth[taxid]["mass"] for taxid in common_order], dtype=float)
    truth_values /= truth_values.sum()
    vectors["truth"] = truth_values
    common_mass = math.fsum(truth[taxid]["mass"] for taxid in common_order)

    for profile in profiles:
        values = profile["values"]
        unresolved = profile["unresolved"]
        unresolved_mass = float(values.get(unresolved, 0.0)) if unresolved else 0.0
        resolved = {feature: float(value) for feature, value in values.items() if feature != unresolved and float(value) > 0}
        resolved_total = math.fsum(resolved.values())
        native_total = math.fsum(float(value) for value in values.values())
        normalized = {feature: value / resolved_total for feature, value in resolved.items()} if resolved_total > 0 else {}
        positive_features = set(resolved)
        threshold_features = {feature for feature, value in normalized.items() if value >= DETECTION_THRESHOLD}
        positive_truth: set[int] = set()
        threshold_truth: set[int] = set()
        for feature in positive_features:
            positive_truth.update(profile["associations"].get(feature, set()))
        for feature in threshold_features:
            threshold_truth.update(profile["associations"].get(feature, set()))
        reference_species = set(profile["reference_species"])
        recovered_raw = positive_truth & reference_species
        recovered_threshold = threshold_truth & reference_species
        nontruth_raw = {feature for feature in positive_features if not profile["associations"].get(feature, set())}
        nontruth_threshold = {feature for feature in threshold_features if not profile["associations"].get(feature, set())}
        ambiguous = {feature for feature in positive_features if len(profile["associations"].get(feature, set())) > 1}
        nontruth_mass = math.fsum(resolved[feature] for feature in nontruth_raw)

        species_raw: defaultdict[int, float] = defaultdict(float)
        species_features: defaultdict[int, list[str]] = defaultdict(list)
        for feature, value in resolved.items():
            associated = profile["associations"].get(feature, set())
            if len(associated) != 1:
                continue
            species = next(iter(associated))
            if species in common_domain:
                species_raw[species] += value
                species_features[species].append(feature)
        common_total = math.fsum(species_raw.values())
        if common_total <= 0:
            raise ValueError(f"No mass in strict common domain for {profile['id']}")
        vector = np.array([species_raw[taxid] / common_total for taxid in common_order], dtype=float)
        vectors[profile["id"]] = vector
        for index, taxid in enumerate(common_order):
            harmonized.append(
                {
                    "ProfileID": profile["id"],
                    "Tool": profile["tool"],
                    "BranchRole": profile["role"],
                    "TruthSpeciesTaxID": taxid,
                    "TruthSpeciesName": truth[taxid]["name"],
                    "ExpectedGenomeFraction": truth_values[index],
                    "EstimatedFraction": vector[index],
                    "NativeFeatureIDs": ";".join(sorted(species_features.get(taxid, []))),
                    "NativeDetectedAt0.01Pct": "Yes" if any(normalized.get(feature, 0.0) >= DETECTION_THRESHOLD for feature in species_features.get(taxid, [])) else "No",
                }
            )

        denom = len(reference_species)
        performance.append(
            {
                "ProfileID": profile["id"],
                "Tool": profile["tool"],
                "BranchRole": profile["role"],
                "FeatureSpace": profile["feature_space"],
                "Database": profile["database"],
                "Parameters": profile["parameters"],
                "NativeUnit": profile["unit"],
                "ReferenceTruthSpecies": denom,
                "ReferenceCoveredExpectedPercent": math.fsum(truth[taxid]["mass"] for taxid in reference_species),
                "NativePositiveFeatures": len(positive_features),
                "NativePositiveTruthSpecies": len(recovered_raw),
                "ThresholdedTruthSpecies": len(recovered_threshold),
                "ReferenceAwareRecoveryRaw": len(recovered_raw) / denom if denom else math.nan,
                "ReferenceAwareRecoveryAt0.01Pct": len(recovered_threshold) / denom if denom else math.nan,
                "NonTruthFeaturesRaw": len(nontruth_raw),
                "NonTruthFeaturesAt0.01Pct": len(nontruth_threshold),
                "NonTruthNativeMassPct": 100.0 * nontruth_mass / native_total if native_total else math.nan,
                "NonTruthResolvedMassPct": 100.0 * nontruth_mass / resolved_total if resolved_total else math.nan,
                "AmbiguousTruthFeatures": len(ambiguous),
                "UnresolvedFeature": unresolved,
                "UnresolvedNativeMass": unresolved_mass,
                "UnresolvedNativeMassPct": 100.0 * unresolved_mass / native_total if native_total else 0.0,
                "NativeMassTotal": native_total,
                "ResolvedMassTotal": resolved_total,
                "CommonDomainSpecies": len(common_domain),
                "CommonDomainExpectedPercent": common_mass,
                "CommonDomainNativeMassPct": 100.0 * common_total / native_total if native_total else math.nan,
                "CommonDomainSpearmanVsTruth": spearman_rho(vector, truth_values),
                "CommonDomainTVVsTruth": total_variation(vector, truth_values),
            }
        )

    pairwise: list[dict[str, Any]] = []
    all_ids = ["truth", *PROFILE_ORDER]
    for left_index, left in enumerate(all_ids):
        for right in all_ids[left_index + 1 :]:
            pairwise.append(
                {
                    "ProfileA": left,
                    "ProfileB": right,
                    "CommonDomainSpecies": len(common_domain),
                    "SpearmanRho": spearman_rho(vectors[left], vectors[right]),
                    "TotalVariation": total_variation(vectors[left], vectors[right]),
                }
            )
    checks.add("performance", "six-performance-rows", len(performance) == 6, len(performance))
    checks.add("performance", "harmonized-row-count", len(harmonized) == 6 * len(common_domain), len(harmonized))
    checks.add("performance", "pairwise-row-count", len(pairwise) == math.comb(7, 2), len(pairwise))
    checks.add("performance", "closed-harmonized-profiles", all(abs(math.fsum(float(row["EstimatedFraction"]) for row in harmonized if row["ProfileID"] == profile_id) - 1.0) < 1e-12 for profile_id in PROFILE_ORDER), "all six")
    return performance, harmonized, pairwise


def parse_elapsed(value: str) -> float:
    parts = value.strip().split(":")
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    return float(parts[0])


def parse_gnu_time(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8")
    elapsed_match = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([^\n]+)", text)
    rss_match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    cpu_match = re.search(r"Percent of CPU this job got:\s*([0-9.]+)%", text)
    exit_match = re.search(r"Exit status:\s*(\d+)", text)
    if not elapsed_match or not rss_match or not exit_match:
        raise ValueError(f"Incomplete GNU time log: {path}")
    return {
        "wall_seconds": parse_elapsed(elapsed_match.group(1)),
        "max_rss_mib": int(rss_match.group(1)) / 1024.0,
        "cpu_percent": float(cpu_match.group(1)) if cpu_match else math.nan,
        "exit_status": int(exit_match.group(1)),
    }


def build_resource_table(
    project_root: Path,
    work_dir: Path,
    previous: dict[str, Any],
    motus_installed_bytes: int,
    checks: Checks,
) -> list[dict[str, Any]]:
    article15 = previous["metaphlan"]
    article16 = previous["kraken_native"]
    article17_resources = read_tsv(project_root / "data/small/17-kraken-database-confidence-frozen/resource-usage.tsv")
    control_rows = [
        row for row in article17_resources
        if row["BranchID"] == "standard16-mock-c010-h2" and row["Tool"] in {"Kraken2", "Bracken"}
    ]
    if len(control_rows) != 2:
        raise ValueError("Missing Article 17 control-aware resource rows")
    standard16_bytes = next(
        int(row["InstalledBytes"])
        for row in previous["kraken_control"]["database_audit"]
        if row["DatabaseID"] == "standard16"
    )

    motus_stages = {
        "map_tax": parse_gnu_time(work_dir / "resources/motus-map-tax.resources.txt"),
        "calc_mgc": parse_gnu_time(work_dir / "resources/motus-calc-mgc.resources.txt"),
        "calc_motu_g1": parse_gnu_time(work_dir / "resources/motus-calc-motu-g1.resources.txt"),
        "calc_motu_g3": parse_gnu_time(work_dir / "resources/motus-calc-motu-g3.resources.txt"),
        "calc_motu_g6": parse_gnu_time(work_dir / "resources/motus-calc-motu-g6.resources.txt"),
    }
    checks.add("resource", "motus-stage-exit-status", all(stage["exit_status"] == 0 for stage in motus_stages.values()), [stage["exit_status"] for stage in motus_stages.values()])

    rows = [
        {
            "ProfileID": "metaphlan-default",
            "Tool": "MetaPhlAn 4",
            "PipelineStages": "MetaPhlAn full profile",
            "WallSeconds": parse_elapsed(article15["resource_summary"]["full"]["elapsed"]),
            "MaxRSSMiB": int(article15["resource_summary"]["full"]["maximum_rss_kb"]) / 1024.0,
            "Threads": 8,
            "InstalledDatabaseBytes": int(article15["database_files"]["total_bytes"]),
            "CacheState": "local filesystem; warm/cold cache not randomized",
            "ComparisonBoundary": "full marker database; read-level estimated counts are not Kraken classifications",
        },
        {
            "ProfileID": "kraken-bracken-native",
            "Tool": "Kraken2 + Bracken",
            "PipelineStages": "Kraken2 Standard-8 + Bracken species",
            "WallSeconds": parse_elapsed(article16["resource_summary"]["kraken2_full"]["elapsed"]) + parse_elapsed(article16["resource_summary"]["bracken_species_r150_t10"]["elapsed"]),
            "MaxRSSMiB": max(int(article16["resource_summary"]["kraken2_full"]["maximum_rss_kb"]) / 1024.0, int(article16["resource_summary"]["bracken_species_r150_t10"]["maximum_rss_kb"]) / 1024.0),
            "Threads": 8,
            "InstalledDatabaseBytes": int(article16["database_files"]["total_bytes"]),
            "CacheState": "local filesystem; warm/cold cache not randomized",
            "ComparisonBoundary": "8 GB capped Standard index; Bracken is serial post-processing",
        },
        {
            "ProfileID": "kraken-bracken-control-aware",
            "Tool": "Kraken2 + Bracken",
            "PipelineStages": "Kraken2 Standard-16 c0.10 + Bracken species",
            "WallSeconds": math.fsum(float(row["WallSeconds"]) for row in control_rows),
            "MaxRSSMiB": max(float(row["MaxRSSMiB"]) for row in control_rows),
            "Threads": 8,
            "InstalledDatabaseBytes": standard16_bytes,
            "CacheState": "local filesystem; warm/cold cache not randomized",
            "ComparisonBoundary": "16 GB capped Standard index; sensitivity branch, not native default",
        },
        {
            "ProfileID": "motus-g3-default",
            "Tool": "mOTUs 4",
            "PipelineStages": "map_tax + calc_mgc + calc_motu g3",
            "WallSeconds": motus_stages["map_tax"]["wall_seconds"] + motus_stages["calc_mgc"]["wall_seconds"] + motus_stages["calc_motu_g3"]["wall_seconds"],
            "MaxRSSMiB": max(motus_stages[key]["max_rss_mib"] for key in ("map_tax", "calc_mgc", "calc_motu_g3")),
            "Threads": 8,
            "InstalledDatabaseBytes": motus_installed_bytes,
            "CacheState": "local filesystem; warm/cold cache not randomized",
            "ComparisonBoundary": "full marker-gene index; one map reused across g=1/3/6",
        },
    ]
    checks.add("resource", "four-comparable-pipelines", len(rows) == 4, len(rows))
    return rows


def normalize_frozen_paths(
    frozen_dir: Path,
    replacements: dict[str, str],
    checks: Checks,
) -> dict[str, Any]:
    ordered = sorted(
        ((key, value) for key, value in replacements.items() if key and not key.startswith("__")),
        key=lambda item: len(item[0]),
        reverse=True,
    )
    changed = 0
    count = 0
    for path in sorted(frozen_dir.rglob("*")):
        if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8")
        updated = text
        for source, target in ordered:
            occurrences = updated.count(source)
            if occurrences:
                updated = updated.replace(source, target)
                count += occurrences
        if updated != text:
            path.write_text(updated, encoding="utf-8")
            changed += 1
    residual = []
    host_root = replacements.get("__HOST_ROOT__", "")
    if host_root:
        for path in sorted(frozen_dir.rglob("*")):
            if path.is_file() and path.suffix in TEXT_SUFFIXES and host_root in path.read_text(encoding="utf-8"):
                residual.append(str(path.relative_to(frozen_dir)))
    checks.add("path", "no-host-specific-paths", not residual, residual)
    return {"changed_files": changed, "replacement_count": count, "remaining_hits": len(residual)}


def write_checksum_manifest(frozen_dir: Path) -> int:
    target = frozen_dir / "file-checksums.sha256"
    rows = []
    for path in sorted(frozen_dir.rglob("*")):
        if not path.is_file() or path == target:
            continue
        rows.append(
            {
                "SHA256": hash_file(path),
                "Bytes": path.stat().st_size,
                "File": str(path.relative_to(frozen_dir)),
            }
        )
    write_tsv(target, rows, ["SHA256", "Bytes", "File"])
    return len(rows)


def verify_checksum_manifest(frozen_dir: Path, checks: Checks) -> int:
    rows = read_tsv(frozen_dir / "file-checksums.sha256")
    failures: list[str] = []
    for row in rows:
        path = frozen_dir / row["File"]
        if not path.is_file():
            failures.append(f"missing:{row['File']}")
            continue
        if path.stat().st_size != int(row["Bytes"]):
            failures.append(f"bytes:{row['File']}")
            continue
        if hash_file(path) != row["SHA256"]:
            failures.append(f"sha256:{row['File']}")
    expected_files = {
        str(path.relative_to(frozen_dir))
        for path in frozen_dir.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    }
    manifest_files = {row["File"] for row in rows}
    if expected_files != manifest_files:
        failures.append(f"coverage:missing={sorted(expected_files - manifest_files)};extra={sorted(manifest_files - expected_files)}")
    checks.add("checksum", "frozen-payload", not failures, failures if failures else len(rows))
    return len(rows)


def audit_motus_database_files(database_root: Path, checks: Checks) -> tuple[list[dict[str, Any]], int]:
    db_dir = database_root / "db_mOTU"
    header_text = (db_dir / "mOTUsv4.1.db").read_text(encoding="utf-8")
    expected_md5 = {
        match.group(1): match.group(2)
        for match in re.finditer(r"^([^:\n]+):\s*([0-9a-f]{32})$", header_text, re.M)
    }
    rows: list[dict[str, Any]] = []
    installed_bytes = 0
    md5_failures: list[str] = []
    for path in sorted(db_dir.iterdir()):
        if not path.is_file():
            continue
        summary = hash_file_multi(path)
        installed_bytes += summary["bytes"]
        official = expected_md5.get(path.name, "")
        status = "NOT_LISTED_IN_VERSION_FILE"
        if official:
            status = "PASS" if summary["md5"] == official else "FAIL"
            if status == "FAIL":
                md5_failures.append(path.name)
        rows.append(
            {
                "SHA256": summary["sha256"],
                "Bytes": summary["bytes"],
                "MD5": summary["md5"],
                "VersionFileMD5": official,
                "VersionFileMD5Status": status,
                "File": path.name,
            }
        )
    checks.add("motus-database", "internal-md5", not md5_failures, md5_failures if md5_failures else len(expected_md5))
    checks.add("motus-database", "installed-bytes", installed_bytes == EXPECTED_MOTUS_INSTALLED_BYTES, installed_bytes)
    return rows, installed_bytes


def profiler_design_rows() -> list[dict[str, Any]]:
    return [
        {
            "ProfileID": "metaphlan-default",
            "Tool": "MetaPhlAn 4.2.5",
            "BranchRole": "Native default",
            "Database": "mpa_vJan26_CHOCOPhlAnSGB_202605",
            "NativeFeature": "SGB; species aggregate",
            "NativeUnit": "relative abundance percent",
            "PrimaryParameters": "default; nproc=8",
            "InputUnit": "199,929 reads >=70 bp from 99,991 pairs",
            "ReferenceEvidence": "NCBI species taxid in compressed metadata pickle",
            "UnresolvedSemantics": "terminal-clade estimated reads; not a fragment classification rate",
        },
        {
            "ProfileID": "kraken-bracken-native",
            "Tool": "Kraken2 2.17.1 + Bracken 3.1p1",
            "BranchRole": "Native confidence 0",
            "Database": "Standard-8-20260626",
            "NativeFeature": "NCBI species taxon",
            "NativeUnit": "Bracken estimated paired fragments / species-model total",
            "PrimaryParameters": "confidence=0; hit-groups=2; S/r150/t10",
            "InputUnit": "99,991 paired fragments",
            "ReferenceEvidence": "sequence accession and NCBI species ancestor",
            "UnresolvedSemantics": "Kraken unclassified paired fragments; outside Bracken table",
        },
        {
            "ProfileID": "kraken-bracken-control-aware",
            "Tool": "Kraken2 2.17.1 + Bracken 3.1p1",
            "BranchRole": "Control-aware sensitivity",
            "Database": "Standard-16-20260626",
            "NativeFeature": "NCBI species taxon",
            "NativeUnit": "Bracken estimated paired fragments / species-model total",
            "PrimaryParameters": "confidence=0.10; hit-groups=2; S/r150/t10",
            "InputUnit": "99,991 paired fragments",
            "ReferenceEvidence": "sequence accession and NCBI species ancestor",
            "UnresolvedSemantics": "Kraken unclassified/ancestor calls; outside Bracken table",
        },
        *[
            {
                "ProfileID": f"motus-g{g}" if g != 3 else "motus-g3-default",
                "Tool": "mOTUs 4.1.0",
                "BranchRole": "Default" if g == 3 else ("Recall sensitivity" if g == 1 else "Precision sensitivity"),
                "Database": "marker-gene database 4.1 / GTDB R226",
                "NativeFeature": "mOTU",
                "NativeUnit": "INSERT_SCALED marker-gene insert counts",
                "PrimaryParameters": f"alignment-length=75; g={g}; INSERT_SCALED",
                "InputUnit": "199,982 reads / 99,991 pairs; insert evidence after mapping",
                "ReferenceEvidence": "exact GCA/GCF assembly accession encoded in genome catalog",
                "UnresolvedSemantics": "mOTU unassigned row; marker-aligned inserts only",
            }
            for g in (1, 3, 6)
        ],
    ]


def initialize_frozen(args: argparse.Namespace) -> None:
    required = {
        "work-dir": args.work_dir,
        "motus-database-root": args.motus_database_root,
        "motus-database-archive": args.motus_database_archive,
        "metaphlan-database-pkl": args.metaphlan_database_pkl,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise SystemExit(f"Initialization requires: {', '.join(missing)}")

    project_root = args.project_root.resolve()
    prefix = args.environment_prefix.resolve()
    frozen_dir = args.frozen_dir.resolve()
    work_dir = args.work_dir.resolve()
    database_root = args.motus_database_root.resolve()
    archive = args.motus_database_archive.resolve()
    metaphlan_pkl = args.metaphlan_database_pkl.resolve()
    frozen_dir.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    versions = environment_checks(project_root, prefix, checks)
    source_manifest, previous = static_source_checks(project_root, checks)
    truth_ordered, truth = truth_species_table(project_root, checks)

    database_manifest = read_tsv(project_root / "data/small/18-database-manifest.tsv")
    checks.add("database", "one-motus-manifest-row", len(database_manifest) == 1, len(database_manifest))
    manifest_row = database_manifest[0]
    checks.add("database", "immutable-zenodo-record", manifest_row["archive_url"] == "https://zenodo.org/records/20322482/files/db_mOTU.tar.gz", manifest_row["archive_url"])
    checks.add("database", "manifest-taxonomy-r226", "GTDB-R226" in manifest_row["notes"], manifest_row["notes"])
    archive_summary = hash_file_multi(archive)
    for key, expected in EXPECTED_MOTUS_ARCHIVE.items():
        checks.add("motus-database", f"archive-{key}", archive_summary[key] == expected, archive_summary[key])

    db_file_rows, installed_bytes = audit_motus_database_files(database_root, checks)
    write_tsv(
        frozen_dir / "database-files.sha256",
        db_file_rows,
        ["SHA256", "Bytes", "MD5", "VersionFileMD5", "VersionFileMD5Status", "File"],
    )

    meta_species_to_features, meta_feature_to_species, meta_summary = build_metaphlan_mapping(
        metaphlan_pkl, set(truth), checks
    )
    (
        motus_species_to_features,
        motus_feature_to_species,
        motus_name_candidates,
        motus_taxonomy,
        motus_summary,
    ) = build_motus_mapping(project_root, database_root, truth, checks)
    crosswalk, common_domain = build_crosswalk(
        truth_ordered,
        meta_species_to_features,
        meta_feature_to_species,
        motus_species_to_features,
        motus_feature_to_species,
        motus_name_candidates,
        motus_summary,
        checks,
    )
    profiles = profile_definitions(
        project_root,
        frozen_dir,
        meta_feature_to_species,
        motus_feature_to_species,
        motus_name_candidates,
        motus_taxonomy,
        truth,
        checks,
    )
    performance, harmonized, pairwise = summarize_profiles(profiles, truth, common_domain, checks)
    resources = build_resource_table(project_root, work_dir, previous, installed_bytes, checks)

    map_log = (work_dir / "logs/motus-map-tax.log").read_text(encoding="utf-8")
    mgc_log = (work_dir / "logs/motus-calc-mgc.log").read_text(encoding="utf-8")
    aligned_matches = re.findall(r"Total reads:\s*(\d+), Total aligned reads\s*(\d+),\s*([0-9.]+)%(?: aligned)?", map_log)
    insert_match = re.search(r"Read\s+(\d+) aligned inserts of which\s+([0-9.]+)% are multimappers", mgc_log)
    checks.add("motus-run", "combined-alignment-ledger", bool(aligned_matches) and tuple(map(int, aligned_matches[-1][:2])) == (199_982, 1_122), aligned_matches[-1] if aligned_matches else "missing")
    checks.add("motus-run", "aligned-inserts", bool(insert_match) and int(insert_match.group(1)) == 673, insert_match.group(1) if insert_match else "missing")
    aligned_reads = int(aligned_matches[-1][1])
    aligned_inserts = int(insert_match.group(1))
    multimapper_pct = float(insert_match.group(2))

    crosswalk_fields = [
        "TruthSpeciesTaxID", "TruthSpeciesName", "ExpectedGenomePercent",
        "TruthAssemblies", "AssemblyAccessions", "MetaPhlAnStatus",
        "MetaPhlAnFeatureCount", "MetaPhlAnFeatures", "KrakenBrackenStatus",
        "KrakenReferenceEvidence", "mOTUsStatus", "mOTUsFeatureCount",
        "mOTUsFeatures", "mOTUsExactAssemblies", "mOTUsExactAssemblyAccessions",
        "mOTUsExactNameCandidates", "NameOnlyEstablishesReference", "CommonOneToOne",
    ]
    write_tsv(frozen_dir / "truth-feature-crosswalk.tsv", crosswalk, crosswalk_fields)

    performance_fields = list(performance[0])
    write_tsv(frozen_dir / "benchmark-performance.tsv", performance, performance_fields)
    native_fields = [
        "ProfileID", "Tool", "BranchRole", "FeatureSpace", "Database", "Parameters",
        "NativeUnit", "NativePositiveFeatures", "NativeMassTotal", "ResolvedMassTotal",
        "UnresolvedFeature", "UnresolvedNativeMass", "UnresolvedNativeMassPct",
        "NonTruthFeaturesRaw", "NonTruthNativeMassPct",
    ]
    write_tsv(
        frozen_dir / "native-profile-summary.tsv",
        [{field: row[field] for field in native_fields} for row in performance],
        native_fields,
    )
    write_tsv(frozen_dir / "harmonized-profile.tsv", harmonized, list(harmonized[0]))
    write_tsv(frozen_dir / "pairwise-agreement.tsv", pairwise, list(pairwise[0]))
    write_tsv(frozen_dir / "resource-usage.tsv", resources, list(resources[0]))
    design = profiler_design_rows()
    write_tsv(frozen_dir / "profiler-design.tsv", design, list(design[0]))

    source_audit = [
        {
            "EvidenceID": f"ERR9765746-{row['Mate']}",
            "Source": row["SourceStage"],
            "Observed": f"{row['Records']} records; {row['CompressedSHA256']}",
            "Expected": f"99991 records; {EXPECTED_INPUT[row['Mate']]['sha256']}",
            "Status": "PASS",
            "Boundary": "same checksum-locked Article 13 clean FASTQ; FASTQ not committed",
        }
        for row in source_manifest
    ]
    source_audit.extend(
        [
            {
                "EvidenceID": "normalized-pair-id-hash",
                "Source": "Article 13 audit",
                "Observed": EXPECTED_INPUT["pair_hash"],
                "Expected": EXPECTED_INPUT["pair_hash"],
                "Status": "PASS",
                "Boundary": "99,991 synchronized pair identifiers",
            },
            {
                "EvidenceID": "truth-definition",
                "Source": "Meslier Table S3 + Article 17 NCBI snapshot",
                "Observed": "71 assemblies; 69 current NCBI species",
                "Expected": "71 assemblies; 69 current NCBI species",
                "Status": "PASS",
                "Boundary": "printed expected percentages preserved without silent renormalization",
            },
        ]
    )
    write_tsv(frozen_dir / "source-audit.tsv", source_audit, list(source_audit[0]))

    database_audit = [
        {
            "Tool": "MetaPhlAn 4.2.5",
            "DatabaseID": previous["metaphlan"]["database_release"],
            "ReleaseDate": "2026-05",
            "Taxonomy": "MetaPhlAn NCBI-taxid lineages + SGB",
            "ArchiveBytes": sum(item["bytes"] for item in previous["metaphlan"]["archive_audits"].values()),
            "ArchiveMD5": ";".join(item["md5"] for item in previous["metaphlan"]["archive_audits"].values()),
            "ArchiveSHA256": ";".join(item["sha256"] for item in previous["metaphlan"]["archive_audits"].values()),
            "InstalledBytes": previous["metaphlan"]["database_files"]["total_bytes"],
            "CatalogFeatures": meta_summary["taxonomy_entries"],
            "RuntimeCatalogFeatures": meta_summary["taxonomy_entries"],
            "ValidationStatus": "REUSED_CHECKSUM_LOCKED_ARTICLE15",
            "Boundary": "full SGB marker database; taxonomy pickle read once for crosswalk",
        },
        {
            "Tool": "Kraken2 2.17.1 + Bracken 3.1p1",
            "DatabaseID": "Standard-8-20260626",
            "ReleaseDate": "2026-06-26",
            "Taxonomy": "NCBI taxonomy",
            "ArchiveBytes": previous["kraken_native"]["database_archive_bytes"],
            "ArchiveMD5": previous["kraken_native"]["database_archive_md5"],
            "ArchiveSHA256": previous["kraken_native"]["database_archive_sha256"],
            "InstalledBytes": previous["kraken_native"]["database_files"]["total_bytes"],
            "CatalogFeatures": previous["kraken_native"]["database_files"]["file_count"],
            "RuntimeCatalogFeatures": previous["kraken_native"]["kraken_species_rows"],
            "ValidationStatus": "REUSED_CHECKSUM_LOCKED_ARTICLE16",
            "Boundary": "8 GB capped Standard minimizer table; native branch",
        },
        {
            "Tool": "Kraken2 2.17.1 + Bracken 3.1p1",
            "DatabaseID": "Standard-16-20260626",
            "ReleaseDate": "2026-06-26",
            "Taxonomy": "NCBI taxonomy",
            "ArchiveBytes": previous["kraken_control"]["database_archives"]["standard16"]["bytes"],
            "ArchiveMD5": previous["kraken_control"]["database_archives"]["standard16"]["md5"],
            "ArchiveSHA256": previous["kraken_control"]["database_archives"]["standard16"]["sha256"],
            "InstalledBytes": next(row["InstalledBytes"] for row in previous["kraken_control"]["database_audit"] if row["DatabaseID"] == "standard16"),
            "CatalogFeatures": next(row["ReferenceTaxIDs"] for row in previous["kraken_control"]["database_audit"] if row["DatabaseID"] == "standard16"),
            "RuntimeCatalogFeatures": "not a native-feature catalog count",
            "ValidationStatus": "REUSED_CHECKSUM_LOCKED_ARTICLE17",
            "Boundary": "16 GB capped Standard minimizer table; control-aware sensitivity branch",
        },
        {
            "Tool": "mOTUs 4.1.0",
            "DatabaseID": "mOTUs-marker-gene-database-4.1",
            "ReleaseDate": motus_summary["header"]["database_date"],
            "Taxonomy": "GTDB R226",
            "ArchiveBytes": archive_summary["bytes"],
            "ArchiveMD5": archive_summary["md5"],
            "ArchiveSHA256": archive_summary["sha256"],
            "InstalledBytes": installed_bytes,
            "CatalogFeatures": motus_summary["header"]["motus"],
            "RuntimeCatalogFeatures": int(motus_summary["header"]["assigned_motus"]) + int(motus_summary["header"]["unassigned_motus"]),
            "ValidationStatus": "VERIFIED; UPSTREAM_HEADER_COUNT_DISCREPANCY_RECORDED",
            "Boundary": "version field says 124,295 mOTUs while assigned + unassigned and runtime loader say 124,300",
        },
    ]
    write_tsv(frozen_dir / "database-audit.tsv", database_audit, list(database_audit[0]))

    path_audit = normalize_frozen_paths(
        frozen_dir,
        {
            str(metaphlan_pkl): "${METAPHLAN_DB_PKL}",
            str(archive): "${MOTUS_DB_ARCHIVE}",
            str(database_root): "${MOTUS_DB_ROOT}",
            str(work_dir): "${WORK_DIR}",
            str(prefix): "${MOTUS_ENV_PREFIX}",
            str(project_root): "${PROJECT_ROOT}",
            str(Path.home()): "${HOME}",
            "__HOST_ROOT__": str(Path.home()),
        },
        checks,
    )

    checks.add("output", "crosswalk-written", (frozen_dir / "truth-feature-crosswalk.tsv").is_file(), len(crosswalk))
    checks.add("output", "performance-written", (frozen_dir / "benchmark-performance.tsv").is_file(), len(performance))
    checks.add("output", "resources-written", (frozen_dir / "resource-usage.tsv").is_file(), len(resources))
    write_tsv(frozen_dir / "initialization-audit.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])
    if checks.failed:
        raise SystemExit(f"Article 18 initialization failed {checks.failed} checks")

    status_counts = {
        tool: {
            status: sum(row[field] == status for row in crosswalk)
            for status in ("One-to-one", "Split", "Merge", "Split + merge", "Absent")
        }
        for tool, field in {
            "MetaPhlAn": "MetaPhlAnStatus",
            "KrakenBracken": "KrakenBrackenStatus",
            "mOTUs": "mOTUsStatus",
        }.items()
    }
    perf_by_id = {row["ProfileID"]: row for row in performance}
    summary = {
        "status": "passed",
        "evidence_date": "2026-07-22",
        "source_project": "PRJEB52977",
        "source_sample": "SAMEA14435832",
        "run_accession": "ERR9765746",
        "input_pairs": EXPECTED_INPUT["pairs"],
        "input_reads": 199_982,
        "pair_id_sha256": EXPECTED_INPUT["pair_hash"],
        "truth_assemblies": 71,
        "truth_species": 69,
        "truth_expected_percent": math.fsum(item["mass"] for item in truth_ordered),
        "detection_threshold_fraction": DETECTION_THRESHOLD,
        "detection_threshold_percent": 100 * DETECTION_THRESHOLD,
        "common_one_to_one_species": len(common_domain),
        "common_one_to_one_expected_percent": math.fsum(truth[taxid]["mass"] for taxid in common_domain),
        "crosswalk_status_counts": status_counts,
        "metaphlan_version": "4.2.5",
        "metaphlan_database": previous["metaphlan"]["database_release"],
        "kraken2_version": "2.17.1",
        "bracken_cli_version": "3.0.1",
        "bracken_package_version": "3.1p1",
        "kraken_native_database": "Standard-8-20260626",
        "kraken_control_database": "Standard-16-20260626",
        "motus_version": versions["mOTUs"],
        "motus_database_version": motus_summary["header"]["database_version"],
        "motus_database_date": motus_summary["header"]["database_date"],
        "motus_taxonomy_release": "R226",
        "motus_archive": archive_summary,
        "motus_installed_bytes": installed_bytes,
        "motus_database_header_motus": int(motus_summary["header"]["motus"]),
        "motus_runtime_catalog_motus": int(motus_summary["header"]["assigned_motus"]) + int(motus_summary["header"]["unassigned_motus"]),
        "motus_exact_truth_assemblies": motus_summary["exact_assemblies"],
        "motus_exact_truth_species": motus_summary["exact_assembly_species"],
        "motus_aligned_reads": aligned_reads,
        "motus_aligned_read_fraction": aligned_reads / 199_982,
        "motus_aligned_inserts": aligned_inserts,
        "motus_multimapper_percent": multimapper_pct,
        "profile_headlines": {
            profile_id: {
                "reference_truth_species": int(perf_by_id[profile_id]["ReferenceTruthSpecies"]),
                "recovered_at_0.01pct": int(perf_by_id[profile_id]["ThresholdedTruthSpecies"]),
                "nontruth_features_at_0.01pct": int(perf_by_id[profile_id]["NonTruthFeaturesAt0.01Pct"]),
                "nontruth_native_mass_pct": float(perf_by_id[profile_id]["NonTruthNativeMassPct"]),
                "common_domain_spearman_vs_truth": float(perf_by_id[profile_id]["CommonDomainSpearmanVsTruth"]),
                "common_domain_tv_vs_truth": float(perf_by_id[profile_id]["CommonDomainTVVsTruth"]),
            }
            for profile_id in PROFILE_ORDER
        },
        "resource_rows": resources,
        "metaphlan_catalog": meta_summary,
        "motus_catalog": {key: value for key, value in motus_summary.items() if key != "species_exact_assemblies"},
        "path_normalization": path_audit,
        "raw_fastq_committed": False,
        "large_database_committed": False,
        "bam_committed": False,
        "mgc_committed": False,
        "one_time_network_access": True,
        "qa_network_access": False,
        "initialization_checks_passed": checks.passed,
        "initialization_checks_failed": checks.failed,
        "frozen_checksum_entries": 0,
    }
    payload_before_summary = sum(1 for path in frozen_dir.rglob("*") if path.is_file() and path.name != "file-checksums.sha256")
    summary["frozen_checksum_entries"] = payload_before_summary + (0 if (frozen_dir / "run-summary.json").exists() else 1)
    write_json(frozen_dir / "run-summary.json", summary)
    entries = write_checksum_manifest(frozen_dir)
    if entries != summary["frozen_checksum_entries"]:
        raise SystemExit(f"Checksum entry count changed unexpectedly: {entries} != {summary['frozen_checksum_entries']}")
    print(f"PASS Article 18 initialization: {checks.passed} checks; {entries} frozen payloads")


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
    fig.savefig(figure_dir / f"{stem}.png", dpi=200, bbox_inches="tight")
    fig.savefig(
        figure_dir / f"{stem}.tiff",
        dpi=350,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def render_crosswalk_figure(crosswalk: list[dict[str, str]], figure_dir: Path) -> None:
    configure_plot_style()
    tools = [
        ("MetaPhlAn", "MetaPhlAnStatus"),
        ("Kraken/Bracken", "KrakenBrackenStatus"),
        ("mOTUs", "mOTUsStatus"),
    ]
    categories = ["One-to-one", "Split", "Merge", "Split + merge", "Absent"]
    colors = {
        "One-to-one": "#2A9D8F",
        "Split": "#E9C46A",
        "Merge": "#F4A261",
        "Split + merge": "#E76F51",
        "Absent": "#7A7A7A",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0), gridspec_kw={"wspace": 0.28})
    x = np.arange(len(tools))
    bottoms_count = np.zeros(len(tools))
    bottoms_mass = np.zeros(len(tools))
    for category in categories:
        counts = np.array([sum(row[field] == category for row in crosswalk) for _, field in tools], dtype=float)
        masses = np.array(
            [math.fsum(float(row["ExpectedGenomePercent"]) for row in crosswalk if row[field] == category) for _, field in tools],
            dtype=float,
        )
        axes[0].bar(x, counts, bottom=bottoms_count, color=colors[category], width=0.68, label=category)
        axes[1].bar(x, masses, bottom=bottoms_mass, color=colors[category], width=0.68, label=category)
        for index, value in enumerate(counts):
            if value >= 3:
                axes[0].text(index, bottoms_count[index] + value / 2, f"{int(value)}", ha="center", va="center", fontsize=8)
        bottoms_count += counts
        bottoms_mass += masses
    axes[0].set_title("A  Truth species by feature relation", loc="left", fontweight="bold")
    axes[0].set_ylabel("Truth species (n)")
    axes[0].set_ylim(0, 72)
    axes[1].set_title("B  Expected mass by feature relation", loc="left", fontweight="bold")
    axes[1].set_ylabel("Expected genome percentage")
    axes[1].axhline(100.01215654459679, color="#333333", linewidth=0.8, linestyle="--")
    axes[1].set_ylim(0, 104)
    for axis in axes:
        axis.set_xticks(x, [label for label, _ in tools])
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.6, zorder=0)
        axis.set_axisbelow(True)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.03))
    fig.suptitle("Database-aware crosswalk: NCBI species, SGBs and mOTUs", y=1.02, fontsize=12, fontweight="bold")
    save_figure(fig, figure_dir, "18-feature-space-crosswalk")


def profile_display(profile_id: str) -> str:
    return {
        "metaphlan-default": "MetaPhlAn\ndefault",
        "kraken-bracken-native": "Kraken/Bracken\nc0",
        "kraken-bracken-control-aware": "Kraken/Bracken\nc0.10",
        "motus-g1": "mOTUs\ng1",
        "motus-g3-default": "mOTUs\ng3",
        "motus-g6": "mOTUs\ng6",
    }[profile_id]


def profile_color(profile_id: str) -> str:
    if profile_id.startswith("metaphlan"):
        return "#3A86FF"
    if profile_id == "kraken-bracken-native":
        return "#E76F51"
    if profile_id == "kraken-bracken-control-aware":
        return "#F4A261"
    return {"motus-g1": "#8AC926", "motus-g3-default": "#2A9D8F", "motus-g6": "#146B62"}[profile_id]


def render_recovery_figure(performance: list[dict[str, str]], figure_dir: Path) -> None:
    configure_plot_style()
    by_id = {row["ProfileID"]: row for row in performance}
    labels = [profile_display(profile_id).replace("\n", " ") for profile_id in PROFILE_ORDER]
    colors = [profile_color(profile_id) for profile_id in PROFILE_ORDER]
    recovery = [100 * float(by_id[profile_id]["ReferenceAwareRecoveryAt0.01Pct"]) for profile_id in PROFILE_ORDER]
    nontruth_count = [int(by_id[profile_id]["NonTruthFeaturesAt0.01Pct"]) for profile_id in PROFILE_ORDER]
    nontruth_mass = [float(by_id[profile_id]["NonTruthNativeMassPct"]) for profile_id in PROFILE_ORDER]
    denominators = [int(by_id[profile_id]["ReferenceTruthSpecies"]) for profile_id in PROFILE_ORDER]
    y = np.arange(len(PROFILE_ORDER))
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 5.0), gridspec_kw={"wspace": 0.28})
    panels = [
        (recovery, "A  Reference-aware recovery", "Recovered reference-present species (%)"),
        (nontruth_count, "B  Non-truth features", "Positive non-truth features (n)"),
        (nontruth_mass, "C  Non-truth native mass", "Non-truth mass (% native total)"),
    ]
    for axis, (values, title, ylabel) in zip(axes, panels):
        bars = axis.barh(y, values, color=colors, height=0.68)
        axis.set_title(title, loc="left", fontweight="bold")
        axis.set_xlabel(ylabel)
        axis.set_yticks(y, labels if axis is axes[0] else ["" for _ in labels])
        axis.invert_yaxis()
        axis.grid(axis="x", color="#DDDDDD", linewidth=0.6)
        axis.set_axisbelow(True)
        offset = max(values) * 0.025 if max(values) > 0 else 0.1
        for bar, value in zip(bars, values):
            label = f"{value:.1f}" if isinstance(value, float) else str(value)
            axis.text(bar.get_width() + offset, bar.get_y() + bar.get_height() / 2, label, ha="left", va="center", fontsize=8)
        axis.set_xlim(0, max(values) * 1.16 if max(values) > 0 else 1)
    axes[0].set_xlim(0, 108)
    for index, denom in enumerate(denominators):
        axes[0].text(2, index, f"denom={denom}", ha="left", va="center", fontsize=7, color="white" if recovery[index] > 15 else "#333333")
    fig.suptitle("Native outputs retained; detection uses a common 0.01% post-profile threshold", y=1.02, fontsize=12, fontweight="bold")
    save_figure(fig, figure_dir, "18-recovery-nontruth")


def pairwise_lookup(pairwise: list[dict[str, str]], left: str, right: str, metric: str) -> float:
    if left == right:
        return 1.0 if metric == "SpearmanRho" else 0.0
    for row in pairwise:
        if {row["ProfileA"], row["ProfileB"]} == {left, right}:
            return float(row[metric])
    raise KeyError((left, right, metric))


def render_agreement_figure(
    harmonized: list[dict[str, str]],
    pairwise: list[dict[str, str]],
    figure_dir: Path,
) -> None:
    configure_plot_style()
    matrix_ids = ["truth", *PRIMARY_PROFILES]
    matrix_labels = ["Truth", "MetaPhlAn", "Kraken/Bracken", "mOTUs"]
    rho = np.array([[pairwise_lookup(pairwise, left, right, "SpearmanRho") for right in matrix_ids] for left in matrix_ids])
    tv = np.array([[pairwise_lookup(pairwise, left, right, "TotalVariation") for right in matrix_ids] for left in matrix_ids])
    fig = plt.figure(figsize=(10.7, 8.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], hspace=0.54, wspace=0.42)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])]
    axes[0].imshow(rho, vmin=-1, vmax=1, cmap="RdBu_r")
    axes[1].imshow(tv, vmin=0, vmax=max(0.5, float(np.max(tv))), cmap="YlOrRd")
    for axis, matrix, title, fmt in (
        (axes[0], rho, "A  Spearman rho", ".2f"),
        (axes[1], tv, "B  Total variation", ".2f"),
    ):
        axis.set_title(title, loc="left", fontweight="bold")
        axis.set_xticks(range(4), matrix_labels, rotation=45, ha="right")
        axis.set_yticks(range(4), matrix_labels)
        for row in range(4):
            for column in range(4):
                value = matrix[row, column]
                axis.text(column, row, format(value, fmt), ha="center", va="center", fontsize=8, color="white" if (title.startswith("A") and abs(value) > 0.65) or (title.startswith("B") and value > 0.32) else "#222222")
    truth_rows = [row for row in harmonized if row["ProfileID"] == "metaphlan-default"]
    top = sorted(truth_rows, key=lambda row: float(row["ExpectedGenomeFraction"]), reverse=True)[:8]
    top_taxids = [row["TruthSpeciesTaxID"] for row in top]
    names = [row["TruthSpeciesName"].replace("Candidatus ", "Ca. ") for row in top]
    value_lookup = {(row["ProfileID"], row["TruthSpeciesTaxID"]): float(row["EstimatedFraction"]) for row in harmonized}
    truth_lookup = {row["TruthSpeciesTaxID"]: float(row["ExpectedGenomeFraction"]) for row in truth_rows}
    series = [
        ("Truth", [truth_lookup[taxid] for taxid in top_taxids], "#333333"),
        ("MetaPhlAn", [value_lookup[("metaphlan-default", taxid)] for taxid in top_taxids], profile_color("metaphlan-default")),
        ("Kraken/Bracken", [value_lookup[("kraken-bracken-native", taxid)] for taxid in top_taxids], profile_color("kraken-bracken-native")),
        ("mOTUs", [value_lookup[("motus-g3-default", taxid)] for taxid in top_taxids], profile_color("motus-g3-default")),
    ]
    y = np.arange(len(top_taxids))
    height = 0.18
    for index, (label, values, color) in enumerate(series):
        axes[2].barh(y + (index - 1.5) * height, values, height=height, label=label, color=color)
    axes[2].set_yticks(y, names, fontsize=7.5)
    axes[2].invert_yaxis()
    axes[2].set_xlabel("Closed fraction within strict common domain")
    axes[2].set_title("C  Highest-mass common species", loc="left", fontweight="bold")
    axes[2].grid(axis="x", color="#DDDDDD", linewidth=0.6)
    axes[2].legend(fontsize=8, loc="lower right")
    fig.suptitle("Composition agreement only after one-to-one mapping and common-domain closure", y=0.985, fontsize=12, fontweight="bold")
    save_figure(fig, figure_dir, "18-composition-agreement")


def render_resource_figure(resources: list[dict[str, str]], figure_dir: Path) -> None:
    configure_plot_style()
    labels = [
        "MetaPhlAn default",
        "Kraken/Bracken c0",
        "Kraken/Bracken c0.10",
        "mOTUs g3",
    ]
    colors = [profile_color(row["ProfileID"]) for row in resources]
    wall = [float(row["WallSeconds"]) for row in resources]
    rss = [float(row["MaxRSSMiB"]) / 1024 for row in resources]
    disk = [float(row["InstalledDatabaseBytes"]) / 1024**3 for row in resources]
    y = np.arange(len(labels))
    fig, axes = plt.subplots(1, 3, figsize=(12.7, 3.9), gridspec_kw={"wspace": 0.42})
    for axis, values, title, xlabel, digits in (
        (axes[0], wall, "A  Wall time", "Seconds (serial stage sum)", 1),
        (axes[1], rss, "B  Peak memory", "Maximum RSS (GiB)", 1),
        (axes[2], disk, "C  Installed database", "Disk footprint (GiB)", 1),
    ):
        bars = axis.barh(y, values, color=colors, height=0.66)
        axis.set_title(title, loc="left", fontweight="bold")
        axis.set_yticks(y, labels if axis is axes[0] else ["" for _ in labels])
        axis.invert_yaxis()
        axis.set_xlabel(xlabel)
        axis.grid(axis="x", color="#DDDDDD", linewidth=0.6)
        axis.set_axisbelow(True)
        offset = max(values) * 0.02
        for bar, value in zip(bars, values):
            axis.text(value + offset, bar.get_y() + bar.get_height() / 2, f"{value:.{digits}f}", va="center", fontsize=8)
        axis.set_xlim(0, max(values) * 1.2)
    fig.suptitle("Same server and input; full versus capped databases remain a design difference", y=1.02, fontsize=12, fontweight="bold")
    save_figure(fig, figure_dir, "18-resource-footprint")


def validate_figures(figure_dir: Path, checks: Checks) -> None:
    for stem in FIGURE_STEMS:
        pdf = figure_dir / f"{stem}.pdf"
        png = figure_dir / f"{stem}.png"
        tiff = figure_dir / f"{stem}.tiff"
        checks.add("figure", f"{stem}-pdf", pdf.is_file() and pdf.stat().st_size > 10_000 and pdf.read_bytes()[:4] == b"%PDF", pdf.stat().st_size if pdf.exists() else "missing")
        checks.add("figure", f"{stem}-png", png.is_file() and png.stat().st_size > 20_000, png.stat().st_size if png.exists() else "missing")
        checks.add("figure", f"{stem}-tiff", tiff.is_file() and tiff.stat().st_size > 20_000, tiff.stat().st_size if tiff.exists() else "missing")
        if png.is_file():
            with Image.open(png) as image:
                checks.add("figure", f"{stem}-png-dimensions", image.width >= 1500 and image.height >= 600, f"{image.width}x{image.height}")
        if tiff.is_file():
            with Image.open(tiff) as image:
                dpi = image.info.get("dpi", (0, 0))
                compression = image.info.get("compression", "")
                checks.add("figure", f"{stem}-tiff-dpi", min(dpi) >= 349, dpi)
                checks.add("figure", f"{stem}-tiff-lzw", str(compression).lower() in {"tiff_lzw", "lzw", "5"}, compression)


def aggregate_crosswalk(crosswalk: list[dict[str, str]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for tool, field in (
        ("MetaPhlAn", "MetaPhlAnStatus"),
        ("Kraken/Bracken", "KrakenBrackenStatus"),
        ("mOTUs", "mOTUsStatus"),
    ):
        for status in ("One-to-one", "Split", "Merge", "Split + merge", "Absent"):
            rows = [row for row in crosswalk if row[field] == status]
            output.append(
                {
                    "Tool": tool,
                    "RelationStatus": status,
                    "TruthSpecies": len(rows),
                    "ExpectedGenomePercent": math.fsum(float(row["ExpectedGenomePercent"]) for row in rows),
                }
            )
    return output


def routine_validate(args: argparse.Namespace) -> None:
    if args.output_dir is None or args.figure_dir is None:
        raise SystemExit("Routine QA requires --output-dir and --figure-dir")
    project_root = args.project_root.resolve()
    prefix = args.environment_prefix.resolve()
    frozen_dir = args.frozen_dir.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    versions = environment_checks(project_root, prefix, checks)
    source_manifest, previous = static_source_checks(project_root, checks)
    truth_ordered, truth = truth_species_table(project_root, checks)
    checksum_entries = verify_checksum_manifest(frozen_dir, checks)
    summary = json.loads((frozen_dir / "run-summary.json").read_text(encoding="utf-8"))
    checks.add("summary", "status", summary["status"] == "passed", summary["status"])
    checks.add("summary", "input-pairs", summary["input_pairs"] == EXPECTED_INPUT["pairs"], summary["input_pairs"])
    checks.add("summary", "truth-species", summary["truth_species"] == 69, summary["truth_species"])
    checks.add("summary", "motus-version", summary["motus_version"] == "4.1.0", summary["motus_version"])
    checks.add("summary", "motus-database-version", summary["motus_database_version"] == "4.1", summary["motus_database_version"])
    checks.add("summary", "motus-taxonomy-release", summary["motus_taxonomy_release"] == "R226", summary["motus_taxonomy_release"])
    checks.add("summary", "checksum-entry-count", summary["frozen_checksum_entries"] == checksum_entries, f"{summary['frozen_checksum_entries']}={checksum_entries}")
    checks.add("summary", "qa-network-free", summary["qa_network_access"] is False, summary["qa_network_access"])

    design = read_tsv(frozen_dir / "profiler-design.tsv")
    crosswalk = read_tsv(frozen_dir / "truth-feature-crosswalk.tsv")
    native = read_tsv(frozen_dir / "native-profile-summary.tsv")
    performance = read_tsv(frozen_dir / "benchmark-performance.tsv")
    harmonized = read_tsv(frozen_dir / "harmonized-profile.tsv")
    pairwise = read_tsv(frozen_dir / "pairwise-agreement.tsv")
    resources = read_tsv(frozen_dir / "resource-usage.tsv")
    database = read_tsv(frozen_dir / "database-audit.tsv")
    source_audit = read_tsv(frozen_dir / "source-audit.tsv")

    checks.add("design", "six-branches", [row["ProfileID"] for row in design] == list(PROFILE_ORDER), [row["ProfileID"] for row in design])
    checks.add("crosswalk", "sixty-nine-rows", len(crosswalk) == 69, len(crosswalk))
    common_rows = [row for row in crosswalk if row["CommonOneToOne"] == "Yes"]
    checks.add("crosswalk", "summary-common-count", len(common_rows) == summary["common_one_to_one_species"], len(common_rows))
    common_mass = math.fsum(float(row["ExpectedGenomePercent"]) for row in common_rows)
    checks.add("crosswalk", "summary-common-mass", abs(common_mass - summary["common_one_to_one_expected_percent"]) < 1e-10, common_mass)
    checks.add("crosswalk", "valid-status-vocabulary", all(row[field] in {"One-to-one", "Split", "Merge", "Split + merge", "Absent"} for row in crosswalk for field in ("MetaPhlAnStatus", "KrakenBrackenStatus", "mOTUsStatus")), "all rows")
    checks.add("crosswalk", "name-only-boundary", all(row["NameOnlyEstablishesReference"] == "No" for row in crosswalk), "all No")

    checks.add("profile", "native-summary-six", len(native) == 6, len(native))
    checks.add("profile", "performance-six", [row["ProfileID"] for row in performance] == list(PROFILE_ORDER), [row["ProfileID"] for row in performance])
    for row in performance:
        recovery = float(row["ReferenceAwareRecoveryAt0.01Pct"])
        nontruth_mass = float(row["NonTruthNativeMassPct"])
        rho = float(row["CommonDomainSpearmanVsTruth"])
        tv = float(row["CommonDomainTVVsTruth"])
        checks.add("performance", f"{row['ProfileID']}-metric-ranges", 0 <= recovery <= 1 and 0 <= nontruth_mass <= 100 and -1 <= rho <= 1 and 0 <= tv <= 1, f"recovery={recovery};mass={nontruth_mass};rho={rho};tv={tv}")
        closed = math.fsum(float(item["EstimatedFraction"]) for item in harmonized if item["ProfileID"] == row["ProfileID"])
        checks.add("performance", f"{row['ProfileID']}-closed", abs(closed - 1.0) < 1e-12, repr(closed))
    checks.add("agreement", "pairwise-21", len(pairwise) == 21, len(pairwise))
    checks.add("agreement", "pairwise-ranges", all(-1 <= float(row["SpearmanRho"]) <= 1 and 0 <= float(row["TotalVariation"]) <= 1 for row in pairwise), "all rows")
    checks.add("resource", "four-pipelines", len(resources) == 4, len(resources))
    checks.add("resource", "positive-values", all(float(row["WallSeconds"]) > 0 and float(row["MaxRSSMiB"]) > 0 and int(row["InstalledDatabaseBytes"]) > 0 for row in resources), "all rows")
    checks.add("database", "four-database-rows", len(database) == 4, len(database))
    motus_db = next(row for row in database if row["Tool"] == "mOTUs 4.1.0")
    checks.add("database", "motus-r226", motus_db["Taxonomy"] == "GTDB R226", motus_db["Taxonomy"])
    checks.add("database", "motus-count-discrepancy-recorded", "DISCREPANCY_RECORDED" in motus_db["ValidationStatus"] and motus_db["CatalogFeatures"] == "124295" and motus_db["RuntimeCatalogFeatures"] == "124300", f"{motus_db['CatalogFeatures']};{motus_db['RuntimeCatalogFeatures']}")

    render_crosswalk_figure(crosswalk, figure_dir)
    render_recovery_figure(performance, figure_dir)
    render_agreement_figure(harmonized, pairwise, figure_dir)
    render_resource_figure(resources, figure_dir)
    validate_figures(figure_dir, checks)

    crosswalk_audit = aggregate_crosswalk(crosswalk)
    write_tsv(output_dir / "source-audit.tsv", source_audit, list(source_audit[0]))
    write_tsv(output_dir / "database-audit.tsv", database, list(database[0]))
    write_tsv(output_dir / "crosswalk-audit.tsv", crosswalk_audit, list(crosswalk_audit[0]))
    write_tsv(output_dir / "native-profile-audit.tsv", native, list(native[0]))
    write_tsv(output_dir / "performance-audit.tsv", performance, list(performance[0]))
    write_tsv(output_dir / "agreement-audit.tsv", pairwise, list(pairwise[0]))
    write_tsv(output_dir / "resource-audit.tsv", resources, list(resources[0]))
    write_tsv(output_dir / "validation-audit.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])
    validation_summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "input_pairs": summary["input_pairs"],
        "truth_assemblies": summary["truth_assemblies"],
        "truth_species": summary["truth_species"],
        "common_one_to_one_species": summary["common_one_to_one_species"],
        "common_one_to_one_expected_percent": summary["common_one_to_one_expected_percent"],
        "metaphlan_version": summary["metaphlan_version"],
        "kraken2_version": summary["kraken2_version"],
        "bracken_package_version": summary["bracken_package_version"],
        "motus_version": summary["motus_version"],
        "motus_database_version": summary["motus_database_version"],
        "motus_taxonomy_release": summary["motus_taxonomy_release"],
        "motus_aligned_reads": summary["motus_aligned_reads"],
        "motus_aligned_inserts": summary["motus_aligned_inserts"],
        "profile_headlines": summary["profile_headlines"],
        "frozen_checksum_entries": checksum_entries,
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "qa_network_access": False,
    }
    write_json(output_dir / "validation-summary.json", validation_summary)
    (output_dir / "validation.log").write_text(
        f"Article 18 routine QA: {validation_summary['status']}\n"
        f"Checks: {checks.passed} passed, {checks.failed} failed\n"
        f"Frozen payloads: {checksum_entries}\n"
        f"Strict common domain: {summary['common_one_to_one_species']} species; "
        f"{summary['common_one_to_one_expected_percent']:.6f}% expected mass\n",
        encoding="utf-8",
    )
    if checks.failed:
        raise SystemExit(f"Article 18 routine QA failed {checks.failed} checks")
    print(f"PASS Article 18 routine QA: {checks.passed} checks; {checksum_entries} frozen payloads")


def main() -> None:
    args = parse_args()
    if args.initialize_frozen:
        initialize_frozen(args)
    else:
        routine_validate(args)


if __name__ == "__main__":
    main()
