#!/usr/bin/env python3
"""Offline acceptance tests for Article 73's public-resource packet."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


FIGURES = (
    "73-resource-layer-map",
    "73-gem-biome-balance",
    "73-gem-geographic-coverage",
    "73-gem-quality-audit",
    "73-tara-sampling-frame",
    "73-tara-metadata-completeness",
    "73-accession-crosswalk",
    "73-catalogue-and-gtdb-releases",
)

SOURCE_LOCKS = {
    "gem-figure1-original.png": (428457, "7e36a5130f753401362c08d86892c7ce53977642318619f0cba27249dd87ec11"),
    "gem-genome-metadata.tsv": (12350919, "fd0ad382e4ec9dbc07915333b6c2e4b53257f6d3a9f47aad7da1d2cad6e83e37"),
    "gem-readme.md": (2850, "29fcef8f24cdf6ae2d3ee0905d3a5af18cfd3abcb478042985625bba75a1ad2a"),
    "gtdb-r226-release-notes.txt": (1809, "c93851a64cc75c425e73fb9b0a472b3cd9f2546c034fb1e4c4ff1a2e7e4ae1f4"),
    "gtdb-r226-version.txt": (29, "f1a34e5c882e437196f2823c2c34281100596f273660c8682da63b529f7d62de"),
    "gtdb-r232-md5.txt": (5482, "a5c2cc52b7d319e70bb678bbe52acc6d4a697dde8f6e130b2ce39e706ba8939d"),
    "gtdb-r232-release-notes.txt": (1363, "c6fb891abcbbec1ac753d1f9d8bc920a8adb060ad0f51b41b99094daee62f2f7"),
    "gtdb-r232-version.txt": (29, "7aaa7dca8b101daaab5635cad0895051ff427223292b04913f8791ad8c53d591"),
    "mgnify-analyses-page1.json": (60273, "cc2525849f194c1eecbfec3e654fc7c378a93d84f06e0f1dc1311b67f13c85c7"),
    "mgnify-analyses-page2.json": (60174, "5e2261673ef284c24bc7d8f315bbaeaad32d01f30d0cd3bfb73b73b6d461cdad"),
    "mgnify-analyses-page3.json": (29460, "8ec9c246e980253b201619e6538c87fa974a72f3aa1f0c91a8406895bc82fee8"),
    "mgnify-catalogues.json": (29351, "e2571384392ee5686ca807b9d02d1f9f9e8187e8b2319a1f599891ffd83bef40"),
    "mgnify-marine-v2.json": (1463, "2d190ed7eab56b4c85fd17f3e37266d0966696295b482543fe8a3e5cdffee5ec"),
    "mgnify-samples-page1.json": (474961, "6d7dec16fde6bc83e9eccf3ea55d795d17b15fb9fa79439fe9e018d951b6e1e4"),
    "mgnify-samples-page2.json": (169616, "b09c27ae61c1c19273f201889cb9b19be158d5e4a8ac240f6555472d759e7030"),
    "mgnify-study.json": (3088, "610ae775d4a7958f44285df1a6d4fcbac29b47023376d4533ec1d8e3b092e385"),
}

ANCHOR_SHA256 = SOURCE_LOCKS["gem-figure1-original.png"][1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--qa-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pixel_sha(path: Path) -> str:
    with Image.open(path) as image:
        normalized = image.convert("RGBA")
        digest = hashlib.sha256()
        digest.update(f"{normalized.width}x{normalized.height}".encode())
        digest.update(normalized.tobytes())
        return digest.hexdigest()


def near(value: object, expected: float, tolerance: float = 1e-9) -> bool:
    try:
        return bool(np.isclose(float(value), expected, rtol=tolerance, atol=tolerance))
    except (TypeError, ValueError):
        return False


@dataclass
class Check:
    category: str
    check: str
    status: bool
    detail: str


class Audit:
    def __init__(self) -> None:
        self.rows: list[Check] = []

    def add(self, category: str, check: str, status: bool, detail: object = "") -> None:
        self.rows.append(Check(category, check, bool(status), str(detail)))

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "Category": row.category,
                    "Check": row.check,
                    "Status": "PASS" if row.status else "FAIL",
                    "Detail": row.detail,
                }
                for row in self.rows
            ]
        )


def audit_checksums(frozen: Path, audit: Audit) -> None:
    checksum_file = frozen / "file-checksums.sha256"
    audit.add("Bundle", "checksum-file", checksum_file.is_file(), checksum_file)
    if not checksum_file.is_file():
        return
    lines = [line for line in checksum_file.read_text().splitlines() if line]
    audit.add("Bundle", "checksum-count", len(lines) == 43, len(lines))
    for line in lines:
        parts = line.split("  ", 1)
        valid = len(parts) == 2 and bool(re.fullmatch(r"[0-9a-f]{64}", parts[0]))
        audit.add("Checksum", f"format-{len(audit.rows):03d}", valid, line)
        if not valid:
            continue
        digest, relative = parts
        path = frozen / relative
        audit.add("Checksum", relative, path.is_file() and sha256(path) == digest, digest)

    manifest_path = frozen / "bundle-manifest.json"
    audit.add("Bundle", "manifest-file", manifest_path.is_file(), manifest_path)
    if not manifest_path.is_file():
        return
    manifest = json.loads(manifest_path.read_text())
    expected = {
        "article": 73,
        "payload_files": 35,
        "script_files": 5,
        "environment_files": 2,
    }
    for field, value in expected.items():
        audit.add("Bundle", field, manifest.get(field) == value, manifest.get(field))
    contract = manifest.get("contract", "")
    audit.add("Bundle", "contract-lineage", "sample/run/analysis lineage" in contract, contract)
    audit.add("Bundle", "contract-release", "GTDB R232" in contract, contract)


def audit_sources(frozen: Path, audit: Audit) -> None:
    manifest_path = frozen / "source-manifest.json"
    audit.add("Source", "manifest-file", manifest_path.is_file(), manifest_path)
    if not manifest_path.is_file():
        return
    manifest = json.loads(manifest_path.read_text())
    expected = {
        "article": 73,
        "snapshot_date": "2026-08-23",
        "mgnify_api": "v2",
        "mgnify_study": "MGYS00000410",
        "ena_project": "PRJEB1787",
        "gem_paper_doi": "10.1038/s41587-020-0718-6",
        "gtdb_release": "R11-RS232",
        "gtdb_release_date": "2026-04-15",
        "resource_count": 16,
    }
    for field, value in expected.items():
        audit.add("Source", field, manifest.get(field) == value, manifest.get(field))
    resources = manifest.get("resources", {})
    audit.add("Source", "resource-set", set(resources) == set(SOURCE_LOCKS), sorted(resources))
    for name, (size, digest) in SOURCE_LOCKS.items():
        record = resources.get(name, {})
        audit.add("Source record", f"{name}-bytes", record.get("bytes") == size, record.get("bytes"))
        audit.add("Source record", f"{name}-sha256", record.get("sha256") == digest, record.get("sha256"))
        audit.add("Source record", f"{name}-https", str(record.get("url", "")).startswith("https://"), record.get("url"))
        path = frozen / "source" / name
        audit.add("Source file", name, path.is_file() and path.stat().st_size == size and sha256(path) == digest, path)

    anchor = resources.get("gem-figure1-original.png", {})
    audit.add("Anchor", "width", anchor.get("width") == 2116, anchor.get("width"))
    audit.add("Anchor", "height", anchor.get("height") == 918, anchor.get("height"))
    flat = frozen / "gem-figure1-original.png"
    source = frozen / "source" / "gem-figure1-original.png"
    audit.add("Anchor", "flat-source-identity", flat.is_file() and source.is_file() and sha256(flat) == sha256(source) == ANCHOR_SHA256, flat)


def audit_mgnify(frozen: Path, audit: Audit) -> None:
    registry = pd.read_csv(frozen / "resource-registry.tsv", sep="\t", dtype=str).fillna("")
    audit.add("Registry", "row-count", len(registry) == 5, len(registry))
    audit.add("Registry", "layer-unique", registry["Layer"].is_unique, registry["Layer"].tolist())
    audit.add("Registry", "resource-set", set(registry["Resource"]) == {"ENA", "MGnify API v2", "GEM", "MGnify Genomes Marine v2.0", "GTDB"}, registry["Resource"].tolist())
    for row in registry.itertuples(index=False):
        audit.add("Registry row", f"{row.Resource}-identity", bool(row.StableIdentity), row.StableIdentity)
        audit.add("Registry row", f"{row.Resource}-unit", len(row.Unit) > 15, row.Unit)
        audit.add("Registry row", f"{row.Resource}-boundary", len(row.Boundary) > 25, row.Boundary)
        audit.add("Registry row", f"{row.Resource}-url", row.URL.startswith("https://"), row.URL)

    crosswalk = pd.read_csv(frozen / "accession-crosswalk.tsv", sep="\t")
    expected_counts = {
        "MGnify study": 1,
        "ENA project IDs": 2,
        "Related BioSamples": 136,
        "Sequencing runs": 249,
        "MGnify analyses": 249,
        "DNA-labelled samples": 89,
        "RNA-labelled samples": 47,
        "Study summary downloads": 6,
    }
    audit.add("Crosswalk", "row-count", len(crosswalk) == 8, len(crosswalk))
    audit.add("Crosswalk", "object-unique", crosswalk["Object"].is_unique, crosswalk["Object"].tolist())
    indexed = crosswalk.set_index("Object")
    for name, count in expected_counts.items():
        audit.add("Crosswalk", name, name in indexed.index and int(indexed.loc[name, "Count"]) == count, indexed.loc[name, "Count"] if name in indexed.index else "MISSING")

    samples = pd.read_csv(frozen / "tara-samples.tsv", sep="\t")
    audit.add("Tara samples", "row-count", len(samples) == 136, len(samples))
    audit.add("Tara samples", "accession-unique", samples["SampleAccession"].is_unique, samples["SampleAccession"].nunique())
    audit.add("Tara samples", "sample-prefix", samples["SampleAccession"].str.startswith("SAMEA").all(), "SAMEA")
    audit.add("Tara samples", "dna-count", int(samples["NucleicAcid"].eq("DNA").sum()) == 89, samples["NucleicAcid"].value_counts().to_dict())
    audit.add("Tara samples", "rna-count", int(samples["NucleicAcid"].eq("RNA").sum()) == 47, samples["NucleicAcid"].value_counts().to_dict())
    audit.add("Tara samples", "protocol-dna", samples.loc[samples["NucleicAcid"].eq("DNA"), "ProtocolLabel"].str.contains("NUC-DNA", regex=False).all(), "NUC-DNA")
    audit.add("Tara samples", "protocol-rna", samples.loc[samples["NucleicAcid"].eq("RNA"), "ProtocolLabel"].str.contains("NUC-RNA", regex=False).all(), "NUC-RNA")
    audit.add("Tara samples", "valid-coordinates", int(samples["ValidCoordinates"].sum()) == 136, samples["ValidCoordinates"].sum())
    audit.add("Tara samples", "depth-annotation-complete", samples["DepthAnnotation"].notna().all(), samples["DepthAnnotation"].isna().sum())
    audit.add("Tara samples", "numeric-depth-count", int(samples["DepthM"].notna().sum()) == 135, samples["DepthM"].notna().sum())
    range_depth = samples.loc[samples["DepthM"].isna()]
    audit.add("Tara samples", "range-depth-row", len(range_depth) == 1 and range_depth.iloc[0]["SampleAccession"] == "SAMEA2623919" and str(range_depth.iloc[0]["DepthAnnotation"]) == "5-160", range_depth[["SampleAccession", "DepthAnnotation"]].to_dict("records"))
    audit.add("Tara samples", "site-complete", samples["SamplingSite"].notna().all(), samples["SamplingSite"].isna().sum())
    audit.add("Tara samples", "date-complete", samples["CollectionDateStart"].notna().all(), samples["CollectionDateStart"].isna().sum())
    audit.add("Tara samples", "temperature-missing", samples["Temperature"].isna().all(), samples["Temperature"].notna().sum())
    audit.add("Tara samples", "salinity-missing", samples["Salinity"].isna().all(), samples["Salinity"].notna().sum())
    for row in samples.itertuples(index=False):
        audit.add("Tara sample row", f"{row.SampleAccession}-coordinate", -90 <= float(row.Latitude) <= 90 and -180 <= float(row.Longitude) <= 180, f"{row.Latitude},{row.Longitude}")

    analyses = pd.read_csv(frozen / "tara-analyses.tsv", sep="\t")
    audit.add("Tara analyses", "row-count", len(analyses) == 249, len(analyses))
    audit.add("Tara analyses", "analysis-unique", analyses["AnalysisAccession"].is_unique, analyses["AnalysisAccession"].nunique())
    audit.add("Tara analyses", "run-unique", analyses["RunAccession"].is_unique, analyses["RunAccession"].nunique())
    audit.add("Tara analyses", "sample-count", analyses["SampleAccession"].nunique() == 136, analyses["SampleAccession"].nunique())
    audit.add("Tara analyses", "analysis-prefix", analyses["AnalysisAccession"].str.startswith("MGYA").all(), "MGYA")
    audit.add("Tara analyses", "run-prefix", analyses["RunAccession"].str.startswith("ERR").all(), "ERR")
    audit.add("Tara analyses", "study", analyses["StudyAccession"].eq("MGYS00000410").all(), analyses["StudyAccession"].unique())
    audit.add("Tara analyses", "dna-count", int(analyses["NucleicAcid"].eq("DNA").sum()) == 149, analyses["NucleicAcid"].value_counts().to_dict())
    audit.add("Tara analyses", "rna-count", int(analyses["NucleicAcid"].eq("RNA").sum()) == 100, analyses["NucleicAcid"].value_counts().to_dict())
    audit.add("Tara analyses", "v2-count", int(analyses["PipelineVersion"].eq("V2").sum()) == 248, analyses["PipelineVersion"].value_counts().to_dict())
    audit.add("Tara analyses", "v3-count", int(analyses["PipelineVersion"].eq("V3").sum()) == 1, analyses["PipelineVersion"].value_counts().to_dict())
    audit.add("Tara analyses", "samples-covered", set(analyses["SampleAccession"]) == set(samples["SampleAccession"]), "136/136")

    lineage = pd.read_csv(frozen / "tara-analysis-lineage.tsv", sep="\t")
    audit.add("Lineage", "rows", len(lineage) == 2, len(lineage))
    audit.add("Lineage", "versions", lineage["PipelineVersion"].tolist() == ["V2", "V3"], lineage["PipelineVersion"].tolist())
    audit.add("Lineage", "analyses", lineage["Analyses"].astype(int).tolist() == [248, 1], lineage["Analyses"].tolist())
    audit.add("Lineage", "runs", lineage["Runs"].astype(int).tolist() == [248, 1], lineage["Runs"].tolist())
    audit.add("Lineage", "dna", lineage["DNAAnalyses"].astype(int).tolist() == [148, 1], lineage["DNAAnalyses"].tolist())
    audit.add("Lineage", "rna", lineage["RNAAnalyses"].astype(int).tolist() == [100, 0], lineage["RNAAnalyses"].tolist())

    completeness = pd.read_csv(frozen / "tara-metadata-completeness.tsv", sep="\t")
    audit.add("Metadata", "row-count", len(completeness) == 9, len(completeness))
    expected_fields = {
        "Coordinates": (136, 100.0),
        "Numeric depth": (135, 100 * 135 / 136),
        "Sampling site": (136, 100.0),
        "Protocol label": (136, 100.0),
        "Collection date": (136, 100.0),
        "Environment feature": (136, 100.0),
        "Environment material": (136, 100.0),
        "Temperature": (0, 0.0),
        "Salinity": (0, 0.0),
    }
    audit.add("Metadata", "field-set", set(completeness["Field"]) == set(expected_fields), completeness["Field"].tolist())
    for row in completeness.itertuples(index=False):
        available, percent = expected_fields[row.Field]
        audit.add("Metadata field", row.Field, int(row.Available) == available and int(row.Total) == 136 and near(row.CompletenessPct, percent), f"{row.Available}/{row.Total}={row.CompletenessPct}")

    downloads = pd.read_csv(frozen / "tara-study-downloads.tsv", sep="\t", dtype=str).fillna("")
    audit.add("Downloads", "row-count", len(downloads) == 6, len(downloads))
    for index, row in enumerate(downloads.itertuples(index=False), start=1):
        audit.add("Download row", f"{index}-name", bool(getattr(row, "DownloadName", "") or getattr(row, "Description", "") or str(row)), str(row))


def audit_catalogues_and_gtdb(frozen: Path, audit: Audit) -> None:
    catalogues = pd.read_csv(frozen / "mgnify-catalogues.tsv", sep="\t")
    audit.add("Catalogues", "row-count", len(catalogues) == 19, len(catalogues))
    audit.add("Catalogues", "id-unique", catalogues["CatalogueID"].is_unique, catalogues["CatalogueID"].nunique())
    audit.add("Catalogues", "published", catalogues["Status"].eq("published").all(), catalogues["Status"].value_counts().to_dict())
    audit.add("Catalogues", "positive-input", catalogues["InputGenomes"].gt(0).all(), catalogues["InputGenomes"].min())
    audit.add("Catalogues", "representatives-bounded", catalogues["RepresentativeClusters"].le(catalogues["InputGenomes"]).all(), "representatives <= inputs")
    for row in catalogues.itertuples(index=False):
        audit.add("Catalogue row", f"{row.CatalogueID}-version", bool(str(row.Version)), row.Version)
        audit.add("Catalogue row", f"{row.CatalogueID}-pipeline", bool(str(row.PipelineVersion)), row.PipelineVersion)
        audit.add("Catalogue row", f"{row.CatalogueID}-url", str(row.FTPURL).startswith(("http://", "https://")), row.FTPURL)
        audit.add("Catalogue row", f"{row.CatalogueID}-boundary", int(row.RepresentativeClusters) <= int(row.InputGenomes), f"{row.RepresentativeClusters}/{row.InputGenomes}")
    indexed = catalogues.set_index("CatalogueID")
    exact = {
        "soil-v1-0": (19472, 20908),
        "marine-v2-0": (13223, 50866),
        "marine-sediment-v1-0": (6158, 10126),
    }
    for catalogue_id, (representatives, inputs) in exact.items():
        row = indexed.loc[catalogue_id]
        audit.add("Catalogue identity", f"{catalogue_id}-counts", (int(row.RepresentativeClusters), int(row.InputGenomes)) == (representatives, inputs), (row.RepresentativeClusters, row.InputGenomes))
    marine = indexed.loc["marine-v2-0"]
    audit.add("Catalogue identity", "marine-pipeline", marine.PipelineVersion == "v2.3.0", marine.PipelineVersion)

    history = pd.read_csv(frozen / "gtdb-release-history.tsv", sep="\t")
    expected = [("R10-RS226", "2025-04-16", 732475, 143614), ("R11-RS232", "2026-04-15", 901341, 199923)]
    audit.add("GTDB history", "rows", len(history) == 2, len(history))
    for index, values in enumerate(expected):
        row = history.iloc[index]
        observed = (row.Release, row.ReleaseDate, int(row.Genomes), int(row.SpeciesClusters))
        audit.add("GTDB release", values[0], observed == values, observed)
    audit.add("GTDB history", "genome-increase", int(history.iloc[1].Genomes) > int(history.iloc[0].Genomes), history["Genomes"].tolist())
    audit.add("GTDB history", "cluster-increase", int(history.iloc[1].SpeciesClusters) > int(history.iloc[0].SpeciesClusters), history["SpeciesClusters"].tolist())

    files = pd.read_csv(frozen / "gtdb-selected-files.tsv", sep="\t")
    audit.add("GTDB files", "rows", len(files) == 5, len(files))
    audit.add("GTDB files", "purpose-unique", files["Purpose"].is_unique, files["Purpose"].tolist())
    for row in files.itertuples(index=False):
        audit.add("GTDB file", f"{row.Purpose}-bytes", int(row.Bytes) > 0, row.Bytes)
        audit.add("GTDB file", f"{row.Purpose}-md5", bool(re.fullmatch(r"[0-9a-f]{32}", row.MD5)), row.MD5)
        audit.add("GTDB file", f"{row.Purpose}-url", row.URL.startswith("https://data.gtdb.ecogenomic.org/releases/release232/232.0/"), row.URL)
    package = files.set_index("Purpose").loc["GTDB-Tk full reference package"]
    audit.add("GTDB files", "package-bytes", int(package.Bytes) == 60806405195, package.Bytes)
    audit.add("GTDB files", "package-md5", package.MD5 == "25a59e0352b1fd150c589f56559767d4", package.MD5)


def audit_gem(frozen: Path, audit: Audit) -> None:
    raw = pd.read_csv(frozen / "source" / "gem-genome-metadata.tsv", sep="\t", low_memory=False)
    audit.add("GEM", "row-count", len(raw) == 52515, len(raw))
    audit.add("GEM", "genome-unique", raw["genome_id"].is_unique, raw["genome_id"].nunique())
    audit.add("GEM", "metagenomes-with-mags", raw["metagenome_id"].nunique() == 7304, raw["metagenome_id"].nunique())
    audit.add("GEM", "species-otus", raw["otu_id"].nunique() == 18028, raw["otu_id"].nunique())
    audit.add("GEM", "hq", int(raw["mimag_quality"].eq("HQ").sum()) == 9143, raw["mimag_quality"].value_counts().to_dict())
    audit.add("GEM", "mq", int(raw["mimag_quality"].eq("MQ").sum()) == 43372, raw["mimag_quality"].value_counts().to_dict())
    audit.add("GEM", "quality-labels", set(raw["mimag_quality"]) == {"HQ", "MQ"}, raw["mimag_quality"].value_counts().to_dict())
    audit.add("GEM", "mean-completeness", near(raw["completeness"].mean(), 83.03630258021516), raw["completeness"].mean())
    audit.add("GEM", "mean-contamination", near(raw["contamination"].mean(), 1.2672501190136154), raw["contamination"].mean())
    audit.add("GEM", "quality-score-identity", np.allclose(raw["quality_score"], raw["completeness"] - 5 * raw["contamination"], atol=1e-8), "completeness - 5*contamination")
    audit.add("GEM", "minimum-completeness", raw["completeness"].ge(50).all(), raw["completeness"].min())
    audit.add("GEM", "maximum-contamination", raw["contamination"].le(5).all(), raw["contamination"].max())
    georef = raw["longitude"].notna() & raw["latitude"].notna() & raw["longitude"].between(-180, 180) & raw["latitude"].between(-90, 90)
    audit.add("GEM", "georeferenced", int(georef.sum()) == 49521, georef.sum())

    summary = pd.read_csv(frozen / "gem-biome-summary.tsv", sep="\t")
    audit.add("GEM summary", "rows", len(summary) == 30, len(summary))
    audit.add("GEM summary", "category-unique", summary["EcosystemCategory"].is_unique, summary["EcosystemCategory"].tolist())
    audit.add("GEM summary", "mags-total", int(summary["MAGs"].sum()) == 52515, summary["MAGs"].sum())
    audit.add("GEM summary", "hq-total", int(summary["HighQuality"].sum()) == 9143, summary["HighQuality"].sum())
    for row in summary.itertuples(index=False):
        audit.add("GEM biome", f"{row.EcosystemCategory}-mags", int(row.MAGs) > 0, row.MAGs)
        audit.add("GEM biome", f"{row.EcosystemCategory}-hq-bounded", 0 <= int(row.HighQuality) <= int(row.MAGs), f"{row.HighQuality}/{row.MAGs}")
        audit.add("GEM biome", f"{row.EcosystemCategory}-hq-pct", near(row.HighQualityPct, 100 * row.HighQuality / row.MAGs), row.HighQualityPct)
    top = summary.set_index("EcosystemCategory")
    audit.add("GEM summary", "aquatic", int(top.loc["Aquatic", "MAGs"]) == 19300, top.loc["Aquatic", "MAGs"])
    audit.add("GEM summary", "human", int(top.loc["Human", "MAGs"]) == 16441, top.loc["Human", "MAGs"])

    preview = pd.read_csv(frozen / "gem-metadata-preview.tsv", sep="\t")
    audit.add("GEM preview", "rows", len(preview) == 12, len(preview))
    audit.add("GEM preview", "columns", len(preview.columns) == 20, preview.columns.tolist())
    audit.add("GEM preview", "subset", set(preview["genome_id"]).issubset(set(raw["genome_id"])), "12/12")

    quality = pd.read_csv(frozen / "gem-quality-visualization-sample.tsv", sep="\t")
    audit.add("GEM plot sample", "quality-rows", len(quality) == 5600, len(quality))
    audit.add("GEM plot sample", "quality-unique", quality["GenomeID"].is_unique, quality["GenomeID"].nunique())
    audit.add("GEM plot sample", "quality-categories", quality["EcosystemCategory"].nunique() == 8, quality["EcosystemCategory"].value_counts().to_dict())
    audit.add("GEM plot sample", "quality-700-each", quality["EcosystemCategory"].value_counts().eq(700).all(), quality["EcosystemCategory"].value_counts().to_dict())
    audit.add("GEM plot sample", "quality-identity", np.allclose(quality["QualityScore"], quality["Completeness"] - 5 * quality["Contamination"], atol=1e-8), "5,600 rows")

    geo = pd.read_csv(frozen / "gem-map-visualization-sample.tsv", sep="\t")
    audit.add("GEM plot sample", "map-rows", len(geo) == 12000, len(geo))
    audit.add("GEM plot sample", "map-unique", geo["GenomeID"].is_unique, geo["GenomeID"].nunique())
    audit.add("GEM plot sample", "map-longitude", geo["Longitude"].between(-180, 180).all(), (geo["Longitude"].min(), geo["Longitude"].max()))
    audit.add("GEM plot sample", "map-latitude", geo["Latitude"].between(-90, 90).all(), (geo["Latitude"].min(), geo["Latitude"].max()))


def audit_contracts_and_environment(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "analysis-metrics.json").read_text())
    expected = {
        "article": 73,
        "snapshot_date": "2026-08-23",
        "analysis_seed": 73001,
        "plot_seed": 20260773,
        "mgnify_api": "v2",
        "mgnify_study": "MGYS00000410",
        "mgnify_samples": 136,
        "mgnify_dna_samples": 89,
        "mgnify_rna_samples": 47,
        "mgnify_runs": 249,
        "mgnify_analyses": 249,
        "mgnify_v2_analyses": 248,
        "mgnify_v3_analyses": 1,
        "mgnify_catalogues": 19,
        "marine_input_genomes": 50866,
        "marine_representative_clusters": 13223,
        "gem_mags": 52515,
        "gem_source_metagenomes_with_mags": 7304,
        "gem_species_otus": 18028,
        "gem_high_quality": 9143,
        "gem_medium_quality": 43372,
        "gem_georeferenced": 49521,
        "gtdb_release": "R11-RS232",
        "gtdb_genomes": 901341,
        "gtdb_species_clusters": 199923,
        "gtdbtk_reference_bytes": 60806405195,
    }
    for field, value in expected.items():
        audit.add("Metrics", field, metrics.get(field) == value, metrics.get(field))
    for field in ("python", "pandas", "numpy"):
        audit.add("Metrics", f"{field}-version", bool(re.fullmatch(r"\d+\.\d+\.\d+", str(metrics.get(field, "")))), metrics.get(field))
    audit.add("Metrics", "mean-completeness", near(metrics.get("gem_mean_completeness"), 83.03630258021516), metrics.get("gem_mean_completeness"))
    audit.add("Metrics", "mean-contamination", near(metrics.get("gem_mean_contamination"), 1.2672501190136154), metrics.get("gem_mean_contamination"))

    contract = json.loads((frozen / "methods-contract.json").read_text())
    for field, value in {"article": 73, "analysis_seed": 73001, "plot_seed": 20260773, "snapshot_date": "2026-08-23"}.items():
        audit.add("Methods contract", field, contract.get(field) == value, contract.get(field))
    units = contract.get("unit_contract", {})
    audit.add("Methods contract", "unit-keys", set(units) == {"archive", "analysis", "gem", "mgnify_catalogue", "gtdb"}, units)
    for key, value in units.items():
        audit.add("Methods unit", key, len(value) > 15, value)
    audit.add("Methods contract", "dna-eligibility", "NUC-DNA" in contract.get("eligibility_rule", "") and "study title" in contract.get("eligibility_rule", ""), contract.get("eligibility_rule"))
    audit.add("Methods contract", "catalogue-limit", "not abundance" in contract.get("catalogue_limit", "") and "not globally dereplicated" in contract.get("catalogue_limit", ""), contract.get("catalogue_limit"))
    audit.add("Methods contract", "taxonomy-limit", "release-bound" in contract.get("taxonomy_limit", "") and "side by side" in contract.get("taxonomy_limit", ""), contract.get("taxonomy_limit"))

    decisions = pd.read_csv(frozen / "download-decisions.tsv", sep="\t", dtype=str).fillna("")
    audit.add("Download plan", "rows", len(decisions) == 7, len(decisions))
    for index, row in enumerate(decisions.itertuples(index=False), start=1):
        audit.add("Download plan", f"{index}-question", len(row.Question) > 8, row.Question)
        audit.add("Download plan", f"{index}-minimum", len(row.MinimumDownload) > 8, row.MinimumDownload)
        audit.add("Download plan", f"{index}-audit", len(row.AuditRequirement) > 12, row.AuditRequirement)

    python_env = (frozen / "env" / "multiomics-python.yml").read_text()
    for pin in ("python=3.10.19", "numpy=2.2.6", "pandas=2.3.3", "pillow=12.2.0", "matplotlib==3.10.8"):
        audit.add("Environment", pin, pin in python_env, pin)
    r_env = pd.read_csv(frozen / "env" / "multiomics-r-packages.tsv", sep="\t", dtype=str)
    for package in ("R", "ggplot2", "jsonlite"):
        audit.add("Environment", f"r-{package}", package in set(r_env["Package"]), package)


def audit_chapter(root: Path, audit: Audit) -> None:
    chapter = root / "chapters" / "73-public-metagenome-resources.qmd"
    audit.add("Chapter", "exists", chapter.is_file(), chapter)
    if not chapter.is_file():
        return
    text = chapter.read_text(encoding="utf-8")
    audit.add("Chapter", "published", "draft: false" in text, "draft: false")
    audit.add("Chapter", "eval-true", "eval: true" in text, "eval: true")
    audit.add("Chapter", "freeze-auto", "freeze: auto" in text, "freeze: auto")
    audit.add("Chapter", "image-count", text.count("](../figures/73-") == 9, text.count("](../figures/73-"))
    audit.add("Chapter", "native-fences", text.count("```{r}") >= 9 and text.count("```{bash}") >= 3 and "~~~{" not in text, "9 R + 3 bash minimum")
    for heading in (
        "对应论文里的哪张图",
        "理论：",
        "准备工作",
        "可复制代码",
        "审计与升级",
        "出版级美化",
        "常见坑",
        "这段 Methods 怎么写",
        "换成你自己的数据怎么做",
        "参考",
    ):
        audit.add("Chapter structure", heading, heading in text, heading)
    phrases = (
        "89 个 `NUC-DNA` 与 47 个 `NUC-RNA`",
        "149 仍是 DNA analyses",
        "不能相加成“32,695 个全球物种”",
        "不表示真实环境数值为零",
        "not recovered in catalogue version X under criterion Y",
        "Taxonomic names are release-bound",
        "Catalogue representatives were not interpreted as sample abundance",
        "60.81 GB",
        "7,304",
        "18,028",
    )
    for index, phrase in enumerate(phrases, start=1):
        audit.add("Chapter boundary", f"boundary-{index}", phrase in text, phrase)
    audit.add("Chapter", "seeds", "73001" in text and "set.seed(20260773)" in text, "73001 / 20260773")
    audit.add("Chapter", "runtime-versions", all(value in text for value in ("Python 3.13.2", "pandas 3.0.3", "NumPy 2.5.0")), "runtime versions")
    audit.add("Chapter", "inline-theme", all(value in text for value in ("theme_pub <- function", "save_pub <- function", "scale_color_pub <- function", "scale_fill_pub  <- function")), "inline plotting functions")
    audit.add("Chapter", "no-source-theme", 'source("R/theme_pub.R")' not in text, "no source dependency")
    audit.add("Chapter", "frozen-input", "data/small/73-public-resources-frozen" in text, "frozen bundle")
    audit.add("Chapter", "checksum-command", "sha256sum -c" in text, "checksum verification")
    forbidden = (
        "本篇可独立跑通",
        "这体现全系列",
        "接口只学一次",
        "作者代码通常长这样",
        "（即本文）",
        "/media/desk16",
        "/tmp/article73",
        "qa_report",
        "review draft",
    )
    for phrase in forbidden:
        audit.add("Chapter prose", f"forbidden-{phrase}", phrase not in text, phrase)
    for stem in FIGURES:
        audit.add("Chapter figure", stem, f"../figures/{stem}.png" in text, stem)
    audit.add("Chapter figure", "anchor", "../figures/73-gem-figure1-original.png" in text, "official anchor")
    for key in ("sunagawa2015ocean", "richardson2023mgnify", "nayfach2021gem", "gurbich2023mgnifygenomes", "parks2022gtdb", "parks2025gtdbr10", "bowers2017mimag"):
        audit.add("Chapter citation", key, f"@{key}" in text, key)


def audit_figures(root: Path, frozen: Path, audit: Audit) -> None:
    figures = root / "figures"
    for stem in FIGURES:
        for suffix in ("pdf", "png", "tiff"):
            path = figures / f"{stem}.{suffix}"
            audit.add("Figure file", f"{stem}.{suffix}", path.is_file() and path.stat().st_size > 10_000, path.stat().st_size if path.is_file() else "MISSING")
        png = figures / f"{stem}.png"
        tiff = figures / f"{stem}.tiff"
        pdf = figures / f"{stem}.pdf"
        if png.is_file():
            with Image.open(png) as image:
                dpi = image.info.get("dpi", (0, 0))
                audit.add("Figure raster", f"{stem}-png-size", image.width >= 1800 and image.height >= 1100, image.size)
                audit.add("Figure raster", f"{stem}-png-dpi", min(dpi) >= 300, dpi)
        if tiff.is_file():
            with Image.open(tiff) as image:
                dpi = image.info.get("dpi", (0, 0))
                audit.add("Figure raster", f"{stem}-tiff-dpi", min(dpi) >= 300, dpi)
                audit.add("Figure raster", f"{stem}-tiff-lzw", image.tag_v2.get(259) == 5, image.tag_v2.get(259))
        if pdf.is_file():
            try:
                result = subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, text=True, timeout=30)
                audit.add("Figure text", f"{stem}-pdf-text", result.returncode == 0 and len(result.stdout.strip()) > 10, len(result.stdout))
                audit.add("Figure text", f"{stem}-english-only", not bool(re.search(r"[\u3400-\u9fff]", result.stdout)), "no CJK glyphs")
            except FileNotFoundError:
                audit.add("Figure text", f"{stem}-pdftotext", False, "pdftotext missing")

    anchor = figures / "73-gem-figure1-original.png"
    frozen_anchor = frozen / "gem-figure1-original.png"
    audit.add("Figure file", "anchor", anchor.is_file() and frozen_anchor.is_file() and sha256(anchor) == sha256(frozen_anchor) == ANCHOR_SHA256, sha256(anchor) if anchor.is_file() else "MISSING")
    if anchor.is_file():
        with Image.open(anchor) as image:
            audit.add("Figure raster", "anchor-dimensions", image.size == (2116, 918), image.size)

    with tempfile.TemporaryDirectory(prefix="article73-replot-") as temporary:
        staged = Path(temporary) / "figures"
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(Path(temporary) / "mpl")
        result = subprocess.run(
            [
                sys.executable,
                str(frozen / "scripts" / "plot_article73_resources.py"),
                "--input-dir",
                str(frozen),
                "--figure-dir",
                str(staged),
            ],
            capture_output=True,
            text=True,
            timeout=240,
            env=environment,
        )
        audit.add("Reanalysis", "plot-script-exit", result.returncode == 0, result.stdout + result.stderr)
        for stem in FIGURES:
            staged_png = staged / f"{stem}.png"
            published_png = figures / f"{stem}.png"
            status = staged_png.is_file() and published_png.is_file() and pixel_sha(staged_png) == pixel_sha(published_png)
            audit.add("Reanalysis", f"{stem}-pixel-identical", status, pixel_sha(staged_png) if staged_png.is_file() else "MISSING")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    qa = args.qa_dir.resolve()
    qa.mkdir(parents=True, exist_ok=True)
    audit = Audit()
    audit_checksums(frozen, audit)
    audit_sources(frozen, audit)
    audit_mgnify(frozen, audit)
    audit_catalogues_and_gtdb(frozen, audit)
    audit_gem(frozen, audit)
    audit_contracts_and_environment(frozen, audit)
    audit_chapter(root, audit)
    audit_figures(root, frozen, audit)
    report = audit.frame()
    report.to_csv(qa / "qa_report.tsv", sep="\t", index=False)
    passed = int(report["Status"].eq("PASS").sum())
    failed = int(report["Status"].eq("FAIL").sum())
    payload = {
        "article": 73,
        "status": "passed" if failed == 0 else "failed",
        "checks": len(report),
        "passed": passed,
        "failed": failed,
        "failed_checks": report.loc[report["Status"].eq("FAIL"), ["Category", "Check", "Detail"]].to_dict("records"),
    }
    (qa / "qa_report.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
