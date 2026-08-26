#!/usr/bin/env python3
"""Offline acceptance tests for Article 72's causal-evidence packet."""

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
    "72-evidence-ladder-framework",
    "72-claim-contracts",
    "72-evidence-coverage",
    "72-human-intervention",
    "72-rcdi-evidence-braid",
    "72-fusobacterium-gap",
    "72-claim-downgrade",
    "72-evidence-packet",
)

PUBLICATIONS = {
    "vannood2013fmt": ("10.1056/nejmoa1205037", 2013),
    "feuerstadt2022ser109": ("10.1056/nejmoa2106516", 2022),
    "smillie2018engraftment": ("10.1016/j.chom.2018.01.003", 2018),
    "li2016fmtstrains": ("10.1126/science.aad8852", 2016),
    "buffie2015cscindens": ("10.1038/nature13828", 2015),
    "theriot2014metabolome": ("10.1038/ncomms4114", 2014),
    "wirbel2019crc": ("10.1038/s41591-019-0406-6", 2019),
    "kostic2013fusobacterium": ("10.1016/j.chom.2013.07.007", 2013),
    "rubinstein2013fada": ("10.1016/j.chom.2013.07.012", 2013),
    "bullman2017fusobacterium": ("10.1126/science.aal5240", 2017),
    "neville2018commensalkoch": ("10.1016/j.mib.2017.10.001", 2018),
}

RUNGS = (
    "Replicated association",
    "Temporality and reversibility",
    "Human causal bridge",
    "Randomized human intervention",
    "Isolation or defined product",
    "Host transfer or perturbation",
    "Molecular perturbation and rescue",
)

ANCHOR_BYTES = 143_655
ANCHOR_SHA256 = (
    "aaa42762d0f3c9287a266a2681774a90b908b64d9614dbf8ec1d1b93588792c9"
)
METADATA_BYTES = 23_575
METADATA_SHA256 = (
    "c48e3ac9ef4ecf0da7c1f92e8dbb59f93b42f3841a79960284de5a9cc2825a53"
)


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


def near(value: object, expected: float, tolerance: float = 1e-10) -> bool:
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

    def add(
        self, category: str, check: str, status: bool, detail: object = ""
    ) -> None:
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
    audit.add("Bundle", "checksum-count", len(lines) == 20, len(lines))
    for line in lines:
        parts = line.split("  ", 1)
        well_formed = len(parts) == 2 and bool(re.fullmatch(r"[0-9a-f]{64}", parts[0]))
        audit.add("Checksum", f"format-{len(audit.rows):03d}", well_formed, line)
        if not well_formed:
            continue
        digest, relative = parts
        path = frozen / relative
        audit.add(
            "Checksum",
            relative,
            path.is_file() and sha256(path) == digest,
            digest,
        )

    manifest_path = frozen / "bundle-manifest.json"
    audit.add("Bundle", "manifest-file", manifest_path.is_file(), manifest_path)
    if not manifest_path.is_file():
        return
    manifest = json.loads(manifest_path.read_text())
    expected = {
        "article": 72,
        "payload_files": 12,
        "script_files": 5,
        "environment_files": 2,
    }
    for field, value in expected.items():
        audit.add("Bundle", field, manifest.get(field) == value, manifest.get(field))
    audit.add(
        "Bundle",
        "contract-non-additive",
        "non-additive" in manifest.get("contract", "")
        and "official anchor" in manifest.get("contract", ""),
        manifest.get("contract"),
    )


def audit_sources(frozen: Path, audit: Audit) -> None:
    path = frozen / "source-manifest.json"
    audit.add("Source", "manifest-file", path.is_file(), path)
    if not path.is_file():
        return
    manifest = json.loads(path.read_text())
    expected = {
        "article": 72,
        "metadata_source": "Crossref REST API",
        "publication_count": 11,
        "publication_rows": 11,
        "evidence_rows": 10,
        "claim_count": 3,
        "rung_count": 7,
        "intervention_arm_rows": 5,
        "framework_citation": "10.1016/j.mib.2017.10.001",
    }
    for field, value in expected.items():
        audit.add("Source", field, manifest.get(field) == value, manifest.get(field))

    observed_dois = set(manifest.get("doi_set", []))
    expected_dois = {value[0] for value in PUBLICATIONS.values()}
    audit.add("Source", "doi-set", observed_dois == expected_dois, sorted(observed_dois))
    for key, (doi, _) in PUBLICATIONS.items():
        audit.add("Source DOI", key, doi in observed_dois, doi)

    resource = manifest.get("resources", {}).get("publication-metadata.json", {})
    audit.add("Source", "metadata-bytes", resource.get("bytes") == METADATA_BYTES, resource)
    audit.add("Source", "metadata-sha256", resource.get("sha256") == METADATA_SHA256, resource)
    audit.add(
        "Source",
        "metadata-url",
        resource.get("url_template") == "https://api.crossref.org/works/{doi}",
        resource.get("url_template"),
    )

    anchor = manifest.get("anchor", {})
    anchor_expected = {
        "name": "buffie-figure4-original.jpg",
        "bytes": ANCHOR_BYTES,
        "sha256": ANCHOR_SHA256,
        "width": 1578,
        "height": 707,
        "paper_doi": "10.1038/nature13828",
        "figure": "Figure 4",
    }
    for field, value in anchor_expected.items():
        audit.add("Anchor", field, anchor.get(field) == value, anchor.get(field))
    audit.add(
        "Anchor",
        "official-https-url",
        str(anchor.get("url", "")).startswith("https://media.springernature.com/"),
        anchor.get("url"),
    )
    frozen_anchor = frozen / "buffie-figure4-original.jpg"
    audit.add(
        "Anchor",
        "frozen-identity",
        frozen_anchor.is_file()
        and frozen_anchor.stat().st_size == ANCHOR_BYTES
        and sha256(frozen_anchor) == ANCHOR_SHA256,
        sha256(frozen_anchor) if frozen_anchor.is_file() else "MISSING",
    )
    if frozen_anchor.is_file():
        with Image.open(frozen_anchor) as image:
            audit.add("Anchor", "frozen-dimensions", image.size == (1578, 707), image.size)


def audit_publications(frozen: Path, audit: Audit) -> pd.DataFrame:
    publications = pd.read_csv(frozen / "publication-metadata.tsv", sep="\t", dtype=str).fillna("")
    audit.add("Publications", "row-count", len(publications) == 11, len(publications))
    audit.add("Publications", "citation-unique", publications["CitationKey"].is_unique, "unique")
    audit.add("Publications", "doi-unique", publications["DOI"].is_unique, "unique")
    audit.add("Publications", "key-set", set(publications["CitationKey"]) == set(PUBLICATIONS), sorted(publications["CitationKey"]))
    indexed = publications.set_index("CitationKey")
    for key, (doi, year) in PUBLICATIONS.items():
        present = key in indexed.index
        audit.add("Publication", f"{key}-present", present, key)
        if not present:
            continue
        row = indexed.loc[key]
        audit.add("Publication", f"{key}-doi", row["DOI"].lower() == doi, row["DOI"])
        audit.add("Publication", f"{key}-year", int(row["Year"]) == year, row["Year"])
        audit.add("Publication", f"{key}-doi-url", row["DOIURL"] == f"https://doi.org/{doi}", row["DOIURL"])
        audit.add("Publication", f"{key}-title", len(row["Title"].strip()) > 12, row["Title"])
        audit.add("Publication", f"{key}-journal", len(row["Journal"].strip()) > 2, row["Journal"])
        audit.add("Publication", f"{key}-author", bool(row["FirstAuthor"].strip()), row["FirstAuthor"])
        dirty = bool(re.search(r"<[^>]+>|&(?:[A-Za-z]+|#\d+);|[\r\n]", row["Title"]))
        audit.add("Publication", f"{key}-clean-title", not dirty, row["Title"])
    return publications


def audit_contracts(frozen: Path, audit: Audit) -> tuple[pd.DataFrame, pd.DataFrame]:
    claims = pd.read_csv(frozen / "claim-contracts.tsv", sep="\t", dtype=str).fillna("")
    audit.add("Claims", "row-count", len(claims) == 3, len(claims))
    audit.add("Claims", "ids", claims["ClaimID"].tolist() == ["C1", "C2", "C3"], claims["ClaimID"].tolist())
    required = (
        "Case",
        "ExactClaim",
        "TargetPopulation",
        "InterventionOrExposure",
        "Comparator",
        "Outcome",
        "TimeHorizon",
        "AllowedConclusion",
        "ForbiddenLeap",
    )
    for row in claims.itertuples(index=False):
        for field in required:
            value = getattr(row, field)
            audit.add("Claim contract", f"{row.ClaimID}-{field}", bool(value.strip()), value)
        audit.add(
            "Claim contract",
            f"{row.ClaimID}-bounded",
            row.AllowedConclusion != row.ForbiddenLeap and len(row.ForbiddenLeap) > 40,
            row.ForbiddenLeap,
        )

    rungs = pd.read_csv(frozen / "rung-definitions.tsv", sep="\t", dtype=str).fillna("")
    audit.add("Rungs", "row-count", len(rungs) == 7, len(rungs))
    audit.add("Rungs", "order", rungs["RungOrder"].astype(int).tolist() == list(range(1, 8)), rungs["RungOrder"].tolist())
    audit.add("Rungs", "names", tuple(rungs["Rung"]) == RUNGS, rungs["Rung"].tolist())
    for row in rungs.itertuples(index=False):
        audit.add("Rung", f"{row.RungOrder}-question", len(row.Question) > 15, row.Question)
        audit.add("Rung", f"{row.RungOrder}-minimum-design", len(row.MinimumDesign) > 20, row.MinimumDesign)
        audit.add("Rung", f"{row.RungOrder}-addition", len(row.WhatItAdds) > 15, row.WhatItAdds)
        audit.add("Rung", f"{row.RungOrder}-boundary", len(row.DoesNotProve) > 15, row.DoesNotProve)

    methods = json.loads((frozen / "methods-contract.json").read_text())
    expected = {
        "article": 72,
        "analysis_seed": 72001,
        "plot_seed": 20260772,
        "unit": "one primary publication by one prespecified claim",
    }
    for field, value in expected.items():
        audit.add("Methods contract", field, methods.get(field) == value, methods.get(field))
    audit.add("Methods contract", "claims", methods.get("claims") == ["C1", "C2", "C3"], methods.get("claims"))
    audit.add("Methods contract", "rungs", tuple(methods.get("rungs", [])) == RUNGS, methods.get("rungs"))
    audit.add(
        "Methods contract",
        "coverage-codes",
        set(methods.get("coverage_codes", [])) == {"Not addressed", "Supporting", "Direct"},
        methods.get("coverage_codes"),
    )
    audit.add(
        "Methods contract",
        "non-additive",
        "without an additive causal score" in methods.get("synthesis", "")
        and "not a quality-weighted meta-analysis" in methods.get("interpretation_limit", ""),
        methods,
    )
    return claims, rungs


def audit_evidence(
    frozen: Path, publications: pd.DataFrame, rungs: pd.DataFrame, audit: Audit
) -> None:
    evidence = pd.read_csv(frozen / "evidence-ledger.tsv", sep="\t")
    audit.add("Evidence", "row-count", len(evidence) == 10, len(evidence))
    audit.add("Evidence", "ids", evidence["EvidenceID"].tolist() == [f"E{i:02d}" for i in range(1, 11)], evidence["EvidenceID"].tolist())
    audit.add("Evidence", "id-unique", evidence["EvidenceID"].is_unique, "unique")
    audit.add("Evidence", "publication-unique", evidence["CitationKey"].is_unique, "one row per primary publication")
    audit.add("Evidence", "claim-counts", evidence["ClaimID"].value_counts().to_dict() == {"C1": 4, "C3": 4, "C2": 2}, evidence["ClaimID"].value_counts().to_dict())
    audit.add("Evidence", "citation-set", set(evidence["CitationKey"]) == set(PUBLICATIONS) - {"neville2018commensalkoch"}, sorted(evidence["CitationKey"]))
    audit.add("Evidence", "roles", set(evidence["EvidenceRole"]) == {"Direct", "Supporting"}, sorted(evidence["EvidenceRole"].unique()))
    audit.add("Evidence", "shotgun-count", int(evidence["ShotgunStrainResolved"].sum()) == 3, evidence["ShotgunStrainResolved"].sum())
    audit.add("Evidence", "replication-count", int(evidence["IndependentReplication"].sum()) == 8, evidence["IndependentReplication"].sum())

    pubs = publications.set_index("CitationKey")
    for row in evidence.itertuples(index=False):
        audit.add("Evidence row", f"{row.EvidenceID}-key", row.CitationKey in pubs.index, row.CitationKey)
        if row.CitationKey in pubs.index:
            publication = pubs.loc[row.CitationKey]
            audit.add("Evidence row", f"{row.EvidenceID}-doi", row.DOI == publication.DOI, row.DOI)
            audit.add("Evidence row", f"{row.EvidenceID}-year", int(row.Year) == int(publication.Year), row.Year)
            audit.add("Evidence row", f"{row.EvidenceID}-title", row.Title == publication.Title, row.Title)
            audit.add("Evidence row", f"{row.EvidenceID}-journal", row.Journal == publication.Journal, row.Journal)
        audit.add("Evidence row", f"{row.EvidenceID}-design", len(str(row.StudyDesign)) > 20, row.StudyDesign)
        audit.add("Evidence row", f"{row.EvidenceID}-population", len(str(row.PopulationModel)) > 15, row.PopulationModel)
        audit.add("Evidence row", f"{row.EvidenceID}-observation", len(str(row.MainObservation)) > 40, row.MainObservation)
        audit.add("Evidence row", f"{row.EvidenceID}-boundary", len(str(row.Boundary)) > 40, row.Boundary)
        for rung in RUNGS:
            value = getattr(row, rung.replace(" ", "_")) if False else evidence.loc[evidence["EvidenceID"].eq(row.EvidenceID), rung].iloc[0]
            audit.add("Evidence coverage", f"{row.EvidenceID}-{rungs.loc[rungs['Rung'].eq(rung), 'RungOrder'].iloc[0]}", value in {"Direct", "Supporting", "Not addressed"}, value)

    indexed = evidence.set_index("EvidenceID")
    exact = {
        "E01": ("C1", "vannood2013fmt", 0, "Direct"),
        "E02": ("C1", "feuerstadt2022ser109", 0, "Direct"),
        "E03": ("C1", "smillie2018engraftment", 1, "Supporting"),
        "E04": ("C1", "li2016fmtstrains", 1, "Supporting"),
        "E05": ("C2", "buffie2015cscindens", 0, "Direct"),
        "E06": ("C2", "theriot2014metabolome", 0, "Supporting"),
        "E07": ("C3", "wirbel2019crc", 0, "Supporting"),
        "E08": ("C3", "kostic2013fusobacterium", 0, "Direct"),
        "E09": ("C3", "rubinstein2013fada", 0, "Direct"),
        "E10": ("C3", "bullman2017fusobacterium", 1, "Direct"),
    }
    for evidence_id, values in exact.items():
        row = indexed.loc[evidence_id]
        observed = (row.ClaimID, row.CitationKey, int(row.ShotgunStrainResolved), row.EvidenceRole)
        audit.add("Evidence identity", evidence_id, observed == values, observed)
    audit.add("Evidence boundary", "human-rct-c1", (indexed.loc[["E01", "E02"], "Randomized human intervention"] == "Direct").all(), indexed.loc[["E01", "E02"], "Randomized human intervention"].tolist())
    audit.add("Evidence boundary", "no-crc-human-rct", (indexed.loc[["E07", "E08", "E09", "E10"], "Randomized human intervention"] == "Not addressed").all(), "C3")
    audit.add("Evidence boundary", "no-human-bridge-borrowing", (evidence["Human causal bridge"] == "Not addressed").all(), evidence["Human causal bridge"].value_counts().to_dict())

    coverage = pd.read_csv(frozen / "study-domain-coverage.tsv", sep="\t")
    audit.add("Coverage", "row-count", len(coverage) == 70, len(coverage))
    audit.add("Coverage", "cell-unique", not coverage.duplicated(["EvidenceID", "Rung"]).any(), "unique")
    audit.add("Coverage", "all-cells", set(zip(coverage["EvidenceID"], coverage["Rung"])) == {(e, r) for e in evidence["EvidenceID"] for r in RUNGS}, "10 x 7")
    mapping = {"Not addressed": 0, "Supporting": 1, "Direct": 2}
    audit.add("Coverage", "codes", set(coverage["Coverage"]) == set(mapping), coverage["Coverage"].value_counts().to_dict())
    audit.add("Coverage", "score-map", coverage.apply(lambda row: mapping.get(row.Coverage) == int(row.CoverageScore), axis=1).all(), "0/1/2")
    audit.add("Coverage", "order-map", coverage.merge(rungs[["Rung", "RungOrder"]].astype({"RungOrder": int}), on="Rung", suffixes=("", "_expected")).eval("RungOrder == RungOrder_expected").all(), "1..7")
    counts = coverage["Coverage"].value_counts().to_dict()
    audit.add("Coverage", "direct-count", counts.get("Direct") == 25, counts)
    audit.add("Coverage", "supporting-count", counts.get("Supporting") == 13, counts)
    audit.add("Coverage", "not-addressed-count", counts.get("Not addressed") == 32, counts)
    for row in coverage.itertuples(index=False):
        expected_value = evidence.set_index("EvidenceID").loc[row.EvidenceID, row.Rung]
        audit.add("Coverage cell", f"{row.EvidenceID}-{row.RungOrder}", row.Coverage == expected_value, row.Coverage)


def audit_outcomes_and_packets(frozen: Path, audit: Audit) -> None:
    outcomes = pd.read_csv(frozen / "human-intervention-outcomes.tsv", sep="\t")
    audit.add("Trials", "row-count", len(outcomes) == 5, len(outcomes))
    expected = [
        ("van Nood 2013", "Donor FMT", 15, 16, 10),
        ("van Nood 2013", "Vancomycin", 4, 13, 10),
        ("van Nood 2013", "Vancomycin + lavage", 3, 13, 10),
        ("Feuerstadt 2022", "SER-109", 78, 89, 8),
        ("Feuerstadt 2022", "Placebo", 56, 93, 8),
    ]
    for index, (study, arm, events, total, weeks) in enumerate(expected):
        row = outcomes.iloc[index]
        audit.add("Trial arm", f"{index + 1}-identity", (row.Study, row.Arm) == (study, arm), (row.Study, row.Arm))
        audit.add("Trial arm", f"{index + 1}-counts", (int(row.FavorableEvents), int(row.Total)) == (events, total), (row.FavorableEvents, row.Total))
        audit.add("Trial arm", f"{index + 1}-time", int(row.TimeWeeks) == weeks, row.TimeWeeks)
        audit.add("Trial arm", f"{index + 1}-rate", near(row.FavorableRate, events / total), row.FavorableRate)
        audit.add("Trial arm", f"{index + 1}-valid", 0 <= events <= total, f"{events}/{total}")
    audit.add("Trials", "ser109-recurrence", (89 - int(outcomes.iloc[3].FavorableEvents), 93 - int(outcomes.iloc[4].FavorableEvents)) == (11, 37), "11/89 vs 37/93")
    raw_rr = (11 / 89) / (37 / 93)
    audit.add("Trials", "ser109-raw-rr-not-published-model", near(raw_rr, 0.3106589736, 1e-8) and not near(raw_rr, 0.32, 1e-4), raw_rr)
    audit.add("Trials", "van-nood-overall", outcomes.iloc[:3]["FavorableEvents"].astype(int).tolist() == [15, 4, 3] and outcomes.iloc[:3]["Total"].astype(int).tolist() == [16, 13, 13], "15/16 vs 4/13 and 3/13")

    downgrades = pd.read_csv(frozen / "claim-downgrade-examples.tsv", sep="\t", dtype=str).fillna("")
    audit.add("Downgrades", "row-count", len(downgrades) == 5, len(downgrades))
    for index, row in enumerate(downgrades.itertuples(index=False), start=1):
        audit.add("Downgrade", f"{index}-overclaim", len(row.Overclaim) > 25, row.Overclaim)
        audit.add("Downgrade", f"{index}-supported", len(row.SupportedClaim) > 70, row.SupportedClaim)
        audit.add("Downgrade", f"{index}-reason", len(row.WhyDowngraded) > 45, row.WhyDowngraded)

    packet = pd.read_csv(frozen / "evidence-packet.tsv", sep="\t", dtype=str).fillna("")
    audit.add("Packet", "row-count", len(packet) == 8, len(packet))
    audit.add("Packet", "order", packet["Order"].astype(int).tolist() == list(range(1, 9)), packet["Order"].tolist())
    expected_packets = (
        "Claim contract",
        "Human lineage",
        "Bias ledger",
        "Perturbation ledger",
        "Entity ledger",
        "Transport ledger",
        "Contradiction ledger",
        "Reproducibility bundle",
    )
    audit.add("Packet", "names", tuple(packet["Packet"]) == expected_packets, packet["Packet"].tolist())
    for row in packet.itertuples(index=False):
        audit.add("Packet item", f"{row.Order}-fields", len(row.RequiredFields) > 30, row.RequiredFields)
        audit.add("Packet item", f"{row.Order}-use", len(row.UsedFor) > 10, row.UsedFor)


def audit_metrics_and_environment(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "analysis-metrics.json").read_text())
    expected = {
        "article": 72,
        "analysis_seed": 72001,
        "plot_seed": 20260772,
        "claim_count": 3,
        "evidence_rows": 10,
        "publication_count": 11,
        "rung_count": 7,
        "direct_cells": 25,
        "supporting_cells": 13,
        "not_addressed_cells": 32,
        "shotgun_strain_resolved_studies": 3,
        "human_randomized_studies": 2,
        "crc_human_randomized_targeted_interventions": 0,
    }
    for field, value in expected.items():
        audit.add("Metrics", field, metrics.get(field) == value, metrics.get(field))
    audit.add("Metrics", "python-recorded", bool(re.fullmatch(r"\d+\.\d+\.\d+", str(metrics.get("python", "")))), metrics.get("python"))
    audit.add("Metrics", "pandas-recorded", bool(re.fullmatch(r"\d+\.\d+\.\d+", str(metrics.get("pandas", "")))), metrics.get("pandas"))

    python_env = (frozen / "env" / "multiomics-python.yml").read_text()
    pins = (
        "python=3.10.19",
        "numpy=2.2.6",
        "pandas=2.3.3",
        "pillow=12.2.0",
        "matplotlib==3.10.8",
    )
    for pin in pins:
        audit.add("Environment", pin, pin in python_env, pin)
    r_env = pd.read_csv(frozen / "env" / "multiomics-r-packages.tsv", sep="\t", dtype=str)
    for package in ("R", "ggplot2", "jsonlite"):
        audit.add("Environment", f"r-{package}", package in set(r_env["Package"]), package)


def audit_chapter(root: Path, audit: Audit) -> None:
    chapter = root / "chapters" / "72-causal-evidence-ladder.qmd"
    audit.add("Chapter", "exists", chapter.is_file(), chapter)
    if not chapter.is_file():
        return
    text = chapter.read_text(encoding="utf-8")
    audit.add("Chapter", "published", "draft: false" in text, "draft: false")
    audit.add("Chapter", "eval-true", "eval: true" in text, "eval: true")
    audit.add("Chapter", "freeze-auto", "freeze: auto" in text, "freeze: auto")
    audit.add("Chapter", "native-fences", text.count("```{r}") >= 7 and text.count("```{bash}") >= 2 and "~~~{" not in text, "7 R + 2 bash minimum")
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
        "不是穷尽检索的 systematic review",
        "不能拼成一个未被任何研究直接测试的 treatment effect",
        "CoverageScore=0/1/2",
        "绝不能求和",
        "随机化的是 FMT/product assignment",
        "0 targeted human CRC randomized trials",
        "can promote CRC phenotypes in specified models",
        "not a systematic review",
        "not summed or interpreted as quality weights",
    )
    for index, phrase in enumerate(phrases, start=1):
        audit.add("Chapter boundary", f"boundary-{index}", phrase in text, phrase)
    audit.add("Chapter", "seeds", "72001" in text and "set.seed(20260772)" in text, "72001 / 20260772")
    audit.add("Chapter", "inline-theme", "theme_pub <- function" in text and "save_pub <- function" in text and "scale_color_pub <- function" in text and "scale_fill_pub  <- function" in text, "inline plotting functions")
    audit.add("Chapter", "no-source-theme", 'source("R/theme_pub.R")' not in text, "no source dependency")
    audit.add("Chapter", "frozen-input", "data/small/72-causal-evidence-frozen" in text, "frozen bundle")
    audit.add("Chapter", "checksum-command", "sha256sum -c" in text, "checksum verification")
    audit.add("Chapter", "versions", all(value in text for value in ("Python 3.13.2", "pandas 3.0.3", "plot seed 20260772")), "runtime versions")
    forbidden = (
        "本篇可独立跑通",
        "这体现全系列",
        "接口只学一次",
        "作者代码通常长这样",
        "（即本文）",
        "/media/desk16",
        "/tmp/article72",
        "results/article72/evidence/qa_report",
    )
    for phrase in forbidden:
        audit.add("Chapter prose", f"forbidden-{phrase}", phrase not in text, phrase)
    for stem in FIGURES:
        audit.add("Chapter figure", stem, f"../figures/{stem}.png" in text, stem)
    audit.add("Chapter figure", "anchor", "../figures/72-buffie-figure4-original.jpg" in text, "official anchor")
    for key in tuple(PUBLICATIONS) + ("imai2010general", "sanderson2022mr"):
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

    anchor = figures / "72-buffie-figure4-original.jpg"
    frozen_anchor = frozen / "buffie-figure4-original.jpg"
    audit.add("Figure file", "anchor", anchor.is_file() and frozen_anchor.is_file() and sha256(anchor) == sha256(frozen_anchor) == ANCHOR_SHA256, sha256(anchor) if anchor.is_file() else "MISSING")
    if anchor.is_file():
        with Image.open(anchor) as image:
            audit.add("Figure raster", "anchor-dimensions", image.size == (1578, 707), image.size)

    with tempfile.TemporaryDirectory(prefix="article72-replot-") as temporary:
        staged = Path(temporary) / "figures"
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(Path(temporary) / "mpl")
        result = subprocess.run(
            [
                sys.executable,
                str(frozen / "scripts" / "plot_article72_evidence.py"),
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
    publications = audit_publications(frozen, audit)
    _, rungs = audit_contracts(frozen, audit)
    audit_evidence(frozen, publications, rungs, audit)
    audit_outcomes_and_packets(frozen, audit)
    audit_metrics_and_environment(frozen, audit)
    audit_chapter(root, audit)
    audit_figures(root, frozen, audit)
    report = audit.frame()
    report.to_csv(qa / "qa_report.tsv", sep="\t", index=False)
    passed = int(report["Status"].eq("PASS").sum())
    failed = int(report["Status"].eq("FAIL").sum())
    payload = {
        "article": 72,
        "status": "passed" if failed == 0 else "failed",
        "checks": len(report),
        "passed": passed,
        "failed": failed,
        "failed_checks": report.loc[
            report["Status"].eq("FAIL"), ["Category", "Check", "Detail"]
        ].to_dict("records"),
    }
    (qa / "qa_report.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
