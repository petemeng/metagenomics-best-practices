#!/usr/bin/env python3
"""Fail-closed validation for Article 66 MAG-metabolite evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from article42_44_validation_utils import Audit, audit_checksums, audit_figures, finish


FIGURES = (
    "66-mag-quality-landscape",
    "66-input-normalization-audit",
    "66-individual-gem-gapfill",
    "66-metabolite-overlap",
    "66-validation-boundary",
    "66-pathway-robustness",
)
EXPECTED_SOURCES = {
    "paper.xml": (
        122_202,
        "6b3fcd4db9bbb0a3c6dd7409db3a70ee7fa5154d04a06d2112626281a34eccec",
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11406951/fullTextXML",
    ),
    "supplementary-files.zip": (
        955_771,
        "49c674c1bad788ea15e3d48ed8b2909b0d67c22e318276103abe4b0d307f28ae",
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11406951/supplementaryFiles",
    ),
}
EXPECTED_MEMBERS = {
    "supplementary-tables.xlsx": (
        101_735,
        "c740876ed3dce897ac9105a5434a88792afa8bc6bc072ae588e5cc1173f6d416",
    ),
    "majzoub-fig2-original.jpg": (
        41_669,
        "e96f19623d9337a6cb41c1ea0294b13ea69426df1166b88cf19fa55c661a284b",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stage-dir", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    return parser.parse_args()


def near(value: object, expected: float, tolerance: float = 1e-10) -> bool:
    try:
        return math.isclose(float(value), expected, rel_tol=tolerance, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pixel_sha(path: Path) -> str:
    with Image.open(path) as image:
        normalized = image.convert("RGBA")
        payload = (
            normalized.size[0].to_bytes(8, "little")
            + normalized.size[1].to_bytes(8, "little")
            + normalized.tobytes()
        )
    return hashlib.sha256(payload).hexdigest()


def audit_chapter(chapter: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    sections = (
        "这一步对应论文里的哪张图",
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
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "bibliography": "bibliography: ../references.bib" in text,
        "analysis-seed": "66001" in text,
        "plot-seed": "20260766" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "hardware-contract": all(token in text for token in ("RAM", "CPU", "磁盘", "耗时")),
        "real-study": all(token in text for token in ("Majzoub", "10.1128/msystems.00746-24", "PMC11406951", "PRJEB50699")),
        "mag-contract": all(token in text for token in ("170", "56", "67", "34", "37", "MIMAG", "rRNA", "tRNA")),
        "metabolomics-contract": all(token in text for token in ("1,348", "1,190", "523", "488", "presence/absence")),
        "overlap-contract": all(token in text for token in ("1,475", "1,530", "0.941", "0.950", "0.957")),
        "pathway-contract": all(token in text for token in ("141", "147", "propanoate", "SMP00016")),
        "phenotype-boundary": all(token in text for token in ("100%", "36.4%", "15", "0.026", "两个供体")),
        "validation-boundary": all(token in text for token in ("sensitivity", "specificity", "accuracy", "flux", "causal")),
        "software-contract": all(token in text for token in ("MetaWRAP v1.3.2", "dRep v2.3.2", "GTDB-Tk v1.5.1", "GTDB r207", "RASTtk v1.073")),
        "methods-template": "**MAG-guided metabolic modelling and metabolomics validation.**" in text,
        "results-template": "**Results template.**" in text,
        "frozen-input": "data/small/66-mag-metabolite-frozen" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "anchor": "../figures/66-majzoub-fig2-original.jpg" in text,
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-meta-prose": not any(token in text for token in ("本篇可独立", "本文可独立", "全系列约定", "接口只学一次", "作者代码通常长这样", "（即本文）", "无头服务器")),
    }
    for section in sections:
        checks[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def stage_figures(root: Path, frozen: Path, stage: Path, python: Path, final: Path, audit: Audit) -> None:
    if stage.exists():
        shutil.rmtree(stage)
    figures = stage / "figures"
    environment = os.environ.copy()
    environment.update(
        {
            "MPLCONFIGDIR": str(stage / "matplotlib"),
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    result = subprocess.run(
        [
            str(python),
            str(frozen / "scripts/plot_article66_mag_metabolite.py"),
            "--input-dir",
            str(frozen),
            "--figure-dir",
            str(figures),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        timeout=240,
    )
    audit.add("Reanalysis", "plot-script-exit", result.returncode == 0, result.stdout + result.stderr)
    for stem in FIGURES:
        staged = figures / f"{stem}.png"
        published = final / f"{stem}.png"
        status = staged.is_file() and published.is_file() and pixel_sha(staged) == pixel_sha(published)
        audit.add("Reanalysis", f"{stem}-pixel-identical", status, pixel_sha(staged) if staged.is_file() else "MISSING")
    staged_anchor = figures / "66-majzoub-fig2-original.jpg"
    final_anchor = final / "66-majzoub-fig2-original.jpg"
    audit.add("Reanalysis", "anchor-byte-identical", staged_anchor.is_file() and final_anchor.is_file() and sha256(staged_anchor) == sha256(final_anchor), sha256(staged_anchor) if staged_anchor.is_file() else "MISSING")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_checksums(frozen, audit)

    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    expected_metrics = {
        "article": 66,
        "mag_total": 170,
        "paper_high_total": 123,
        "paper_medium_total": 47,
        "normalized_mag_total": 71,
        "species_resolved_total": 166,
        "individual_model_rows": 142,
        "community_model_configurations": 14,
        "metabolomics_longitudinal_samples": 11,
        "within_donor_pathway_records": 302,
        "between_donor_pathway_records": 18,
        "combined_pathway_records": 297,
        "combined_donor_specific_pathways": 5,
        "negative_gapfill_deltas": 1,
        "missing_blocked_reaction_counts": 2,
        "seed": 66001,
        "plot_seed": 20260766,
    }
    for key, expected in expected_metrics.items():
        audit.add("Metric", key, metrics.get(key) == expected, metrics.get(key))
    expected_donors = {
        "Donor 1": (86, 56, 30, 34, 1500, 1493, 28583, 28172),
        "Donor 2": (84, 67, 17, 37, 1562, 1545, 33661, 31902),
    }
    for donor, expected in expected_donors.items():
        record = metrics.get("donors", {}).get(donor, {})
        observed = tuple(record.get(key) for key in ("mags", "paper_high", "paper_medium", "normalized_mags", "reference_unique_normalized", "mag_unique_normalized", "reference_total_normalized", "mag_total_normalized"))
        audit.add("Metric", donor, observed == expected, observed)
    expected_jaccard = {
        "Donor 1 · Predicted": 1475 / 1568,
        "Donor 2 · Predicted": 1530 / 1611,
        "Donor 1 · Confirmed": 180 / 188,
        "Donor 2 · Confirmed": 176 / 184,
    }
    for key, expected in expected_jaccard.items():
        value = metrics.get("overlap_jaccard", {}).get(key)
        audit.add("Metric", key, near(value, expected), value)

    manifest = json.loads((frozen / "source-manifest.json").read_text(encoding="utf-8"))
    audit.add("Source", "doi", manifest.get("doi") == "10.1128/msystems.00746-24", manifest.get("doi"))
    audit.add("Source", "pmcid", manifest.get("pmcid") == "PMC11406951", manifest.get("pmcid"))
    audit.add("Source", "ena", manifest.get("human_metagenomics_accession") == "PRJEB50699", manifest.get("human_metagenomics_accession"))
    for name, (size, digest, url) in EXPECTED_SOURCES.items():
        record = manifest.get("resources", {}).get(name, {})
        audit.add("Source", f"{name}-bytes", record.get("bytes") == size, record.get("bytes"))
        audit.add("Source", f"{name}-sha256", record.get("sha256") == digest, record.get("sha256"))
        audit.add("Source", f"{name}-official-url", record.get("url") == url, record.get("url"))
    for name, (size, digest) in EXPECTED_MEMBERS.items():
        record = manifest.get("extracted_members", {}).get(name, {})
        audit.add("Source", f"{name}-bytes", record.get("bytes") == size, record.get("bytes"))
        audit.add("Source", f"{name}-sha256", record.get("sha256") == digest, record.get("sha256"))

    ledger = pd.read_csv(frozen / "mag-ledger.tsv", sep="\t")
    audit.add("MAG", "rows", len(ledger) == 170, len(ledger))
    audit.add("MAG", "unique-within-donor", not ledger.duplicated(["Donor", "MAGID"]).any(), int(ledger.duplicated(["Donor", "MAGID"]).sum()))
    audit.add("MAG", "all-paper-gated", set(ledger["PaperQuality"]) == {"Paper high", "Paper medium"}, sorted(ledger["PaperQuality"].unique()))
    audit.add("MAG", "mimag-not-assessable", ledger["MIMAGHighQualityStatus"].eq("Not assessable: rRNA/tRNA evidence absent").all(), ledger["MIMAGHighQualityStatus"].value_counts().to_dict())
    audit.add("MAG", "normalized-high-only", ledger.loc[ledger["NormalizedInput"].astype(bool), "PaperQuality"].eq("Paper high").all(), int(ledger["NormalizedInput"].astype(bool).sum()))

    models = pd.read_csv(frozen / "individual-model-audit.tsv", sep="\t")
    audit.add("Model", "rows", len(models) == 142, len(models))
    audit.add("Model", "finite-core", np.isfinite(models[["CompoundsBeforeGapfill", "CompoundsAfterGapfill", "TotalReactions", "BlockedReactionPercent", "GapfillAddedCompounds"]].to_numpy(float)).all(), "finite")
    invalid_gapfill = ~models["GapfillDeltaValid"].astype(bool)
    audit.add(
        "Model",
        "one-published-negative-gapfill-delta",
        invalid_gapfill.sum() == 1 and near(models["GapfillAddedCompounds"].min(), -59),
        models.loc[invalid_gapfill, ["Donor", "ModelID", "GapfillAddedCompounds"]].to_dict("records"),
    )
    audit.add(
        "Model",
        "valid-gapfill-deltas-nonnegative",
        models.loc[~invalid_gapfill, "GapfillAddedCompounds"].ge(0).all(),
        float(models.loc[~invalid_gapfill, "GapfillAddedCompounds"].min()),
    )
    audit.add(
        "Model",
        "two-published-missing-blocked-counts",
        models["BlockedReactions"].isna().sum() == 2,
        int(models["BlockedReactions"].isna().sum()),
    )
    anomalies = pd.read_csv(frozen / "source-anomaly-ledger.tsv", sep="\t")
    audit.add(
        "Model",
        "source-anomaly-ledger",
        len(anomalies) == 3 and anomalies["ModelID"].nunique() == 3,
        anomalies.to_dict("records"),
    )
    audit.add("Model", "mag-quality-mapped", models.loc[models["Approach"].eq("MAG-guided"), "Completeness"].notna().all(), int(models.loc[models["Approach"].eq("MAG-guided"), "Completeness"].notna().sum()))

    community = pd.read_csv(frozen / "community-model-audit.tsv", sep="\t")
    audit.add("Community", "rows", len(community) == 14, len(community))
    audit.add("Community", "all-positive", (community.select_dtypes(include=[np.number]) > 0).all().all(), community.min(numeric_only=True).to_dict())
    overlap = pd.read_csv(frozen / "metabolite-overlap-audit.tsv", sep="\t")
    audit.add("Overlap", "four-panels", len(overlap) == 4, len(overlap))
    audit.add("Overlap", "arithmetic", ((overlap["ReferenceOnly"] + overlap["Shared"] == overlap["ReferenceTotal"]) & (overlap["MAGOnly"] + overlap["Shared"] == overlap["MAGTotal"]) & (overlap[["ReferenceOnly", "Shared", "MAGOnly"]].sum(axis=1) == overlap["Union"])).all(), overlap.to_dict("records"))

    coverage = pd.read_csv(frozen / "metabolomics-coverage.tsv", sep="\t")
    observed_coverage = [
        tuple(int(value) for value in row)
        for row in coverage[
            [
                "LongitudinalSamples",
                "DetectedMetabolites",
                "MetabolitesWithKEGGID",
                "KEGGIDsSearched",
            ]
        ].to_numpy()
    ]
    audit.add("Metabolomics", "coverage", observed_coverage == [(4, 1348, 489, 523), (7, 1190, 453, 488)], observed_coverage)
    phenotype = pd.read_csv(frozen / "phenotype-evidence.tsv", sep="\t")
    audit.add("Phenotype", "two-donor-boundary", len(phenotype) == 2 and phenotype["IndependentDonorUnits"].eq(2).all() and phenotype["RecipientsAcrossComparison"].eq(15).all(), phenotype.to_dict("records"))

    counts = pd.read_csv(frozen / "pathway-counts.tsv", sep="\t")
    expected_count_rows = {
        ("Donor 1", "Common", "Within donor · separate approaches"): 141,
        ("Donor 2", "Common", "Within donor · separate approaches"): 147,
        ("Reference-guided", "Present in donor 1", "Between donors · separate approaches"): 5,
        ("Reference-guided", "Present in donor 2", "Between donors · separate approaches"): 7,
        ("MAG-guided", "Present in donor 1", "Between donors · separate approaches"): 1,
        ("MAG-guided", "Present in donor 2", "Between donors · separate approaches"): 5,
    }
    indexed = counts.set_index(["Parent", "Category", "Comparison"])["Count"]
    for key, expected in expected_count_rows.items():
        value = int(indexed.loc[key]) if key in indexed.index else -1
        audit.add("Pathway", " · ".join(key), value == expected, value)
    donor_pathways = pd.read_csv(frozen / "combined-donor-pathways.tsv", sep="\t")
    propanoate = donor_pathways[donor_pathways["PathwayID"].eq("SMP00016")]
    audit.add("Pathway", "propanoate", len(propanoate) == 1 and propanoate.iloc[0]["Donor"] == "Donor 2" and near(propanoate.iloc[0]["QValue"], 0.0071), propanoate.to_dict("records"))

    anchor = frozen / "majzoub-fig2-original.jpg"
    audit.add("Anchor", "source-checksum", anchor.is_file() and sha256(anchor) == EXPECTED_MEMBERS["majzoub-fig2-original.jpg"][1], sha256(anchor) if anchor.is_file() else "MISSING")
    audit_chapter(args.chapter.resolve(), audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    stage_figures(root, frozen, args.stage_dir.resolve(), args.python.resolve(), args.figure_dir.resolve(), audit)
    return_code = finish(
        article=66,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "study": "Majzoub et al. mSystems 2024",
            "mag_rows": len(ledger),
            "individual_model_rows": len(models),
            "independent_donor_units": 2,
            "figure_sets": len(FIGURES),
        },
    )
    raise SystemExit(return_code)


if __name__ == "__main__":
    main()
