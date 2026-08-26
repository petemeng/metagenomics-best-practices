#!/usr/bin/env python3
"""Prepare MAG-to-metabolite evidence from Majzoub et al. for Article 66."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import re
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr


SEED = 66_001
PLOT_SEED = 20_260_766
APPROACHES = (
    "Reference ≥5%",
    "Reference ≥2.5%",
    "Reference ≥1%",
    "Reference ≥0.5%",
    "MAG high",
    "MAG high + medium",
    "MAG normalized",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def clean_text(element: ET.Element) -> str:
    return " ".join("".join(element.itertext()).split())


def number(value: str) -> float:
    return float(value.replace(",", ""))


def parse_table1(root: ET.Element) -> pd.DataFrame:
    table = root.find(".//table-wrap[@id='T1']/table/tbody")
    if table is None:
        raise RuntimeError("Article Table 1 was not found in the full-text XML")
    records: dict[tuple[str, str], dict[str, float]] = {}
    donor = ""
    for row in table.findall("./tr"):
        cells = [clean_text(cell) for cell in row.findall("./*")]
        if cells and cells[0].startswith("Donor"):
            donor = cells.pop(0)
        if len(cells) != 8 or not donor:
            raise RuntimeError(f"Unexpected Table 1 row: {cells}")
        parameter, *values = cells
        for approach, value in zip(APPROACHES, values, strict=True):
            records.setdefault((donor, approach), {})[parameter] = number(value)
    frame = pd.DataFrame(
        [
            {"Donor": donor, "Approach": approach, **values}
            for (donor, approach), values in records.items()
        ]
    )
    expected = {
        "Taxa",
        "Total predicted metabolites",
        "Total extracellular metabolites",
        "Unique metabolites",
        "Confirmed unique metabolites (total)",
        "Data loss (%)",
        "Confirmed unique metabolites (%)",
    }
    if len(frame) != 14 or not expected.issubset(frame.columns):
        raise RuntimeError("Table 1 parsing contract failed")
    return frame


def parse_table2(root: ET.Element) -> pd.DataFrame:
    table = root.find(".//table-wrap[@id='T2']/table/tbody")
    if table is None:
        raise RuntimeError("Article Table 2 was not found in the full-text XML")
    records: list[dict[str, object]] = []
    donor = ""
    for row in table.findall("./tr"):
        cells = [clean_text(cell) for cell in row.findall("./*")]
        if cells and cells[0].startswith("Donor"):
            donor = cells.pop(0)
        if len(cells) != 4 or not donor:
            raise RuntimeError(f"Unexpected Table 2 row: {cells}")
        pathway_id, pathway_name, source, q_value = cells
        records.append(
            {
                "Donor": donor,
                "PathwayID": pathway_id,
                "PathwayName": pathway_name,
                "Source": source,
                "QValue": float(q_value),
            }
        )
    return pd.DataFrame(records)


def taxonomy_rank(taxonomy: str, prefix: str) -> str:
    match = re.search(rf"(?:^|;){re.escape(prefix)}([^;]*)", taxonomy)
    return match.group(1) if match else ""


def load_mag_ledger(workbook: Path, donor: int, sheet: str) -> pd.DataFrame:
    frame = pd.read_excel(workbook, sheet_name=sheet, header=4)
    frame = frame.dropna(subset=["Bin Id"]).copy()
    frame = frame.rename(
        columns={
            "Bin Id": "MAGID",
            "GTDBv207_classification": "GTDBr207Taxonomy",
        }
    )
    frame.insert(0, "Donor", f"Donor {donor}")
    frame["Completeness"] = pd.to_numeric(frame["Completeness"], errors="raise")
    frame["Contamination"] = pd.to_numeric(frame["Contamination"], errors="raise")
    high = (frame["Completeness"] > 90) & (frame["Contamination"] < 5)
    medium = (frame["Completeness"] >= 50) & (frame["Contamination"] < 10)
    if not medium.all():
        raise RuntimeError(f"Published donor {donor} MAG sheet contains an out-of-gate bin")
    frame["PaperQuality"] = np.where(high, "Paper high", "Paper medium")
    for column, prefix in (
        ("Domain", "d__"),
        ("Phylum", "p__"),
        ("Class", "c__"),
        ("Order", "o__"),
        ("Family", "f__"),
        ("Genus", "g__"),
        ("Species", "s__"),
    ):
        frame[column] = frame["GTDBr207Taxonomy"].astype(str).map(
            lambda value, rank=prefix: taxonomy_rank(value, rank)
        )
    frame["SpeciesResolved"] = frame["Species"].ne("")
    frame["MIMAGHighQualityStatus"] = "Not assessable: rRNA/tRNA evidence absent"
    return frame


def extract_individual_models(workbook: Path, mag_ledger: pd.DataFrame) -> pd.DataFrame:
    raw = pd.read_excel(workbook, sheet_name="4", header=None)
    records: list[dict[str, object]] = []
    for donor, rows in ((1, range(6, 40)), (2, range(45, 82))):
        for row_index in rows:
            for approach, offset in (("Reference-guided", 0), ("MAG-guided", 7)):
                values = raw.iloc[row_index, offset : offset + 6].tolist()
                records.append(
                    {
                        "Donor": f"Donor {donor}",
                        "Approach": approach,
                        "ModelID": str(values[0]),
                        "CompoundsBeforeGapfill": float(values[1]),
                        "CompoundsAfterGapfill": float(values[2]),
                        "BlockedReactions": pd.to_numeric(values[3], errors="coerce"),
                        "TotalReactions": float(values[4]),
                        "BlockedReactionPercent": float(values[5]),
                    }
                )
    frame = pd.DataFrame(records)
    frame["GapfillAddedCompounds"] = (
        frame["CompoundsAfterGapfill"] - frame["CompoundsBeforeGapfill"]
    )
    frame["GapfillDeltaValid"] = frame["GapfillAddedCompounds"].ge(0)
    frame["GapfillAddedFraction"] = np.where(
        frame["GapfillDeltaValid"],
        frame["GapfillAddedCompounds"] / frame["CompoundsAfterGapfill"],
        np.nan,
    )
    frame["SourceAnomaly"] = ""
    frame.loc[~frame["GapfillDeltaValid"], "SourceAnomaly"] = (
        "Published after-gapfill compound count is lower than before-gapfill count"
    )
    missing_blocked = frame["BlockedReactions"].isna()
    blocked_note = "Published blocked-reaction count is missing while percent is reported"
    frame.loc[missing_blocked, "SourceAnomaly"] = frame.loc[
        missing_blocked, "SourceAnomaly"
    ].map(lambda value: f"{value}; {blocked_note}".strip("; "))
    annotations = mag_ledger[
        ["Donor", "MAGID", "Completeness", "Contamination", "Species", "PaperQuality"]
    ].rename(columns={"MAGID": "ModelID"})
    frame = frame.merge(annotations, how="left", on=["Donor", "ModelID"])
    mag_rows = frame["Approach"].eq("MAG-guided")
    if frame.loc[mag_rows, "Completeness"].isna().any():
        raise RuntimeError("One or more normalized MAG model IDs did not map to the MAG ledger")
    return frame


def model_summary(frame: pd.DataFrame) -> pd.DataFrame:
    metrics = (
        "CompoundsBeforeGapfill",
        "CompoundsAfterGapfill",
        "GapfillAddedCompounds",
        "GapfillAddedFraction",
        "BlockedReactionPercent",
    )
    records: list[dict[str, object]] = []
    for (donor, approach), group in frame.groupby(["Donor", "Approach"], sort=False):
        record: dict[str, object] = {
            "Donor": donor,
            "Approach": approach,
            "Models": len(group),
        }
        for metric in metrics:
            if metric.startswith("Gapfill"):
                values = group.loc[group["GapfillDeltaValid"], metric].dropna()
            else:
                values = group[metric].dropna()
            record[f"{metric}Median"] = float(values.median())
            record[f"{metric}Q1"] = float(values.quantile(0.25))
            record[f"{metric}Q3"] = float(values.quantile(0.75))
        records.append(record)
    return pd.DataFrame(records)


def descriptive_model_comparison(summary: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    metrics = (
        "CompoundsBeforeGapfill",
        "CompoundsAfterGapfill",
        "GapfillAddedCompounds",
        "BlockedReactionPercent",
    )
    for donor in ("Donor 1", "Donor 2"):
        subset = summary[summary["Donor"].eq(donor)].set_index("Approach")
        for metric in metrics:
            reference = float(subset.loc["Reference-guided", f"{metric}Median"])
            mag = float(subset.loc["MAG-guided", f"{metric}Median"])
            records.append(
                {
                    "Donor": donor,
                    "Metric": metric,
                    "ReferenceMedian": reference,
                    "MAGMedian": mag,
                    "MAGMinusReference": mag - reference,
                    "MAGToReferenceRatio": mag / reference if reference else np.nan,
                    "Inference": "Descriptive only; models are nested within two donors",
                }
            )
    return pd.DataFrame(records)


def quality_function_correlations(frame: pd.DataFrame) -> pd.DataFrame:
    mag = frame[frame["Approach"].eq("MAG-guided")].copy()
    records: list[dict[str, object]] = []
    for donor, group in mag.groupby("Donor", sort=False):
        for metric in (
            "CompoundsBeforeGapfill",
            "GapfillAddedCompounds",
            "BlockedReactionPercent",
        ):
            selected = (
                group[group["GapfillDeltaValid"]]
                if metric == "GapfillAddedCompounds"
                else group
            )
            result = spearmanr(selected["Completeness"], selected[metric])
            records.append(
                {
                    "Donor": donor,
                    "Models": len(selected),
                    "Metric": metric,
                    "SpearmanRho": float(result.statistic),
                    "PValueNotReported": True,
                    "Interpretation": "Range-restricted, selected MAG set; descriptive association only",
                }
            )
    return pd.DataFrame(records)


def extract_gapfill_summary(workbook: Path) -> pd.DataFrame:
    raw = pd.read_excel(workbook, sheet_name="12", header=None)
    records: list[dict[str, object]] = []
    donor = ""
    for row_index in range(6, 14):
        if pd.notna(raw.iat[row_index, 0]):
            donor = str(raw.iat[row_index, 0])
        metric = str(raw.iat[row_index, 1])
        for approach, column in (("Reference ≥0.5%", 2), ("MAG normalized", 3)):
            records.append(
                {
                    "Donor": donor,
                    "Approach": approach,
                    "Metric": metric,
                    "Value": float(raw.iat[row_index, column]),
                }
            )
    return pd.DataFrame(records)


def extract_pathways(
    workbook: Path,
    sheet: str,
    left_label: str,
    right_label: str,
) -> pd.DataFrame:
    raw = pd.read_excel(workbook, sheet_name=sheet, header=None)
    records: list[dict[str, object]] = []
    for parent, offset in ((left_label, 0), (right_label, 5)):
        category = ""
        for row_index in range(5, len(raw)):
            marker = raw.iat[row_index, offset]
            pathway_id = raw.iat[row_index, offset + 1]
            if pd.notna(marker):
                category = str(marker).strip()
            if pd.isna(pathway_id):
                continue
            records.append(
                {
                    "Parent": parent,
                    "Category": category,
                    "PathwayID": str(pathway_id).strip(),
                    "PathwayName": str(raw.iat[row_index, offset + 2]).strip(),
                    "Source": str(raw.iat[row_index, offset + 3]).strip(),
                }
            )
    return pd.DataFrame(records)


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((cache / "download-manifest.json").read_text(encoding="utf-8"))
    if manifest.get("article") != 66:
        raise RuntimeError("Article 66 source manifest is missing or incompatible")
    shutil.copy2(cache / "download-manifest.json", output / "source-manifest.json")
    shutil.copy2(cache / "majzoub-fig2-original.jpg", output / "majzoub-fig2-original.jpg")

    workbook = cache / "supplementary-tables.xlsx"
    root = ET.parse(cache / "paper.xml").getroot()
    community = parse_table1(root)
    combined_donor_pathways = parse_table2(root)

    mag_ledger = pd.concat(
        [load_mag_ledger(workbook, 1, "2"), load_mag_ledger(workbook, 2, "3")],
        ignore_index=True,
    )
    individual = extract_individual_models(workbook, mag_ledger)
    normalized_ids = set(
        zip(
            individual.loc[individual["Approach"].eq("MAG-guided"), "Donor"],
            individual.loc[individual["Approach"].eq("MAG-guided"), "ModelID"],
        )
    )
    mag_ledger["NormalizedInput"] = [
        (donor, mag_id) in normalized_ids
        for donor, mag_id in zip(mag_ledger["Donor"], mag_ledger["MAGID"], strict=True)
    ]
    if not mag_ledger.loc[mag_ledger["NormalizedInput"], "PaperQuality"].eq("Paper high").all():
        raise RuntimeError("Normalized input unexpectedly contains a paper-medium MAG")

    quality_summary = (
        mag_ledger.groupby("Donor", sort=False)
        .agg(
            MAGs=("MAGID", "size"),
            PaperHigh=("PaperQuality", lambda x: int((x == "Paper high").sum())),
            PaperMedium=("PaperQuality", lambda x: int((x == "Paper medium").sum())),
            NormalizedInput=("NormalizedInput", "sum"),
            SpeciesResolved=("SpeciesResolved", "sum"),
            CompletenessMedian=("Completeness", "median"),
            CompletenessMinimum=("Completeness", "min"),
            ContaminationMedian=("Contamination", "median"),
            ContaminationMaximum=("Contamination", "max"),
        )
        .reset_index()
    )

    individual_summary = model_summary(individual)
    model_comparison = descriptive_model_comparison(individual_summary)
    correlations = quality_function_correlations(individual)
    gapfill = extract_gapfill_summary(workbook)

    overlaps = pd.DataFrame(
        [
            ("Donor 1", "Predicted", 25, 1475, 68),
            ("Donor 2", "Predicted", 32, 1530, 49),
            ("Donor 1", "Confirmed", 3, 180, 5),
            ("Donor 2", "Confirmed", 4, 176, 4),
        ],
        columns=["Donor", "Evidence", "ReferenceOnly", "Shared", "MAGOnly"],
    )
    overlaps["ReferenceTotal"] = overlaps["ReferenceOnly"] + overlaps["Shared"]
    overlaps["MAGTotal"] = overlaps["MAGOnly"] + overlaps["Shared"]
    overlaps["Union"] = overlaps[["ReferenceOnly", "Shared", "MAGOnly"]].sum(axis=1)
    overlaps["Jaccard"] = overlaps["Shared"] / overlaps["Union"]
    expected_overlap_totals = {
        ("Donor 1", "Predicted"): (1500, 1543),
        ("Donor 2", "Predicted"): (1562, 1579),
        ("Donor 1", "Confirmed"): (183, 185),
        ("Donor 2", "Confirmed"): (180, 180),
    }
    for row in overlaps.itertuples(index=False):
        if (row.ReferenceTotal, row.MAGTotal) != expected_overlap_totals[(row.Donor, row.Evidence)]:
            raise RuntimeError("Figure 2 overlap arithmetic failed")

    metabolomics = pd.DataFrame(
        [
            ("Donor 1", 4, "24;26;34;40", 1348, 489, 34, 523),
            ("Donor 2", 7, "20;22;31;33;41;42;70", 1190, 453, 35, 488),
        ],
        columns=[
            "Donor",
            "LongitudinalSamples",
            "Weeks",
            "DetectedMetabolites",
            "MetabolitesWithKEGGID",
            "AdditionalSynonymousKEGGIDs",
            "KEGGIDsSearched",
        ],
    )
    metabolomics["KEGGAnnotatedFraction"] = (
        metabolomics["MetabolitesWithKEGGID"] / metabolomics["DetectedMetabolites"]
    )
    metabolomics["ValidationTransform"] = "Presence/absence per donor before imputation"
    metabolomics["BiologicalOriginBoundary"] = "Microbial + host + diet + other biota"

    pathways_within = extract_pathways(workbook, "5", "Donor 1", "Donor 2")
    pathways_between = extract_pathways(
        workbook, "6", "Reference-guided", "MAG-guided"
    )
    pathways_combined = extract_pathways(workbook, "7", "Donor 1", "Donor 2")
    pathway_counts = pd.concat(
        [
            pathways_within.groupby(["Parent", "Category"]).size().rename("Count").reset_index().assign(Comparison="Within donor · separate approaches"),
            pathways_between.groupby(["Parent", "Category"]).size().rename("Count").reset_index().assign(Comparison="Between donors · separate approaches"),
            pathways_combined.groupby(["Parent", "Category"]).size().rename("Count").reset_index().assign(Comparison="Within donor · combined output"),
        ],
        ignore_index=True,
    )
    phenotype = pd.DataFrame(
        [
            ("Donor 1 microbiota", 100.0),
            ("Donor 2 microbiota", 36.4),
        ],
        columns=["DonorExposure", "ReportedUCRecipientEfficacyPercent"],
    )
    phenotype["RecipientsAcrossComparison"] = 15
    phenotype["PublishedFishersP"] = 0.026
    phenotype["IndependentDonorUnits"] = 2
    phenotype["InferenceBoundary"] = (
        "Descriptive donor contrast; insufficient donor units for generalizable MAG-pathway-phenotype inference"
    )

    methods_contract = {
        "human_reads": "ENA PRJEB50699",
        "binning": "MetaBAT2, MaxBin2 and CONCOCT in MetaWRAP 1.3.2",
        "mag_quality": "CheckM; paper high >90% completeness and <5% contamination; paper medium >=50% completeness and <10% contamination",
        "dereplication": "dRep 2.3.2 at 99% ANI",
        "taxonomy": "GTDB-Tk 1.5.1 with GTDB r207",
        "annotation": "RASTtk 1.073, similarity e-value cutoff 1e-6",
        "gem": "KBase ModelSEED/OMEGGA applications",
        "gapfill": "MS2 Improved Gapfill Metabolic Models with OMEGGA",
        "metabolomics": "Metabolon Precision Metabolomics global LC-MS; presence/absence per donor before imputation",
        "enrichment": "MBROLE 2.0; SMPDB, KEGG and UniPathway; Benjamini-Hochberg q<0.05",
        "mimag_boundary": "The published high-quality label cannot be upgraded to MIMAG high-quality without rRNA/tRNA evidence",
    }
    (output / "methods-contract.json").write_text(
        json.dumps(methods_contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    outputs = {
        "community-model-audit.tsv": community,
        "combined-donor-pathways.tsv": combined_donor_pathways,
        "mag-ledger.tsv": mag_ledger,
        "mag-quality-summary.tsv": quality_summary,
        "individual-model-audit.tsv": individual,
        "source-anomaly-ledger.tsv": individual.loc[
            individual["SourceAnomaly"].ne(""),
            [
                "Donor",
                "Approach",
                "ModelID",
                "CompoundsBeforeGapfill",
                "CompoundsAfterGapfill",
                "GapfillAddedCompounds",
                "BlockedReactions",
                "BlockedReactionPercent",
                "SourceAnomaly",
            ],
        ],
        "individual-model-summary.tsv": individual_summary,
        "model-comparison-summary.tsv": model_comparison,
        "quality-function-correlations.tsv": correlations,
        "community-gapfill-audit.tsv": gapfill,
        "metabolite-overlap-audit.tsv": overlaps,
        "metabolomics-coverage.tsv": metabolomics,
        "pathways-within-donor.tsv": pathways_within,
        "pathways-between-donors.tsv": pathways_between,
        "pathways-combined.tsv": pathways_combined,
        "pathway-counts.tsv": pathway_counts,
        "phenotype-evidence.tsv": phenotype,
    }
    for filename, frame in outputs.items():
        frame.to_csv(output / filename, sep="\t", index=False)

    normalized = community[community["Approach"].isin(["Reference ≥0.5%", "MAG normalized"])]
    normalized_pivot = normalized.pivot(index="Donor", columns="Approach")
    metrics = {
        "article": 66,
        "seed": SEED,
        "plot_seed": PLOT_SEED,
        "mag_total": int(len(mag_ledger)),
        "paper_high_total": int(mag_ledger["PaperQuality"].eq("Paper high").sum()),
        "paper_medium_total": int(mag_ledger["PaperQuality"].eq("Paper medium").sum()),
        "normalized_mag_total": int(mag_ledger["NormalizedInput"].sum()),
        "species_resolved_total": int(mag_ledger["SpeciesResolved"].sum()),
        "individual_model_rows": int(len(individual)),
        "community_model_configurations": int(len(community)),
        "metabolomics_longitudinal_samples": int(metabolomics["LongitudinalSamples"].sum()),
        "within_donor_pathway_records": int(len(pathways_within)),
        "between_donor_pathway_records": int(len(pathways_between)),
        "combined_pathway_records": int(len(pathways_combined)),
        "combined_donor_specific_pathways": int(len(combined_donor_pathways)),
        "negative_gapfill_deltas": int((~individual["GapfillDeltaValid"]).sum()),
        "missing_blocked_reaction_counts": int(individual["BlockedReactions"].isna().sum()),
        "donors": {},
        "overlap_jaccard": {},
    }
    for donor in ("Donor 1", "Donor 2"):
        q = quality_summary.set_index("Donor").loc[donor]
        metrics["donors"][donor] = {
            "mags": int(q["MAGs"]),
            "paper_high": int(q["PaperHigh"]),
            "paper_medium": int(q["PaperMedium"]),
            "normalized_mags": int(q["NormalizedInput"]),
            "reference_unique_normalized": int(normalized_pivot.loc[donor, ("Unique metabolites", "Reference ≥0.5%")]),
            "mag_unique_normalized": int(normalized_pivot.loc[donor, ("Unique metabolites", "MAG normalized")]),
            "reference_total_normalized": int(normalized_pivot.loc[donor, ("Total predicted metabolites", "Reference ≥0.5%")]),
            "mag_total_normalized": int(normalized_pivot.loc[donor, ("Total predicted metabolites", "MAG normalized")]),
        }
    for row in overlaps.itertuples(index=False):
        metrics["overlap_jaccard"][f"{row.Donor} · {row.Evidence}"] = float(row.Jaccard)
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    versions = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "openpyxl": importlib.metadata.version("openpyxl"),
    }
    (output / "software-versions.json").write_text(
        json.dumps(versions, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
