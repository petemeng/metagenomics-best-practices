#!/usr/bin/env python3
"""Create compact, auditable summaries from the Article 61 model outputs."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


POSITIVE_GROWTH = 1e-6
PRIMARY_TRADEOFF = 0.5
DISPLAY_NAMES = {
    "ac": "Acetate",
    "but": "Butyrate",
    "ppa": "Propionate",
    "succ": "Succinate",
    "lac_L": "L-lactate",
    "lac_D": "D-lactate",
    "for": "Formate",
    "etoh": "Ethanol",
    "co2": "Carbon dioxide",
    "h2": "Hydrogen",
    "nh4": "Ammonium",
    "glc__D": "D-glucose",
    "fru": "D-fructose",
    "gal": "D-galactose",
    "o2": "Oxygen",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--summary-dir", type=Path, required=True)
    return parser.parse_args()


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def normalize_compound(value: object) -> str:
    text = str(value)
    text = re.sub(r"^R_", "", text)
    text = re.sub(r"^EX_", "", text)
    text = re.sub(r"^M_", "", text)
    text = re.sub(r"\[(?:e|m)\]$", "", text)
    text = re.sub(r"_(?:e|m)(?:_pool)?$", "", text)
    text = re.sub(r"\((?:e|m)\)$", "", text)
    return text


def display_name(compound: str) -> str:
    return DISPLAY_NAMES.get(compound, compound.replace("__", "-").replace("_", " "))


def growth_summary(frame: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    data = frame.copy()
    data["PositiveGrowth"] = data.growth_rate.gt(POSITIVE_GROWTH)
    data["WeightedGrowth"] = data.growth_rate * data.abundance
    data["PositiveAbundance"] = data.abundance.where(data.PositiveGrowth, 0.0)
    summary = (
        data.groupby(group_columns, dropna=False, sort=True)
        .agg(
            ModeledTaxa=("taxon", "nunique"),
            PositiveTaxa=("PositiveGrowth", "sum"),
            FractionGrowing=("PositiveGrowth", "mean"),
            ModeledAbundance=("abundance", "sum"),
            PositiveAbundance=("PositiveAbundance", "sum"),
            CommunityGrowth=("WeightedGrowth", "sum"),
            MedianTaxonGrowth=("growth_rate", "median"),
        )
        .reset_index()
    )
    summary["AbundanceWeightedFractionGrowing"] = (
        summary.PositiveAbundance / summary.ModeledAbundance
    )
    return summary


def main() -> None:
    args = parse_args()
    work = args.work_dir.resolve()
    out = args.summary_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    selected = read_tsv(work / "selected-samples.tsv")
    sample_map = selected.set_index("SampleID").SubjectID.to_dict()
    coverage = read_tsv(work / "model-coverage.tsv")
    write_tsv(out / "sample-selection-audit.tsv", read_tsv(work / "sample-selection-audit.tsv"))
    write_tsv(out / "selected-samples.tsv", selected)
    write_tsv(out / "model-coverage.tsv", coverage)
    write_tsv(out / "taxon-match-audit.tsv", read_tsv(work / "taxon-match-audit.tsv"))

    tradeoff = read_tsv(work / "micom-tradeoff-growth.tsv")
    tradeoff = tradeoff.loc[tradeoff.tradeoff.notna()].copy()
    tradeoff["SubjectID"] = tradeoff.sample_id.map(sample_map)
    tradeoff_summary = growth_summary(tradeoff, ["sample_id", "SubjectID", "tradeoff"])
    write_tsv(out / "tradeoff-summary.tsv", tradeoff_summary)

    primary_growth = read_tsv(work / "micom-primary-growth.tsv")
    primary_growth["SubjectID"] = primary_growth.sample_id.map(sample_map)
    primary_growth["PositiveGrowth"] = primary_growth.growth_rate.gt(POSITIVE_GROWTH)
    primary_growth["WeightedGrowth"] = primary_growth.growth_rate * primary_growth.abundance
    write_tsv(out / "primary-growth.tsv", primary_growth)
    primary_summary = growth_summary(primary_growth, ["sample_id", "SubjectID"])
    write_tsv(out / "primary-growth-summary.tsv", primary_summary)

    medium = read_tsv(work / "micom-medium-sensitivity.tsv")
    medium["SubjectID"] = medium.sample_id.map(sample_map)
    medium_summary = growth_summary(medium, ["sample_id", "SubjectID", "MediumScale"])
    write_tsv(out / "medium-sensitivity-summary.tsv", medium_summary)

    equal = read_tsv(work / "micom-equal-abundance-growth.tsv")
    equal["SubjectID"] = equal.sample_id.map(sample_map)
    equal_summary = growth_summary(equal, ["sample_id", "SubjectID"])
    observed = tradeoff.loc[np.isclose(tradeoff.tradeoff, PRIMARY_TRADEOFF)].copy()
    observed_summary = growth_summary(observed, ["sample_id", "SubjectID"])
    comparison = observed_summary.merge(
        equal_summary,
        on=["sample_id", "SubjectID"],
        suffixes=("Observed", "Equal"),
        validate="one_to_one",
    )
    comparison["CommunityGrowthRatioEqualToObserved"] = (
        comparison.CommunityGrowthEqual / comparison.CommunityGrowthObserved.replace(0, np.nan)
    )
    comparison["GrowingFractionDifference"] = (
        comparison.FractionGrowingEqual - comparison.FractionGrowingObserved
    )
    write_tsv(out / "abundance-sensitivity-summary.tsv", comparison)

    exchanges = read_tsv(work / "micom-primary-exchanges.tsv")
    exchanges["SubjectID"] = exchanges.sample_id.map(sample_map)
    exchanges["Compound"] = exchanges.reaction.map(normalize_compound)
    exchanges["CompoundName"] = exchanges.Compound.map(display_name)
    exchanges["CommunityScaledFlux"] = np.where(
        exchanges.taxon.eq("medium"),
        exchanges.flux,
        exchanges.flux * exchanges.abundance,
    )
    write_tsv(out / "primary-exchanges.tsv", exchanges)

    net = exchanges.loc[exchanges.taxon.eq("medium")].copy()
    net["NetDirection"] = np.where(net.CommunityScaledFlux.gt(0), "Export", "Import")
    net["AbsoluteFlux"] = net.CommunityScaledFlux.abs()
    net = net[
        [
            "sample_id", "SubjectID", "Compound", "CompoundName", "reaction",
            "CommunityScaledFlux", "AbsoluteFlux", "NetDirection",
        ]
    ].sort_values(["sample_id", "AbsoluteFlux"], ascending=[True, False])
    write_tsv(out / "net-community-flux.tsv", net)

    taxa_flux = exchanges.loc[~exchanges.taxon.eq("medium")].copy()
    crossfeed_rows: list[dict[str, object]] = []
    for (sample_id, compound), group in taxa_flux.groupby(["sample_id", "Compound"], sort=True):
        producers = group.loc[group.CommunityScaledFlux.gt(POSITIVE_GROWTH)]
        consumers = group.loc[group.CommunityScaledFlux.lt(-POSITIVE_GROWTH)]
        if producers.empty or consumers.empty:
            continue
        export_flux = float(producers.CommunityScaledFlux.sum())
        import_flux = float(-consumers.CommunityScaledFlux.sum())
        crossfeed_rows.append(
            {
                "SampleID": sample_id,
                "SubjectID": sample_map[sample_id],
                "Compound": compound,
                "CompoundName": display_name(compound),
                "Producers": producers.taxon.nunique(),
                "Consumers": consumers.taxon.nunique(),
                "ExportFlux": export_flux,
                "ImportFlux": import_flux,
                "PotentialTurnover": min(export_flux, import_flux),
            }
        )
    crossfeed_columns = [
        "SampleID", "SubjectID", "Compound", "CompoundName", "Producers",
        "Consumers", "ExportFlux", "ImportFlux", "PotentialTurnover",
    ]
    crossfeed = pd.DataFrame(crossfeed_rows, columns=crossfeed_columns)
    if not crossfeed.empty:
        crossfeed.sort_values(
            ["SampleID", "PotentialTurnover"], ascending=[True, False], inplace=True
        )
    write_tsv(out / "micom-crossfeeding-potential.tsv", crossfeed)

    subcommunity = read_tsv(work / "smetana-subcommunity.tsv")
    focal_sample = subcommunity.SampleID.iloc[0]
    focal_taxa = set(subcommunity.ModelID)
    focal = taxa_flux.loc[
        taxa_flux.sample_id.eq(focal_sample) & taxa_flux.taxon.isin(focal_taxa)
    ].copy()
    focal_matrix = focal[
        ["taxon", "Compound", "CompoundName", "CommunityScaledFlux"]
    ].sort_values(["Compound", "taxon"])
    write_tsv(out / "focal-micom-fluxes.tsv", focal_matrix)

    edge_rows: list[dict[str, object]] = []
    for compound, group in focal.groupby("Compound", sort=True):
        donors = group.loc[group.CommunityScaledFlux.gt(POSITIVE_GROWTH)]
        receivers = group.loc[group.CommunityScaledFlux.lt(-POSITIVE_GROWTH)]
        for donor in donors.itertuples(index=False):
            for receiver in receivers.itertuples(index=False):
                edge_rows.append(
                    {
                        "Donor": donor.taxon,
                        "Receiver": receiver.taxon,
                        "Compound": compound,
                        "CompoundName": display_name(compound),
                        "DonorExport": donor.CommunityScaledFlux,
                        "ReceiverImport": -receiver.CommunityScaledFlux,
                        "PotentialFlux": min(donor.CommunityScaledFlux, -receiver.CommunityScaledFlux),
                    }
                )
    edge_columns = [
        "Donor", "Receiver", "Compound", "CompoundName", "DonorExport",
        "ReceiverImport", "PotentialFlux",
    ]
    micom_edges = pd.DataFrame(edge_rows, columns=edge_columns)
    if not micom_edges.empty:
        micom_edges.sort_values("PotentialFlux", ascending=False, inplace=True)
    write_tsv(out / "focal-micom-potential-edges.tsv", micom_edges)

    smetana_global = pd.read_csv(
        work / "smetana-global.tsv", sep="\t", keep_default_na=False
    )
    smetana_detailed = read_tsv(work / "smetana-detailed.tsv")
    compatibility_audit = read_tsv(work / "smetana-compatibility-audit.tsv")
    compatibility_summary = json.loads(
        (work / "smetana-compatibility-summary.json").read_text(encoding="utf-8")
    )
    for column in ("scs", "mus", "mps", "smetana"):
        smetana_detailed[column] = pd.to_numeric(smetana_detailed[column], errors="coerce")
    smetana_detailed["Compound"] = smetana_detailed.compound.map(normalize_compound)
    smetana_detailed["CompoundName"] = smetana_detailed.Compound.map(display_name)
    smetana_detailed["PositiveSMETANA"] = smetana_detailed.smetana.fillna(0).gt(0)
    write_tsv(out / "smetana-global.tsv", smetana_global)
    write_tsv(out / "smetana-detailed.tsv", smetana_detailed)
    write_tsv(out / "smetana-compatibility-audit.tsv", compatibility_audit)
    component_summary = pd.DataFrame(
        [
            {
                "Component": "SCS",
                "PositiveRows": compatibility_summary["positive_scs_rows"],
            },
            {
                "Component": "MUS",
                "PositiveRows": compatibility_summary["positive_mus_rows"],
            },
            {
                "Component": "MPS",
                "PositiveRows": compatibility_summary["positive_mps_rows"],
            },
            {
                "Component": "Composite",
                "PositiveRows": compatibility_summary["positive_smetana_rows"],
            },
        ]
    )
    write_tsv(out / "smetana-component-summary.tsv", component_summary)

    positive_smetana = smetana_detailed.loc[smetana_detailed.PositiveSMETANA].copy()
    if positive_smetana.empty:
        pair_scores = pd.DataFrame(
            columns=["receiver", "donor", "TotalSMETANA", "MaxSMETANA", "PositiveCompounds"]
        )
    else:
        pair_scores = (
            positive_smetana.groupby(["receiver", "donor"], sort=True)
            .agg(
                TotalSMETANA=("smetana", "sum"),
                MaxSMETANA=("smetana", "max"),
                PositiveCompounds=("Compound", "nunique"),
            )
            .reset_index()
        )
    write_tsv(out / "smetana-pair-summary.tsv", pair_scores)

    micom_compare = micom_edges.rename(
        columns={"Donor": "donor", "Receiver": "receiver", "PotentialFlux": "MICOMPotentialFlux"}
    )[["receiver", "donor", "Compound", "MICOMPotentialFlux"]] if not micom_edges.empty else pd.DataFrame(
        columns=["receiver", "donor", "Compound", "MICOMPotentialFlux"]
    )
    comparison_limited = (
        compatibility_summary["interpretation"]
        == "software_model_interface_limitation_not_biological_absence"
    )
    if comparison_limited:
        concordance = micom_compare.copy()
        for column in ("smetana", "scs", "mus", "mps"):
            concordance[column] = np.nan
        concordance["MICOMPositive"] = True
        concordance["SMETANAPositive"] = False
        concordance["EvidenceClass"] = "MICOM candidate; SMETANA unavailable"
    else:
        smetana_compare = positive_smetana[
            ["receiver", "donor", "Compound", "smetana", "scs", "mus", "mps"]
        ]
        concordance = micom_compare.merge(
            smetana_compare,
            on=["receiver", "donor", "Compound"],
            how="outer",
        )
        concordance["MICOMPositive"] = concordance.MICOMPotentialFlux.fillna(0).gt(0)
        concordance["SMETANAPositive"] = concordance.smetana.fillna(0).gt(0)
        concordance["EvidenceClass"] = np.select(
            [
                concordance.MICOMPositive & concordance.SMETANAPositive,
                concordance.MICOMPositive,
                concordance.SMETANAPositive,
            ],
            ["Both", "MICOM only", "SMETANA only"],
            default="Neither",
        )
    concordance.sort_values(
        ["EvidenceClass", "MICOMPotentialFlux", "smetana"],
        ascending=[True, False, False],
        inplace=True,
    )
    write_tsv(out / "cross-method-concordance.tsv", concordance)

    ledger = read_tsv(work / "run-ledger.tsv")
    write_tsv(out / "run-ledger.tsv", ledger)
    metrics = {
        "article": 61,
        "selected_samples": len(selected),
        "independent_subjects": selected.SubjectID.nunique(),
        "minimum_reads": int(selected.Reads.min()),
        "minimum_modeled_abundance": float(coverage.ModeledAbundance.min()),
        "primary_tradeoff": PRIMARY_TRADEOFF,
        "primary_positive_taxa": int(primary_growth.PositiveGrowth.sum()),
        "primary_total_taxa": len(primary_growth),
        "primary_fraction_growing_median": float(primary_summary.FractionGrowing.median()),
        "net_flux_rows": len(net),
        "crossfeeding_compound_rows": len(crossfeed),
        "smetana_subcommunity_size": len(subcommunity),
        "smetana_detailed_rows": len(smetana_detailed),
        "smetana_positive_rows": int(smetana_detailed.PositiveSMETANA.sum()),
        "smetana_positive_scs_rows": compatibility_summary["positive_scs_rows"],
        "smetana_positive_mus_rows": compatibility_summary["positive_mus_rows"],
        "smetana_positive_mps_rows": compatibility_summary["positive_mps_rows"],
        "smetana_global_estimable": False,
        "smetana_cross_method_comparison_estimable": not comparison_limited,
        "smetana_medium_matches": compatibility_summary["matched_pooled_exchanges"],
        "smetana_standalone_positive_models": compatibility_summary["standalone_positive_models"],
        "smetana_interacting_positive_models": compatibility_summary["interacting_positive_models"],
        "smetana_legacy_noninteracting_positive_models": compatibility_summary["legacy_noninteracting_positive_models"],
        "cross_method_both_rows": int(concordance.EvidenceClass.eq("Both").sum()),
        "completed_steps": int(ledger.Status.eq("passed").sum()),
        "passed_with_limitation_steps": int(
            ledger.Status.eq("passed_with_limitation").sum()
        ),
        "not_estimable_steps": int(ledger.Status.eq("not_estimable").sum()),
        "elapsed_seconds": float(ledger.ElapsedSeconds.sum()),
        "maximum_recorded_rss_kb": int(ledger.CumulativeMaxRSSKB.max()),
    }
    (out / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
