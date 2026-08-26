#!/usr/bin/env python3
"""Prepare threshold, provenance, and MGX/MTX/MPX evidence for Article 65."""

from __future__ import annotations

import argparse
import gzip
import html
import json
import platform
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


SEED = 65_001
PLOT_SEED = 20_260_765
BOOTSTRAPS = 2_000
PREVALENCE_GATE = 0.10
DIAGNOSIS_ORDER = ("Control", "CD", "UC")
THRESHOLDS = {
    "1 peptide · 1% FDR": "mpx-1pep-1pct.tsv.gz",
    "1 peptide · 5% FDR": "mpx-1pep-5pct.tsv.gz",
    "2 peptides · 1% FDR": "mpx-2pep-1pct.tsv.gz",
    "2 peptides · 5% FDR": "mpx-2pep-5pct.tsv.gz",
}
PAIR_ORDER = ("DNA–RNA", "DNA–Protein", "RNA–Protein")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def feature_parts(value: str) -> tuple[str, str]:
    identifier, separator, description = value.partition(":")
    description = html.unescape(description.strip()) if separator else identifier
    return identifier.strip(), description


def stream_unstratified(source: Path, target: Path) -> int:
    """Copy only HUMAnN unstratified rows without loading the large table."""
    count = 0
    with gzip.open(source, "rt", encoding="utf-8") as input_handle, gzip.open(
        target, "wt", encoding="utf-8", compresslevel=6
    ) as output_handle:
        header = input_handle.readline()
        if not header:
            raise RuntimeError(f"Empty HUMAnN table: {source}")
        output_handle.write(header)
        for line in input_handle:
            feature = line.split("\t", 1)[0]
            if "|" not in feature:
                output_handle.write(line)
                count += 1
    return count


def load_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    feature = frame.columns[0]
    if frame[feature].duplicated().any():
        raise RuntimeError(f"Duplicated feature identifiers in {path}")
    return frame.set_index(feature).astype(float)


def write_matrix(frame: pd.DataFrame, path: Path, index_name: str) -> None:
    output = frame.copy()
    output.index.name = index_name
    output.reset_index().to_csv(path, sep="\t", index=False, compression="gzip")


def bh(values: pd.Series) -> np.ndarray:
    output = np.full(len(values), np.nan, dtype=float)
    keep = values.notna().to_numpy()
    if keep.any():
        output[keep] = multipletests(values.to_numpy()[keep], method="fdr_bh")[1]
    return output


def bootstrap_median(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan
    indices = rng.integers(0, values.size, size=(BOOTSTRAPS, values.size))
    estimates = np.median(values[indices], axis=1)
    return tuple(np.quantile(estimates, [0.025, 0.975]).tolist())


def namespace(identifier: str) -> str:
    if identifier.startswith("HR:"):
        return "Explicit host namespace"
    if identifier.startswith(("CNTM:", "Contaminant_")):
        return "Contaminant namespace"
    if ":" in identifier:
        return "Taxon-prefixed reference"
    return "Generic accession"


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((cache / "download-manifest.json").read_text(encoding="utf-8"))
    if manifest.get("article") != 65:
        raise RuntimeError("Article 65 source manifest is missing or incompatible")

    # Four published protein-filter tables: quantify sensitivity before selecting one.
    threshold_records: list[dict[str, object]] = []
    richness: dict[str, pd.Series] = {}
    identifier_sets: dict[str, set[str]] = {}
    primary: pd.DataFrame | None = None
    for label, filename in THRESHOLDS.items():
        table = pd.read_csv(cache / filename, sep="\t", low_memory=False)
        identifier_column = table.columns[0]
        identifiers = table[identifier_column].astype(str)
        values = table.iloc[:, 1:].astype(float)
        detected_per_protein = values.gt(0).any(axis=1)
        detected_per_sample = values.gt(0).sum(axis=0)
        threshold_records.append(
            {
                "Threshold": label,
                "ProteinTableRows": len(table),
                "Profiles": values.shape[1],
                "DetectedProteinIDs": int(detected_per_protein.sum()),
                "MedianDetectedPerSample": float(detected_per_sample.median()),
                "IQRLowDetectedPerSample": float(detected_per_sample.quantile(0.25)),
                "IQRHighDetectedPerSample": float(detected_per_sample.quantile(0.75)),
                "TotalReportedCounts": float(values.to_numpy().sum()),
            }
        )
        richness[label] = detected_per_sample.rename_axis("SampleID")
        identifier_sets[label] = set(identifiers)
        if label == "2 peptides · 1% FDR":
            primary = table
        else:
            del table, values
    if primary is None:
        raise RuntimeError("Primary 2-peptide/1%-FDR table was not loaded")
    threshold_audit = pd.DataFrame(threshold_records)

    richness_table = pd.concat(richness, axis=1).reset_index()
    richness_correlations: list[dict[str, object]] = []
    for first, second in combinations(THRESHOLDS, 2):
        paired = richness_table[[first, second]].dropna()
        result = spearmanr(paired[first], paired[second])
        richness_correlations.append(
            {
                "ThresholdA": first,
                "ThresholdB": second,
                "SharedProfiles": len(paired),
                "SpearmanRho": float(result.statistic),
                "PValue": float(result.pvalue),
            }
        )
    richness_correlation = pd.DataFrame(richness_correlations)

    overlap_records: list[dict[str, object]] = []
    for first, second in combinations(THRESHOLDS, 2):
        a, b = identifier_sets[first], identifier_sets[second]
        overlap_records.append(
            {
                "ThresholdA": first,
                "ThresholdB": second,
                "Intersection": len(a & b),
                "Union": len(a | b),
                "OnlyA": len(a - b),
                "OnlyB": len(b - a),
                "Jaccard": len(a & b) / len(a | b),
            }
        )
    identifier_overlap = pd.DataFrame(overlap_records)

    primary_ids = primary.iloc[:, 0].astype(str)
    primary_values = primary.iloc[:, 1:].astype(float)
    namespace_audit = pd.DataFrame(
        {
            "Namespace": primary_ids.map(namespace),
            "ReportedCounts": primary_values.sum(axis=1).to_numpy(float),
            "Detected": primary_values.gt(0).any(axis=1).to_numpy(),
            "DetectedEntries": primary_values.gt(0).sum(axis=1).to_numpy(int),
        }
    )
    namespace_audit = (
        namespace_audit.groupby("Namespace", as_index=False, sort=False)
        .agg(
            ProteinIDs=("Detected", "size"),
            DetectedProteinIDs=("Detected", "sum"),
            ReportedCounts=("ReportedCounts", "sum"),
            DetectedEntries=("DetectedEntries", "sum"),
        )
    )
    namespace_audit["ReportedCountShare"] = (
        namespace_audit["ReportedCounts"] / namespace_audit["ReportedCounts"].sum()
    )
    del primary, primary_values

    # Stream the very large HUMAnN EC tables down to their unstratified rows.
    mgx_small = output / "mgx-ec-unstratified.tsv.gz"
    mtx_small = output / "mtx-ec-unstratified.tsv.gz"
    mgx_rows = stream_unstratified(cache / "mgx-ecs-rela.tsv.gz", mgx_small)
    mtx_rows = stream_unstratified(cache / "mtx-ecs-rela.tsv.gz", mtx_small)
    dna_source = load_matrix(mgx_small)
    rna_source = load_matrix(mtx_small)
    protein_source = load_matrix(cache / "mpx-ecs.tsv.gz")
    if (mgx_rows, mtx_rows) != (len(dna_source), len(rna_source)):
        raise RuntimeError("Streamed EC row count mismatch")

    triple_raw = sorted(set(dna_source.columns) & set(rna_source.columns) & set(protein_source.columns))
    technical_replicates = sorted(sample for sample in triple_raw if sample.endswith("_TR"))
    triple_raw = [sample for sample in triple_raw if not sample.endswith("_TR")]
    shared_ec = sorted(set(dna_source.index) & set(rna_source.index) & set(protein_source.index))
    dna = dna_source.loc[shared_ec, triple_raw].T
    rna = rna_source.loc[shared_ec, triple_raw].T
    protein = protein_source.loc[shared_ec, triple_raw].T
    complete = dna.sum(axis=1).gt(0) & rna.sum(axis=1).gt(0) & protein.sum(axis=1).gt(0)
    excluded_zero = sorted(dna.index[~complete].tolist())
    dna, rna, protein = dna.loc[complete].copy(), rna.loc[complete].copy(), protein.loc[complete].copy()
    dna_relative = dna.div(dna.sum(axis=1), axis=0)
    rna_relative = rna.div(rna.sum(axis=1), axis=0)
    protein_relative = protein.div(protein.sum(axis=1), axis=0)

    metadata_source = pd.read_csv(cache / "hmp2-metadata.csv", low_memory=False)
    meta = metadata_source.loc[
        metadata_source["data_type"].eq("proteomics")
        & metadata_source["External ID"].isin(dna_relative.index)
    ].copy()
    if meta["External ID"].duplicated().any() or set(meta["External ID"]) != set(dna_relative.index):
        raise RuntimeError("Proteomics metadata does not map one-to-one to triple-assay samples")
    meta = meta.set_index("External ID").loc[dna_relative.index]
    diagnosis = meta["diagnosis"].replace({"nonIBD": "Control"})
    if set(diagnosis) != set(DIAGNOSIS_ORDER):
        raise RuntimeError(f"Unexpected diagnosis labels: {sorted(set(diagnosis))}")
    sample_metadata = pd.DataFrame(
        {
            "SampleID": dna_relative.index,
            "SubjectID": meta["Participant ID"].astype(str).to_numpy(),
            "CollectionID": meta["site_sub_coll"].astype(str).to_numpy(),
            "Week": pd.to_numeric(meta["week_num"], errors="coerce").to_numpy(),
            "Diagnosis": diagnosis.to_numpy(),
            "DNAECSum": dna.sum(axis=1).to_numpy(),
            "RNAECSum": rna.sum(axis=1).to_numpy(),
            "ProteinECSum": protein.sum(axis=1).to_numpy(),
        }
    )
    if sample_metadata.groupby("SubjectID")["Diagnosis"].nunique().max() != 1:
        raise RuntimeError("Diagnosis changes within a participant")

    metabolomics_ids = set(
        metadata_source.loc[metadata_source["data_type"].eq("metabolomics"), "External ID"].astype(str)
    )
    raw_four_layer = sorted(set(triple_raw) & metabolomics_ids)
    complete_four_layer = sorted(set(dna_relative.index) & metabolomics_ids)
    sample_metadata["MetabolomicsAvailable"] = sample_metadata["SampleID"].isin(metabolomics_ids)

    prevalence = pd.DataFrame(
        {
            "Feature": shared_ec,
            "DNAPrevalence": dna_relative.gt(0).mean(axis=0).reindex(shared_ec).to_numpy(),
            "RNAPrevalence": rna_relative.gt(0).mean(axis=0).reindex(shared_ec).to_numpy(),
            "ProteinPrevalence": protein_relative.gt(0).mean(axis=0).reindex(shared_ec).to_numpy(),
        }
    )
    parts = prevalence["Feature"].map(feature_parts)
    prevalence["EC"] = parts.map(lambda value: value[0])
    prevalence["Description"] = parts.map(lambda value: value[1])
    prevalence["MinimumPrevalence"] = prevalence[
        ["DNAPrevalence", "RNAPrevalence", "ProteinPrevalence"]
    ].min(axis=1)
    prevalence["Selected"] = prevalence["MinimumPrevalence"].ge(PREVALENCE_GATE)
    selected = prevalence.loc[prevalence["Selected"], "Feature"].tolist()
    dna_selected = dna_relative.loc[:, selected]
    rna_selected = rna_relative.loc[:, selected]
    protein_selected = protein_relative.loc[:, selected]

    layers = {"DNA": dna_selected, "RNA": rna_selected, "Protein": protein_selected}
    pairs = (("DNA", "RNA"), ("DNA", "Protein"), ("RNA", "Protein"))
    sample_records: list[dict[str, object]] = []
    sample_lookup = sample_metadata.set_index("SampleID")
    for sample in dna_selected.index:
        for first, second in pairs:
            result = spearmanr(layers[first].loc[sample], layers[second].loc[sample])
            sample_records.append(
                {
                    "SampleID": sample,
                    "SubjectID": sample_lookup.loc[sample, "SubjectID"],
                    "Diagnosis": sample_lookup.loc[sample, "Diagnosis"],
                    "LayerPair": f"{first}–{second}",
                    "SpearmanRho": float(result.statistic),
                }
            )
    sample_concordance = pd.DataFrame(sample_records)
    subject_concordance = (
        sample_concordance.groupby(["SubjectID", "Diagnosis", "LayerPair"], as_index=False)["SpearmanRho"]
        .median()
        .rename(columns={"SpearmanRho": "MedianSpearmanRho"})
    )
    rng = np.random.default_rng(SEED)
    concordance_summary_records: list[dict[str, object]] = []
    for pair in PAIR_ORDER:
        values = subject_concordance.loc[
            subject_concordance["LayerPair"].eq(pair), "MedianSpearmanRho"
        ].to_numpy(float)
        ci_low, ci_high = bootstrap_median(values, rng)
        concordance_summary_records.append(
            {
                "LayerPair": pair,
                "Subjects": values.size,
                "MedianSpearmanRho": float(np.median(values)),
                "CILow": ci_low,
                "CIHigh": ci_high,
            }
        )
    concordance_summary = pd.DataFrame(concordance_summary_records)

    # Collapse repeated samples before EC-wise correlations; BH spans all 3 x 263 tests.
    subjects = sample_metadata.set_index("SampleID")["SubjectID"]
    subject_layers: dict[str, pd.DataFrame] = {}
    for layer, frame in layers.items():
        work = frame.copy()
        work.insert(0, "SubjectID", subjects.loc[work.index])
        subject_layers[layer] = work.groupby("SubjectID").median(numeric_only=True)
    correlation_records: list[dict[str, object]] = []
    for first, second in pairs:
        for feature in selected:
            result = spearmanr(subject_layers[first][feature], subject_layers[second][feature])
            identifier, description = feature_parts(feature)
            correlation_records.append(
                {
                    "Feature": feature,
                    "EC": identifier,
                    "Description": description,
                    "LayerPair": f"{first}–{second}",
                    "Subjects": subject_layers[first].shape[0],
                    "SpearmanRho": float(result.statistic),
                    "PValue": float(result.pvalue),
                }
            )
    ec_correlations = pd.DataFrame(correlation_records)
    ec_correlations["QValue"] = bh(ec_correlations["PValue"])
    ec_correlations["BH05"] = ec_correlations["QValue"].lt(0.05)

    mpx_profiles = protein_source.shape[1]
    mpx_mgx = sorted(set(protein_source.columns) & set(dna_source.columns))
    attrition = pd.DataFrame(
        {
            "Stage": [
                "MPX profiles",
                "Exact MPX + MGX",
                "Exact MPX + MGX + MTX",
                "Three-layer EC-complete",
                "Triple profiles with MBX product",
                "Four-layer & EC-complete",
            ],
            "Count": [
                mpx_profiles,
                len(mpx_mgx),
                len(triple_raw),
                len(dna_relative),
                len(raw_four_layer),
                len(complete_four_layer),
            ],
            "Unit": ["profiles", "samples", "samples", "samples", "samples", "samples"],
        }
    )

    subject_diagnosis = sample_metadata.drop_duplicates("SubjectID").set_index("SubjectID")["Diagnosis"]
    threshold_audit.to_csv(output / "threshold-audit.tsv", sep="\t", index=False)
    richness_table.to_csv(output / "sample-protein-richness.tsv", sep="\t", index=False)
    richness_correlation.to_csv(output / "richness-correlations.tsv", sep="\t", index=False)
    identifier_overlap.to_csv(output / "protein-id-overlap.tsv", sep="\t", index=False)
    namespace_audit.to_csv(output / "protein-namespace-audit.tsv", sep="\t", index=False)
    sample_metadata.to_csv(output / "sample-metadata.tsv", sep="\t", index=False)
    prevalence.to_csv(output / "ec-feature-audit.tsv", sep="\t", index=False)
    sample_concordance.to_csv(output / "sample-concordance.tsv", sep="\t", index=False)
    subject_concordance.to_csv(output / "subject-concordance.tsv", sep="\t", index=False)
    concordance_summary.to_csv(output / "concordance-summary.tsv", sep="\t", index=False)
    ec_correlations.to_csv(output / "ec-correlations.tsv", sep="\t", index=False)
    attrition.to_csv(output / "sample-attrition.tsv", sep="\t", index=False)
    write_matrix(dna_selected, output / "dna-ec-relative.tsv.gz", "SampleID")
    write_matrix(rna_selected, output / "rna-ec-relative.tsv.gz", "SampleID")
    write_matrix(protein_selected, output / "protein-ec-relative.tsv.gz", "SampleID")
    for layer, frame in subject_layers.items():
        write_matrix(frame, output / f"subject-{layer.lower()}-ec-relative.tsv.gz", "SubjectID")

    pair_metrics = {}
    for pair in PAIR_ORDER:
        section = ec_correlations.loc[ec_correlations["LayerPair"].eq(pair)]
        pair_metrics[pair] = {
            "median_subject_sample_spearman": float(
                concordance_summary.loc[concordance_summary["LayerPair"].eq(pair), "MedianSpearmanRho"].iloc[0]
            ),
            "median_ec_cross_subject_spearman": float(section["SpearmanRho"].median()),
            "ec_bh05": int(section["BH05"].sum()),
        }
    metrics = {
        "article": 65,
        "seed": SEED,
        "plot_seed": PLOT_SEED,
        "bootstrap_replicates": BOOTSTRAPS,
        "prevalence_gate": PREVALENCE_GATE,
        "source_mgx_profiles": int(dna_source.shape[1]),
        "source_mtx_profiles": int(rna_source.shape[1]),
        "source_mpx_profiles": int(protein_source.shape[1]),
        "source_mgx_unstratified_ecs": int(dna_source.shape[0]),
        "source_mtx_unstratified_ecs": int(rna_source.shape[0]),
        "source_mpx_ec_rows": int(protein_source.shape[0]),
        "shared_ecs": len(shared_ec),
        "selected_ecs": len(selected),
        "raw_triple_samples": len(triple_raw),
        "technical_replicates_excluded": len(technical_replicates),
        "zero_layer_samples_excluded": len(excluded_zero),
        "analysis_samples": len(dna_relative),
        "independent_subjects": int(sample_metadata["SubjectID"].nunique()),
        "diagnosis_subjects": subject_diagnosis.value_counts().sort_index().to_dict(),
        "raw_four_layer_samples": len(raw_four_layer),
        "complete_four_layer_samples": len(complete_four_layer),
        "ec_correlation_tests": len(ec_correlations),
        "ec_correlation_bh05_total": int(ec_correlations["BH05"].sum()),
        "pair_metrics": pair_metrics,
    }
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    software = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "seed": SEED,
    }
    (output / "software-versions.json").write_text(
        json.dumps(software, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "source-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "exclusion-ledger.json").write_text(
        json.dumps(
            {
                "technical_replicates": technical_replicates,
                "zero_layer_samples": excluded_zero,
                "four_layer_exact_sample_ids": complete_four_layer,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
