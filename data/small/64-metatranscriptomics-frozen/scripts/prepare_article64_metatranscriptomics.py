#!/usr/bin/env python3
"""Prepare the paired MGX/MTX analysis and frozen evidence for Article 64."""

from __future__ import annotations

import argparse
import html
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import kruskal, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests


SEED = 64_001
PLOT_SEED = 20_260_764
PSEUDOCOUNT = 1e-6
CODETECTION_GATE = 0.20
BOOTSTRAPS = 2_000
DIAGNOSIS_ORDER = ("Control", "CD", "UC")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def load_unstratified(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    feature = frame.columns[0]
    frame[feature] = frame[feature].astype(str).str.strip()
    frame = frame.loc[~frame[feature].str.contains("|", regex=False)].copy()
    if frame[feature].duplicated().any():
        raise RuntimeError(f"Duplicated unstratified pathways in {path}")
    return frame.set_index(feature)


def feature_parts(value: str) -> tuple[str, str]:
    identifier, separator, description = value.partition(":")
    description = html.unescape(description.strip()) if separator else identifier
    return identifier.strip(), description


def bootstrap_median(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan
    indices = rng.integers(0, values.size, size=(BOOTSTRAPS, values.size))
    estimates = np.median(values[indices], axis=1)
    return tuple(np.quantile(estimates, [0.025, 0.975]).tolist())


def bh(values: pd.Series) -> np.ndarray:
    output = np.full(len(values), np.nan, dtype=float)
    keep = values.notna().to_numpy()
    if keep.any():
        output[keep] = multipletests(values.to_numpy()[keep], method="fdr_bh")[1]
    return output


def write_matrix(frame: pd.DataFrame, path: Path, index_name: str) -> None:
    output = frame.copy()
    output.index.name = index_name
    output.reset_index().to_csv(path, sep="\t", index=False, compression="gzip")


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((cache / "download-manifest.json").read_text(encoding="utf-8"))
    if manifest.get("article") != 64:
        raise RuntimeError("Article 64 source manifest is missing or incompatible")

    dna_source = load_unstratified(cache / "mgx-pathabundance-rela.tsv.gz")
    rna_source = load_unstratified(cache / "mtx-pathabundance-rela.tsv.gz")
    metadata_source = pd.read_csv(cache / "hmp2-metadata.csv", low_memory=False)

    paired = sorted(set(dna_source.columns).intersection(rna_source.columns))
    technical_replicates = sorted(sample for sample in paired if sample.endswith("_TR"))
    paired = [sample for sample in paired if not sample.endswith("_TR")]
    shared_features = sorted(set(dna_source.index).intersection(rna_source.index))

    dna = dna_source.loc[shared_features, paired].T.astype(float)
    rna = rna_source.loc[shared_features, paired].T.astype(float)
    dna_total = dna.sum(axis=1)
    rna_total = rna.sum(axis=1)
    assay_complete = (dna_total > 0) & (rna_total > 0)
    excluded_zero = sorted(dna.index[~assay_complete].tolist())
    dna = dna.loc[assay_complete].copy()
    rna = rna.loc[assay_complete].copy()

    meta = metadata_source.loc[
        metadata_source["data_type"].eq("metagenomics")
        & metadata_source["External ID"].isin(dna.index)
    ].copy()
    if meta["External ID"].duplicated().any() or set(meta["External ID"]) != set(dna.index):
        raise RuntimeError("MGX metadata does not map one-to-one to paired biospecimens")
    meta = meta.set_index("External ID").loc[dna.index]
    diagnosis = meta["diagnosis"].replace({"nonIBD": "Control"})
    if set(diagnosis) != set(DIAGNOSIS_ORDER):
        raise RuntimeError(f"Unexpected diagnosis labels: {sorted(set(diagnosis))}")

    sample_metadata = pd.DataFrame(
        {
            "SampleID": dna.index,
            "SubjectID": meta["Participant ID"].astype(str).to_numpy(),
            "CollectionID": meta["site_sub_coll"].astype(str).to_numpy(),
            "Week": pd.to_numeric(meta["week_num"], errors="coerce").to_numpy(),
            "Diagnosis": diagnosis.to_numpy(),
            "DNAPathwaySum": dna.sum(axis=1).to_numpy(),
            "RNAPathwaySum": rna.sum(axis=1).to_numpy(),
        }
    )
    if sample_metadata["SubjectID"].nunique() != 104:
        raise RuntimeError("Unexpected number of independent participants")
    subject_diagnosis = sample_metadata.drop_duplicates("SubjectID").set_index("SubjectID")["Diagnosis"]
    if sample_metadata.groupby("SubjectID")["Diagnosis"].nunique().max() != 1:
        raise RuntimeError("Diagnosis changes within a participant")

    detected_dna = dna.gt(0)
    detected_rna = rna.gt(0)
    codetected = detected_dna & detected_rna
    feature_audit = pd.DataFrame(
        {
            "Feature": shared_features,
            "DNAPrevalence": detected_dna.mean(axis=0).reindex(shared_features).to_numpy(),
            "RNAPrevalence": detected_rna.mean(axis=0).reindex(shared_features).to_numpy(),
            "CoDetection": codetected.mean(axis=0).reindex(shared_features).to_numpy(),
            "MeanDNA": dna.mean(axis=0).reindex(shared_features).to_numpy(),
            "MeanRNA": rna.mean(axis=0).reindex(shared_features).to_numpy(),
        }
    )
    parts = feature_audit["Feature"].map(feature_parts)
    feature_audit["PathwayID"] = parts.map(lambda item: item[0])
    feature_audit["Description"] = parts.map(lambda item: item[1])
    feature_audit["Selected"] = feature_audit["CoDetection"] >= CODETECTION_GATE
    selected = feature_audit.loc[feature_audit["Selected"], "Feature"].tolist()

    dna_selected = dna.loc[:, selected]
    rna_selected = rna.loc[:, selected]
    log_dna = np.log2(dna_selected + PSEUDOCOUNT)
    log_rna = np.log2(rna_selected + PSEUDOCOUNT)

    concordance_records: list[dict[str, object]] = []
    for sample in dna_selected.index:
        statistic = spearmanr(log_dna.loc[sample], log_rna.loc[sample]).statistic
        concordance_records.append({"SampleID": sample, "SpearmanRho": statistic})
    concordance = pd.DataFrame(concordance_records).merge(sample_metadata, on="SampleID", validate="one_to_one")
    subject_concordance = (
        concordance.groupby(["SubjectID", "Diagnosis"], as_index=False)["SpearmanRho"]
        .median()
        .rename(columns={"SpearmanRho": "MedianSpearmanRho"})
    )

    activity_values = np.full(rna_selected.shape, np.nan, dtype=float)
    both = codetected.loc[:, selected].to_numpy()
    ratio = rna_selected.to_numpy() / np.where(dna_selected.to_numpy() > 0, dna_selected.to_numpy(), 1.0)
    activity_values[both] = np.log2(ratio[both])
    activity = pd.DataFrame(activity_values, index=dna_selected.index, columns=selected)
    activity_with_subject = activity.copy()
    activity_with_subject.insert(0, "SubjectID", sample_metadata.set_index("SampleID").loc[activity.index, "SubjectID"])
    subject_activity = activity_with_subject.groupby("SubjectID").median(numeric_only=True)

    rng = np.random.default_rng(SEED)
    activity_records: list[dict[str, object]] = []
    for feature in selected:
        values = subject_activity[feature].dropna().to_numpy(float)
        if values.size < 30:
            p_value = np.nan
        else:
            p_value = float(wilcoxon(values, zero_method="wilcox", alternative="two-sided").pvalue)
        ci_low, ci_high = bootstrap_median(values, rng)
        identifier, description = feature_parts(feature)
        activity_records.append(
            {
                "Feature": feature,
                "PathwayID": identifier,
                "Description": description,
                "Subjects": values.size,
                "MedianLog2RNAoverDNA": float(np.median(values)),
                "MeanLog2RNAoverDNA": float(np.mean(values)),
                "CILow": ci_low,
                "CIHigh": ci_high,
                "PValue": p_value,
            }
        )
    activity_results = pd.DataFrame(activity_records)
    activity_results["QValue"] = bh(activity_results["PValue"])
    activity_results["BH05"] = activity_results["QValue"] < 0.05

    diagnosis_records: list[dict[str, object]] = []
    for feature in selected:
        groups = {
            group: subject_activity.loc[subject_diagnosis.eq(group), feature].dropna().to_numpy(float)
            for group in DIAGNOSIS_ORDER
        }
        group_sizes = {group: values.size for group, values in groups.items()}
        medians = {
            group: (float(np.median(values)) if values.size else np.nan)
            for group, values in groups.items()
        }
        if min(group_sizes.values()) >= 15:
            statistic, p_value = kruskal(*(groups[group] for group in DIAGNOSIS_ORDER))
        else:
            statistic, p_value = np.nan, np.nan
        diagnosis_records.append(
            {
                "Feature": feature,
                "PathwayID": feature_parts(feature)[0],
                "ControlN": group_sizes["Control"],
                "CDN": group_sizes["CD"],
                "UCN": group_sizes["UC"],
                "ControlMedian": medians["Control"],
                "CDMedian": medians["CD"],
                "UCMedian": medians["UC"],
                "MedianRange": float(np.nanmax(list(medians.values())) - np.nanmin(list(medians.values()))),
                "KruskalH": statistic,
                "PValue": p_value,
            }
        )
    diagnosis_results = pd.DataFrame(diagnosis_records)
    diagnosis_results["QValue"] = bh(diagnosis_results["PValue"])
    diagnosis_results["BH05"] = diagnosis_results["QValue"] < 0.05

    concordance_summary_records: list[dict[str, object]] = []
    for group in ("All",) + DIAGNOSIS_ORDER:
        values = subject_concordance.loc[
            subject_concordance["Diagnosis"].eq(group) if group != "All" else np.ones(len(subject_concordance), dtype=bool),
            "MedianSpearmanRho",
        ].to_numpy(float)
        ci_low, ci_high = bootstrap_median(values, rng)
        concordance_summary_records.append(
            {
                "Diagnosis": group,
                "Subjects": values.size,
                "MedianSpearmanRho": float(np.median(values)),
                "CILow": ci_low,
                "CIHigh": ci_high,
            }
        )
    concordance_summary = pd.DataFrame(concordance_summary_records)

    sensitivity_records: list[dict[str, object]] = []
    primary_median = subject_activity.median(axis=0)
    for delta in (1e-7, 1e-6, 1e-5):
        alternate = np.log2((rna_selected + delta) / (dna_selected + delta))
        alternate.insert(0, "SubjectID", sample_metadata.set_index("SampleID").loc[alternate.index, "SubjectID"])
        alternate_subject = alternate.groupby("SubjectID").median(numeric_only=True)
        alternate_median = alternate_subject.median(axis=0)
        rank_rho = spearmanr(primary_median, alternate_median).statistic
        sensitivity_records.append(
            {
                "Analysis": "Pseudocount",
                "Setting": f"{delta:.0e}",
                "Features": len(selected),
                "RankSpearman": rank_rho,
                "MedianAbsoluteShift": float(np.median(np.abs(alternate_median - primary_median))),
            }
        )
    for gate in (0.10, 0.20, 0.30, 0.50):
        sensitivity_records.append(
            {
                "Analysis": "Co-detection gate",
                "Setting": f"{gate:.0%}",
                "Features": int((feature_audit["CoDetection"] >= gate).sum()),
                "RankSpearman": 1.0 if gate == CODETECTION_GATE else np.nan,
                "MedianAbsoluteShift": 0.0 if gate == CODETECTION_GATE else np.nan,
            }
        )
    sensitivity = pd.DataFrame(sensitivity_records)

    attrition = pd.DataFrame(
        {
            "Stage": [
                "MGX profiles",
                "MTX profiles",
                "Raw paired columns",
                "Unique biospecimens",
                "Assay-complete pairs",
                "Independent subjects",
            ],
            "Count": [
                dna_source.shape[1],
                rna_source.shape[1],
                len(paired) + len(technical_replicates),
                len(paired),
                len(dna),
                sample_metadata["SubjectID"].nunique(),
            ],
            "Unit": ["profiles", "profiles", "columns", "samples", "samples", "subjects"],
        }
    )

    sample_metadata.to_csv(output / "sample-metadata.tsv", sep="\t", index=False)
    feature_audit.to_csv(output / "feature-audit.tsv", sep="\t", index=False)
    concordance.to_csv(output / "sample-concordance.tsv", sep="\t", index=False)
    subject_concordance.to_csv(output / "subject-concordance.tsv", sep="\t", index=False)
    concordance_summary.to_csv(output / "concordance-summary.tsv", sep="\t", index=False)
    activity_results.to_csv(output / "activity-results.tsv", sep="\t", index=False)
    diagnosis_results.to_csv(output / "diagnosis-results.tsv", sep="\t", index=False)
    sensitivity.to_csv(output / "sensitivity-summary.tsv", sep="\t", index=False)
    attrition.to_csv(output / "sample-attrition.tsv", sep="\t", index=False)
    write_matrix(dna_selected, output / "dna-relative.tsv.gz", "SampleID")
    write_matrix(rna_selected, output / "rna-relative.tsv.gz", "SampleID")
    write_matrix(activity, output / "activity-log2-rna-dna.tsv.gz", "SampleID")
    write_matrix(subject_activity, output / "subject-activity.tsv.gz", "SubjectID")

    metrics = {
        "article": 64,
        "seed": SEED,
        "plot_seed": PLOT_SEED,
        "source_mgx_profiles": int(dna_source.shape[1]),
        "source_mtx_profiles": int(rna_source.shape[1]),
        "source_mgx_unstratified_pathways": int(dna_source.shape[0]),
        "source_mtx_unstratified_pathways": int(rna_source.shape[0]),
        "raw_paired_columns": len(paired) + len(technical_replicates),
        "technical_replicates_excluded": len(technical_replicates),
        "unique_paired_biospecimens": len(paired),
        "zero_layer_samples_excluded": len(excluded_zero),
        "zero_layer_sample_ids": excluded_zero,
        "analysis_samples": len(dna),
        "independent_subjects": int(sample_metadata["SubjectID"].nunique()),
        "diagnosis_subjects": subject_diagnosis.value_counts().sort_index().to_dict(),
        "shared_unstratified_pathways": len(shared_features),
        "selected_pathways": len(selected),
        "codetection_gate": CODETECTION_GATE,
        "pseudocount": PSEUDOCOUNT,
        "bootstrap_replicates": BOOTSTRAPS,
        "median_sample_spearman": float(concordance["SpearmanRho"].median()),
        "median_subject_spearman": float(subject_concordance["MedianSpearmanRho"].median()),
        "relative_activity_bh05": int(activity_results["BH05"].sum()),
        "relative_activity_bh05_positive": int((activity_results["BH05"] & activity_results["MedianLog2RNAoverDNA"].gt(0)).sum()),
        "relative_activity_bh05_negative": int((activity_results["BH05"] & activity_results["MedianLog2RNAoverDNA"].lt(0)).sum()),
        "diagnosis_tests": int(diagnosis_results["PValue"].notna().sum()),
        "diagnosis_bh05": int(diagnosis_results["BH05"].sum()),
        "minimum_diagnosis_q": float(diagnosis_results["QValue"].min()),
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
