#!/usr/bin/env python3
"""Prepare subject-aware longitudinal HMP2 evidence for Article 67."""

from __future__ import annotations

import argparse
import importlib.metadata
import itertools
import json
import platform
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import braycurtis, cdist


SEED = 67_001
PLOT_SEED = 20_260_767
READ_GATE = 1_000_000
MIN_VISITS = 4
DETECTION = 1e-4
PSEUDOCOUNT = 1e-6
SHIFT_THRESHOLD = 0.54
DIAGNOSIS_ORDER = ("Control", "CD", "UC")
LAG_BINS = (-np.inf, 2, 4, 8, 16, 32, 60)
LAG_LABELS = ("0–2", ">2–4", ">4–8", ">8–16", ">16–32", ">32–60")
ABUNDANCE_LABELS = ("0.01–0.1%", "0.1–1%", "≥1%")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    return parser.parse_args()


def bootstrap_ci(
    values: np.ndarray,
    rng: np.random.Generator,
    iterations: int,
    statistic: str = "median",
) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, np.nan
    function = np.median if statistic == "median" else np.mean
    estimate = float(function(values))
    draws = np.empty(iterations, dtype=float)
    for index in range(iterations):
        draws[index] = function(rng.choice(values, size=len(values), replace=True))
    lower, upper = np.quantile(draws, [0.025, 0.975])
    return estimate, float(lower), float(upper)


def safe_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def load_sources(cache: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    manifest = json.loads((cache / "download-manifest.json").read_text(encoding="utf-8"))
    if manifest.get("article") != 67:
        raise RuntimeError("Article identity mismatch in download manifest")
    profile = pd.read_csv(cache / "taxonomic_profiles.tsv.gz", sep="\t")
    metadata = pd.read_csv(cache / "hmp2-metadata.csv", low_memory=False)
    return profile, metadata, manifest


def select_profiles(
    profile: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_ids = pd.Index(profile.columns[1:].astype(str), name="SampleID")
    meta = metadata.loc[metadata["data_type"].eq("metagenomics")].copy()
    meta["External ID"] = meta["External ID"].astype(str)
    meta = meta.loc[meta["External ID"].isin(sample_ids)].copy()
    if len(meta) != len(sample_ids) or meta["External ID"].duplicated().any():
        raise RuntimeError("HMP2 profile-to-metadata mapping is not one-to-one")
    if set(meta["External ID"]) != set(sample_ids):
        raise RuntimeError("One or more profile columns lack metagenomics metadata")

    meta["TechnicalRerun"] = meta["External ID"].str.endswith("_TR")
    meta["PreferredP"] = meta["External ID"].str.endswith("_P")
    meta["reads_filtered"] = pd.to_numeric(meta["reads_filtered"], errors="raise")
    ranked = meta.sort_values(
        ["site_sub_coll", "TechnicalRerun", "PreferredP", "reads_filtered", "External ID"],
        ascending=[True, True, False, False, True],
    ).copy()
    ranked["SelectionRank"] = ranked.groupby("site_sub_coll").cumcount() + 1
    ranked["PrimaryProfile"] = ranked["SelectionRank"].eq(1)
    ranked["SelectionReason"] = np.where(
        ranked["PrimaryProfile"],
        "Primary: non-TR, then _P, then highest filtered reads, then SampleID",
        "Excluded technical/processing duplicate",
    )
    selection = ranked[
        [
            "External ID",
            "site_sub_coll",
            "Participant ID",
            "diagnosis",
            "week_num",
            "reads_filtered",
            "TechnicalRerun",
            "PreferredP",
            "SelectionRank",
            "PrimaryProfile",
            "SelectionReason",
        ]
    ].rename(columns={"External ID": "SampleID", "Participant ID": "SubjectID"})

    primary = ranked.loc[ranked["PrimaryProfile"]].copy()
    depth_pass = primary.loc[primary["reads_filtered"].ge(READ_GATE)].copy()
    visit_counts = depth_pass.groupby("Participant ID")["External ID"].size()
    eligible_subjects = visit_counts.loc[visit_counts.ge(MIN_VISITS)].index
    eligible = depth_pass.loc[depth_pass["Participant ID"].isin(eligible_subjects)].copy()
    eligible["Diagnosis"] = eligible["diagnosis"].replace({"nonIBD": "Control"})
    eligible["SampleID"] = eligible["External ID"].astype(str)
    eligible["SubjectID"] = eligible["Participant ID"].astype(str)
    eligible["CollectionID"] = eligible["site_sub_coll"].astype(str)
    eligible["Week"] = pd.to_numeric(eligible["week_num"], errors="raise")
    eligible["FilteredReads"] = eligible["reads_filtered"].astype(float)
    eligible["Antibiotics"] = eligible["Antibiotics"].astype(str)
    eligible = eligible.sort_values(["SubjectID", "Week", "SampleID"])

    stages = []
    for order, (stage, frame) in enumerate(
        (
            ("Published MetaPhlAn profiles", meta),
            ("One primary profile per biospecimen", primary),
            ("Filtered reads ≥1,000,000", depth_pass),
            ("Subjects with ≥4 eligible visits", eligible),
        ),
        start=1,
    ):
        subject_col = "Participant ID"
        stages.append(
            {
                "StageOrder": order,
                "Stage": stage,
                "Profiles": len(frame),
                "Subjects": frame[subject_col].nunique(),
                "CD": int((frame["diagnosis"] == "CD").sum()),
                "UC": int((frame["diagnosis"] == "UC").sum()),
                "Control": int((frame["diagnosis"] == "nonIBD").sum()),
            }
        )
    return eligible, selection, pd.DataFrame(stages)


def species_matrix(
    profile: pd.DataFrame, eligible: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    taxa = profile["#SampleID"].astype(str)
    terminal_species = taxa.str.contains(r"(?:^|\|)s__", regex=True) & ~taxa.str.contains(
        r"\|t__", regex=True
    )
    species_source = profile.loc[terminal_species].copy()
    values = species_source.iloc[:, 1:].apply(pd.to_numeric, errors="raise")
    ever_observed = values.gt(0).any(axis=1).to_numpy()
    species_source = species_source.loc[ever_observed].reset_index(drop=True)
    values = values.loc[ever_observed].reset_index(drop=True)
    taxa = species_source["#SampleID"].astype(str)
    species_names = taxa.str.rsplit("|", n=1).str[-1]
    if species_names.duplicated().any():
        raise RuntimeError("Terminal species labels are not unique")

    sample_ids = eligible["SampleID"].tolist()
    matrix = values.loc[:, sample_ids].T
    matrix.columns = species_names
    raw_sums = matrix.sum(axis=1)
    if raw_sums.le(0).any():
        raise RuntimeError("An eligible sample has no terminal-species abundance")
    relative = matrix.div(raw_sums, axis=0)
    relative.index.name = "SampleID"
    feature_map = pd.DataFrame(
        {
            "FeatureID": [f"S{index:03d}" for index in range(1, len(taxa) + 1)],
            "Taxon": taxa,
            "Species": species_names,
        }
    )
    relative.columns = feature_map["FeatureID"]
    normalization = pd.DataFrame(
        {
            "SampleID": sample_ids,
            "TerminalSpeciesSumBeforeRenormalization": raw_sums.to_numpy(float),
            "TerminalSpeciesSumAfterRenormalization": relative.sum(axis=1).to_numpy(float),
        }
    )
    return relative, feature_map, normalization


def dysbiosis_analysis(
    eligible: pd.DataFrame,
    relative: pd.DataFrame,
    rng: np.random.Generator,
    iterations: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample = eligible[
        ["SampleID", "SubjectID", "CollectionID", "Diagnosis", "Week", "FilteredReads", "Antibiotics"]
    ].copy()
    sample = sample.set_index("SampleID").loc[relative.index].reset_index()
    reference_mask = sample["Diagnosis"].eq("Control") & sample["Week"].ge(20)
    reference_ids = sample.loc[reference_mask, "SampleID"].tolist()
    distances = cdist(relative.to_numpy(float), relative.loc[reference_ids].to_numpy(float), metric="braycurtis")
    reference_lookup = {sample_id: index for index, sample_id in enumerate(reference_ids)}
    scores = np.empty(len(sample), dtype=float)
    for row, sample_id in enumerate(sample["SampleID"]):
        values = distances[row].copy()
        if sample_id in reference_lookup:
            values[reference_lookup[sample_id]] = np.nan
        scores[row] = np.nanmedian(values)
    sample["DysbiosisScore"] = scores
    threshold = float(np.quantile(sample.loc[sample["Diagnosis"].eq("Control"), "DysbiosisScore"], 0.90))
    sample["Dysbiotic"] = sample["DysbiosisScore"].gt(threshold)
    sample["WeekYearCentered"] = (sample["Week"] - 26.0) / 52.0
    sample["Log10ReadsCentered"] = np.log10(sample["FilteredReads"])
    sample["Log10ReadsCentered"] -= sample["Log10ReadsCentered"].mean()

    summary_records = []
    for diagnosis in DIAGNOSIS_ORDER:
        group = sample.loc[sample["Diagnosis"].eq(diagnosis)]
        subject_scores = group.groupby("SubjectID")["DysbiosisScore"].median().to_numpy(float)
        median, lower, upper = bootstrap_ci(subject_scores, rng, iterations, "median")
        summary_records.append(
            {
                "Diagnosis": diagnosis,
                "Samples": len(group),
                "Subjects": group["SubjectID"].nunique(),
                "MedianDysbiosisScore": float(group["DysbiosisScore"].median()),
                "SubjectMedianScore": median,
                "SubjectMedianCILower": lower,
                "SubjectMedianCIUpper": upper,
                "DysbioticSamples": int(group["Dysbiotic"].sum()),
                "DysbioticSamplePercent": float(group["Dysbiotic"].mean() * 100),
            }
        )
    reference_audit = pd.DataFrame(
        [
            {
                "ReferenceDefinition": "Control samples collected at week ≥20",
                "ReferenceSamples": len(reference_ids),
                "ReferenceSubjects": sample.loc[reference_mask, "SubjectID"].nunique(),
                "ScoreDefinition": "Median Bray-Curtis distance to reference samples; self-distance excluded",
                "ThresholdDefinition": "90th percentile of all eligible control-sample scores",
                "DysbiosisThreshold": threshold,
            }
        ]
    )
    return sample, pd.DataFrame(summary_records), reference_audit


def longitudinal_distances(
    sample: pd.DataFrame,
    relative: pd.DataFrame,
    rng: np.random.Generator,
    iterations: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    consecutive_rows = []
    sample_lookup = sample.set_index("SampleID")
    for subject, visits in sample.sort_values(["SubjectID", "Week"]).groupby("SubjectID", sort=True):
        visits = visits.reset_index(drop=True)
        for first, second in itertools.combinations(range(len(visits)), 2):
            left = visits.iloc[first]
            right = visits.iloc[second]
            gap = float(right["Week"] - left["Week"])
            distance = float(braycurtis(relative.loc[left["SampleID"]], relative.loc[right["SampleID"]]))
            rows.append(
                {
                    "SubjectID": subject,
                    "Diagnosis": left["Diagnosis"],
                    "Sample1": left["SampleID"],
                    "Sample2": right["SampleID"],
                    "Week1": left["Week"],
                    "Week2": right["Week"],
                    "LagWeeks": gap,
                    "BrayCurtis": distance,
                }
            )
        for first in range(len(visits) - 1):
            left = visits.iloc[first]
            right = visits.iloc[first + 1]
            gap = float(right["Week"] - left["Week"])
            distance = float(braycurtis(relative.loc[left["SampleID"]], relative.loc[right["SampleID"]]))
            consecutive_rows.append(
                {
                    "SubjectID": subject,
                    "Diagnosis": left["Diagnosis"],
                    "Sample1": left["SampleID"],
                    "Sample2": right["SampleID"],
                    "Week1": left["Week"],
                    "Week2": right["Week"],
                    "LagWeeks": gap,
                    "BrayCurtis": distance,
                    "ShortInterval": 1 <= gap <= 3,
                    "Shift054": (1 <= gap <= 3) and distance > SHIFT_THRESHOLD,
                    "AntibioticsAtSecondVisit": sample_lookup.loc[right["SampleID"], "Antibiotics"],
                }
            )
    pairwise = pd.DataFrame(rows)
    # Different biospecimens occasionally share a rounded study week. They do not
    # provide a positive time lag, so they are retained in the sample ledger but
    # excluded from the lag-response estimand.
    pairwise = pairwise.loc[pairwise["LagWeeks"].gt(0)].copy()
    pairwise["LagBin"] = pd.cut(
        pairwise["LagWeeks"], bins=LAG_BINS, labels=LAG_LABELS, right=True
    ).astype("object")
    pairwise = pairwise.loc[pairwise["LagBin"].notna()].copy()
    subject_lag = (
        pairwise.groupby(["SubjectID", "Diagnosis", "LagBin"], observed=True)["BrayCurtis"]
        .agg(MedianBrayCurtis="median", PairCount="size")
        .reset_index()
    )
    lag_summary_records = []
    for diagnosis in DIAGNOSIS_ORDER:
        for lag in LAG_LABELS:
            values = subject_lag.loc[
                subject_lag["Diagnosis"].eq(diagnosis) & subject_lag["LagBin"].eq(lag),
                "MedianBrayCurtis",
            ].to_numpy(float)
            median, lower, upper = bootstrap_ci(values, rng, iterations, "median")
            lag_summary_records.append(
                {
                    "Diagnosis": diagnosis,
                    "LagBin": lag,
                    "Subjects": len(values),
                    "SubjectMedianBrayCurtis": median,
                    "CILower": lower,
                    "CIUpper": upper,
                }
            )
    consecutive = pd.DataFrame(consecutive_rows)
    short = consecutive.loc[consecutive["ShortInterval"]].copy()
    subject_shifts = (
        short.groupby(["SubjectID", "Diagnosis"])
        .agg(
            ShortIntervals=("Shift054", "size"),
            ShiftEvents=("Shift054", "sum"),
            MedianBrayCurtis=("BrayCurtis", "median"),
        )
        .reset_index()
    )
    subject_shifts["ShiftFraction"] = subject_shifts["ShiftEvents"] / subject_shifts["ShortIntervals"]
    shift_summary_records = []
    for diagnosis in DIAGNOSIS_ORDER:
        group = subject_shifts.loc[subject_shifts["Diagnosis"].eq(diagnosis)]
        mean, lower, upper = bootstrap_ci(group["ShiftFraction"].to_numpy(float), rng, iterations, "mean")
        shift_summary_records.append(
            {
                "Diagnosis": diagnosis,
                "Subjects": len(group),
                "ShortIntervals": int(group["ShortIntervals"].sum()),
                "ShiftEvents": int(group["ShiftEvents"].sum()),
                "PooledShiftPercent": float(group["ShiftEvents"].sum() / group["ShortIntervals"].sum() * 100),
                "MeanSubjectShiftFraction": mean,
                "MeanSubjectCILower": lower,
                "MeanSubjectCIUpper": upper,
                "MedianSubjectShiftFraction": float(group["ShiftFraction"].median()),
                "MedianSubjectBrayCurtis": float(group["MedianBrayCurtis"].median()),
            }
        )
    return (
        pairwise,
        pd.DataFrame(lag_summary_records),
        consecutive,
        subject_shifts,
        pd.DataFrame(shift_summary_records),
    )


def retention_analysis(
    consecutive: pd.DataFrame,
    relative: pd.DataFrame,
    rng: np.random.Generator,
    iterations: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts: dict[tuple[str, str, str], list[int]] = {}
    for row in consecutive.loc[consecutive["ShortInterval"]].itertuples(index=False):
        baseline = relative.loc[row.Sample1].to_numpy(float)
        followup = relative.loc[row.Sample2].to_numpy(float)
        bands = (
            (baseline >= DETECTION) & (baseline < 1e-3),
            (baseline >= 1e-3) & (baseline < 1e-2),
            baseline >= 1e-2,
        )
        for label, mask in zip(ABUNDANCE_LABELS, bands, strict=True):
            key = (row.SubjectID, row.Diagnosis, label)
            present, retained = counts.setdefault(key, [0, 0])
            counts[key][0] = present + int(mask.sum())
            counts[key][1] = retained + int(((followup >= DETECTION) & mask).sum())
    subject = pd.DataFrame(
        [
            {
                "SubjectID": subject_id,
                "Diagnosis": diagnosis,
                "BaselineAbundance": band,
                "BaselineDetections": present,
                "RetainedDetections": retained,
                "RetentionFraction": retained / present if present else np.nan,
            }
            for (subject_id, diagnosis, band), (present, retained) in counts.items()
        ]
    ).dropna(subset=["RetentionFraction"])
    records = []
    for diagnosis in DIAGNOSIS_ORDER:
        for band in ABUNDANCE_LABELS:
            values = subject.loc[
                subject["Diagnosis"].eq(diagnosis) & subject["BaselineAbundance"].eq(band),
                "RetentionFraction",
            ].to_numpy(float)
            median, lower, upper = bootstrap_ci(values, rng, iterations, "median")
            records.append(
                {
                    "Diagnosis": diagnosis,
                    "BaselineAbundance": band,
                    "Subjects": len(values),
                    "SubjectMedianRetention": median,
                    "CILower": lower,
                    "CIUpper": upper,
                }
            )
    return subject, pd.DataFrame(records)


def feature_analysis(
    sample: pd.DataFrame,
    relative: pd.DataFrame,
    feature_map: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prevalence = (relative >= DETECTION).mean(axis=0)
    mean = relative.mean(axis=0)
    audit = feature_map.copy()
    audit["Prevalence"] = audit["FeatureID"].map(prevalence)
    audit["MeanRelativeAbundance"] = audit["FeatureID"].map(mean)
    audit["SelectedForMixedModel"] = audit["Prevalence"].ge(0.20) & audit[
        "MeanRelativeAbundance"
    ].ge(0.001)
    selected_ids = audit.loc[audit["SelectedForMixedModel"], "FeatureID"].tolist()
    selected = relative.loc[:, selected_ids]
    closed = (selected + PSEUDOCOUNT).div((selected + PSEUDOCOUNT).sum(axis=1), axis=0)
    clr = np.log(closed).sub(np.log(closed).mean(axis=1), axis=0)
    clr.insert(0, "SampleID", clr.index)

    pcopri = audit.loc[audit["Species"].eq("s__Prevotella_copri")]
    if len(pcopri) != 1:
        raise RuntimeError("Prevotella copri feature contract failed")
    pcopri_id = pcopri.iloc[0]["FeatureID"]
    trajectory = sample[["SampleID", "SubjectID", "Diagnosis", "Week", "Antibiotics"]].copy()
    trajectory["RelativeAbundance"] = trajectory["SampleID"].map(relative[pcopri_id])
    subject_audit = (
        trajectory.groupby(["SubjectID", "Diagnosis"])["RelativeAbundance"]
        .agg(Visits="size", Prevalence=lambda x: float((x >= DETECTION).mean()), SD="std")
        .reset_index()
    )
    subject_audit = subject_audit.loc[
        subject_audit["Visits"].ge(8) & subject_audit["Prevalence"].ge(0.20)
    ].copy()
    subject_audit["SelectionRank"] = subject_audit.groupby("Diagnosis")["SD"].rank(
        method="first", ascending=False
    )
    selected_subjects = subject_audit.loc[subject_audit["SelectionRank"].le(4), "SubjectID"]
    trajectory = trajectory.loc[trajectory["SubjectID"].isin(selected_subjects)].merge(
        subject_audit[["SubjectID", "SelectionRank"]], on="SubjectID", how="left"
    )
    return audit, clr, trajectory


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    profile, metadata, manifest = load_sources(cache)
    eligible, selection, attrition = select_profiles(profile, metadata)
    relative, feature_map, normalization = species_matrix(profile, eligible)
    sample, dysbiosis_summary, reference_audit = dysbiosis_analysis(
        eligible, relative, rng, args.bootstrap_iterations
    )
    pairwise, lag_summary, consecutive, subject_shifts, shift_summary = longitudinal_distances(
        sample, relative, rng, args.bootstrap_iterations
    )
    retention_subject, retention_summary = retention_analysis(
        consecutive, relative, rng, args.bootstrap_iterations
    )
    feature_audit, clr, pcopri = feature_analysis(sample, relative, feature_map)

    species_export = relative.reset_index()
    outputs = {
        "profile-selection-ledger.tsv": selection,
        "sample-attrition.tsv": attrition,
        "normalization-audit.tsv": normalization,
        "species-feature-audit.tsv": feature_audit,
        "species-relative.tsv.gz": species_export,
        "sample-ledger.tsv": sample,
        "dysbiosis-summary.tsv": dysbiosis_summary,
        "dysbiosis-reference-audit.tsv": reference_audit,
        "within-subject-pairs.tsv.gz": pairwise,
        "lag-summary.tsv": lag_summary,
        "consecutive-intervals.tsv.gz": consecutive,
        "subject-shift-summary.tsv": subject_shifts,
        "shift-summary.tsv": shift_summary,
        "retention-subject-summary.tsv": retention_subject,
        "retention-summary.tsv": retention_summary,
        "selected-species-clr.tsv.gz": clr,
        "prevotella-selected-trajectories.tsv": pcopri,
    }
    for name, frame in outputs.items():
        write_tsv(frame, output / name)

    shutil.copy2(cache / "lloyd-price-fig3-original.png", output / "lloyd-price-fig3-original.png")
    (output / "source-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    software = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }
    (output / "software-versions-python.json").write_text(
        json.dumps(software, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    reference = reference_audit.iloc[0]
    metrics = {
        "article": 67,
        "seed": SEED,
        "plot_seed": PLOT_SEED,
        "bootstrap_iterations": args.bootstrap_iterations,
        "published_profiles": int(len(profile.columns) - 1),
        "published_subjects": int(selection["SubjectID"].nunique()),
        "primary_biospecimens": int(selection["PrimaryProfile"].sum()),
        "excluded_duplicate_profiles": int((~selection["PrimaryProfile"]).sum()),
        "eligible_samples": int(len(sample)),
        "eligible_subjects": int(sample["SubjectID"].nunique()),
        "eligible_by_diagnosis": {
            diagnosis: {
                "samples": int((sample["Diagnosis"] == diagnosis).sum()),
                "subjects": int(sample.loc[sample["Diagnosis"] == diagnosis, "SubjectID"].nunique()),
            }
            for diagnosis in DIAGNOSIS_ORDER
        },
        "terminal_species_rows": int(
            (
                profile["#SampleID"].astype(str).str.contains(r"(?:^|\|)s__", regex=True)
                & ~profile["#SampleID"].astype(str).str.contains(r"\|t__", regex=True)
            ).sum()
        ),
        "ever_observed_species": int(len(feature_map)),
        "selected_species": int(feature_audit["SelectedForMixedModel"].sum()),
        "dysbiosis_reference_samples": int(reference["ReferenceSamples"]),
        "dysbiosis_reference_subjects": int(reference["ReferenceSubjects"]),
        "dysbiosis_threshold": float(reference["DysbiosisThreshold"]),
        "short_intervals": int(consecutive["ShortInterval"].sum()),
        "shift_events": int(consecutive["Shift054"].sum()),
        "shift_threshold": SHIFT_THRESHOLD,
        "pcopri_selected_subjects": int(pcopri["SubjectID"].nunique()),
        "read_gate": READ_GATE,
        "minimum_visits": MIN_VISITS,
        "detection_threshold": DETECTION,
        "pseudocount": PSEUDOCOUNT,
    }
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    contract = {
        "primary_profile_rule": "within site_sub_coll: non-_TR first, then _P, then highest reads_filtered, then lexicographic SampleID",
        "eligibility": "reads_filtered >= 1,000,000 and subject has >=4 eligible primary visits",
        "composition": "terminal MetaPhlAn species rows, renormalized within sample",
        "dysbiosis_score": "median Bray-Curtis distance to week>=20 control reference profiles, excluding self",
        "dysbiosis_threshold": "90th percentile of all eligible control-sample dysbiosis scores",
        "short_interval": "consecutive eligible visits 1-3 weeks apart",
        "shift_sensitivity": "Bray-Curtis >0.54; declared sensitivity analysis, not the paper's lag-dependent shift call",
        "retention": "baseline species >=1e-4 remains >=1e-4 at next short-interval visit",
        "uncertainty_unit": "subject-level bootstrap",
        "feature_screen": "prevalence >=20% and mean relative abundance >=0.1%; screen declared before mixed modelling",
        "interpretation": "within-person association and detection stability; not proof of permanent colonization or causality",
    }
    (output / "methods-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
