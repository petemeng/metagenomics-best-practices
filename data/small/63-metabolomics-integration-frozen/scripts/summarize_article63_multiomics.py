#!/usr/bin/env python3
"""Summarize HAllA replication, global concordance, and DIABLO for Article 63."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def write_tsv(path: Path, frame: pd.DataFrame, compressed: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {"method": "gzip", "mtime": 0} if compressed else None
    frame.to_csv(
        path, sep="\t", index=False, lineterminator="\n", compression=compression
    )


def bh(values: np.ndarray) -> np.ndarray:
    return multipletests(values, alpha=0.05, method="fdr_bh")[1]


def load_feature_matrix(path: Path) -> pd.DataFrame:
    frame = read_tsv(path, index_col=0)
    if frame.index.duplicated().any() or frame.columns.duplicated().any():
        raise RuntimeError(f"Duplicate identifiers in {path}")
    return frame.astype(float)


def block_summary(branch: str, path: Path) -> pd.DataFrame:
    blocks = read_tsv(path)
    if blocks.empty:
        return pd.DataFrame(
            [{
                "Branch": branch, "Blocks": 0, "MedianMicrobes": 0,
                "MedianMetabolites": 0, "MaximumPairCells": 0,
                "BestAdjustedP": np.nan,
            }]
        )
    microbes = blocks.cluster_X.str.split(";").str.len()
    metabolites = blocks.cluster_Y.str.split(";").str.len()
    return pd.DataFrame(
        [{
            "Branch": branch,
            "Blocks": len(blocks),
            "MedianMicrobes": float(microbes.median()),
            "MedianMetabolites": float(metabolites.median()),
            "MaximumPairCells": int((microbes * metabolites).max()),
            "BestAdjustedP": float(blocks.best_adjusted_pvalue.min()),
        }]
    )


def main() -> None:
    args = parse_args()
    source = args.input_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    features = read_tsv(source / "feature-audit.tsv")
    feature_by_id = features.set_index("FeatureID", verify_integrity=True)
    adjusted = read_tsv(source / "halla-results/adjusted/all_associations.txt")
    raw = read_tsv(source / "halla-results/raw/all_associations.txt")
    for branch, frame in (("adjusted", adjusted), ("raw", raw)):
        if len(frame) != 166 * 153:
            raise RuntimeError(f"Unexpected {branch} HAllA pair count")
        if frame[["X_features", "Y_features"]].duplicated().any():
            raise RuntimeError(f"Duplicate HAllA pairs in {branch}")

    adjusted_pairs = set(
        zip(
            adjusted.loc[adjusted["q-values"].lt(0.05), "X_features"],
            adjusted.loc[adjusted["q-values"].lt(0.05), "Y_features"],
        )
    )
    raw_pairs = set(
        zip(
            raw.loc[raw["q-values"].lt(0.05), "X_features"],
            raw.loc[raw["q-values"].lt(0.05), "Y_features"],
        )
    )
    overlap = pd.DataFrame(
        [
            {"EvidenceClass": "Both", "Pairs": len(adjusted_pairs & raw_pairs)},
            {"EvidenceClass": "Adjusted only", "Pairs": len(adjusted_pairs - raw_pairs)},
            {"EvidenceClass": "Raw only", "Pairs": len(raw_pairs - adjusted_pairs)},
            {
                "EvidenceClass": "Neither",
                "Pairs": 166 * 153 - len(adjusted_pairs | raw_pairs),
            },
        ]
    )
    write_tsv(output / "halla-branch-overlap.tsv", overlap)

    validation_microbes = load_feature_matrix(
        source / "halla/validation-microbiome-adjusted.tsv"
    )
    validation_metabolites = load_feature_matrix(
        source / "halla/validation-metabolome-adjusted.tsv"
    )
    if not validation_microbes.columns.equals(validation_metabolites.columns):
        raise RuntimeError("Validation HAllA matrices are not sample aligned")

    discovery = adjusted.loc[adjusted["q-values"].lt(0.05)].copy()
    validation_rho = np.empty(len(discovery), dtype=float)
    validation_p = np.empty(len(discovery), dtype=float)
    for position, row in enumerate(discovery.itertuples(index=False)):
        result = spearmanr(
            validation_microbes.loc[row.X_features].to_numpy(float),
            validation_metabolites.loc[row.Y_features].to_numpy(float),
        )
        validation_rho[position] = result.statistic
        validation_p[position] = result.pvalue
    discovery = discovery.rename(
        columns={
            "X_features": "MicrobeID",
            "Y_features": "MetaboliteID",
            "association": "DiscoveryRho",
            "p-values": "DiscoveryP",
            "q-values": "DiscoveryQ",
        }
    )
    discovery["ValidationRho"] = validation_rho
    discovery["ValidationP"] = validation_p
    discovery["ValidationQ"] = bh(validation_p)
    discovery["SameDirection"] = (
        np.sign(discovery.DiscoveryRho) == np.sign(discovery.ValidationRho)
    )
    discovery["Replicated"] = discovery.SameDirection & discovery.ValidationQ.lt(0.05)

    microbe_labels = feature_by_id.loc[
        discovery.MicrobeID,
        ["DisplayName", "Phylum", "RawFeature"],
    ].reset_index(drop=True)
    metabolite_labels = feature_by_id.loc[
        discovery.MetaboliteID,
        ["DisplayName", "HMDB", "ChemicalClass", "RawFeature"],
    ].reset_index(drop=True)
    discovery["Microbe"] = microbe_labels.DisplayName.to_numpy()
    discovery["Phylum"] = microbe_labels.Phylum.to_numpy()
    discovery["MicrobeRawFeature"] = microbe_labels.RawFeature.to_numpy()
    discovery["Metabolite"] = metabolite_labels.DisplayName.to_numpy()
    discovery["HMDB"] = metabolite_labels.HMDB.to_numpy()
    discovery["ChemicalClass"] = metabolite_labels.ChemicalClass.to_numpy()
    discovery["MetaboliteRawFeature"] = metabolite_labels.RawFeature.to_numpy()
    discovery["CombinedMagnitude"] = (
        discovery.DiscoveryRho.abs() + discovery.ValidationRho.abs()
    ) / 2
    discovery = discovery.sort_values(
        ["Replicated", "ValidationQ", "CombinedMagnitude", "MicrobeID", "MetaboliteID"],
        ascending=[False, True, False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    write_tsv(output / "halla-pair-validation.tsv.gz", discovery, compressed=True)
    write_tsv(output / "top-replicated-pairs.tsv", discovery.loc[discovery.Replicated].head(30))

    replication_summary = pd.DataFrame(
        [
            {"Stage": "PRISM discovery BH q < 0.05", "Pairs": len(discovery)},
            {"Stage": "Same direction in Validation", "Pairs": int(discovery.SameDirection.sum())},
            {"Stage": "Validation nominal p < 0.05", "Pairs": int(discovery.ValidationP.lt(0.05).sum())},
            {"Stage": "Validation BH q < 0.05", "Pairs": int(discovery.ValidationQ.lt(0.05).sum())},
            {"Stage": "Same direction + Validation BH", "Pairs": int(discovery.Replicated.sum())},
        ]
    )
    write_tsv(output / "halla-replication-summary.tsv", replication_summary)

    block_table = pd.concat(
        [
            block_summary(
                "Covariate-adjusted",
                source / "halla-results/adjusted/sig_clusters.txt",
            ),
            block_summary(
                "Raw sensitivity",
                source / "halla-results/raw/sig_clusters.txt",
            ),
        ],
        ignore_index=True,
    )
    block_table["MarginalSignificantPairs"] = [len(adjusted_pairs), len(raw_pairs)]
    write_tsv(output / "halla-branch-summary.tsv", block_table)

    diablo_dir = source / "diablo"
    selected = read_tsv(diablo_dir / "final-selected-features.tsv")
    stability = read_tsv(diablo_dir / "bootstrap-feature-stability.tsv")
    selected_stability = selected.merge(
        stability[
            ["FeatureID", "Block", "Component", "SelectedBootstraps",
             "Bootstraps", "SelectionFrequency"]
        ],
        on=["FeatureID", "Block", "Component"],
        how="left",
        validate="one_to_one",
    )
    selected_stability["Stable70"] = selected_stability.SelectionFrequency.ge(0.70)
    selected_stability = selected_stability.sort_values(
        ["Block", "Component", "SelectionFrequency", "FeatureID"],
        ascending=[True, True, False, True],
    )
    write_tsv(output / "diablo-selected-stability.tsv", selected_stability)

    global_tests = read_tsv(source / "global/global-concordance.tsv")
    external_metrics = read_tsv(diablo_dir / "external-metrics.tsv")
    null_summary = read_tsv(diablo_dir / "label-permutation-summary.tsv")
    tuning = read_tsv(diablo_dir / "tuning-summary.tsv")
    latent = read_tsv(diablo_dir / "latent-correlations.tsv")
    metric_by_name = external_metrics.set_index("Metric").to_dict("index")
    metrics = {
        "article": 63,
        "samples": 220,
        "independent_subjects": 220,
        "discovery_samples": 155,
        "external_validation_samples": 65,
        "selected_microbes": 166,
        "selected_metabolites": 153,
        "halla_adjusted_significant_pairs": len(adjusted_pairs),
        "halla_raw_significant_pairs": len(raw_pairs),
        "halla_replicated_pairs": int(discovery.Replicated.sum()),
        "halla_same_direction_pairs": int(discovery.SameDirection.sum()),
        "diablo_external_balanced_accuracy": float(
            metric_by_name["BalancedAccuracy"]["Estimate"]
        ),
        "diablo_external_balanced_accuracy_low": float(
            metric_by_name["BalancedAccuracy"]["Low"]
        ),
        "diablo_external_balanced_accuracy_high": float(
            metric_by_name["BalancedAccuracy"]["High"]
        ),
        "diablo_external_macro_f1": float(metric_by_name["MacroF1"]["Estimate"]),
        "diablo_label_permutation_p": float(null_summary.EmpiricalP.iloc[0]),
        "diablo_stable_selected_features": int(selected_stability.Stable70.sum()),
        "diablo_final_selected_rows": len(selected_stability),
        "global_tests": global_tests.to_dict("records"),
        "diablo_tuning": tuning.to_dict("records"),
        "diablo_latent_correlations": latent.to_dict("records"),
        "versions": {
            package: importlib.metadata.version(package)
            for package in ("numpy", "pandas", "scipy", "statsmodels")
        },
    }
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
