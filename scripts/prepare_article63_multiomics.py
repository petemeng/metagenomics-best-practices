#!/usr/bin/env python3
"""Prepare paired, leakage-audited microbiome-metabolome inputs for Article 63."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 63001
PSEUDOCOUNT = 1e-6
MICROBE_PREVALENCE = 0.20
MICROBE_MEAN_ABUNDANCE = 5e-4
EXPECTED_GROUPS = {
    ("PRISM", "CD"): 68,
    ("PRISM", "UC"): 53,
    ("PRISM", "Control"): 34,
    ("Validation", "CD"): 20,
    ("Validation", "UC"): 23,
    ("Validation", "Control"): 22,
}
MEDICATIONS = ("antibiotic", "immunosuppressant", "mesalamine", "steroids")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def write_tsv(path: Path, frame: pd.DataFrame, index: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(
        path,
        sep="\t",
        index=index,
        lineterminator="\n",
        compression=compression,
    )


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def cohort_from_sample(sample: str) -> str:
    cohort = str(sample).split(".", 1)[0]
    if cohort not in {"PRISM", "Validation"}:
        raise ValueError(f"Unexpected sample prefix: {sample}")
    return cohort


def taxonomy_label(taxonomy: str) -> tuple[str, str, str]:
    parts = {part[:3]: part[3:] for part in str(taxonomy).split(";") if len(part) >= 3}
    genus = parts.get("g__", "")
    family = parts.get("f__", "")
    phylum = parts.get("p__", "")
    label = genus or family or phylum or "Unclassified"
    return genus or "Unclassified", phylum or "Unclassified", label


def centered_log_ratio(frame: pd.DataFrame, pseudocount: float) -> pd.DataFrame:
    logged = np.log(frame.astype(float) + pseudocount)
    return logged.sub(logged.mean(axis=1), axis=0)


def design_matrix(metadata: pd.DataFrame) -> pd.DataFrame:
    design = pd.DataFrame(index=metadata.index)
    design["Intercept"] = 1.0
    age = pd.to_numeric(metadata["Age"], errors="coerce")
    age = age.fillna(age.median())
    scale = age.std(ddof=0)
    design["Age_z"] = (age - age.mean()) / (scale if scale > 0 else 1.0)
    group = pd.Categorical(metadata["Study.Group"], categories=["Control", "CD", "UC"])
    for level in ("CD", "UC"):
        design[f"Group_{level}"] = (group == level).astype(float)
    for medication in MEDICATIONS:
        values = metadata[medication].fillna("Unknown").astype(str)
        for level in ("Yes", "Unknown"):
            column = (values == level).astype(float)
            if column.nunique() > 1:
                design[f"{medication}_{level}"] = column
    return design.astype(float)


def residualize(frame: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    design = design_matrix(metadata)
    values = frame.loc[metadata.index].to_numpy(float)
    coefficients, *_ = np.linalg.lstsq(design.to_numpy(float), values, rcond=None)
    residuals = values - design.to_numpy(float) @ coefficients
    scale = residuals.std(axis=0, ddof=1)
    if np.any(scale <= 0) or not np.isfinite(scale).all():
        bad = frame.columns[(scale <= 0) | ~np.isfinite(scale)].tolist()
        raise RuntimeError(f"Residualization produced constant features: {bad[:10]}")
    residuals = (residuals - residuals.mean(axis=0)) / scale
    result = pd.DataFrame(residuals, index=metadata.index, columns=frame.columns)
    return result, design


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(cache / "metadata.tsv", sep="\t", dtype={"Sample": str})
    microbes = pd.read_csv(cache / "genera.tsv", sep="\t", dtype={"Sample": str})
    metabolites = pd.read_csv(cache / "mtb.tsv", sep="\t", dtype={"Sample": str})
    metabolite_map = pd.read_csv(cache / "mtb.map.tsv", sep="\t")

    for name, frame in {
        "metadata": metadata,
        "microbiome": microbes,
        "metabolome": metabolites,
    }.items():
        if frame["Sample"].duplicated().any():
            raise RuntimeError(f"Duplicate sample IDs in {name}")
    sample_order = metadata["Sample"].tolist()
    if microbes["Sample"].tolist() != sample_order or metabolites["Sample"].tolist() != sample_order:
        raise RuntimeError("The three source tables are not identically sample aligned")
    if metadata["Subject"].astype(str).duplicated().any():
        raise RuntimeError("Article 63 requires one sample per independent subject")

    metadata = metadata.set_index("Sample", drop=False)
    metadata["Cohort"] = [cohort_from_sample(sample) for sample in metadata.index]
    observed_groups = metadata.groupby(["Cohort", "Study.Group"], observed=True).size().to_dict()
    if observed_groups != EXPECTED_GROUPS:
        raise RuntimeError(f"Unexpected cohort/group counts: {observed_groups}")

    microbes = microbes.set_index("Sample").astype(float)
    metabolites = metabolites.set_index("Sample").astype(float)
    if (microbes < 0).any().any() or (metabolites < 0).any().any():
        raise RuntimeError("Negative source values are not permitted")
    closure = microbes.sum(axis=1)
    if not np.allclose(closure, 1.0, atol=1e-10, rtol=1e-10):
        raise RuntimeError("The curated genus table does not close to one per sample")
    if set(metabolite_map["Compound"].astype(str)) != set(metabolites.columns):
        raise RuntimeError("Metabolite mapping and abundance feature universes differ")

    training_samples = metadata.index[metadata["Cohort"].eq("PRISM")]
    training_microbes = microbes.loc[training_samples]
    microbe_prevalence = training_microbes.gt(0).mean(axis=0)
    microbe_mean = training_microbes.mean(axis=0)
    selected_microbes = sorted(
        microbe_prevalence.index[
            microbe_prevalence.ge(MICROBE_PREVALENCE)
            & microbe_mean.ge(MICROBE_MEAN_ABUNDANCE)
        ]
    )
    if not selected_microbes:
        raise RuntimeError("No microbes passed the prespecified discovery-cohort filter")

    metabolite_map = metabolite_map.set_index("Compound", drop=False)
    high_confidence = (
        metabolite_map["High.Confidence.Annotation"]
        .astype(str)
        .str.strip()
        .str.lower()
        .eq("true")
    )
    selected_metabolites = sorted(
        metabolite_map.index[high_confidence & metabolite_map["HMDB"].notna()].astype(str)
    )
    if len(selected_metabolites) != 153:
        raise RuntimeError(
            f"Expected 153 high-confidence HMDB features, observed {len(selected_metabolites)}"
        )
    if metabolite_map.loc[selected_metabolites, "HMDB"].duplicated().any():
        raise RuntimeError("Primary metabolite features must have unique HMDB identifiers")

    microbe_id = {feature: f"MB{index:03d}" for index, feature in enumerate(selected_microbes, 1)}
    metabolite_id = {
        feature: f"MT{index:03d}" for index, feature in enumerate(selected_metabolites, 1)
    }
    microbiome_relative = microbes[selected_microbes].rename(columns=microbe_id)
    microbiome_clr = centered_log_ratio(microbiome_relative, PSEUDOCOUNT)
    metabolome_intensity = metabolites[selected_metabolites].rename(columns=metabolite_id)
    metabolome_log1p = np.log1p(metabolome_intensity)
    if not np.allclose(microbiome_clr.sum(axis=1), 0.0, atol=1e-10):
        raise RuntimeError("CLR row sums do not close to zero")

    feature_rows: list[dict[str, object]] = []
    for raw in selected_microbes:
        genus, phylum, label = taxonomy_label(raw)
        feature_rows.append(
            {
                "Modality": "Microbiome",
                "FeatureID": microbe_id[raw],
                "DisplayName": label,
                "Genus": genus,
                "Phylum": phylum,
                "HMDB": "",
                "KEGG": "",
                "ChemicalClass": "",
                "RawFeature": raw,
                "TrainingPrevalence": microbe_prevalence[raw],
                "TrainingMean": microbe_mean[raw],
                "TrainingVariance": training_microbes[raw].var(ddof=1),
                "HighConfidenceAnnotation": True,
            }
        )
    for raw in selected_metabolites:
        row = metabolite_map.loc[raw]
        feature_rows.append(
            {
                "Modality": "Metabolome",
                "FeatureID": metabolite_id[raw],
                "DisplayName": str(row["Compound.Name"]),
                "Genus": "",
                "Phylum": "",
                "HMDB": str(row["HMDB"]),
                "KEGG": "" if pd.isna(row["KEGG"]) else str(row["KEGG"]),
                "ChemicalClass": "" if pd.isna(row["Putative.Chemical.Class"]) else str(row["Putative.Chemical.Class"]),
                "RawFeature": raw,
                "TrainingPrevalence": metabolites.loc[training_samples, raw].gt(0).mean(),
                "TrainingMean": metabolites.loc[training_samples, raw].mean(),
                "TrainingVariance": np.log1p(metabolites.loc[training_samples, raw]).var(ddof=1),
                "HighConfidenceAnnotation": bool(row["High.Confidence.Annotation"]),
            }
        )
    features = pd.DataFrame(feature_rows)

    selected_metadata = metadata[
        [
            "Sample", "Subject", "Cohort", "Study.Group", "Age",
            "Fecal.Calprotectin", *MEDICATIONS,
        ]
    ].reset_index(drop=True)
    write_tsv(output / "sample-metadata.tsv", selected_metadata)
    write_tsv(output / "feature-audit.tsv", features)
    for name, frame in {
        "microbiome-relative.tsv.gz": microbiome_relative,
        "microbiome-clr.tsv.gz": microbiome_clr,
        "metabolome-intensity.tsv.gz": metabolome_intensity,
        "metabolome-log1p.tsv.gz": metabolome_log1p,
    }.items():
        export = frame.copy()
        export.insert(0, "Sample", export.index)
        write_tsv(output / name, export.reset_index(drop=True))

    design_rows: list[pd.DataFrame] = []
    for cohort in ("PRISM", "Validation"):
        cohort_metadata = metadata.loc[metadata["Cohort"].eq(cohort)]
        mb_adjusted, design = residualize(microbiome_clr, cohort_metadata)
        mt_adjusted, _ = residualize(metabolome_log1p, cohort_metadata)
        mb_raw = microbiome_clr.loc[cohort_metadata.index]
        mt_raw = metabolome_log1p.loc[cohort_metadata.index]
        for branch, mb_frame, mt_frame in (
            ("raw", mb_raw, mt_raw),
            ("adjusted", mb_adjusted, mt_adjusted),
        ):
            write_tsv(
                output / "halla" / f"{cohort.lower()}-microbiome-{branch}.tsv",
                mb_frame.T,
                index=True,
            )
            write_tsv(
                output / "halla" / f"{cohort.lower()}-metabolome-{branch}.tsv",
                mt_frame.T,
                index=True,
            )
        current_design = design.copy()
        current_design.insert(0, "Sample", current_design.index)
        current_design.insert(1, "Cohort", cohort)
        design_rows.append(current_design.reset_index(drop=True))
    write_tsv(output / "covariate-design.tsv", pd.concat(design_rows, ignore_index=True))

    attrition = pd.DataFrame(
        [
            {"Modality": "Microbiome", "Stage": "Source features", "Features": microbes.shape[1]},
            {
                "Modality": "Microbiome",
                "Stage": "PRISM prevalence ≥20%",
                "Features": int(microbe_prevalence.ge(MICROBE_PREVALENCE).sum()),
            },
            {
                "Modality": "Microbiome",
                "Stage": "Mean abundance ≥0.05%",
                "Features": len(selected_microbes),
            },
            {"Modality": "Metabolome", "Stage": "Source peaks", "Features": metabolites.shape[1]},
            {
                "Modality": "Metabolome",
                "Stage": "Named compounds",
                "Features": int(metabolite_map["Compound.Name"].notna().sum()),
            },
            {
                "Modality": "Metabolome",
                "Stage": "High-confidence HMDB",
                "Features": len(selected_metabolites),
            },
        ]
    )
    write_tsv(output / "feature-attrition.tsv", attrition)

    contract = {
        "article": 63,
        "seed": SEED,
        "samples": len(metadata),
        "independent_subjects": metadata["Subject"].nunique(),
        "discovery_samples": int(metadata["Cohort"].eq("PRISM").sum()),
        "external_validation_samples": int(metadata["Cohort"].eq("Validation").sum()),
        "source_microbe_features": microbes.shape[1],
        "source_metabolite_features": metabolites.shape[1],
        "selected_microbes": len(selected_microbes),
        "selected_metabolites": len(selected_metabolites),
        "microbe_prevalence_threshold": MICROBE_PREVALENCE,
        "microbe_mean_abundance_threshold": MICROBE_MEAN_ABUNDANCE,
        "microbiome_transform": f"CLR after fixed pseudocount {PSEUDOCOUNT:g}",
        "metabolome_transform": "log1p of author-processed nonnegative intensity",
        "metabolite_annotation_gate": "High.Confidence.Annotation == TRUE and HMDB non-missing",
        "feature_filter_fit_on": "PRISM discovery cohort only",
        "halla_primary": "within-cohort covariate residuals",
        "halla_covariates": ["Study.Group", "Age", *MEDICATIONS],
        "diablo_training": "PRISM discovery cohort only",
        "diablo_external_test": "Validation cohort, untouched until final model",
    }
    (output / "analysis-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    versions = {
        package: package_version(package)
        for package in ("pandas", "numpy", "scipy")
    }
    versions["seed"] = SEED
    (output / "prepare-software-versions.json").write_text(
        json.dumps(versions, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
