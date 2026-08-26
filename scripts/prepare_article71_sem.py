#!/usr/bin/env python3
"""Prepare the real Franzosa antibiotic–microbiome–calprotectin SEM cohort."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


ARTICLE = 71
ANALYSIS_SEED = 71_001
PLOT_SEED = 20_260_771
PSEUDOCOUNT = 1e-6
EXPECTED_COMMIT = "89a519d8c832008fbc6e650453e83e2f04858d02"
FAECALIBACTERIUM_SUFFIX = ";g__Faecalibacterium"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
        float_format="%.15g",
    )


def z_from_reference(
    values: pd.Series, reference: pd.Series
) -> tuple[pd.Series, float, float]:
    mean = float(reference.mean())
    sd = float(reference.std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Reference standard deviation must be positive")
    return (values - mean) / sd, mean, sd


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    manifest = json.loads(
        (cache / "download-manifest.json").read_text(encoding="utf-8")
    )
    if manifest["article"] != ARTICLE:
        raise ValueError("Unexpected download manifest article")
    if manifest["repository_commit"] != EXPECTED_COMMIT:
        raise ValueError("Unexpected curated-resource commit")
    for name, record in manifest["resources"].items():
        path = cache / name
        if path.stat().st_size != record["bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"Checksum mismatch for {name}")

    metadata = pd.read_csv(cache / "metadata.tsv", sep="\t")
    genera = pd.read_csv(cache / "genera.tsv", sep="\t")
    if len(metadata) != 220 or len(genera) != 220:
        raise ValueError("Expected 220 paired subjects")
    if not metadata["Sample"].is_unique or not metadata["Subject"].is_unique:
        raise ValueError("Samples and subjects must be unique")
    if not genera["Sample"].is_unique:
        raise ValueError("Genus sample identifiers must be unique")
    feature_columns = [column for column in genera.columns if column != "Sample"]
    if len(feature_columns) != 11_720:
        raise ValueError("Unexpected genus feature count")
    faec_columns = [
        column
        for column in feature_columns
        if column.endswith(FAECALIBACTERIUM_SUFFIX)
    ]
    if len(faec_columns) != 1:
        raise ValueError(f"Expected one Faecalibacterium column, got {faec_columns}")

    abundance = genera[feature_columns].to_numpy(dtype=float)
    if not np.isfinite(abundance).all() or (abundance < 0).any():
        raise ValueError("Genus matrix must be finite and nonnegative")
    row_sums = abundance.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-10):
        raise ValueError("Genus profiles do not sum to one")
    positive = abundance > 0
    safe = np.where(positive, abundance, 1.0)
    shannon = -(np.where(positive, abundance * np.log(safe), 0.0)).sum(axis=1)
    richness = positive.sum(axis=1)
    faec_index = feature_columns.index(faec_columns[0])
    faecalibacterium = abundance[:, faec_index]

    microbial = pd.DataFrame(
        {
            "Sample": genera["Sample"],
            "ProfileSum": row_sums,
            "DetectedGenera": richness,
            "Shannon": shannon,
            "FaecalibacteriumRelative": faecalibacterium,
            "Log10Faecalibacterium": np.log10(
                faecalibacterium + PSEUDOCOUNT
            ),
        }
    )
    cohort = metadata.merge(microbial, on="Sample", how="inner", validate="one_to_one")
    if len(cohort) != 220:
        raise ValueError("Metadata–profile join lost samples")
    cohort["Cohort"] = cohort["Sample"].str.split(".").str[0]
    cohort["Diagnosis"] = cohort["Study.Group"]
    cohort["LogCalprotectin"] = np.log1p(cohort["Fecal.Calprotectin"])
    for source, target in (
        ("antibiotic", "Antibiotic"),
        ("immunosuppressant", "Immunosuppressant"),
        ("mesalamine", "Mesalamine"),
        ("steroids", "Steroids"),
    ):
        cohort[target] = cohort[source].map({"No": 0, "Yes": 1})
    cohort["CD"] = cohort["Diagnosis"].eq("CD").astype(int)
    cohort["UC"] = cohort["Diagnosis"].eq("UC").astype(int)

    required = [
        "Fecal.Calprotectin",
        "Antibiotic",
        "Immunosuppressant",
        "Mesalamine",
        "Steroids",
        "Age",
    ]
    primary_mask = cohort["Cohort"].eq("PRISM") & cohort[required].notna().all(axis=1)
    validation_mask = cohort["Cohort"].eq("Validation") & cohort[required].notna().all(axis=1)
    validation_broad_mask = (
        cohort["Cohort"].eq("Validation")
        & cohort[["Fecal.Calprotectin", "Age"]].notna().all(axis=1)
    )
    primary = cohort.loc[primary_mask].copy()
    validation = cohort.loc[validation_mask].copy()
    if len(primary) != 90 or int(primary["Antibiotic"].sum()) != 13:
        raise ValueError("Unexpected PRISM primary cohort")
    if len(validation) != 38 or int(validation["Antibiotic"].sum()) != 0:
        raise ValueError("Unexpected strict Validation cohort")
    if int(validation_broad_mask.sum()) != 60:
        raise ValueError("Unexpected broad Validation phenotype cohort")

    standardization: list[dict[str, object]] = []
    for source, target in (
        ("Age", "AgeZ"),
        ("Shannon", "ShannonZ"),
        ("LogCalprotectin", "LogCalprotectinZ"),
        ("Log10Faecalibacterium", "LogFaecalibacteriumZ"),
    ):
        reference = primary[source]
        cohort[target], mean, sd = z_from_reference(cohort[source], reference)
        standardization.append(
            {
                "Variable": source,
                "ReferenceCohort": "PRISM primary complete cases",
                "Mean": mean,
                "SD": sd,
                "TransformedColumn": target,
            }
        )
    primary = cohort.loc[primary_mask].copy()
    validation = cohort.loc[validation_mask].copy()

    columns = [
        "Sample",
        "Subject",
        "Cohort",
        "Diagnosis",
        "Age",
        "AgeZ",
        "Antibiotic",
        "Immunosuppressant",
        "Mesalamine",
        "Steroids",
        "CD",
        "UC",
        "Fecal.Calprotectin",
        "LogCalprotectin",
        "LogCalprotectinZ",
        "ProfileSum",
        "DetectedGenera",
        "Shannon",
        "ShannonZ",
        "FaecalibacteriumRelative",
        "Log10Faecalibacterium",
        "LogFaecalibacteriumZ",
    ]
    write_tsv(primary[columns], output / "sem-primary-cohort.tsv")
    write_tsv(validation[columns], output / "sem-validation-cohort.tsv")
    write_tsv(
        cohort[
            [
                "Sample",
                "Subject",
                "Cohort",
                "Diagnosis",
                "Fecal.Calprotectin",
                "antibiotic",
                "Shannon",
                "DetectedGenera",
                "FaecalibacteriumRelative",
            ]
        ],
        output / "all-sample-metrics.tsv",
    )
    write_tsv(
        pd.DataFrame(standardization),
        output / "standardization-parameters.tsv",
    )

    attrition = pd.DataFrame(
        [
            (1, "Public paired profiles and metadata", 220),
            (
                2,
                "Fecal calprotectin available",
                int(cohort["Fecal.Calprotectin"].notna().sum()),
            ),
            (
                3,
                "Calprotectin and antibiotic available",
                int(
                    cohort[
                        ["Fecal.Calprotectin", "Antibiotic"]
                    ].notna().all(axis=1).sum()
                ),
            ),
            (
                4,
                "All prespecified variables available",
                int(cohort[required].notna().all(axis=1).sum()),
            ),
            (5, "PRISM primary model cohort", len(primary)),
            (6, "Validation same-variable complete cases", len(validation)),
        ],
        columns=["StageOrder", "Stage", "Subjects"],
    )
    write_tsv(attrition, output / "sample-attrition.tsv")

    exposure = (
        cohort.groupby(["Cohort", "Diagnosis"], observed=True)
        .agg(
            Subjects=("Subject", "size"),
            AntibioticKnown=("Antibiotic", lambda x: int(x.notna().sum())),
            AntibioticExposed=("Antibiotic", lambda x: int(x.fillna(0).sum())),
            CalprotectinKnown=("Fecal.Calprotectin", lambda x: int(x.notna().sum())),
        )
        .reset_index()
    )
    write_tsv(exposure, output / "exposure-by-cohort-diagnosis.tsv")

    variable_rows: list[dict[str, object]] = []
    for label, frame in (("PRISM primary", primary), ("Validation strict", validation)):
        for variable in (
            "Age",
            "Fecal.Calprotectin",
            "LogCalprotectin",
            "Shannon",
            "DetectedGenera",
            "FaecalibacteriumRelative",
        ):
            values = frame[variable].dropna().astype(float)
            variable_rows.append(
                {
                    "Cohort": label,
                    "Variable": variable,
                    "N": len(values),
                    "Minimum": values.min(),
                    "Median": values.median(),
                    "Mean": values.mean(),
                    "SD": values.std(ddof=1),
                    "Maximum": values.max(),
                }
            )
    write_tsv(pd.DataFrame(variable_rows), output / "variable-summary.tsv")

    feature_contract = pd.DataFrame(
        [
            {
                "Node": "Environment",
                "Variable": "Antibiotic",
                "Definition": "metadata Yes versus No at the profiled visit",
                "Scale": "binary 1 versus 0",
                "TimingLimit": "collection timing relative to treatment initiation is unavailable",
            },
            {
                "Node": "Microbiome",
                "Variable": "ShannonZ",
                "Definition": "Shannon entropy over all 11,720 genus-profile columns",
                "Scale": "PRISM-primary standard deviations",
                "TimingLimit": "same stool/sample window as phenotype metadata",
            },
            {
                "Node": "Phenotype",
                "Variable": "LogCalprotectinZ",
                "Definition": "log1p fecal calprotectin",
                "Scale": "PRISM-primary standard deviations",
                "TimingLimit": "same cross-sectional visit; no temporal mediation ordering",
            },
        ]
    )
    write_tsv(feature_contract, output / "node-contract.tsv")

    shutil.copy2(
        cache / "franzosa-fig1-original.png",
        output / "franzosa-fig1-original.png",
    )
    source_manifest = dict(manifest)
    source_manifest.update(
        {
            "profile_rows": len(genera),
            "genus_features": len(feature_columns),
            "independent_subjects": 220,
            "primary_subjects": len(primary),
            "validation_strict_subjects": len(validation),
            "validation_broad_subjects": int(validation_broad_mask.sum()),
            "faecalibacterium_feature": faec_columns[0],
        }
    )
    (output / "source-manifest.json").write_text(
        json.dumps(source_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    methods_contract = {
        "article": ARTICLE,
        "design": "cross-sectional piecewise structural equation model",
        "independent_unit": "one subject and one paired stool profile per row",
        "primary_cohort": "PRISM complete cases",
        "environment": "contemporaneous antibiotic metadata, Yes versus No",
        "microbiome": "Shannon entropy over the complete genus relative-abundance profile",
        "phenotype": "log1p fecal calprotectin",
        "adjustment": [
            "Crohn disease indicator",
            "ulcerative colitis indicator",
            "age z-score",
            "immunosuppressant",
            "mesalamine",
            "steroids",
        ],
        "primary_graph": (
            "antibiotic -> Shannon; Shannon -> calprotectin; "
            "antibiotic -> calprotectin; covariates -> both endogenous nodes"
        ),
        "bootstrap": 5000,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "interpretation_limit": (
            "Exposure, stool microbiome, and calprotectin lack a verified temporal "
            "ordering; path coefficients are conditional associations, not identified "
            "causal mediation effects."
        ),
    }
    (output / "methods-contract.json").write_text(
        json.dumps(methods_contract, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    metrics = {
        "article": ARTICLE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "public_subjects": 220,
        "genus_features": 11_720,
        "calprotectin_available": int(
            cohort["Fecal.Calprotectin"].notna().sum()
        ),
        "primary_subjects": len(primary),
        "primary_antibiotic_exposed": int(primary["Antibiotic"].sum()),
        "primary_controls": int(primary["Diagnosis"].eq("Control").sum()),
        "primary_cd": int(primary["Diagnosis"].eq("CD").sum()),
        "primary_uc": int(primary["Diagnosis"].eq("UC").sum()),
        "primary_exposed_controls": int(
            (
                primary["Diagnosis"].eq("Control")
                & primary["Antibiotic"].eq(1)
            ).sum()
        ),
        "validation_strict_subjects": len(validation),
        "validation_antibiotic_exposed": int(validation["Antibiotic"].sum()),
        "pseudocount": PSEUDOCOUNT,
    }
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    software = {
        "Python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    (output / "software-versions-python.json").write_text(
        json.dumps(software, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
