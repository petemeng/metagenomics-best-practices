#!/usr/bin/env python3
"""Prepare a leakage-auditable shotgun microbiome survival cohort for Article 68."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


SEED = 68_001
PLOT_SEED = 20_260_768
PSEUDOCOUNT_PPM = 25.0
QUALITY_GB_GATE = 0.5
QUALITY_ASSEMBLED_GATE = 50.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        compression="gzip" if path.suffix == ".gz" else None,
    )


def safe_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace({"N_A": np.nan, "": np.nan}), errors="coerce")


def extract_anchor(paper: Path, target: Path) -> None:
    """Rasterize page 7 and crop the published Figure 3A Kaplan-Meier panel."""
    with tempfile.TemporaryDirectory(prefix="article68-anchor-") as temporary:
        prefix = Path(temporary) / "page"
        subprocess.run(
            [
                "pdftoppm",
                "-f",
                "7",
                "-l",
                "7",
                "-singlefile",
                "-png",
                "-r",
                "240",
                str(paper),
                str(prefix),
            ],
            check=True,
            capture_output=True,
        )
        page = Image.open(prefix.with_suffix(".png")).convert("RGB")
        width, height = page.size
        crop = page.crop(
            (
                round(width * 0.055),
                round(height * 0.045),
                round(width * 0.452),
                round(height * 0.315),
            )
        )
        crop.save(target, format="PNG", optimize=True)


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((cache / "download-manifest.json").read_text(encoding="utf-8"))
    if manifest.get("article") != 68:
        raise RuntimeError("Article identity mismatch in download manifest")

    workbook = cache / "spencer-human-wgs.xlsx"
    metadata = pd.read_excel(workbook, sheet_name="Metadata", dtype=str, keep_default_na=False)
    ppm = pd.read_excel(workbook, sheet_name="LKT_PPM").set_index("row.names")
    feature = pd.read_excel(workbook, sheet_name="LKT_featuretable").set_index("row.names")
    ppm = ppm.apply(pd.to_numeric, errors="raise")
    feature = feature.reindex(ppm.index)

    if len(metadata) != 167 or metadata["SRA_accession"].nunique() != 167:
        raise RuntimeError("Unexpected Spencer WGS metadata dimensions")
    if ppm.shape != (225, 167) or set(ppm.columns) != set(metadata["SRA_accession"]):
        raise RuntimeError("Unexpected LKT matrix dimensions or sample identifiers")
    if feature.isna().all(axis=1).any():
        raise RuntimeError("One or more filtered LKT rows lack taxonomy")

    pfs_days = numeric(metadata["pfs_d"])
    event = numeric(metadata["pfsevent"])
    complete_pfs = pfs_days.notna() & event.isin([0, 1]) & pfs_days.gt(0)
    eligible_meta = metadata.loc[complete_pfs].copy()
    eligible_meta["PFS_days"] = pfs_days.loc[complete_pfs].astype(int)
    eligible_meta["Event"] = event.loc[complete_pfs].astype(int)
    if len(eligible_meta) != 110 or int(eligible_meta["Event"].sum()) != 61:
        raise RuntimeError("Unexpected complete-PFS WGS cohort")

    sample_ids = eligible_meta["SRA_accession"].tolist()
    eligible_ppm = ppm.loc[:, sample_ids].copy()
    totals = eligible_ppm.sum(axis=0)
    if totals.le(0).any():
        raise RuntimeError("A survival sample has no reported LKT abundance")
    relative = eligible_ppm.div(totals, axis=1)
    positive = eligible_ppm.gt(0).sum(axis=0)
    richness_50 = eligible_ppm.ge(50).sum(axis=0)
    shannon = -(relative.where(relative.gt(0)) * np.log(relative.where(relative.gt(0)))).sum(axis=0)
    inverse_simpson = 1.0 / relative.pow(2).sum(axis=0)

    faec_mask = feature["Genus"].eq("g__Faecalibacterium")
    faec_rows = feature.index[faec_mask]
    if len(faec_rows) != 3:
        raise RuntimeError(f"Expected three Faecalibacterium LKT rows, observed {len(faec_rows)}")
    faecalibacterium = eligible_ppm.loc[faec_rows].sum(axis=0)
    faec_log2 = np.log2(faecalibacterium + PSEUDOCOUNT_PPM)

    meta_index = eligible_meta.set_index("SRA_accession")
    cohort = pd.DataFrame(index=sample_ids)
    cohort.index.name = "SampleID"
    cohort["PFS_days"] = meta_index.loc[sample_ids, "PFS_days"]
    cohort["PFS_months"] = cohort["PFS_days"] / 30.4375
    cohort["Event"] = meta_index.loc[sample_ids, "Event"]
    cohort["Age"] = numeric(meta_index.loc[sample_ids, "Age"])
    cohort["Sex"] = meta_index.loc[sample_ids, "Sex"]
    cohort["BMI"] = numeric(meta_index.loc[sample_ids, "BMI"])
    cohort["PrimarySubtype"] = meta_index.loc[sample_ids, "primary_cat"]
    cohort["AdvancedSubstage"] = meta_index.loc[sample_ids, "adv_substage"]
    cohort["LDH"] = meta_index.loc[sample_ids, "LDH"]
    cohort["Treatment"] = meta_index.loc[sample_ids, "treatment"]
    cohort["TreatmentNaive"] = meta_index.loc[sample_ids, "Treatment_naive"]
    cohort["Response"] = meta_index.loc[sample_ids, "response"]
    cohort["FiberCategory"] = meta_index.loc[sample_ids, "fiber_cat"]
    cohort["Probiotics"] = meta_index.loc[sample_ids, "probiotics"]
    cohort["Antibiotics"] = meta_index.loc[sample_ids, "antibiotics"]
    cohort["GbNonHuman"] = numeric(meta_index.loc[sample_ids, "GbNAHS"])
    cohort["PercentAssembled"] = numeric(meta_index.loc[sample_ids, "PctAss"])
    cohort["ReportedLKTPPM"] = totals.loc[sample_ids]
    cohort["PositiveLKTFeatures"] = positive.loc[sample_ids]
    cohort["Richness50PPM"] = richness_50.loc[sample_ids]
    cohort["ShannonReportedLKT"] = shannon.loc[sample_ids]
    cohort["InverseSimpsonReportedLKT"] = inverse_simpson.loc[sample_ids]
    cohort["FaecalibacteriumPPM"] = faecalibacterium.loc[sample_ids]
    cohort["FaecalibacteriumLog2"] = faec_log2.loc[sample_ids]
    cohort["QualitySensitivityPass"] = (
        cohort["GbNonHuman"].ge(QUALITY_GB_GATE)
        & cohort["PercentAssembled"].ge(QUALITY_ASSEMBLED_GATE)
    )
    cohort = cohort.reset_index()

    primary_covariates = ["BMI", "PrimarySubtype", "AdvancedSubstage", "LDH"]
    if cohort[primary_covariates].isna().any().any():
        raise RuntimeError("Primary paper-aligned adjustment variables must be complete")
    expected_levels = {
        "PrimarySubtype": {"Cutaneous_or_unknown", "Mucosal_or_acral"},
        "AdvancedSubstage": {"Stage_M1C", "Stage_M1D"},
        "LDH": {"No", "Yes"},
    }
    for column, expected in expected_levels.items():
        observed = set(cohort[column])
        if observed != expected:
            raise RuntimeError(f"Unexpected {column} levels: {observed}")

    attrition = pd.DataFrame(
        [
            {
                "StageOrder": 1,
                "Stage": "Public baseline WGS profiles",
                "Samples": len(metadata),
                "Events": int(event.fillna(0).sum()),
            },
            {
                "StageOrder": 2,
                "Stage": "Complete positive PFS time and event indicator",
                "Samples": len(cohort),
                "Events": int(cohort["Event"].sum()),
            },
            {
                "StageOrder": 3,
                "Stage": "Paper-aligned clinical covariates complete",
                "Samples": int(cohort[primary_covariates].notna().all(axis=1).sum()),
                "Events": int(cohort.loc[cohort[primary_covariates].notna().all(axis=1), "Event"].sum()),
            },
            {
                "StageOrder": 4,
                "Stage": "Sequencing-quality sensitivity subset",
                "Samples": int(cohort["QualitySensitivityPass"].sum()),
                "Events": int(cohort.loc[cohort["QualitySensitivityPass"], "Event"].sum()),
            },
        ]
    )

    feature_audit = feature.loc[faec_rows].reset_index().rename(columns={"row.names": "LKTFeature"})
    feature_audit["SurvivalCohortPrevalence"] = [
        float(eligible_ppm.loc[row].gt(0).mean()) for row in faec_rows
    ]
    feature_audit["SurvivalCohortMedianPPM"] = [
        float(eligible_ppm.loc[row].median()) for row in faec_rows
    ]

    variable_rows = []
    for source, output_name in (
        ("pfs_d", "PFS_days"),
        ("pfsevent", "Event"),
        ("BMI", "BMI"),
        ("primary_cat", "PrimarySubtype"),
        ("adv_substage", "AdvancedSubstage"),
        ("LDH", "LDH"),
        ("GbNAHS", "GbNonHuman"),
        ("PctAss", "PercentAssembled"),
    ):
        values = metadata[source].replace({"N_A": np.nan, "": np.nan})
        variable_rows.append(
            {
                "SourceVariable": source,
                "AnalysisVariable": output_name,
                "PublicRows": len(metadata),
                "PublicNonMissing": int(values.notna().sum()),
                "SurvivalRows": len(cohort),
                "SurvivalNonMissing": int(cohort[output_name].notna().sum()),
            }
        )
    variable_audit = pd.DataFrame(variable_rows)

    composition_audit = cohort[
        [
            "SampleID",
            "ReportedLKTPPM",
            "PositiveLKTFeatures",
            "Richness50PPM",
            "ShannonReportedLKT",
            "InverseSimpsonReportedLKT",
            "FaecalibacteriumPPM",
            "GbNonHuman",
            "PercentAssembled",
            "QualitySensitivityPass",
        ]
    ].copy()

    matrix_out = eligible_ppm.reset_index().rename(columns={"row.names": "LKTFeature"})
    feature_out = feature.loc[eligible_ppm.index].reset_index().rename(columns={"row.names": "LKTFeature"})
    methods = {
        "article": 68,
        "analysis_seed": SEED,
        "plot_seed": PLOT_SEED,
        "pfs_definition": "Days from treatment start to progression/death or last vital assessment",
        "event_definition": "1=progression or death; 0=right-censored",
        "primary_microbiome_feature": "Sum of all filtered LKT rows assigned to g__Faecalibacterium",
        "feature_transform": f"log2(PPM + {PSEUDOCOUNT_PPM:g})",
        "primary_adjustment": ["primary melanoma subtype", "advanced substage", "LDH category", "BMI z-score"],
        "quality_sensitivity": f"Gb non-human >= {QUALITY_GB_GATE:g} and percent assembled >= {QUALITY_ASSEMBLED_GATE:g}",
        "independent_unit": "patient; one baseline stool WGS sample per patient",
        "cutoff_rule": "optimize only inside training data; apply unchanged to held-out test data",
        "cv_rule": "20 repeated event-stratified 5-fold splits; every prediction is out-of-fold",
    }
    metrics = {
        "article": 68,
        "public_wgs_samples": 167,
        "complete_pfs_samples": len(cohort),
        "pfs_events": int(cohort["Event"].sum()),
        "pfs_censored": int((1 - cohort["Event"]).sum()),
        "reported_lkt_features": len(eligible_ppm),
        "faecalibacterium_lkt_rows": len(faec_rows),
        "faecalibacterium_detected_samples": int(cohort["FaecalibacteriumPPM"].gt(0).sum()),
        "quality_sensitivity_samples": int(cohort["QualitySensitivityPass"].sum()),
        "quality_sensitivity_events": int(cohort.loc[cohort["QualitySensitivityPass"], "Event"].sum()),
        "analysis_seed": SEED,
        "plot_seed": PLOT_SEED,
        "pseudocount_ppm": PSEUDOCOUNT_PPM,
    }

    shutil.copy2(cache / "download-manifest.json", output / "source-manifest.json")
    shutil.copy2(cache / "spencer-data-dictionary.xlsx", output / "spencer-data-dictionary.xlsx")
    extract_anchor(cache / "spencer-paper.pdf", output / "spencer-fig3a-original.png")
    write_tsv(cohort, output / "survival-cohort.tsv")
    write_tsv(attrition, output / "cohort-attrition.tsv")
    write_tsv(feature_audit, output / "faecalibacterium-feature-audit.tsv")
    write_tsv(variable_audit, output / "metadata-variable-audit.tsv")
    write_tsv(composition_audit, output / "composition-audit.tsv")
    write_tsv(matrix_out, output / "lkt-survival-ppm.tsv.gz")
    write_tsv(feature_out, output / "lkt-feature-map.tsv")
    (output / "methods-contract.json").write_text(
        json.dumps(methods, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    versions = {
        "python": platform.python_version(),
        "numpy": safe_version("numpy"),
        "pandas": safe_version("pandas"),
        "openpyxl": safe_version("openpyxl"),
        "pillow": safe_version("Pillow"),
    }
    (output / "software-versions-python.json").write_text(
        json.dumps(versions, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
