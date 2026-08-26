#!/usr/bin/env python3
"""Prepare a timing- and assumption-auditable mediation cohort for Article 69."""

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


SEED = 69_001
PLOT_SEED = 20_260_769
PSEUDOCOUNT_PPM = 25.0
FIBER_THRESHOLD = 20.0
QUALITY_GB_GATE = 0.5
QUALITY_ASSEMBLED_GATE = 50.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace({"N_A": np.nan, "": np.nan}), errors="coerce")


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="gzip" if path.suffix == ".gz" else None)


def safe_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def extract_anchor(paper: Path, target: Path) -> None:
    """Rasterize page 7 and crop the published Figure 3A/B patient panels."""
    with tempfile.TemporaryDirectory(prefix="article69-anchor-") as temporary:
        prefix = Path(temporary) / "page"
        subprocess.run(
            ["pdftoppm", "-f", "7", "-l", "7", "-singlefile", "-png", "-r", "240", str(paper), str(prefix)],
            check=True,
            capture_output=True,
        )
        page = Image.open(prefix.with_suffix(".png")).convert("RGB")
        width, height = page.size
        crop = page.crop((round(width * 0.055), round(height * 0.045), round(width * 0.950), round(height * 0.315)))
        crop.save(target, format="PNG", optimize=True)


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((cache / "download-manifest.json").read_text(encoding="utf-8"))
    if manifest.get("article") != 69:
        raise RuntimeError("Article identity mismatch in download manifest")

    workbook = cache / "spencer-human-wgs.xlsx"
    metadata = pd.read_excel(workbook, sheet_name="Metadata", dtype=str, keep_default_na=False)
    ppm = pd.read_excel(workbook, sheet_name="LKT_PPM").set_index("row.names").apply(pd.to_numeric, errors="raise")
    feature = pd.read_excel(workbook, sheet_name="LKT_featuretable").set_index("row.names").reindex(ppm.index)
    if len(metadata) != 167 or metadata["SRA_accession"].nunique() != 167 or ppm.shape != (225, 167):
        raise RuntimeError("Unexpected Spencer public WGS dimensions")
    if set(ppm.columns) != set(metadata["SRA_accession"]):
        raise RuntimeError("WGS profile columns do not match metadata")

    fiber = numeric(metadata["DSQfib"])
    category_known = metadata["fiber_cat"].isin(["Insufficient_intake", "Sufficient_intake"])
    response_known = metadata["response"].isin(["Responder", "Non_Responder"])
    eligible = metadata.loc[category_known & response_known & fiber.notna()].copy()
    eligible["FiberGrams"] = fiber.loc[eligible.index]
    eligible["ExposureSufficient"] = eligible["fiber_cat"].eq("Sufficient_intake").astype(int)
    eligible["OutcomeResponder"] = eligible["response"].eq("Responder").astype(int)
    if len(eligible) != 94:
        raise RuntimeError(f"Expected 94 WGS patients with fiber and response, observed {len(eligible)}")
    if not (
        eligible.loc[eligible["ExposureSufficient"].eq(1), "FiberGrams"].ge(FIBER_THRESHOLD).all()
        and eligible.loc[eligible["ExposureSufficient"].eq(0), "FiberGrams"].lt(FIBER_THRESHOLD).all()
    ):
        raise RuntimeError("Published fiber category does not match the locked 20 g/day threshold")

    sample_ids = eligible["SRA_accession"].tolist()
    eligible_ppm = ppm.loc[:, sample_ids]
    faec_rows = feature.index[feature["Genus"].eq("g__Faecalibacterium")]
    rumi_rows = feature.index[feature["Family"].eq("f__Ruminococcaceae")]
    if len(faec_rows) != 3 or len(rumi_rows) != 26:
        raise RuntimeError(f"Unexpected prespecified feature rows: Faec={len(faec_rows)}, Ruminococcaceae={len(rumi_rows)}")
    faec = eligible_ppm.loc[faec_rows].sum(axis=0)
    rumi = eligible_ppm.loc[rumi_rows].sum(axis=0)

    index = eligible.set_index("SRA_accession")
    cohort = pd.DataFrame(index=sample_ids)
    cohort.index.name = "SampleID"
    cohort["ExposureSufficient"] = index.loc[sample_ids, "ExposureSufficient"]
    cohort["FiberCategory"] = index.loc[sample_ids, "fiber_cat"]
    cohort["FiberGrams"] = index.loc[sample_ids, "FiberGrams"]
    cohort["OutcomeResponder"] = index.loc[sample_ids, "OutcomeResponder"]
    cohort["Response"] = index.loc[sample_ids, "response"]
    cohort["FaecalibacteriumPPM"] = faec.loc[sample_ids]
    cohort["FaecalibacteriumLog2"] = np.log2(cohort["FaecalibacteriumPPM"] + PSEUDOCOUNT_PPM)
    cohort["RuminococcaceaePPM"] = rumi.loc[sample_ids]
    cohort["RuminococcaceaeLog2"] = np.log2(cohort["RuminococcaceaePPM"] + PSEUDOCOUNT_PPM)
    cohort["BMI"] = numeric(index.loc[sample_ids, "BMI"])
    cohort["PrimarySubtype"] = index.loc[sample_ids, "primary_cat"]
    cohort["AdvancedSubstage"] = index.loc[sample_ids, "adv_substage"]
    cohort["LDH"] = index.loc[sample_ids, "LDH"]
    cohort["Treatment"] = index.loc[sample_ids, "treatment"]
    cohort["TreatmentNaive"] = index.loc[sample_ids, "Treatment_naive"]
    cohort["Antibiotics"] = index.loc[sample_ids, "antibiotics"]
    cohort["Probiotics"] = index.loc[sample_ids, "probiotics"]
    cohort["GbNonHuman"] = numeric(index.loc[sample_ids, "GbNAHS"])
    cohort["PercentAssembled"] = numeric(index.loc[sample_ids, "PctAss"])
    cohort["QualitySensitivityPass"] = cohort["GbNonHuman"].ge(QUALITY_GB_GATE) & cohort["PercentAssembled"].ge(QUALITY_ASSEMBLED_GATE)
    cohort = cohort.reset_index()

    covariates = ["BMI", "PrimarySubtype", "AdvancedSubstage", "LDH"]
    if cohort[covariates].isna().any().any():
        raise RuntimeError("Paper-aligned confounder set must be complete")
    if cohort["ExposureSufficient"].value_counts().to_dict() != {0: 71, 1: 23}:
        raise RuntimeError("Unexpected fiber exposure counts")
    if cohort["OutcomeResponder"].value_counts().to_dict() != {1: 60, 0: 34}:
        raise RuntimeError("Unexpected response counts")

    attrition = pd.DataFrame(
        [
            {"StageOrder": 1, "Stage": "Public baseline WGS profiles", "Samples": 167, "Responders": int(metadata["response"].eq("Responder").sum())},
            {"StageOrder": 2, "Stage": "Response status available", "Samples": int(response_known.sum()), "Responders": int(metadata.loc[response_known, "response"].eq("Responder").sum())},
            {"StageOrder": 3, "Stage": "Published fiber category and grams available", "Samples": len(cohort), "Responders": int(cohort["OutcomeResponder"].sum())},
            {"StageOrder": 4, "Stage": "Paper-aligned confounders complete", "Samples": int(cohort[covariates].notna().all(axis=1).sum()), "Responders": int(cohort.loc[cohort[covariates].notna().all(axis=1), "OutcomeResponder"].sum())},
            {"StageOrder": 5, "Stage": "Sequencing-quality sensitivity subset", "Samples": int(cohort["QualitySensitivityPass"].sum()), "Responders": int(cohort.loc[cohort["QualitySensitivityPass"], "OutcomeResponder"].sum())},
        ]
    )
    group_summary = (
        cohort.groupby(["FiberCategory", "ExposureSufficient"], as_index=False)
        .agg(
            Patients=("SampleID", "size"),
            Responders=("OutcomeResponder", "sum"),
            ResponseRate=("OutcomeResponder", "mean"),
            MedianFiberGrams=("FiberGrams", "median"),
            MedianFaecalibacteriumPPM=("FaecalibacteriumPPM", "median"),
        )
    )
    feature_audit = pd.concat(
        [
            feature.loc[faec_rows].assign(MediatorDefinition="Primary: genus Faecalibacterium"),
            feature.loc[rumi_rows].assign(MediatorDefinition="Sensitivity: family Ruminococcaceae"),
        ]
    ).reset_index().rename(columns={"row.names": "LKTFeature"})
    feature_audit["CohortPrevalence"] = [float(eligible_ppm.loc[row].gt(0).mean()) for row in feature_audit["LKTFeature"]]

    methods = {
        "article": 69,
        "analysis_seed": SEED,
        "plot_seed": PLOT_SEED,
        "exposure": "habitual dietary fiber sufficient (>=20 g/day) versus insufficient (<20 g/day)",
        "mediator": "log2(sum of three g__Faecalibacterium LKT rows in PPM + 25)",
        "outcome": "binary ICB responder by the study RECIST-based definition",
        "adjustment": ["primary melanoma subtype", "advanced substage", "LDH category", "BMI z-score"],
        "primary_estimand": "standardized interventional-analog risk-difference decomposition under parametric mediator and outcome models",
        "timing_limit": "habitual-fiber questionnaire and baseline stool were contemporaneous; exposure-before-mediator ordering is not experimentally established",
        "independent_unit": "patient; one baseline stool WGS sample per patient",
    }
    metrics = {
        "article": 69,
        "public_wgs_samples": 167,
        "mediation_samples": len(cohort),
        "sufficient_fiber": int(cohort["ExposureSufficient"].sum()),
        "insufficient_fiber": int((1 - cohort["ExposureSufficient"]).sum()),
        "responders": int(cohort["OutcomeResponder"].sum()),
        "nonresponders": int((1 - cohort["OutcomeResponder"]).sum()),
        "faecalibacterium_lkt_rows": len(faec_rows),
        "ruminococcaceae_lkt_rows": len(rumi_rows),
        "faecalibacterium_detected": int(cohort["FaecalibacteriumPPM"].gt(0).sum()),
        "quality_sensitivity_samples": int(cohort["QualitySensitivityPass"].sum()),
        "analysis_seed": SEED,
        "plot_seed": PLOT_SEED,
        "fiber_threshold_g_day": FIBER_THRESHOLD,
        "pseudocount_ppm": PSEUDOCOUNT_PPM,
    }

    shutil.copy2(cache / "download-manifest.json", output / "source-manifest.json")
    shutil.copy2(cache / "spencer-data-dictionary.xlsx", output / "spencer-data-dictionary.xlsx")
    extract_anchor(cache / "spencer-paper.pdf", output / "spencer-fig3ab-original.png")
    write_tsv(cohort, output / "mediation-cohort.tsv")
    write_tsv(attrition, output / "cohort-attrition.tsv")
    write_tsv(group_summary, output / "exposure-outcome-summary.tsv")
    write_tsv(feature_audit, output / "mediator-feature-audit.tsv")
    write_tsv(eligible_ppm.reset_index().rename(columns={"row.names": "LKTFeature"}), output / "lkt-mediation-ppm.tsv.gz")
    write_tsv(feature.loc[eligible_ppm.index].reset_index().rename(columns={"row.names": "LKTFeature"}), output / "lkt-feature-map.tsv")
    (output / "methods-contract.json").write_text(json.dumps(methods, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "analysis-metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    versions = {
        "python": platform.python_version(),
        "numpy": safe_version("numpy"),
        "pandas": safe_version("pandas"),
        "openpyxl": safe_version("openpyxl"),
        "pillow": safe_version("Pillow"),
    }
    (output / "software-versions-python.json").write_text(json.dumps(versions, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
