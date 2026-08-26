#!/usr/bin/env python3
"""Create a checksum-covered Article 67 longitudinal evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd


WORK_FILES = (
    "analysis-metrics.json",
    "consecutive-intervals.tsv.gz",
    "dysbiosis-reference-audit.tsv",
    "dysbiosis-summary.tsv",
    "lag-summary.tsv",
    "lloyd-price-fig3-original.png",
    "methods-contract.json",
    "normalization-audit.tsv",
    "prevotella-selected-trajectories.tsv",
    "primary-dysbiosis-lmm.rds",
    "primary-fixed-effects.tsv",
    "primary-marginal-predictions.tsv",
    "primary-model-audit.tsv",
    "primary-model-diagnostics.tsv",
    "primary-random-effects.tsv",
    "primary-type3-anova.tsv",
    "profile-selection-ledger.tsv",
    "r-session-info.txt",
    "random-slope-model-comparison.tsv",
    "retention-subject-summary.tsv",
    "retention-summary.tsv",
    "sample-attrition.tsv",
    "sample-ledger.tsv",
    "selected-species-clr.tsv.gz",
    "shift-summary.tsv",
    "software-versions-python.json",
    "software-versions-r.tsv",
    "source-manifest.json",
    "species-feature-audit.tsv",
    "species-mixed-model-results.tsv",
    "species-model-audit.tsv",
    "species-relative.tsv.gz",
    "subject-shift-summary.tsv",
    "within-subject-pairs.tsv.gz",
)

SCRIPT_FILES = (
    "download_article67_longitudinal_data.py",
    "prepare_article67_longitudinal.py",
    "run_article67_longitudinal_models.R",
    "plot_article67_longitudinal.py",
    "freeze_article67_longitudinal.py",
    "validate_article67_longitudinal.py",
    "article42_44_validation_utils.py",
)

ENV_FILES = ("multiomics-python.yml", "multiomics-r-packages.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy(source: Path, target: Path) -> None:
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    target = args.output_dir.resolve()
    metrics = json.loads((work / "analysis-metrics.json").read_text(encoding="utf-8"))
    if metrics.get("article") != 67:
        raise RuntimeError("Article identity mismatch")

    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    for name in WORK_FILES:
        copy(work / name, target / name)
    for name in SCRIPT_FILES:
        copy(root / "scripts" / name, target / "scripts" / name)
    for name in ENV_FILES:
        copy(root / "env" / name, target / "env" / name)

    fixed = pd.read_csv(work / "primary-fixed-effects.tsv", sep="\t")
    lookup = fixed.set_index("Term")
    comparison = pd.read_csv(
        work / "random-slope-model-comparison.tsv", sep="\t"
    )
    contract = {
        "article": 67,
        "study": "Lloyd-Price et al. Nature 2019",
        "doi": "10.1038/s41586-019-1237-9",
        "pmcid": "PMC6650278",
        "cohort": "HMP2 / IBDMDB",
        "source_profiles_included": False,
        "source_metadata_included": False,
        "source_figure_included": True,
        "source_manifest_and_checksums_included": True,
        "eligible_samples": metrics["eligible_samples"],
        "eligible_subjects": metrics["eligible_subjects"],
        "selected_species": metrics["selected_species"],
        "primary_model": "REML LMM with subject random intercept and time slope; Kenward-Roger Type III tests",
        "cd_midpoint_estimate": float(lookup.loc["DiagnosisCD", "Estimate"]),
        "antibiotics_estimate": float(lookup.loc["AntibioticsYes", "Estimate"]),
        "cd_time_interaction_estimate": float(
            lookup.loc["DiagnosisCD:WeekYearCentered", "Estimate"]
        ),
        "random_slope_lrt_p": float(comparison.iloc[-1]["Pr(>Chisq)"]),
        "uncertainty_unit": "subject",
        "seed": metrics["seed"],
        "plot_seed": metrics["plot_seed"],
        "interpretation": "within-person association and detection stability; not proof of permanent colonization, treatment effect, or causality",
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    notice = """Article 67 frozen evidence bundle

The source study is Lloyd-Price et al., Nature 2019 (PMCID PMC6650278,
DOI 10.1038/s41586-019-1237-9). Exact official HMP2 profile, metadata and
publisher-figure URLs, byte counts and SHA-256 hashes are recorded in
source-manifest.json. The large source profile and metadata files are fetched
by the included downloader and are not duplicated here; Figure 3 is retained
as the visual anchor.

The bundle contains the deterministic profile-selection ledger, eligibility
attrition, renormalized terminal-species matrix, subject-aware Bray-Curtis and
retention summaries, primary random-slope mixed model, exploratory 63-species
mixed-model screen, full diagnostics and exact software versions. Bootstrap
resampling is performed at subject level. The fixed Bray-Curtis >0.54 event
definition is explicitly a sensitivity analysis, not a reconstruction of the
paper's lag-dependent shift algorithm. Detection at two visits does not prove
permanent colonization, and longitudinal association does not establish a
treatment effect or causality.
"""
    (target / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    payloads = sorted(
        path
        for path in target.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    (target / "file-checksums.sha256").write_text(
        "\n".join(
            f"{sha256(path)}  {path.relative_to(target).as_posix()}" for path in payloads
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Frozen Article 67 bundle: {len(payloads)} payload files in {target}")


if __name__ == "__main__":
    main()
