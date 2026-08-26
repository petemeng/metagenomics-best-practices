#!/usr/bin/env python3
"""Validate Article 20 functional-profile normalization and interpretation."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article20-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


SPECIAL_FEATURES = ("UNMAPPED", "UNINTEGRATED", "UNGROUPED")
FEATURE_SPACES = {
    "Gene families": "genefamilies-rpk.tsv",
    "Regrouped reactions": "reactions-rpk.tsv",
    "Pathway abundance": "pathabundance-rpk.tsv",
}
FIGURE_STEMS = (
    "20-normalization-denominators",
    "20-special-feature-budget",
    "20-pathway-contributions",
    "20-prevalence-zero-sensitivity",
)
EXPECTED_RAW = {
    "pathway_abundance": {
        "experimenthub_id": "EH7089",
        "bytes": 367_712,
        "sha256": "ead7c78c075fec92a7d641b731594e068b2ba2a47479151d081c338f615af121",
    },
    "pathway_coverage": {
        "experimenthub_id": "EH7090",
        "bytes": 262_967,
        "sha256": "73a1b77b70f88e9028e8707ba3e99b93f0ff99cd91401a3966c4f7a31dbfc3a1",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--environment-prefix", required=True, type=Path)
    parser.add_argument("--article19-dir", required=True, type=Path)
    parser.add_argument("--cohort-dir", required=True, type=Path)
    parser.add_argument("--frozen-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--figure-dir", required=True, type=Path)
    return parser.parse_args()


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict[str, str]] = []

    def add(self, category: str, check_id: str, passed: bool, detail: Any) -> None:
        self.rows.append(
            {
                "Category": category,
                "CheckID": check_id,
                "Status": "PASS" if passed else "FAIL",
                "Detail": str(detail),
            }
        )

    @property
    def passed(self) -> int:
        return sum(row["Status"] == "PASS" for row in self.rows)

    @property
    def failed(self) -> int:
        return sum(row["Status"] == "FAIL" for row in self.rows)


def feature_base(feature: str) -> str:
    return feature.split("|", 1)[0]


def feature_id(feature: str) -> str:
    return feature_base(feature).split(": ", 1)[0]


def stratum(feature: str) -> str:
    return feature.split("|", 1)[1] if "|" in feature else ""


def is_special(feature: str) -> bool:
    return feature_id(feature) in SPECIAL_FEATURES


def read_humann_table(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("#"):
        raise ValueError(f"Unexpected HUMAnN header in {path}")
    header = lines[0].split("\t")
    if len(header) != 2:
        raise ValueError(f"Article 20 expects a one-sample HUMAnN table: {path}")
    rows: list[dict[str, Any]] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError(f"Malformed HUMAnN row: {line[:100]}")
        rows.append(
            {
                "Feature": fields[0],
                "BaseFeature": feature_base(fields[0]),
                "FeatureID": feature_id(fields[0]),
                "Stratum": stratum(fields[0]),
                "Level": len(fields[0].split("|")),
                "Special": is_special(fields[0]),
                "Value": float(fields[1]),
            }
        )
    return {"header": header, "rows": rows}


def verify_checksum_manifest(
    directory: Path, checks: Checks, category: str
) -> int:
    manifest = directory / "file-checksums.sha256"
    entries: dict[str, str] = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        relative = relative.strip()
        path = directory / relative
        observed = hash_file(path) if path.is_file() else "missing"
        checks.add(category, f"sha256-{relative}", observed == expected, observed)
        entries[relative] = expected
    payloads = {
        path.relative_to(directory).as_posix()
        for path in directory.rglob("*")
        if path.is_file() and path != manifest
    }
    checks.add(
        category,
        "manifest-complete",
        payloads == set(entries),
        f"payloads={len(payloads)}; entries={len(entries)}",
    )
    return len(entries)


def normalize_table(
    rows: list[dict[str, Any]], unit: str, mode: str, special: str
) -> tuple[list[dict[str, Any]], dict[int, float]]:
    kept = [row for row in rows if special == "y" or not row["Special"]]
    totals: defaultdict[int, float] = defaultdict(float)
    for row in kept:
        totals[row["Level"]] += row["Value"]
    safe_totals = {level: total if total != 0 else 1.0 for level, total in totals.items()}
    target = 1_000_000.0 if unit == "cpm" else 1.0
    output: list[dict[str, Any]] = []
    for row in kept:
        denominator = safe_totals[row["Level"]] if mode == "levelwise" else safe_totals[1]
        normalized = row["Value"] / denominator * target
        output.append({**row, "Normalized": normalized})
    return output, safe_totals


def build_normalization_audits(
    tables: dict[str, dict[str, Any]], checks: Checks
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], int]:
    feature_rows: list[dict[str, Any]] = []
    branch_rows: list[dict[str, Any]] = []
    special_rows: list[dict[str, Any]] = []
    closure_failures = 0
    for space, table in tables.items():
        rows = table["rows"]
        level1 = [row for row in rows if row["Level"] == 1]
        ordinary = [row for row in level1 if not row["Special"]]
        special = [row for row in level1 if row["Special"]]
        all_total = math.fsum(row["Value"] for row in level1)
        ordinary_total = math.fsum(row["Value"] for row in ordinary)
        strata_total = math.fsum(row["Value"] for row in rows if row["Level"] == 2)
        feature_rows.append(
            {
                "FeatureSpace": space,
                "NativeUnit": "RPK" if space != "Pathway abundance" else "RPK-derived pathway abundance",
                "Rows": len(rows),
                "CommunityRows": len(level1),
                "StratifiedRows": len(rows) - len(level1),
                "OrdinaryCommunityFeatures": len(ordinary),
                "SpecialCommunityRows": len(special),
                "CommunityTotalWithSpecial": all_total,
                "CommunityTotalWithoutSpecial": ordinary_total,
                "StrataTotalWithSpecial": strata_total,
                "CanAddStrataToCommunity": "No" if space == "Pathway abundance" else "Yes",
            }
        )
        values = defaultdict(float)
        values["Ordinary"] = ordinary_total
        for row in special:
            values[row["FeatureID"]] += row["Value"]
        for category in ("Ordinary", *SPECIAL_FEATURES):
            value = values[category]
            special_rows.append(
                {
                    "FeatureSpace": space,
                    "Category": category,
                    "NativeValue": value,
                    "FractionOfCommunityTotal": value / all_total if all_total else 0.0,
                    "IncludedWhenSpecialN": "Yes" if category == "Ordinary" else "No",
                }
            )
        for unit in ("cpm", "relab"):
            for mode in ("community", "levelwise"):
                for special_choice in ("y", "n"):
                    normalized, totals = normalize_table(rows, unit, mode, special_choice)
                    target = 1_000_000.0 if unit == "cpm" else 1.0
                    sums = {
                        level: math.fsum(
                            row["Normalized"] for row in normalized if row["Level"] == level
                        )
                        for level in sorted(totals)
                    }
                    levels_to_close = sorted(totals) if mode == "levelwise" else [1]
                    failures = sum(
                        not math.isclose(sums[level], target, rel_tol=1e-10, abs_tol=1e-7)
                        for level in levels_to_close
                    )
                    closure_failures += failures
                    branch_rows.append(
                        {
                            "FeatureSpace": space,
                            "Unit": unit,
                            "Mode": mode,
                            "Special": special_choice,
                            "InputRows": len(rows),
                            "OutputRows": len(normalized),
                            "CommunityDenominator": totals.get(1, 0.0),
                            "StrataDenominator": totals.get(2, 0.0),
                            "CommunityOutputSum": sums.get(1, 0.0),
                            "StrataOutputSum": sums.get(2, 0.0),
                            "Target": target,
                            "RequiredClosureLevels": ";".join(map(str, levels_to_close)),
                            "ClosureFailures": failures,
                            "ActualToolVerified": "Pending" if space == "Pathway abundance" and unit == "relab" else "Not selected",
                        }
                    )
    checks.add("normalization", "three-feature-spaces", len(tables) == 3, len(tables))
    checks.add("normalization", "twenty-four-branches", len(branch_rows) == 24, len(branch_rows))
    checks.add("normalization", "closure-failures", closure_failures == 0, closure_failures)
    return feature_rows, branch_rows, special_rows, closure_failures


def compare_actual_outputs(
    path_rows: list[dict[str, Any]], frozen_dir: Path, branch_rows: list[dict[str, Any]], checks: Checks
) -> tuple[float, int]:
    maximum_relative_error = 0.0
    verified = 0
    for mode in ("community", "levelwise"):
        for special in ("y", "n"):
            expected, _ = normalize_table(path_rows, "relab", mode, special)
            actual = read_humann_table(
                frozen_dir / f"pathabundance-relab-{mode}-special-{special}.tsv"
            )["rows"]
            same_features = [row["Feature"] for row in expected] == [row["Feature"] for row in actual]
            errors: list[float] = []
            for exp, obs in zip(expected, actual):
                denominator = max(abs(exp["Normalized"]), 1e-15)
                errors.append(abs(exp["Normalized"] - obs["Value"]) / denominator)
            max_error = max(errors, default=0.0)
            maximum_relative_error = max(maximum_relative_error, max_error)
            passed = same_features and len(expected) == len(actual) and max_error <= 5e-6
            checks.add(
                "actual-tool",
                f"path-relab-{mode}-special-{special}",
                passed,
                f"rows={len(actual)};max_relative_error={max_error:.9g}",
            )
            if passed:
                verified += 1
            for row in branch_rows:
                if (
                    row["FeatureSpace"] == "Pathway abundance"
                    and row["Unit"] == "relab"
                    and row["Mode"] == mode
                    and row["Special"] == special
                ):
                    row["ActualToolVerified"] = "Yes" if passed else "No"
    return maximum_relative_error, verified


def read_matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        features: list[str] = []
        values: list[list[float]] = []
        for fields in reader:
            features.append(fields[0])
            values.append([float(value) for value in fields[1:]])
    return features, header[1:], np.asarray(values, dtype=float)


def cohort_audits(
    cohort_dir: Path, checks: Checks
) -> tuple[
    list[dict[str, Any]], list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]
]:
    manifest = read_tsv(cohort_dir / "resource-manifest.tsv")
    manifest_by_resource = {row["resource"]: row for row in manifest}
    lineage: list[dict[str, Any]] = []
    for resource, expected in EXPECTED_RAW.items():
        row = manifest_by_resource[resource]
        raw_path = cohort_dir / row["raw_file"]
        checks.add("cohort-source", f"{resource}-eh-id", row["experimenthub_id"] == expected["experimenthub_id"], row["experimenthub_id"])
        checks.add("cohort-source", f"{resource}-bytes", raw_path.stat().st_size == expected["bytes"] == int(row["raw_bytes"]), raw_path.stat().st_size)
        checks.add("cohort-source", f"{resource}-sha256", hash_file(raw_path) == expected["sha256"] == row["raw_sha256"], hash_file(raw_path))
        checks.add("cohort-source", f"{resource}-package", row["package_version"] == "3.12.0", row["package_version"])
        lineage.append(
            {
                "Source": "curatedMetagenomicData",
                "Resource": resource,
                "Identifier": row["experimenthub_id"],
                "Release": "2021-10-14 / package 3.12.0",
                "Rows": int(row["rows"]),
                "Samples": int(row["samples"]),
                "SHA256": row["raw_sha256"],
                "Use": "Prevalence and zero-value method audit",
            }
        )

    abundance_features, abundance_samples, abundance = read_matrix(cohort_dir / "pathway-abundance.tsv.gz")
    coverage_features, coverage_samples, coverage = read_matrix(cohort_dir / "pathway-coverage.tsv.gz")
    metadata = read_tsv(cohort_dir / "sample-metadata.tsv")
    metadata_samples = [row["sample_id"] for row in metadata]
    subjects = [row["subject_id"] for row in metadata]
    checks.add("cohort", "abundance-dimensions", abundance.shape == (11173, 24), abundance.shape)
    checks.add("cohort", "coverage-dimensions", coverage.shape == (11173, 24), coverage.shape)
    checks.add("cohort", "feature-order-aligned", abundance_features == coverage_features, len(abundance_features))
    checks.add("cohort", "sample-order-aligned", abundance_samples == coverage_samples == metadata_samples, len(metadata_samples))
    checks.add("cohort", "fifteen-subjects", len(set(subjects)) == 15, len(set(subjects)))
    checks.add("cohort", "no-missing-values", not np.isnan(abundance).any() and not np.isnan(coverage).any(), "abundance and coverage")
    out_of_range = int(np.sum((coverage < 0) | (coverage > 1)))
    checks.add("coverage", "zero-to-one", out_of_range == 0, out_of_range)
    checks.add("coverage", "not-renormalized", not np.allclose(coverage.sum(axis=0), 1.0), tuple(np.round([coverage.sum(axis=0).min(), coverage.sum(axis=0).max()], 6)))

    keep = np.asarray(
        ["|" not in feature and not is_special(feature) for feature in abundance_features],
        dtype=bool,
    )
    ordinary_features = [feature for feature, flag in zip(abundance_features, keep) if flag]
    ordinary_abundance = abundance[keep]
    ordinary_coverage = coverage[keep]
    subject_order = list(dict.fromkeys(subjects))
    subject_indices = [np.asarray([subject == item for subject in subjects]) for item in subject_order]
    detection_rules = {
        "Abundance > 0": ordinary_abundance > 0,
        "Coverage > 0": ordinary_coverage > 0,
        "Coverage >= 0.5": ordinary_coverage >= 0.5,
    }
    prevalence_rows: list[dict[str, Any]] = []
    for rule, detected in detection_rules.items():
        subject_detected = np.column_stack(
            [np.any(detected[:, index], axis=1) for index in subject_indices]
        )
        for unit, matrix, denominator in (
            ("Profile", detected, len(abundance_samples)),
            ("Subject (any visit)", subject_detected, len(subject_order)),
        ):
            prevalence = matrix.mean(axis=1)
            for threshold in (0.10, 0.25, 0.50):
                retained = int(np.sum(prevalence >= threshold))
                prevalence_rows.append(
                    {
                        "AnalysisUnit": unit,
                        "DetectionRule": rule,
                        "PrevalenceThreshold": threshold,
                        "Units": denominator,
                        "OrdinaryUnstratifiedPathways": len(ordinary_features),
                        "RetainedPathways": retained,
                        "RemovedPathways": len(ordinary_features) - retained,
                    }
                )

    presence = np.mean(ordinary_abundance > 0, axis=1)
    positive_median = np.asarray([
        np.median(row[row > 0]) if np.any(row > 0) else 0.0
        for row in ordinary_abundance
    ])
    denominator_candidates = np.flatnonzero(presence == np.max(presence))
    denominator_index = int(denominator_candidates[np.argmax(positive_median[denominator_candidates])])
    numerator_candidates = np.flatnonzero(
        (presence >= 0.25)
        & (presence <= 0.50)
        & (np.arange(len(presence)) != denominator_index)
    )
    numerator_score = ordinary_abundance[numerator_candidates].mean(axis=1)
    numerator_index = int(numerator_candidates[np.argmax(numerator_score)])
    numerator = ordinary_abundance[numerator_index]
    denominator = ordinary_abundance[denominator_index]
    pseudocount_rows: list[dict[str, Any]] = []
    for pseudocount in (1e-6, 1e-5, 1e-4):
        log_ratio = np.log2((numerator + pseudocount) / (denominator + pseudocount))
        pseudocount_rows.append(
            {
                "NumeratorPathway": ordinary_features[numerator_index],
                "DenominatorPathway": ordinary_features[denominator_index],
                "NumeratorZeroProfiles": int(np.sum(numerator == 0)),
                "DenominatorZeroProfiles": int(np.sum(denominator == 0)),
                "Pseudocount": pseudocount,
                "MinimumLog2Ratio": float(np.min(log_ratio)),
                "FirstQuartileLog2Ratio": float(np.quantile(log_ratio, 0.25)),
                "MedianLog2Ratio": float(np.median(log_ratio)),
                "ThirdQuartileLog2Ratio": float(np.quantile(log_ratio, 0.75)),
                "MaximumLog2Ratio": float(np.max(log_ratio)),
                "RangeLog2Ratio": float(np.ptp(log_ratio)),
            }
        )
    checks.add("prevalence", "ordinary-unstratified-nonempty", len(ordinary_features) > 100, len(ordinary_features))
    checks.add("prevalence", "eighteen-sensitivity-rows", len(prevalence_rows) == 18, len(prevalence_rows))
    checks.add("zero", "numerator-has-zero-and-positive", 0 < int(np.sum(numerator == 0)) < 24, int(np.sum(numerator == 0)))
    checks.add(
        "zero",
        "denominator-more-prevalent-than-numerator",
        int(np.sum(denominator == 0)) < int(np.sum(numerator == 0)),
        f"denominator_zeros={int(np.sum(denominator == 0))};numerator_zeros={int(np.sum(numerator == 0))}",
    )
    checks.add("zero", "three-pseudocounts", len(pseudocount_rows) == 3, [row["Pseudocount"] for row in pseudocount_rows])
    cohort_summary = {
        "rows": len(abundance_features),
        "samples": len(abundance_samples),
        "subjects": len(set(subjects)),
        "ordinary_unstratified_pathways": len(ordinary_features),
        "coverage_out_of_range": out_of_range,
        "numerator_pathway": ordinary_features[numerator_index],
        "denominator_pathway": ordinary_features[denominator_index],
    }
    plot_context = [
        {
            "features": ordinary_features,
            "abundance": ordinary_abundance,
            "coverage": ordinary_coverage,
        }
    ]
    return prevalence_rows, pseudocount_rows, cohort_summary, lineage + plot_context


def pathway_contribution_audit(
    path_rows: list[dict[str, Any]], checks: Checks
) -> list[dict[str, Any]]:
    community = {
        row["BaseFeature"]: row
        for row in path_rows
        if row["Level"] == 1 and not row["Special"]
    }
    strata: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in path_rows:
        if row["Level"] == 2 and not row["Special"]:
            strata[row["BaseFeature"]].append(row)
    selected = sorted(
        (row for base, row in community.items() if math.fsum(item["Value"] for item in strata[base]) > 0),
        key=lambda row: row["Value"],
        reverse=True,
    )[:6]
    output: list[dict[str, Any]] = []
    for row in selected:
        members = strata[row["BaseFeature"]]
        strata_sum = math.fsum(item["Value"] for item in members)
        dominant = max(members, key=lambda item: item["Value"])
        unclassified = math.fsum(item["Value"] for item in members if item["Stratum"] == "unclassified")
        dominant_value = dominant["Value"] if dominant["Stratum"] != "unclassified" else 0.0
        other_classified = strata_sum - unclassified - dominant_value
        for category, value in (
            ("Dominant taxon", dominant_value),
            ("Other classified", other_classified),
            ("Unclassified", unclassified),
        ):
            output.append(
                {
                    "Pathway": row["BaseFeature"],
                    "PathwayID": row["FeatureID"],
                    "CommunityAbundance": row["Value"],
                    "StrataAbundanceSum": strata_sum,
                    "StrataToCommunityRatio": strata_sum / row["Value"] if row["Value"] else 0.0,
                    "Category": category,
                    "Taxon": dominant["Stratum"] if category == "Dominant taxon" else category,
                    "StratifiedAbundance": value,
                    "ShareWithinStratifiedEvidence": value / strata_sum if strata_sum else 0.0,
                }
            )
    per_path = defaultdict(float)
    for row in output:
        per_path[row["Pathway"]] += row["ShareWithinStratifiedEvidence"]
    maximum_error = max((abs(value - 1.0) for value in per_path.values()), default=0.0)
    checks.add("contribution", "six-pathways", len(per_path) == 6, len(per_path))
    checks.add("contribution", "within-strata-closure", maximum_error < 1e-10, maximum_error)
    checks.add("contribution", "separate-community-ratio", all("StrataToCommunityRatio" in row for row in output), len(output))
    return output


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / f"{stem}.png", dpi=350, bbox_inches="tight")
    fig.savefig(
        figure_dir / f"{stem}.tiff",
        dpi=350,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def render_normalization_denominators(
    branches: list[dict[str, Any]], figure_dir: Path
) -> None:
    configure_plot_style()
    spaces = list(FEATURE_SPACES)
    colors = {"y": "#457B9D", "n": "#E76F51"}
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.9))
    width = 0.34
    x = np.arange(len(spaces))
    for offset, special in ((-width / 2, "y"), (width / 2, "n")):
        selected = [
            next(row for row in branches if row["FeatureSpace"] == space and row["Unit"] == "relab" and row["Mode"] == "community" and row["Special"] == special)
            for space in spaces
        ]
        axes[0].bar(
            x + offset,
            [row["CommunityDenominator"] for row in selected],
            width,
            label=f"Special {special}",
            color=colors[special],
        )
    axes[0].set_yscale("log")
    axes[0].set_xticks(x, ["Gene families", "Reactions", "Pathways"])
    axes[0].set_ylabel("Community denominator (native scale)")
    axes[0].set_title("A  Removing special rows changes the denominator", loc="left", fontweight="bold")
    axes[0].legend()
    axes[0].grid(axis="y", color="#DDDDDD", linewidth=0.6)

    group_x = np.arange(len(spaces) * 2)
    labels: list[str] = []
    community_values: list[float] = []
    levelwise_values: list[float] = []
    for space in spaces:
        for special in ("y", "n"):
            row = next(item for item in branches if item["FeatureSpace"] == space and item["Unit"] == "relab" and item["Mode"] == "community" and item["Special"] == special)
            labels.append(("Gene" if space == "Gene families" else "Reaction" if space == "Regrouped reactions" else "Pathway") + f"\nspecial {special}")
            community_values.append(row["StrataOutputSum"])
            levelwise_values.append(1.0)
    axes[1].bar(group_x - width / 2, community_values, width, color="#2A9D8F", label="Community mode")
    axes[1].bar(group_x + width / 2, levelwise_values, width, color="#F4A261", label="Levelwise mode")
    axes[1].axhline(1.0, color="#555555", linestyle="--", linewidth=0.8)
    axes[1].set_xticks(group_x, labels)
    axes[1].set_ylabel("Sum of normalized stratified rows")
    axes[1].set_title("B  Only levelwise mode closes each level", loc="left", fontweight="bold")
    axes[1].legend()
    axes[1].grid(axis="y", color="#DDDDDD", linewidth=0.6)
    fig.tight_layout()
    save_figure(fig, figure_dir, "20-normalization-denominators")


def render_special_budget(rows: list[dict[str, Any]], figure_dir: Path) -> None:
    configure_plot_style()
    spaces = list(FEATURE_SPACES)
    categories = ["Ordinary", "UNMAPPED", "UNINTEGRATED", "UNGROUPED"]
    colors = ["#2A9D8F", "#B7B7A4", "#E9C46A", "#E76F51"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    bottom = np.zeros(len(spaces))
    for category, color in zip(categories, colors):
        values = [
            next(row["FractionOfCommunityTotal"] for row in rows if row["FeatureSpace"] == space and row["Category"] == category)
            for space in spaces
        ]
        ax.bar(spaces, values, bottom=bottom, color=color, label=category, width=0.65)
        bottom += np.asarray(values)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of community denominator")
    ax.set_title("Special features are part of the normalization budget", loc="left", fontweight="bold")
    ax.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    fig.tight_layout()
    save_figure(fig, figure_dir, "20-special-feature-budget")


def render_pathway_contributions(rows: list[dict[str, Any]], figure_dir: Path) -> None:
    configure_plot_style()
    pathways = list(dict.fromkeys(row["PathwayID"] for row in rows))
    categories = ["Dominant taxon", "Other classified", "Unclassified"]
    colors = ["#457B9D", "#2A9D8F", "#B7B7A4"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), gridspec_kw={"width_ratios": [1.7, 1]})
    y = np.arange(len(pathways))
    left = np.zeros(len(pathways))
    for category, color in zip(categories, colors):
        values = [
            next(row["ShareWithinStratifiedEvidence"] for row in rows if row["PathwayID"] == pathway and row["Category"] == category)
            for pathway in pathways
        ]
        axes[0].barh(y, values, left=left, color=color, label=category, height=0.62)
        left += np.asarray(values)
    axes[0].set_yticks(y, pathways)
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("Share within stratified evidence")
    axes[0].set_title("A  Taxon contributions use the strata denominator", loc="left", fontweight="bold")
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3)
    axes[0].grid(axis="x", color="#DDDDDD", linewidth=0.6)
    ratios = [
        next(row["StrataToCommunityRatio"] for row in rows if row["PathwayID"] == pathway)
        for pathway in pathways
    ]
    axes[1].barh(y, ratios, color="#E76F51", height=0.62)
    axes[1].axvline(1.0, color="#555555", linestyle="--", linewidth=0.8)
    axes[1].set_yticks(y, [""] * len(pathways))
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Strata sum / community abundance")
    axes[1].set_title("B  Community abundance is separate", loc="left", fontweight="bold")
    axes[1].grid(axis="x", color="#DDDDDD", linewidth=0.6)
    fig.tight_layout()
    save_figure(fig, figure_dir, "20-pathway-contributions")


def render_prevalence_zero(
    prevalence: list[dict[str, Any]], pseudocounts: list[dict[str, Any]], figure_dir: Path
) -> None:
    configure_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    profile_rows = [row for row in prevalence if row["AnalysisUnit"] == "Profile"]
    colors = {"Abundance > 0": "#457B9D", "Coverage > 0": "#2A9D8F", "Coverage >= 0.5": "#E76F51"}
    for rule, color in colors.items():
        selected = sorted(
            (row for row in profile_rows if row["DetectionRule"] == rule),
            key=lambda row: row["PrevalenceThreshold"],
        )
        axes[0].plot(
            [row["PrevalenceThreshold"] * 100 for row in selected],
            [row["RetainedPathways"] for row in selected],
            marker="o",
            linewidth=2,
            color=color,
            label=rule,
        )
    axes[0].set_xlabel("Profile prevalence threshold (%)")
    axes[0].set_ylabel("Retained pathways")
    axes[0].set_title("A  Filtering depends on the detection rule", loc="left", fontweight="bold")
    axes[0].legend()
    axes[0].grid(color="#DDDDDD", linewidth=0.6)

    x = np.arange(len(pseudocounts))
    medians = [row["MedianLog2Ratio"] for row in pseudocounts]
    low = [row["FirstQuartileLog2Ratio"] for row in pseudocounts]
    high = [row["ThirdQuartileLog2Ratio"] for row in pseudocounts]
    axes[1].errorbar(
        x,
        medians,
        yerr=[np.asarray(medians) - np.asarray(low), np.asarray(high) - np.asarray(medians)],
        fmt="o",
        color="#6A4C93",
        capsize=5,
        linewidth=1.5,
    )
    axes[1].set_xticks(x, [f"{row['Pseudocount']:.0e}" for row in pseudocounts])
    axes[1].set_xlabel("Pseudocount")
    axes[1].set_ylabel("Log2 pathway ratio (median and IQR)")
    axes[1].set_title("B  Zero replacement changes log-ratios", loc="left", fontweight="bold")
    axes[1].grid(axis="y", color="#DDDDDD", linewidth=0.6)
    fig.tight_layout()
    save_figure(fig, figure_dir, "20-prevalence-zero-sensitivity")


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    prefix = args.environment_prefix.resolve()
    article19_dir = args.article19_dir.resolve()
    cohort_dir = args.cohort_dir.resolve()
    frozen_dir = args.frozen_dir.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    article19_entries = verify_checksum_manifest(article19_dir, checks, "article19-frozen")
    cohort_entries = verify_checksum_manifest(cohort_dir, checks, "cohort-frozen")
    normalization_entries = verify_checksum_manifest(frozen_dir, checks, "normalization-frozen")
    checks.add("frozen", "article19-entry-count", article19_entries == 34, article19_entries)
    checks.add("frozen", "cohort-entry-count", cohort_entries == 7, cohort_entries)
    checks.add("frozen", "normalization-entry-count", normalization_entries == 11, normalization_entries)

    version_result = subprocess.run(
        [str(prefix / "bin/humann"), "--version"],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONHASHSEED": "0"},
    )
    version_text = (version_result.stdout + version_result.stderr).strip()
    checks.add("environment", "humann-3.9", version_result.returncode == 0 and "3.9" in version_text, version_text)
    run_summary = json.loads((frozen_dir / "run-summary.json").read_text(encoding="utf-8"))
    checks.add("actual-tool", "four-frozen-branches", run_summary["actual_branches"] == 4, run_summary["actual_branches"])
    checks.add("actual-tool", "offline-run", run_summary["network_access"] is False, run_summary["network_access"])

    tables = {
        space: read_humann_table(article19_dir / filename)
        for space, filename in FEATURE_SPACES.items()
    }
    feature_rows, branch_rows, special_rows, closure_failures = build_normalization_audits(tables, checks)
    actual_error, actual_verified = compare_actual_outputs(
        tables["Pathway abundance"]["rows"], frozen_dir, branch_rows, checks
    )
    checks.add("actual-tool", "all-four-verified", actual_verified == 4, actual_verified)

    coverage_native = read_humann_table(article19_dir / "pathcoverage.tsv")
    mock_coverage_out = sum(
        not 0 <= row["Value"] <= 1
        for row in coverage_native["rows"]
        if not row["Special"]
    )
    checks.add("coverage", "mock-zero-to-one", mock_coverage_out == 0, mock_coverage_out)
    contribution_rows = pathway_contribution_audit(tables["Pathway abundance"]["rows"], checks)
    prevalence_rows, pseudocount_rows, cohort_summary, cohort_context = cohort_audits(cohort_dir, checks)
    cohort_lineage = cohort_context[:-1]
    matrix_context = cohort_context[-1]

    data_lineage = [
        {
            "Source": "Article 19 frozen HUMAnN output",
            "Resource": space,
            "Identifier": filename,
            "Release": "HUMAnN 3.9 / UniRef90 v201901b / MetaCyc",
            "Rows": len(tables[space]["rows"]),
            "Samples": 1,
            "SHA256": hash_file(article19_dir / filename),
            "Use": "Normalization denominator and contribution audit",
        }
        for space, filename in FEATURE_SPACES.items()
    ] + cohort_lineage

    write_tsv(output_dir / "data-lineage.tsv", data_lineage, list(data_lineage[0]))
    write_tsv(output_dir / "feature-space-audit.tsv", feature_rows, list(feature_rows[0]))
    write_tsv(output_dir / "normalization-branch-audit.tsv", branch_rows, list(branch_rows[0]))
    write_tsv(output_dir / "special-feature-audit.tsv", special_rows, list(special_rows[0]))
    write_tsv(output_dir / "pathway-contribution-audit.tsv", contribution_rows, list(contribution_rows[0]))
    write_tsv(output_dir / "prevalence-filter-audit.tsv", prevalence_rows, list(prevalence_rows[0]))
    write_tsv(output_dir / "zero-pseudocount-audit.tsv", pseudocount_rows, list(pseudocount_rows[0]))

    render_normalization_denominators(branch_rows, figure_dir)
    render_special_budget(special_rows, figure_dir)
    render_pathway_contributions(contribution_rows, figure_dir)
    render_prevalence_zero(prevalence_rows, pseudocount_rows, figure_dir)
    for stem in FIGURE_STEMS:
        for suffix in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            checks.add("figure", f"{stem}-{suffix}", path.is_file() and path.stat().st_size > 0, path.stat().st_size if path.is_file() else "missing")
        with Image.open(figure_dir / f"{stem}.tiff") as image:
            compression = image.info.get("compression", "")
            checks.add("figure", f"{stem}-tiff-lzw", str(compression).lower() in {"tiff_lzw", "lzw", "5"}, compression)

    summary = {
        "status": "passed" if checks.failed == 0 else "failed",
        "qa_network_access": False,
        "humann_version": "3.9",
        "normalization_feature_spaces": len(FEATURE_SPACES),
        "normalization_branches": len(branch_rows),
        "normalization_closure_failures": closure_failures,
        "actual_renorm_outputs_verified": actual_verified,
        "actual_renorm_max_relative_error": actual_error,
        "coverage_renormalized": False,
        "coverage_out_of_range": mock_coverage_out + cohort_summary["coverage_out_of_range"],
        "cohort_resource": "AsnicarF_2017 pathway_abundance/pathway_coverage",
        "cohort_package_version": "3.12.0",
        "cohort_pathway_rows": cohort_summary["rows"],
        "cohort_samples": cohort_summary["samples"],
        "cohort_subjects": cohort_summary["subjects"],
        "cohort_ordinary_unstratified_pathways": cohort_summary["ordinary_unstratified_pathways"],
        "pseudocounts": [1e-6, 1e-5, 1e-4],
        "pseudocount_numerator_pathway": cohort_summary["numerator_pathway"],
        "pseudocount_denominator_pathway": cohort_summary["denominator_pathway"],
        "biological_group_tests": 0,
        "article19_checksum_entries": article19_entries,
        "cohort_checksum_entries": cohort_entries,
        "normalization_checksum_entries": normalization_entries,
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "figures": list(FIGURE_STEMS),
    }
    write_tsv(
        output_dir / "validation-audit.tsv",
        checks.rows,
        ["Category", "CheckID", "Status", "Detail"],
    )
    write_json(output_dir / "validation-summary.json", summary)
    (output_dir / "validation.log").write_text(
        "Article 20 functional-profile normalization validation\n"
        f"Status: {summary['status']}\n"
        f"Checks passed: {checks.passed}\n"
        f"Checks failed: {checks.failed}\n"
        "Network access: false\n"
        "Biological group tests: 0\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    if checks.failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
