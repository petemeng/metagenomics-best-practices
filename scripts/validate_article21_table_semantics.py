#!/usr/bin/env python3
"""Validate Article 21 metagenomic table semantics and redraw its figures."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article21-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


EXPECTED_MICROBECENSUS_COMMIT = "dfc42d356bfd7943633cde6c0fbfc0b116f29ae2"
EXPECTED_TOTAL_BASES = 29_809_773
EXPECTED_SOURCE_READS = 199_982
EXPECTED_SOURCE_PAIRS = 99_991
SPECIAL_FEATURES = ("UNMAPPED", "UNINTEGRATED", "UNGROUPED")
FIGURE_STEMS = (
    "21-table-unit-map",
    "21-denominator-closure",
    "21-genome-equivalent-calibration",
    "21-zero-strata-semantics",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--article13-dir", type=Path, required=True)
    parser.add_argument("--article15-dir", type=Path, required=True)
    parser.add_argument("--article16-dir", type=Path, required=True)
    parser.add_argument("--article19-dir", type=Path, required=True)
    parser.add_argument("--cohort-dir", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
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


def verify_checksum_manifest(directory: Path, checks: Checks, category: str) -> int:
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
        f"payloads={len(payloads)};entries={len(entries)}",
    )
    return len(entries)


def feature_base(feature: str) -> str:
    return feature.split("|", 1)[0]


def feature_id(feature: str) -> str:
    return feature_base(feature).split(": ", 1)[0]


def is_special(feature: str) -> bool:
    return feature_id(feature) in SPECIAL_FEATURES


def read_humann_table(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("#"):
        raise ValueError(f"Unexpected HUMAnN header: {path}")
    for line in lines[1:]:
        if not line.strip():
            continue
        feature, value = line.split("\t")
        rows.append(
            {
                "Feature": feature,
                "FeatureID": feature_id(feature),
                "Level": len(feature.split("|")),
                "Special": is_special(feature),
                "Value": float(value),
            }
        )
    return rows


def read_metaphlan(path: Path) -> list[dict[str, Any]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header_index = next(
        index for index, line in enumerate(lines) if line.startswith("#clade_name\t")
    )
    header = [item.lstrip("#") for item in lines[header_index].split("\t")]
    rows: list[dict[str, Any]] = []
    for line in lines[header_index + 1 :]:
        if not line.strip():
            continue
        values = line.split("\t")
        row = dict(zip(header, values))
        last = row["clade_name"].split("|")[-1]
        row["rank"] = last.split("__", 1)[0] if "__" in last else last
        row["relative_abundance"] = float(row["relative_abundance"])
        row["estimated_number_of_reads_from_the_clade"] = int(
            row["estimated_number_of_reads_from_the_clade"]
        )
        rows.append(row)
    return rows


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


def build_contract() -> list[dict[str, Any]]:
    rows = [
        ("Kraken report", "direct reads", "read assignment", "taxonomy node", "paired fragments", "No", "all input fragments", "No", "0 = no direct assignment", "rank-aware count audit"),
        ("Kraken report", "clade reads", "subtree assignment", "taxonomy subtree", "paired fragments", "No", "all input fragments", "No", "0 = no subtree assignment", "classification ledger"),
        ("Bracken", "kraken_assigned_reads", "rank input", "species", "paired fragments", "No", "eligible rank model", "No", "0 = no eligible direct assignment", "model input audit"),
        ("Bracken", "added_reads", "redistributed mass", "species", "estimated paired fragments", "No", "database/read-length model", "No", "0 = no redistributed mass", "redistribution audit"),
        ("Bracken", "new_est_reads", "rank abundance estimate", "species", "estimated paired fragments", "No", "reported species estimates", "No", "0 = not estimated at rank", "conditional count-like analysis"),
        ("Bracken", "fraction_total_reads", "rank composition", "species", "fraction", "No", "sum of reported rank estimates", "Yes: 1", "0 = zero estimated fraction", "composition analysis"),
        ("MetaPhlAn", "relative_abundance", "marker-model composition", "one taxonomic rank", "percent", "Genome/marker model", "terminal clades at declared rank", "Yes: 100", "0 = not detected/modelled", "composition analysis"),
        ("MetaPhlAn", "estimated_number_of_reads_from_the_clade", "read-equivalent model", "taxonomic clade", "estimated reads", "Genome/marker model", "per-rank model ledger", "No raw-count closure", "0 = no modelled read equivalents", "model diagnostics"),
        ("HUMAnN", "Abundance-RPKs", "gene-family abundance", "gene family", "reads per kilobase", "Gene length", "none across samples", "No", "0 = no retained alignment support", "length-corrected continuous analysis"),
        ("HUMAnN", "Abundance-CPM", "gene-family composition", "gene family", "copies per million (CoPM)", "Gene length", "community RPK sum", "Yes: 1e6", "0 = no retained support", "TSS/compositional analysis"),
        ("HUMAnN", "relative abundance", "gene-family composition", "gene family", "fraction", "Gene length", "community RPK sum", "Yes: 1", "0 = no retained support", "TSS/compositional analysis"),
        ("HUMAnN", "pathway abundance", "reconstructed pathway potential", "MetaCyc pathway", "RPK-derived abundance", "Gene/reaction evidence", "none before renorm", "No", "0 = not reconstructed", "continuous functional analysis"),
        ("HUMAnN", "pathway coverage", "pathway detection/completeness evidence", "MetaCyc pathway", "bounded score", "Pathway rule", "none", "No; bounded 0-1", "0 = no sufficient pathway evidence", "detection/coverage analysis"),
        ("MicrobeCensus", "average_genome_size", "community genome-size estimate", "sample", "bp per genome", "Universal marker model", "marker-hit model", "No", "NA when marker evidence fails", "library calibration"),
        ("MicrobeCensus", "genome_equivalents", "aggregate sequenced genome coverage", "sample", "genome coverage", "Average genome size", "total bp / AGS", "No", "0 only for empty sequence input", "gene-abundance denominator"),
        ("Derived", "RPKG", "gene abundance per genome coverage", "gene family", "reads/kb/genome equivalent", "Gene length and AGS", "same sequence universe GE", "No", "0 = no retained gene support", "genome-calibrated continuous analysis"),
    ]
    fields = [
        "SourceTable",
        "Column",
        "MeasurementObject",
        "FeatureSpace",
        "NativeUnit",
        "Correction",
        "Denominator",
        "ClosureTarget",
        "ZeroMeaning",
        "AllowedUse",
    ]
    return [dict(zip(fields, row)) for row in rows]


def transformation_contract() -> list[dict[str, str]]:
    rows = [
        ("Observed fragment counts", "Sum within one nonoverlapping rank", "Allowed", "Keep sample and rank fixed; taxonomy hierarchy is not additive across ranks."),
        ("Observed fragment counts", "Count likelihood", "Conditional", "Use biological replicates and a library-size offset; database assignment is still uncertain."),
        ("Bracken estimates", "Count likelihood", "Conditional", "Values are redistributed estimates; a count model ignores re-estimation uncertainty."),
        ("Bracken estimates", "Treat as raw observations", "Forbidden", "new_est_reads includes model-added fragments."),
        ("MetaPhlAn relative abundance", "Close within one rank", "Allowed", "Filter to the terminal rank before checking 100%."),
        ("MetaPhlAn relative abundance", "Raw count likelihood", "Forbidden", "Percentages are marker/genome-model compositions."),
        ("MetaPhlAn estimated reads", "Read-level event audit", "Forbidden", "Read equivalents can exceed processed reads and use a balancing correction."),
        ("HUMAnN RPK", "Compare gene lengths within sample", "Allowed", "RPK corrects feature length but not sample sequencing depth."),
        ("HUMAnN RPK", "Compare samples without another denominator", "Forbidden", "Sequencing depth remains in the scale."),
        ("HUMAnN CoPM", "Interpret as counts per million", "Forbidden", "HUMAnN CPM means copies per million after RPK TSS; label it CoPM."),
        ("HUMAnN CoPM/relative", "Log-ratio or suitable continuous model", "Conditional", "Declare special rows, zero handling, and repeated-subject design."),
        ("Pathway coverage", "Sum or TSS across pathways", "Forbidden", "Coverage is an independent bounded score, not mass."),
        ("Pathway coverage", "Detection/prevalence analysis", "Conditional", "Predeclare the coverage threshold and analysis unit."),
        ("Genome equivalents", "Interpret as cells or biomass", "Forbidden", "GE is total sequenced bp divided by estimated AGS."),
        ("Genome equivalents", "Normalize matching-library gene RPK", "Allowed", "RPKG requires the exact same filtered sequence universe."),
        ("RPKG", "Interpret as copies per gram", "Forbidden", "Absolute concentration needs spike-in, qPCR/flow cytometry, mass, or volume information."),
        ("Stratified gene families", "Sum strata to community", "Allowed", "Only after verifying the additive gene-family contract."),
        ("Stratified pathways", "Force strata to community closure", "Forbidden", "Community and taxon pathways are reconstructed independently."),
    ]
    return [
        {"TableFamily": table, "Operation": operation, "Verdict": verdict, "Condition": condition}
        for table, operation, verdict, condition in rows
    ]


def render_unit_map(figure_dir: Path) -> None:
    configure_plot_style()
    labels = [
        "Kraken fragments",
        "Bracken estimates",
        "MetaPhlAn abundance",
        "HUMAnN RPK",
        "HUMAnN CoPM",
        "Pathway coverage",
        "RPKG",
    ]
    properties = [
        "Observed\nevents",
        "Model\nderived",
        "Length\ncorrected",
        "TSS\nclosed",
        "Bounded\n0-1",
        "Genome\ncalibrated",
    ]
    values = np.asarray(
        [
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 0],
            [0, 1, 1, 1, 0, 0],
            [0, 1, 1, 0, 0, 0],
            [0, 1, 1, 1, 0, 0],
            [0, 1, 0, 0, 1, 0],
            [0, 1, 1, 0, 0, 1],
        ],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(8.3, 4.6))
    ax.imshow(values, cmap=matplotlib.colors.ListedColormap(["#F1F1ED", "#2A9D8F"]), vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(properties)), properties)
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            ax.text(column, row, "Yes" if values[row, column] else "No", ha="center", va="center", color="white" if values[row, column] else "#555555", fontsize=8, fontweight="bold" if values[row, column] else "normal")
    ax.set_xticks(np.arange(-0.5, len(properties), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("A table name does not define its statistical meaning", loc="left", fontweight="bold", pad=18)
    fig.tight_layout()
    save_figure(fig, figure_dir, "21-table-unit-map")


def render_denominator_closure(
    closure_rows: list[dict[str, Any]], metaphlan_summary: dict[str, int], figure_dir: Path
) -> None:
    configure_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0), gridspec_kw={"width_ratios": [1.15, 1]})
    selected_names = ["Bracken fraction", "MetaPhlAn species %", "HUMAnN CoPM", "HUMAnN relative"]
    selected = [next(row for row in closure_rows if row["Metric"] == name) for name in selected_names]
    ratios = [row["ObservedSum"] / row["ClosureTarget"] for row in selected]
    colors = ["#457B9D", "#6A4C93", "#2A9D8F", "#E9C46A"]
    axes[0].barh(np.arange(len(selected)), ratios, color=colors)
    axes[0].axvline(1.0, color="#555555", linestyle="--", linewidth=0.9)
    axes[0].set_yticks(np.arange(len(selected)), selected_names)
    axes[0].invert_yaxis()
    axes[0].set_xlim(0.985, 1.006)
    axes[0].set_xlabel("Observed sum / declared closure target")
    axes[0].set_title("A  Closure applies only after fixing a denominator", loc="left", fontweight="bold")
    axes[0].grid(axis="x", color="#DDDDDD", linewidth=0.6)
    for index, ratio in enumerate(ratios):
        axes[0].text(ratio, index, f" {ratio:.5f}", va="center", fontsize=8)

    names = ["Processed\nreads", "Known-clade\nestimate", "UNCLASSIFIED\ncorrection", "Reconciled\ntotal"]
    values = [
        metaphlan_summary["processed_reads"],
        metaphlan_summary["known_estimated_reads"],
        metaphlan_summary["unclassified_correction"],
        metaphlan_summary["reconciled_reads"],
    ]
    axes[1].bar(np.arange(4), values, color=["#457B9D", "#E76F51", "#B7B7A4", "#2A9D8F"])
    axes[1].axhline(0, color="#555555", linewidth=0.8)
    axes[1].set_xticks(np.arange(4), names)
    axes[1].tick_params(axis="x", labelsize=8)
    axes[1].set_ylabel("Read-equivalent value")
    axes[1].set_title("B  MetaPhlAn read equivalents are a model ledger", loc="left", fontweight="bold")
    axes[1].grid(axis="y", color="#DDDDDD", linewidth=0.6)
    for index, value in enumerate(values):
        axes[1].text(index, value + (5000 if value >= 0 else -5000), f"{value:,}", ha="center", va="bottom" if value >= 0 else "top", fontsize=8)
    fig.tight_layout()
    save_figure(fig, figure_dir, "21-denominator-closure")


def render_genome_equivalent(
    genome_rows: list[dict[str, Any]], gene_rows: list[dict[str, Any]], figure_dir: Path
) -> None:
    configure_plot_style()
    fig = plt.figure(figsize=(11.6, 4.8))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.2, 1])
    ax0 = fig.add_subplot(grid[0, 0])
    branch_labels = [f"{row['ReadLengthBp']} bp" for row in genome_rows]
    x = np.arange(len(genome_rows))
    ags = [row["AverageGenomeSizeBp"] / 1e6 for row in genome_rows]
    ge = [row["GenomeEquivalents"] for row in genome_rows]
    bars = ax0.bar(x - 0.17, ags, width=0.34, color="#457B9D", label="Average genome size")
    ax0.set_xticks(x, branch_labels)
    ax0.set_ylabel("Average genome size (Mbp)", color="#457B9D")
    ax0.tick_params(axis="y", labelcolor="#457B9D")
    ax0.set_ylim(0, max(ags) * 1.28)
    ax0b = ax0.twinx()
    ax0b.spines["top"].set_visible(False)
    ax0b.bar(x + 0.17, ge, width=0.34, color="#E76F51", label="Genome equivalents")
    ax0b.set_ylabel("Genome equivalents", color="#E76F51")
    ax0b.tick_params(axis="y", labelcolor="#E76F51")
    ax0b.set_ylim(0, max(ge) * 1.28)
    ax0.set_title("A  Read-length sensitivity", loc="left", fontweight="bold")
    ax0.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    ax0.text(0.5, -0.22, "Total sequence: 29,809,773 bp", transform=ax0.transAxes, ha="center", fontsize=8)

    genes = [row["GeneFamily"] for row in gene_rows]
    y = np.arange(len(genes))
    ax1 = fig.add_subplot(grid[0, 1])
    ax1.barh(y, [row["CoPM"] for row in gene_rows], color="#2A9D8F")
    ax1.set_yticks(y, genes)
    ax1.invert_yaxis()
    ax1.set_xlabel("Copies per million")
    ax1.set_title("B  TSS denominator", loc="left", fontweight="bold")
    ax1.grid(axis="x", color="#DDDDDD", linewidth=0.6)

    ax2 = fig.add_subplot(grid[0, 2])
    ax2.barh(y, [row["RPKG"] for row in gene_rows], color="#F4A261")
    ax2.set_yticks(y, [""] * len(y))
    ax2.invert_yaxis()
    ax2.set_xlabel("RPK per genome equivalent")
    ax2.set_title("C  Genome-equivalent denominator", loc="left", fontweight="bold")
    ax2.grid(axis="x", color="#DDDDDD", linewidth=0.6)
    fig.tight_layout()
    save_figure(fig, figure_dir, "21-genome-equivalent-calibration")


def render_zero_strata(
    zero_rows: list[dict[str, Any]], strata_rows: list[dict[str, Any]], figure_dir: Path
) -> None:
    configure_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4.1))
    colors = ["#B7B7A4", "#E9C46A", "#6A4C93", "#2A9D8F"]
    x = np.arange(len(zero_rows))
    counts = [row["CellCount"] for row in zero_rows]
    axes[0].bar(x, counts, color=colors)
    axes[0].set_yscale("symlog", linthresh=1)
    axes[0].set_xticks(x, [row["StateLabel"] for row in zero_rows], rotation=18, ha="right")
    axes[0].set_ylabel("Pathway-profile cells")
    axes[0].set_title("A  Abundance and coverage zeros are not interchangeable", loc="left", fontweight="bold")
    axes[0].grid(axis="y", color="#DDDDDD", linewidth=0.6)
    for index, count in enumerate(counts):
        axes[0].text(index, max(count, 1) * 1.12, f"{count:,}", ha="center", va="bottom", fontsize=8)

    labels = [row["Output"] for row in strata_rows]
    features = [row["Features"] for row in strata_rows]
    additive = [row["Semantics"] == "Additive" for row in strata_rows]
    axes[1].bar(np.arange(len(labels)), features, color=["#2A9D8F" if item else "#E76F51" for item in additive])
    axes[1].set_yscale("log")
    axes[1].set_xticks(np.arange(len(labels)), labels)
    axes[1].set_ylabel("Ordinary community features")
    axes[1].set_ylim(40, 20000)
    axes[1].set_title("B  Stratification semantics depend on output type", loc="left", fontweight="bold", pad=12)
    axes[1].grid(axis="y", color="#DDDDDD", linewidth=0.6)
    for index, row in enumerate(strata_rows):
        label = "Additive" if row["Semantics"] == "Additive" else "Independent\nreconstruction"
        axes[1].text(
            index,
            max(55, row["Features"] / 2.2),
            f"{row['Features']:,}\n{label}",
            ha="center",
            va="center",
            fontsize=7.5,
            color="white",
            fontweight="bold",
        )
    fig.tight_layout()
    save_figure(fig, figure_dir, "21-zero-strata-semantics")


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    article13_dir = args.article13_dir.resolve()
    article15_dir = args.article15_dir.resolve()
    article16_dir = args.article16_dir.resolve()
    article19_dir = args.article19_dir.resolve()
    cohort_dir = args.cohort_dir.resolve()
    frozen_dir = args.frozen_dir.resolve()
    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    checksum_counts = {
        "article13": verify_checksum_manifest(article13_dir, checks, "article13-frozen"),
        "article15": verify_checksum_manifest(article15_dir, checks, "article15-frozen"),
        "article16": verify_checksum_manifest(article16_dir, checks, "article16-frozen"),
        "article19": verify_checksum_manifest(article19_dir, checks, "article19-frozen"),
        "cohort": verify_checksum_manifest(cohort_dir, checks, "cohort-frozen"),
        "article21": verify_checksum_manifest(frozen_dir, checks, "article21-frozen"),
    }

    article13 = json.loads((article13_dir / "run-summary.json").read_text(encoding="utf-8"))
    article15 = json.loads((article15_dir / "run-summary.json").read_text(encoding="utf-8"))
    article16 = json.loads((article16_dir / "run-summary.json").read_text(encoding="utf-8"))
    article19 = json.loads((article19_dir / "run-summary.json").read_text(encoding="utf-8"))
    article21 = json.loads((frozen_dir / "run-summary.json").read_text(encoding="utf-8"))
    checks.add("lineage", "same-source-pairs", all(item["input_pairs"] == EXPECTED_SOURCE_PAIRS for item in (article15, article16, article19)), [item["input_pairs"] for item in (article15, article16, article19)])
    checks.add("lineage", "article13-pairs", article13["retained_pairs"] == EXPECTED_SOURCE_PAIRS, article13["retained_pairs"])
    checks.add("lineage", "source-reads", article21["source_reads"] == EXPECTED_SOURCE_READS, article21["source_reads"])
    checks.add("lineage", "source-total-bases", article21["source_total_bases"] == EXPECTED_TOTAL_BASES, article21["source_total_bases"])
    checks.add("lineage", "raw-fastq-not-committed", article21["raw_fastq_committed"] is False, article21["raw_fastq_committed"])
    checks.add("lineage", "routine-qa-no-fastq", article21["qa_reads_fastq"] is False, article21["qa_reads_fastq"])
    checks.add("lineage", "routine-qa-no-network", article21["qa_network_access"] is False, article21["qa_network_access"])

    source_manifest = read_tsv(frozen_dir / "source-manifest.tsv")
    versions = {row["Component"]: row["Version"] for row in read_tsv(frozen_dir / "tool-versions.tsv")}
    commands = (frozen_dir / "commands.sh").read_text(encoding="utf-8")
    checks.add("microbecensus", "source-manifest-rows", len(source_manifest) == 3, len(source_manifest))
    checks.add("microbecensus", "source-tag", versions.get("MicrobeCensus source tag") == "v1.1.1", versions.get("MicrobeCensus source tag"))
    checks.add("microbecensus", "source-commit", versions.get("MicrobeCensus source commit") == EXPECTED_MICROBECENSUS_COMMIT, versions.get("MicrobeCensus source commit"))
    checks.add("microbecensus", "internal-version", versions.get("MicrobeCensus internal version") == "1.1.0", versions.get("MicrobeCensus internal version"))
    checks.add("microbecensus", "rapsearch-version", versions.get("RAPsearch2") == "2.15", versions.get("RAPsearch2"))
    checks.add("microbecensus", "shim-scope", versions.get("Compatibility shim") == "RAPsearch2 preflight bytes-to-text only", versions.get("Compatibility shim"))
    checks.add("microbecensus", "command-wrapper", "run_article21_microbecensus.sh" in commands, "run_article21_microbecensus.sh")

    genome_rows_raw = read_tsv(frozen_dir / "microbecensus-read-length-sensitivity.tsv")
    genome_rows: list[dict[str, Any]] = []
    for row in genome_rows_raw:
        parsed = {
            **row,
            "ReadLengthBp": int(row["ReadLengthBp"]),
            "ReadsInSequenceUniverse": int(row["ReadsInSequenceUniverse"]),
            "ReadsSampledForAGS": int(row["ReadsSampledForAGS"]),
            "ReadsTooShort": int(row["ReadsTooShort"]),
            "MarkerHits": int(row["MarkerHits"]),
            "AssignedMarkerReads": int(row["AssignedMarkerReads"]),
            "AverageGenomeSizeBp": float(row["AverageGenomeSizeBp"]),
            "TotalBases": int(row["TotalBases"]),
            "GenomeEquivalents": float(row["GenomeEquivalents"]),
            "NReadsCeiling": int(row["NReadsCeiling"]),
            "Threads": int(row["Threads"]),
        }
        genome_rows.append(parsed)
        length = parsed["ReadLengthBp"]
        checks.add("microbecensus", f"total-bases-{length}", parsed["TotalBases"] == EXPECTED_TOTAL_BASES, parsed["TotalBases"])
        checks.add("microbecensus", f"source-reads-{length}", parsed["ReadsInSequenceUniverse"] == EXPECTED_SOURCE_READS, parsed["ReadsInSequenceUniverse"])
        checks.add("microbecensus", f"all-reads-ceiling-{length}", parsed["NReadsCeiling"] > EXPECTED_SOURCE_READS, parsed["NReadsCeiling"])
        checks.add("microbecensus", f"ge-formula-{length}", math.isclose(parsed["GenomeEquivalents"], parsed["TotalBases"] / parsed["AverageGenomeSizeBp"], rel_tol=1e-12), parsed["GenomeEquivalents"])
        checks.add("microbecensus", f"no-qc-{length}", row["QualityFiltering"] == "No" and row["DuplicateFiltering"] == "No", f"quality={row['QualityFiltering']};duplicate={row['DuplicateFiltering']}")
        checks.add("microbecensus", f"marker-evidence-{length}", parsed["AssignedMarkerReads"] > 0 and parsed["MarkerHits"] >= parsed["AssignedMarkerReads"], f"hits={parsed['MarkerHits']};assigned={parsed['AssignedMarkerReads']}")
    checks.add("microbecensus", "two-read-length-branches", [row["ReadLengthBp"] for row in genome_rows] == [150, 100], [row["ReadLengthBp"] for row in genome_rows])
    primary = genome_rows[0]

    bracken = read_tsv(article16_dir / "bracken-species-r150-t10.tsv")
    for row in bracken:
        for column in ("kraken_assigned_reads", "added_reads", "new_est_reads"):
            row[column] = int(row[column])
        row["fraction_total_reads"] = float(row["fraction_total_reads"])
    bracken_identity_failures = sum(row["new_est_reads"] != row["kraken_assigned_reads"] + row["added_reads"] for row in bracken)
    bracken_assigned = sum(row["kraken_assigned_reads"] for row in bracken)
    bracken_added = sum(row["added_reads"] for row in bracken)
    bracken_estimated = sum(row["new_est_reads"] for row in bracken)
    bracken_fraction_sum = math.fsum(row["fraction_total_reads"] for row in bracken)
    checks.add("bracken", "species-features", len(bracken) == 65, len(bracken))
    checks.add("bracken", "row-identities", bracken_identity_failures == 0, bracken_identity_failures)
    checks.add("bracken", "assigned-plus-added", bracken_assigned + bracken_added == bracken_estimated == 84_147, f"{bracken_assigned}+{bracken_added}={bracken_estimated}")
    checks.add("bracken", "rounded-fraction-closure", abs(bracken_fraction_sum - 1) < 5e-4, bracken_fraction_sum)

    metaphlan = read_metaphlan(article15_dir / "profile-all.tsv")
    species = [row for row in metaphlan if row["rank"] == "s"]
    unclassified = next(row for row in metaphlan if row["clade_name"] == "UNCLASSIFIED")
    metaphlan_percent_sum = math.fsum(row["relative_abundance"] for row in species)
    metaphlan_estimated_sum = sum(row["estimated_number_of_reads_from_the_clade"] for row in species)
    processed_reads = int(article15["profiling_reads"])
    unclassified_correction = int(unclassified["estimated_number_of_reads_from_the_clade"])
    reconciled_reads = metaphlan_estimated_sum + unclassified_correction
    checks.add("metaphlan", "species-features", len(species) == 31, len(species))
    checks.add("metaphlan", "species-percent-closure", abs(metaphlan_percent_sum - 100) < 1e-3, metaphlan_percent_sum)
    checks.add("metaphlan", "known-estimate", metaphlan_estimated_sum == 261_721, metaphlan_estimated_sum)
    checks.add("metaphlan", "negative-unclassified-correction", unclassified_correction == -61_792, unclassified_correction)
    checks.add("metaphlan", "read-equivalent-reconciliation", reconciled_reads == processed_reads == 199_929, reconciled_reads)
    checks.add("metaphlan", "estimate-exceeds-input", metaphlan_estimated_sum > processed_reads, f"estimate={metaphlan_estimated_sum};input={processed_reads}")

    gene_rpk = read_humann_table(article19_dir / "genefamilies-rpk.tsv")
    gene_copm = read_humann_table(article19_dir / "genefamilies-cpm.tsv")
    gene_relab = read_humann_table(article19_dir / "genefamilies-relab.tsv")
    pathway_rpk = read_humann_table(article19_dir / "pathabundance-rpk.tsv")
    pathway_coverage = read_humann_table(article19_dir / "pathcoverage.tsv")
    checks.add("humann", "aligned-feature-order-copm", [row["Feature"] for row in gene_rpk] == [row["Feature"] for row in gene_copm], len(gene_rpk))
    checks.add("humann", "aligned-feature-order-relab", [row["Feature"] for row in gene_rpk] == [row["Feature"] for row in gene_relab], len(gene_rpk))
    community_copm = math.fsum(row["Value"] for row in gene_copm if row["Level"] == 1)
    community_relab = math.fsum(row["Value"] for row in gene_relab if row["Level"] == 1)
    community_rpk = math.fsum(row["Value"] for row in gene_rpk if row["Level"] == 1)
    copm_relab_error = max(abs(c["Value"] - r["Value"] * 1_000_000) for c, r in zip(gene_copm, gene_relab))
    coverage_ordinary = [row for row in pathway_coverage if row["Level"] == 1 and not row["Special"]]
    coverage_out_of_range = sum(not 0 <= row["Value"] <= 1 for row in coverage_ordinary)
    coverage_sum = math.fsum(row["Value"] for row in coverage_ordinary)
    checks.add("humann", "community-copm-closure", abs(community_copm - 1_000_000) < 1, community_copm)
    checks.add("humann", "community-relab-closure", abs(community_relab - 1) < 1e-6, community_relab)
    checks.add("humann", "copm-relab-proportional", copm_relab_error < 1e-8, copm_relab_error)
    checks.add("humann", "coverage-bounded", coverage_out_of_range == 0, coverage_out_of_range)
    checks.add("humann", "coverage-not-closed", not math.isclose(coverage_sum, 1, rel_tol=1e-3), coverage_sum)
    pathway_rpk_features = [row["Feature"] for row in pathway_rpk]
    pathway_coverage_features = [row["Feature"] for row in pathway_coverage]
    pathway_order_differences = sum(
        left != right
        for left, right in zip(pathway_rpk_features, pathway_coverage_features)
    )
    checks.add(
        "humann",
        "pathway-feature-sets-aligned",
        set(pathway_rpk_features) == set(pathway_coverage_features),
        len(pathway_rpk_features),
    )
    checks.add(
        "humann",
        "pathway-order-requires-id-alignment",
        pathway_order_differences > 0,
        pathway_order_differences,
    )

    ordinary_gene_indices = [
        index
        for index, row in enumerate(gene_rpk)
        if row["Level"] == 1 and not row["Special"]
    ]
    top_indices = sorted(ordinary_gene_indices, key=lambda index: gene_rpk[index]["Value"], reverse=True)[:10]
    gene_rows: list[dict[str, Any]] = []
    for index in top_indices:
        rpk = gene_rpk[index]["Value"]
        gene_rows.append(
            {
                "GeneFamily": gene_rpk[index]["FeatureID"],
                "NativeRPK": rpk,
                "CoPM": gene_copm[index]["Value"],
                "RelativeAbundance": gene_relab[index]["Value"],
                "PrimaryGenomeEquivalents": primary["GenomeEquivalents"],
                "RPKG": rpk / primary["GenomeEquivalents"],
                "Formula": "RPK / genome equivalents",
                "SequenceUniverse": "ERR9765746_MOCK1 199982 clean reads",
            }
        )
    checks.add("rpkg", "top-ten-ordinary", len(gene_rows) == 10 and all(not row["GeneFamily"].startswith("UN") for row in gene_rows), len(gene_rows))
    checks.add("rpkg", "formula", all(math.isclose(row["RPKG"], row["NativeRPK"] / row["PrimaryGenomeEquivalents"], rel_tol=1e-12) for row in gene_rows), "RPK/GE")

    abundance_features, abundance_samples, abundance = read_matrix(cohort_dir / "pathway-abundance.tsv.gz")
    coverage_features, coverage_samples, cohort_coverage = read_matrix(cohort_dir / "pathway-coverage.tsv.gz")
    metadata = read_tsv(cohort_dir / "sample-metadata.tsv")
    metadata_samples = [row["sample_id"] for row in metadata]
    subjects = {row["subject_id"] for row in metadata}
    checks.add("cohort", "abundance-dimensions", abundance.shape == (11_173, 24), abundance.shape)
    checks.add("cohort", "coverage-dimensions", cohort_coverage.shape == (11_173, 24), cohort_coverage.shape)
    checks.add("cohort", "feature-order", abundance_features == coverage_features, len(abundance_features))
    checks.add("cohort", "sample-order", abundance_samples == coverage_samples == metadata_samples, len(metadata_samples))
    checks.add("cohort", "subject-count", len(subjects) == 15, len(subjects))
    checks.add("cohort", "coverage-range", int(np.sum((cohort_coverage < 0) | (cohort_coverage > 1))) == 0, (float(cohort_coverage.min()), float(cohort_coverage.max())))
    ordinary_keep = np.asarray(["|" not in feature and not is_special(feature) for feature in abundance_features], dtype=bool)
    ordinary_abundance = abundance[ordinary_keep]
    ordinary_coverage = cohort_coverage[ordinary_keep]
    checks.add("cohort", "ordinary-unstratified-pathways", ordinary_abundance.shape == (445, 24), ordinary_abundance.shape)
    states = [
        ("Abundance = 0\nCoverage = 0", (ordinary_abundance == 0) & (ordinary_coverage == 0), "No abundance or coverage evidence"),
        ("Abundance > 0\nCoverage = 0", (ordinary_abundance > 0) & (ordinary_coverage == 0), "Abundance reconstructed below coverage evidence"),
        ("Abundance = 0\nCoverage > 0", (ordinary_abundance == 0) & (ordinary_coverage > 0), "Coverage evidence without abundance"),
        ("Abundance > 0\nCoverage > 0", (ordinary_abundance > 0) & (ordinary_coverage > 0), "Both evidence types positive"),
    ]
    zero_rows = [
        {
            "StateLabel": label,
            "CellCount": int(np.sum(mask)),
            "FractionOfOrdinaryCells": float(np.mean(mask)),
            "Interpretation": interpretation,
            "BiologicalAbsenceClaim": "No",
        }
        for label, mask, interpretation in states
    ]
    checks.add("zero", "state-partition", sum(row["CellCount"] for row in zero_rows) == 445 * 24, sum(row["CellCount"] for row in zero_rows))
    checks.add("zero", "abundance-positive-coverage-zero-exists", zero_rows[1]["CellCount"] > 0, zero_rows[1]["CellCount"])
    checks.add("zero", "no-biological-absence-claim", all(row["BiologicalAbsenceClaim"] == "No" for row in zero_rows), "all states")

    strata_source = read_tsv(article19_dir / "stratification-audit.tsv")
    strata_rows: list[dict[str, Any]] = []
    for row in strata_source:
        output = row["Output"]
        strata_rows.append(
            {
                "Output": output,
                "Features": int(row["Features"]),
                "Violations": int(row["Violations"]),
                "NonadditiveFeatures": int(row["NonadditiveFeatures"]),
                "MedianStrataToCommunityRatio": float(row["MedianStrataToCommunityRatio"]),
                "ExpectedRelationship": row["ExpectedRelationship"],
                "Semantics": "Independent reconstruction" if output == "Pathway abundance" else "Additive",
                "CanForceClosure": "No" if output == "Pathway abundance" else "Yes after audit",
            }
        )
    checks.add("strata", "three-output-types", len(strata_rows) == 3, len(strata_rows))
    checks.add("strata", "gene-reaction-additive", all(row["Violations"] == 0 and row["Semantics"] == "Additive" for row in strata_rows[:2]), strata_rows[:2])
    checks.add("strata", "pathway-nonadditive", strata_rows[2]["NonadditiveFeatures"] == strata_rows[2]["Features"] == 147 and strata_rows[2]["CanForceClosure"] == "No", strata_rows[2])

    closure_rows = [
        {"Metric": "Bracken fraction", "FeatureLevel": "Species", "ObservedSum": bracken_fraction_sum, "ClosureTarget": 1.0, "ClosureStatus": "Rounded pass", "Denominator": "Sum of reported species estimates", "MayMixRanks": "No"},
        {"Metric": "MetaPhlAn species %", "FeatureLevel": "Species", "ObservedSum": metaphlan_percent_sum, "ClosureTarget": 100.0, "ClosureStatus": "Pass", "Denominator": "Species-level marker model", "MayMixRanks": "No"},
        {"Metric": "HUMAnN CoPM", "FeatureLevel": "Community gene families", "ObservedSum": community_copm, "ClosureTarget": 1_000_000.0, "ClosureStatus": "Pass", "Denominator": "Community gene-family RPK sum", "MayMixRanks": "No"},
        {"Metric": "HUMAnN relative", "FeatureLevel": "Community gene families", "ObservedSum": community_relab, "ClosureTarget": 1.0, "ClosureStatus": "Pass", "Denominator": "Community gene-family RPK sum", "MayMixRanks": "No"},
        {"Metric": "Pathway coverage", "FeatureLevel": "Ordinary community pathways", "ObservedSum": coverage_sum, "ClosureTarget": "NA", "ClosureStatus": "Not applicable", "Denominator": "Independent bounded pathway evidence", "MayMixRanks": "No"},
    ]
    native_rows = [
        {"Source": "Bracken", "Table": "bracken-species-r150-t10.tsv", "FeatureLevel": "Species", "NativeUnit": "Estimated paired fragments", "Rows": len(bracken), "NativeSum": bracken_estimated, "ObservedOrModel": "Model-derived", "LengthCorrected": "No", "DepthNormalized": "No", "Closed": "Only fraction_total_reads"},
        {"Source": "MetaPhlAn", "Table": "profile-all.tsv", "FeatureLevel": "Species", "NativeUnit": "Percent", "Rows": len(species), "NativeSum": metaphlan_percent_sum, "ObservedOrModel": "Model-derived", "LengthCorrected": "Marker/genome model", "DepthNormalized": "Yes: rank composition", "Closed": "Yes within rank"},
        {"Source": "HUMAnN", "Table": "genefamilies-rpk.tsv", "FeatureLevel": "Community gene families", "NativeUnit": "RPK", "Rows": len([row for row in gene_rpk if row["Level"] == 1]), "NativeSum": community_rpk, "ObservedOrModel": "Alignment-derived", "LengthCorrected": "Yes", "DepthNormalized": "No", "Closed": "No"},
        {"Source": "HUMAnN", "Table": "genefamilies-cpm.tsv", "FeatureLevel": "Community gene families", "NativeUnit": "Copies per million (CoPM)", "Rows": len([row for row in gene_copm if row["Level"] == 1]), "NativeSum": community_copm, "ObservedOrModel": "RPK-derived", "LengthCorrected": "Yes", "DepthNormalized": "TSS", "Closed": "1,000,000"},
        {"Source": "HUMAnN", "Table": "pathcoverage.tsv", "FeatureLevel": "Ordinary community pathways", "NativeUnit": "0-1 coverage score", "Rows": len(coverage_ordinary), "NativeSum": coverage_sum, "ObservedOrModel": "Pathway-model evidence", "LengthCorrected": "No", "DepthNormalized": "No TSS", "Closed": "Not applicable"},
        {"Source": "MicrobeCensus", "Table": "microbecensus-primary.tsv", "FeatureLevel": "Sample", "NativeUnit": "Genome equivalents", "Rows": 1, "NativeSum": primary["GenomeEquivalents"], "ObservedOrModel": "AGS-model-derived", "LengthCorrected": "Total bp / AGS", "DepthNormalized": "Library calibration", "Closed": "No"},
    ]

    contract_rows = build_contract()
    transform_rows = transformation_contract()
    checks.add("contract", "sixteen-unit-rows", len(contract_rows) == 16, len(contract_rows))
    checks.add("contract", "all-required-fields", all(all(str(row[field]).strip() for field in row) for row in contract_rows), len(contract_rows))
    checks.add("contract", "copm-expanded", any(row["NativeUnit"] == "copies per million (CoPM)" for row in contract_rows), "CoPM")
    checks.add("contract", "coverage-no-closure", any(row["Column"] == "pathway coverage" and str(row["ClosureTarget"]).startswith("No") for row in contract_rows), "pathway coverage")
    checks.add("contract", "transformation-decisions", len(transform_rows) >= 18, len(transform_rows))
    checks.add("contract", "forbidden-cases", sum(row["Verdict"] == "Forbidden" for row in transform_rows) >= 8, sum(row["Verdict"] == "Forbidden" for row in transform_rows))

    data_lineage = [
        {"Source": "Article 13 clean MOCK1", "Identifier": "PRJEB52977/SAMEA14435832/ERR9765746", "Release": "first 100,000 pairs; 99,991 retained", "RowsOrReads": EXPECTED_SOURCE_READS, "Samples": 1, "Unit": "reads", "Use": "shared sequence universe; FASTQ excluded from routine QA"},
        {"Source": "Article 15 MetaPhlAn", "Identifier": "profile-all.tsv", "Release": "MetaPhlAn 4.2.5 / vJan26", "RowsOrReads": len(metaphlan), "Samples": 1, "Unit": "percent and estimated reads", "Use": "model-ledger semantics"},
        {"Source": "Article 16 Bracken", "Identifier": "bracken-species-r150-t10.tsv", "Release": "Kraken2 2.17.1 / Bracken 3.1p1", "RowsOrReads": len(bracken), "Samples": 1, "Unit": "estimated paired fragments", "Use": "observed-versus-estimated audit"},
        {"Source": "Article 19 HUMAnN", "Identifier": "gene/pathway tables", "Release": "HUMAnN 3.9 / UniRef90 v201901b", "RowsOrReads": len(gene_rpk), "Samples": 1, "Unit": "RPK/CoPM/relative/coverage", "Use": "length, TSS, coverage and strata semantics"},
        {"Source": "Article 20 cMD", "Identifier": "AsnicarF_2017 EH7089/EH7090", "Release": "curatedMetagenomicData 3.12.0", "RowsOrReads": len(abundance_features), "Samples": len(abundance_samples), "Unit": "pathway abundance and coverage", "Use": "multi-profile zero-state audit"},
        {"Source": "Article 21 MicrobeCensus", "Identifier": EXPECTED_MICROBECENSUS_COMMIT, "Release": "tag v1.1.1; internal 1.1.0", "RowsOrReads": EXPECTED_SOURCE_READS, "Samples": 1, "Unit": "average genome size and genome equivalents", "Use": "actual matching-library calibration"},
    ]

    write_tsv(output_dir / "data-lineage.tsv", data_lineage, list(data_lineage[0]))
    write_tsv(output_dir / "table-semantics-contract.tsv", contract_rows, list(contract_rows[0]))
    write_tsv(output_dir / "native-unit-audit.tsv", native_rows, list(native_rows[0]))
    write_tsv(output_dir / "closure-denominator-audit.tsv", closure_rows, list(closure_rows[0]))
    write_tsv(output_dir / "genome-equivalent-audit.tsv", genome_rows, list(genome_rows[0]))
    write_tsv(output_dir / "gene-rpkg-audit.tsv", gene_rows, list(gene_rows[0]))
    write_tsv(output_dir / "zero-state-audit.tsv", zero_rows, list(zero_rows[0]))
    write_tsv(output_dir / "stratification-semantics-audit.tsv", strata_rows, list(strata_rows[0]))
    write_tsv(output_dir / "transformation-legality.tsv", transform_rows, list(transform_rows[0]))

    metaphlan_summary = {
        "processed_reads": processed_reads,
        "known_estimated_reads": metaphlan_estimated_sum,
        "unclassified_correction": unclassified_correction,
        "reconciled_reads": reconciled_reads,
    }
    render_unit_map(figure_dir)
    render_denominator_closure(closure_rows, metaphlan_summary, figure_dir)
    render_genome_equivalent(genome_rows, gene_rows, figure_dir)
    render_zero_strata(zero_rows, strata_rows, figure_dir)
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
        "qa_reads_fastq": False,
        "source_pairs": EXPECTED_SOURCE_PAIRS,
        "source_reads": EXPECTED_SOURCE_READS,
        "primary_total_bases": int(primary["TotalBases"]),
        "microbecensus_source_tag": versions["MicrobeCensus source tag"],
        "microbecensus_source_commit": versions["MicrobeCensus source commit"],
        "microbecensus_internal_version": versions["MicrobeCensus internal version"],
        "microbecensus_branches": len(genome_rows),
        "primary_read_length_bp": int(primary["ReadLengthBp"]),
        "primary_reads_sampled_for_ags": int(primary["ReadsSampledForAGS"]),
        "primary_average_genome_size_bp": float(primary["AverageGenomeSizeBp"]),
        "primary_genome_equivalents": float(primary["GenomeEquivalents"]),
        "sensitivity_average_genome_size_bp": float(genome_rows[1]["AverageGenomeSizeBp"]),
        "sensitivity_genome_equivalents": float(genome_rows[1]["GenomeEquivalents"]),
        "bracken_species_features": len(bracken),
        "bracken_assigned_fragments": bracken_assigned,
        "bracken_added_fragments": bracken_added,
        "bracken_estimated_fragments": bracken_estimated,
        "bracken_fraction_sum": bracken_fraction_sum,
        "metaphlan_species_features": len(species),
        "metaphlan_species_percent_sum": metaphlan_percent_sum,
        "metaphlan_processed_reads": processed_reads,
        "metaphlan_known_estimated_reads": metaphlan_estimated_sum,
        "metaphlan_unclassified_correction": unclassified_correction,
        "humann_gene_rpk_community_sum": community_rpk,
        "humann_gene_copm_community_sum": community_copm,
        "humann_gene_relative_community_sum": community_relab,
        "humann_copm_relative_max_error": copm_relab_error,
        "pathway_coverage_sum_no_closure_target": coverage_sum,
        "coverage_out_of_range": coverage_out_of_range,
        "cohort_pathway_rows": len(abundance_features),
        "cohort_profiles": len(abundance_samples),
        "cohort_subjects": len(subjects),
        "cohort_ordinary_unstratified_pathways": int(ordinary_abundance.shape[0]),
        "zero_state_cells": int(ordinary_abundance.size),
        "biological_group_tests": 0,
        "table_semantics_rows": len(contract_rows),
        "transformation_decisions": len(transform_rows),
        "checksum_entries": checksum_counts,
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "figures": list(FIGURE_STEMS),
    }
    write_tsv(output_dir / "validation-audit.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])
    write_json(output_dir / "validation-summary.json", summary)
    (output_dir / "validation.log").write_text(
        "Article 21 metagenomic table semantics validation\n"
        f"Status: {summary['status']}\n"
        f"Checks passed: {checks.passed}\n"
        f"Checks failed: {checks.failed}\n"
        "Network access: false\n"
        "FASTQ access: false\n"
        "Biological group tests: 0\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    if checks.failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
