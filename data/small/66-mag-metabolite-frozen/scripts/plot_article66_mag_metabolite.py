#!/usr/bin/env python3
"""Create publication-ready Article 66 figures from frozen evidence tables."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PLOT_SEED = 20_260_766
REFERENCE = "#0072B2"
MAG = "#D55E00"
SHARED = "#8C8C8C"
DONOR_COLORS = {"Donor 1": "#009E73", "Donor 2": "#CC79A7"}
FIGURES = (
    "66-mag-quality-landscape",
    "66-input-normalization-audit",
    "66-individual-gem-gapfill",
    "66-metabolite-overlap",
    "66-validation-boundary",
    "66-pathway-robustness",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    return parser.parse_args()


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def configure() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def save(fig: plt.Figure, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(target.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(
        target.with_suffix(".tiff"),
        dpi=350,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def panel(ax: plt.Axes, letter: str, title: str) -> None:
    ax.set_title(f"{letter}. {title}", loc="left", pad=8)


def plot_quality(input_dir: Path, figure_dir: Path) -> None:
    data = read(input_dir / "mag-ledger.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), sharex=True, sharey=True)
    for ax, donor, letter in zip(axes, ("Donor 1", "Donor 2"), ("A", "B"), strict=True):
        subset = data[data["Donor"].eq(donor)]
        for quality, color in (("Paper high", "#0072B2"), ("Paper medium", "#E69F00")):
            group = subset[subset["PaperQuality"].eq(quality)]
            ax.scatter(
                group["Completeness"],
                group["Contamination"],
                s=32,
                color=color,
                alpha=0.72,
                edgecolor="white",
                linewidth=0.35,
                label=quality,
            )
        selected = subset[subset["NormalizedInput"].astype(bool)]
        ax.scatter(
            selected["Completeness"],
            selected["Contamination"],
            s=64,
            facecolor="none",
            edgecolor="black",
            linewidth=0.8,
            label="Normalized input",
        )
        ax.axvline(90, color="#444444", linestyle="--", linewidth=0.9)
        ax.axhline(5, color="#444444", linestyle="--", linewidth=0.9)
        ax.set_xlim(74, 101)
        ax.set_ylim(-0.15, 5.55)
        ax.set_xlabel("CheckM completeness (%)")
        ax.set_ylabel("CheckM contamination (%)" if ax is axes[0] else "")
        ax.grid(alpha=0.22)
        high = int(subset["PaperQuality"].eq("Paper high").sum())
        medium = int(subset["PaperQuality"].eq("Paper medium").sum())
        normalized = int(subset["NormalizedInput"].astype(bool).sum())
        ax.text(
            0.03,
            0.96,
            f"High: {high}  ·  Medium: {medium}\nNormalized: {normalized}",
            transform=ax.transAxes,
            va="top",
            fontsize=8.2,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D7DE"},
        )
        panel(ax, letter, f"{donor} MAG quality and selected inputs")
    handles, labels = axes[1].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.055))
    fig.text(
        0.5,
        0.015,
        "Published high-quality labels use completeness/contamination only; MIMAG rRNA/tRNA evidence is unavailable.",
        ha="center",
        fontsize=8.2,
        color="#37474F",
    )
    fig.tight_layout(rect=(0, 0.145, 1, 1), w_pad=1.6)
    save(fig, figure_dir / "66-mag-quality-landscape")


def plot_normalization(input_dir: Path, figure_dir: Path) -> None:
    data = read(input_dir / "community-model-audit.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(11.1, 4.55), sharey=True)
    for ax, donor, letter in zip(axes, ("Donor 1", "Donor 2"), ("A", "B"), strict=True):
        subset = data[data["Donor"].eq(donor)].copy()
        ref = subset[subset["Approach"].str.startswith("Reference")]
        mag = subset[subset["Approach"].str.startswith("MAG")]
        ax.plot(ref["Taxa"], ref["Unique metabolites"], color=REFERENCE, marker="o", linewidth=1.6, label="Reference-guided")
        ax.plot(mag["Taxa"], mag["Unique metabolites"], color=MAG, marker="s", linewidth=1.2, linestyle="--", label="MAG-guided")
        normalized = subset[subset["Approach"].isin(["Reference ≥0.5%", "MAG normalized"])].sort_values("Unique metabolites")
        ax.scatter(normalized["Taxa"], normalized["Unique metabolites"], s=95, facecolor="none", edgecolor="black", linewidth=1.1, zorder=5)
        values = normalized.set_index("Approach")["Unique metabolites"]
        difference = int(values["MAG normalized"] - values["Reference ≥0.5%"])
        x = float(normalized["Taxa"].iloc[0])
        y = float(normalized["Unique metabolites"].max()) + 20
        ax.text(x + 1.5, y, f"Matched inputs: Δ = {difference:+d}", fontsize=8.2, ha="left")
        ax.set_xlabel("Input genomes or taxa")
        ax.set_ylabel("Non-redundant predicted metabolites" if ax is axes[0] else "")
        ax.set_ylim(1020, 1625)
        ax.grid(alpha=0.22)
        panel(ax, letter, f"{donor}: predictions approach saturation")
    axes[0].legend(frameon=False, loc="lower right")
    fig.tight_layout(w_pad=1.8)
    save(fig, figure_dir / "66-input-normalization-audit")


def plot_gapfill(input_dir: Path, figure_dir: Path) -> None:
    data = read(input_dir / "individual-model-audit.tsv")
    palette = {"Reference-guided": REFERENCE, "MAG-guided": MAG}
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), sharex=True, sharey=True)
    for ax, donor, letter in zip(axes, ("Donor 1", "Donor 2"), ("A", "B"), strict=True):
        subset = data[data["Donor"].eq(donor)]
        for approach in ("Reference-guided", "MAG-guided"):
            group = subset[subset["Approach"].eq(approach)]
            valid = group[group["GapfillDeltaValid"].astype(bool)]
            ax.scatter(
                valid["CompoundsBeforeGapfill"],
                valid["GapfillAddedCompounds"],
                s=31,
                alpha=0.72,
                color=palette[approach],
                edgecolor="white",
                linewidth=0.35,
                label=approach,
            )
            ax.scatter(
                [valid["CompoundsBeforeGapfill"].median()],
                [valid["GapfillAddedCompounds"].median()],
                marker="X",
                s=95,
                color=palette[approach],
                edgecolor="black",
                linewidth=0.7,
                zorder=5,
            )
        invalid = subset[~subset["GapfillDeltaValid"].astype(bool)]
        if len(invalid):
            ax.scatter(
                invalid["CompoundsBeforeGapfill"],
                invalid["GapfillAddedCompounds"],
                marker="v",
                s=80,
                color="#CC79A7",
                edgecolor="black",
                linewidth=0.7,
                zorder=6,
                label="Published count anomaly",
            )
            row = invalid.iloc[0]
            ax.annotate(
                "Published after < before",
                (row["CompoundsBeforeGapfill"], row["GapfillAddedCompounds"]),
                xytext=(18, 8),
                textcoords="offset points",
                fontsize=7.8,
                arrowprops={"arrowstyle": "-", "color": "#555555", "linewidth": 0.7},
            )
        ax.set_xlim(120, 1535)
        ax.set_ylim(-70, 250)
        ax.set_xlabel("Compounds before gap filling")
        ax.set_ylabel("Compounds added by gap filling" if ax is axes[0] else "")
        ax.grid(alpha=0.22)
        ax.text(0.98, 0.04, "X = group median", transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
        panel(ax, letter, f"{donor}: gap-filling burden is model-specific")
    axes[0].legend(frameon=False, loc="upper right")
    fig.tight_layout(w_pad=1.7)
    save(fig, figure_dir / "66-individual-gem-gapfill")


def plot_overlap(input_dir: Path, figure_dir: Path) -> None:
    data = read(input_dir / "metabolite-overlap-audit.tsv")
    data["Label"] = data["Donor"] + " · " + data["Evidence"]
    order = ["Donor 1 · Predicted", "Donor 2 · Predicted", "Donor 1 · Confirmed", "Donor 2 · Confirmed"]
    data = data.set_index("Label").loc[order].reset_index()
    fig, ax = plt.subplots(figsize=(10.6, 4.8))
    y = np.arange(len(data))
    left = np.zeros(len(data))
    for column, label, color in (
        ("ReferenceOnly", "Reference only", REFERENCE),
        ("Shared", "Shared", SHARED),
        ("MAGOnly", "MAG only", MAG),
    ):
        values = data[column].to_numpy(float)
        bars = ax.barh(y, values, left=left, color=color, label=label, edgecolor="white", height=0.62)
        for bar, value in zip(bars, values, strict=True):
            if value >= 3:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_y() + bar.get_height() / 2, f"{int(value):,}", ha="center", va="center", color="white" if column != "ReferenceOnly" else "white", fontsize=8.3, fontweight="bold")
        left += values
    ax.set_yticks(y, order)
    ax.invert_yaxis()
    ax.set_xlabel("Union of non-redundant metabolites")
    ax.set_ylabel("")
    ax.grid(axis="x", alpha=0.22)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, data["Union"].max() * 1.13)
    for i, row in data.iterrows():
        ax.text(row["Union"] + data["Union"].max() * 0.015, i, f"J = {row['Jaccard']:.3f}", va="center", fontsize=8.5)
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.25))
    panel(ax, "A", "High overlap does not make either input strategy complete")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, figure_dir / "66-metabolite-overlap")


def plot_validation(input_dir: Path, figure_dir: Path) -> None:
    data = read(input_dir / "community-model-audit.tsv").copy()
    data["Family"] = np.where(data["Approach"].str.startswith("Reference"), "Reference-guided", "MAG-guided")
    palette = {"Reference-guided": REFERENCE, "MAG-guided": MAG}
    markers = {"Donor 1": "o", "Donor 2": "s"}
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.65))
    for donor, marker in markers.items():
        subset = data[data["Donor"].eq(donor)]
        for family in ("Reference-guided", "MAG-guided"):
            group = subset[subset["Family"].eq(family)]
            axes[0].scatter(group["Taxa"], group["Data loss (%)"], color=palette[family], marker=marker, s=47, alpha=0.8, edgecolor="white", linewidth=0.4)
            axes[1].scatter(group["Taxa"], group["Confirmed unique metabolites (%)"], color=palette[family], marker=marker, s=47, alpha=0.8, edgecolor="white", linewidth=0.4)
    axes[0].set_xlabel("Input genomes or taxa")
    axes[0].set_ylabel("Unmatched identifiers / naming loss (%)")
    axes[0].set_ylim(26.5, 32.2)
    panel(axes[0], "A", "About 30% of predictions were not directly searchable")
    axes[1].set_xlabel("Input genomes or taxa")
    axes[1].set_ylabel("Confirmed unique metabolites (%)")
    axes[1].set_ylim(15.8, 20.2)
    panel(axes[1], "B", "Confirmed fraction is not sensitivity or accuracy")
    for ax in axes:
        ax.grid(alpha=0.22)
    family_handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", color=color, label=label)
        for label, color in palette.items()
    ]
    donor_handles = [
        mpl.lines.Line2D([], [], marker=marker, linestyle="", color="#555555", label=donor)
        for donor, marker in markers.items()
    ]
    fig.legend(family_handles + donor_handles, list(palette) + list(markers), frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.015))
    fig.tight_layout(rect=(0, 0.08, 1, 1), w_pad=1.9)
    save(fig, figure_dir / "66-validation-boundary")


def plot_pathways(input_dir: Path, figure_dir: Path) -> None:
    counts = read(input_dir / "pathway-counts.tsv")
    within = counts[counts["Comparison"].eq("Within donor · separate approaches")]
    categories = ["Common", "Only reference-guided", "Only MAG-guided"]
    colors = [SHARED, REFERENCE, MAG]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.7), gridspec_kw={"width_ratios": [1.05, 1.25]})
    ax = axes[0]
    donors = ["Donor 1", "Donor 2"]
    bottom = np.zeros(2)
    for category, color in zip(categories, colors, strict=True):
        values = []
        for donor in donors:
            selected = within[(within["Parent"].eq(donor)) & (within["Category"].eq(category))]
            values.append(int(selected["Count"].iloc[0]) if len(selected) else 0)
        bars = ax.bar(donors, values, bottom=bottom, color=color, edgecolor="white", label=category.replace("Only ", ""))
        for bar, value, base in zip(bars, values, bottom, strict=True):
            if value:
                ax.text(bar.get_x() + bar.get_width() / 2, base + value / 2, str(value), ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
        bottom += values
    ax.set_ylabel("Enriched pathways (q < 0.05)")
    ax.set_xlabel("")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    panel(ax, "A", "Within-donor calls depend on genome input")

    ax = axes[1]
    approaches = ["Reference-guided", "MAG-guided", "Combined"]
    donor1 = [5, 1, 2]
    donor2 = [7, 5, 3]
    x = np.arange(3)
    width = 0.34
    ax.bar(x - width / 2, donor1, width, color=DONOR_COLORS["Donor 1"], label="Donor 1-specific")
    ax.bar(x + width / 2, donor2, width, color=DONOR_COLORS["Donor 2"], label="Donor 2-specific")
    ax.set_xticks(x, approaches)
    ax.set_ylabel("Donor-specific enriched pathways")
    ax.set_ylim(0, 8.6)
    ax.legend(frameon=False, loc="upper center")
    for xpos, value in zip(x - width / 2, donor1, strict=True):
        ax.text(xpos, value + 0.15, str(value), ha="center", fontsize=8.5)
    for xpos, value in zip(x + width / 2, donor2, strict=True):
        ax.text(xpos, value + 0.15, str(value), ha="center", fontsize=8.5)
    ax.text(
        0.03,
        0.96,
        "Only propanoate metabolism\nremained Donor 2-specific\nin all three constructions.",
        transform=ax.transAxes,
        va="top",
        fontsize=8.2,
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#D0D7DE"},
    )
    panel(ax, "B", "Between-donor phenotype stories are construction-sensitive")
    fig.tight_layout(w_pad=2.0)
    save(fig, figure_dir / "66-pathway-robustness")


def main() -> None:
    args = parse_args()
    configure()
    np.random.seed(PLOT_SEED)
    input_dir = args.input_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    plot_quality(input_dir, figure_dir)
    plot_normalization(input_dir, figure_dir)
    plot_gapfill(input_dir, figure_dir)
    plot_overlap(input_dir, figure_dir)
    plot_validation(input_dir, figure_dir)
    plot_pathways(input_dir, figure_dir)
    anchor = input_dir / "majzoub-fig2-original.jpg"
    if not anchor.is_file():
        raise FileNotFoundError(anchor)
    shutil.copy2(anchor, figure_dir / "66-majzoub-fig2-original.jpg")
    print(f"Rendered {len(FIGURES)} Article 66 figure sets in {figure_dir}")


if __name__ == "__main__":
    main()
