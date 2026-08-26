#!/usr/bin/env python3
"""Create publication-ready Article 64 figures from prepared evidence tables."""

from __future__ import annotations

import argparse
import json
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PLOT_SEED = 20_260_764
PALETTE = {
    "Control": "#009E73",
    "CD": "#0072B2",
    "UC": "#D55E00",
    "DNA": "#0072B2",
    "RNA": "#CC79A7",
    "Higher RNA share": "#D55E00",
    "Higher DNA share": "#0072B2",
    "Not significant": "#9AA0A6",
}
FIGURES = (
    "64-pairing-audit",
    "64-dna-rna-concordance",
    "64-sample-concordance",
    "64-relative-activity",
    "64-diagnosis-audit",
    "64-ratio-sensitivity",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
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


def clipped_label(identifier: str, description: str, width: int = 45) -> str:
    clean = " ".join(str(description).split())
    return f"{identifier} · {textwrap.shorten(clean, width=width, placeholder='…')}"


def plot_pairing_audit(input_dir: Path, figure_dir: Path) -> None:
    attrition = read_tsv(input_dir / "sample-attrition.tsv")
    metadata = read_tsv(input_dir / "sample-metadata.tsv")
    metrics = json.loads((input_dir / "analysis-metrics.json").read_text(encoding="utf-8"))

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.5), gridspec_kw={"width_ratios": [1.55, 1]})
    ax = axes[0]
    bars = ax.barh(
        attrition["Stage"],
        attrition["Count"],
        color=["#4C78A8", "#B279A2", "#72B7B2", "#59A14F", "#F28E2B", "#E15759"],
        edgecolor="white",
    )
    ax.invert_yaxis()
    ax.set_xlabel("Count")
    ax.set_ylabel("")
    ax.grid(axis="x", alpha=0.25)
    ax.grid(axis="y", visible=False)
    for bar, unit in zip(bars, attrition["Unit"]):
        ax.text(
            bar.get_width() + max(attrition["Count"]) * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"{int(bar.get_width()):,} {unit}",
            va="center",
            fontsize=8,
        )
    ax.set_xlim(0, max(attrition["Count"]) * 1.23)
    panel(ax, "A", "From assay profiles to independent subjects")

    ax = axes[1]
    subjects = metadata.drop_duplicates("SubjectID")["Diagnosis"].value_counts().reindex(["Control", "CD", "UC"])
    bars = ax.bar(subjects.index, subjects.values, color=[PALETTE[x] for x in subjects.index], width=0.65)
    ax.set_ylabel("Independent subjects")
    ax.set_xlabel("")
    ax.grid(axis="y", alpha=0.25)
    ax.grid(axis="x", visible=False)
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8, f"{int(bar.get_height())}", ha="center")
    ax.set_ylim(0, max(subjects.values) * 1.25)
    ax.text(
        0.02,
        0.96,
        f"Excluded before analysis\n{metrics['technical_replicates_excluded']} paired technical replicates\n{metrics['zero_layer_samples_excluded']} zero-layer pairs",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#D0D7DE"},
    )
    panel(ax, "B", "Diagnosis balance after subject de-duplication")
    fig.tight_layout(w_pad=2.0)
    save(fig, figure_dir / "64-pairing-audit")


def plot_dna_rna(input_dir: Path, figure_dir: Path) -> None:
    audit = read_tsv(input_dir / "feature-audit.tsv")
    activity = read_tsv(input_dir / "activity-results.tsv")
    data = audit.loc[audit["Selected"]].merge(
        activity[["Feature", "MedianLog2RNAoverDNA"]], on="Feature", validate="one_to_one"
    )
    x = np.log10(data["MeanDNA"] + 1e-8)
    y = np.log10(data["MeanRNA"] + 1e-8)
    color = data["MedianLog2RNAoverDNA"]
    limit = float(np.nanmax(np.abs(color)))

    fig, ax = plt.subplots(figsize=(7.1, 5.8))
    points = ax.scatter(
        x,
        y,
        c=color,
        cmap="coolwarm",
        vmin=-limit,
        vmax=limit,
        s=34,
        alpha=0.82,
        edgecolor="white",
        linewidth=0.35,
    )
    lo = min(x.min(), y.min()) - 0.15
    hi = max(x.max(), y.max()) + 0.15
    ax.plot([lo, hi], [lo, hi], color="#37474F", linestyle="--", linewidth=1)
    extremes = pd.concat([data.nsmallest(3, "MedianLog2RNAoverDNA"), data.nlargest(3, "MedianLog2RNAoverDNA")])
    offsets = ((5, 5), (5, 12), (5, -14), (5, 5), (5, 12), (5, -12))
    for (_, row), offset in zip(extremes.iterrows(), offsets):
        i = data.index[data["Feature"].eq(row["Feature"])][0]
        ax.annotate(
            row["PathwayID"],
            (x.loc[i], y.loc[i]),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.5,
        )
    colorbar = fig.colorbar(points, ax=ax, pad=0.02)
    colorbar.set_label("Median log2 relative RNA / DNA")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Mean DNA pathway relative abundance (log10)")
    ax.set_ylabel("Mean RNA pathway relative abundance (log10)")
    ax.grid(alpha=0.22)
    panel(ax, "A", "Pathway potential and transcript share are related but not identical")
    ax.text(
        0.02,
        0.98,
        f"{len(data)} pathways · co-detection ≥20%",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#D0D7DE"},
    )
    fig.tight_layout()
    save(fig, figure_dir / "64-dna-rna-concordance")


def plot_sample_concordance(input_dir: Path, figure_dir: Path) -> None:
    data = read_tsv(input_dir / "subject-concordance.tsv")
    summary = read_tsv(input_dir / "concordance-summary.tsv")
    order = ["Control", "CD", "UC"]
    rng = np.random.default_rng(PLOT_SEED)

    fig, ax = plt.subplots(figsize=(6.9, 5.0))
    sns.violinplot(
        data=data,
        x="Diagnosis",
        y="MedianSpearmanRho",
        order=order,
        palette=PALETTE,
        inner=None,
        cut=0,
        linewidth=0.8,
        alpha=0.35,
        ax=ax,
        hue="Diagnosis",
        legend=False,
    )
    sns.boxplot(
        data=data,
        x="Diagnosis",
        y="MedianSpearmanRho",
        order=order,
        width=0.24,
        showfliers=False,
        boxprops={"facecolor": "white", "alpha": 0.85},
        medianprops={"color": "black", "linewidth": 1.3},
        ax=ax,
    )
    for position, group in enumerate(order):
        values = data.loc[data["Diagnosis"].eq(group), "MedianSpearmanRho"].to_numpy()
        jitter = rng.uniform(-0.12, 0.12, size=len(values))
        ax.scatter(position + jitter, values, s=17, color=PALETTE[group], alpha=0.72, edgecolor="white", linewidth=0.25)
    overall = summary.loc[summary["Diagnosis"].eq("All")].iloc[0]
    ax.axhline(overall["MedianSpearmanRho"], color="#37474F", linestyle="--", linewidth=1)
    ax.text(
        0.02,
        0.04,
        f"Overall subject median = {overall['MedianSpearmanRho']:.3f}\n95% subject-bootstrap CI {overall['CILow']:.3f}–{overall['CIHigh']:.3f}",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#D0D7DE"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("Within-sample DNA–RNA Spearman rho\n(subject median)")
    ax.grid(axis="y", alpha=0.23)
    ax.grid(axis="x", visible=False)
    panel(ax, "A", "DNA and RNA pathway profiles remain concordant within subjects")
    fig.tight_layout()
    save(fig, figure_dir / "64-sample-concordance")


def plot_activity(input_dir: Path, figure_dir: Path) -> None:
    data = read_tsv(input_dir / "activity-results.tsv")
    display = pd.concat(
        [data.nsmallest(6, "MedianLog2RNAoverDNA"), data.nlargest(6, "MedianLog2RNAoverDNA")]
    ).sort_values("MedianLog2RNAoverDNA")
    display = display.copy()
    display["Label"] = [
        clipped_label(identifier, description) for identifier, description in zip(display["PathwayID"], display["Description"])
    ]
    y = np.arange(len(display))
    colors = [
        PALETTE["Higher RNA share"] if value > 0 else PALETTE["Higher DNA share"]
        for value in display["MedianLog2RNAoverDNA"]
    ]

    fig, ax = plt.subplots(figsize=(9.3, 5.9))
    ax.hlines(y, display["CILow"], display["CIHigh"], color=colors, linewidth=2)
    ax.scatter(display["MedianLog2RNAoverDNA"], y, color=colors, s=48, edgecolor="white", linewidth=0.6, zorder=3)
    ax.axvline(0, color="#37474F", linewidth=1, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(display["Label"])
    ax.set_xlabel("Median log2 relative RNA / DNA (95% subject-bootstrap CI)")
    ax.set_ylabel("")
    ax.grid(axis="x", alpha=0.23)
    ax.grid(axis="y", visible=False)
    panel(ax, "A", "Pathways with the largest relative transcriptional allocation shifts")
    ax.text(
        0.01,
        0.01,
        "Ratios use co-detected pairs only; they are relative allocation, not absolute transcription rates.",
        transform=ax.transAxes,
        fontsize=8,
        color="#455A64",
    )
    fig.tight_layout()
    save(fig, figure_dir / "64-relative-activity")


def plot_diagnosis_audit(input_dir: Path, figure_dir: Path) -> None:
    data = read_tsv(input_dir / "diagnosis-results.tsv").dropna(subset=["PValue", "QValue"]).copy()
    data["MinusLog10P"] = -np.log10(data["PValue"].clip(lower=np.finfo(float).tiny))
    top = data.nsmallest(6, "PValue")

    fig, ax = plt.subplots(figsize=(7.4, 5.5))
    points = ax.scatter(
        data["MedianRange"],
        data["MinusLog10P"],
        c=data["QValue"],
        cmap="viridis_r",
        vmin=0,
        vmax=1,
        s=32,
        alpha=0.8,
        edgecolor="white",
        linewidth=0.35,
    )
    ax.axhline(-np.log10(0.05), color="#37474F", linestyle="--", linewidth=1, label="Raw P = 0.05")
    for _, row in top.iterrows():
        ax.annotate(row["PathwayID"], (row["MedianRange"], row["MinusLog10P"]), xytext=(4, 4), textcoords="offset points", fontsize=7.2)
    cbar = fig.colorbar(points, ax=ax, pad=0.02)
    cbar.set_label("BH q-value")
    ax.set_xlabel("Range of diagnosis-specific median log2 RNA / DNA")
    ax.set_ylabel("−log10 raw Kruskal–Wallis P")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(alpha=0.22)
    panel(ax, "A", "Nominal diagnosis contrasts do not survive pathway-wide FDR")
    ax.text(
        0.02,
        0.97,
        f"0 / {len(data)} pathways at BH q < 0.05\nminimum BH q = {data['QValue'].min():.3f}",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#D0D7DE"},
    )
    fig.tight_layout()
    save(fig, figure_dir / "64-diagnosis-audit")


def plot_sensitivity(input_dir: Path, figure_dir: Path) -> None:
    data = read_tsv(input_dir / "sensitivity-summary.tsv")
    gates = data.loc[data["Analysis"].eq("Co-detection gate")]
    pseudo = data.loc[data["Analysis"].eq("Pseudocount")]

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.4))
    ax = axes[0]
    bars = ax.bar(gates["Setting"], gates["Features"], color=["#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00"])
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 4, str(int(bar.get_height())), ha="center", fontsize=8)
    ax.set_ylim(0, gates["Features"].max() * 1.18)
    ax.set_xlabel("Required DNA–RNA co-detection")
    ax.set_ylabel("Retained pathways")
    ax.grid(axis="y", alpha=0.22)
    ax.grid(axis="x", visible=False)
    panel(ax, "A", "Feature count depends on the detection gate")

    ax = axes[1]
    bars = ax.bar(pseudo["Setting"], pseudo["RankSpearman"], color=["#80B1D3", "#4C78A8", "#2F4B7C"])
    for bar, shift in zip(bars, pseudo["MedianAbsoluteShift"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.015,
            f"rho {bar.get_height():.2f}\nshift {shift:.2f}",
            ha="center",
            fontsize=7.5,
        )
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Pseudocount in all-pair ratio sensitivity")
    ax.set_ylabel("Rank correlation with co-detected primary analysis")
    ax.grid(axis="y", alpha=0.22)
    ax.grid(axis="x", visible=False)
    panel(ax, "B", "Zero handling can reorder pathway activity")
    fig.tight_layout(w_pad=2.0)
    save(fig, figure_dir / "64-ratio-sensitivity")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    configure()
    plot_pairing_audit(input_dir, figure_dir)
    plot_dna_rna(input_dir, figure_dir)
    plot_sample_concordance(input_dir, figure_dir)
    plot_activity(input_dir, figure_dir)
    plot_diagnosis_audit(input_dir, figure_dir)
    plot_sensitivity(input_dir, figure_dir)
    missing = [stem for stem in FIGURES if not (figure_dir / f"{stem}.png").is_file()]
    if missing:
        raise RuntimeError(f"Missing figures: {missing}")
    print(f"Created {len(FIGURES)} Article 64 figure sets in {figure_dir}")


if __name__ == "__main__":
    main()
