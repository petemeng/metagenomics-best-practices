#!/usr/bin/env python3
"""Create publication-ready Article 65 figures from frozen evidence tables."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PLOT_SEED = 20_260_765
PAIR_ORDER = ("DNA–RNA", "DNA–Protein", "RNA–Protein")
PAIR_COLORS = {
    "DNA–RNA": "#0072B2",
    "DNA–Protein": "#E69F00",
    "RNA–Protein": "#CC79A7",
}
THRESHOLD_ORDER = (
    "1 peptide · 1% FDR",
    "1 peptide · 5% FDR",
    "2 peptides · 1% FDR",
    "2 peptides · 5% FDR",
)
THRESHOLD_LABELS = {
    "1 peptide · 1% FDR": "≥1 peptide\n1% protein FDR",
    "1 peptide · 5% FDR": "≥1 peptide\n5% protein FDR",
    "2 peptides · 1% FDR": "≥2 peptides\n1% protein FDR",
    "2 peptides · 5% FDR": "≥2 peptides\n5% protein FDR",
}
THRESHOLD_COLORS = ("#4C78A8", "#72B7B2", "#F28E2B", "#E15759")
FIGURES = (
    "65-evidence-threshold-audit",
    "65-protein-namespace-audit",
    "65-multiomic-sample-alignment",
    "65-three-layer-concordance",
    "65-ec-cross-subject-correlation",
    "65-threshold-richness-sensitivity",
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


def clipped_ec(row: pd.Series, width: int = 42) -> str:
    description = " ".join(str(row["Description"]).split())
    return f"{row['EC']} · {textwrap.shorten(description, width=width, placeholder='…')}"


def plot_threshold_audit(input_dir: Path, figure_dir: Path) -> None:
    audit = read_tsv(input_dir / "threshold-audit.tsv").set_index("Threshold").loc[list(THRESHOLD_ORDER)].reset_index()
    richness = read_tsv(input_dir / "sample-protein-richness.tsv")
    long = richness.melt(id_vars="SampleID", var_name="Threshold", value_name="Detected proteins").dropna()
    long["Display"] = long["Threshold"].map(THRESHOLD_LABELS)
    display_order = [THRESHOLD_LABELS[x] for x in THRESHOLD_ORDER]

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8), gridspec_kw={"width_ratios": [1.05, 1.35]})
    ax = axes[0]
    positions = np.arange(len(audit))
    bars = ax.barh(positions, audit["DetectedProteinIDs"], color=THRESHOLD_COLORS, edgecolor="white")
    ax.set_yticks(positions, display_order)
    ax.invert_yaxis()
    ax.set_xlabel("Protein IDs detected in ≥1 profile")
    ax.set_ylabel("")
    ax.grid(axis="x", alpha=0.25)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, audit["DetectedProteinIDs"].max() * 1.18)
    for bar, (_, row) in zip(bars, audit.iterrows()):
        ax.text(
            bar.get_width() + audit["DetectedProteinIDs"].max() * 0.02,
            bar.get_y() + bar.get_height() / 2,
            f"{int(row['DetectedProteinIDs']):,}",
            va="center",
            fontsize=8.5,
        )
    panel(ax, "A", "Protein evidence gates change the reported universe")

    ax = axes[1]
    sns.boxplot(
        data=long,
        x="Display",
        y="Detected proteins",
        order=display_order,
        hue="Display",
        palette=dict(zip(display_order, THRESHOLD_COLORS)),
        legend=False,
        showfliers=False,
        width=0.62,
        ax=ax,
    )
    audit_by_threshold = audit.set_index("Threshold")
    for i, threshold in enumerate(THRESHOLD_ORDER):
        value = audit_by_threshold.loc[threshold, "MedianDetectedPerSample"]
        label_y = audit_by_threshold.loc[threshold, "IQRHighDetectedPerSample"] + 70
        ax.text(
            i,
            label_y,
            f"median {value:,.1f}",
            ha="center",
            va="bottom",
            fontsize=7.8,
            bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
        )
    ax.set_xlabel("")
    ax.set_ylabel("Detected protein IDs per profile")
    ax.tick_params(axis="x", rotation=0)
    ax.grid(axis="y", alpha=0.25)
    ax.grid(axis="x", visible=False)
    panel(ax, "B", "Per-sample coverage is threshold-dependent")
    fig.tight_layout(w_pad=2.2)
    save(fig, figure_dir / "65-evidence-threshold-audit")


def plot_namespace_audit(input_dir: Path, figure_dir: Path) -> None:
    audit = read_tsv(input_dir / "protein-namespace-audit.tsv")
    order = [
        "Generic accession",
        "Taxon-prefixed reference",
        "Explicit host namespace",
        "Contaminant namespace",
    ]
    audit = audit.set_index("Namespace").loc[order].reset_index()
    labels = [value.replace(" namespace", "\nnamespace").replace(" reference", "\nreference") for value in order]
    colors = ["#0072B2", "#009E73", "#CC79A7", "#D55E00"]

    fig, axes = plt.subplots(1, 2, figsize=(10.9, 4.7), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    x = np.arange(len(audit))
    width = 0.36
    ax.bar(x - width / 2, audit["ProteinIDs"], width, color="#B0BEC5", label="Rows in table")
    ax.bar(x + width / 2, audit["DetectedProteinIDs"], width, color=colors, label="Detected IDs")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Protein IDs (symlog scale)")
    ax.set_xlabel("")
    ax.grid(axis="y", alpha=0.23)
    ax.grid(axis="x", visible=False)
    ax.legend(frameon=False, loc="upper right")
    for i, row in audit.iterrows():
        ax.text(i - width / 2, row["ProteinIDs"] * 1.22 if row["ProteinIDs"] else 0.3, f"{int(row['ProteinIDs']):,}", ha="center", fontsize=7.5)
        detected_y = row["DetectedProteinIDs"] * 1.22 if row["DetectedProteinIDs"] else 0.35
        ax.text(i + width / 2, detected_y, f"{int(row['DetectedProteinIDs']):,}", ha="center", fontsize=7.5)
    panel(ax, "A", "Identifier syntax is not a taxonomic classifier")

    ax = axes[1]
    shares = audit["ReportedCountShare"] * 100
    bars = ax.barh(labels, shares, color=colors, edgecolor="white")
    ax.invert_yaxis()
    ax.set_xlabel("Share of reported spectral counts (%)")
    ax.set_ylabel("")
    ax.grid(axis="x", alpha=0.25)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 104)
    for bar, share in zip(bars, shares):
        x_text = share + 1 if share > 0 else 0.8
        ax.text(x_text, bar.get_y() + bar.get_height() / 2, f"{share:.1f}%", va="center", fontsize=8)
    ax.text(
        0.33,
        0.36,
        "Rows with zero reported counts\nremain database candidates, not detections.",
        transform=ax.transAxes,
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#D0D7DE"},
    )
    panel(ax, "B", "Detection evidence must accompany the namespace")
    fig.tight_layout(w_pad=2.0)
    save(fig, figure_dir / "65-protein-namespace-audit")


def plot_sample_alignment(input_dir: Path, figure_dir: Path) -> None:
    attrition = read_tsv(input_dir / "sample-attrition.tsv").set_index("Stage")
    metadata = read_tsv(input_dir / "sample-metadata.tsv")
    counts = {stage: int(attrition.loc[stage, "Count"]) for stage in attrition.index}

    fig, axes = plt.subplots(1, 2, figsize=(11.3, 4.8), gridspec_kw={"width_ratios": [1.65, 1]})
    ax = axes[0]
    ax.axis("off")
    nodes = {
        "MPX profiles": (0.08, 0.72),
        "Exact MPX + MGX": (0.33, 0.72),
        "Exact MPX + MGX + MTX": (0.60, 0.72),
        "Three-layer EC-complete": (0.87, 0.72),
        "Triple profiles with MBX product": (0.60, 0.26),
        "Four-layer & EC-complete": (0.87, 0.26),
    }
    short = {
        "MPX profiles": "MPX",
        "Exact MPX + MGX": "MPX + MGX",
        "Exact MPX + MGX + MTX": "MPX + MGX + MTX",
        "Three-layer EC-complete": "EC-complete triple",
        "Triple profiles with MBX product": "+ MBX product",
        "Four-layer & EC-complete": "Four-layer complete",
    }
    for stage, (x, y) in nodes.items():
        ax.text(
            x,
            y,
            f"{short[stage]}\n{counts[stage]:,} samples",
            ha="center",
            va="center",
            fontsize=8.5,
            bbox={"boxstyle": "round,pad=0.45", "facecolor": "#F7FAFC", "edgecolor": "#607D8B", "linewidth": 1.1},
            transform=ax.transAxes,
        )
    arrows = [
        ("MPX profiles", "Exact MPX + MGX"),
        ("Exact MPX + MGX", "Exact MPX + MGX + MTX"),
        ("Exact MPX + MGX + MTX", "Three-layer EC-complete"),
        ("Exact MPX + MGX + MTX", "Triple profiles with MBX product"),
        ("Triple profiles with MBX product", "Four-layer & EC-complete"),
    ]
    for first, second in arrows:
        x1, y1 = nodes[first]
        x2, y2 = nodes[second]
        ax.annotate(
            "",
            xy=(x2 - (0.09 if x2 > x1 else 0), y2),
            xytext=(x1 + (0.09 if x2 > x1 else 0), y1),
            xycoords=ax.transAxes,
            textcoords=ax.transAxes,
            arrowprops={"arrowstyle": "->", "color": "#455A64", "lw": 1.2, "connectionstyle": "arc3,rad=0"},
        )
    ax.text(0.43, 0.47, "exact biospecimen ID", transform=ax.transAxes, ha="center", fontsize=7.5, color="#455A64")
    panel(ax, "A", "Four assays align only after an explicit sample ledger")

    ax = axes[1]
    subjects = metadata.drop_duplicates("SubjectID")["Diagnosis"].value_counts().reindex(["Control", "CD", "UC"])
    colors = ["#009E73", "#0072B2", "#D55E00"]
    bars = ax.bar(subjects.index, subjects.values, color=colors, width=0.64)
    ax.set_ylabel("Independent subjects")
    ax.set_xlabel("")
    ax.grid(axis="y", alpha=0.25)
    ax.grid(axis="x", visible=False)
    ax.set_ylim(0, max(subjects) * 1.26)
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.7, f"{int(bar.get_height())}", ha="center")
    ax.text(
        0.04,
        0.95,
        "186 EC-complete samples\n76 independent subjects",
        transform=ax.transAxes,
        va="top",
        fontsize=8.3,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#D0D7DE"},
    )
    panel(ax, "B", "Repeated profiles collapse to the participant level")
    fig.tight_layout(w_pad=1.6)
    save(fig, figure_dir / "65-multiomic-sample-alignment")


def plot_three_layer_concordance(input_dir: Path, figure_dir: Path) -> None:
    data = read_tsv(input_dir / "subject-concordance.tsv")
    summary = read_tsv(input_dir / "concordance-summary.tsv").set_index("LayerPair").loc[list(PAIR_ORDER)].reset_index()
    rng = np.random.default_rng(PLOT_SEED)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.9), gridspec_kw={"width_ratios": [1.45, 1]})
    ax = axes[0]
    sns.violinplot(
        data=data,
        x="LayerPair",
        y="MedianSpearmanRho",
        order=list(PAIR_ORDER),
        hue="LayerPair",
        palette=PAIR_COLORS,
        legend=False,
        inner=None,
        cut=0,
        linewidth=0.8,
        ax=ax,
    )
    sns.boxplot(
        data=data,
        x="LayerPair",
        y="MedianSpearmanRho",
        order=list(PAIR_ORDER),
        width=0.24,
        showfliers=False,
        boxprops={"facecolor": "white", "alpha": 0.9},
        medianprops={"color": "black", "linewidth": 1.3},
        ax=ax,
    )
    for position, pair in enumerate(PAIR_ORDER):
        values = data.loc[data["LayerPair"].eq(pair), "MedianSpearmanRho"].to_numpy(float)
        jitter = rng.uniform(-0.12, 0.12, size=len(values))
        ax.scatter(position + jitter, values, s=15, color=PAIR_COLORS[pair], alpha=0.62, edgecolor="white", linewidth=0.2)
    ax.axhline(0, color="#455A64", linestyle="--", linewidth=0.9)
    ax.set_xlabel("")
    ax.set_ylabel("Within-sample EC Spearman rho\n(subject median)")
    ax.grid(axis="y", alpha=0.23)
    ax.grid(axis="x", visible=False)
    panel(ax, "A", "RNA resembles protein more than DNA does")

    ax = axes[1]
    y = np.arange(len(summary))[::-1]
    for position, (_, row) in zip(y, summary.iterrows()):
        pair = row["LayerPair"]
        ax.errorbar(
            row["MedianSpearmanRho"],
            position,
            xerr=[[row["MedianSpearmanRho"] - row["CILow"]], [row["CIHigh"] - row["MedianSpearmanRho"]]],
            fmt="o",
            color=PAIR_COLORS[pair],
            capsize=3,
            markersize=6,
        )
        ax.text(row["CIHigh"] + 0.025, position, f"{row['MedianSpearmanRho']:.3f}", va="center", fontsize=8.5)
    ax.set_yticks(y, summary["LayerPair"])
    ax.axvline(0, color="#455A64", linestyle="--", linewidth=0.9)
    ax.set_xlim(-0.05, 0.86)
    ax.set_xlabel("Median rho (95% subject-bootstrap CI)")
    ax.set_ylabel("")
    ax.grid(axis="x", alpha=0.23)
    ax.grid(axis="y", visible=False)
    panel(ax, "B", "Effect sizes retain the independent sampling unit")
    fig.tight_layout(w_pad=2.0)
    save(fig, figure_dir / "65-three-layer-concordance")


def plot_ec_correlations(input_dir: Path, figure_dir: Path) -> None:
    data = read_tsv(input_dir / "ec-correlations.tsv")
    data["BH05"] = data["BH05"].astype(str).str.lower().isin(["true", "1", "t", "yes"])
    top = data.loc[data["LayerPair"].eq("RNA–Protein")].nlargest(10, "SpearmanRho").copy()
    top["Label"] = top.apply(clipped_ec, axis=1)
    top = top.sort_values("SpearmanRho")

    fig, axes = plt.subplots(1, 2, figsize=(12.1, 5.1), gridspec_kw={"width_ratios": [1, 1.55]})
    ax = axes[0]
    sns.violinplot(
        data=data,
        x="LayerPair",
        y="SpearmanRho",
        order=list(PAIR_ORDER),
        hue="LayerPair",
        palette=PAIR_COLORS,
        legend=False,
        inner="box",
        cut=0,
        linewidth=0.8,
        ax=ax,
    )
    ax.axhline(0, color="#455A64", linestyle="--", linewidth=0.9)
    ax.set_xlabel("")
    ax.set_ylabel("Cross-subject EC Spearman rho")
    ax.grid(axis="y", alpha=0.23)
    ax.grid(axis="x", visible=False)
    for position, pair in enumerate(PAIR_ORDER):
        section = data.loc[data["LayerPair"].eq(pair)]
        ax.text(
            position,
            0.94,
            f"median {section['SpearmanRho'].median():.3f}\n{int(section['BH05'].sum())}/263 BH q<0.05",
            ha="center",
            va="top",
            fontsize=7.5,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "#D0D7DE"},
        )
    ax.set_ylim(min(-0.55, data["SpearmanRho"].min() - 0.05), 1.0)
    panel(ax, "A", "Functional concordance varies by molecular layer")

    ax = axes[1]
    y = np.arange(len(top))
    ax.hlines(y, 0, top["SpearmanRho"], color="#CC79A7", alpha=0.7, linewidth=1.5)
    ax.scatter(top["SpearmanRho"], y, color="#CC79A7", s=38, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_yticks(y, top["Label"])
    ax.axvline(0, color="#455A64", linestyle="--", linewidth=0.9)
    ax.set_xlabel("RNA–protein Spearman rho across 76 subjects")
    ax.set_ylabel("")
    ax.grid(axis="x", alpha=0.23)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, max(0.62, top["SpearmanRho"].max() + 0.08))
    for position, (_, row) in zip(y, top.iterrows()):
        ax.text(row["SpearmanRho"] + 0.012, position, f"q={row['QValue']:.2g}", va="center", fontsize=7.2)
    panel(ax, "B", "Top RNA–protein EC associations remain hypotheses")
    fig.tight_layout(w_pad=2.0)
    save(fig, figure_dir / "65-ec-cross-subject-correlation")


def plot_threshold_sensitivity(input_dir: Path, figure_dir: Path) -> None:
    richness = read_tsv(input_dir / "sample-protein-richness.tsv")
    correlations = read_tsv(input_dir / "richness-correlations.tsv")
    first = "1 peptide · 5% FDR"
    second = "2 peptides · 1% FDR"
    paired = richness[[first, second]].dropna()
    pair_rho = correlations.loc[
        ((correlations["ThresholdA"].eq(first)) & (correlations["ThresholdB"].eq(second)))
        | ((correlations["ThresholdA"].eq(second)) & (correlations["ThresholdB"].eq(first))),
        "SpearmanRho",
    ].iloc[0]

    matrix = pd.DataFrame(np.eye(len(THRESHOLD_ORDER)), index=THRESHOLD_ORDER, columns=THRESHOLD_ORDER)
    for _, row in correlations.iterrows():
        matrix.loc[row["ThresholdA"], row["ThresholdB"]] = row["SpearmanRho"]
        matrix.loc[row["ThresholdB"], row["ThresholdA"]] = row["SpearmanRho"]
    short = ["1 pep / 1%", "1 pep / 5%", "2 pep / 1%", "2 pep / 5%"]
    matrix.index = short
    matrix.columns = short

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    ax.scatter(paired[first], paired[second], s=20, alpha=0.58, color="#0072B2", edgecolor="white", linewidth=0.25)
    lo = min(paired[first].min(), paired[second].min())
    hi = max(paired[first].max(), paired[second].max())
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="#455A64", linewidth=1)
    ax.set_xlabel("Detected proteins: ≥1 peptide, 5% FDR")
    ax.set_ylabel("Detected proteins: ≥2 peptides, 1% FDR")
    ax.text(
        0.04,
        0.94,
        f"Shared profiles = {len(paired)}\nSpearman rho = {pair_rho:.3f}",
        transform=ax.transAxes,
        va="top",
        fontsize=8.3,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#D0D7DE"},
    )
    ax.grid(alpha=0.23)
    panel(ax, "A", "Sample rankings can be stable while counts shift")

    ax = axes[1]
    sns.heatmap(
        matrix,
        annot=True,
        fmt=".3f",
        cmap="YlGnBu",
        vmin=0.95,
        vmax=1.0,
        square=True,
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": "Spearman rho", "shrink": 0.78},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=35)
    ax.tick_params(axis="y", rotation=0)
    panel(ax, "B", "All four richness rankings are highly concordant")
    fig.tight_layout(w_pad=2.2)
    save(fig, figure_dir / "65-threshold-richness-sensitivity")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    configure()
    plot_threshold_audit(input_dir, figure_dir)
    plot_namespace_audit(input_dir, figure_dir)
    plot_sample_alignment(input_dir, figure_dir)
    plot_three_layer_concordance(input_dir, figure_dir)
    plot_ec_correlations(input_dir, figure_dir)
    plot_threshold_sensitivity(input_dir, figure_dir)
    print(f"Created {len(FIGURES)} Article 65 figure sets in {figure_dir}")


if __name__ == "__main__":
    main()
