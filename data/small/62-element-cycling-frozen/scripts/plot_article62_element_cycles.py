#!/usr/bin/env python3
"""Create publication-ready, English-only figures for Article 62."""

from __future__ import annotations

import argparse
import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/article62-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D


BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
YELLOW = "#E69F00"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
GRAY = "#6B7478"
LIGHT = "#EDF2F4"
DARK = "#263238"
WHITE = "#FFFFFF"
RED = "#B2182B"
ELEMENT_COLORS = {
    "Carbon": GREEN,
    "Carbon/Nitrogen": PURPLE,
    "Nitrogen": BLUE,
    "Sulfur": ORANGE,
}
REGIME_ORDER = ["below_30", "30_50", "50_70", "70_80", "80_100"]
REGIME_LABELS = {
    "below_30": "<30 °C",
    "30_50": "30–<50 °C",
    "50_70": "50–<70 °C",
    "70_80": "70–<80 °C",
    "80_100": "80–100 °C",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    if frame.empty:
        raise ValueError(f"Empty plotting table: {path}")
    return frame


def truth(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    values = series.astype(str).str.lower()
    if not values.isin(["true", "false"]).all():
        raise ValueError(f"Unexpected logical values: {sorted(values.unique())}")
    return values.eq("true")


def wrap(value: str, width: int = 28) -> str:
    return "\n".join(
        textwrap.wrap(str(value), width=width, break_long_words=False)
    )


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.5,
            "axes.edgecolor": "#455A64",
            "axes.linewidth": 0.7,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "figure.facecolor": WHITE,
            "axes.facecolor": WHITE,
            "savefig.facecolor": WHITE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_pub(fig: plt.Figure, directory: Path, stem: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(directory / f"{stem}.png", dpi=350, bbox_inches="tight")
    fig.savefig(
        directory / f"{stem}.tiff",
        dpi=350,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def rule_order(analysis: Path) -> tuple[pd.DataFrame, list[str], dict[str, str]]:
    rules = read_tsv(analysis / "process-rules.tsv")
    order = rules.ProcessID.tolist()
    labels = dict(zip(rules.ProcessID, rules.Process))
    return rules, order, labels


def plot_landscape(analysis: Path, figures: Path) -> None:
    rules, order, labels = rule_order(analysis)
    regime = read_tsv(analysis / "temperature-regime-summary.tsv")
    matrix = regime.pivot_table(
        index="ProcessID", columns="temperature_regime", values="MedianIndex"
    ).reindex(index=order, columns=REGIME_ORDER)

    standardized = matrix.astype(float).copy()
    for process_id, values in matrix.iterrows():
        positive = values[values > 0]
        floor = float(positive.min() / 2) if len(positive) else 1e-12
        logged = np.log10(values.to_numpy(float) + floor)
        spread = float(np.std(logged, ddof=0))
        standardized.loc[process_id] = (
            logged - float(np.mean(logged))
        ) / (spread if spread > 0 else 1.0)
    standardized.index = [wrap(labels[value], 25) for value in standardized.index]
    standardized.columns = [REGIME_LABELS[value] for value in standardized.columns]

    fig, ax = plt.subplots(figsize=(9.2, 9.4))
    cmap = LinearSegmentedColormap.from_list("within_process", [BLUE, WHITE, RED])
    sns.heatmap(
        standardized,
        cmap=cmap,
        center=0,
        linewidths=0.5,
        linecolor=WHITE,
        cbar_kws={"label": "Within-process z score of log10 median KO index"},
        ax=ax,
    )
    ax.set_xlabel("Sample temperature regime")
    ax.set_ylabel("Marker-defined biogeochemical process")
    ax.set_title("Community marker profiles shift across temperature regimes")
    ax.tick_params(axis="x", rotation=0)
    element_boundaries = np.flatnonzero(
        rules.Element.to_numpy()[1:] != rules.Element.to_numpy()[:-1]
    ) + 1
    for boundary in element_boundaries:
        ax.axhline(boundary, color=DARK, linewidth=1.3)
    all_zero = matrix.fillna(0).abs().sum(axis=1).eq(0).to_numpy()
    for row_index in np.flatnonzero(all_zero):
        for column_index in range(matrix.shape[1]):
            ax.text(
                column_index + 0.5,
                row_index + 0.5,
                "ND",
                ha="center",
                va="center",
                fontsize=7,
                color=GRAY,
            )
    fig.text(
        0.01,
        0.01,
        "Each row is standardized independently; colors must not be compared as absolute abundance across processes.",
        fontsize=8,
        color=GRAY,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    save_pub(fig, figures, "62-process-landscape")


def plot_environment(analysis: Path, figures: Path) -> None:
    rules, order, labels = rule_order(analysis)
    associations = read_tsv(analysis / "environment-associations.tsv")
    display_order = order[::-1]
    y = np.arange(len(display_order))
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 9.0), sharey=True)
    for ax, variable, title in zip(
        axes,
        ("Temperature", "pH"),
        ("A. Temperature", "B. pH"),
    ):
        frame = (
            associations.loc[associations.Variable.eq(variable)]
            .set_index("ProcessID")
            .reindex(display_order)
        )
        for index, row in enumerate(frame.itertuples()):
            if not np.isfinite(row.SpearmanRho):
                continue
            significant = bool(np.isfinite(row.FDR) and row.FDR < 0.05)
            color = ELEMENT_COLORS[row.Element] if significant else "#B0BEC5"
            ax.plot(
                [row.BootstrapLow95, row.BootstrapHigh95],
                [index, index],
                color=color,
                linewidth=1.5,
                alpha=1.0 if significant else 0.8,
            )
            ax.scatter(
                row.SpearmanRho,
                index,
                s=43 if significant else 28,
                color=color,
                edgecolor=WHITE,
                linewidth=0.55,
                zorder=3,
            )
            if significant:
                ax.text(0.95, index, "q<0.05", ha="right", va="center", fontsize=7, color=color)
        ax.axvline(0, color=DARK, linewidth=0.8, linestyle=":")
        ax.set_xlim(-1, 1)
        ax.set_xlabel("Spring-level Spearman correlation (95% bootstrap CI)")
        ax.set_title(title)
    axes[0].set_yticks(y, [wrap(labels[value], 25) for value in display_order])
    axes[0].set_ylabel("Marker-defined process")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color, label=element, markersize=7)
        for element, color in ELEMENT_COLORS.items()
    ]
    axes[1].legend(handles=handles, title="FDR-significant element", loc="lower right")
    sns.despine(fig=fig)
    fig.suptitle(
        "Environmental associations use 56 hot-spring medians as independent units",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    save_pub(fig, figures, "62-environment-associations")


def plot_carrier_map(analysis: Path, figures: Path) -> None:
    _, order, labels = rule_order(analysis)
    frame = read_tsv(analysis / "carrier-phylum-summary.tsv")
    top_phyla = (
        frame.groupby("Phylum")["CarrierMAGs"]
        .sum()
        .sort_values(ascending=False)
        .head(12)
        .index.tolist()
    )
    subset = frame.loc[frame.Phylum.isin(top_phyla)].copy()
    x_lookup = {value: index for index, value in enumerate(top_phyla)}
    display_order = order[::-1]
    y_lookup = {value: index for index, value in enumerate(display_order)}
    subset["x"] = subset.Phylum.map(x_lookup)
    subset["y"] = subset.ProcessID.map(y_lookup)
    abundance = subset.MeanRecoveredMAGAbundance.clip(lower=1e-8)
    color_value = np.log10(abundance)
    norm = Normalize(vmin=float(color_value.min()), vmax=float(color_value.max()))

    fig, ax = plt.subplots(figsize=(13.2, 9.3))
    points = ax.scatter(
        subset.x,
        subset.y,
        s=18 + 25 * np.sqrt(subset.CarrierMAGs),
        c=color_value,
        cmap="viridis",
        norm=norm,
        edgecolor=WHITE,
        linewidth=0.55,
        alpha=0.9,
    )
    ax.set_xticks(np.arange(len(top_phyla)), top_phyla, rotation=38, ha="right")
    ax.set_yticks(
        np.arange(len(display_order)),
        [wrap(labels[value], 25) for value in display_order],
    )
    ax.set_xlim(-0.6, len(top_phyla) - 0.4)
    ax.set_ylim(-0.6, len(display_order) - 0.4)
    ax.set_xlabel("Top phyla among strict carrier MAGs")
    ax.set_ylabel("Marker-defined process")
    ax.set_title("Strict carrier calls are distributed across recovered MAG lineages")
    ax.grid(color="#E3E8EA", linewidth=0.55)
    colorbar = fig.colorbar(points, ax=ax, pad=0.015)
    colorbar.set_label("log10 mean fraction of the recovered-MAG pool")
    size_handles = [
        ax.scatter([], [], s=18 + 25 * np.sqrt(value), color=GRAY, edgecolor=WHITE, label=str(value))
        for value in (1, 10, 50, 100)
    ]
    ax.legend(
        handles=size_handles,
        title="Strict carrier MAGs",
        loc="lower left",
        bbox_to_anchor=(1.15, 0.0),
    )
    fig.tight_layout()
    save_pub(fig, figures, "62-mag-carrier-map")


def plot_concordance(analysis: Path, figures: Path) -> None:
    concordance = read_tsv(analysis / "community-mag-concordance.tsv")
    carriers = read_tsv(analysis / "mag-carrier-summary.tsv")
    eligible = concordance.merge(
        carriers[["ProcessID", "StrictCarrierMAGs"]],
        on="ProcessID",
        validate="one_to_one",
    )
    selected = (
        eligible.loc[eligible.FDR.lt(0.05) & eligible.StrictCarrierMAGs.ge(10)]
        .sort_values("SpearmanRho", ascending=False)
        .head(4)
    )
    if len(selected) != 4:
        raise RuntimeError("Expected four prespecified concordance panels")
    community = read_tsv(analysis / "spring-process-index.tsv")
    mag = read_tsv(analysis / "spring-carrier-fraction.tsv")
    joined = community.merge(
        mag,
        on=["hotspring", "ProcessID", "Element", "Process"],
        validate="one_to_one",
        suffixes=("Community", "MAG"),
    )

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.2))
    points = None
    for ax, row in zip(axes.flat, selected.itertuples(index=False)):
        frame = joined.loc[joined.ProcessID.eq(row.ProcessID)].copy()
        positive = frame.CommunityIndex.loc[frame.CommunityIndex.gt(0)]
        floor = float(positive.min() / 2) if len(positive) else 1e-12
        frame["LogCommunityIndex"] = np.log10(frame.CommunityIndex + floor)
        points = ax.scatter(
            frame.LogCommunityIndex,
            100 * frame.RecoveredMAGCarrierFraction,
            c=100 * frame.RecruitmentRate,
            s=30 + 5 * np.sqrt(frame.SamplesCommunity),
            cmap="viridis",
            vmin=0,
            vmax=80,
            edgecolor=WHITE,
            linewidth=0.55,
            alpha=0.88,
        )
        ax.text(
            0.04,
            0.94,
            f"Spearman ρ = {row.SpearmanRho:.2f}\nFDR = {row.FDR:.3g}\nStrict MAGs = {int(row.StrictCarrierMAGs)}",
            transform=ax.transAxes,
            va="top",
            fontsize=8,
            bbox={"facecolor": WHITE, "edgecolor": "none", "alpha": 0.82},
        )
        ax.set_xlabel("Community KO process index (log10 + process floor)")
        ax.set_ylabel("Carrier fraction within recovered MAGs (%)")
        ax.set_title(wrap(row.Process, 34))
    sns.despine(fig=fig)
    fig.suptitle(
        "Community markers and recovered-MAG carriers agree for selected processes",
        fontsize=12,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.08, right=0.86, bottom=0.08, top=0.91, hspace=0.34, wspace=0.28)
    if points is not None:
        color_axis = fig.add_axes([0.89, 0.18, 0.018, 0.64])
        colorbar = fig.colorbar(points, cax=color_axis)
        colorbar.set_label("Median reads recruited to MAGs (%)")
    save_pub(fig, figures, "62-community-mag-concordance")


def plot_nitrogen_steps(analysis: Path, figures: Path) -> None:
    process_ids = ["nitrate_to_nitrite", "nitrite_to_no", "no_to_n2o", "n2o_to_n2"]
    short_labels = ["NO₃⁻ → NO₂⁻", "NO₂⁻ → NO", "NO → N₂O", "N₂O → N₂"]
    summary = (
        read_tsv(analysis / "nitrogen-step-summary.tsv")
        .set_index("ProcessID")
        .loc[process_ids]
    )
    evidence = read_tsv(analysis / "mag-process-evidence.tsv.gz")
    evidence["StrictCarrier"] = truth(evidence.StrictCarrier)
    matrix = (
        evidence.loc[evidence.ProcessID.isin(process_ids)]
        .pivot(index="MAG", columns="ProcessID", values="StrictCarrier")
        .reindex(columns=process_ids, fill_value=False)
        .fillna(False)
        .astype(bool)
    )
    pattern = matrix.astype(int).astype(str).agg("".join, axis=1)
    patterns = pattern.value_counts().drop(labels="0000", errors="ignore").head(10)

    fig = plt.figure(figsize=(12.8, 9.0))
    grid = fig.add_gridspec(2, 3, height_ratios=[0.85, 1.15], width_ratios=[1, 1, 0.72])
    ax_counts = fig.add_subplot(grid[0, 0])
    ax_fraction = fig.add_subplot(grid[0, 1])
    ax_note = fig.add_subplot(grid[0, 2])
    ax_matrix = fig.add_subplot(grid[1, :2])
    ax_patterns = fig.add_subplot(grid[1, 2], sharey=ax_matrix)

    positions = np.arange(4)
    ax_counts.barh(positions, summary.StrictCarrierMAGs, color=BLUE)
    ax_counts.set_yticks(positions, short_labels)
    ax_counts.invert_yaxis()
    ax_counts.set_xlabel("Strict carrier MAGs")
    ax_counts.set_title("A. Step-specific genomic potential")
    for index, value in enumerate(summary.StrictCarrierMAGs):
        ax_counts.text(value + 2, index, str(int(value)), va="center", fontsize=8)

    ax_fraction.barh(
        positions,
        100 * summary.MedianConditionalFraction,
        color=GREEN,
    )
    ax_fraction.set_yticks(positions, short_labels)
    ax_fraction.invert_yaxis()
    ax_fraction.set_xlabel("Median fraction within recovered MAGs (%)")
    ax_fraction.set_title("B. Conditional abundance")

    complete = int(summary.CompleteChainMAGs.max())
    ax_note.axis("off")
    ax_note.text(
        0.5,
        0.58,
        f"{complete}",
        ha="center",
        va="center",
        fontsize=34,
        fontweight="bold",
        color=ORANGE,
    )
    ax_note.text(
        0.5,
        0.30,
        "MAG carries all\nfour marker steps",
        ha="center",
        va="center",
        fontsize=10,
    )

    pattern_labels = patterns.index.tolist()[::-1]
    y = np.arange(len(pattern_labels))
    for row_index, code in enumerate(pattern_labels):
        present = np.array([char == "1" for char in code])
        ax_matrix.plot(positions, np.repeat(row_index, 4), color="#CFD8DC", linewidth=1.1)
        ax_matrix.scatter(positions[~present], np.repeat(row_index, (~present).sum()), s=32, facecolor=WHITE, edgecolor="#90A4AE")
        ax_matrix.scatter(positions[present], np.repeat(row_index, present.sum()), s=46, color=PURPLE, edgecolor=WHITE, linewidth=0.5)
    ax_matrix.set_xticks(positions, short_labels)
    ax_matrix.set_yticks(y, [f"Pattern {index + 1}" for index in range(len(y))])
    ax_matrix.set_xlim(-0.4, 3.4)
    ax_matrix.set_ylim(-0.6, len(y) - 0.4)
    ax_matrix.set_title("C. Most frequent step co-occurrence patterns across MAGs")
    ax_matrix.set_xlabel("Marker-defined denitrification step")
    ax_matrix.set_ylabel("Presence pattern")
    ax_patterns.barh(y, patterns.loc[pattern_labels], color=PURPLE)
    ax_patterns.set_xlabel("MAGs")
    ax_patterns.tick_params(axis="y", labelleft=False)
    ax_patterns.set_title("Pattern count")
    sns.despine(fig=fig)
    fig.suptitle(
        "Distributed nitrogen-step genes suggest co-occurrence, not demonstrated metabolic handoff",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    save_pub(fig, figures, "62-nitrogen-step-cooccurrence")


def plot_recovery_ceiling(analysis: Path, figures: Path) -> None:
    recovery = read_tsv(analysis / "mag-recovery-ceiling.tsv")
    values = 100 * recovery.MedianRecruitment
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4), gridspec_kw={"width_ratios": [1.25, 0.85]})
    points = axes[0].scatter(
        recovery.MedianTemperature,
        values,
        c=recovery.MedianpH,
        s=28 + 7 * np.sqrt(recovery.Samples),
        cmap="plasma",
        edgecolor=WHITE,
        linewidth=0.6,
        alpha=0.88,
    )
    axes[0].axhline(50, color=ORANGE, linestyle="--", linewidth=1, label="50% recruitment")
    axes[0].set_xlabel("Median spring temperature (°C)")
    axes[0].set_ylabel("Median reads recruited to the 780 MAGs (%)")
    axes[0].set_title("A. Recovery varies among 56 hot springs")
    axes[0].legend(loc="upper left")
    colorbar = fig.colorbar(points, ax=axes[0], pad=0.015)
    colorbar.set_label("Median spring pH")

    ordered = np.sort(values)
    cumulative = np.arange(1, len(ordered) + 1) / len(ordered)
    axes[1].step(ordered, cumulative, where="post", color=BLUE, linewidth=2.2)
    median = float(np.median(values))
    axes[1].axvline(median, color=GREEN, linestyle="--", linewidth=1.2, label=f"Median = {median:.1f}%")
    axes[1].axvline(50, color=ORANGE, linestyle="--", linewidth=1.0)
    axes[1].set_xlabel("Median reads recruited to MAGs (%)")
    axes[1].set_ylabel("Cumulative fraction of hot springs")
    axes[1].set_ylim(0, 1.02)
    axes[1].set_title("B. Empirical recovery ceiling")
    axes[1].legend(loc="lower right")
    sns.despine(fig=fig)
    fig.suptitle(
        "MAG carrier fractions describe the recovered-genome pool, not the whole community",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    save_pub(fig, figures, "62-mag-recovery-ceiling")


def main() -> None:
    args = parse_args()
    setup_style()
    plot_landscape(args.analysis_dir, args.figure_dir)
    plot_environment(args.analysis_dir, args.figure_dir)
    plot_carrier_map(args.analysis_dir, args.figure_dir)
    plot_concordance(args.analysis_dir, args.figure_dir)
    plot_nitrogen_steps(args.analysis_dir, args.figure_dir)
    plot_recovery_ceiling(args.analysis_dir, args.figure_dir)
    print("Article 62 publication figures: PASS")


if __name__ == "__main__":
    main()
