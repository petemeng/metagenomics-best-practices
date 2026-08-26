#!/usr/bin/env python3
"""Create publication-ready, English-only Article 60 figures."""

from __future__ import annotations

import argparse
import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/article60-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


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
TOOL_COLORS = {"gapseq": BLUE, "CarveMe": ORANGE}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if frame.empty:
        raise ValueError(f"Empty plotting table: {path}")
    return frame


def numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    result = frame.copy()
    for column in columns:
        result[column] = pd.to_numeric(result[column], errors="raise")
    return result


def truth(series: pd.Series) -> pd.Series:
    values = series.str.lower()
    if not values.isin(["true", "false"]).all():
        raise ValueError(f"Unexpected logical values: {sorted(values.unique())}")
    return values.eq("true")


def species_short(name: str) -> str:
    words = name.split()
    return f"{words[0][0]}. {' '.join(words[1:])}" if len(words) > 1 else name


def wrap(value: str, width: int = 25) -> str:
    return "\n".join(textwrap.wrap(value, width=width, break_long_words=False))


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
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


def plot_model_size(summary: Path, figures: Path) -> None:
    frame = numeric(
        read_tsv(summary / "model-structure-summary.tsv"),
        ["Reactions", "Metabolites", "Genes"],
    )
    frame = frame.loc[frame["Genome"].str.startswith("SGB_")].copy()
    meta = frame[["Genome", "Species"]].drop_duplicates().set_index("Genome")
    order = (
        frame.loc[frame["Stage"].eq("Draft")]
        .groupby("Genome")["Reactions"]
        .mean()
        .sort_values()
        .index.tolist()
    )
    labels = [f"{genome} · {species_short(meta.loc[genome, 'Species'])}" for genome in order]
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 6.2), sharey=True)
    for ax, tool in zip(axes, ("gapseq", "CarveMe")):
        subset = frame.loc[frame["Tool"].eq(tool)]
        draft = subset.loc[subset["Stage"].eq("Draft")].set_index("Genome").loc[order]
        filled = subset.loc[subset["Stage"].eq("Gap-filled")].set_index("Genome").loc[order]
        y = np.arange(len(order))
        ax.hlines(y, draft["Reactions"], filled["Reactions"], color="#B0BEC5", linewidth=2.0)
        ax.scatter(draft["Reactions"], y, s=46, color=GRAY, label="Draft", zorder=3)
        ax.scatter(filled["Reactions"], y, s=52, color=TOOL_COLORS[tool], label="Gap-filled", zorder=3)
        for index, (left, right) in enumerate(zip(draft["Reactions"], filled["Reactions"])):
            ax.text(max(left, right) + 12, index, f"+{right-left}", va="center", fontsize=7, color=DARK)
        ax.set_title(tool)
        ax.set_xlabel("Model reactions")
        ax.set_yticks(y, labels)
        ax.legend(loc="lower right")
        ax.grid(axis="y", visible=False)
    axes[0].set_ylabel("Quality-audited real MAG")
    fig.suptitle("Gap filling changes model size, but the added reactions are hypotheses", fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_pub(fig, figures, "60-model-size-gapfill")


def plot_gapfill_burden(summary: Path, figures: Path) -> None:
    gaps = numeric(
        read_tsv(summary / "gapfill-burden.tsv"),
        ["AddedFractionOfFilledPct", "AddedReactions", "AddedWithoutGPR"],
    )
    inputs = numeric(
        read_tsv(summary / "input-mag-ledger.tsv"),
        ["RetentionObservedPct"],
    )
    primary = gaps.loc[gaps["Genome"].str.startswith("SGB_")].merge(
        inputs[["Genome", "CompletenessPct", "ContaminationPct"]], on="Genome"
    )
    primary["CompletenessPct"] = pd.to_numeric(primary["CompletenessPct"], errors="raise")
    primary["ContaminationPct"] = pd.to_numeric(primary["ContaminationPct"], errors="raise")
    trunc = gaps.loc[gaps["Genome"].str.startswith("TRUNC_")].merge(
        inputs[["Genome", "RetentionObservedPct"]], on="Genome"
    )
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.3))
    for tool, group in primary.groupby("Tool"):
        axes[0].scatter(
            group["CompletenessPct"], group["AddedFractionOfFilledPct"],
            s=45 + 8 * group["ContaminationPct"], color=TOOL_COLORS[tool],
            edgecolor=WHITE, linewidth=0.7, label=tool,
        )
    audit_labels = primary.loc[
        primary["Genome"].isin(["SGB_008", "SGB_018"])
        | (primary["Genome"].eq("SGB_006") & primary["Tool"].eq("gapseq"))
    ]
    for _, row in audit_labels.iterrows():
        axes[0].annotate(row["Genome"], (row["CompletenessPct"], row["AddedFractionOfFilledPct"]), xytext=(3, 3), textcoords="offset points", fontsize=7)
    axes[0].axvline(90, color=GRAY, linestyle=":", linewidth=1)
    axes[0].axhline(10, color=ORANGE, linestyle="--", linewidth=0.9)
    axes[0].set_xlabel("CheckM2 completeness (%)")
    axes[0].set_ylabel("Added reactions / filled model (%)")
    axes[0].set_title("A. Real MAG quality versus repair burden")
    axes[0].legend(title="Tool")

    for tool, group in trunc.groupby("Tool"):
        group = group.sort_values("RetentionObservedPct")
        axes[1].plot(
            group["RetentionObservedPct"], group["AddedFractionOfFilledPct"],
            marker="o", linewidth=2, color=TOOL_COLORS[tool], label=tool,
        )
        for _, row in group.iterrows():
            axes[1].annotate(str(int(row["AddedReactions"])), (row["RetentionObservedPct"], row["AddedFractionOfFilledPct"]), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=7)
    axes[1].axhline(10, color=ORANGE, linestyle="--", linewidth=0.9, label="10% audit flag")
    axes[1].set_xlabel("Retained parent genome (%)")
    axes[1].set_ylabel("Added reactions / filled model (%)")
    axes[1].set_title("B. Controlled incompleteness sensitivity")
    axes[1].legend(title="Tool / audit")
    sns.despine(fig=fig)
    fig.tight_layout()
    save_pub(fig, figures, "60-gapfill-burden")


def plot_medium_feasibility(summary: Path, figures: Path) -> None:
    frame = read_tsv(summary / "medium-feasibility.tsv")
    frame = frame.loc[
        frame["Stage"].eq("Gap-filled") & ~frame["Medium"].eq("No uptake audit")
    ].copy()
    frame["Growth"] = truth(frame["GrowthAbove1e-6"]).astype(int)
    frame["Column"] = frame["Tool"] + "\n" + frame["Medium"].map(lambda x: wrap(x, 20))
    order = (
        frame[["Genome", "Species"]].drop_duplicates().sort_values("Genome")["Genome"].tolist()
    )
    columns = list(dict.fromkeys(frame["Column"]))
    matrix = frame.pivot(index="Genome", columns="Column", values="Growth").loc[order, columns]
    species = frame[["Genome", "Species"]].drop_duplicates().set_index("Genome")["Species"]
    cmap = LinearSegmentedColormap.from_list("feasibility", ["#ECEFF1", GREEN])
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    sns.heatmap(
        matrix, cmap=cmap, vmin=0, vmax=1, cbar=False, linewidths=0.7,
        linecolor=WHITE, ax=ax, annot=matrix.replace({0: "No", 1: "Yes"}), fmt="",
    )
    ax.set_xticklabels(columns, rotation=35, ha="right")
    ax.set_yticklabels(
        [f"{genome} · {species_short(species[genome])}" for genome in order], rotation=0
    )
    ax.set_xlabel("Constraint set used only for feasibility testing")
    ax.set_ylabel("Real MAG or deterministic truncation")
    ax.set_title("Feasible biomass depends on both gap filling and medium assumptions")
    fig.tight_layout()
    save_pub(fig, figures, "60-medium-feasibility")


def plot_truncation(summary: Path, figures: Path) -> None:
    frame = numeric(
        read_tsv(summary / "truncation-sensitivity.tsv"),
        ["RetentionObservedPct", "Genes", "DraftReactions", "GapfillAddedFractionPct"],
    )
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.7))
    gene = frame.groupby("RetentionObservedPct", as_index=False)["Genes"].first().sort_values("RetentionObservedPct")
    axes[0].plot(gene["RetentionObservedPct"], gene["Genes"], color=PURPLE, marker="o", linewidth=2)
    axes[0].set_xlabel("Retained parent genome (%)")
    axes[0].set_ylabel("Predicted proteins")
    axes[0].set_title("A. Input evidence")
    for tool, group in frame.groupby("Tool"):
        group = group.sort_values("RetentionObservedPct")
        axes[1].plot(group["RetentionObservedPct"], group["DraftReactions"], color=TOOL_COLORS[tool], marker="o", linewidth=2, label=tool)
        axes[2].plot(group["RetentionObservedPct"], group["GapfillAddedFractionPct"], color=TOOL_COLORS[tool], marker="o", linewidth=2, label=tool)
    axes[1].set_xlabel("Retained parent genome (%)")
    axes[1].set_ylabel("Draft reactions")
    axes[1].set_title("B. Evidence-supported draft")
    axes[2].set_xlabel("Retained parent genome (%)")
    axes[2].set_ylabel("Gap-fill burden (%)")
    axes[2].set_title("C. Optimization repair")
    axes[2].axhline(10, color=ORANGE, linestyle="--", linewidth=0.9)
    axes[1].legend(title="Tool")
    sns.despine(fig=fig)
    fig.suptitle("Genome incompleteness propagates into both network evidence and repair", fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_pub(fig, figures, "60-truncation-sensitivity")


def plot_audit(summary: Path, figures: Path) -> None:
    structures = numeric(
        read_tsv(summary / "model-structure-summary.tsv"),
        ["TopologicalDeadEndMetabolites", "Genes"],
    )
    structures = structures.loc[structures["Stage"].eq("Gap-filled")]
    gaps = numeric(read_tsv(summary / "gapfill-burden.tsv"), ["AddedFractionOfFilledPct"])
    frame = structures.merge(gaps[["Genome", "Tool", "AddedFractionOfFilledPct"]], on=["Genome", "Tool"])
    leaks = read_tsv(summary / "medium-feasibility.tsv")
    leaks = leaks.loc[
        leaks["Stage"].eq("Gap-filled") & leaks["Medium"].eq("No uptake audit")
    ][["Genome", "Tool", "GrowthAbove1e-6"]]
    frame = frame.merge(leaks, on=["Genome", "Tool"])
    frame["Leak"] = truth(frame["GrowthAbove1e-6"])
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    for tool, group in frame.groupby("Tool"):
        for leak, subgroup in group.groupby("Leak"):
            ax.scatter(
                subgroup["AddedFractionOfFilledPct"],
                subgroup["TopologicalDeadEndMetabolites"],
                s=np.clip(subgroup["Genes"] / 10, 30, 240),
                color=TOOL_COLORS[tool], marker="X" if leak else "o",
                alpha=0.82, edgecolor=WHITE, linewidth=0.7,
                label=tool if not leak else f"{tool} · no-uptake flag",
            )
    # Label a small, prespecified audit set.  Restricting the crowded low-burden
    # CarveMe cluster to its biological stress case keeps names legible while
    # retaining the points themselves.
    label_keys = {
        ("SGB_006", "gapseq"),
        ("SGB_008", "gapseq"),
        ("SGB_018", "gapseq"),
        ("SGB_021", "gapseq"),
        ("TRUNC_050", "gapseq"),
        ("SGB_018", "CarveMe"),
    }
    label_rows = frame.loc[
        frame.apply(lambda row: (row["Genome"], row["Tool"]) in label_keys, axis=1)
    ]
    for _, row in label_rows.iterrows():
        ax.annotate(
            row["Genome"],
            (row["AddedFractionOfFilledPct"], row["TopologicalDeadEndMetabolites"]),
            xytext=(4, 4), textcoords="offset points", fontsize=7,
        )
    ax.axvline(10, color=ORANGE, linestyle="--", linewidth=1, label="10% burden audit flag")
    ax.set_xlabel("Added reactions / filled model (%)")
    ax.set_ylabel("Topological dead-end metabolites")
    ax.set_title("Model quality is multidimensional after a solver returns a network")
    ax.legend(title="Tool / audit")
    sns.despine(fig=fig)
    fig.tight_layout()
    save_pub(fig, figures, "60-model-audit")


def plot_evidence_ladder(summary: Path, figures: Path) -> None:
    frame = numeric(read_tsv(summary / "evidence-ladder.tsv"), ["Level"])
    fig, ax = plt.subplots(figsize=(12.0, 7.0))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, len(frame) + 1)
    ax.axis("off")
    colors = ["#E3F2FD", "#E8F5E9", "#FFF8E1", "#FCE4EC", "#EDE7F6", "#ECEFF1"]
    for index, row in frame.sort_values("Level").iterrows():
        level = int(row["Level"])
        y = len(frame) - level + 0.7
        box = FancyBboxPatch(
            (0.3, y), 11.4, 0.78, boxstyle="round,pad=0.03,rounding_size=0.06",
            facecolor=colors[(level - 1) % len(colors)], edgecolor="#90A4AE", linewidth=0.8,
        )
        ax.add_patch(box)
        ax.text(0.55, y + 0.39, str(level), ha="center", va="center", fontsize=11, fontweight="bold", color=DARK)
        ax.text(1.0, y + 0.50, row["Evidence"], va="center", fontsize=9.4, fontweight="bold", color=DARK)
        ax.text(1.0, y + 0.22, f"Supports: {row['SupportedClaim']}", va="center", fontsize=8.2, color=GREEN)
        ax.text(6.4, y + 0.22, f"Does not prove: {row['UnsupportedClaim']}", va="center", fontsize=8.2, color=ORANGE)
    ax.text(6, len(frame) + 0.55, "Evidence ladder for MAG-derived genome-scale models", ha="center", va="center", fontsize=13, fontweight="bold", color=DARK)
    ax.annotate(
        "More external validation",
        xy=(11.55, 1.2), xytext=(11.55, 6.0),
        arrowprops={"arrowstyle": "->", "color": GRAY, "lw": 1.5},
        rotation=90, ha="center", va="center", color=GRAY, fontsize=8,
    )
    save_pub(fig, figures, "60-evidence-ladder")


def main() -> None:
    args = parse_args()
    summary = args.summary_dir.resolve()
    figures = args.figure_dir.resolve()
    setup_style()
    plot_model_size(summary, figures)
    plot_gapfill_burden(summary, figures)
    plot_medium_feasibility(summary, figures)
    plot_truncation(summary, figures)
    plot_audit(summary, figures)
    plot_evidence_ladder(summary, figures)
    print("PASS\tArticle 60 figures\t6 figure families")


if __name__ == "__main__":
    main()
