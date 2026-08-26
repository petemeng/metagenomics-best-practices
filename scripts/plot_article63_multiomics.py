#!/usr/bin/env python3
"""Create publication-ready, English-only figures for Article 63."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import textwrap

os.environ.setdefault("MPLCONFIGDIR", "/tmp/article63-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
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
DIAGNOSIS_COLORS = {"Control": GREEN, "CD": BLUE, "UC": ORANGE}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--summary-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


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
        directory / f"{stem}.tiff", dpi=350, bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def wrap(value: object, width: int = 25) -> str:
    return "\n".join(
        textwrap.wrap(str(value), width=width, break_long_words=False)
    )


def plot_audit(source: Path, figures: Path) -> None:
    metadata = read_tsv(source / "sample-metadata.tsv")
    attrition = read_tsv(source / "feature-attrition.tsv")
    counts = (
        metadata.groupby(["Cohort", "Study.Group"], observed=True)
        .size().unstack(fill_value=0)
        .reindex(index=["PRISM", "Validation"], columns=["Control", "CD", "UC"])
    )
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.8),
                             gridspec_kw={"width_ratios": [1.05, 1.45, 0.9]})
    bottom = np.zeros(len(counts))
    x = np.arange(len(counts))
    for diagnosis in counts.columns:
        values = counts[diagnosis].to_numpy()
        axes[0].bar(x, values, bottom=bottom, color=DIAGNOSIS_COLORS[diagnosis],
                    width=0.62, label=diagnosis)
        for position, (base, value) in enumerate(zip(bottom, values)):
            axes[0].text(position, base + value / 2, str(int(value)), ha="center",
                         va="center", color=WHITE, fontweight="bold", fontsize=8)
        bottom += values
    axes[0].set_xticks(x, ["PRISM\ndiscovery", "Independent\nValidation"])
    axes[0].set_ylabel("Independent subjects")
    axes[0].set_title("A. Locked cohort split")
    axes[0].legend(title="Diagnosis", loc="upper right")

    modalities = ["Microbiome", "Metabolome"]
    colors = {"Microbiome": BLUE, "Metabolome": PURPLE}
    offsets = {"Microbiome": -0.16, "Metabolome": 0.16}
    stages = []
    for modality in modalities:
        current = attrition.loc[attrition.Modality.eq(modality)]
        for order, row in enumerate(current.itertuples(index=False)):
            label = wrap(row.Stage, 22)
            stages.append((modality, order, label, row.Features))
    max_stages = max(item[1] for item in stages) + 1
    for modality in modalities:
        current = [item for item in stages if item[0] == modality]
        positions = np.arange(max_stages) + offsets[modality]
        values = np.array([item[3] for item in current])
        axes[1].plot(positions, values, marker="o", linewidth=2.2,
                     color=colors[modality], label=modality)
        for position, value in zip(positions, values):
            axes[1].text(position, value * 1.13, f"{int(value):,}", ha="center",
                         color=colors[modality], fontsize=7.5)
    stage_labels = [wrap(value, 17) for value in (
        "Source table", "Annotation or prevalence gate", "Primary analysis set"
    )]
    axes[1].set_xticks(np.arange(max_stages), stage_labels)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Features (log scale)")
    axes[1].set_title("B. Prespecified feature attrition")
    axes[1].legend()

    audit_lines = [
        ("220", "paired samples"),
        ("220", "unique subjects"),
        ("166", "CLR genera"),
        ("153", "high-confidence HMDB metabolites"),
    ]
    axes[2].set_axis_off()
    axes[2].add_patch(plt.Rectangle((0.04, 0.04), 0.92, 0.92,
                                    transform=axes[2].transAxes,
                                    facecolor=LIGHT, edgecolor="#CFD8DC"))
    for index, (number, label) in enumerate(audit_lines):
        y = 0.82 - index * 0.21
        axes[2].text(0.13, y, number, transform=axes[2].transAxes,
                     fontsize=18, fontweight="bold", color=DARK, va="center")
        axes[2].text(0.13, y - 0.08, wrap(label, 27), transform=axes[2].transAxes,
                     fontsize=8.5, color=GRAY, va="center")
    axes[2].set_title("C. Pairing and annotation gates")
    sns.despine(fig=fig)
    fig.suptitle("Paired Franzosa microbiome–metabolome data with an untouched external cohort",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_pub(fig, figures, "63-paired-data-audit")


def plot_procrustes(source: Path, figures: Path) -> None:
    scores = read_tsv(source / "global/procrustes-scores.tsv")
    tests = read_tsv(source / "global/global-concordance.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.8))
    for ax, cohort in zip(axes, ["PRISM", "Validation"]):
        current = scores.loc[scores.Cohort.eq(cohort)]
        for row in current.itertuples(index=False):
            color = DIAGNOSIS_COLORS[row.Diagnosis]
            ax.plot([row.MicrobiomeAxis1, row.MetabolomeAxis1],
                    [row.MicrobiomeAxis2, row.MetabolomeAxis2],
                    color=color, alpha=0.17, linewidth=0.55, zorder=1)
        for diagnosis in ["Control", "CD", "UC"]:
            subset = current.loc[current.Diagnosis.eq(diagnosis)]
            ax.scatter(subset.MicrobiomeAxis1, subset.MicrobiomeAxis2,
                       s=14, color=DIAGNOSIS_COLORS[diagnosis], alpha=0.66,
                       marker="o", linewidth=0, zorder=2)
            ax.scatter(subset.MetabolomeAxis1, subset.MetabolomeAxis2,
                       s=18, facecolor=WHITE, edgecolor=DIAGNOSIS_COLORS[diagnosis],
                       alpha=0.72, marker="^", linewidth=0.7, zorder=3)
        result = tests.loc[
            tests.Cohort.eq(cohort) & tests.Restriction.eq("Within diagnosis")
        ].iloc[0]
        title = "A. PRISM discovery" if cohort == "PRISM" else "B. Independent Validation"
        ax.set_title(
            f"{title}\nProcrustes r = {result.ProcrustesR:.3f}, "
            f"P = {result.ProcrustesP:.4f}"
        )
        ax.set_xlabel("Symmetric Procrustes axis 1")
        ax.set_ylabel("Symmetric Procrustes axis 2")
        ax.set_aspect("equal", adjustable="datalim")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DARK,
               label="Microbiome", markersize=5),
        Line2D([0], [0], marker="^", color=DARK, markerfacecolor=WHITE,
               label="Metabolome", markersize=5),
    ] + [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=DIAGNOSIS_COLORS[group], label=group, markersize=5)
        for group in ["Control", "CD", "UC"]
    ]
    axes[1].legend(handles=handles, loc="best", ncol=2)
    sns.despine(fig=fig)
    fig.suptitle("Global sample geometry replicates after projection through PRISM-fitted axes",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_pub(fig, figures, "63-procrustes-concordance")


def plot_global_tests(source: Path, figures: Path) -> None:
    tests = read_tsv(source / "global/global-concordance.tsv")
    restricted = tests.loc[tests.Restriction.eq("Within diagnosis")].copy()
    colors = {"PRISM": BLUE, "Validation": ORANGE}
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.5))
    y = np.arange(2)
    for index, row in restricted.reset_index(drop=True).iterrows():
        axes[0].errorbar(row.ProcrustesR, index,
                         xerr=[[row.ProcrustesR - row.ProcrustesBootLow],
                               [row.ProcrustesBootHigh - row.ProcrustesR]],
                         fmt="o", color=colors[row.Cohort], capsize=4,
                         markersize=7, linewidth=1.8)
        axes[0].text(row.ProcrustesBootHigh + 0.015, index,
                     f"P = {row.ProcrustesP:.4f}", va="center", fontsize=8)
        axes[1].scatter(row.MantelRho, index, s=58, color=colors[row.Cohort])
        axes[1].text(row.MantelRho + 0.015, index,
                     f"P = {row.MantelP:.4f}", va="center", fontsize=8)
    labels = ["PRISM discovery", "Independent Validation"]
    for ax in axes:
        ax.set_yticks(y, labels)
        ax.set_xlim(0, 0.86)
        ax.axvline(0, color=GRAY, linewidth=0.8)
    axes[0].set_xlabel("Procrustes correlation (95% subject bootstrap CI)")
    axes[0].set_title("A. Concordance of ten-dimensional ordinations")
    axes[1].set_xlabel("Spearman Mantel correlation")
    axes[1].set_title("B. Concordance of pairwise distances")
    fig.text(0.5, 0.01,
             "P values use 9,999 permutations restricted within diagnosis; statistics are global, not feature-level mechanisms.",
             ha="center", color=GRAY, fontsize=8.5)
    sns.despine(fig=fig)
    fig.suptitle("Two global tests agree in discovery and independent validation",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    save_pub(fig, figures, "63-global-concordance")


def plot_halla(source: Path, summary: Path, figures: Path) -> None:
    pairs = read_tsv(summary / "halla-pair-validation.tsv.gz")
    stages = read_tsv(summary / "halla-replication-summary.tsv")
    overlap = read_tsv(summary / "halla-branch-overlap.tsv")
    category = np.where(
        pairs.Replicated,
        "Replicated after BH",
        np.where(pairs.SameDirection, "Same direction, not BH", "Direction discordant"),
    )
    pairs = pairs.assign(Category=category)
    palette = {
        "Replicated after BH": BLUE,
        "Same direction, not BH": "#B0BEC5",
        "Direction discordant": ORANGE,
    }
    order = ["Same direction, not BH", "Direction discordant", "Replicated after BH"]
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.2),
                             gridspec_kw={"width_ratios": [1.45, 1.0, 0.9]})
    for name in order:
        current = pairs.loc[pairs.Category.eq(name)]
        axes[0].scatter(current.DiscoveryRho, current.ValidationRho,
                        s=10 if name != "Replicated after BH" else 13,
                        color=palette[name], alpha=0.45 if name != "Replicated after BH" else 0.62,
                        linewidth=0, label=f"{name} ({len(current):,})")
    axes[0].axhline(0, color=GRAY, linewidth=0.7)
    axes[0].axvline(0, color=GRAY, linewidth=0.7)
    axes[0].plot([-0.8, 0.8], [-0.8, 0.8], color=DARK, linestyle=":", linewidth=0.9)
    axes[0].set_xlim(-0.78, 0.78)
    axes[0].set_ylim(-0.78, 0.78)
    axes[0].set_xlabel("PRISM covariate-adjusted Spearman rho")
    axes[0].set_ylabel("Validation covariate-adjusted Spearman rho")
    axes[0].set_title("A. Prespecified pairwise replication")
    axes[0].legend(loc="lower right", fontsize=7.2)

    stage_labels = [
        "Discovery BH",
        "Same direction",
        "Validation nominal",
        "Validation BH",
        "Direction + BH",
    ]
    values = stages.Pairs.to_numpy()
    positions = np.arange(len(stages))[::-1]
    axes[1].barh(positions, values, color=[PURPLE, SKY, YELLOW, GREEN, BLUE])
    axes[1].set_yticks(positions, stage_labels)
    axes[1].set_xlabel("Microbe–metabolite pairs")
    axes[1].set_title("B. Replication gates")
    for position, value in zip(positions, values):
        axes[1].text(value + 65, position, f"{int(value):,}", va="center", fontsize=8)
    axes[1].set_xlim(0, values.max() * 1.22)

    overlap_order = ["Both", "Adjusted only", "Raw only", "Neither"]
    overlap = overlap.set_index("EvidenceClass").loc[overlap_order]
    axes[2].bar(np.arange(4), overlap.Pairs,
                color=[PURPLE, BLUE, ORANGE, "#CFD8DC"])
    axes[2].set_xticks(np.arange(4), ["Both", "Adjusted\nonly", "Raw\nonly", "Neither"])
    axes[2].set_ylabel("All 25,398 tested pairs")
    axes[2].set_title("C. Covariate sensitivity")
    for position, value in enumerate(overlap.Pairs):
        axes[2].text(position, value + 280, f"{int(value):,}", ha="center", fontsize=7.5)
    sns.despine(fig=fig)
    fig.suptitle("HAllA discovery contracts sharply after covariate adjustment and external BH control",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_pub(fig, figures, "63-halla-discovery-replication")


def plot_diablo_validation(source: Path, figures: Path) -> None:
    diablo = source / "diablo"
    confusion = read_tsv(diablo / "external-confusion.tsv")
    metrics = read_tsv(diablo / "external-metrics.tsv")
    null = read_tsv(diablo / "label-permutation-null.tsv")
    null_summary = read_tsv(diablo / "label-permutation-summary.tsv").iloc[0]
    classes = ["Control", "CD", "UC"]
    matrix = confusion.pivot(index="Truth", columns="Predicted", values="Samples").reindex(
        index=classes, columns=classes, fill_value=0
    )
    row_percent = matrix.div(matrix.sum(axis=1), axis=0) * 100
    annotations = matrix.astype(str) + "\n(" + row_percent.round(1).astype(str) + "%)"
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.8),
                             gridspec_kw={"width_ratios": [1.1, 1.05, 1.15]})
    sns.heatmap(row_percent, annot=annotations, fmt="", cmap="Blues", vmin=0, vmax=100,
                linewidths=0.6, linecolor=WHITE, cbar_kws={"label": "Row percentage"},
                ax=axes[0])
    axes[0].set_xlabel("Predicted diagnosis")
    axes[0].set_ylabel("Observed diagnosis")
    axes[0].set_title("A. Untouched 65-subject cohort")

    metric_order = ["Accuracy", "BalancedAccuracy", "MacroF1"]
    display = {"Accuracy": "Accuracy", "BalancedAccuracy": "Balanced accuracy",
               "MacroF1": "Macro F1"}
    current = metrics.set_index("Metric").loc[metric_order].reset_index()
    y = np.arange(len(current))[::-1]
    axes[1].errorbar(current.Estimate, y,
                     xerr=np.vstack([current.Estimate - current.Low,
                                     current.High - current.Estimate]),
                     fmt="o", color=BLUE, capsize=4, linewidth=1.8, markersize=7)
    axes[1].axvline(1 / 3, color=ORANGE, linestyle="--", linewidth=1,
                    label="Three-class chance")
    axes[1].set_yticks(y, [display[value] for value in current.Metric])
    axes[1].set_xlim(0.2, 0.88)
    axes[1].set_xlabel("Estimate (95% subject bootstrap CI)")
    axes[1].set_title("B. External performance")
    axes[1].legend(loc="lower right")

    axes[2].hist(null.BalancedAccuracy, bins=20, color="#B0BEC5", edgecolor=WHITE)
    axes[2].axvline(null_summary.ObservedBalancedAccuracy, color=RED,
                    linewidth=2.2, label="Observed external model")
    axes[2].axvline(null_summary.NullQ975, color=GRAY, linestyle=":",
                    linewidth=1.2, label="Null 97.5th percentile")
    axes[2].set_xlabel("External balanced accuracy")
    axes[2].set_ylabel("Label-permuted training fits")
    axes[2].set_title(f"C. Fixed-complexity null\nEmpirical P = {null_summary.EmpiricalP:.4f}")
    axes[2].legend(loc="upper left")
    sns.despine(fig=fig)
    fig.suptitle("DIABLO tuning stays inside PRISM before one independent external evaluation",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_pub(fig, figures, "63-diablo-external-validation")


def plot_diablo_stability(source: Path, summary: Path, figures: Path) -> None:
    selected = read_tsv(summary / "diablo-selected-stability.tsv")
    latent = read_tsv(source / "diablo/latent-scores.tsv")
    top_rows = []
    for block in ["microbiome", "metabolome"]:
        current = selected.loc[selected.Block.eq(block)].sort_values(
            ["SelectionFrequency", "Loading"], ascending=[False, False]
        ).head(10)
        top_rows.append(current)
    top = pd.concat(top_rows, ignore_index=True)
    top["Label"] = top.apply(
        lambda row: f"C{int(row.Component)} · {wrap(row.DisplayName, 18)}", axis=1
    )
    top = top.sort_values(["Block", "SelectionFrequency"], ascending=[True, True])
    fig, axes = plt.subplots(1, 3, figsize=(16.6, 6.0),
                             gridspec_kw={"width_ratios": [1.3, 1.0, 1.0]})
    y = np.arange(len(top))
    bar_colors = [BLUE if value == "microbiome" else PURPLE for value in top.Block]
    axes[0].barh(y, top.SelectionFrequency, color=bar_colors)
    axes[0].axvline(0.70, color=ORANGE, linestyle="--", linewidth=1,
                    label="70% stability gate")
    axes[0].set_yticks(y, top.Label)
    axes[0].set_xlim(0, 1.05)
    axes[0].set_xlabel("Selection frequency across 100 PRISM bootstraps")
    axes[0].set_title("A. Most stable final-model features")
    axes[0].legend(loc="lower right")

    validation = latent.loc[latent.Cohort.eq("Validation")]
    for ax, component in zip(axes[1:], [1, 2]):
        current = validation.loc[validation.Component.eq(component)]
        for diagnosis in ["Control", "CD", "UC"]:
            subset = current.loc[current.Diagnosis.eq(diagnosis)]
            ax.scatter(subset.MicrobiomeScore, subset.MetabolomeScore,
                       color=DIAGNOSIS_COLORS[diagnosis], s=30, alpha=0.75,
                       edgecolor=WHITE, linewidth=0.5, label=diagnosis)
        rho = current[["MicrobiomeScore", "MetabolomeScore"]].corr(
            method="spearman"
        ).iloc[0, 1]
        ax.set_xlabel(f"Microbiome latent component {component}")
        ax.set_ylabel(f"Metabolome latent component {component}")
        ax.set_title(f"{'B' if component == 1 else 'C'}. Validation component {component}\nSpearman rho = {rho:.3f}")
    axes[2].legend(title="Diagnosis", loc="best")
    sns.despine(fig=fig)
    fig.suptitle("Sparse features and cross-block latent structure are audited separately",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_pub(fig, figures, "63-diablo-feature-stability")


def main() -> None:
    args = parse_args()
    setup_style()
    plot_audit(args.input_dir, args.figure_dir)
    plot_procrustes(args.input_dir, args.figure_dir)
    plot_global_tests(args.input_dir, args.figure_dir)
    plot_halla(args.input_dir, args.summary_dir, args.figure_dir)
    plot_diablo_validation(args.input_dir, args.figure_dir)
    plot_diablo_stability(args.input_dir, args.summary_dir, args.figure_dir)
    print("Article 63 publication figures: PASS")


if __name__ == "__main__":
    main()
