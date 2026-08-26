#!/usr/bin/env python3
"""Create deterministic publication figures for Article 71."""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image


PLOT_SEED = 20_260_771
COLORS = {
    "Exposure": "#D55E00",
    "Microbiome": "#2A9D8F",
    "Phenotype": "#4C78A8",
    "CD": "#E45756",
    "UC": "#59A14F",
    "Control": "#4C78A8",
    "Neutral": "#9AA5B1",
    "Warning": "#F2C14E",
    "Purple": "#B279A2",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    return parser.parse_args()


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.6,
            "axes.titlesize": 10,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.7,
            "ytick.labelsize": 7.7,
            "legend.fontsize": 7.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.10,
        1.06,
        label,
        transform=axis.transAxes,
        fontweight="bold",
        fontsize=10,
        va="top",
    )


def save(fig: plt.Figure, output: Path, stem: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output / f"{stem}.png", dpi=360, bbox_inches="tight")
    fig.savefig(
        output / f"{stem}.tiff",
        dpi=360,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)
    print(stem)


def data_positivity(input_dir: Path, output: Path) -> None:
    attrition = pd.read_csv(input_dir / "sample-attrition.tsv", sep="\t")
    overlap = pd.read_csv(
        input_dir / "antibiotic-overlap-by-diagnosis.tsv", sep="\t"
    )
    cohort = pd.read_csv(input_dir / "sem-primary-cohort.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(11.3, 4.15))

    axis = axes[0]
    display = attrition.iloc[::-1]
    y = np.arange(len(display))
    bars = axis.barh(y, display["Subjects"], color="#8FB3C9", height=0.62)
    axis.set_yticks(y, display["Stage"].str.replace(" available", "", regex=False))
    axis.set_xlabel("Subjects")
    axis.set_title("Complete-case losses are visible", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    for bar, value in zip(bars, display["Subjects"], strict=True):
        axis.text(value + 3, bar.get_y() + bar.get_height() / 2, str(value), va="center")
    panel_label(axis, "A")

    axis = axes[1]
    x = np.arange(len(overlap))
    axis.bar(x, overlap["Unexposed"], color="#B9C2C9", label="Unexposed")
    axis.bar(
        x,
        overlap["Exposed"],
        bottom=overlap["Unexposed"],
        color=COLORS["Exposure"],
        label="Antibiotic exposed",
    )
    axis.set_xticks(x, overlap["Diagnosis"])
    axis.set_ylabel("PRISM complete-case subjects")
    axis.set_title("No exposed control subject", loc="left", fontweight="bold")
    axis.legend(frameon=False)
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    for xx, row in overlap.iterrows():
        axis.text(
            xx,
            row["Unexposed"] + row["Exposed"] + 1,
            f"{int(row['Exposed'])}/{int(row['Unexposed'] + row['Exposed'])}",
            ha="center",
            fontsize=7.5,
        )
    panel_label(axis, "B")

    axis = axes[2]
    for diagnosis in ("Control", "CD", "UC"):
        subset = cohort.loc[cohort["Diagnosis"].eq(diagnosis)]
        axis.scatter(
            subset["Shannon"],
            subset["LogCalprotectin"],
            s=np.where(subset["Antibiotic"].eq(1), 55, 28),
            marker="o",
            color=COLORS[diagnosis],
            edgecolor=np.where(subset["Antibiotic"].eq(1), "#D55E00", "white"),
            linewidth=np.where(subset["Antibiotic"].eq(1), 1.5, 0.5),
            alpha=0.78,
            label=diagnosis,
        )
    axis.set_xlabel("Genus-profile Shannon entropy")
    axis.set_ylabel("log(1 + fecal calprotectin)")
    axis.set_title("All three nodes share one visit", loc="left", fontweight="bold")
    axis.legend(frameon=False, title="Diagnosis")
    axis.grid(color="#E8E8E8", linewidth=0.55)
    axis.text(
        0.02,
        0.02,
        "Orange outline: antibiotic exposed",
        transform=axis.transAxes,
        fontsize=7.3,
        color="#6B4B00",
    )
    panel_label(axis, "C")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "71-data-positivity")


def add_box(
    axis: plt.Axes,
    xy: tuple[float, float],
    text: str,
    color: str,
    width: float = 0.23,
    height: float = 0.16,
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.018",
        facecolor=color,
        edgecolor="#59636A",
        linewidth=0.9,
    )
    axis.add_patch(patch)
    axis.text(x, y, text, ha="center", va="center", fontsize=8, fontweight="bold")


def add_arrow(
    axis: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = "#56636B",
    style: str = "-",
    rad: float = 0.0,
) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=11,
        linewidth=1.35,
        color=color,
        linestyle=style,
        connectionstyle=f"arc3,rad={rad}",
    )
    axis.add_patch(arrow)


def prespecified_dag(input_dir: Path, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.15))
    axis = axes[0]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    add_box(axis, (0.17, 0.62), "Antibiotic\nmetadata", "#F7D6CE")
    add_box(axis, (0.50, 0.62), "Shannon\ndiversity", "#D4ECE7")
    add_box(axis, (0.83, 0.62), "log(1 + fecal\ncalprotectin)", "#D7E4F2", 0.25)
    add_box(
        axis,
        (0.50, 0.20),
        "Diagnosis + age +\nthree medications",
        "#ECEFF1",
        0.34,
        0.15,
    )
    add_arrow(axis, (0.29, 0.62), (0.38, 0.62), COLORS["Exposure"])
    add_arrow(axis, (0.62, 0.62), (0.70, 0.62), COLORS["Microbiome"])
    add_arrow(axis, (0.26, 0.68), (0.71, 0.68), COLORS["Exposure"], rad=-0.24)
    add_arrow(axis, (0.43, 0.27), (0.25, 0.52), COLORS["Neutral"])
    add_arrow(axis, (0.50, 0.28), (0.50, 0.51), COLORS["Neutral"])
    add_arrow(axis, (0.57, 0.27), (0.76, 0.52), COLORS["Neutral"])
    axis.text(0.335, 0.66, "a", color=COLORS["Exposure"], fontweight="bold")
    axis.text(0.665, 0.66, "b", color=COLORS["Microbiome"], fontweight="bold")
    axis.text(0.49, 0.84, "direct path", color=COLORS["Exposure"], ha="center")
    axis.set_title("Prespecified conditional-association graph", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    axis.set_xlim(0, 10)
    axis.set_ylim(0, 4)
    axis.axis("off")
    axis.plot([0.8, 9.2], [2.0, 2.0], color="#707B82", linewidth=1.4)
    events = [
        (1.4, "Treatment\ninitiation", "unknown"),
        (4.6, "Antibiotic\nmetadata", "recorded"),
        (6.0, "Stool genus\nprofile", "recorded"),
        (7.5, "Fecal\ncalprotectin", "recorded"),
    ]
    for x, label, status in events:
        color = COLORS["Warning"] if status == "unknown" else "#D7E4F2"
        axis.scatter(x, 2.0, s=100, color=color, edgecolor="#59636A", zorder=3)
        axis.text(x, 2.42, label, ha="center", va="bottom", fontweight="bold")
        axis.text(x, 1.62, status, ha="center", va="top", color="#59636A")
    axis.annotate(
        "",
        xy=(7.9, 0.72),
        xytext=(4.1, 0.72),
        arrowprops={"arrowstyle": "<->", "color": "#D55E00", "linestyle": "--"},
    )
    axis.text(
        6.0,
        0.48,
        "Ordering inside the visit is not identified",
        ha="center",
        color="#8B3D00",
        fontweight="bold",
    )
    axis.set_title("A path arrow is not a time stamp", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "71-prespecified-dag")


def local_paths(input_dir: Path, output: Path) -> None:
    coefficients = pd.read_csv(
        input_dir / "local-path-coefficients-hc3.tsv", sep="\t"
    )
    effects = pd.read_csv(input_dir / "path-effect-summary.tsv", sep="\t")
    primary = effects.loc[effects["Model"].eq("Primary Shannon path")].set_index("Effect")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), gridspec_kw={"width_ratios": [0.95, 1.25]})

    axis = axes[0]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    add_box(axis, (0.16, 0.62), "Antibiotic", "#F7D6CE", 0.22, 0.14)
    add_box(axis, (0.50, 0.62), "Shannon", "#D4ECE7", 0.22, 0.14)
    add_box(axis, (0.84, 0.62), "Calprotectin", "#D7E4F2", 0.24, 0.14)
    add_arrow(axis, (0.27, 0.62), (0.39, 0.62), COLORS["Exposure"])
    add_arrow(axis, (0.61, 0.62), (0.71, 0.62), COLORS["Microbiome"])
    add_arrow(axis, (0.24, 0.68), (0.74, 0.68), COLORS["Exposure"], rad=-0.28)
    axis.text(
        0.30,
        0.47,
        f"a = {primary.loc['A', 'Estimate']:.3f}\n95% boot CI {primary.loc['A', 'CILower']:.2f}, {primary.loc['A', 'CIUpper']:.2f}",
        ha="center",
        va="top",
        color="#8B3D00",
        fontsize=7.2,
    )
    axis.text(
        0.70,
        0.35,
        f"b = {primary.loc['B', 'Estimate']:.3f}\n95% boot CI {primary.loc['B', 'CILower']:.2f}, {primary.loc['B', 'CIUpper']:.2f}",
        ha="center",
        va="top",
        color="#176B5B",
        fontsize=7.2,
    )
    axis.text(
        0.50,
        0.85,
        f"direct = {primary.loc['Direct', 'Estimate']:.3f}\n95% boot CI {primary.loc['Direct', 'CILower']:.2f}, {primary.loc['Direct', 'CIUpper']:.2f}",
        ha="center",
        color="#8B3D00",
    )
    axis.text(
        0.50,
        0.10,
        "All continuous nodes use PRISM SD units.\nBinary exposure remains a 0-to-1 contrast.",
        ha="center",
        color="#59636A",
    )
    axis.set_title("Only the exposure-to-diversity path is precise", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    selected = coefficients.loc[
        (
            coefficients["Model"].eq("Microbiome node")
            & coefficients["Term"].isin(["Antibiotic", "CD", "UC"])
        )
        | (
            coefficients["Model"].eq("Phenotype node")
            & coefficients["Term"].isin(["Antibiotic", "ShannonZ", "CD", "UC"])
        )
    ].copy()
    selected["Label"] = (
        selected["Model"].str.replace(" node", "", regex=False)
        + ": "
        + selected["Term"].replace({"ShannonZ": "Shannon"})
    )
    selected = selected.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(selected))
    colors = [
        COLORS["Exposure"] if term == "Antibiotic" else COLORS["Microbiome"]
        if term == "ShannonZ" else COLORS["Neutral"]
        for term in selected["Term"]
    ]
    for yy, row, color in zip(y, selected.itertuples(), colors, strict=True):
        axis.errorbar(
            row.Estimate,
            yy,
            xerr=[[row.Estimate - row.CILower], [row.CIUpper - row.Estimate]],
            fmt="o",
            color=color,
            ecolor=color,
            markersize=5,
            elinewidth=1.5,
            capsize=2.5,
            zorder=3,
        )
    axis.axvline(0, color="#555555", linestyle="--", linewidth=0.9)
    axis.set_yticks(y, selected["Label"])
    axis.set_xlabel("HC3 coefficient (standardized continuous outcome)")
    axis.set_title("Local equations, not one global causal coefficient", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "71-local-paths")


def path_decomposition(input_dir: Path, output: Path) -> None:
    summary = pd.read_csv(input_dir / "path-effect-summary.tsv", sep="\t")
    summary = summary.loc[
        summary["Model"].eq("Primary Shannon path")
        & summary["Effect"].isin(["Direct", "Indirect", "Total"])
    ].set_index("Effect").loc[["Direct", "Indirect", "Total"]].reset_index()
    draws = pd.read_csv(input_dir / "sem-path-bootstrap.tsv.gz", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.35))

    axis = axes[0]
    y = np.arange(len(summary))[::-1]
    palette = [COLORS["Exposure"], COLORS["Microbiome"], COLORS["Purple"]]
    for yy, row, color in zip(y, summary.itertuples(), palette, strict=True):
        axis.errorbar(
            row.Estimate,
            yy,
            xerr=[[row.Estimate - row.CILower], [row.CIUpper - row.Estimate]],
            fmt="o",
            color=color,
            ecolor=color,
            markersize=6,
            elinewidth=2,
            capsize=3,
            zorder=3,
        )
    axis.axvline(0, color="#555555", linestyle="--", linewidth=0.9)
    axis.set_yticks(y, summary["Effect"])
    axis.set_xlabel("Path effect (PRISM SD units)")
    axis.set_title("Indirect effect is uncertain", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    for yy, row in zip(y, summary.itertuples(), strict=True):
        axis.text(
            max(summary["CIUpper"]) + 0.04,
            yy,
            f"{row.Estimate:.2f} [{row.CILower:.2f}, {row.CIUpper:.2f}]",
            va="center",
            fontsize=7.5,
        )
    panel_label(axis, "A")

    axis = axes[1]
    bins = np.linspace(
        draws[["Direct", "Indirect", "Total"]].min().min(),
        draws[["Direct", "Indirect", "Total"]].max().max(),
        45,
    )
    for effect, color in zip(("Direct", "Indirect", "Total"), palette, strict=True):
        axis.hist(
            draws[effect],
            bins=bins,
            density=True,
            histtype="step",
            linewidth=1.7,
            color=color,
            label=effect,
        )
    axis.axvline(0, color="#555555", linestyle="--", linewidth=0.9)
    axis.set_xlabel("Bootstrap path effect")
    axis.set_ylabel("Density")
    axis.set_title("5,000 full path refits", loc="left", fontweight="bold")
    axis.legend(frameon=False)
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    axis.text(
        0.02,
        0.03,
        "Opposite direct and indirect signs\nmake a mediated proportion unstable.",
        transform=axis.transAxes,
        fontsize=7.5,
        color="#59636A",
    )
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.1)
    save(fig, output, "71-path-decomposition")


def model_fit_dsep(input_dir: Path, output: Path) -> None:
    models = pd.read_csv(input_dir / "sem-fit-comparison.tsv", sep="\t")
    dsep = pd.read_csv(input_dir / "directed-separation-claims.tsv", sep="\t").iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))

    axis = axes[0]
    names = ["Partial", "Microbiome-only", "Reverse"]
    colors = [COLORS["Exposure"], COLORS["Microbiome"], COLORS["Purple"]]
    bars = axis.bar(names, models["DeltaAIC"], color=colors, alpha=0.85)
    axis.set_ylabel("Delta AIC")
    axis.set_title("Direction is not selected by AIC here", loc="left", fontweight="bold")
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    for bar, row in zip(bars, models.itertuples(), strict=True):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.04,
            f"AIC {row.AIC:.1f}",
            ha="center",
            fontsize=7.6,
        )
    axis.set_ylim(0, max(models["DeltaAIC"].max() + 0.35, 1.5))
    panel_label(axis, "A")

    axis = axes[1]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    add_box(axis, (0.25, 0.74), "Partial path\n0 d-sep claims", "#F7D6CE", 0.32, 0.16)
    add_box(
        axis,
        (0.75, 0.74),
        "Microbiome-only\n1 d-sep claim",
        "#D4ECE7",
        0.32,
        0.16,
    )
    axis.text(
        0.25,
        0.48,
        "Saturated\nFisher C unavailable\n(df = 0)",
        ha="center",
        va="center",
        color="#8B3D00",
        fontweight="bold",
    )
    axis.text(
        0.75,
        0.48,
        f"Fisher C = 4.66\nP = 0.098\nomitted direct path P = {dsep['PValue']:.3f}",
        ha="center",
        va="center",
        color="#176B5B",
        fontweight="bold",
    )
    axis.text(
        0.50,
        0.17,
        "Not rejected is not validated:\n13 exposed subjects give little power.",
        ha="center",
        va="center",
        color="#59636A",
    )
    axis.set_title("Global fit exists only when the graph omits a path", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.1)
    save(fig, output, "71-model-fit-dsep")


def overlap_influence(input_dir: Path, output: Path) -> None:
    overlap = pd.read_csv(
        input_dir / "antibiotic-overlap-by-diagnosis.tsv", sep="\t"
    )
    leave = pd.read_csv(input_dir / "leave-one-out-paths.tsv", sep="\t")
    metrics = pd.read_json(input_dir / "model-metrics.json", typ="series")
    fig, axes = plt.subplots(1, 2, figsize=(10.7, 4.25))

    axis = axes[0]
    bars = axis.bar(
        overlap["Diagnosis"],
        100 * overlap["ExposureFraction"],
        color=[COLORS.get(x, COLORS["Neutral"]) for x in overlap["Diagnosis"]],
    )
    axis.set_ylabel("Antibiotic-exposed subjects (%)")
    axis.set_title("Structural positivity fails in controls", loc="left", fontweight="bold")
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    for bar, row in zip(bars, overlap.itertuples(), strict=True):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.8,
            f"{int(row.Exposed)}/{int(row.Exposed + row.Unexposed)}",
            ha="center",
        )
    axis.text(
        0.02,
        0.96,
        "Validation exposed: 0/38\nPropensity max |coef|: 19.4",
        transform=axis.transAxes,
        va="top",
        fontsize=7.5,
        color="#8B3D00",
    )
    panel_label(axis, "A")

    axis = axes[1]
    ordered = leave.sort_values("Indirect").reset_index(drop=True)
    axis.plot(
        np.arange(len(ordered)),
        ordered["Indirect"],
        color=COLORS["Microbiome"],
        linewidth=1.4,
    )
    axis.scatter(
        np.arange(len(ordered)),
        ordered["Indirect"],
        s=np.where(ordered["OmittedAntibiotic"].eq(1), 34, 16),
        color=np.where(
            ordered["OmittedAntibiotic"].eq(1),
            COLORS["Exposure"],
            COLORS["Microbiome"],
        ),
        alpha=0.82,
    )
    axis.axhline(
        metrics["shannon_indirect"],
        color="#333333",
        linestyle="--",
        linewidth=0.9,
        label="All-subject indirect path",
    )
    axis.axhline(0, color="#777777", linewidth=0.7)
    axis.set_xlabel("Omitted-subject analyses, ordered by estimate")
    axis.set_ylabel("Indirect path effect")
    axis.set_title("Leave-one-out stability is not positivity", loc="left", fontweight="bold")
    axis.legend(frameon=False)
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    axis.text(
        0.02,
        0.03,
        "Orange: omitted exposed subject",
        transform=axis.transAxes,
        fontsize=7.4,
        color="#8B3D00",
    )
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.1)
    save(fig, output, "71-overlap-influence")


def transport_sensitivity(input_dir: Path, output: Path) -> None:
    transport = pd.read_csv(input_dir / "outcome-path-transport.tsv", sep="\t")
    effects = pd.read_csv(input_dir / "path-effect-summary.tsv", sep="\t")
    indirect = effects.loc[effects["Effect"].eq("Indirect")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))

    axis = axes[0]
    y = np.arange(len(transport))[::-1]
    labels = ["PRISM (n=90)", "Validation (n=38)"]
    transport_colors = [COLORS["Exposure"], COLORS["Phenotype"]]
    for yy, row, color in zip(y, transport.itertuples(), transport_colors, strict=True):
        axis.errorbar(
            row.Estimate,
            yy,
            xerr=[[row.Estimate - row.CILower], [row.CIUpper - row.Estimate]],
            fmt="o",
            color=color,
            ecolor=color,
            markersize=6,
            elinewidth=1.8,
            capsize=3,
            zorder=3,
        )
    axis.axvline(0, color="#555555", linestyle="--", linewidth=0.9)
    axis.set_yticks(y, labels)
    axis.set_xlabel("Shannon-to-calprotectin coefficient")
    axis.set_title("Outcome path does not transport in direction", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    axis.text(
        0.02,
        0.50,
        "Same covariate set and PRISM scaling;\nValidation has no exposed subject.",
        transform=axis.transAxes,
        fontsize=7.4,
        color="#59636A",
        va="center",
    )
    panel_label(axis, "A")

    axis = axes[1]
    y = np.arange(len(indirect))[::-1]
    colors = [COLORS["Microbiome"], COLORS["Purple"]]
    labels = ["Shannon", "Faecalibacterium"]
    for yy, row, color in zip(y, indirect.itertuples(), colors, strict=True):
        axis.errorbar(
            row.Estimate,
            yy,
            xerr=[[row.Estimate - row.CILower], [row.CIUpper - row.Estimate]],
            fmt="o",
            color=color,
            ecolor=color,
            markersize=6,
            elinewidth=1.8,
            capsize=3,
            zorder=3,
        )
    axis.axvline(0, color="#555555", linestyle="--", linewidth=0.9)
    axis.set_yticks(y, labels)
    axis.set_xlabel("Antibiotic-to-microbiome-to-calprotectin path")
    axis.set_title("Mediator definition changes the number", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.1)
    save(fig, output, "71-transport-sensitivity")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output = args.figure_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    np.random.seed(PLOT_SEED)
    style()
    anchor = input_dir / "franzosa-fig1-original.png"
    shutil.copy2(anchor, output / "71-franzosa-fig1-original.png")
    data_positivity(input_dir, output)
    prespecified_dag(input_dir, output)
    local_paths(input_dir, output)
    path_decomposition(input_dir, output)
    model_fit_dsep(input_dir, output)
    overlap_influence(input_dir, output)
    transport_sensitivity(input_dir, output)


if __name__ == "__main__":
    main()
