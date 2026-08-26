#!/usr/bin/env python3
"""Create publication-ready, English-only figures for Article 61."""

from __future__ import annotations

import argparse
import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/article61-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Patch


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


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


def short_taxon(value: str) -> str:
    words = str(value).replace("_", " ").split()
    if len(words) < 2:
        return " ".join(words)
    return f"{words[0][0]}. {' '.join(words[1:])}"


def wrap(value: str, width: int = 24) -> str:
    return "\n".join(textwrap.wrap(str(value), width=width, break_long_words=False))


def signed_log(values: pd.DataFrame | pd.Series) -> pd.DataFrame | pd.Series:
    return np.sign(values) * np.log10(1.0 + np.abs(values))


def plot_coverage(summary: Path, figures: Path) -> None:
    coverage = (
        read_tsv(summary / "model-coverage.tsv")
        .sort_values("SubjectID")
        .reset_index(drop=True)
    )
    categories = pd.DataFrame(
        {
            "Modeled abundance": coverage["ModeledAbundance"].to_numpy(),
            "Matched below cutoff": (
                coverage["AllMatchedAbundance"] - coverage["ModeledAbundance"]
            ).clip(lower=0).to_numpy(),
            "Unmatched species": (
                coverage["SpeciesResolvedAbundance"] - coverage["AllMatchedAbundance"]
            ).clip(lower=0).to_numpy(),
            "Not species-resolved": (
                1.0 - coverage["SpeciesResolvedAbundance"]
            ).clip(lower=0).to_numpy(),
        },
        index=coverage["SubjectID"].to_numpy(),
    )
    colors = [GREEN, SKY, ORANGE, GRAY]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2), gridspec_kw={"width_ratios": [1.45, 1]})
    left = np.zeros(len(categories))
    y = np.arange(len(categories))
    for column, color in zip(categories.columns, colors):
        values = categories[column].to_numpy() * 100
        axes[0].barh(y, values, left=left, height=0.62, color=color, label=column)
        left += values
    axes[0].axvline(50, color=DARK, linestyle=":", linewidth=1)
    axes[0].set_yticks(y, categories.index)
    axes[0].set_xlim(0, 100)
    axes[0].set_xlabel("Whole-profile relative abundance (%)")
    axes[0].set_ylabel("Independent adult subject")
    axes[0].set_title("A. Model coverage retains the original denominator")
    axes[0].legend(loc="lower center", bbox_to_anchor=(0.5, -0.34), ncol=2)

    reads = coverage["Reads"].to_numpy() / 1e6
    axes[1].hlines(y, 0, reads, color="#B0BEC5", linewidth=2.2)
    axes[1].scatter(reads, y, color=BLUE, s=48, zorder=3)
    axes[1].axvline(1, color=ORANGE, linestyle="--", linewidth=1, label="1 M read gate")
    for index, value in enumerate(reads):
        axes[1].text(value + 0.45, index, f"{value:.1f}", va="center", fontsize=7.5)
    axes[1].set_yticks(y, categories.index)
    axes[1].set_xlabel("Reported whole-metagenome reads (millions)")
    axes[1].set_ylabel("")
    axes[1].set_title("B. Prespecified read-depth gate")
    axes[1].legend(loc="lower right")
    sns.despine(fig=fig)
    fig.suptitle(
        "Six independent AsnicarF_2017 adult profiles enter personalized community models",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    save_pub(fig, figures, "61-model-coverage")


def plot_tradeoff(summary: Path, figures: Path) -> None:
    frame = read_tsv(summary / "tradeoff-summary.tsv").sort_values(["SubjectID", "tradeoff"])
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.0))
    for subject, group in frame.groupby("SubjectID", sort=True):
        axes[0].plot(group["tradeoff"], group["FractionGrowing"], color=SKY, alpha=0.45, linewidth=1.1)
        axes[1].plot(group["tradeoff"], group["CommunityGrowth"], color=PURPLE, alpha=0.45, linewidth=1.1)
    median = frame.groupby("tradeoff", as_index=False)[["FractionGrowing", "CommunityGrowth"]].median()
    axes[0].plot(median["tradeoff"], median["FractionGrowing"], color=BLUE, marker="o", linewidth=2.4, label="Median")
    axes[1].plot(median["tradeoff"], median["CommunityGrowth"], color=PURPLE, marker="o", linewidth=2.4, label="Median")
    for ax in axes:
        ax.axvline(0.5, color=ORANGE, linestyle="--", linewidth=1.1, label="Primary trade-off = 0.5")
        ax.set_xlabel("Cooperative trade-off")
        ax.set_xlim(0.08, 0.92)
        ax.legend(loc="best")
    axes[0].set_ylabel("Fraction of modeled taxa with growth > 10⁻⁶")
    axes[0].set_ylim(-0.03, 1.03)
    axes[0].set_title("A. Taxon coexistence changes with the objective")
    axes[1].set_ylabel("Abundance-weighted biomass objective")
    axes[1].set_title("B. Community objective retained")
    sns.despine(fig=fig)
    fig.suptitle("Trade-off selection is a modeling decision, not a fitted biological constant", fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_pub(fig, figures, "61-tradeoff-curves")


def plot_sensitivity(summary: Path, figures: Path) -> None:
    medium = read_tsv(summary / "medium-sensitivity-summary.tsv")
    baseline = medium.loc[np.isclose(medium["MediumScale"], 1.0), ["sample_id", "CommunityGrowth", "FractionGrowing"]].rename(
        columns={"CommunityGrowth": "BaselineGrowth", "FractionGrowing": "BaselineFraction"}
    )
    medium = medium.merge(baseline, on="sample_id", validate="many_to_one")
    medium["GrowthRatio"] = medium["CommunityGrowth"] / medium["BaselineGrowth"].replace(0, np.nan)
    abundance = read_tsv(summary / "abundance-sensitivity-summary.tsv")

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.0))
    for subject, group in medium.groupby("SubjectID", sort=True):
        group = group.sort_values("MediumScale")
        axes[0].plot(group["MediumScale"], group["GrowthRatio"], color=SKY, alpha=0.65, marker="o", linewidth=1.3)
    med = medium.groupby("MediumScale", as_index=False)["GrowthRatio"].median()
    axes[0].plot(med["MediumScale"], med["GrowthRatio"], color=BLUE, marker="o", linewidth=2.5, label="Median")
    axes[0].axhline(1, color=GRAY, linestyle=":", linewidth=1)
    axes[0].set_xticks([0.5, 1.0, 2.0], ["0.5×", "1×", "2×"])
    axes[0].set_xlabel("All Western-diet uptake bounds scaled together")
    axes[0].set_ylabel("Community objective / 1× objective")
    axes[0].set_title("A. Medium-bound sensitivity")
    axes[0].legend()

    minimum = min(abundance["CommunityGrowthObserved"].min(), abundance["CommunityGrowthEqual"].min())
    maximum = max(abundance["CommunityGrowthObserved"].max(), abundance["CommunityGrowthEqual"].max())
    axes[1].plot([minimum, maximum], [minimum, maximum], color=GRAY, linestyle="--", linewidth=1, label="No change")
    axes[1].scatter(
        abundance["CommunityGrowthObserved"], abundance["CommunityGrowthEqual"],
        s=55, color=GREEN, edgecolor=WHITE, linewidth=0.7,
    )
    for row in abundance.itertuples(index=False):
        axes[1].annotate(row.SubjectID, (row.CommunityGrowthObserved, row.CommunityGrowthEqual), xytext=(3, 3), textcoords="offset points", fontsize=7)
    axes[1].set_xlabel("Observed-abundance community objective")
    axes[1].set_ylabel("Equal-abundance community objective")
    axes[1].set_title("B. Composition-weight sensitivity")
    axes[1].legend()
    sns.despine(fig=fig)
    fig.suptitle("Predictions are conditional on both environmental and abundance constraints", fontsize=12, fontweight="bold")
    fig.tight_layout()
    save_pub(fig, figures, "61-constraint-sensitivity")


def plot_net_flux(summary: Path, figures: Path) -> None:
    frame = read_tsv(summary / "net-community-flux.tsv")
    all_subjects = sorted(frame.SubjectID.unique())
    all_fluxes = frame.pivot_table(
        index="CompoundName",
        columns="SubjectID",
        values="CommunityScaledFlux",
        aggfunc="sum",
        fill_value=0,
    ).reindex(columns=all_subjects, fill_value=0)
    ranking = (
        all_fluxes.abs().median(axis=1).sort_values(ascending=False).head(16).index
    )
    matrix = all_fluxes.loc[list(ranking)[::-1]]
    transformed = signed_log(matrix)
    limit = float(np.nanmax(np.abs(transformed.to_numpy()))) or 1.0
    cmap = LinearSegmentedColormap.from_list("signed_flux", [BLUE, WHITE, RED])
    fig, ax = plt.subplots(figsize=(10.8, 7.0))
    sns.heatmap(
        transformed,
        cmap=cmap,
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
        linewidths=0.5,
        linecolor=WHITE,
        cbar_kws={"label": "Signed log10(1 + |flux|): import ← 0 → export"},
        ax=ax,
    )
    ax.set_xlabel("Independent adult subject")
    ax.set_ylabel("Modeled extracellular metabolite")
    ax.set_title("Net community exchange under the locked Western-diet constraints")
    fig.tight_layout()
    save_pub(fig, figures, "61-net-community-flux")


def plot_micom_crossfeeding(summary: Path, figures: Path) -> None:
    flux = read_tsv(summary / "focal-micom-fluxes.tsv")
    edges = read_tsv(summary / "focal-micom-potential-edges.tsv")
    if edges.empty:
        ranking = flux.groupby("CompoundName")["CommunityScaledFlux"].apply(lambda x: x.abs().sum()).sort_values(ascending=False).head(16).index
    else:
        ranking = edges.groupby("CompoundName")["PotentialFlux"].sum().sort_values(ascending=False).head(16).index
    subset = flux.loc[flux["CompoundName"].isin(ranking)].copy()
    subset["TaxonLabel"] = subset["taxon"].map(short_taxon)
    matrix = subset.pivot_table(index="CompoundName", columns="TaxonLabel", values="CommunityScaledFlux", aggfunc="sum", fill_value=0)
    matrix = matrix.reindex(list(ranking)[::-1]).fillna(0)
    transformed = signed_log(matrix)
    limit = float(np.nanmax(np.abs(transformed.to_numpy()))) or 1.0
    cmap = LinearSegmentedColormap.from_list("taxon_flux", [BLUE, WHITE, RED])
    fig, ax = plt.subplots(figsize=(11.8, 7.0))
    sns.heatmap(
        transformed,
        cmap=cmap,
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
        linewidths=0.45,
        linecolor=WHITE,
        cbar_kws={"label": "Abundance-scaled signed log10(1 + |flux|): import ← 0 → export"},
        ax=ax,
    )
    ax.set_xlabel("Focal community member")
    ax.set_ylabel("Candidate exchanged metabolite")
    ax.set_title("MICOM pFBA: simultaneous export and import indicate cross-feeding potential")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right")
    fig.tight_layout()
    save_pub(fig, figures, "61-micom-crossfeeding")


def plot_smetana(summary: Path, figures: Path) -> None:
    audit = read_tsv(summary / "smetana-compatibility-audit.tsv")
    components = read_tsv(summary / "smetana-component-summary.tsv")
    import json
    metrics = json.loads((summary / "analysis-metrics.json").read_text(encoding="utf-8"))

    modes = {
        "StandaloneCompleteGrowth": "Standalone model",
        "InteractingCompleteGrowth": "Interacting merge",
        "NoninteractingCompleteGrowth": "Legacy non-interacting merge",
    }
    long = audit.melt(
        id_vars="ModelID",
        value_vars=list(modes),
        var_name="Mode",
        value_name="MaximumGrowth",
    )
    long["Mode"] = long.Mode.map(modes)
    long["Taxon"] = long.ModelID.map(short_taxon)
    long["LogGrowth"] = np.log10(1.0 + long.MaximumGrowth.clip(lower=0))

    fig, axes = plt.subplots(
        1, 3, figsize=(16.0, 5.4), gridspec_kw={"width_ratios": [1.45, 0.85, 1.05]}
    )
    sns.barplot(
        data=long,
        x="Taxon",
        y="LogGrowth",
        hue="Mode",
        palette=[GREEN, SKY, ORANGE],
        ax=axes[0],
    )
    axes[0].set_xlabel("Exported AGORA2 species model")
    axes[0].set_ylabel("log10(1 + maximum growth)")
    axes[0].set_title("A. Complete-environment structural control")
    axes[0].tick_params(axis="x", rotation=40)
    for label in axes[0].get_xticklabels():
        label.set_horizontalalignment("right")
    axes[0].legend(title="Model context", loc="upper right")

    colors = [GRAY, BLUE, GREEN, ORANGE]
    axes[1].bar(
        np.arange(len(components)), components.PositiveRows, color=colors
    )
    axes[1].set_xticks(np.arange(len(components)), components.Component)
    axes[1].set_ylabel("Detailed rows with score > 0")
    axes[1].set_title("B. Detailed components returned")
    for index, value in enumerate(components.PositiveRows):
        axes[1].text(index, value, f"{int(value):,}", ha="center", va="bottom", fontsize=8)

    axes[2].set_axis_off()
    status_rows = [
        ("SBML identity convention", "6 / 6 passed", GREEN),
        ("Medium-to-pool match", f"{metrics['smetana_medium_matches']} / 159", GREEN),
        ("Detailed enumeration", f"{metrics['smetana_detailed_rows']:,} rows", BLUE),
        ("Global MIP / MRO", "Not estimable", ORANGE),
        ("Cross-method overlap", "Not estimable", ORANGE),
    ]
    axes[2].text(0.0, 1.0, "C. Fail-closed interpretation", fontsize=10.5, fontweight="bold", va="top")
    for index, (label, value, color) in enumerate(status_rows):
        y = 0.84 - index * 0.14
        axes[2].text(0.0, y, label, color=DARK, va="center", fontsize=9)
        axes[2].text(
            0.98, y, value, color=color, va="center", ha="right",
            fontsize=9, fontweight="bold",
        )
        axes[2].plot([0, 1], [y - 0.06, y - 0.06], color="#CFD8DC", linewidth=0.7)
    axes[2].text(
        0.0, 0.05,
        "Software/model-interface limitation;\nnot evidence that biological exchange is absent.",
        color=RED, fontsize=9, fontweight="bold", va="bottom",
    )
    sns.despine(fig=fig)
    fig.suptitle(
        "SMETANA 1.2.1 × exported AGORA2 requires an explicit compatibility audit",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    save_pub(fig, figures, "61-smetana-concordance")


def main() -> None:
    args = parse_args()
    setup_style()
    plot_coverage(args.summary_dir, args.figure_dir)
    plot_tradeoff(args.summary_dir, args.figure_dir)
    plot_sensitivity(args.summary_dir, args.figure_dir)
    plot_net_flux(args.summary_dir, args.figure_dir)
    plot_micom_crossfeeding(args.summary_dir, args.figure_dir)
    plot_smetana(args.summary_dir, args.figure_dir)
    print("Article 61 publication figures: PASS")


if __name__ == "__main__":
    main()
