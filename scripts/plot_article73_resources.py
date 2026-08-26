#!/usr/bin/env python3
"""Create deterministic publication figures for Article 73."""

from __future__ import annotations

import argparse
import shutil
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PLOT_SEED = 20_260_773
COLORS = {
    "Archive": "#4C78A8",
    "Analysis": "#76A5AF",
    "Catalogue": "#59A14F",
    "Taxonomy": "#B279A2",
    "DNA": "#2A9D8F",
    "RNA": "#D55E00",
    "HQ": "#2A9D8F",
    "MQ": "#9AA5B1",
    "Risk": "#C44E52",
    "Ink": "#33434A",
    "Grid": "#E6EAED",
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
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.08,
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


def add_box(
    axis: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    face: str,
    edge: str = "#66737A",
    fontsize: float = 7.5,
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.015",
        facecolor=face,
        edgecolor=edge,
        linewidth=0.9,
    )
    axis.add_patch(patch)
    axis.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
    )


def resource_layers(input_dir: Path, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.3), gridspec_kw={"width_ratios": [1.25, 0.75]})
    axis = axes[0]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    boxes = [
        (0.03, "ENA\nStudy · BioSample\nRun · FASTQ", COLORS["Archive"]),
        (0.28, "MGnify\nMGYS · MGYA\nversioned output", COLORS["Analysis"]),
        (0.53, "GEM / biome catalogue\nMAG · representative\nprotein catalogue", COLORS["Catalogue"]),
        (0.78, "GTDB R232\nrelease-bound\ngenome taxonomy", COLORS["Taxonomy"]),
    ]
    for index, (x, label, color) in enumerate(boxes):
        add_box(axis, x, 0.57, 0.19, 0.25, label, color, fontsize=8.0)
        if index < len(boxes) - 1:
            axis.add_patch(
                FancyArrowPatch(
                    (x + 0.195, 0.695),
                    (boxes[index + 1][0] - 0.008, 0.695),
                    arrowstyle="-|>",
                    mutation_scale=11,
                    linewidth=1.1,
                    color="#59636A",
                )
            )
    notes = [
        "Archive identity",
        "Pipeline identity",
        "Catalogue identity",
        "Taxonomy identity",
    ]
    for (x, _, color), note in zip(boxes, notes, strict=True):
        axis.text(x + 0.095, 0.47, note, ha="center", fontweight="bold", color=color)
    axis.text(
        0.50,
        0.27,
        "Never replace an accession crosswalk with a filename.\nNever replace a release with the word ‘latest’.",
        ha="center",
        color=COLORS["Risk"],
        fontweight="bold",
        fontsize=9,
    )
    axis.text(
        0.50,
        0.10,
        "The four layers answer different questions and can change independently.",
        ha="center",
        color=COLORS["Ink"],
    )
    axis.set_title("Public metagenome resources form a lineage—not one database", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    decisions = [
        ("Reanalyse reads", "ENA runs + FASTQ", COLORS["Archive"]),
        ("Reuse annotations", "Exact MGYA download", COLORS["Analysis"]),
        ("Screen novelty", "Biome-matched representatives", COLORS["Catalogue"]),
        ("Name your MAGs", "Pinned GTDB release", COLORS["Taxonomy"]),
    ]
    for index, (question, answer, color) in enumerate(decisions):
        y = 0.80 - index * 0.20
        axis.text(0.04, y + 0.055, question, fontweight="bold", color=color)
        add_box(axis, 0.04, y - 0.055, 0.91, 0.10, answer, "#F1F4F5", edge=color, fontsize=7.8)
    axis.text(0.50, 0.03, "Download only the layer needed for the estimand.", ha="center", color=COLORS["Risk"], fontweight="bold")
    axis.set_title("Start from the question", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "73-resource-layer-map")


def gem_biome_balance(input_dir: Path, output: Path) -> None:
    data = pd.read_csv(input_dir / "gem-biome-summary.tsv", sep="\t").head(12).copy()
    data = data.sort_values("MAGs")
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.6), gridspec_kw={"width_ratios": [1.15, 0.85]})
    y = np.arange(len(data))
    axis = axes[0]
    colors = [COLORS["Catalogue"] if name in {"Aquatic", "Terrestrial"} else "#9AA5B1" for name in data["EcosystemCategory"]]
    bars = axis.barh(y, data["MAGs"], color=colors, height=0.68)
    axis.set_yticks(y, data["EcosystemCategory"])
    axis.set_xlabel("GEM MAGs")
    axis.grid(axis="x", color=COLORS["Grid"], linewidth=0.55)
    for bar, value in zip(bars, data["MAGs"], strict=True):
        axis.text(bar.get_width() + 180, bar.get_y() + bar.get_height() / 2, f"{value:,}", va="center", fontsize=7.2)
    axis.set_title("The fixed GEM snapshot is not biome-balanced", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    axis.scatter(data["MeanCompleteness"], data["HighQualityPct"], s=np.sqrt(data["MAGs"]) * 7, c=COLORS["HQ"], alpha=0.78, edgecolor="white", linewidth=0.6)
    for row in data.itertuples():
        if row.EcosystemCategory in {"Aquatic", "Human", "Terrestrial", "Built environment", "Wastewater"}:
            axis.text(row.MeanCompleteness + 0.18, row.HighQualityPct + 0.3, row.EcosystemCategory, fontsize=7.0)
    axis.set_xlabel("Mean completeness (%)")
    axis.set_ylabel("High-quality MAGs (%)")
    axis.grid(color=COLORS["Grid"], linewidth=0.55)
    axis.set_title("Recovery quality also differs by source", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.suptitle("52,515 published GEM MAGs require an explicit sampling-frame caveat", x=0.02, ha="left", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=2.2)
    save(fig, output, "73-gem-biome-balance")


def geo_axes(axis: plt.Axes) -> None:
    axis.set_xlim(-180, 180)
    axis.set_ylim(-90, 90)
    axis.set_xticks(np.arange(-180, 181, 60))
    axis.set_yticks(np.arange(-90, 91, 30))
    axis.set_xlabel("Longitude")
    axis.set_ylabel("Latitude")
    axis.grid(color="#DDE3E6", linewidth=0.45)
    axis.axhline(0, color="#AAB3B8", linewidth=0.65)
    axis.axvline(0, color="#AAB3B8", linewidth=0.65)


def gem_world_map(input_dir: Path, output: Path) -> None:
    data = pd.read_csv(input_dir / "gem-map-visualization-sample.tsv", sep="\t")
    top = data["EcosystemCategory"].value_counts().head(5).index
    palette = {
        "Aquatic": "#4C78A8",
        "Human": "#E76F51",
        "Terrestrial": "#59A14F",
        "Built environment": "#B279A2",
        "Wastewater": "#F2C14E",
        "Other": "#A9B0B4",
    }
    data["DisplayCategory"] = data["EcosystemCategory"].where(data["EcosystemCategory"].isin(top), "Other")
    fig, axes = plt.subplots(2, 1, figsize=(11.7, 7.2), gridspec_kw={"height_ratios": [1.25, 0.75]})
    axis = axes[0]
    for category in ["Other"] + list(top[::-1]):
        subset = data.loc[data["DisplayCategory"].eq(category)]
        axis.scatter(subset["Longitude"], subset["Latitude"], s=5, alpha=0.35 if category == "Other" else 0.58, label=category, color=palette.get(category, "#A9B0B4"), linewidth=0)
    geo_axes(axis)
    axis.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.27))
    axis.set_title("Georeferenced recovery is globally broad—but visibly uneven", loc="left", fontweight="bold")
    axis.text(0.99, 0.03, "12,000-row deterministic display subset", transform=axis.transAxes, ha="right", color="#59636A")
    panel_label(axis, "A")

    axis = axes[1]
    counts = data.groupby(["DisplayCategory", "MIMAGQuality"]).size().unstack(fill_value=0)
    order = counts.sum(axis=1).sort_values().index
    counts = counts.loc[order]
    y = np.arange(len(counts))
    mq = counts.get("MQ", pd.Series(0, index=counts.index))
    hq = counts.get("HQ", pd.Series(0, index=counts.index))
    axis.barh(y, mq, color=COLORS["MQ"], label="Medium quality")
    axis.barh(y, hq, left=mq, color=COLORS["HQ"], label="High quality")
    axis.set_yticks(y, counts.index)
    axis.set_xlabel("Displayed georeferenced MAGs")
    axis.grid(axis="x", color=COLORS["Grid"], linewidth=0.55)
    axis.legend(frameon=False, loc="lower right")
    axis.set_title("Map density mixes sampling effort and recovery success", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(h_pad=2.2)
    save(fig, output, "73-gem-geographic-coverage")


def gem_quality(input_dir: Path, output: Path) -> None:
    data = pd.read_csv(input_dir / "gem-quality-visualization-sample.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.0))
    axis = axes[0]
    for label, color, alpha in (("MQ", COLORS["MQ"], 0.28), ("HQ", COLORS["HQ"], 0.55)):
        subset = data.loc[data["MIMAGQuality"].eq(label)]
        axis.scatter(subset["Contamination"], subset["Completeness"], s=8, alpha=alpha, color=color, label={"MQ": "Medium quality", "HQ": "High quality"}[label], linewidth=0)
    x = np.linspace(0, 5, 100)
    axis.plot(x, 50 + 5 * x, color=COLORS["Risk"], linewidth=1.2, linestyle="--", label="Quality score = 50")
    axis.axvline(5, color=COLORS["Risk"], linewidth=0.8, linestyle=":")
    axis.set_xlim(-0.1, 5.15)
    axis.set_ylim(48, 101)
    axis.set_xlabel("Contamination (%)")
    axis.set_ylabel("Completeness (%)")
    axis.grid(color=COLORS["Grid"], linewidth=0.5)
    axis.legend(frameon=False, loc="lower right")
    axis.set_title("All released GEM MAGs pass the catalogue floor", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    order = data["EcosystemCategory"].value_counts().index.tolist()
    values = [data.loc[data["EcosystemCategory"].eq(category), "QualityScore"].to_numpy() for category in order]
    parts = axis.violinplot(values, positions=np.arange(len(order)), showmeans=False, showmedians=True, widths=0.8)
    for body in parts["bodies"]:
        body.set_facecolor(COLORS["Catalogue"])
        body.set_edgecolor("white")
        body.set_alpha(0.66)
    parts["cmedians"].set_color(COLORS["Ink"])
    for key in ("cbars", "cmins", "cmaxes"):
        parts[key].set_color("#78868D")
    axis.set_xticks(np.arange(len(order)), order, rotation=30, ha="right")
    axis.set_ylabel("Quality score = completeness − 5 × contamination")
    axis.grid(axis="y", color=COLORS["Grid"], linewidth=0.5)
    axis.set_title("Passing one floor does not erase biome structure", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.1)
    save(fig, output, "73-gem-quality-audit")


def tara_sampling(input_dir: Path, output: Path) -> None:
    samples = pd.read_csv(input_dir / "tara-samples.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw={"width_ratios": [1.25, 0.75]})
    axis = axes[0]
    for label in ("RNA", "DNA"):
        subset = samples.loc[samples["NucleicAcid"].eq(label)]
        axis.scatter(subset["Longitude"], subset["Latitude"], s=18, alpha=0.72, color=COLORS[label], label=f"{label}-labelled samples (n={len(subset)})", edgecolor="white", linewidth=0.3)
    geo_axes(axis)
    axis.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=2)
    axis.set_title("MGYS00000410 relationship metadata spans DNA and RNA protocols", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    bins = [0, 10, 50, 200, 500, 1000, 4000]
    for label in ("DNA", "RNA"):
        subset = samples.loc[samples["NucleicAcid"].eq(label), "DepthM"]
        axis.hist(subset, bins=bins, alpha=0.58, color=COLORS[label], label=label, edgecolor="white")
    axis.set_xscale("symlog", linthresh=10)
    axis.set_xlabel("Sampling depth (m; symlog scale)")
    axis.set_ylabel("Samples")
    axis.grid(axis="y", color=COLORS["Grid"], linewidth=0.5)
    axis.legend(frameon=False)
    axis.set_title("Depth and protocol define eligibility", loc="left", fontweight="bold")
    axis.text(0.02, 0.04, "Study title alone is not a cohort filter", transform=axis.transAxes, color=COLORS["Risk"], fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "73-tara-sampling-frame")


def metadata_completeness(input_dir: Path, output: Path) -> None:
    data = pd.read_csv(input_dir / "tara-metadata-completeness.tsv", sep="\t").sort_values("CompletenessPct")
    fig, axis = plt.subplots(figsize=(10.8, 5.2))
    colors = [COLORS["Risk"] if value == 0 else COLORS["Analysis"] for value in data["CompletenessPct"]]
    bars = axis.barh(np.arange(len(data)), data["CompletenessPct"], color=colors, height=0.64)
    axis.set_yticks(np.arange(len(data)), data["Field"])
    axis.set_xlim(0, 108)
    axis.set_xlabel("Metadata completeness (%)")
    axis.grid(axis="x", color=COLORS["Grid"], linewidth=0.55)
    for bar, row in zip(bars, data.itertuples(), strict=True):
        axis.text(max(bar.get_width() + 1.2, 2), bar.get_y() + bar.get_height() / 2, f"{row.Available}/{row.Total}", va="center", fontweight="bold", fontsize=7.6)
    axis.text(0.99, 0.04, "Missing in this API snapshot does not mean unmeasured in the original expedition.", transform=axis.transAxes, ha="right", color=COLORS["Risk"], fontweight="bold")
    axis.set_title("Environmental covariates must be audited before biological filtering", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, output, "73-tara-metadata-completeness")


def accession_crosswalk(input_dir: Path, output: Path) -> None:
    crosswalk = pd.read_csv(input_dir / "accession-crosswalk.tsv", sep="\t")
    counts = dict(zip(crosswalk["Object"], crosswalk["Count"], strict=True))
    fig, axis = plt.subplots(figsize=(11.8, 5.4))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    steps = [
        (0.03, 0.63, 0.17, "MGnify study\nMGYS00000410\n1", COLORS["Analysis"]),
        (0.25, 0.63, 0.17, "ENA projects\nPRJEB1787 / ERP001736\n2 IDs", COLORS["Archive"]),
        (0.47, 0.63, 0.17, f"BioSamples\nSAMEA / ERS\n{counts['Related BioSamples']}", "#DCE7F2"),
        (0.69, 0.63, 0.13, f"Runs\nERR\n{counts['Sequencing runs']}", "#DCE7F2"),
        (0.86, 0.63, 0.11, f"Analyses\nMGYA\n{counts['MGnify analyses']}", COLORS["Analysis"]),
    ]
    for index, (x, y, width, label, color) in enumerate(steps):
        add_box(axis, x, y, width, 0.20, label, color, fontsize=7.6)
        if index < len(steps) - 1:
            next_x = steps[index + 1][0]
            axis.add_patch(FancyArrowPatch((x + width + 0.004, y + 0.10), (next_x - 0.008, y + 0.10), arrowstyle="-|>", mutation_scale=10, color="#59636A"))
    add_box(axis, 0.37, 0.26, 0.23, 0.16, f"DNA-labelled\n{counts['DNA-labelled samples']} samples", "#DDEFEA", edge=COLORS["DNA"], fontsize=8.0)
    add_box(axis, 0.63, 0.26, 0.23, 0.16, f"RNA-labelled\n{counts['RNA-labelled samples']} samples", "#F8E4D8", edge=COLORS["RNA"], fontsize=8.0)
    axis.add_patch(FancyArrowPatch((0.555, 0.63), (0.485, 0.43), arrowstyle="-|>", mutation_scale=10, color=COLORS["DNA"]))
    axis.add_patch(FancyArrowPatch((0.575, 0.63), (0.745, 0.43), arrowstyle="-|>", mutation_scale=10, color=COLORS["RNA"]))
    axis.text(0.50, 0.12, "One study title · 136 samples · 249 runs/analyses · two nucleic-acid protocols", ha="center", fontweight="bold", color=COLORS["Risk"], fontsize=9)
    axis.text(0.50, 0.05, "The biological cohort is created by explicit eligibility rules, not by accession prefix.", ha="center", color=COLORS["Ink"])
    axis.set_title("Keep the accession lineage until the final sample sheet", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, output, "73-accession-crosswalk")


def catalogue_landscape(input_dir: Path, output: Path) -> None:
    catalogues = pd.read_csv(input_dir / "mgnify-catalogues.tsv", sep="\t").head(10).copy().sort_values("InputGenomes")
    history = pd.read_csv(input_dir / "gtdb-release-history.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.6), gridspec_kw={"width_ratios": [1.25, 0.75]})
    axis = axes[0]
    y = np.arange(len(catalogues))
    axis.barh(y, catalogues["InputGenomes"], color="#DDE3E6", label="Input genomes")
    axis.barh(y, catalogues["RepresentativeClusters"], color=COLORS["Catalogue"], label="Representative clusters")
    axis.set_yticks(y, catalogues["Name"])
    axis.set_xscale("log")
    axis.set_xlabel("Catalogue rows (log scale)")
    axis.grid(axis="x", color=COLORS["Grid"], linewidth=0.55)
    axis.legend(frameon=False, loc="lower right")
    axis.set_title("MGnify catalogues have distinct inputs and dereplication scopes", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    x = np.arange(len(history))
    width = 0.34
    axis.bar(x - width / 2, history["Genomes"] / 1e3, width, color=COLORS["Taxonomy"], label="Genomes")
    axis.bar(x + width / 2, history["SpeciesClusters"] / 1e3, width, color=COLORS["Analysis"], label="Species clusters")
    axis.set_xticks(x, history["Release"])
    axis.set_ylabel("Count (thousands)")
    axis.grid(axis="y", color=COLORS["Grid"], linewidth=0.55)
    axis.legend(frameon=False)
    for xx, row in zip(x, history.itertuples(), strict=True):
        axis.text(xx - width / 2, row.Genomes / 1e3 + 18, f"{row.Genomes:,}", ha="center", rotation=90, fontsize=7.0)
        axis.text(xx + width / 2, row.SpeciesClusters / 1e3 + 18, f"{row.SpeciesClusters:,}", ha="center", rotation=90, fontsize=7.0)
    axis.set_ylim(0, 1040)
    axis.set_title("GTDB labels are release-bound", loc="left", fontweight="bold")
    axis.text(0.50, 0.04, "R232 released 15 Apr 2026", transform=axis.transAxes, ha="center", color=COLORS["Risk"], fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "73-catalogue-and-gtdb-releases")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output = args.figure_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    np.random.seed(PLOT_SEED)
    style()
    shutil.copy2(input_dir / "gem-figure1-original.png", output / "73-gem-figure1-original.png")
    resource_layers(input_dir, output)
    gem_biome_balance(input_dir, output)
    gem_world_map(input_dir, output)
    gem_quality(input_dir, output)
    tara_sampling(input_dir, output)
    metadata_completeness(input_dir, output)
    accession_crosswalk(input_dir, output)
    catalogue_landscape(input_dir, output)


if __name__ == "__main__":
    main()
