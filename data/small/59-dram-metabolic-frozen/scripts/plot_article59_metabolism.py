#!/usr/bin/env python3
"""Create publication-ready, English-only Article 59 figures."""

from __future__ import annotations

import argparse
import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/article59-matplotlib")
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
SKY = "#56B4E9"
YELLOW = "#E69F00"
PURPLE = "#CC79A7"
GRAY = "#7A7A7A"
DARK = "#253238"
LIGHT = "#EEF3F5"
WHITE = "#FFFFFF"


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
    if len(words) >= 2:
        return f"{words[0][0]}. {' '.join(words[1:])}"
    return name


def wrap_label(value: str, width: int = 24) -> str:
    return "\n".join(textwrap.wrap(value, width=width, break_long_words=False))


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "axes.titleweight": "bold",
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


def save_pub(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
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


def plot_pathway_matrix(summary: Path, figures: Path) -> None:
    frame = read_tsv(summary / "pathway-evidence-matrix.tsv")
    frame = numeric(frame, ["DRAMCoverage"])
    frame["AgreementBool"] = truth(frame["Agreement"])
    genome_order = (
        frame[["Genome", "Phylum", "Species"]]
        .drop_duplicates()
        .sort_values(["Phylum", "Species", "Genome"])["Genome"]
        .tolist()
    )
    module_meta = frame[["ModuleID", "DRAMModuleName"]].drop_duplicates()
    module_order = module_meta["ModuleID"].tolist()
    coverage = frame.pivot(index="Genome", columns="ModuleID", values="DRAMCoverage").loc[genome_order, module_order]
    agreement = frame.pivot(index="Genome", columns="ModuleID", values="AgreementBool").loc[genome_order, module_order]
    species = frame[["Genome", "Species"]].drop_duplicates().set_index("Genome")["Species"]
    module_name = module_meta.set_index("ModuleID")["DRAMModuleName"]
    xlabels = [f"{module}\n{wrap_label(module_name[module], 18)}" for module in module_order]
    ylabels = [f"{genome} · {species_short(species[genome])}" for genome in genome_order]

    fig, ax = plt.subplots(figsize=(16.5, 9.0), constrained_layout=True)
    sns.heatmap(
        coverage,
        ax=ax,
        cmap="viridis",
        vmin=0,
        vmax=1,
        linewidths=0.35,
        linecolor=WHITE,
        cbar_kws={"label": "DRAM module-step coverage", "shrink": 0.75},
    )
    discordant = np.argwhere(~agreement.to_numpy(dtype=bool))
    for row, column in discordant:
        ax.scatter(column + 0.5, row + 0.5, marker="x", s=23, linewidth=1.0, color="black")
    ax.set_xticklabels(xlabels, rotation=55, ha="right")
    ax.set_yticklabels(ylabels, rotation=0)
    ax.set_xlabel("Shared KEGG module (DRAM coverage; × marks METABOLIC disagreement at 0.75)")
    ax.set_ylabel("Quality-audited representative MAG")
    ax.set_title("Module coverage is continuous; pathway presence depends on an explicit cutoff")
    save_pub(fig, figures, "59-pathway-module-heatmap")


def plot_key_processes(summary: Path, figures: Path) -> None:
    frame = read_tsv(summary / "key-process-evidence.tsv")
    frame = frame.loc[frame["AnalysisSet"].eq("Primary real MAG")].copy()
    frame["PresentBool"] = truth(frame["Present"])
    prevalence = frame.groupby(["Category", "Function"], sort=False)["PresentBool"].sum().reset_index(name="PositiveMAGs")
    selected = prevalence.loc[prevalence["PositiveMAGs"].gt(0)].sort_values(
        ["Category", "PositiveMAGs", "Function"], ascending=[True, False, True]
    )
    if selected.empty:
        raise RuntimeError("No curated key-process rule was positive")
    functions = selected["Function"].tolist()
    genome_order = (
        frame[["Genome", "Phylum", "Species"]]
        .drop_duplicates()
        .sort_values(["Phylum", "Species", "Genome"])["Genome"]
        .tolist()
    )
    matrix = (
        frame.loc[frame["Function"].isin(functions)]
        .pivot(index="Genome", columns="Function", values="PresentBool")
        .loc[genome_order, functions]
        .astype(int)
    )
    species = frame[["Genome", "Species"]].drop_duplicates().set_index("Genome")["Species"]
    cmap = LinearSegmentedColormap.from_list("binary_pub", ["#F2F4F5", BLUE])
    fig, (ax_top, ax) = plt.subplots(
        2,
        1,
        figsize=(max(12.5, 0.7 * len(functions)), 9.5),
        gridspec_kw={"height_ratios": [1.0, 6.0], "hspace": 0.05},
    )
    counts = matrix.sum(axis=0)
    ax_top.bar(np.arange(len(functions)) + 0.5, counts, width=0.78, color=GREEN)
    ax_top.set_xlim(0, len(functions))
    ax_top.set_ylim(0, max(counts.max() * 1.18, 1))
    ax_top.set_ylabel("Positive\nMAGs")
    ax_top.set_xticks([])
    ax_top.spines[["top", "right", "bottom"]].set_visible(False)
    ax_top.grid(axis="x", visible=False)
    sns.heatmap(
        matrix,
        ax=ax,
        cmap=cmap,
        vmin=0,
        vmax=1,
        cbar=False,
        linewidths=0.35,
        linecolor=WHITE,
    )
    ax.set_xticklabels([wrap_label(name, 18) for name in functions], rotation=55, ha="right")
    ax.set_yticklabels([f"{genome} · {species_short(species[genome])}" for genome in genome_order], rotation=0)
    ax.set_xlabel("METABOLIC curated multi-gene or motif rule")
    ax.set_ylabel("Quality-audited representative MAG")
    fig.suptitle("Diagnostic metabolic traits are evidence rules, not measurements of activity", fontsize=12, fontweight="bold")
    fig.subplots_adjust(left=0.21, right=0.99, top=0.93, bottom=0.30)
    save_pub(fig, figures, "59-key-process-evidence")


def plot_tool_concordance(summary: Path, figures: Path) -> None:
    counts = read_tsv(summary / "ko-counts-by-genome.tsv")
    counts = counts.loc[counts["AnalysisSet"].eq("Primary real MAG")]
    counts = numeric(counts, ["DRAMComparableKOCount", "METABOLICKOCount", "Jaccard", "GenomeBp"])
    modules = numeric(read_tsv(summary / "module-agreement-summary.tsv"), ["AgreementRate", "Agreements", "PrimaryGenomes"])
    modules = modules.sort_values("AgreementRate")
    module_labels = modules["ModuleID"].tolist()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6), gridspec_kw={"width_ratios": [1.15, 1.0]})
    marker_map = {"Bacteria": "o", "Archaea": "^"}
    scatter = None
    for domain, group in counts.groupby("Domain"):
        scatter = ax1.scatter(
            group["DRAMComparableKOCount"],
            group["METABOLICKOCount"],
            c=group["Jaccard"],
            cmap="cividis",
            vmin=0,
            vmax=1,
            s=np.clip(group["GenomeBp"] / 35_000, 35, 180),
            marker=marker_map.get(domain, "s"),
            edgecolor=WHITE,
            linewidth=0.7,
            label=domain,
        )
    low = min(counts["DRAMComparableKOCount"].min(), counts["METABOLICKOCount"].min()) * 0.96
    high = max(counts["DRAMComparableKOCount"].max(), counts["METABOLICKOCount"].max()) * 1.04
    ax1.plot([low, high], [low, high], linestyle="--", color=GRAY, linewidth=1)
    counts = counts.assign(Difference=(counts["DRAMComparableKOCount"] - counts["METABOLICKOCount"]).abs())
    for _, row in counts.nlargest(4, "Difference").iterrows():
        ax1.annotate(row["Genome"], (row["DRAMComparableKOCount"], row["METABOLICKOCount"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax1.set_xlim(low, high)
    ax1.set_ylim(low, high)
    ax1.set_xlabel("Distinct KOs detected by DRAM (shared 2,678-KO grid)")
    ax1.set_ylabel("Distinct KOs detected by METABOLIC")
    ax1.set_title("A. Workflow-level KO concordance")
    ax1.legend(title="Domain", loc="lower right")
    if scatter is None:
        raise RuntimeError("No KO concordance points")
    colorbar = fig.colorbar(scatter, ax=ax1, fraction=0.046, pad=0.03)
    colorbar.set_label("KO-set Jaccard")

    colors = [GREEN if value >= 0.9 else YELLOW if value >= 0.75 else ORANGE for value in modules["AgreementRate"]]
    ax2.barh(module_labels, 100 * modules["AgreementRate"], color=colors)
    for index, row in enumerate(modules.itertuples(index=False)):
        ax2.text(min(100 * row.AgreementRate + 1.0, 97), index, f"{int(row.Agreements)}/{int(row.PrimaryGenomes)}", va="center", fontsize=7)
    ax2.axvline(75, color=GRAY, linewidth=0.8, linestyle=":")
    ax2.set_xlim(0, 105)
    ax2.set_xlabel("Agreement at the 0.75 presence cutoff (%)")
    ax2.set_ylabel("Shared KEGG module")
    ax2.set_title("B. Module calls also include tool-specific rules")
    sns.despine(fig=fig)
    fig.tight_layout()
    save_pub(fig, figures, "59-tool-concordance")


def plot_absence_audit(summary: Path, figures: Path) -> None:
    frame = numeric(
        read_tsv(summary / "completeness-absence-audit.tsv"),
        ["CompletenessPct", "ContaminationPct", "DRAMComparableKOCount", "METABOLICKOCount", "KOJaccard"],
    )
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.4))
    for row in frame.itertuples(index=False):
        ax1.plot(
            [row.CompletenessPct, row.CompletenessPct],
            [row.DRAMComparableKOCount, row.METABOLICKOCount],
            color="#B0BEC5",
            linewidth=0.8,
            zorder=1,
        )
    ax1.scatter(frame["CompletenessPct"], frame["DRAMComparableKOCount"], color=BLUE, s=42, label="DRAM", zorder=2)
    ax1.scatter(frame["CompletenessPct"], frame["METABOLICKOCount"], color=ORANGE, marker="s", s=38, label="METABOLIC", zorder=2)
    for _, row in frame.loc[(frame["CompletenessPct"] < 70) | (frame["ContaminationPct"] >= 5)].iterrows():
        right_edge = row["CompletenessPct"] > 95
        ax1.annotate(
            row["Genome"],
            (row["CompletenessPct"], max(row["DRAMComparableKOCount"], row["METABOLICKOCount"])),
            xytext=(-4 if right_edge else 3, 4),
            textcoords="offset points",
            ha="right" if right_edge else "left",
            fontsize=7,
        )
    ax1.axvspan(0, 70, color=ORANGE, alpha=0.06)
    ax1.axvline(90, color=GRAY, linestyle=":", linewidth=0.9)
    ax1.set_xlabel("CheckM2 completeness (%)")
    ax1.set_ylabel("Distinct KOs detected in the shared grid")
    ax1.set_title("A. Low completeness reduces observable gene content")
    ax1.legend()

    points = ax2.scatter(
        frame["CompletenessPct"],
        frame["KOJaccard"],
        c=frame["ContaminationPct"],
        cmap="magma_r",
        s=65,
        edgecolor=WHITE,
        linewidth=0.7,
    )
    for _, row in frame.loc[(frame["CompletenessPct"] < 70) | (frame["ContaminationPct"] >= 5)].iterrows():
        right_edge = row["CompletenessPct"] > 95
        ax2.annotate(
            row["Genome"],
            (row["CompletenessPct"], row["KOJaccard"]),
            xytext=(-4 if right_edge else 3, 3),
            textcoords="offset points",
            ha="right" if right_edge else "left",
            fontsize=7,
        )
    ax2.axvspan(0, 70, color=ORANGE, alpha=0.06)
    ax2.axvline(90, color=GRAY, linestyle=":", linewidth=0.9)
    ax2.set_ylim(0, 1.03)
    ax2.set_xlabel("CheckM2 completeness (%)")
    ax2.set_ylabel("DRAM–METABOLIC KO-set Jaccard")
    ax2.set_title("B. Contamination can create composite pathway evidence")
    colorbar = fig.colorbar(points, ax=ax2, fraction=0.046, pad=0.03)
    colorbar.set_label("CheckM2 contamination (%)")
    sns.despine(fig=fig)
    fig.tight_layout()
    save_pub(fig, figures, "59-completeness-absence-audit")


def plot_truncation(summary: Path, figures: Path) -> None:
    frame = numeric(
        read_tsv(summary / "truncation-sensitivity.tsv"),
        [
            "RetentionObservedPct", "DRAMFullKORetentionPct", "METABOLICKORetentionPct",
            "DRAMParentModules", "DRAMRetainedParentModules", "METABOLICParentModules",
            "METABOLICRetainedParentModules",
        ],
    ).sort_values("RetentionObservedPct")
    frame["DRAMModuleRetentionPct"] = np.where(
        frame["DRAMParentModules"] > 0,
        100 * frame["DRAMRetainedParentModules"] / frame["DRAMParentModules"],
        np.nan,
    )
    frame["METABOLICModuleRetentionPct"] = np.where(
        frame["METABOLICParentModules"] > 0,
        100 * frame["METABOLICRetainedParentModules"] / frame["METABOLICParentModules"],
        np.nan,
    )
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.8))
    for column, label, color, marker in (
        ("DRAMFullKORetentionPct", "DRAM full KOfam", BLUE, "o"),
        ("METABOLICKORetentionPct", "METABOLIC", ORANGE, "s"),
    ):
        ax1.plot(frame["RetentionObservedPct"], frame[column], color=color, marker=marker, linewidth=1.8, label=label)
    ax1.plot([45, 102], [45, 102], color=GRAY, linestyle="--", linewidth=0.9, label="Sequence retained")
    ax1.set_xlim(45, 103)
    ax1.set_ylim(40, 103)
    ax1.set_xlabel("Parent genome sequence retained (%)")
    ax1.set_ylabel("Parent KO set retained (%)")
    ax1.set_title("A. Missing sequence creates false metabolic absences")
    ax1.legend()

    for column, label, color, marker in (
        ("DRAMModuleRetentionPct", "DRAM", BLUE, "o"),
        ("METABOLICModuleRetentionPct", "METABOLIC", ORANGE, "s"),
    ):
        if frame[column].notna().any():
            ax2.plot(frame["RetentionObservedPct"], frame[column], color=color, marker=marker, linewidth=1.8, label=label)
    ax2.set_xlim(45, 103)
    ax2.set_ylim(-3, 103)
    ax2.set_xlabel("Parent genome sequence retained (%)")
    ax2.set_ylabel("Parent positive modules retained (%)")
    ax2.set_title("B. Binary pathway calls hide continuous information loss")
    ax2.legend()
    sns.despine(fig=fig)
    fig.tight_layout()
    save_pub(fig, figures, "59-truncation-sensitivity")


def plot_evidence_ladder(summary: Path, figures: Path) -> None:
    frame = numeric(read_tsv(summary / "evidence-ladder.tsv"), ["Rank"]).sort_values("Rank")
    colors = ["#DDEBF3", "#B9DCEB", "#8BC7D9", "#4FAEA5", "#1B7F68"]
    fig, ax = plt.subplots(figsize=(13.2, 5.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0.25, 5.75)
    ax.axis("off")
    for index, row in enumerate(frame.itertuples(index=False), start=1):
        y = 6 - index
        width = 4.6 + 0.85 * index
        x = (10 - width) / 2
        box = FancyBboxPatch(
            (x, y - 0.37),
            width,
            0.74,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            linewidth=0.8,
            edgecolor=WHITE,
            facecolor=colors[index - 1],
        )
        ax.add_patch(box)
        ax.text(x + 0.18, y, f"{int(row.Rank)}", va="center", ha="left", fontsize=12, fontweight="bold", color=DARK)
        ax.text(x + 0.62, y + 0.12, row.Evidence, va="center", ha="left", fontsize=9, fontweight="bold", color=DARK)
        ax.text(x + 0.62, y - 0.15, f"Claim ceiling: {row.ClaimCeiling}", va="center", ha="left", fontsize=8, color=DARK)
    ax.text(5, 5.62, "Evidence must rise before the biological claim rises", ha="center", fontsize=13, fontweight="bold", color=DARK)
    ax.annotate("", xy=(9.55, 0.75), xytext=(9.55, 5.1), arrowprops={"arrowstyle": "-|>", "color": GRAY, "lw": 1.2})
    ax.text(9.72, 2.9, "Increasing biological support", rotation=90, va="center", fontsize=8, color=GRAY)
    save_pub(fig, figures, "59-metabolism-evidence-ladder")


def main() -> None:
    args = parse_args()
    summary = args.summary_dir.resolve()
    figures = args.figure_dir.resolve()
    required = (
        "pathway-evidence-matrix.tsv",
        "key-process-evidence.tsv",
        "ko-counts-by-genome.tsv",
        "module-agreement-summary.tsv",
        "completeness-absence-audit.tsv",
        "truncation-sensitivity.tsv",
        "evidence-ladder.tsv",
    )
    missing = [name for name in required if not (summary / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing Article 59 summary tables: {missing}")
    setup_style()
    plot_pathway_matrix(summary, figures)
    plot_key_processes(summary, figures)
    plot_tool_concordance(summary, figures)
    plot_absence_audit(summary, figures)
    plot_truncation(summary, figures)
    plot_evidence_ladder(summary, figures)
    print(f"Created 6 Article 59 figure families in {figures}")


if __name__ == "__main__":
    main()
