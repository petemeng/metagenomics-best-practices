#!/usr/bin/env python3
"""Render Article 76 publication figures from the frozen reporting ledgers."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch, Patch


PLOT_SEED = 20_260_776
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PINK = "#CC79A7"
GOLD = "#E69F00"
SKY = "#56B4E9"
INK = "#263238"
GREY = "#90A4AE"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 14,
            "axes.labelsize": 10,
            "axes.titleweight": "bold",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
        }
    )


def save(fig: plt.Figure, output: Path, stem: str) -> list[Path]:
    output.mkdir(parents=True, exist_ok=True)
    paths = [output / f"{stem}.png", output / f"{stem}.svg"]
    fig.savefig(paths[0], dpi=360, bbox_inches="tight", pad_inches=0.14)
    fig.savefig(paths[1], bbox_inches="tight", pad_inches=0.14)
    plt.close(fig)
    return paths


def selection_plot(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "standard-selection-matrix.tsv", sep="\t")
    scenarios = frame["Scenario"].drop_duplicates().tolist()
    standards = ["STORMS", "STREAMS", "MIMAG", "MIUViG"]
    matrix = frame.pivot(index="Scenario", columns="Standard", values="RoleCode").loc[scenarios, standards]
    fig, ax = plt.subplots(figsize=(9.4, 6.4))
    ax.imshow(matrix, aspect="auto", vmin=0, vmax=2, cmap=ListedColormap(["#EDF1F2", "#E9C46A", "#2A9D8F"]))
    ax.set_xticks(np.arange(4), standards, fontsize=11, weight="bold")
    ax.set_yticks(np.arange(len(scenarios)), scenarios, fontsize=9)
    labels = {0: "—", 1: "ADD", 2: "CORE"}
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix.iloc[i, j])
            ax.text(j, i, labels[value], ha="center", va="center", color="white" if value == 2 else INK, weight="bold")
    ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(scenarios), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("Choose one study checklist, then add entity standards", loc="left", pad=12)
    fig.text(0.01, 0.01, "STORMS: human studies · STREAMS: environmental/non-human/synthetic · MIMAG/MIUViG: entity add-ons", fontsize=8.5, color="#51636B")
    return save(fig, output, "76-standard-selection")


def layer_stack(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "reporting-layer-map.tsv", sep="\t")
    layers = frame["Layer"].drop_duplicates().tolist()
    colors = [BLUE, SKY, GOLD, GREEN]
    fig, ax = plt.subplots(figsize=(11.2, 6.8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7.1)
    for index, layer in enumerate(layers):
        y = 5.6 - index * 1.35
        color = colors[index]
        ax.add_patch(FancyBboxPatch((0.35, y - 0.45), 2.5, 0.9, boxstyle="round,pad=0.02,rounding_size=0.05", facecolor=color, edgecolor="none"))
        ax.text(1.6, y, layer.replace(" and ", "\nand "), ha="center", va="center", color="white", weight="bold")
        subset = frame[frame["Layer"].eq(layer)]
        for j, row in enumerate(subset.itertuples(index=False)):
            x = 3.2 + j * 2.1
            active = "Core" in row.Role
            face = color if active else "#E8EEF0"
            text_color = "white" if active else INK
            ax.add_patch(FancyBboxPatch((x, y - 0.36), 1.75, 0.72, boxstyle="round,pad=0.015,rounding_size=0.04", facecolor=face, edgecolor="white", linewidth=1.5))
            ax.text(x + 0.875, y + 0.10, row.Standard, ha="center", va="center", color=text_color, weight="bold", fontsize=9.5)
            ax.text(x + 0.875, y - 0.16, row.Role, ha="center", va="center", color=text_color, fontsize=7.4)
        if index < len(layers) - 1:
            ax.annotate("", xy=(1.6, y - 0.83), xytext=(1.6, y - 0.52), arrowprops={"arrowstyle": "->", "lw": 1.8, "color": GREY})
    ax.text(3.2, 6.55, "Study-level standards", color=BLUE, weight="bold")
    ax.text(7.4, 6.55, "Entity-level add-ons", color=GREEN, weight="bold")
    ax.axis("off")
    ax.set_title("Reporting is a linked record stack, not a single checklist score", loc="left", pad=10)
    return save(fig, output, "76-reporting-layer-stack")


def section_counts(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "checklist-section-counts.tsv", sep="\t")
    sections = ["Abstract", "Introduction", "Methods", "Results", "Discussion", "Other information"]
    storms = frame[frame["Standard"].eq("STORMS")].set_index("ManuscriptSection").loc[sections, "ExpandedRecommendations"]
    streams = frame[frame["Standard"].eq("STREAMS")].set_index("ManuscriptSection").loc[sections, "ExpandedRecommendations"]
    y = np.arange(len(sections))
    fig, ax = plt.subplots(figsize=(9.4, 6.1))
    ax.barh(y + 0.18, storms, height=0.34, color=BLUE, label="STORMS v1.03 · 69 expanded rows")
    ax.barh(y - 0.18, streams, height=0.34, color=GREEN, label="STREAMS v1.0 · 67 recommendations")
    for values, offset in ((storms, 0.18), (streams, -0.18)):
        for i, value in enumerate(values):
            ax.text(value + 0.7, i + offset, str(int(value)), va="center", fontsize=9, color=INK)
    ax.set_yticks(y, sections)
    ax.invert_yaxis()
    ax.set_xlabel("Expanded recommendations")
    ax.set_xlim(0, max(max(storms), max(streams)) + 8)
    ax.grid(axis="x", color="#EDF1F2")
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("Most reporting obligations accumulate in Methods", loc="left")
    fig.text(0.01, 0.01, "STORMS has 17 top-level items; its 69 expanded rows include numbered subitems.", fontsize=8.5, color="#51636B")
    return save(fig, output, "76-checklist-sections")


def crosswalk_plot(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "standards-crosswalk.tsv", sep="\t")
    domains = frame["ReportingDomain"].drop_duplicates().tolist()
    standards = ["STORMS", "STREAMS", "MIMAG", "MIUViG"]
    matrix = frame.pivot(index="ReportingDomain", columns="Standard", values="CoverageCode").loc[domains, standards]
    fig, ax = plt.subplots(figsize=(9.6, 7.5))
    ax.imshow(matrix, aspect="auto", vmin=0, vmax=2, cmap=ListedColormap(["#ECEFF1", "#E9C46A", "#2A9D8F"]))
    ax.set_xticks(np.arange(4), standards, fontsize=10.5, weight="bold")
    ax.set_yticks(np.arange(len(domains)), [textwrap.fill(x, 34) for x in domains], fontsize=8.6)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix.iloc[i, j])
            coverage = frame[(frame["ReportingDomain"].eq(domains[i])) & (frame["Standard"].eq(standards[j]))]["Coverage"].iloc[0]
            ax.text(j, i, coverage.replace("Not primary", "not primary").replace("Context-specific", "context"), ha="center", va="center", fontsize=7.5, color="white" if value == 2 else INK, weight="bold" if value == 2 else "normal")
    ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(domains), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("The four standards cover different reporting domains", loc="left", pad=12)
    fig.text(0.01, 0.01, "Interpretive crosswalk based on the version-locked source tables; it is not an official equivalence map.", fontsize=8.3, color="#51636B")
    return save(fig, output, "76-standards-crosswalk")


def mimag_plot(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "article44-mimag-compliance.tsv", sep="\t")
    colors = frame["Article44ExtendedTier"].map({"High quality": GREEN, "Medium quality": BLUE, "Low/failed": ORANGE})
    marker = frame["CompleteRRNASet"].astype(str).str.lower().map({"true": "o", "false": "X"})
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    for shape, label in (("o", "Complete 5S/16S/23S"), ("X", "Incomplete rRNA set")):
        subset = marker.eq(shape)
        ax.scatter(frame.loc[subset, "CheckM2Completeness"], frame.loc[subset, "CheckM2Contamination"], s=55 + frame.loc[subset, "TRNAIsotypes"] * 3, c=colors[subset], marker=shape, edgecolor="white", linewidth=0.8, alpha=0.9, label=label)
    ax.axvline(90, color="#455A64", linestyle="--", linewidth=1.2)
    ax.axhline(5, color="#455A64", linestyle="--", linewidth=1.2)
    ax.axvline(50, color=GREY, linestyle=":", linewidth=1.2)
    ax.axhline(10, color=GREY, linestyle=":", linewidth=1.2)
    ax.text(90.4, 8.8, "HQ completeness gate", fontsize=8, color="#455A64", rotation=90, va="top")
    ax.text(51, 9.35, "MQ boundary", fontsize=8, color=GREY)
    ax.set_xlabel("CheckM2 completeness (%)")
    ax.set_ylabel("CheckM2 contamination (%)")
    ax.set_xlim(48, 101.5)
    ax.set_ylim(-0.35, 10.7)
    ax.grid(color="#EDF1F2")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Completeness and contamination alone do not confer MIMAG high quality", loc="left")
    fig.text(0.01, 0.01, "n = 23 MAGs · point size reflects tRNA isotype count · color is the Article 44 MIMAG + GUNC tier", fontsize=8.3, color="#51636B")
    return save(fig, output, "76-mimag-compliance")


def miuvig_plot(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "article54-miuvig-compliance.tsv", sep="\t")
    colors = frame["Status"].map({"Complete": GREEN, "Partial": GOLD, "Missing": ORANGE})
    y = np.arange(len(frame))
    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    ax.barh(y, np.ones(len(frame)), color=colors, height=0.62)
    ax.set_yticks(y, [textwrap.fill(x, 34) for x in frame["MandatoryMetadata"]], fontsize=9)
    ax.invert_yaxis()
    for i, status in enumerate(frame["Status"]):
        ax.text(0.5, i, status.upper(), ha="center", va="center", color="white" if status != "Partial" else INK, weight="bold", fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_title("A viral quality label is only one of eight MIUViG records", loc="left")
    ax.legend(handles=[Patch(facecolor=GREEN, label="Complete"), Patch(facecolor=GOLD, label="Partial"), Patch(facecolor=ORANGE, label="Missing")], frameon=False, ncol=3, loc="lower right")
    fig.text(0.01, 0.01, "Audit of the Article 54 public CheckV fixture: 4 complete · 2 partial · 2 missing fields.", fontsize=8.3, color="#51636B")
    return save(fig, output, "76-miuvig-readiness")


def owner_timeline(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "field-responsibility-matrix.tsv", sep="\t")
    owners = frame["Owner"].drop_duplicates().tolist()
    palette = [
        BLUE,
        GREEN,
        ORANGE,
        PINK,
        GOLD,
        SKY,
        "#264653",
        "#8D99AE",
        "#6A4C93",
        "#457B9D",
    ]
    colors = dict(zip(owners, palette[: len(owners)], strict=True))
    fig, ax = plt.subplots(figsize=(11.3, 6.7))
    ax.plot(frame["Order"], np.zeros(len(frame)), color="#CFD8DC", linewidth=5, zorder=1)
    for row in frame.itertuples(index=False):
        up = row.Order % 2 == 1
        y = 1.05 if up else -1.05
        ax.scatter(row.Order, 0, s=330, color=colors[row.Owner], edgecolor="white", linewidth=1.4, zorder=3)
        ax.text(row.Order, 0, str(row.Order), ha="center", va="center", color="white", weight="bold", fontsize=8.5, zorder=4)
        ax.plot([row.Order, row.Order], [0.2 if up else -0.2, y * 0.72], color=colors[row.Owner], linewidth=1.2)
        ax.text(row.Order, y, textwrap.fill(row.Milestone, 16), ha="center", va="center", weight="bold", color=INK, fontsize=8)
        ax.text(row.Order, y - (0.31 if up else -0.31), textwrap.fill(row.Owner, 18), ha="center", va="center", color=colors[row.Owner], fontsize=7.2)
    ax.set_xlim(0.3, 12.7)
    ax.set_ylim(-1.75, 1.75)
    ax.axis("off")
    ax.set_title("Reporting fields need owners before the manuscript is written", loc="left", pad=8)
    fig.text(0.01, 0.01, "Collection-time metadata cannot be reconstructed reliably at submission; freeze each hand-off artifact at its due stage.", fontsize=8.4, color="#51636B")
    return save(fig, output, "76-owner-timeline")


def na_plot(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "not-applicable-ledger.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(10.7, 6.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.7, len(frame) - 0.2)
    for i, row in enumerate(frame.itertuples(index=False)):
        y = len(frame) - i - 1
        color = GREEN if row.Disposition == "Report" else SKY
        ax.add_patch(FancyBboxPatch((0.15, y - 0.34), 9.5, 0.68, boxstyle="round,pad=0.01,rounding_size=0.025", facecolor=color, edgecolor="none", alpha=0.10))
        ax.text(0.35, y + 0.08, textwrap.fill(row.Field, 28), va="center", weight="bold", color=INK, fontsize=8.3)
        ax.text(3.15, y + 0.08, row.StandardItem, va="center", color=color, weight="bold", fontsize=8)
        ax.text(4.75, y + 0.08, row.Disposition.upper(), va="center", color=color, weight="bold", fontsize=8)
        ax.text(6.35, y + 0.08, row.Approver, va="center", color=INK, fontsize=8)
        ax.text(0.35, y - 0.19, textwrap.fill(row.Reason, 110), va="center", color="#52636B", fontsize=7.1)
    ax.text(0.35, len(frame) - 0.15, "Field", weight="bold")
    ax.text(3.15, len(frame) - 0.15, "Standard", weight="bold")
    ax.text(4.75, len(frame) - 0.15, "Decision", weight="bold")
    ax.text(6.35, len(frame) - 0.15, "Approver", weight="bold")
    ax.axis("off")
    ax.set_title("Not applicable is an adjudicated record, not an empty cell", loc="left", pad=12)
    return save(fig, output, "76-na-ledger")


def readiness_plot(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "submission-readiness.tsv", sep="\t")
    order = ["Complete", "Partial", "Pending", "Missing", "Not selected"]
    colors = {"Complete": GREEN, "Partial": GOLD, "Pending": PINK, "Missing": ORANGE, "Not selected": GREY}
    groups = ["STREAMS", "MIMAG", "MIUViG", "All four", "Other"]
    frame["Group"] = np.select(
        [
            frame["Standard"].str.contains("STREAMS", regex=False),
            frame["Standard"].str.contains("MIMAG", regex=False),
            frame["Standard"].str.contains("MIUViG", regex=False),
            frame["Standard"].eq("All four"),
        ],
        ["STREAMS", "MIMAG", "MIUViG", "All four"],
        default="Other",
    )
    counts = frame.groupby(["Group", "Status"]).size().unstack(fill_value=0).reindex(index=groups, columns=order, fill_value=0)
    fig, ax = plt.subplots(figsize=(9.5, 5.9))
    left = np.zeros(len(counts))
    for status in order:
        values = counts[status].to_numpy()
        bars = ax.barh(counts.index, values, left=left, color=colors[status], label=status)
        for bar, value, base in zip(bars, values, left, strict=True):
            if value:
                ax.text(base + value / 2, bar.get_y() + bar.get_height() / 2, str(int(value)), ha="center", va="center", color="white" if status != "Partial" else INK, weight="bold")
        left += values
    ax.set_xlabel("Readiness records")
    ax.grid(axis="x", color="#EDF1F2")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    ax.set_title("The worked evidence bundle is not yet submission-ready", loc="left")
    fig.text(0.01, 0.01, "Missing or partial records trigger a named next action; checklist completion is never converted to a quality score.", fontsize=8.3, color="#51636B")
    return save(fig, output, "76-submission-readiness")


def main() -> None:
    args = parse_args()
    data = args.input_dir.resolve()
    output = args.figures_dir.resolve()
    np.random.seed(PLOT_SEED)
    style()
    generated: list[Path] = []
    generated += selection_plot(data, output)
    generated += layer_stack(data, output)
    generated += section_counts(data, output)
    generated += crosswalk_plot(data, output)
    generated += mimag_plot(data, output)
    generated += miuvig_plot(data, output)
    generated += owner_timeline(data, output)
    generated += na_plot(data, output)
    generated += readiness_plot(data, output)

    anchor = output / "76-streams-figure1-original.png"
    shutil.copy2(data / "streams-figure1-original.png", anchor)
    generated.append(anchor)
    records = [
        {
            "file": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "kind": "source-anchor" if path == anchor else path.suffix.lstrip("."),
        }
        for path in sorted(generated)
    ]
    manifest = {"article": 76, "plot_seed": PLOT_SEED, "figures": records}
    (data / "figure-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"plotted\t{output}\t{len(generated)} files")


if __name__ == "__main__":
    main()
