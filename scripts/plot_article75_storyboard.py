#!/usr/bin/env python3
"""Render Article 75 publication figures from the prepared storyboard tables."""

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


PLOT_SEED = 20_260_775
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#264653", "#8D99AE"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def save(fig: plt.Figure, output: Path, name: str) -> list[Path]:
    output.mkdir(parents=True, exist_ok=True)
    png = output / f"{name}.png"
    svg = output / f"{name}.svg"
    fig.savefig(png, dpi=360, bbox_inches="tight", pad_inches=0.14)
    fig.savefig(svg, bbox_inches="tight", pad_inches=0.14)
    plt.close(fig)
    return [png, svg]


def box(ax: plt.Axes, x: float, y: float, width: float, height: float, color: str, title: str, body: str = "") -> None:
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.05",
        facecolor=color, edgecolor="none",
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height * 0.64, title, ha="center", va="center", color="white", weight="bold", fontsize=10)
    if body:
        ax.text(x + width / 2, y + height * 0.28, body, ha="center", va="center", color="white", fontsize=8.2)


def paper_arc(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "wirbel-main-figure-ledger.tsv", sep="\t")
    short = ["Core markers", "CRC subgroups", "Model transfer", "Function + qPCR", "External cohorts"]
    fig, ax = plt.subplots(figsize=(11.0, 5.4))
    ax.set_xlim(0.4, 5.6)
    ax.set_ylim(-0.35, 3.4)
    ax.plot(frame["Figure"], np.repeat(1.55, len(frame)), color="#CBD5DA", lw=5, zorder=1)
    for idx, row in frame.iterrows():
        number = int(row["Figure"])
        color = PALETTE[idx]
        ax.scatter(number, 1.55, s=760, color=color, edgecolor="white", lw=2, zorder=3)
        ax.text(number, 1.55, str(number), ha="center", va="center", color="white", weight="bold", fontsize=14)
        ypos = 2.22 if number % 2 else 0.5
        ax.text(number, ypos + 0.35, short[idx], ha="center", va="center", weight="bold", fontsize=10)
        ax.text(number, ypos, f"{int(row['Panels'])} panels", ha="center", va="center", fontsize=8.5, color="#455A64")
        ax.plot([number, number], [1.85 if number % 2 else 1.25, ypos + (0.16 if number % 2 else 0.18)], color=color, lw=1.3)
    ax.text(0.48, 3.08, "Discovery", weight="bold", color=PALETTE[0], fontsize=10)
    ax.annotate("", xy=(5.48, 3.08), xytext=(1.08, 3.08), arrowprops={"arrowstyle": "->", "lw": 1.7, "color": "#264653"})
    ax.text(5.5, 3.08, "Validation", ha="right", weight="bold", color=PALETTE[3], fontsize=10)
    ax.axis("off")
    ax.set_title("A strong paper advances one evidence chain across its main figures", loc="left", pad=10)
    fig.text(0.01, 0.01, "Wirbel et al. 2019: discovery → heterogeneity → transfer → orthogonal assay → independent populations", fontsize=8.5, color="#51636B")
    return save(fig, output, "75-wirbel-figure-arc")


def five_figure_storyboard(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "main-figure-storyboard.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(11.2, 7.0))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7.1)
    ax.text(1.35, 6.42, "Narrative role", weight="bold", color="#263238")
    ax.text(4.75, 6.42, "Question", weight="bold", color="#263238")
    ax.text(8.25, 6.42, "Validation gate", weight="bold", color="#263238")
    for idx, row in enumerate(frame.itertuples(index=False)):
        y = 5.78 - idx * 1.08
        color = PALETTE[idx]
        ax.add_patch(FancyBboxPatch((0.32, y - 0.38), 11.2, 0.78, boxstyle="round,pad=0.01,rounding_size=0.03", facecolor=color, edgecolor="none", alpha=0.10))
        ax.scatter(0.78, y, s=520, color=color, edgecolor="white", linewidth=1.5, zorder=3)
        ax.text(0.78, y, str(row.Figure), ha="center", va="center", color="white", weight="bold", fontsize=12)
        ax.text(1.35, y + 0.13, textwrap.fill(row.ShortTitle, width=27), va="center", weight="bold", fontsize=8.5, color=color, linespacing=0.95)
        ax.text(1.35, y - 0.21, f"Unit: {row.PrimaryUnit}", va="center", fontsize=7.2, color="#607078")
        ax.text(4.75, y, textwrap.fill(row.Question, width=31), va="center", fontsize=8.3, color="#263238")
        ax.text(8.25, y, textwrap.fill(row.ValidationGate, width=38), va="center", fontsize=8.5, color="#263238", weight="bold")
        if idx < len(frame) - 1:
            ax.annotate("", xy=(0.78, y - 0.66), xytext=(0.78, y - 0.45), arrowprops={"arrowstyle": "->", "lw": 1.5, "color": "#90A4AE"})
    ax.axhspan(0.05, 0.55, color="#F2F5F6")
    ax.text(0.35, 0.30, "STOP RULE", va="center", weight="bold", color="#D55E00")
    ax.text(1.65, 0.30, "A failed validation gate narrows the claim and remains visible in the supplement.", va="center", color="#37474F")
    ax.axis("off")
    ax.set_title("Five figures are five claim gates—not five software outputs", loc="left", pad=10)
    return save(fig, output, "75-five-figure-storyboard")


def panel_budget(data: Path, output: Path) -> list[Path]:
    panels = pd.read_csv(data / "panel-register.tsv", sep="\t")
    roles = ["Context", "Diagnostic", "Primary", "Secondary", "Validation", "Sensitivity"]
    colors = dict(zip(roles, ["#B0BEC5", "#56B4E9", "#0072B2", "#E69F00", "#009E73", "#CC79A7"], strict=True))
    counts = panels.groupby(["Figure", "PanelRole"]).size().unstack(fill_value=0).reindex(columns=roles, fill_value=0)
    fig, ax = plt.subplots(figsize=(10.0, 5.5))
    bottom = np.zeros(len(counts))
    for role in roles:
        values = counts[role].to_numpy()
        bars = ax.bar(counts.index.astype(str), values, bottom=bottom, color=colors[role], width=0.65, label=role)
        for bar, value, base in zip(bars, values, bottom, strict=True):
            if value:
                ax.text(bar.get_x() + bar.get_width() / 2, base + value / 2, str(int(value)), ha="center", va="center", fontsize=9, color="white" if role in {"Primary", "Validation"} else "#263238", weight="bold")
        bottom += values
    ax.set_xlabel("Main figure")
    ax.set_ylabel("Planned panels")
    ax.set_ylim(0, max(bottom) + 1.4)
    ax.set_title("Panel budgets protect the primary result and its validation", loc="left")
    ax.legend(frameon=False, ncol=3, loc="upper center")
    ax.grid(axis="y", color="#EDF1F2")
    fig.text(0.01, 0.01, "Context and diagnostics establish denominators; validation and sensitivity defend the claim.", fontsize=8.5, color="#51636B")
    return save(fig, output, "75-panel-budget")


def main_supplement(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "main-supplement-map.tsv", sep="\t")
    matrix = frame[["MainSpace", "SupplementDetail"]].to_numpy()
    labels = [
        "Study flow", "Read QC + controls", "All feature effects", "Compositional sensitivity",
        "External performance", "Assembly diagnostics", "MAG acceptance", "Strain robustness",
        "Virus + vOTU audit", "Virus-host negatives", "Model + gap-fill", "Medium sensitivity",
        "Version ledger", "Null and failed analyses",
    ]
    fig, ax = plt.subplots(figsize=(8.8, 7.2))
    cmap = ListedColormap(["#F4F6F7", "#E9C46A", "#2A9D8F"])
    ax.imshow(matrix, aspect="auto", vmin=0, vmax=2, cmap=cmap)
    ax.set_xticks([0, 1], ["Main text", "Supplement"])
    ax.set_yticks(np.arange(len(labels)), labels)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix[i, j])
            ax.text(j, i, ["—", "summary", "full"][value], ha="center", va="center", color="white" if value == 2 else "#37474F", weight="bold", fontsize=8.5)
    ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("The supplement stores the denominator, not the leftovers", loc="left", pad=12)
    fig.text(0.01, 0.01, "full = complete audit payload; summary = one claim-facing panel; — = no decorative main-text panel", fontsize=8, color="#51636B")
    return save(fig, output, "75-main-supplement")


def claim_evidence(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "claim-evidence-matrix.tsv", sep="\t")
    claims = frame[["ClaimID", "Claim"]].drop_duplicates()
    layers = frame["EvidenceLayer"].drop_duplicates().tolist()
    matrix = frame.pivot(index="Claim", columns="EvidenceLayer", values="RequirementCode").loc[claims["Claim"], layers].to_numpy()
    cmap = ListedColormap(["#ECEFF1", "#D9E2E6", "#E9C46A", "#2A9D8F"])
    fig, ax = plt.subplots(figsize=(11.2, 5.2))
    ax.imshow(matrix + 1, aspect="auto", vmin=0, vmax=3, cmap=cmap)
    ax.set_xticks(np.arange(len(layers)), [text.replace(" and ", "\n") for text in layers], fontsize=8.5)
    ax.set_yticks(np.arange(len(claims)), claims["Claim"], fontsize=9)
    symbol = {-1: "N/A", 0: "alone ≠ enough", 1: "support", 2: "required"}
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix[i, j])
            ax.text(j, i, symbol[value], ha="center", va="center", fontsize=7.7, color="white" if value == 2 else "#37474F", weight="bold" if value == 2 else "normal")
    ax.set_xticks(np.arange(-0.5, len(layers), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(claims), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("Claim strength determines the minimum evidence package", loc="left", pad=12)
    return save(fig, output, "75-claim-evidence")


def traceability_audit(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "result-traceability-ledger.tsv", sep="\t")
    columns = ["UnitRecorded", "StatisticRecorded", "ValidationRecorded", "ChecksumRequired"]
    matrix = frame[columns].astype(bool).astype(int).to_numpy()
    fig, ax = plt.subplots(figsize=(8.6, 8.0))
    ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap=ListedColormap(["#E76F51", "#2A9D8F"]))
    ax.set_xticks(np.arange(4), ["Unit", "Statistic", "Validation gate", "Checksum"], fontsize=9)
    ax.set_yticks(np.arange(len(frame)), frame["ResultID"], fontsize=7.5)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, "✓" if matrix[i, j] else "!", ha="center", va="center", color="white", weight="bold", fontsize=8)
    for boundary in np.cumsum(frame.groupby("Figure").size())[:-1] - 0.5:
        ax.axhline(boundary, color="white", lw=3)
    ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
    ax.grid(which="minor", axis="x", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False)
    ax.set_title("Every plotted result needs a four-field traceability record", loc="left", pad=12)
    fig.text(0.01, 0.01, "A panel is not final until its independent unit, statistic, gate, and byte identity are recorded.", fontsize=8.3, color="#51636B")
    return save(fig, output, "75-traceability-audit")


def sensitivity_matrix(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "sensitivity-matrix.tsv", sep="\t")
    axes = frame["SensitivityAxis"].drop_duplicates().tolist()
    matrix = frame.pivot(index="SensitivityAxis", columns="Figure", values="RequirementCode").loc[axes, [1, 2, 3, 4, 5]].to_numpy()
    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    ax.imshow(matrix, aspect="auto", vmin=0, vmax=2, cmap=ListedColormap(["#ECEFF1", "#E9C46A", "#2A9D8F"]))
    ax.set_xticks(np.arange(5), [f"Figure {number}" for number in range(1, 6)])
    ax.set_yticks(np.arange(len(axes)), axes, fontsize=8.5)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix[i, j])
            ax.text(j, i, ["—", "conditional", "required"][value], ha="center", va="center", fontsize=7.5, color="white" if value == 2 else "#37474F", weight="bold" if value == 2 else "normal")
    ax.set_xticks(np.arange(-0.5, 5, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(axes), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("Sensitivity analyses follow the claim they can overturn", loc="left", pad=12)
    return save(fig, output, "75-sensitivity-matrix")


def style_contract(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "figure-style-contract.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(10.0, 6.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.8, len(frame) - 0.2)
    marker_map = {"circle": "o", "diamond": "D", "cross": "X", "none": None}
    for idx, row in enumerate(frame.itertuples(index=False)):
        y = len(frame) - idx - 1
        marker = marker_map[row.Marker]
        if marker is not None:
            ax.scatter(0.75, y, s=220, marker=marker, color=row.Hex, edgecolor="white", lw=1.3)
        else:
            ax.plot([0.3, 1.2], [y, y], color=row.Hex, lw=3, linestyle="--" if row.LineType == "dashed" else "-")
        ax.text(1.55, y, row.Semantic, va="center", weight="bold", fontsize=9.3)
        ax.text(4.55, y, row.Channel, va="center", color="#455A64", fontsize=8.7)
        ax.text(6.1, y, row.DoNotReuseFor, va="center", color="#455A64", fontsize=8.5)
    ax.text(1.55, len(frame) - 0.25, "Meaning", weight="bold", color="#263238")
    ax.text(4.55, len(frame) - 0.25, "Channel", weight="bold", color="#263238")
    ax.text(6.1, len(frame) - 0.25, "Reserved use", weight="bold", color="#263238")
    for y in np.arange(0.5, len(frame), 1):
        ax.axhline(y, color="#EDF1F2", lw=0.8)
    ax.axis("off")
    ax.set_title("Color encodes biology; shape and outline encode study role", loc="left", pad=12)
    fig.text(0.01, 0.01, "Keep the same semantic mapping across every panel, supplement, graphical abstract, and slide.", fontsize=8.3, color="#51636B")
    return save(fig, output, "75-style-contract")


def reviewer_attack_map(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "reviewer-attack-map.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(9.5, 6.4))
    colors = [PALETTE[1] if impact == 5 else PALETTE[4] for impact in frame["Impact"]]
    jitter = {
        "Outcome leakage": (0.08, 0.08),
        "Cohort confounding": (-0.02, -0.02),
        "Missing external validation": (-0.12, -0.12),
        "Compositional artifact": (-0.10, 0.05),
        "Incomplete MAG QC": (-0.10, 0.08),
        "Strain threshold drift": (0.00, 0.00),
        "Virus false positive": (0.00, -0.02),
        "Unsupported host link": (0.10, -0.12),
        "Gap-fill dependence": (0.05, -0.05),
        "Version drift": (0.15, -0.15),
    }
    label_positions = {
        "Outcome leakage": (4.62, 5.23),
        "Cohort confounding": (4.62, 5.02),
        "Missing external validation": (4.40, 4.72),
        "Compositional artifact": (3.12, 4.19),
        "Incomplete MAG QC": (3.52, 5.23),
        "Strain threshold drift": (3.10, 3.76),
        "Virus false positive": (3.52, 5.00),
        "Unsupported host link": (3.50, 4.73),
        "Gap-fill dependence": (4.18, 4.18),
        "Version drift": (4.18, 3.76),
    }
    points = []
    for row in frame.itertuples(index=False):
        dx, dy = jitter[row.Attack]
        points.append((row.Likelihood + dx, row.Impact + dy))
    ax.scatter([point[0] for point in points], [point[1] for point in points], s=170, c=colors, edgecolor="white", lw=1.2, zorder=3)
    for row, point in zip(frame.itertuples(index=False), points, strict=True):
        text_x, text_y = label_positions[row.Attack]
        ax.annotate(row.Attack, xy=point, xytext=(text_x, text_y), fontsize=8.0, color="#263238", arrowprops={"arrowstyle": "-", "color": "#90A4AE", "lw": 0.8}, bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.82})
    ax.axvspan(4.5, 5.35, color="#FDECEC", zorder=0)
    ax.axhspan(4.5, 5.35, color="#FDECEC", zorder=0)
    ax.set_xlim(2.8, 5.25)
    ax.set_ylim(3.55, 5.35)
    ax.set_xticks([3, 4, 5], ["Medium", "High", "Very high"])
    ax.set_yticks([4, 5], ["Major", "Critical"])
    ax.set_xlabel("Likelihood of reviewer challenge")
    ax.set_ylabel("Impact if unresolved")
    ax.grid(color="#EDF1F2", zorder=0)
    ax.set_title("Design the supplement around predictable reviewer attacks", loc="left")
    fig.text(0.01, 0.01, "Red points can overturn the central claim; gold points mainly narrow interpretation.", fontsize=8.3, color="#51636B")
    return save(fig, output, "75-reviewer-attack-map")


def main() -> None:
    args = parse_args()
    data = args.input_dir.resolve()
    output = args.figures_dir.resolve()
    np.random.seed(PLOT_SEED)
    style()
    generated: list[Path] = []
    for function in (
        paper_arc,
        five_figure_storyboard,
        panel_budget,
        main_supplement,
        claim_evidence,
        traceability_audit,
        sensitivity_matrix,
        style_contract,
        reviewer_attack_map,
    ):
        generated.extend(function(data, output))

    anchor = output / "75-wirbel-figure1-original.jpg"
    shutil.copy2(data / "wirbel-figure1-original.jpg", anchor)
    generated.append(anchor)
    manifest = {
        "article": 75,
        "plot_seed": PLOT_SEED,
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in sorted(generated)
        },
    }
    (data / "figure-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"plotted\t{len(generated)} files\t{output}")


if __name__ == "__main__":
    main()
