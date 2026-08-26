#!/usr/bin/env python3
"""Create deterministic publication figures for Article 72."""

from __future__ import annotations

import argparse
import shutil
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PLOT_SEED = 20_260_772
COLORS = {
    "Association": "#9AA5B1",
    "Temporal": "#76A5AF",
    "Causal": "#B279A2",
    "Human": "#4C78A8",
    "Culture": "#F2C14E",
    "Transfer": "#59A14F",
    "Mechanism": "#D55E00",
    "Direct": "#2A9D8F",
    "Supporting": "#F2C14E",
    "Absent": "#F2F4F5",
    "Risk": "#C44E52",
    "Ink": "#33434A",
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
            "legend.fontsize": 7.4,
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


def evidence_ladder(input_dir: Path, output: Path) -> None:
    rungs = pd.read_csv(input_dir / "rung-definitions.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.9), gridspec_kw={"width_ratios": [1.15, 0.85]})

    axis = axes[0]
    axis.set_xlim(0, 8.6)
    axis.set_ylim(0, 8.3)
    axis.axis("off")
    short = [
        "Replicated\nassociation",
        "Temporality &\nreversibility",
        "Human causal\nbridge",
        "Randomized human\nintervention",
        "Isolation /\ndefined product",
        "Host transfer /\nperturbation",
        "Molecular perturbation\n& rescue",
    ]
    palette = [
        COLORS["Association"], COLORS["Temporal"], COLORS["Causal"],
        COLORS["Human"], COLORS["Culture"], COLORS["Transfer"], COLORS["Mechanism"],
    ]
    for index, (label, color) in enumerate(zip(short, palette, strict=True), start=1):
        x = 0.22 + (index - 1) * 0.92
        y = 0.55 + (index - 1) * 0.92
        width = 2.15
        add_box(axis, x, y, width, 0.72, label, color, fontsize=7.2)
        axis.text(x + 0.12, y + 0.55, str(index), fontweight="bold", color="#243238")
    axis.annotate(
        "",
        xy=(7.72, 7.70),
        xytext=(0.35, 0.40),
        arrowprops={"arrowstyle": "->", "color": "#59636A", "linewidth": 1.2},
    )
    axis.text(4.0, 0.06, "Increasing intervention and entity specificity", ha="center", color="#59636A")
    axis.set_title("Seven questions—not seven points", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    domains = [
        ("Human relevance", 0.84, COLORS["Human"]),
        ("Time ordering", 0.70, COLORS["Temporal"]),
        ("Intervention", 0.56, COLORS["Causal"]),
        ("Entity specificity", 0.42, COLORS["Culture"]),
        ("Host transport", 0.28, COLORS["Transfer"]),
        ("Mechanism", 0.14, COLORS["Mechanism"]),
    ]
    for label, y, color in domains:
        add_box(axis, 0.04, y - 0.045, 0.34, 0.09, label, color, fontsize=7.3)
        axis.add_patch(
            FancyArrowPatch(
                (0.39, y), (0.62, 0.50), arrowstyle="-|>",
                mutation_scale=10, linewidth=1.0, color=color,
            )
        )
    add_box(axis, 0.63, 0.37, 0.31, 0.26, "Claim contract\nPopulation\nIntervention / exposure\nOutcome + horizon", "#E8EEF1", fontsize=7.6)
    axis.text(
        0.50,
        0.02,
        "A strong packet may skip a rung, but it may not\nborrow evidence from a different claim.",
        ha="center",
        color=COLORS["Risk"],
        fontweight="bold",
    )
    axis.set_title("Evidence is a braid anchored to one claim", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "72-evidence-ladder-framework")


def claim_contracts(input_dir: Path, output: Path) -> None:
    claims = pd.read_csv(input_dir / "claim-contracts.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.7))
    colors = [COLORS["Human"], COLORS["Transfer"], COLORS["Mechanism"]]
    short_claims = [
        "Tested microbiota interventions\ncan reduce short-term recurrent CDI\nin trial-eligible adults.",
        "C. scindens can increase\nC. difficile resistance in specified\nantibiotic-perturbed model systems.",
        "F. nucleatum can promote CRC\nphenotypes in specified cell, mouse,\nand xenograft systems.",
    ]
    for axis, row, color, short_claim in zip(axes, claims.itertuples(), colors, short_claims, strict=True):
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.axis("off")
        axis.text(0.04, 0.95, row.ClaimID, fontsize=13, fontweight="bold", color=color, va="top")
        axis.text(0.16, 0.95, row.Case, fontweight="bold", va="top")
        add_box(axis, 0.05, 0.64, 0.90, 0.20, short_claim, "#EEF3F5", fontsize=8.1)
        axis.text(0.05, 0.56, "Permitted", color=COLORS["Direct"], fontweight="bold")
        axis.text(
            0.05,
            0.51,
            "\n".join(textwrap.wrap(row.AllowedConclusion, width=45)),
            va="top",
            fontsize=7.3,
        )
        axis.text(0.05, 0.25, "Forbidden leap", color=COLORS["Risk"], fontweight="bold")
        axis.text(
            0.05,
            0.20,
            "\n".join(textwrap.wrap(row.ForbiddenLeap, width=45)),
            va="top",
            fontsize=7.3,
        )
    fig.suptitle("Three neighboring claims require three separate evidence packets", x=0.02, ha="left", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=1.7)
    save(fig, output, "72-claim-contracts")


def evidence_coverage(input_dir: Path, output: Path) -> None:
    coverage = pd.read_csv(input_dir / "study-domain-coverage.tsv", sep="\t")
    publications = pd.read_csv(input_dir / "publication-metadata.tsv", sep="\t")
    rungs = pd.read_csv(input_dir / "rung-definitions.tsv", sep="\t")
    evidence = pd.read_csv(input_dir / "evidence-ledger.tsv", sep="\t")
    study_order = evidence["EvidenceID"].tolist()
    rung_order = rungs["Rung"].tolist()
    matrix = coverage.pivot(index="EvidenceID", columns="Rung", values="CoverageScore").loc[study_order, rung_order]
    metadata = evidence.set_index("EvidenceID")
    author = publications.set_index("CitationKey")["FirstAuthor"]
    labels = [
        f"{metadata.loc[key, 'Year']} {author.loc[metadata.loc[key, 'CitationKey']]} ({metadata.loc[key, 'ClaimID']})"
        for key in study_order
    ]
    fig, axis = plt.subplots(figsize=(11.8, 6.0))
    cmap = ListedColormap([COLORS["Absent"], COLORS["Supporting"], COLORS["Direct"]])
    image = axis.imshow(matrix.to_numpy(), cmap=cmap, vmin=-0.5, vmax=2.5, aspect="auto")
    axis.set_xticks(
        np.arange(len(rung_order)),
        ["Association", "Time", "Mediation / MR", "Human RCT", "Defined entity", "Host transfer", "Mechanism"],
        rotation=25,
        ha="right",
    )
    axis.set_yticks(np.arange(len(labels)), labels)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = int(matrix.iloc[i, j])
            axis.text(j, i, {0: "", 1: "S", 2: "D"}[value], ha="center", va="center", fontweight="bold", color="#243238")
    for boundary in (3.5, 5.5):
        axis.axhline(boundary, color="white", linewidth=4)
        axis.axhline(boundary, color="#67747B", linewidth=0.8)
    axis.text(6.64, 1.5, "C1 clinical restoration", rotation=90, va="center", color=COLORS["Human"], fontweight="bold")
    axis.text(6.64, 4.5, "C2 bile-acid resistance", rotation=90, va="center", color=COLORS["Transfer"], fontweight="bold")
    axis.text(6.64, 7.5, "C3 CRC phenotypes", rotation=90, va="center", color=COLORS["Mechanism"], fontweight="bold")
    axis.set_title("Coverage is claim-specific and non-additive", loc="left", fontweight="bold")
    axis.text(0, -0.16, "D = direct for the prespecified claim; S = supporting; blank = not addressed", transform=axis.transAxes, color="#59636A")
    fig.tight_layout()
    save(fig, output, "72-evidence-coverage")


def human_intervention(input_dir: Path, output: Path) -> None:
    data = pd.read_csv(input_dir / "human-intervention-outcomes.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4))
    for axis, study, panel in zip(axes, data["Study"].unique(), ("A", "B"), strict=True):
        subset = data.loc[data["Study"].eq(study)].copy()
        colors = [COLORS["Human"]] + [COLORS["Association"]] * (len(subset) - 1)
        bars = axis.barh(np.arange(len(subset))[::-1], 100 * subset["FavorableRate"], color=colors, height=0.62)
        axis.set_yticks(np.arange(len(subset))[::-1], subset["Arm"])
        axis.set_xlim(0, 105)
        axis.set_xlabel(f"{subset.iloc[0]['FavorableEndpoint']} (%)")
        axis.set_title(f"{study}: randomized clinical endpoint", loc="left", fontweight="bold")
        axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
        for bar, row in zip(bars, subset.itertuples(), strict=True):
            axis.text(
                min(bar.get_width() + 2, 98),
                bar.get_y() + bar.get_height() / 2,
                f"{row.FavorableEvents}/{row.Total}",
                va="center",
                fontsize=7.7,
                fontweight="bold",
            )
        axis.text(
            0.02,
            0.03,
            f"Endpoint at {int(subset.iloc[0]['TimeWeeks'])} weeks",
            transform=axis.transAxes,
            color="#59636A",
        )
        panel_label(axis, panel)
    axes[1].text(
        0.98,
        0.03,
        "Recurrence RR 0.32\n95% CI 0.18–0.58",
        transform=axes[1].transAxes,
        ha="right",
        color=COLORS["Human"],
        fontweight="bold",
    )
    fig.tight_layout(w_pad=2.5)
    save(fig, output, "72-human-intervention")


def rcdi_braid(input_dir: Path, output: Path) -> None:
    fig, axis = plt.subplots(figsize=(11.7, 5.0))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    axis.text(0.02, 0.94, "C1 · Clinical intervention claim", color=COLORS["Human"], fontweight="bold", fontsize=10)
    axis.text(0.02, 0.45, "C2 · C. scindens mechanism claim", color=COLORS["Transfer"], fontweight="bold", fontsize=10)

    upper = [
        (0.03, "Randomized\nclinical effect\nvan Nood 2013", COLORS["Human"]),
        (0.27, "Randomized defined\nspore product\nFeuerstadt 2022", COLORS["Human"]),
        (0.51, "Shotgun strain\nengraftment\nSmillie 2018", COLORS["Temporal"]),
        (0.75, "C1 permitted claim\nintervention bundle\nreduces recurrence", "#DCE7F2"),
    ]
    lower = [
        (0.03, "Antibiotic ecology\n& bile acids\nTheriot 2014", COLORS["Association"]),
        (0.27, "Cultured\nC. scindens\nBuffie 2015", COLORS["Culture"]),
        (0.51, "Mouse transfer +\nsecondary bile acids\n+ ex-vivo block", COLORS["Transfer"]),
        (0.75, "C2 permitted claim\nmodel-system\ncolonization resistance", "#DDEEDC"),
    ]
    for lane, y in ((upper, 0.67), (lower, 0.18)):
        for index, (x, label, color) in enumerate(lane):
            add_box(axis, x, y, 0.19, 0.17, label, color, fontsize=7.2)
            if index < len(lane) - 1:
                axis.add_patch(FancyArrowPatch((x + 0.195, y + 0.085), (lane[index + 1][0] - 0.005, y + 0.085), arrowstyle="-|>", mutation_scale=10, linewidth=1.1, color="#59636A"))
    axis.plot([0.72, 0.72], [0.42, 0.60], color=COLORS["Risk"], linewidth=2)
    axis.text(
        0.735,
        0.51,
        "Do not merge:\nclinical efficacy ≠\nsingle-strain mechanism",
        va="center",
        color=COLORS["Risk"],
        fontweight="bold",
    )
    axis.set_title("The rCDI evidence braid contains two neighboring—but different—claims", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, output, "72-rcdi-evidence-braid")


def fusobacterium_gap(input_dir: Path, output: Path) -> None:
    evidence = pd.read_csv(input_dir / "evidence-ledger.tsv", sep="\t")
    crc = evidence.loc[evidence["ClaimID"].eq("C3")].copy()
    publications = pd.read_csv(input_dir / "publication-metadata.tsv", sep="\t").set_index("CitationKey")
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), gridspec_kw={"width_ratios": [1.3, 0.7]})

    axis = axes[0]
    x = np.arange(len(crc))
    domains = [
        ("Human association", "Replicated association", COLORS["Association"]),
        ("Defined entity", "Isolation or defined product", COLORS["Culture"]),
        ("Host perturbation", "Host transfer or perturbation", COLORS["Transfer"]),
        ("Mechanism / rescue", "Molecular perturbation and rescue", COLORS["Mechanism"]),
        ("Human targeted RCT", "Randomized human intervention", COLORS["Human"]),
    ]
    bottom = np.zeros(len(crc))
    for label, column, color in domains:
        values = crc[column].map({"Not addressed": 0, "Supporting": 0.5, "Direct": 1}).to_numpy(float)
        axis.bar(x, values, bottom=bottom, label=label, color=color, width=0.68)
        bottom += values
    labels = [f"{row.Year}\n{publications.loc[row.CitationKey, 'FirstAuthor']}" for row in crc.itertuples()]
    axis.set_xticks(x, labels)
    axis.set_ylabel("Evidence-domain coverage\n(direct=1, supporting=0.5; not a causal score)")
    axis.set_title("F. nucleatum evidence spans models, not a human targeted trial", loc="left", fontweight="bold")
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    axis.legend(frameon=False, ncol=2, loc="upper left")
    panel_label(axis, "A")

    axis = axes[1]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    add_box(axis, 0.10, 0.68, 0.80, 0.17, "Replicated human\nassociation", "#E6EAED", fontsize=8.2)
    add_box(axis, 0.10, 0.43, 0.80, 0.17, "Mouse / cell / xenograft\nperturbation and mechanism", "#E2EFE4", fontsize=8.2)
    axis.add_patch(FancyArrowPatch((0.50, 0.68), (0.50, 0.61), arrowstyle="-|>", mutation_scale=11, color="#59636A"))
    axis.add_patch(FancyArrowPatch((0.50, 0.43), (0.50, 0.33), arrowstyle="-|>", mutation_scale=11, color="#59636A"))
    add_box(axis, 0.10, 0.12, 0.80, 0.17, "Missing for a clinical-effect claim:\nrandomized targeted human intervention", "#F8DAD8", edge=COLORS["Risk"], fontsize=8.2)
    axis.text(0.50, 0.04, "Permitted verb: can promote in specified models", ha="center", color=COLORS["Mechanism"], fontweight="bold")
    axis.set_title("The gap controls the verb", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "72-fusobacterium-gap")


def claim_downgrade(input_dir: Path, output: Path) -> None:
    data = pd.read_csv(input_dir / "claim-downgrade-examples.tsv", sep="\t")
    fig, axis = plt.subplots(figsize=(11.8, 6.0))
    axis.set_xlim(0, 1)
    axis.set_ylim(-0.4, len(data) - 0.55)
    axis.axis("off")
    axis.text(0.02, len(data) - 0.62, "Overclaim", color=COLORS["Risk"], fontweight="bold", fontsize=10)
    axis.text(0.53, len(data) - 0.62, "Supported wording", color=COLORS["Direct"], fontweight="bold", fontsize=10)
    for index, row in enumerate(data.itertuples()):
        y = len(data) - 1.3 - index * 0.88
        add_box(axis, 0.02, y - 0.30, 0.40, 0.58, "\n".join(textwrap.wrap(row.Overclaim, width=42)), "#F8DEDC", edge=COLORS["Risk"], fontsize=7.2)
        axis.add_patch(FancyArrowPatch((0.43, y), (0.51, y), arrowstyle="-|>", mutation_scale=11, color="#59636A"))
        add_box(axis, 0.52, y - 0.30, 0.46, 0.58, "\n".join(textwrap.wrap(row.SupportedClaim, width=58)), "#DDEFEA", edge=COLORS["Direct"], fontsize=7.1)
    axis.set_title("The evidence packet determines the verb—not the P value", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, output, "72-claim-downgrade")


def evidence_packet(input_dir: Path, output: Path) -> None:
    packet = pd.read_csv(input_dir / "evidence-packet.tsv", sep="\t")
    fig, axes = plt.subplots(2, 4, figsize=(12.0, 5.7))
    palette = [
        COLORS["Human"], COLORS["Temporal"], COLORS["Association"], COLORS["Causal"],
        COLORS["Culture"], COLORS["Transfer"], COLORS["Mechanism"], "#7A8B93",
    ]
    for axis, row, color in zip(axes.flat, packet.itertuples(), palette, strict=True):
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.axis("off")
        axis.text(0.05, 0.92, f"{row.Order:02d}", fontsize=13, fontweight="bold", color=color, va="top")
        axis.text(0.20, 0.90, row.Packet, fontweight="bold", va="top")
        add_box(axis, 0.05, 0.32, 0.90, 0.43, "\n".join(textwrap.wrap(row.RequiredFields, width=35)), "#F2F5F6", fontsize=7.1)
        axis.text(0.05, 0.15, row.UsedFor, color="#59636A", fontsize=7.0, va="top")
    fig.suptitle("A submission-ready causal evidence packet has eight auditable ledgers", x=0.02, ha="left", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=1.3, h_pad=1.5)
    save(fig, output, "72-evidence-packet")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output = args.figure_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    np.random.seed(PLOT_SEED)
    style()
    shutil.copy2(
        input_dir / "buffie-figure4-original.jpg",
        output / "72-buffie-figure4-original.jpg",
    )
    evidence_ladder(input_dir, output)
    claim_contracts(input_dir, output)
    evidence_coverage(input_dir, output)
    human_intervention(input_dir, output)
    rcdi_braid(input_dir, output)
    fusobacterium_gap(input_dir, output)
    claim_downgrade(input_dir, output)
    evidence_packet(input_dir, output)


if __name__ == "__main__":
    main()
