#!/usr/bin/env python3
"""Create deterministic publication figures for Article 67."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


SEED = 20_260_767
COLORS = {"Control": "#4C78A8", "CD": "#E45756", "UC": "#72B7B2"}
ORDER = ("Control", "CD", "UC")
FIGURES = (
    "67-sampling-design",
    "67-dysbiosis-trajectories",
    "67-lag-stability",
    "67-short-interval-shifts",
    "67-species-retention",
    "67-prevotella-trajectories",
    "67-mixed-model-effects",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    return parser.parse_args()


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.dpi": 120,
            "savefig.dpi": 360,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    sns.set_style("ticks")


def save(fig: plt.Figure, output: Path, stem: str) -> None:
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output / f"{stem}.png", dpi=360, bbox_inches="tight")
    fig.savefig(
        output / f"{stem}.tiff",
        dpi=360,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.13, 1.06, label, transform=axis.transAxes, fontweight="bold", fontsize=11)


def sampling_design(input_dir: Path, output: Path, rng: np.random.Generator) -> None:
    attrition = pd.read_csv(input_dir / "sample-attrition.tsv", sep="\t")
    sample = pd.read_csv(input_dir / "sample-ledger.tsv", sep="\t")
    visits = sample.groupby(["SubjectID", "Diagnosis"]).size().rename("Visits").reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    bars = ax.barh(
        np.arange(len(attrition)),
        attrition["Profiles"],
        color=["#B8C4CE", "#88A6B8", "#5E8EA8", "#2F6F91"],
        height=0.68,
    )
    ax.set_yticks(np.arange(len(attrition)), attrition["Stage"])
    ax.invert_yaxis()
    ax.set_xlabel("Profiles retained")
    ax.set_title("Deterministic sample attrition", loc="left", fontweight="bold")
    ax.set_xlim(0, 1770)
    for bar, profiles, subjects in zip(bars, attrition["Profiles"], attrition["Subjects"], strict=True):
        ax.text(profiles + 22, bar.get_y() + bar.get_height() / 2, f"{profiles:,} / {subjects} subjects", va="center", fontsize=8)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    panel_label(ax, "A")

    ax = axes[1]
    positions = np.arange(len(ORDER))
    for position, diagnosis in enumerate(ORDER):
        values = visits.loc[visits["Diagnosis"].eq(diagnosis), "Visits"].to_numpy(float)
        violin = ax.violinplot(values, positions=[position], widths=0.72, showextrema=False)
        for body in violin["bodies"]:
            body.set_facecolor(COLORS[diagnosis])
            body.set_edgecolor("none")
            body.set_alpha(0.28)
        jitter = rng.uniform(-0.13, 0.13, size=len(values))
        ax.scatter(position + jitter, values, s=17, color=COLORS[diagnosis], alpha=0.75, edgecolor="white", linewidth=0.35)
        ax.plot([position - 0.22, position + 0.22], [np.median(values)] * 2, color="#222222", linewidth=1.6)
        ax.text(position, values.max() + 0.8, f"n={len(values)}", ha="center", fontsize=8)
    ax.set_xticks(positions, ORDER)
    ax.set_ylabel("Eligible visits per subject")
    ax.set_title("Repeated observations are the unit structure", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.6)
    panel_label(ax, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "67-sampling-design")


def dysbiosis_trajectories(input_dir: Path, output: Path) -> None:
    sample = pd.read_csv(input_dir / "sample-ledger.tsv", sep="\t")
    predictions = pd.read_csv(input_dir / "primary-marginal-predictions.tsv", sep="\t")
    threshold = pd.read_csv(input_dir / "dysbiosis-reference-audit.tsv", sep="\t").iloc[0]["DysbiosisThreshold"]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8), sharex=True, sharey=True)
    for index, (axis, diagnosis) in enumerate(zip(axes, ORDER, strict=True)):
        group = sample.loc[sample["Diagnosis"].eq(diagnosis)].sort_values(["SubjectID", "Week"])
        for _, subject in group.groupby("SubjectID"):
            axis.plot(subject["Week"], subject["DysbiosisScore"], color=COLORS[diagnosis], alpha=0.16, linewidth=0.65)
        antibiotic = group.loc[group["Antibiotics"].eq("Yes")]
        axis.scatter(
            antibiotic["Week"], antibiotic["DysbiosisScore"], marker="x", s=13,
            color="#B27900", linewidth=0.75, alpha=0.75, label="Antibiotics: Yes" if index == 0 else None,
        )
        curve = predictions.loc[predictions["Diagnosis"].eq(diagnosis)].sort_values("Week")
        axis.fill_between(curve["Week"], curve["CILower"], curve["CIUpper"], color=COLORS[diagnosis], alpha=0.22, linewidth=0)
        axis.plot(curve["Week"], curve["Predicted"], color=COLORS[diagnosis], linewidth=2.2, label="Fixed-effect prediction" if index == 0 else None)
        axis.axhline(threshold, color="#333333", linestyle="--", linewidth=0.9, label="Control 90th percentile" if index == 0 else None)
        axis.set_title(f"{diagnosis} · {group['SubjectID'].nunique()} subjects", color=COLORS[diagnosis], fontweight="bold")
        axis.set_xlabel("Study week")
        axis.grid(color="#ECECEC", linewidth=0.55)
        panel_label(axis, chr(ord("A") + index))
    axes[0].set_ylabel("Dysbiosis score\n(median Bray–Curtis to control reference)")
    axes[0].legend(frameon=False, loc="upper left")
    fig.tight_layout(w_pad=1.2)
    save(fig, output, "67-dysbiosis-trajectories")


def lag_stability(input_dir: Path, output: Path) -> None:
    summary = pd.read_csv(input_dir / "lag-summary.tsv", sep="\t")
    labels = ["0–2", ">2–4", ">4–8", ">8–16", ">16–32", ">32–60"]
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    x = np.arange(len(labels))
    offsets = {"Control": -0.10, "CD": 0, "UC": 0.10}
    for diagnosis in ORDER:
        group = summary.loc[summary["Diagnosis"].eq(diagnosis)].set_index("LagBin").loc[labels]
        y = group["SubjectMedianBrayCurtis"].to_numpy(float)
        lower = y - group["CILower"].to_numpy(float)
        upper = group["CIUpper"].to_numpy(float) - y
        ax.errorbar(
            x + offsets[diagnosis], y, yerr=np.vstack([lower, upper]),
            color=COLORS[diagnosis], marker="o", markersize=5, linewidth=1.8,
            elinewidth=1.0, capsize=2.5, label=diagnosis,
        )
    ax.set_xticks(x, labels)
    ax.set_xlabel("Within-subject time lag (weeks)")
    ax.set_ylabel("Subject-equal median Bray–Curtis")
    ax.set_title("Community dissimilarity rises with elapsed time", loc="left", fontweight="bold")
    ax.text(0.01, 0.02, "Points: median across subject-level medians · Bars: subject bootstrap 95% CI", transform=ax.transAxes, fontsize=8, color="#555555")
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.6)
    fig.tight_layout()
    save(fig, output, "67-lag-stability")


def short_interval_shifts(input_dir: Path, output: Path, rng: np.random.Generator) -> None:
    subject = pd.read_csv(input_dir / "subject-shift-summary.tsv", sep="\t")
    summary = pd.read_csv(input_dir / "shift-summary.tsv", sep="\t").set_index("Diagnosis")
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    for position, diagnosis in enumerate(ORDER):
        values = subject.loc[subject["Diagnosis"].eq(diagnosis), "ShiftFraction"].to_numpy(float) * 100
        jitter = rng.uniform(-0.16, 0.16, size=len(values))
        ax.scatter(position + jitter, values, s=22, color=COLORS[diagnosis], alpha=0.67, edgecolor="white", linewidth=0.4)
        estimate = summary.loc[diagnosis, "MeanSubjectShiftFraction"] * 100
        lower = summary.loc[diagnosis, "MeanSubjectCILower"] * 100
        upper = summary.loc[diagnosis, "MeanSubjectCIUpper"] * 100
        ax.errorbar(position, estimate, yerr=[[estimate - lower], [upper - estimate]], color="#202020", marker="D", markersize=5, capsize=4, linewidth=1.3)
        ax.text(
            position, 103,
            f"{int(summary.loc[diagnosis, 'ShiftEvents'])}/{int(summary.loc[diagnosis, 'ShortIntervals'])} intervals",
            ha="center", fontsize=8,
        )
    ax.set_xticks(np.arange(len(ORDER)), ORDER)
    ax.set_ylabel("Subject-level shift fraction (%)")
    ax.set_ylim(-3, 110)
    ax.set_title("A fixed 0.54 threshold is a sensitivity analysis", loc="left", fontweight="bold")
    fig.text(0.02, 0.012, "Diamonds: mean across subjects · Bars: subject bootstrap 95% CI", fontsize=8, color="#555555")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.6)
    fig.tight_layout(rect=[0, 0.045, 1, 1])
    save(fig, output, "67-short-interval-shifts")


def species_retention(input_dir: Path, output: Path) -> None:
    summary = pd.read_csv(input_dir / "retention-summary.tsv", sep="\t")
    bands = ["0.01–0.1%", "0.1–1%", "≥1%"]
    x = np.arange(len(bands))
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    offsets = {"Control": -0.08, "CD": 0, "UC": 0.08}
    for diagnosis in ORDER:
        group = summary.loc[summary["Diagnosis"].eq(diagnosis)].set_index("BaselineAbundance").loc[bands]
        y = group["SubjectMedianRetention"].to_numpy(float) * 100
        lower = y - group["CILower"].to_numpy(float) * 100
        upper = group["CIUpper"].to_numpy(float) * 100 - y
        ax.errorbar(
            x + offsets[diagnosis], y, yerr=np.vstack([lower, upper]),
            color=COLORS[diagnosis], marker="o", linewidth=1.8, capsize=3,
            markersize=5, label=diagnosis,
        )
    ax.set_xticks(x, bands)
    ax.set_xlabel("Baseline relative-abundance band")
    ax.set_ylabel("Retained at next 1–3-week visit (%)")
    ax.set_ylim(45, 101)
    ax.set_title("Detection stability depends strongly on starting abundance", loc="left", fontweight="bold")
    fig.text(0.02, 0.012, "Retention: detected at ≥0.01% at both visits · Subject-equal medians and bootstrap 95% CIs", fontsize=8, color="#555555")
    ax.legend(frameon=False, ncol=3, loc="lower right")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.6)
    fig.tight_layout(rect=[0, 0.045, 1, 1])
    save(fig, output, "67-species-retention")


def prevotella_trajectories(input_dir: Path, output: Path) -> None:
    data = pd.read_csv(input_dir / "prevotella-selected-trajectories.tsv", sep="\t")
    fig, axes = plt.subplots(4, 3, figsize=(10.2, 9.0), sharex=True, sharey=True)
    for column, diagnosis in enumerate(ORDER):
        subjects = (
            data.loc[data["Diagnosis"].eq(diagnosis), ["SubjectID", "SelectionRank"]]
            .drop_duplicates()
            .sort_values("SelectionRank")["SubjectID"]
            .tolist()
        )
        for row, subject_id in enumerate(subjects):
            axis = axes[row, column]
            subject = data.loc[data["SubjectID"].eq(subject_id)].sort_values("Week")
            abundance = subject["RelativeAbundance"].to_numpy(float) * 100
            axis.plot(subject["Week"], abundance, color=COLORS[diagnosis], marker="o", markersize=2.8, linewidth=1.2)
            antibiotic = subject["Antibiotics"].eq("Yes").to_numpy()
            if antibiotic.any():
                axis.scatter(subject.loc[antibiotic, "Week"], abundance[antibiotic], marker="x", color="#B27900", s=24, linewidth=1.0)
            axis.set_yscale("symlog", linthresh=0.01, linscale=0.7)
            axis.set_ylim(-0.001, 100)
            axis.set_yticks([0, 0.01, 0.1, 1, 10, 100])
            axis.set_yticklabels(["0", "0.01", "0.1", "1", "10", "100"])
            axis.text(0.03, 0.88, subject_id, transform=axis.transAxes, fontsize=8, fontweight="bold")
            axis.grid(axis="y", color="#ECECEC", linewidth=0.5)
            if row == 0:
                axis.set_title(diagnosis, color=COLORS[diagnosis], fontweight="bold")
            if row == 3:
                axis.set_xlabel("Study week")
            if column == 0:
                axis.set_ylabel("P. copri (%)")
    fig.suptitle("Predeclared descriptive view: four highest-variance eligible subjects per diagnosis", x=0.08, ha="left", y=0.995, fontsize=11, fontweight="bold")
    fig.text(0.99, 0.005, "Eligibility: ≥8 visits and ≥20% detection; selection is descriptive, not inferential", ha="right", fontsize=8, color="#555555")
    fig.tight_layout(rect=[0, 0.02, 1, 0.975])
    save(fig, output, "67-prevotella-trajectories")


def mixed_model_effects(input_dir: Path, output: Path) -> None:
    primary = pd.read_csv(input_dir / "primary-fixed-effects.tsv", sep="\t")
    primary_labels = {
        "DiagnosisCD": "CD vs Control at week 26",
        "DiagnosisUC": "UC vs Control at week 26",
        "WeekYearCentered": "Time in Control (per year)",
        "AntibioticsYes": "Antibiotics: Yes vs No",
        "Log10ReadsCentered": "Log10 filtered reads",
        "DiagnosisCD:WeekYearCentered": "CD × time",
        "DiagnosisUC:WeekYearCentered": "UC × time",
    }
    primary = primary.loc[primary["Term"].isin(primary_labels)].copy()
    primary["Label"] = primary["Term"].map(primary_labels)
    primary = primary.iloc[::-1]
    secondary = pd.read_csv(input_dir / "species-mixed-model-results.tsv", sep="\t").head(12).copy()
    secondary["Label"] = secondary["Species"].str.replace("s__", "", regex=False).str.replace("_", " ", regex=False)
    term_labels = {
        "AntibioticsYes": "Antibiotics",
        "DiagnosisCD:WeekYearCentered": "CD × time",
        "DiagnosisUC:WeekYearCentered": "UC × time",
        "DiagnosisCD": "CD",
        "DiagnosisUC": "UC",
    }
    term_colors = {
        "AntibioticsYes": "#B27900",
        "DiagnosisCD:WeekYearCentered": COLORS["CD"],
        "DiagnosisUC:WeekYearCentered": COLORS["UC"],
        "DiagnosisCD": "#9D3E3B",
        "DiagnosisUC": "#3B827E",
    }
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.3), gridspec_kw={"width_ratios": [0.9, 1.35]})
    ax = axes[0]
    y = np.arange(len(primary))
    ax.axvline(0, color="#888888", linewidth=0.8)
    ax.errorbar(
        primary["Estimate"], y,
        xerr=np.vstack([primary["Estimate"] - primary["CILower"], primary["CIUpper"] - primary["Estimate"]]),
        fmt="o", color="#2F6F91", ecolor="#2F6F91", capsize=3, markersize=5,
    )
    ax.set_yticks(y, primary["Label"])
    ax.set_xlabel("Change in dysbiosis score (95% CI)")
    ax.set_title("Primary random-slope LMM", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    panel_label(ax, "A")

    ax = axes[1]
    secondary = secondary.iloc[::-1]
    y = np.arange(len(secondary))
    ax.axvline(0, color="#888888", linewidth=0.8)
    for index, row in enumerate(secondary.itertuples(index=False)):
        color = term_colors.get(row.Term, "#666666")
        ax.errorbar(row.Estimate, index, xerr=[[row.Estimate - row.CILower], [row.CIUpper - row.Estimate]], fmt="o", color=color, ecolor=color, capsize=2.5, markersize=4.5)
    ax.set_yticks(y, secondary["Label"])
    ax.set_xlabel("CLR coefficient (95% CI)")
    ax.set_title("Top exploratory species coefficients", loc="left", fontweight="bold")
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", color=color, label=term_labels[term])
        for term, color in term_colors.items()
        if term in set(secondary["Term"])
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=2,
        title="Model term",
    )
    ax.text(0.01, 0.01, "Ranked by global BH q-value across 315 tests", transform=ax.transAxes, fontsize=8, color="#555555")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    panel_label(ax, "B")
    fig.tight_layout(w_pad=2.1, rect=[0, 0.08, 1, 1])
    save(fig, output, "67-mixed-model-effects")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output = args.figure_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure()
    rng = np.random.default_rng(SEED)
    sampling_design(input_dir, output, rng)
    dysbiosis_trajectories(input_dir, output)
    lag_stability(input_dir, output)
    short_interval_shifts(input_dir, output, rng)
    species_retention(input_dir, output)
    prevotella_trajectories(input_dir, output)
    mixed_model_effects(input_dir, output)
    shutil.copy2(input_dir / "lloyd-price-fig3-original.png", output / "67-lloyd-price-fig3-original.png")
    print(f"Wrote {len(FIGURES)} figure sets and the article anchor to {output}")


if __name__ == "__main__":
    main()
