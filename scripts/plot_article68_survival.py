#!/usr/bin/env python3
"""Create deterministic publication figures for Article 68."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess


SEED = 20_260_768
COLORS = {
    "Clinical": "#4C78A8",
    "Clinical + Faecalibacterium": "#E45756",
    "Low": "#9AA5B1",
    "High": "#2A9D8F",
    "Event": "#D55E00",
    "Censored": "#4C78A8",
}
FIGURES = (
    "68-cohort-design",
    "68-cox-forest",
    "68-assumption-and-spline",
    "68-cutoff-leakage",
    "68-time-dependent-roc",
    "68-crossvalidated-calibration",
    "68-influence-sensitivity",
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


def cohort_design(input_dir: Path, output: Path, rng: np.random.Generator) -> None:
    attrition = pd.read_csv(input_dir / "cohort-attrition.tsv", sep="\t")
    cohort = pd.read_csv(input_dir / "survival-cohort.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.7, 4.35), gridspec_kw={"width_ratios": [1.05, 1.35]})
    axis = axes[0]
    colors = ["#B7C6D0", "#82A9BC", "#4C86A8", "#2F6680"]
    bars = axis.barh(np.arange(len(attrition)), attrition["Samples"], color=colors, height=0.67)
    axis.set_yticks(np.arange(len(attrition)), attrition["Stage"])
    axis.invert_yaxis()
    axis.set_xlabel("Patients retained")
    axis.set_xlim(0, 184)
    axis.set_title("One baseline WGS sample per patient", loc="left", fontweight="bold")
    for bar, samples, events in zip(bars, attrition["Samples"], attrition["Events"], strict=True):
        axis.text(
            samples + 3,
            bar.get_y() + bar.get_height() / 2,
            f"{samples} / {events} events",
            va="center",
            fontsize=8,
        )
    axis.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    panel_label(axis, "A")

    axis = axes[1]
    status = np.where(cohort["Event"].eq(1), "Event", "Censored")
    for label in ("Censored", "Event"):
        mask = status == label
        jitter = rng.normal(0, 0.025, mask.sum())
        axis.scatter(
            cohort.loc[mask, "FaecalibacteriumLog2"],
            cohort.loc[mask, "PFS_months"] + jitter,
            s=27,
            color=COLORS[label],
            alpha=0.72,
            edgecolor="white",
            linewidth=0.45,
            label=f"{label} (n={mask.sum()})",
        )
    low_quality = ~cohort["QualitySensitivityPass"].astype(bool)
    axis.scatter(
        cohort.loc[low_quality, "FaecalibacteriumLog2"],
        cohort.loc[low_quality, "PFS_months"],
        s=80,
        facecolor="none",
        edgecolor="#111111",
        linewidth=1.2,
        label="Sequencing-QC sensitivity exclusion",
    )
    axis.set_xlabel(r"log$_2$(Faecalibacterium PPM + 25)")
    axis.set_ylabel("Observed PFS or censoring time (months)")
    axis.set_title("Censoring cannot be analyzed as an ordinary scatter", loc="left", fontweight="bold")
    axis.legend(frameon=False, loc="upper right")
    axis.grid(color="#EAEAEA", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "68-cohort-design")


def cox_forest(input_dir: Path, output: Path) -> None:
    estimates = pd.read_csv(input_dir / "cox-model-estimates.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.9, 4.55), gridspec_kw={"width_ratios": [1.0, 1.25]})

    axis = axes[0]
    order = ["Unadjusted", "Adjusted primary", "Sequencing-QC sensitivity", "Anti-PD1 sensitivity"]
    labels = {
        "Unadjusted": "Unadjusted",
        "Adjusted primary": "Adjusted primary",
        "Sequencing-QC sensitivity": "Sequencing-QC sensitivity",
        "Anti-PD1 sensitivity": "Anti-PD1 sensitivity",
    }
    data = estimates.loc[estimates["Term"].eq("FaecalibacteriumLog2")].set_index("Model").loc[order]
    y = np.arange(len(order))[::-1]
    axis.errorbar(
        data["HazardRatio"],
        y,
        xerr=np.vstack([data["HazardRatio"] - data["CILower"], data["CIUpper"] - data["HazardRatio"]]),
        fmt="o",
        color="#2F6F91",
        ecolor="#2F6F91",
        capsize=3,
        linewidth=1.4,
    )
    axis.axvline(1, color="#333333", linestyle="--", linewidth=0.9)
    axis.set_yticks(y, [labels[item] for item in order])
    axis.set_xlabel("Hazard ratio per abundance doubling (95% CI)")
    axis.set_xscale("log")
    axis.set_xlim(0.78, 1.16)
    axis.set_title("Predeclared microbiome feature", loc="left", fontweight="bold")
    for yy, (_, row) in zip(y, data.iterrows(), strict=True):
        axis.text(1.16, yy, f"P={row['PValue']:.3f}\nn={int(row['N'])}", va="center", fontsize=7.5)
    axis.grid(axis="x", color="#E6E6E6", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    primary = estimates.loc[estimates["Model"].eq("Adjusted primary")].copy()
    term_order = [
        "FaecalibacteriumLog2",
        "PrimarySubtypeMucosal_or_acral",
        "AdvancedSubstageStage_M1D",
        "LDHYes",
        "BMIz",
    ]
    term_labels = {
        "FaecalibacteriumLog2": "Faecalibacterium (per doubling)",
        "PrimarySubtypeMucosal_or_acral": "Mucosal/acral vs cutaneous/unknown",
        "AdvancedSubstageStage_M1D": "M1D vs M1C",
        "LDHYes": "LDH above upper limit vs not",
        "BMIz": "BMI (per SD)",
    }
    primary = primary.set_index("Term").loc[term_order]
    y = np.arange(len(term_order))[::-1]
    axis.errorbar(
        primary["HazardRatio"],
        y,
        xerr=np.vstack([primary["HazardRatio"] - primary["CILower"], primary["CIUpper"] - primary["HazardRatio"]]),
        fmt="o",
        color="#E45756",
        ecolor="#E45756",
        capsize=3,
        linewidth=1.4,
    )
    axis.axvline(1, color="#333333", linestyle="--", linewidth=0.9)
    axis.set_yticks(y, [term_labels[item] for item in term_order])
    axis.set_xscale("log")
    axis.set_xlabel("Adjusted hazard ratio (95% CI)")
    axis.set_title("Paper-aligned clinical adjustment", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E6E6E6", linewidth=0.55)
    panel_label(axis, "B")
    fig.text(
        0.02,
        0.01,
        "Primary model: melanoma subtype + advanced substage + LDH category + BMI; one baseline sample per patient.",
        fontsize=8,
        color="#4F5962",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1], w_pad=2.4)
    save(fig, output, "68-cox-forest")


def assumption_and_spline(input_dir: Path, output: Path) -> None:
    residual = pd.read_csv(input_dir / "faecalibacterium-schoenfeld.tsv", sep="\t")
    ph = pd.read_csv(input_dir / "proportional-hazards-tests.tsv", sep="\t")
    spline = pd.read_csv(input_dir / "spline-effect-curve.tsv", sep="\t")
    nonlinear = pd.read_csv(input_dir / "nonlinearity-test.tsv", sep="\t").iloc[0]
    cohort = pd.read_csv(input_dir / "survival-cohort.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.35))
    axis = axes[0]
    x = residual["EventTimeDays"].to_numpy(float) / 30.4375
    y = residual["ScaledSchoenfeld"].to_numpy(float)
    smooth = lowess(y, x, frac=0.45, return_sorted=True)
    axis.scatter(x, y, s=22, color="#7D8A95", alpha=0.6, edgecolor="white", linewidth=0.35)
    axis.plot(smooth[:, 0], smooth[:, 1], color="#D55E00", linewidth=2.0)
    axis.axhline(0, color="#333333", linestyle="--", linewidth=0.8)
    faec_p = ph.loc[ph["Term"].eq("FaecalibacteriumLog2"), "PValue"].iloc[0]
    global_p = ph.loc[ph["Term"].eq("GLOBAL"), "PValue"].iloc[0]
    axis.text(0.03, 0.95, f"Feature PH P={faec_p:.3f}\nGlobal PH P={global_p:.3f}", transform=axis.transAxes, va="top")
    axis.set_xlabel("Event time (months)")
    axis.set_ylabel("Scaled Schoenfeld residual")
    axis.set_title("Proportional-hazards diagnostic", loc="left", fontweight="bold")
    axis.grid(color="#ECECEC", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    axis.fill_between(
        spline["FaecalibacteriumLog2"],
        spline["CILower"],
        spline["CIUpper"],
        color="#4C78A8",
        alpha=0.20,
        linewidth=0,
    )
    axis.plot(spline["FaecalibacteriumLog2"], spline["HazardRatio"], color="#4C78A8", linewidth=2.0)
    axis.axhline(1, color="#333333", linestyle="--", linewidth=0.8)
    rug_y = np.full(len(cohort), axis.get_ylim()[0] * 1.02)
    axis.plot(cohort["FaecalibacteriumLog2"], rug_y, "|", color="#303030", alpha=0.35, markersize=5)
    axis.text(0.03, 0.95, f"Nonlinearity LRT P={nonlinear['PValue']:.3f}", transform=axis.transAxes, va="top")
    axis.set_xlabel(r"log$_2$(Faecalibacterium PPM + 25)")
    axis.set_ylabel("Adjusted hazard ratio vs cohort median")
    axis.set_yscale("log")
    axis.set_title("Do not force a linear effect without checking", loc="left", fontweight="bold")
    axis.grid(axis="y", color="#ECECEC", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "68-assumption-and-spline")


def plot_km(axis: plt.Axes, curves: pd.DataFrame, audit: pd.Series, title: str, leaky: bool) -> None:
    for group in ("Low", "High"):
        data = curves.loc[curves["Group"].eq(group)].sort_values("TimeMonths")
        axis.step(data["TimeMonths"], data["Survival"], where="post", color=COLORS[group], linewidth=2.0, label=f"{group} (n={int(audit[group + 'N'])})")
        axis.fill_between(
            data["TimeMonths"], data["CILower"], data["CIUpper"],
            step="post", color=COLORS[group], alpha=0.13, linewidth=0,
        )
    axis.set_xlim(0, 62)
    axis.set_ylim(0, 1.03)
    axis.set_xlabel("Progression-free survival (months)")
    axis.set_ylabel("Survival probability")
    axis.set_title(title, loc="left", fontweight="bold")
    prefix = "Leaky descriptive" if leaky else "Held-out"
    axis.text(
        0.03,
        0.08,
        f"{prefix} log-rank P={audit['LogRankPValue']:.3f}\nHR={audit['HazardRatioHighVsLow']:.2f}",
        transform=axis.transAxes,
        fontsize=8,
    )
    axis.legend(frameon=False, loc="upper right")
    axis.grid(color="#EEEEEE", linewidth=0.5)


def cutoff_leakage(input_dir: Path, output: Path) -> None:
    train = pd.read_csv(input_dir / "cutoff-search-training.tsv", sep="\t")
    full = pd.read_csv(input_dir / "cutoff-search-full-leaky.tsv", sep="\t")
    audit = pd.read_csv(input_dir / "cutoff-evaluation.tsv", sep="\t")
    curves = pd.read_csv(input_dir / "cutoff-km-curves.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.25))
    axis = axes[0]
    for data, label, color in (
        (train, "Training only", "#2A9D8F"),
        (full, "Full data (leaky)", "#D55E00"),
    ):
        axis.plot(data["CutoffPPM"], -np.log10(data["PValue"]), color=color, linewidth=1.7, label=label)
        selected = data.loc[data["Selected"].astype(str).str.lower().eq("true")]
        axis.scatter(selected["CutoffPPM"], -np.log10(selected["PValue"]), s=45, color=color, edgecolor="white", zorder=3)
    axis.set_xscale("log")
    axis.set_xlabel("Candidate Faecalibacterium cutoff (PPM)")
    axis.set_ylabel(r"$-\log_{10}$(log-rank P)")
    axis.set_title("Outcome-guided cutoff search", loc="left", fontweight="bold")
    axis.legend(frameon=False)
    axis.grid(color="#ECECEC", linewidth=0.55)
    panel_label(axis, "A")

    heldout_audit = audit.loc[audit["EvaluationData"].eq("Held-out test")].iloc[0]
    heldout_curves = curves.loc[curves["EvaluationData"].eq("Held-out test")]
    plot_km(axes[1], heldout_curves, heldout_audit, "Training cutoff applied once", leaky=False)
    panel_label(axes[1], "B")

    full_audit = audit.loc[audit["EvaluationData"].eq("Full data")].iloc[0]
    full_curves = curves.loc[curves["EvaluationData"].eq("Full data")]
    plot_km(axes[2], full_curves, full_audit, "Same data select and evaluate", leaky=True)
    panel_label(axes[2], "C")
    fig.text(
        0.02,
        0.01,
        "The selected threshold happens to coincide here; the held-out P value still changes from 0.002 to 0.312.",
        fontsize=8,
        color="#4F5962",
    )
    fig.tight_layout(rect=[0, 0.045, 1, 1], w_pad=1.8)
    save(fig, output, "68-cutoff-leakage")


def time_dependent_roc(input_dir: Path, output: Path) -> None:
    roc = pd.read_csv(input_dir / "time-dependent-roc-365d.tsv", sep="\t")
    perf = pd.read_csv(input_dir / "cv-performance-summary.tsv", sep="\t")
    auc = perf.loc[perf["Metric"].eq("Cumulative/dynamic AUC") & ~perf["Model"].eq("Increment")].copy()
    delta = perf.loc[perf["Metric"].eq("Cumulative/dynamic AUC") & perf["Model"].eq("Increment")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4.4))
    axis = axes[0]
    for model in ("Clinical", "Clinical + Faecalibacterium"):
        data = roc.loc[roc["Model"].eq(model)].sort_values("FalsePositiveRate")
        auc_value = data["AUC"].iloc[0]
        axis.plot(data["FalsePositiveRate"], data["TruePositiveRate"], color=COLORS[model], linewidth=2.0, label=f"{model} · AUC={auc_value:.3f}")
    axis.plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=0.8)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_xlabel("1 − specificity at 12 months")
    axis.set_ylabel("Sensitivity at 12 months")
    axis.set_title("Cross-fitted time-dependent ROC", loc="left", fontweight="bold")
    axis.legend(frameon=False, loc="lower right")
    axis.grid(color="#ECECEC", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    offsets = {"Clinical": -8, "Clinical + Faecalibacterium": 8}
    for model in ("Clinical", "Clinical + Faecalibacterium"):
        data = auc.loc[auc["Model"].eq(model)].sort_values("HorizonDays")
        x = data["HorizonDays"] / 30.4375 + offsets[model] / 30.4375
        y = data["Estimate"]
        axis.errorbar(
            x,
            y,
            yerr=np.vstack([y - data["CILower"], data["CIUpper"] - y]),
            color=COLORS[model],
            marker="o",
            linewidth=1.8,
            capsize=3,
            label=model,
        )
    axis.axhline(0.5, color="#777777", linestyle="--", linewidth=0.8)
    for _, row in delta.iterrows():
        axis.text(row["HorizonDays"] / 30.4375, 0.49, f"Δ {row['Estimate']:+.03f}", ha="center", va="top", fontsize=7.5)
    axis.set_xticks(np.array([180, 365, 548]) / 30.4375, ["6", "12", "18"])
    axis.set_xlabel("Prediction horizon (months)")
    axis.set_ylabel("Cumulative/dynamic AUC (95% CI)")
    axis.set_ylim(0.42, 0.86)
    axis.set_title("Discrimination depends on time horizon", loc="left", fontweight="bold")
    axis.legend(frameon=False, loc="upper left")
    axis.grid(axis="y", color="#ECECEC", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "68-time-dependent-roc")


def crossvalidated_calibration(input_dir: Path, output: Path) -> None:
    calibration = pd.read_csv(input_dir / "calibration-365d.tsv", sep="\t")
    perf = pd.read_csv(input_dir / "cv-performance-summary.tsv", sep="\t")
    cindex = perf.loc[perf["Metric"].eq("Harrell C-index")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.35))
    axis = axes[0]
    axis.plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=0.8)
    for model in ("Clinical", "Clinical + Faecalibacterium"):
        data = calibration.loc[calibration["Model"].eq(model)].sort_values("MeanPredictedRisk")
        axis.errorbar(
            data["MeanPredictedRisk"],
            data["ObservedRiskKM"],
            yerr=np.vstack([data["ObservedRiskKM"] - data["ObservedRiskLower"], data["ObservedRiskUpper"] - data["ObservedRiskKM"]]),
            color=COLORS[model],
            marker="o",
            linewidth=1.7,
            capsize=3,
            label=model,
        )
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_xlabel("Mean cross-fitted predicted 12-month risk")
    axis.set_ylabel("Kaplan–Meier observed 12-month risk")
    axis.set_title("Calibration is not implied by AUC", loc="left", fontweight="bold")
    axis.legend(frameon=False, loc="upper left")
    axis.grid(color="#ECECEC", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    models = ["Clinical", "Clinical + Faecalibacterium"]
    data = cindex.set_index("Model").loc[models]
    y = np.arange(len(models))[::-1]
    for yy, model in zip(y, models, strict=True):
        row = data.loc[model]
        axis.errorbar(
            row["Estimate"],
            yy,
            xerr=[[row["Estimate"] - row["CILower"]], [row["CIUpper"] - row["Estimate"]]],
            fmt="o",
            color=COLORS[model],
            ecolor="#4F5962",
            capsize=4,
            linewidth=1.5,
            markersize=7,
        )
    axis.set_yticks(y, models)
    axis.set_xlim(0.50, 0.77)
    axis.set_xlabel("Cross-fitted Harrell C-index (95% CI)")
    increment = cindex.loc[cindex["Model"].eq("Increment")].iloc[0]
    axis.text(
        0.03,
        0.12,
        f"Increment {increment['Estimate']:+.3f}\n95% CI {increment['CILower']:+.3f} to {increment['CIUpper']:+.3f}",
        transform=axis.transAxes,
    )
    axis.set_title("Use paired out-of-fold predictions", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#ECECEC", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "68-crossvalidated-calibration")


def influence_sensitivity(input_dir: Path, output: Path) -> None:
    influence = pd.read_csv(input_dir / "leave-one-out-influence.tsv", sep="\t")
    estimates = pd.read_csv(input_dir / "cox-model-estimates.tsv", sep="\t")
    cohort = pd.read_csv(input_dir / "survival-cohort.tsv", sep="\t")
    primary = estimates.loc[
        estimates["Model"].eq("Adjusted primary") & estimates["Term"].eq("FaecalibacteriumLog2")
    ].iloc[0]
    influence = influence.sort_values("HazardRatio").reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.35), gridspec_kw={"width_ratios": [1.25, 1]})
    axis = axes[0]
    x = np.arange(len(influence))
    quality = influence["OmittedQualitySensitivityPass"].astype(str).str.lower().eq("true")
    axis.scatter(x[quality], influence.loc[quality, "HazardRatio"], s=22, color="#6F7D86", alpha=0.75, label="Omit one primary-analysis patient")
    axis.scatter(x[~quality], influence.loc[~quality, "HazardRatio"], s=55, color="#D55E00", edgecolor="white", linewidth=0.5, label="Omit sequencing-QC outlier")
    axis.axhline(primary["HazardRatio"], color="#2F6F91", linewidth=1.6, label="All patients")
    axis.axhspan(primary["CILower"], primary["CIUpper"], color="#2F6F91", alpha=0.10)
    axis.set_xlabel("Leave-one-patient-out refits (ordered by HR)")
    axis.set_ylabel("Adjusted HR per abundance doubling")
    axis.set_title("No single patient should define the conclusion", loc="left", fontweight="bold")
    axis.legend(frameon=False, loc="lower right")
    axis.grid(axis="y", color="#ECECEC", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    quality = cohort["QualitySensitivityPass"].astype(bool)
    axis.scatter(
        cohort.loc[quality, "PercentAssembled"],
        cohort.loc[quality, "FaecalibacteriumPPM"],
        s=26,
        color="#4C78A8",
        alpha=0.62,
        edgecolor="white",
        linewidth=0.4,
        label="Primary cohort",
    )
    axis.scatter(
        cohort.loc[~quality, "PercentAssembled"],
        cohort.loc[~quality, "FaecalibacteriumPPM"],
        s=65,
        color="#D55E00",
        edgecolor="white",
        linewidth=0.6,
        label="QC sensitivity exclusion",
    )
    axis.set_yscale("symlog", linthresh=25)
    axis.set_xlabel("Assembled non-human reads (%)")
    axis.set_ylabel("Faecalibacterium abundance (PPM)")
    axis.set_title("Keep the low-quality sample visible", loc="left", fontweight="bold")
    axis.legend(frameon=False, loc="upper left")
    axis.grid(color="#ECECEC", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "68-influence-sensitivity")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output = args.figure_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure()
    rng = np.random.default_rng(SEED)
    cohort_design(input_dir, output, rng)
    cox_forest(input_dir, output)
    assumption_and_spline(input_dir, output)
    cutoff_leakage(input_dir, output)
    time_dependent_roc(input_dir, output)
    crossvalidated_calibration(input_dir, output)
    influence_sensitivity(input_dir, output)
    shutil.copy2(input_dir / "spencer-fig3a-original.png", output / "68-spencer-fig3a-original.png")
    for stem in FIGURES:
        print(stem)


if __name__ == "__main__":
    main()
