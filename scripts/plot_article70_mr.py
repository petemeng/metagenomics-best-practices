#!/usr/bin/env python3
"""Create deterministic publication figures for Article 70."""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns


SEED = 20_260_770
COLORS = {
    "IVW": "#4C78A8",
    "Median": "#2A9D8F",
    "Egger": "#E45756",
    "Mode": "#B279A2",
    "PRESSO": "#F2A541",
    "Neutral": "#7C8A93",
    "Risk": "#D55E00",
}
FIGURES = (
    "70-study-design-assumptions",
    "70-instrument-harmonisation",
    "70-mr-method-comparison",
    "70-mr-scatter",
    "70-heterogeneity-pleiotropy",
    "70-leave-one-out",
    "70-outlier-sensitivity",
    "70-steiger-directionality",
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
    axis.text(-0.12, 1.06, label, transform=axis.transAxes, fontweight="bold", fontsize=11)


def add_box(
    axis: plt.Axes,
    center: tuple[float, float],
    text: str,
    facecolor: str,
    width: float,
    height: float,
    edgecolor: str = "#3C464D",
    linestyle: str = "-",
) -> None:
    x, y = center
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.018",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=1.0,
        linestyle=linestyle,
    )
    axis.add_patch(patch)
    axis.text(x, y, text, ha="center", va="center", fontsize=8.4, fontweight="bold")


def add_arrow(axis: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.35,
            color="#48545C",
        )
    )


def study_design_assumptions(input_dir: Path, output: Path) -> None:
    metrics = json.loads((input_dir / "analysis-metrics.json").read_text())
    model = json.loads((input_dir / "model-metrics.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.65), gridspec_kw={"width_ratios": [1.15, 1]})
    axis = axes[0]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    add_box(axis, (0.16, 0.70), "Exposure GWAS\nBMI · ieu-a-2\nN up to 339,224", "#DDEBF7", 0.25, 0.23)
    add_box(axis, (0.50, 0.70), "79 preselected\ngenome-wide hits\nP < 5e-8", "#FFF1C7", 0.24, 0.23)
    add_box(axis, (0.84, 0.70), "Outcome GWAS\nCHD · ieu-a-7\n60,801 / 123,504", "#F9DED8", 0.25, 0.23)
    add_arrow(axis, (0.29, 0.70), (0.37, 0.70))
    add_arrow(axis, (0.63, 0.70), (0.71, 0.70))
    axis.text(0.5, 0.48, "Same effect allele · summary associations only", ha="center", color="#45515A")
    add_box(axis, (0.25, 0.25), "No individual-level\nshotgun data", "#ECEFF1", 0.28, 0.18)
    add_box(axis, (0.75, 0.25), "Microbiome MR needs an\nexternal host mGWAS", "#E6F2EA", 0.32, 0.18)
    add_arrow(axis, (0.38, 0.25), (0.58, 0.25))
    axis.set_title("API-free two-sample MR teaching example", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    cards = [
        (0.78, "IV1 · Relevance", f"Supported\nminimum F = {model['minimum_f']:.1f}", "#D5EFE8"),
        (0.50, "IV2 · Exclusion", "Not provable\ncheck pleiotropy and colocalization", "#FFF1C7"),
        (0.22, "IV3 · Independence", "Not provable\naudit ancestry and confounder paths", "#F9DED8"),
    ]
    for y, title, detail, color in cards:
        add_box(axis, (0.5, y), f"{title}\n{detail}", color, 0.76, 0.19)
    axis.text(
        0.5,
        0.04,
        "A non-significant diagnostic does not verify an instrumental-variable assumption.",
        ha="center",
        fontsize=8,
        color="#4F5962",
    )
    axis.set_title("Three assumptions, only one directly measurable here", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "70-study-design-assumptions")


def instrument_harmonisation(input_dir: Path, output: Path) -> None:
    instruments = pd.read_csv(input_dir / "harmonised-instruments.tsv.gz", sep="\t")
    audit = pd.read_csv(input_dir / "harmonisation-audit.tsv", sep="\t").set_index("Quantity")
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.15), gridspec_kw={"width_ratios": [0.8, 1, 1.05]})
    axis = axes[0]
    stages = ["Exposure hits", "Outcome lookup", "Harmonised", "Retained"]
    values = [79, 79, 79, 79]
    colors = ["#BFD3E2", "#9EBFD1", "#6B9DB8", "#2F6B87"]
    bars = axis.barh(np.arange(4), values, color=colors, height=0.62)
    axis.set_yticks(np.arange(4), stages)
    axis.invert_yaxis()
    axis.set_xlim(0, 92)
    axis.set_xlabel("SNPs")
    for bar, value in zip(bars, values, strict=True):
        axis.text(value + 2, bar.get_y() + bar.get_height() / 2, str(value), va="center")
    axis.set_title("No SNP attrition", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    sns.histplot(instruments["FStatistic"], bins=22, color="#4C78A8", alpha=0.42, edgecolor="white", ax=axis)
    axis.axvline(10, color="#D55E00", linestyle="--", linewidth=1.2, label="Weak-instrument threshold")
    axis.axvline(instruments["FStatistic"].min(), color="#2A9D8F", linewidth=1.5, label=f"Minimum = {instruments['FStatistic'].min():.1f}")
    axis.set_xscale("log")
    axis.set_xlabel("Per-SNP F statistic (log scale)")
    axis.set_ylabel("Instruments")
    axis.set_title("Instrument strength", loc="left", fontweight="bold")
    axis.legend(frameon=False, fontsize=7.3)
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")

    axis = axes[2]
    ordered = instruments.sort_values("R2Exposure", ascending=False).reset_index(drop=True)
    axis.plot(np.arange(1, len(ordered) + 1), 100 * ordered["R2Exposure"].cumsum(), color="#2A9D8F", linewidth=2)
    axis.scatter(np.arange(1, len(ordered) + 1), 100 * ordered["R2Exposure"].cumsum(), s=10, color="#2A9D8F")
    total = 100 * ordered["R2Exposure"].sum()
    axis.axhline(total, color="#6F7D86", linestyle="--", linewidth=0.9)
    axis.text(78, total - 0.08, f"Approx. R² = {total:.2f}%", ha="right", va="top", fontsize=8)
    axis.set_xlabel("Instruments ordered by exposure R²")
    axis.set_ylabel("Cumulative BMI variance explained (%)")
    axis.set_title("Strength is not validity", loc="left", fontweight="bold")
    axis.grid(color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "C")
    fig.text(
        0.02,
        0.01,
        f"Palindromic retained: {int(audit.loc['Palindromic SNPs', 'Count'])} · ambiguous: {int(audit.loc['Ambiguous palindromic SNPs', 'Count'])} · outcome beta flips: {int(audit.loc['Outcome beta sign flips', 'Count'])}",
        fontsize=8,
        color="#4F5962",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1], w_pad=2.0)
    save(fig, output, "70-instrument-harmonisation")


def method_comparison(input_dir: Path, output: Path) -> None:
    estimates = pd.read_csv(input_dir / "mr-estimates.tsv", sep="\t")
    presso = pd.read_csv(input_dir / "mr-presso-estimates.tsv", sep="\t")
    method_order = [
        "Inverse variance weighted",
        "Weighted median",
        "MR Egger",
        "Simple mode",
        "Weighted mode",
    ]
    estimates = estimates.set_index("method").loc[method_order].reset_index()
    rows = pd.DataFrame(
        {
            "Method": estimates["method"],
            "OR": estimates["or"],
            "Lower": estimates["or_lci95"],
            "Upper": estimates["or_uci95"],
            "P": estimates["pval"],
            "Color": [COLORS["IVW"], COLORS["Median"], COLORS["Egger"], COLORS["Mode"], COLORS["Mode"]],
        }
    )
    corrected = presso.loc[presso["Analysis"].eq("Outlier-corrected")].iloc[0]
    rows.loc[len(rows)] = [
        "MR-PRESSO corrected",
        corrected["OddsRatio"],
        corrected["OddsRatioLower"],
        corrected["OddsRatioUpper"],
        corrected["PValue"],
        COLORS["PRESSO"],
    ]
    y = np.arange(len(rows))[::-1]
    fig, axis = plt.subplots(figsize=(8.5, 4.8))
    for yy, (_, row) in zip(y, rows.iterrows(), strict=True):
        axis.errorbar(
            row["OR"],
            yy,
            xerr=[[row["OR"] - row["Lower"]], [row["Upper"] - row["OR"]]],
            fmt="o",
            color=row["Color"],
            capsize=4,
            linewidth=1.6,
            markersize=6,
        )
        axis.text(2.42, yy, f"OR {row['OR']:.2f} ({row['Lower']:.2f}–{row['Upper']:.2f})", va="center", fontsize=8)
    axis.axvline(1, color="#333333", linestyle="--", linewidth=0.9)
    axis.set_yticks(y, rows["Method"])
    axis.set_xlim(0.92, 2.80)
    axis.set_xlabel("CHD odds ratio per genetically predicted 1-SD higher BMI")
    axis.set_title("Concordant direction across estimators with different validity assumptions", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    fig.text(
        0.02,
        0.01,
        "Confidence intervals use each estimator's native standard error; concordance does not prove the exclusion restriction.",
        fontsize=8,
        color="#4F5962",
    )
    fig.tight_layout(rect=[0, 0.045, 1, 1])
    save(fig, output, "70-mr-method-comparison")


def mr_scatter(input_dir: Path, output: Path) -> None:
    instruments = pd.read_csv(input_dir / "harmonised-instruments.tsv.gz", sep="\t")
    estimates = pd.read_csv(input_dir / "mr-estimates.tsv", sep="\t").set_index("method")
    egger = pd.read_csv(input_dir / "egger-intercept.tsv", sep="\t").iloc[0]
    outliers = pd.read_csv(input_dir / "mr-presso-outliers.tsv", sep="\t")
    outlier_snps = set(outliers.loc[outliers["DistortionOutlier"].astype(str).str.lower().eq("true"), "SNP"])
    instruments["Outlier"] = instruments["SNP"].isin(outlier_snps)
    weights = 1 / instruments["se.outcome"].pow(2)
    sizes = 16 + 72 * (weights - weights.min()) / (weights.max() - weights.min())
    fig, axis = plt.subplots(figsize=(7.25, 5.7))
    regular = instruments.loc[~instruments["Outlier"]]
    special = instruments.loc[instruments["Outlier"]]
    axis.scatter(
        regular["BetaExposureOriented"],
        regular["BetaOutcomeOriented"],
        s=sizes.loc[regular.index],
        color="#4C78A8",
        alpha=0.65,
        edgecolor="white",
        linewidth=0.5,
        label="Instrument SNP",
    )
    axis.scatter(
        special["BetaExposureOriented"],
        special["BetaOutcomeOriented"],
        s=sizes.loc[special.index] + 35,
        marker="D",
        color="#D55E00",
        edgecolor="white",
        linewidth=0.7,
        label="MR-PRESSO outlier",
        zorder=4,
    )
    for _, row in special.iterrows():
        axis.annotate(row["SNP"], (row["BetaExposureOriented"], row["BetaOutcomeOriented"]), xytext=(5, 5), textcoords="offset points", fontsize=7.5)
    x = np.linspace(0, instruments["BetaExposureOriented"].max() * 1.05, 100)
    ivw = estimates.loc["Inverse variance weighted", "b"]
    egger_beta = estimates.loc["MR Egger", "b"]
    axis.plot(x, ivw * x, color=COLORS["IVW"], linewidth=2, label="IVW")
    axis.plot(x, egger["egger_intercept"] + egger_beta * x, color=COLORS["Egger"], linewidth=2, linestyle="--", label="MR-Egger")
    axis.axhline(0, color="#777777", linewidth=0.7)
    axis.set_xlabel("SNP effect on BMI (oriented positive; SD units)")
    axis.set_ylabel("SNP effect on CHD (oriented log odds)")
    axis.set_title("Allele-oriented summary set", loc="left", fontweight="bold")
    axis.legend(frameon=False)
    axis.grid(color="#E8E8E8", linewidth=0.55)
    fig.text(
        0.02,
        0.01,
        "Point area reflects SNP–outcome precision; orientation preserves each Wald ratio.",
        fontsize=8,
        color="#4F5962",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    save(fig, output, "70-mr-scatter")


def _parse_pvalue(text: object) -> tuple[float, str]:
    label = str(text)
    if label.startswith("<"):
        return float(label[1:]), label
    value = float(label)
    return value, f"{value:.3g}"


def heterogeneity_pleiotropy(input_dir: Path, output: Path) -> None:
    single = pd.read_csv(input_dir / "single-snp-estimates.tsv.gz", sep="\t")
    single = single.loc[~single["SNP"].str.startswith("All -")].copy()
    estimates = pd.read_csv(input_dir / "mr-estimates.tsv", sep="\t").set_index("method")
    heterogeneity = pd.read_csv(input_dir / "mr-heterogeneity.tsv", sep="\t")
    egger = pd.read_csv(input_dir / "egger-intercept.tsv", sep="\t").iloc[0]
    presso = pd.read_csv(input_dir / "mr-presso-tests.tsv", sep="\t").set_index("Test")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.55))
    axis = axes[0]
    precision = 1 / single["se"]
    axis.scatter(single["b"], precision, s=25, color="#6B8DA3", alpha=0.67, edgecolor="white", linewidth=0.4)
    axis.axvline(estimates.loc["Inverse variance weighted", "b"], color=COLORS["IVW"], linewidth=1.8, label="IVW")
    axis.axvline(estimates.loc["MR Egger", "b"], color=COLORS["Egger"], linewidth=1.8, linestyle="--", label="MR-Egger")
    axis.set_xlabel("Single-SNP Wald ratio (log odds / SD BMI)")
    axis.set_ylabel("Precision (1 / SE)")
    axis.set_title("Funnel asymmetry is a diagnostic, not a test alone", loc="left", fontweight="bold")
    axis.legend(frameon=False)
    axis.grid(color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    ivw_q = heterogeneity.loc[heterogeneity["method"].eq("Inverse variance weighted")].iloc[0]
    egger_q = heterogeneity.loc[heterogeneity["method"].eq("MR Egger")].iloc[0]
    global_p, global_label = _parse_pvalue(presso.loc["Global", "PValueText"])
    distortion_p, distortion_label = _parse_pvalue(presso.loc["Distortion", "PValueText"])
    tests = pd.DataFrame(
        {
            "Test": ["IVW heterogeneity", "Egger heterogeneity", "Egger intercept", "MR-PRESSO global", "PRESSO distortion"],
            "P": [ivw_q["Q_pval"], egger_q["Q_pval"], egger["pval"], global_p, distortion_p],
            "Label": [f"{ivw_q['Q_pval']:.2g}", f"{egger_q['Q_pval']:.2g}", f"{egger['pval']:.3f}", global_label, distortion_label],
        }
    )
    y = np.arange(len(tests))[::-1]
    colors = [COLORS["Risk"], COLORS["Risk"], COLORS["Neutral"], COLORS["Risk"], COLORS["Neutral"]]
    values = -np.log10(tests["P"].clip(lower=1e-300))
    axis = axes[1]
    axis.barh(y, values, color=colors, alpha=0.82, height=0.62)
    axis.axvline(-math.log10(0.05), color="#333333", linestyle="--", linewidth=0.9)
    axis.set_yticks(y, tests["Test"])
    axis.set_xlabel("−log10(P)")
    axis.set_title("Different diagnostics answer different questions", loc="left", fontweight="bold")
    for yy, value, label in zip(y, values, tests["Label"], strict=True):
        axis.text(value + 0.08, yy, f"P={label}", va="center", fontsize=7.7)
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "70-heterogeneity-pleiotropy")


def leave_one_out_plot(input_dir: Path, output: Path) -> None:
    leave = pd.read_csv(input_dir / "leave-one-out.tsv", sep="\t")
    leave = leave.loc[leave["SNP"].ne("All")].sort_values("b").reset_index(drop=True)
    estimates = pd.read_csv(input_dir / "mr-estimates.tsv", sep="\t").set_index("method")
    ivw = estimates.loc["Inverse variance weighted", "b"]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.45), gridspec_kw={"width_ratios": [1.35, 0.85]})
    axis = axes[0]
    axis.fill_between(np.arange(len(leave)), leave["CILower"], leave["CIUpper"], color="#BFD3E2", alpha=0.38, linewidth=0)
    axis.plot(np.arange(len(leave)), leave["b"], color="#4C78A8", linewidth=1.5)
    axis.scatter(np.arange(len(leave)), leave["b"], s=18, color="#4C78A8")
    axis.axhline(ivw, color="#333333", linestyle="--", linewidth=0.9, label="All-instrument IVW")
    axis.axhline(0, color="#888888", linewidth=0.7)
    axis.set_xlabel("Omitted SNP analyses, ordered by estimate")
    axis.set_ylabel("IVW estimate (log odds / SD BMI)")
    axis.set_title("No single omission reverses the estimate", loc="left", fontweight="bold")
    axis.legend(frameon=False)
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    extremes = pd.concat([leave.head(4), leave.tail(4)]).sort_values("b")
    axis = axes[1]
    y = np.arange(len(extremes))[::-1]
    axis.errorbar(
        extremes["b"],
        y,
        xerr=np.vstack([extremes["b"] - extremes["CILower"], extremes["CIUpper"] - extremes["b"]]),
        fmt="o",
        color="#2A9D8F",
        capsize=3,
    )
    axis.axvline(ivw, color="#333333", linestyle="--", linewidth=0.9)
    axis.set_yticks(y, [f"Omit {snp}" for snp in extremes["SNP"]])
    axis.set_xlabel("IVW estimate")
    axis.set_title("Most influential omissions", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "70-leave-one-out")


def outlier_sensitivity(input_dir: Path, output: Path) -> None:
    main = pd.read_csv(input_dir / "mr-estimates.tsv", sep="\t")
    corrected = pd.read_csv(input_dir / "presso-outlier-exclusion-estimates.tsv", sep="\t")
    radial = pd.read_csv(input_dir / "radial-ivw-outliers.tsv", sep="\t").sort_values("Q_statistic", ascending=False)
    presso = pd.read_csv(input_dir / "mr-presso-outliers.tsv", sep="\t")
    outlier_snps = set(presso.loc[presso["DistortionOutlier"].astype(str).str.lower().eq("true"), "SNP"])
    selected_methods = ["Inverse variance weighted", "Weighted median", "MR Egger"]
    records = []
    for label, frame in (("All 79 SNPs", main), ("Exclude 2 PRESSO outliers", corrected)):
        subset = frame.set_index("method").loc[selected_methods]
        for method, row in subset.iterrows():
            records.append(
                {
                    "Specification": label,
                    "Method": method,
                    "Estimate": row["b"],
                    "Lower": row["lo_ci"],
                    "Upper": row["up_ci"],
                }
            )
    records = pd.DataFrame(records)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.55))
    axis = axes[0]
    method_colors = {"Inverse variance weighted": COLORS["IVW"], "Weighted median": COLORS["Median"], "MR Egger": COLORS["Egger"]}
    positions = {"All 79 SNPs": 0.08, "Exclude 2 PRESSO outliers": -0.08}
    base_y = {method: 2 - index for index, method in enumerate(selected_methods)}
    for _, row in records.iterrows():
        yy = base_y[row["Method"]] + positions[row["Specification"]]
        marker = "o" if row["Specification"] == "All 79 SNPs" else "D"
        axis.errorbar(
            row["Estimate"],
            yy,
            xerr=[[row["Estimate"] - row["Lower"]], [row["Upper"] - row["Estimate"]]],
            fmt=marker,
            color=method_colors[row["Method"]],
            capsize=3,
            alpha=0.95,
        )
    axis.axvline(0, color="#333333", linestyle="--", linewidth=0.8)
    axis.set_yticks([2, 1, 0], selected_methods)
    axis.set_xlabel("MR estimate (log CHD odds / SD BMI)")
    axis.set_title("Direction survives outlier exclusion", loc="left", fontweight="bold")
    axis.scatter([], [], marker="o", color="#555555", label="All 79 SNPs")
    axis.scatter([], [], marker="D", color="#555555", label="Exclude 2 PRESSO outliers")
    axis.legend(frameon=False, loc="lower right")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    top = radial.head(12).sort_values("Q_statistic")
    colors = [COLORS["Risk"] if snp in outlier_snps else "#7C8A93" for snp in top["SNP"]]
    axis.barh(np.arange(len(top)), top["Q_statistic"], color=colors, height=0.65)
    axis.set_yticks(np.arange(len(top)), top["SNP"])
    axis.set_xlabel("Radial IVW contribution to Q")
    axis.set_title("PRESSO outliers are also heterogeneity leaders", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "70-outlier-sensitivity")


def steiger_directionality(input_dir: Path, output: Path) -> None:
    aggregate = pd.read_csv(input_dir / "steiger-directionality.tsv", sep="\t")
    per_snp = pd.read_csv(input_dir / "steiger-per-snp.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.7, 4.45))
    axis = axes[0]
    axis.plot(aggregate["AssumedCHDPrevalence"] * 100, aggregate["ExposureR2"] * 100, marker="o", color="#4C78A8", linewidth=2, label="BMI exposure")
    axis.plot(aggregate["AssumedCHDPrevalence"] * 100, aggregate["OutcomeLiabilityR2"] * 100, marker="o", color="#E45756", linewidth=2, label="CHD liability")
    axis.set_xlabel("Assumed population CHD prevalence (%)")
    axis.set_ylabel("Sum of SNP R² (%)")
    axis.set_title("Aggregate Steiger direction is prevalence-robust", loc="left", fontweight="bold")
    axis.legend(frameon=False)
    axis.grid(color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    forward = per_snp["ExposureExplainsMore"].astype(str).str.lower().eq("true")
    axis.scatter(
        per_snp.loc[forward, "ExposureR2"],
        per_snp.loc[forward, "OutcomeLiabilityR2"],
        s=25,
        color="#2A9D8F",
        alpha=0.7,
        label=f"Exposure R² larger (n={int(forward.sum())})",
    )
    axis.scatter(
        per_snp.loc[~forward, "ExposureR2"],
        per_snp.loc[~forward, "OutcomeLiabilityR2"],
        s=29,
        color="#D55E00",
        alpha=0.8,
        label=f"Outcome R² larger (n={int((~forward).sum())})",
    )
    lower = min(per_snp[["ExposureR2", "OutcomeLiabilityR2"]].replace(0, np.nan).min()) * 0.7
    upper = max(per_snp[["ExposureR2", "OutcomeLiabilityR2"]].max()) * 1.4
    axis.plot([lower, upper], [lower, upper], color="#555555", linestyle="--", linewidth=0.9)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(lower, upper)
    axis.set_ylim(lower, upper)
    axis.set_xlabel("Per-SNP BMI R²")
    axis.set_ylabel("Per-SNP CHD liability R²")
    axis.set_title("Per-SNP direction is not unanimous", loc="left", fontweight="bold")
    axis.legend(frameon=False, fontsize=7.4)
    axis.grid(color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "70-steiger-directionality")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output = args.figure_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure()
    np.random.seed(SEED)
    study_design_assumptions(input_dir, output)
    instrument_harmonisation(input_dir, output)
    method_comparison(input_dir, output)
    mr_scatter(input_dir, output)
    heterogeneity_pleiotropy(input_dir, output)
    leave_one_out_plot(input_dir, output)
    outlier_sensitivity(input_dir, output)
    steiger_directionality(input_dir, output)
    shutil.copy2(input_dir / "hemani-figure1-original.png", output / "70-hemani-figure1-original.png")
    for stem in FIGURES:
        print(stem)


if __name__ == "__main__":
    main()
