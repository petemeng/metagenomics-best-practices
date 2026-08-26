#!/usr/bin/env python3
"""Create deterministic publication figures for Article 69."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.proportion import proportion_confint


SEED = 20_260_769
COLORS = {
    "Insufficient": "#9AA5B1",
    "Sufficient": "#2A9D8F",
    "Direct": "#4C78A8",
    "Indirect": "#E45756",
    "Total": "#B279A2",
    "Before weighting": "#D55E00",
    "After IPTW": "#4C78A8",
}
FIGURES = (
    "69-cohort-observed-path",
    "69-dag-identification",
    "69-path-models",
    "69-effect-decomposition",
    "69-bootstrap-distributions",
    "69-positivity-balance",
    "69-sensitivity-analysis",
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
    fig.savefig(output / f"{stem}.tiff", dpi=360, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.13, 1.06, label, transform=axis.transAxes, fontweight="bold", fontsize=11)


def cohort_observed_path(input_dir: Path, output: Path, rng: np.random.Generator) -> None:
    attrition = pd.read_csv(input_dir / "cohort-attrition.tsv", sep="\t")
    cohort = pd.read_csv(input_dir / "mediation-cohort.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.25), gridspec_kw={"width_ratios": [1.2, 0.9, 1.05]})
    axis = axes[0]
    colors = ["#C2CDD5", "#A6BBC7", "#78A3B7", "#4F829C", "#2F6179"]
    bars = axis.barh(np.arange(len(attrition)), attrition["Samples"], color=colors, height=0.68)
    axis.set_yticks(np.arange(len(attrition)), attrition["Stage"])
    axis.invert_yaxis()
    axis.set_xlim(0, 184)
    axis.set_xlabel("Patients retained")
    axis.set_title("Real WGS–fiber–response intersection", loc="left", fontweight="bold")
    for bar, n, r in zip(bars, attrition["Samples"], attrition["Responders"], strict=True):
        axis.text(n + 3, bar.get_y() + bar.get_height() / 2, f"{n} / {r} responders", va="center", fontsize=7.5)
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    records = []
    for exposure, label in ((0, "Insufficient"), (1, "Sufficient")):
        group = cohort.loc[cohort["ExposureSufficient"].eq(exposure)]
        count = int(group["OutcomeResponder"].sum())
        n = len(group)
        lower, upper = proportion_confint(count, n, method="wilson")
        records.append((label, n, count / n, lower, upper))
    for position, (label, n, rate, lower, upper) in enumerate(records):
        axis.bar(position, rate * 100, color=COLORS[label], width=0.62, alpha=0.9)
        axis.errorbar(position, rate * 100, yerr=[[rate * 100 - lower * 100], [upper * 100 - rate * 100]], color="#222222", capsize=4, linewidth=1.1)
        axis.text(position, rate * 100 + 9, f"{rate * 100:.1f}%\nn={n}", ha="center", fontsize=8)
    axis.set_xticks([0, 1], ["Insufficient\n<20 g/day", "Sufficient\n≥20 g/day"])
    axis.set_ylim(0, 110)
    axis.set_ylabel("Observed ICB response (%)")
    axis.set_title("Exposure–outcome association", loc="left", fontweight="bold")
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")

    axis = axes[2]
    for position, (exposure, label) in enumerate(((0, "Insufficient"), (1, "Sufficient"))):
        values = cohort.loc[cohort["ExposureSufficient"].eq(exposure), "FaecalibacteriumLog2"].to_numpy(float)
        violin = axis.violinplot(values, positions=[position], widths=0.72, showextrema=False)
        for body in violin["bodies"]:
            body.set_facecolor(COLORS[label])
            body.set_edgecolor("none")
            body.set_alpha(0.28)
        axis.scatter(position + rng.uniform(-0.13, 0.13, len(values)), values, s=22, color=COLORS[label], alpha=0.72, edgecolor="white", linewidth=0.4)
        axis.plot([position - 0.22, position + 0.22], [np.median(values)] * 2, color="#222222", linewidth=1.5)
    axis.set_xticks([0, 1], ["Insufficient", "Sufficient"])
    axis.set_ylabel(r"log$_2$(Faecalibacterium PPM + 25)")
    axis.set_title("Exposure–mediator association", loc="left", fontweight="bold")
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "C")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "69-cohort-observed-path")


def add_box(axis: plt.Axes, xy: tuple[float, float], text: str, color: str, width: float = 0.20, height: float = 0.16, linestyle: str = "-") -> None:
    x, y = xy
    box = FancyBboxPatch((x - width / 2, y - height / 2), width, height, boxstyle="round,pad=0.02", facecolor=color, edgecolor="#30363B", linewidth=1.0, linestyle=linestyle)
    axis.add_patch(box)
    axis.text(x, y, text, ha="center", va="center", fontsize=8.5, fontweight="bold")


def add_arrow(axis: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = "#38434B", linestyle: str = "-", label: str | None = None, curve: float = 0.0) -> None:
    arrow = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12, linewidth=1.35, color=color, linestyle=linestyle, connectionstyle=f"arc3,rad={curve}")
    axis.add_patch(arrow)
    if label:
        middle = ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
        axis.text(middle[0], middle[1] + 0.035, label, color=color, ha="center", fontsize=8)


def dag_identification(input_dir: Path, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), gridspec_kw={"width_ratios": [1.25, 1]})
    axis = axes[0]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    add_box(axis, (0.16, 0.55), "A\nFiber sufficiency", "#D5EFE8", width=0.22)
    add_box(axis, (0.50, 0.55), "M\nFaecalibacterium", "#F9DED8", width=0.23)
    add_box(axis, (0.84, 0.55), "Y\nICB response", "#DFE8F2", width=0.20)
    add_box(axis, (0.50, 0.86), "C\nMeasured baseline\nconfounders", "#ECEFF1", width=0.25, height=0.18)
    add_box(axis, (0.50, 0.18), "U\nUnmeasured diet, health,\ndisease and behavior", "#FFF1C7", width=0.31, height=0.19, linestyle="--")
    add_arrow(axis, (0.28, 0.55), (0.38, 0.55), label="a path")
    add_arrow(axis, (0.62, 0.55), (0.74, 0.55), label="b path")
    add_arrow(axis, (0.25, 0.49), (0.75, 0.49), color="#4C78A8", curve=0.20, label="direct path")
    add_arrow(axis, (0.46, 0.79), (0.22, 0.64), color="#6B747B")
    add_arrow(axis, (0.50, 0.77), (0.50, 0.65), color="#6B747B")
    add_arrow(axis, (0.54, 0.79), (0.78, 0.64), color="#6B747B")
    add_arrow(axis, (0.45, 0.27), (0.43, 0.45), color="#B27900", linestyle="--")
    add_arrow(axis, (0.56, 0.27), (0.76, 0.48), color="#B27900", linestyle="--")
    axis.text(0.03, 0.03, "Identification requires no unmeasured A–M, A–Y, or M–Y confounding\nand no exposure-induced mediator–outcome confounder.", fontsize=8, color="#4F5962")
    axis.set_title("The DAG states assumptions; it does not verify them", loc="left", fontweight="bold")
    panel_label(axis, "A")

    axis = axes[1]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    axis.plot([0.12, 0.88], [0.52, 0.52], color="#505A61", linewidth=2)
    positions = [0.18, 0.48, 0.82]
    labels = [
        ("Habitual fiber\nquestionnaire", "Exposure", "#D5EFE8"),
        ("Baseline stool\nshotgun WGS", "Mediator", "#F9DED8"),
        ("RECIST-based\nICB response", "Outcome", "#DFE8F2"),
    ]
    for position, (title, role, color) in zip(positions, labels, strict=True):
        axis.scatter(position, 0.52, s=80, color="#30363B", zorder=3)
        box = FancyBboxPatch((position - 0.12, 0.61), 0.24, 0.19, boxstyle="round,pad=0.02", facecolor=color, edgecolor="#30363B", linewidth=1)
        axis.add_patch(box)
        axis.text(position, 0.70, title, ha="center", va="center", fontsize=8.5, fontweight="bold")
        axis.text(position, 0.42, role, ha="center", color="#4F5962", fontsize=8)
    axis.annotate("Same baseline window", xy=(0.33, 0.55), xytext=(0.33, 0.26), ha="center", arrowprops={"arrowstyle": "-[", "color": "#D55E00", "lw": 1.2}, color="#D55E00", fontsize=8)
    axis.text(0.5, 0.08, "Exposure-before-mediator ordering is plausible but not experimentally established.", ha="center", fontsize=8, color="#4F5962")
    axis.set_title("Timing is the weakest link in this public example", loc="left", fontweight="bold")
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.0)
    save(fig, output, "69-dag-identification")


def path_models(input_dir: Path, output: Path) -> None:
    path = pd.read_csv(input_dir / "path-model-estimates.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.25))
    axis = axes[0]
    mediator = path.loc[path["Model"].eq("Mediator model") & path["Term"].eq("A")].iloc[0]
    continuous = path.loc[path["Model"].eq("Continuous-fiber mediator model") & path["Term"].eq("I(FiberGrams/5)")].iloc[0]
    rows = [mediator, continuous]
    labels = ["Sufficient vs insufficient", "Per 5 g/day fiber"]
    y = np.array([1, 0])
    for yy, row, color in zip(y, rows, ["#2A9D8F", "#4C78A8"], strict=True):
        axis.errorbar(row["Estimate"], yy, xerr=[[row["Estimate"] - row["CILower"]], [row["CIUpper"] - row["Estimate"]]], fmt="o", color=color, capsize=4, linewidth=1.5)
        axis.text(2.55, yy, f"P={row['PValue']:.3f}", va="center", fontsize=8)
    axis.axvline(0, color="#333333", linestyle="--", linewidth=0.8)
    axis.set_yticks(y, labels)
    axis.set_xlabel(r"Adjusted difference in log$_2$(Faecalibacterium PPM + 25)")
    axis.set_xlim(-0.8, 2.9)
    axis.set_title("Exposure → mediator model", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    total = path.loc[path["Model"].eq("Total-association model") & path["Term"].eq("A")].iloc[0]
    direct = path.loc[path["Model"].eq("Outcome model") & path["Term"].eq("A")].iloc[0]
    mediator_y = path.loc[path["Model"].eq("Outcome model") & path["Term"].eq("M")].iloc[0]
    rows = [total, direct, mediator_y]
    labels = ["Fiber total-association OR", "Fiber direct-model OR", "Mediator OR per doubling"]
    colors = [COLORS["Total"], COLORS["Direct"], COLORS["Indirect"]]
    y = np.arange(3)[::-1]
    for yy, row, color in zip(y, rows, colors, strict=True):
        axis.errorbar(row["OddsRatio"], yy, xerr=[[row["OddsRatio"] - row["OddsRatioLower"]], [row["OddsRatioUpper"] - row["OddsRatio"]]], fmt="o", color=color, capsize=4, linewidth=1.5)
        axis.text(13.0, yy, f"P={row['PValue']:.3f}", va="center", fontsize=8)
    axis.axvline(1, color="#333333", linestyle="--", linewidth=0.8)
    axis.set_xscale("log")
    axis.set_xlim(0.7, 15.8)
    axis.set_yticks(y, labels)
    axis.set_xlabel("Adjusted odds ratio (95% CI)")
    axis.set_title("Mediator and outcome models", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.4)
    save(fig, output, "69-path-models")


def effect_decomposition(input_dir: Path, output: Path) -> None:
    summary = pd.read_csv(input_dir / "gformula-effect-summary.tsv", sep="\t")
    primary = summary.loc[summary["Variant"].eq("Primary")].set_index("Effect")
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.45), gridspec_kw={"width_ratios": [1.2, 1]})
    axis = axes[0]
    scenarios = ["P00", "P10", "P11"]
    labels = ["Low fiber + M(low)", "Sufficient fiber + M(low)", "Sufficient fiber + M(high)"]
    x = np.arange(3)
    data = primary.loc[scenarios]
    axis.errorbar(
        x,
        data["Estimate"] * 100,
        yerr=np.vstack([(data["Estimate"] - data["CILower"]) * 100, (data["CIUpper"] - data["Estimate"]) * 100]),
        fmt="o-",
        color="#384E5C",
        ecolor="#7C8A93",
        capsize=4,
        linewidth=1.8,
        markersize=6,
    )
    axis.annotate("Direct +20.6 pp", xy=(0.5, 69), ha="center", color=COLORS["Direct"], fontsize=8.5, fontweight="bold")
    axis.annotate("Indirect +1.0 pp", xy=(1.5, 81.5), ha="center", color=COLORS["Indirect"], fontsize=8.5, fontweight="bold")
    axis.set_xticks(x, labels, rotation=12, ha="right")
    axis.set_ylabel("Standardized response probability (%)")
    axis.set_ylim(38, 103)
    axis.set_title("G-computation counterfactual scenarios", loc="left", fontweight="bold")
    axis.grid(axis="y", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    effects = ["Direct", "Indirect", "Total"]
    data = primary.loc[effects]
    y = np.arange(3)[::-1]
    for yy, effect in zip(y, effects, strict=True):
        row = data.loc[effect]
        axis.errorbar(row["Estimate"] * 100, yy, xerr=[[(row["Estimate"] - row["CILower"]) * 100], [(row["CIUpper"] - row["Estimate"]) * 100]], fmt="o", color=COLORS[effect], capsize=4, linewidth=1.6, markersize=6)
        axis.text(42, yy, f"{row['Estimate']*100:+.1f} pp\nP={row['BootstrapP']:.3f}", va="center", fontsize=8)
    axis.axvline(0, color="#333333", linestyle="--", linewidth=0.8)
    axis.set_yticks(y, ["Direct", "Indirect through mediator", "Total"])
    axis.set_xlabel("Risk difference (percentage points; 95% bootstrap CI)")
    axis.set_xlim(-8, 49)
    axis.set_title("Model-based decomposition", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "69-effect-decomposition")


def bootstrap_distributions(input_dir: Path, output: Path) -> None:
    bootstrap = pd.read_csv(input_dir / "primary-gformula-bootstrap.tsv.gz", sep="\t")
    summary = pd.read_csv(input_dir / "gformula-effect-summary.tsv", sep="\t")
    primary = summary.loc[summary["Variant"].eq("Primary")].set_index("Effect")
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8), sharex=True)
    for index, effect in enumerate(("Direct", "Indirect", "Total")):
        axis = axes[index]
        values = bootstrap[effect].dropna().to_numpy(float) * 100
        sns.histplot(values, bins=38, stat="density", color=COLORS[effect], alpha=0.35, edgecolor="white", ax=axis)
        sns.kdeplot(values, color=COLORS[effect], linewidth=1.8, ax=axis)
        row = primary.loc[effect]
        axis.axvline(0, color="#333333", linestyle="--", linewidth=0.8)
        axis.axvline(row["Estimate"] * 100, color=COLORS[effect], linewidth=2.0)
        axis.axvspan(row["CILower"] * 100, row["CIUpper"] * 100, color=COLORS[effect], alpha=0.12)
        axis.set_title(f"{effect} effect", loc="left", fontweight="bold")
        axis.set_xlabel("Risk difference (percentage points)")
        axis.set_ylabel("Bootstrap density" if index == 0 else "")
        axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
        panel_label(axis, chr(ord("A") + index))
    fig.text(0.02, 0.01, "5,000 patient-level bootstrap refits · vertical color line: point estimate · shaded span: percentile 95% CI", fontsize=8, color="#4F5962")
    fig.tight_layout(rect=[0, 0.045, 1, 1], w_pad=1.5)
    save(fig, output, "69-bootstrap-distributions")


def positivity_balance(input_dir: Path, output: Path) -> None:
    overlap = pd.read_csv(input_dir / "exposure-overlap.tsv", sep="\t")
    balance = pd.read_csv(input_dir / "exposure-balance.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.25))
    axis = axes[0]
    for exposure, label in ((0, "Insufficient"), (1, "Sufficient")):
        values = overlap.loc[overlap["A"].eq(exposure), "Propensity"]
        sns.kdeplot(values, fill=True, alpha=0.22, linewidth=1.8, color=COLORS[label], label=f"{label} (n={len(values)})", ax=axis, clip=(0, 1))
    axis.set_xlim(0, 0.55)
    axis.set_xlabel("Estimated probability of sufficient fiber")
    axis.set_ylabel("Density")
    axis.set_title("Positivity / overlap audit", loc="left", fontweight="bold")
    axis.text(0.03, 0.93, f"Maximum IPTW = {overlap['IPTW'].max():.2f}", transform=axis.transAxes, va="top")
    axis.legend(frameon=False)
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    labels = {"BMIz": "BMI", "Mucosal": "Mucosal/acral subtype", "StageM1D": "Stage M1D", "LDHHigh": "LDH above upper limit"}
    y = np.arange(len(balance))[::-1]
    axis.scatter(balance["SMDUnweighted"].abs(), y + 0.08, color=COLORS["Before weighting"], s=40, label="Before weighting")
    axis.scatter(balance["SMDIPTW"].abs(), y - 0.08, color=COLORS["After IPTW"], s=40, label="After IPTW diagnostic")
    for yy, before, after in zip(y, balance["SMDUnweighted"].abs(), balance["SMDIPTW"].abs(), strict=True):
        axis.plot([before, after], [yy + 0.08, yy - 0.08], color="#AAB2B8", linewidth=0.8)
    axis.axvline(0.1, color="#333333", linestyle="--", linewidth=0.8)
    axis.set_yticks(y, [labels[item] for item in balance["Covariate"]])
    axis.set_xlabel("Absolute standardized mean difference")
    axis.set_xlim(-0.01, 0.34)
    axis.set_title("Measured balance does not prove exchangeability", loc="left", fontweight="bold")
    axis.legend(frameon=False, loc="lower right")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "69-positivity-balance")


def sensitivity_analysis(input_dir: Path, output: Path) -> None:
    summary = pd.read_csv(input_dir / "gformula-effect-summary.tsv", sep="\t")
    rho = pd.read_csv(input_dir / "residual-correlation-sensitivity.tsv", sep="\t")
    indirect = summary.loc[summary["Effect"].eq("Indirect")].copy()
    order = ["Primary", "Sequencing-QC subset", "Ruminococcaceae mediator", "Exposure-mediator interaction"]
    indirect = indirect.set_index("Variant").loc[order]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.35), gridspec_kw={"width_ratios": [1.05, 1.15]})
    axis = axes[0]
    y = np.arange(len(order))[::-1]
    for yy, variant in zip(y, order, strict=True):
        row = indirect.loc[variant]
        color = "#E45756" if variant == "Primary" else "#6F7D86"
        axis.errorbar(row["Estimate"] * 100, yy, xerr=[[(row["Estimate"] - row["CILower"]) * 100], [(row["CIUpper"] - row["Estimate"]) * 100]], fmt="o", color=color, capsize=4, linewidth=1.5)
        axis.text(22, yy, f"{row['Estimate']*100:+.1f} pp", va="center", fontsize=8)
    axis.axvline(0, color="#333333", linestyle="--", linewidth=0.8)
    axis.set_yticks(y, order)
    axis.set_xlim(-5, 25)
    axis.set_xlabel("Indirect risk difference (95% bootstrap CI)")
    axis.set_title("Mediator definition and interaction matter", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "A")

    axis = axes[1]
    axis.fill_between(rho["Rho"], rho["LowerAverage"] * 100, rho["UpperAverage"] * 100, color="#B279A2", alpha=0.18, linewidth=0)
    axis.plot(rho["Rho"], rho["ACMEAverage"] * 100, color="#B279A2", linewidth=2.0)
    axis.axhline(0, color="#333333", linestyle="--", linewidth=0.8)
    axis.axvline(0, color="#777777", linestyle=":", linewidth=0.8)
    axis.scatter([0.05], [0], color="#D55E00", s=38, zorder=3)
    axis.text(0.07, 0.8, "Zero crossing near ρ=0.05", color="#D55E00", fontsize=8)
    axis.set_xlabel("Assumed mediator/outcome residual correlation (ρ)")
    axis.set_ylabel("Average indirect effect (percentage points)")
    axis.set_title("Sequential ignorability is highly fragile", loc="left", fontweight="bold")
    axis.grid(color="#E8E8E8", linewidth=0.55)
    panel_label(axis, "B")
    fig.tight_layout(w_pad=2.2)
    save(fig, output, "69-sensitivity-analysis")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output = args.figure_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure()
    rng = np.random.default_rng(SEED)
    cohort_observed_path(input_dir, output, rng)
    dag_identification(input_dir, output)
    path_models(input_dir, output)
    effect_decomposition(input_dir, output)
    bootstrap_distributions(input_dir, output)
    positivity_balance(input_dir, output)
    sensitivity_analysis(input_dir, output)
    shutil.copy2(input_dir / "spencer-fig3ab-original.png", output / "69-spencer-fig3ab-original.png")
    for stem in FIGURES:
        print(stem)


if __name__ == "__main__":
    main()
