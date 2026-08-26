#!/usr/bin/env python3
"""Render publication figures for Article 74 from checksum-frozen tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


PLOT_SEED = 20_260_774
PALETTE = ["#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", "#264653", "#6C5CE7"]


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
            "axes.titlesize": 13,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.titleweight": "bold",
            "svg.fonttype": "none",
        }
    )


def save(fig: plt.Figure, output: Path, name: str) -> list[Path]:
    output.mkdir(parents=True, exist_ok=True)
    png = output / f"{name}.png"
    svg = output / f"{name}.svg"
    fig.savefig(png, dpi=360, bbox_inches="tight", pad_inches=0.12)
    fig.savefig(svg, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    return [png, svg]


def rounded_box(ax: plt.Axes, xy: tuple[float, float], width: float, height: float, color: str, text: str, *, text_color: str = "white") -> None:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.04",
        linewidth=0,
        facecolor=color,
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", color=text_color, weight="bold", fontsize=9)


def release_lock(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "release-lock.tsv", sep="\t")
    dates = pd.to_datetime(frame["PublishedUTC"], utc=True)
    origin = pd.Timestamp("2026-04-01", tz="UTC")
    days = (dates - origin).dt.total_seconds() / 86400
    fig, ax = plt.subplots(figsize=(9.0, 4.2))
    y = np.arange(len(frame))[::-1]
    ax.hlines(y, 0, days, color="#DCE4E8", linewidth=4, zorder=1)
    for idx, row in frame.iterrows():
        yi = y[idx]
        ax.scatter(days.iloc[idx], yi, s=180, color=PALETTE[idx], edgecolor="white", linewidth=1.5, zorder=3)
        commit = row["Commit"][:10] if row["Commit"] != "release asset" else "signed asset"
        ax.text(days.iloc[idx] + 3, yi, f"{row['Component']} {row['Release']}\n{commit}", va="center", fontsize=9)
    tick_dates = pd.to_datetime(["2026-04-01", "2026-05-01", "2026-06-01", "2026-07-01", "2026-08-01"], utc=True)
    ticks = (tick_dates - origin).total_seconds() / 86400
    ax.set_xticks(ticks, [date.strftime("%b %Y") for date in tick_dates])
    ax.set_yticks([])
    ax.set_xlim(-3, max(days) + 45)
    ax.set_title("Release locks turn a moving workflow into an auditable object", loc="left")
    ax.set_xlabel("Official release date (UTC)")
    ax.grid(axis="x", color="#EDF1F2", linewidth=0.8)
    fig.text(0.01, 0.01, "Tag + full commit + artifact SHA256 are recorded together", fontsize=8, color="#51636B")
    return save(fig, output, "74-release-lock")


def parameter_precedence(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "parameter-precedence.tsv", sep="\t").sort_values("Rank")
    fig, ax = plt.subplots(figsize=(8.3, 5.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.3)
    widths = [5.4, 5.7, 6.0, 6.3]
    for idx, row in frame.iterrows():
        y = 0.55 + idx * 1.05
        rounded_box(ax, (0.5, y), widths[idx], 0.72, PALETTE[idx], f"{int(row['Rank'])}. {row['Layer']}")
        ax.text(7.15, y + 0.36, row["Role"], ha="left", va="center", fontsize=8.5, color="#37474F")
    ax.annotate("higher precedence", xy=(9.55, 4.75), xytext=(9.55, 0.35), arrowprops={"arrowstyle": "->", "color": "#264653", "lw": 1.6}, ha="center", va="center", rotation=90, color="#264653", fontsize=9)
    ax.axis("off")
    ax.set_title("Keep scientific parameters out of infrastructure config", loc="left", pad=8)
    fig.text(0.08, 0.02, "Emergency CLI overrides must be copied back into the archived parameter file", fontsize=8.5, color="#51636B")
    return save(fig, output, "74-parameter-precedence")


def execution_units(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "execution-unit-matrix.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylim(-0.7, 3.8)
    for i, row in frame.iterrows():
        y = 3 - i
        assembly_color = PALETTE[0] if row["AssemblyUnit"] == "sample" else PALETTE[3]
        coverage_color = PALETTE[4] if row["CoverageUnit"] == "sample" else PALETTE[2]
        rounded_box(ax, (-0.35, y - 0.28), 1.08, 0.56, assembly_color, row["AssemblyUnit"].upper())
        rounded_box(ax, (1.05, y - 0.28), 1.55, 0.56, coverage_color, "OWN" if row["binning_map_mode"] == "own" else str(row["binning_map_mode"]).upper())
        ax.text(2.78, y, row["UseCase"], va="center", fontsize=8.3, wrap=True)
    for y in [0.5, 1.5, 2.5]:
        ax.axhline(y, color="#EDF1F2", lw=0.8)
    ax.text(0.2, 3.55, "Assembly unit", ha="center", weight="bold")
    ax.text(1.82, 3.55, "Coverage mapping", ha="center", weight="bold")
    ax.text(2.78, 3.55, "Design implication", ha="center", weight="bold")
    ax.set_yticks(np.arange(4)[::-1], frame["Strategy"], fontsize=8.5)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("A samplesheet group does not imply co-assembly", loc="left", pad=14)
    return save(fig, output, "74-execution-unit-matrix")


def profile_separation(data: Path, output: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(9.6, 5.3))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    columns = [
        (0.25, PALETTE[0], "Scientific contract", ["params.publication.yml", "assembly / binning", "QC / database releases", "random seeds"]),
        (3.55, PALETTE[5], "Software contract", ["nf-core/mag 5.5.0", "Nextflow 26.04.0", "-profile apptainer", "immutable SIF cache"]),
        (6.85, PALETTE[4], "Infrastructure contract", ["hpc.slurm.config", "queue / account", "resource ceilings", "mounts / scratch"]),
    ]
    for x, color, title, items in columns:
        rounded_box(ax, (x, 4.8), 2.9, 0.72, color, title)
        for j, item in enumerate(items):
            rounded_box(ax, (x + 0.22, 3.85 - j * 0.78), 2.46, 0.52, "#EDF3F4", item, text_color="#263238")
        ax.annotate("", xy=(x + 1.45, 0.62), xytext=(x + 1.45, 1.25), arrowprops={"arrowstyle": "->", "color": color, "lw": 1.7})
    rounded_box(ax, (2.55, 0.08), 4.9, 0.72, PALETTE[3], "trace + report + timeline + DAG + outputs")
    ax.axis("off")
    ax.set_title("Reproducibility needs three independent locks", loc="left", pad=8)
    return save(fig, output, "74-profile-separation")


def resume_audit(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "stub-runtime-trace.tsv", sep="\t")
    counts = frame.groupby(["Run", "Status"]).size().unstack(fill_value=0).reindex(["first-success", "resume"])
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    bottom = np.zeros(len(counts))
    status_colors = {"COMPLETED": PALETTE[0], "CACHED": PALETTE[1]}
    for status in ["COMPLETED", "CACHED"]:
        values = counts[status].to_numpy() if status in counts else np.zeros(len(counts))
        bars = ax.bar(["First successful stub", "Same work + -resume"], values, bottom=bottom, color=status_colors[status], width=0.58, label=status.title())
        for bar, value, base in zip(bars, values, bottom, strict=True):
            if value:
                ax.text(bar.get_x() + bar.get_width() / 2, base + value / 2, str(int(value)), ha="center", va="center", weight="bold", color="#263238")
        bottom += values
    ax.set_ylim(0, 5.8)
    ax.set_ylabel("Tasks in official test_minimal DAG")
    ax.set_title("Resume reused 4/5 tasks; the aggregate report rebuilt", loc="left")
    ax.legend(frameon=False, ncol=2, loc="upper center")
    ax.grid(axis="y", color="#EDF1F2")
    fig.text(0.01, 0.01, "Stub mode validates orchestration and cache behavior—not biological output", fontsize=8, color="#51636B")
    return save(fig, output, "74-resume-audit")


def resource_envelope(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "hardware-envelope.tsv", sep="\t")
    memory = [15, 220, 140, 40, 36]
    labels = ["Local stub", "Queue ceiling", "GTDB-Tk", "MEGAHIT", "CheckM2 / GUNC"]
    fig, ax = plt.subplots(figsize=(8.5, 4.9))
    order = np.arange(len(labels))[::-1]
    bars = ax.barh(order, memory, color=[PALETTE[1], PALETTE[3], PALETTE[4], PALETTE[0], PALETTE[5]], height=0.62)
    ax.set_yticks(order, labels)
    ax.set_xlabel("RAM (GB; request or planning ceiling)")
    ax.set_xlim(0, 240)
    for bar, value in zip(bars, memory, strict=True):
        ax.text(value + 4, bar.get_y() + bar.get_height() / 2, f"{value} GB", va="center", weight="bold", fontsize=9)
    ax.set_title("The taxonomy step—not the launcher—sets the memory floor", loc="left")
    ax.grid(axis="x", color="#EDF1F2")
    fig.text(0.01, 0.01, "Production estimates require a pilot run on the target cohort and filesystem", fontsize=8, color="#51636B")
    return save(fig, output, "74-resource-envelope")


def database_lock(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "database-lock.tsv", sep="\t")
    gb = frame["Bytes"] / 1e9
    labels = ["GTDB R232", "CheckM2 v3", "GUNC ProGenomes 2.1"]
    fig, ax = plt.subplots(figsize=(8.3, 4.7))
    bars = ax.bar(labels, gb, color=[PALETTE[4], PALETTE[0], PALETTE[2]], width=0.58)
    ax.set_yscale("log")
    ax.set_ylim(0.8, float(gb.max()) * 2.4)
    ax.set_ylabel("Download size (GB, log scale)")
    ax.set_title("Database identity is release + file checksum", loc="left")
    for bar, value, checksum in zip(bars, gb, frame["Checksum"], strict=True):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.12, f"{value:.2f} GB\n{checksum[:8]}…", ha="center", va="bottom", fontsize=8.5)
    ax.grid(axis="y", which="both", color="#EDF1F2")
    return save(fig, output, "74-database-lock")


def provenance_bundle(data: Path, output: Path) -> list[Path]:
    frame = pd.read_csv(data / "provenance-bundle.tsv", sep="\t")
    present = frame["PresentInLocalPacket"].astype(bool)
    fig, ax = plt.subplots(figsize=(9.3, 5.5))
    y = np.arange(len(frame))[::-1]
    colors = np.where(present, PALETTE[0], PALETTE[3])
    ax.scatter(np.zeros(len(frame)), y, s=180, color=colors, edgecolor="white", linewidth=1.2)
    for yi, row, yes in zip(y, frame.itertuples(index=False), present, strict=True):
        ax.text(0, yi, "✓" if yes else "!", ha="center", va="center", color="white", weight="bold")
        ax.text(0.15, yi, row.Layer, va="center", weight="bold", fontsize=9)
        ax.text(3.1, yi, row.RequiredArtifact, va="center", fontsize=8.5, color="#455A64")
    ax.set_xlim(-0.2, 9.5)
    ax.set_ylim(-0.7, len(frame) - 0.3)
    ax.axis("off")
    ax.set_title("A successful run is not yet a publication packet", loc="left", pad=8)
    ax.text(7.3, -0.35, "green = frozen locally   red = requires full HPC run", fontsize=8, color="#51636B")
    return save(fig, output, "74-provenance-bundle")


def main() -> None:
    args = parse_args()
    data = args.input_dir.resolve()
    output = args.figures_dir.resolve()
    np.random.seed(PLOT_SEED)
    style()
    generated: list[Path] = []
    for function in (
        release_lock,
        parameter_precedence,
        execution_units,
        profile_separation,
        resume_audit,
        resource_envelope,
        database_lock,
        provenance_bundle,
    ):
        generated.extend(function(data, output))

    for source_name, figure_name in (
        ("mag-metromap-original.png", "74-mag-metromap-original.png"),
        ("funcscan-metromap-original.png", "74-funcscan-metromap-original.png"),
    ):
        destination = output / figure_name
        shutil.copy2(data / source_name, destination)
        generated.append(destination)

    manifest = {
        "article": 74,
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
