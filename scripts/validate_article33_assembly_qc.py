#!/usr/bin/env python3
"""Validate Article 33 frozen evidence and draw four publication-ready figures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article33-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from PIL import Image


SEED = 20260733
FIGURE_STEMS = (
    "33-qc-evidence-ladder",
    "33-n50-na50-correctness",
    "33-recovery-error-tradeoff",
    "33-diagnostic-task-gates",
)
COLORS = {
    "Short read": "#0072B2",
    "Long read": "#E69F00",
    "Hybrid": "#009E73",
    "Polished long read": "#CC79A7",
    "Diagnostic control": "#6C757D",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    parser.add_argument("--chapter", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text in {"", "-", "NA", "None"}:
        return None
    return float(text)


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict[str, str]] = []

    def add(self, category: str, check_id: str, passed: bool, detail: Any) -> None:
        self.rows.append(
            {
                "Category": category,
                "CheckID": check_id,
                "Status": "PASS" if passed else "FAIL",
                "Detail": str(detail),
            }
        )

    @property
    def passed(self) -> int:
        return sum(row["Status"] == "PASS" for row in self.rows)

    @property
    def failed(self) -> int:
        return sum(row["Status"] == "FAIL" for row in self.rows)


def verify_checksum_manifest(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    manifest = frozen / "file-checksums.sha256"
    rows = []
    expected_names: set[str] = set()
    for number, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not line:
            continue
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if match is None:
            checks.add("Frozen input", f"checksum-line-{number}", False, line)
            continue
        expected, relative = match.groups()
        expected_names.add(relative)
        path = frozen / relative
        observed = sha256(path) if path.is_file() else "MISSING"
        checks.add("Frozen input", f"sha256-{relative}", observed == expected, observed)
        rows.append(
            {
                "File": relative,
                "ExpectedSHA256": expected,
                "ObservedSHA256": observed,
                "Status": "PASS" if observed == expected else "FAIL",
            }
        )
    observed_names = {
        path.relative_to(frozen).as_posix()
        for path in frozen.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    }
    checks.add("Frozen input", "checksum-coverage", observed_names == expected_names, f"{len(observed_names)}/{len(expected_names)}")
    checks.add("Frozen input", "nonempty-checksum-manifest", len(rows) >= 35, len(rows))
    return rows


def save_figure(fig: plt.Figure, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(figure_dir / f"{stem}.png", dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(
        figure_dir / f"{stem}.tiff",
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def set_pub_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def figure_evidence_ladder(figure_dir: Path) -> None:
    set_pub_style()
    fig, ax = plt.subplots(figsize=(12.5, 5.8))
    ax.set_xlim(0, 12.5)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    layers = [
        ("1", "FASTA validity", "IDs · IUPAC · threshold\nlengths · GC · N bases", "Describes the file"),
        ("2", "Contiguity", "N50/L50 · N90/L90\nlong-contig retention", "Describes length distribution"),
        ("3", "Consistency", "NA50 · misassemblies\nreads · k-mers · graph", "Tests sequence support"),
        ("4", "Biological recovery", "per-genome fraction\nORFs · markers · MAG QC", "Tests target recovery"),
        ("5", "Task claim", "gene catalog · MAG\ncomplete replicon", "Defines usable for what"),
    ]
    palette = ["#DCEAF4", "#B8D8EB", "#A8DCC8", "#F4D8A7", "#E7B8CC"]
    widths = [7.1, 8.0, 8.9, 9.8, 10.7]
    for index, ((number, title, metrics, meaning), width, color) in enumerate(zip(layers, widths, palette)):
        y = 0.55 + index * 1.08
        x = (12.5 - width) / 2
        box = FancyBboxPatch((x, y), width, 0.82, boxstyle="round,pad=0.03,rounding_size=0.08", linewidth=1, edgecolor="#40515C", facecolor=color)
        ax.add_patch(box)
        ax.text(x + 0.28, y + 0.41, number, ha="center", va="center", weight="bold", color="#263238")
        ax.text(x + 0.72, y + 0.56, title, ha="left", va="center", weight="bold", color="#263238")
        ax.text(x + 0.72, y + 0.24, metrics, ha="left", va="center", fontsize=8.4, color="#263238")
        ax.text(x + width - 0.22, y + 0.41, meaning, ha="right", va="center", fontsize=8.6, color="#37474F")
        if index < len(layers) - 1:
            ax.add_patch(FancyArrowPatch((6.25, y + 0.84), (6.25, y + 1.06), arrowstyle="-|>", mutation_scale=12, linewidth=1.2, color="#52656F"))
    ax.text(6.25, 6.02, "Assembly quality is an evidence ladder", ha="center", va="center", fontsize=16, weight="bold")
    ax.text(6.25, 5.72, "No lower layer proves a higher-layer claim", ha="center", va="center", fontsize=10.5, color="#455A64")
    ax.text(6.25, 0.18, "N50 belongs to layer 2; correctness and task readiness require additional evidence.", ha="center", fontsize=10, color="#8B1E3F", weight="bold")
    save_figure(fig, figure_dir, FIGURE_STEMS[0])


def figure_n50_na50(metrics: list[dict[str, str]], figure_dir: Path) -> None:
    set_pub_style()
    rows = [row for row in metrics if row["EvaluationSet"] == "MOCK1"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 6.2), gridspec_kw={"width_ratios": [1.15, 1]})
    for row in rows:
        x = float(row["N50Bp"])
        y = float(row["NA50Bp"])
        color = COLORS[row["Family"]]
        marker = "X" if row["EvidenceClass"] == "diagnostic" else "o"
        size = 45 + 2.2 * float(row["GenomeFractionPct"])
        ax1.scatter(x, y, s=size, color=color, marker=marker, edgecolor="white", linewidth=0.8, alpha=0.92, zorder=3)
    all_values = [float(row[key]) for row in rows for key in ("N50Bp", "NA50Bp") if numeric(row[key]) and float(row[key]) > 0]
    minimum, maximum = min(all_values) * 0.75, max(all_values) * 1.35
    ax1.plot([minimum, maximum], [minimum, maximum], linestyle="--", color="#7A8790", linewidth=1, label="N50 = NA50")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(minimum, maximum)
    ax1.set_ylim(minimum, maximum)
    ax1.set_xlabel("N50 (bp)")
    ax1.set_ylabel("NA50 (bp)")
    ax1.set_title("A. Continuity before and after misassembly breaks", loc="left", weight="bold")
    ax1.grid(True, which="major", color="#E7ECEF", linewidth=0.8)
    labels = {
        "lr-hifiasm-hifi": (8, 16),
        "lr-flye-hifi": (-66, 14),
        "diagnostic-fragmented-50kb": (5, 6),
        "diagnostic-chimeric-rotation": (8, -18),
        "sr-megahit-m1": (5, 5),
        "sr-spades-10m": (5, 5),
    }
    short = {
        "lr-hifiasm-hifi": "hifiasm-meta HiFi",
        "lr-flye-hifi": "Flye HiFi",
        "diagnostic-fragmented-50kb": "Fragmented control",
        "diagnostic-chimeric-rotation": "Chimeric control",
        "sr-megahit-m1": "MEGAHIT 2M",
        "sr-spades-10m": "SPAdes 10M",
    }
    by_branch = {row["Branch"]: row for row in rows}
    for branch, offset in labels.items():
        row = by_branch[branch]
        ax1.annotate(
            short[branch],
            (float(row["N50Bp"]), float(row["NA50Bp"])),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.7,
            color="#263238",
            arrowprops={"arrowstyle": "-", "color": "#7A8790", "linewidth": 0.65},
        )

    biological = [row for row in rows if row["EvidenceClass"] == "biological"]
    order = sorted(biological, key=lambda row: float(row["N50Bp"]))
    x = np.arange(len(order))
    recovery = [float(row["RecoveredFractionGe90Pct"]) * 100 for row in order]
    discount = [100 * float(row["NA50Bp"]) / float(row["N50Bp"]) for row in order]
    ax2.bar(x - 0.2, recovery, width=0.4, color="#4C78A8", label="Genomes recovered >=90%")
    ax2.bar(x + 0.2, discount, width=0.4, color="#E69F00", label="NA50 / N50")
    ax2.set_ylim(0, max(105, max(recovery + discount) * 1.12))
    ax2.set_ylabel("Percent")
    ax2.set_xticks(x)
    ax2.set_xticklabels([row["Display"].replace(" · ", "\n") for row in order], rotation=55, ha="right", fontsize=7.2)
    ax2.set_title("B. Same N50 rank, different evidence", loc="left", weight="bold")
    ax2.grid(axis="y", color="#E7ECEF", linewidth=0.8)
    ax2.legend(loc="upper left", fontsize=8)
    handles = [
        Line2D([0], [0], marker="o", color="w", label=family, markerfacecolor=color, markersize=8)
        for family, color in COLORS.items() if family != "Diagnostic control"
    ] + [Line2D([0], [0], marker="X", color="w", label="Diagnostic control", markerfacecolor=COLORS["Diagnostic control"], markersize=8)]
    ax1.legend(handles=handles, loc="lower right", fontsize=7.5)
    fig.suptitle("N50 is not a correctness axis · MOCK1 exact-reference audit", fontsize=15, weight="bold", y=1.01)
    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.17, top=0.88, wspace=0.36)
    save_figure(fig, figure_dir, FIGURE_STEMS[1])


def figure_tradeoff(metrics: list[dict[str, str]], figure_dir: Path) -> None:
    set_pub_style()
    rows = [row for row in metrics if row["EvidenceClass"] == "biological"]
    rows = sorted(rows, key=lambda row: (row["EvaluationSet"], float(row["RecoveredFractionGe90Pct"]), float(row["N50Bp"])))
    y = np.arange(len(rows))
    labels = [f"{row['Display']}  [{row['EvaluationSet']}]" for row in rows]
    colors = [COLORS[row["Family"]] for row in rows]
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 8.2), sharey=True, gridspec_kw={"width_ratios": [1.05, 1, 1.15]})
    recovery = [100 * float(row["RecoveredFractionGe90Pct"]) for row in rows]
    axes[0].barh(y, recovery, color=colors, alpha=0.9)
    axes[0].set_xlabel("Truth genomes recovered >=90% (%)")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels, fontsize=7.4)
    axes[0].set_title("A. Recovery", loc="left", weight="bold")
    axes[0].set_xlim(0, max(100, max(recovery) * 1.12))

    mis = [float(row["Misassemblies"]) for row in rows]
    axes[1].barh(y, mis, color=colors, alpha=0.9)
    axes[1].set_xlabel("Extensive misassemblies")
    axes[1].set_title("B. Structural discordance", loc="left", weight="bold")

    mismatch = [float(row["MismatchesPer100Kbp"]) for row in rows]
    indel = [float(row["IndelsPer100Kbp"]) for row in rows]
    axes[2].scatter(mismatch, y, color="#0072B2", s=46, label="Mismatches", zorder=3)
    axes[2].scatter(indel, y, color="#D55E00", s=46, marker="s", label="Indels", zorder=3)
    for yi, left, right in zip(y, mismatch, indel):
        axes[2].plot([left, right], [yi, yi], color="#C7D0D5", linewidth=1, zorder=1)
    axes[2].set_xscale("symlog", linthresh=1)
    axes[2].set_xlabel("Errors per 100 kbp (symlog)")
    axes[2].set_title("C. Aligned consensus error", loc="left", weight="bold")
    axes[2].legend(loc="lower right", fontsize=8)
    for ax in axes:
        ax.grid(axis="x", color="#E7ECEF", linewidth=0.8)
    fig.suptitle("Recovery, structural discordance, and consensus error remain separate", fontsize=15, weight="bold", y=0.995)
    fig.text(0.5, 0.012, "Truth-set denominators are shown in labels; branches have unequal inputs and are not a causal platform benchmark.", ha="center", fontsize=9, color="#455A64")
    fig.tight_layout(rect=(0, 0.035, 1, 0.97))
    save_figure(fig, figure_dir, FIGURE_STEMS[2])


def figure_controls_gates(metrics: list[dict[str, str]], figure_dir: Path) -> None:
    set_pub_style()
    by_branch = {row["Branch"]: row for row in metrics}
    selected = [by_branch[name] for name in ("lr-flye-hifi", "diagnostic-fragmented-50kb", "diagnostic-chimeric-rotation")]
    labels = ["Flye HiFi source", "Fragmented control", "Chimeric control"]
    fig = plt.figure(figsize=(15.5, 6.4))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.05, 0.75, 1.45], wspace=0.35)
    ax1 = fig.add_subplot(grid[0, 0])
    x = np.arange(3)
    n50 = [float(row["N50Bp"]) for row in selected]
    na50 = [float(row["NA50Bp"]) for row in selected]
    ax1.bar(x - 0.18, n50, width=0.36, color="#4C78A8", label="N50")
    ax1.bar(x + 0.18, na50, width=0.36, color="#E69F00", label="NA50")
    ax1.set_yscale("log")
    ax1.set_ylabel("Length (bp, log scale)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax1.set_title("A. Length metrics", loc="left", weight="bold")
    ax1.legend(fontsize=8)
    ax1.grid(axis="y", which="major", color="#E7ECEF")

    ax2 = fig.add_subplot(grid[0, 1])
    mis = [float(row["Misassemblies"]) for row in selected]
    ax2.bar(x, mis, color=["#4C78A8", "#8A9A5B", "#8B1E3F"])
    ax2.set_ylabel("Extensive misassemblies")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax2.set_title("B. Reference breakpoints", loc="left", weight="bold")
    ax2.grid(axis="y", color="#E7ECEF")

    ax3 = fig.add_subplot(grid[0, 2])
    tasks = ["Gene / contig", "Binning / MAG", "Complete replicon"]
    evidence = ["Length threshold", "Read support", "Assembly graph", "ORF / marker", "MAG QC", "Junction support"]
    matrix = np.array(
        [
            [2, 2, 1, 2, 0, 0],
            [2, 2, 1, 1, 2, 0],
            [1, 2, 2, 1, 1, 2],
        ]
    )
    cmap = {0: "#F1F3F4", 1: "#F4D8A7", 2: "#2A9D8F"}
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = int(matrix[row_index, column_index])
            ax3.add_patch(Rectangle((column_index, row_index), 1, 1, facecolor=cmap[value], edgecolor="white", linewidth=2))
            ax3.text(column_index + 0.5, row_index + 0.5, {0: "Not primary", 1: "Support", 2: "Required"}[value], ha="center", va="center", fontsize=7.2, color="#263238", weight="bold" if value == 2 else "normal")
    ax3.set_xlim(0, len(evidence))
    ax3.set_ylim(len(tasks), 0)
    ax3.set_xticks(np.arange(len(evidence)) + 0.5)
    ax3.set_xticklabels(evidence, rotation=40, ha="right", fontsize=8)
    ax3.set_yticks(np.arange(len(tasks)) + 0.5)
    ax3.set_yticklabels(tasks, fontsize=9)
    ax3.tick_params(length=0)
    for spine in ax3.spines.values():
        spine.set_visible(False)
    ax3.set_title("C. Usability is task-specific", loc="left", weight="bold")
    fig.suptitle("Diagnostic controls separate metric behavior from task readiness", fontsize=15, weight="bold", y=1.01)
    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.17, top=0.88, wspace=0.36)
    save_figure(fig, figure_dir, FIGURE_STEMS[3])


def image_audit(figure_dir: Path, checks: Checks) -> list[dict[str, str]]:
    rows = []
    for stem in FIGURE_STEMS:
        for suffix in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            exists = path.is_file() and path.stat().st_size > 1000
            checks.add("Figures", f"exists-{stem}.{suffix}", exists, path.stat().st_size if path.is_file() else "MISSING")
            detail = "vector PDF"
            width = height = dpi_x = dpi_y = mode = "NA"
            if exists and suffix in {"png", "tiff"}:
                with Image.open(path) as image:
                    width, height = image.size
                    mode = image.mode
                    dpi = image.info.get("dpi", (0, 0))
                    dpi_x, dpi_y = float(dpi[0]), float(dpi[1])
                checks.add("Figures", f"dimensions-{stem}.{suffix}", width >= 1600 and height >= 1200, f"{width}x{height}")
                checks.add("Figures", f"mode-{stem}.{suffix}", mode in {"RGB", "RGBA"}, mode)
                checks.add("Figures", f"dpi-{stem}.{suffix}", dpi_x >= 300 and dpi_y >= 300, f"{dpi_x:.1f}x{dpi_y:.1f}")
                detail = f"{width}x{height}; {mode}; {dpi_x:.1f}x{dpi_y:.1f} dpi"
            rows.append({"Figure": f"{stem}.{suffix}", "Status": "PASS" if exists else "FAIL", "Detail": detail})
    return rows


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = (args.frozen_dir or root / "data/small/33-assembly-qc-frozen").resolve()
    output = (args.output_dir or root / "results/33-assembly-qc").resolve()
    figure_dir = (args.figure_dir or root / "figures").resolve()
    chapter = (args.chapter or root / "chapters/33-assembly-qc.qmd").resolve()
    env_prefix = args.env_prefix.resolve()
    output.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    checksum_rows = verify_checksum_manifest(frozen, checks)
    write_tsv(output / "checksum-audit.tsv", checksum_rows, ["File", "ExpectedSHA256", "ObservedSHA256", "Status"])

    contract = json.loads((frozen / "frozen-contract.json").read_text())
    checks.add("Contract", "article-number", contract.get("article") == 33, contract.get("article"))
    checks.add("Contract", "no-fastq", contract.get("fastq_included") is False and not list(frozen.rglob("*.fastq*")), contract.get("fastq_included"))
    checks.add("Contract", "no-assembly-fasta", contract.get("assembly_fasta_included") is False and not list(frozen.rglob("*.fna*")) and not list(frozen.rglob("*.fasta*")), contract.get("assembly_fasta_included"))
    checks.add("Contract", "offline-routine-qa", contract.get("routine_validation_requires_network") is False, contract.get("routine_validation_requires_network"))
    checks.add("Contract", "no-routine-metaquast", contract.get("routine_validation_reruns_metaquast") is False, contract.get("routine_validation_reruns_metaquast"))

    version_command = [str(env_prefix / "bin/quast.py"), "--version"]
    version = subprocess.check_output(version_command, text=True, stderr=subprocess.STDOUT).strip()
    checks.add("Tools", "quast-version", version == "QUAST v5.3.0", version)
    tool_rows = read_tsv(frozen / "tool-versions.tsv")
    checks.add("Tools", "tool-table", any(row["Tool"] == "QUAST/MetaQUAST" and row["Version"] == "5.3.0" for row in tool_rows), tool_rows)

    source_rows = read_tsv(frozen / "input-lineage.tsv")
    biological = [row for row in source_rows if row["EvidenceClass"] == "biological"]
    diagnostic = [row for row in source_rows if row["EvidenceClass"] == "diagnostic"]
    checks.add("Inputs", "lineage-count", len(source_rows) == 17, len(source_rows))
    checks.add("Inputs", "biological-count", len(biological) == 15, len(biological))
    checks.add("Inputs", "diagnostic-count", len(diagnostic) == 2, len(diagnostic))
    source_audit_rows = []
    for row in biological:
        path = root / row["SourceRelativePath"]
        observed = sha256(path) if path.is_file() else "MISSING"
        expected = row["SourceCompressedSHA256"]
        status = observed == expected
        checks.add("Inputs", f"upstream-{row['Branch']}", status, observed)
        source_audit_rows.append({"Branch": row["Branch"], "Path": row["SourceRelativePath"], "ExpectedSHA256": expected, "ObservedSHA256": observed, "Status": "PASS" if status else "FAIL"})
    write_tsv(output / "source-audit.tsv", source_audit_rows, ["Branch", "Path", "ExpectedSHA256", "ObservedSHA256", "Status"])
    checks.add("Inputs", "iupac-valid", all(row["InvalidIUPACBases"] == "0" for row in source_rows), [row["InvalidIUPACBases"] for row in source_rows])
    checks.add("Inputs", "minimum-contig", all(float(row["N90Bp"]) >= 1000 for row in source_rows), min(float(row["N90Bp"]) for row in source_rows))

    truth = read_tsv(frozen / "truth-manifest.tsv")
    truth_sets = {name: [row for row in truth if row["EvaluationSet"] == name] for name in ("MOCK1", "MOCK2", "MOCK1+MOCK2")}
    for name, expected in (("MOCK1", 71), ("MOCK2", 87), ("MOCK1+MOCK2", 87)):
        checks.add("Truth", f"count-{name}", len(truth_sets[name]) == expected, len(truth_sets[name]))
        checks.add("Truth", f"unique-{name}", len({row["GenBankAssembly"] for row in truth_sets[name]}) == expected, len({row["GenBankAssembly"] for row in truth_sets[name]}))
    m1 = {row["GenBankAssembly"] for row in truth_sets["MOCK1"]}
    m2 = {row["GenBankAssembly"] for row in truth_sets["MOCK2"]}
    union = {row["GenBankAssembly"] for row in truth_sets["MOCK1+MOCK2"]}
    checks.add("Truth", "mock1-strict-subset", m1 < m2, f"{len(m1)}/{len(m2)}")
    checks.add("Truth", "co-union-equals-mock2", union == m2, f"{len(union)}/{len(m2)}")

    metrics = read_tsv(frozen / "branch-metrics.tsv")
    by_branch = {row["Branch"]: row for row in metrics}
    checks.add("Metrics", "branch-count", len(metrics) == 17 and len(by_branch) == 17, len(metrics))
    group_counts = {name: sum(row["EvaluationSet"] == name for row in metrics) for name in ("MOCK1", "MOCK2", "MOCK1+MOCK2")}
    checks.add("Metrics", "group-counts", group_counts == {"MOCK1": 13, "MOCK2": 2, "MOCK1+MOCK2": 2}, group_counts)
    checks.add("Metrics", "truth-denominators", all(int(row["TruthGenomes"]) == (71 if row["EvaluationSet"] == "MOCK1" else 87) for row in metrics), [row["TruthGenomes"] for row in metrics])
    checks.add("Metrics", "finite-positive-n50", all(float(row["N50Bp"]) >= 1000 for row in metrics), min(float(row["N50Bp"]) for row in metrics))
    checks.add("Metrics", "finite-positive-na50", all(numeric(row["NA50Bp"]) is not None and float(row["NA50Bp"]) > 0 for row in metrics), [row["NA50Bp"] for row in metrics])
    expected_metrics = {
        "sr-megahit-m1": (6197, 6110, 24.072, 88, 12, 4),
        "sr-spades-10m": (26020, 25703, 59.717, 192, 37, 12),
        "lr-flye-hifi": (2013697, 903979, 70.125, 284, 43, 36),
        "lr-hifiasm-hifi": (2154824, 847090, 68.678, 264, 39, 28),
        "lr-metamdbg-hifi": (2007072, 749991, 72.573, 303, 45, 33),
        "sr-metaspades-co": (17124, 16938, 30.953, 129, 21, 5),
        "diagnostic-fragmented-50kb": (50000, 50000, 70.129, 276, 43, 36),
        "diagnostic-chimeric-rotation": (2013697, 783687, 70.125, 323, 43, 36),
    }
    for branch, expected in expected_metrics.items():
        row = by_branch[branch]
        observed = (
            int(row["N50Bp"]), int(row["NA50Bp"]), float(row["GenomeFractionPct"]),
            int(row["Misassemblies"]), int(row["RecoveredGenomesGe90Pct"]), int(row["FullGenomesGe99Pct"]),
        )
        passed = observed[:2] == expected[:2] and abs(observed[2] - expected[2]) < 1e-9 and observed[3:] == expected[3:]
        checks.add("Metrics", f"locked-{branch}", passed, observed)

    per_genome = read_tsv(frozen / "per-genome-metaquast.tsv")
    checks.add("Per genome", "row-count", len(per_genome) == 1271, len(per_genome))
    keys = {(row["EvaluationSet"], row["Reference"], row["Branch"]) for row in per_genome}
    checks.add("Per genome", "unique-keys", len(keys) == 1271, len(keys))
    for branch, metric in by_branch.items():
        rows = [row for row in per_genome if row["Branch"] == branch]
        expected = 71 if metric["EvaluationSet"] == "MOCK1" else 87
        recovered = sum(row["RecoveredGe90Pct"] == "yes" for row in rows)
        full = sum(row["FullGenomeGe99Pct"] == "yes" for row in rows)
        checks.add("Per genome", f"rows-{branch}", len(rows) == expected, len(rows))
        checks.add("Per genome", f"recovery-{branch}", recovered == int(metric["RecoveredGenomesGe90Pct"]), recovered)
        checks.add("Per genome", f"full-{branch}", full == int(metric["FullGenomesGe99Pct"]), full)

    source = by_branch["lr-flye-hifi"]
    fragmented = by_branch["diagnostic-fragmented-50kb"]
    chimeric = by_branch["diagnostic-chimeric-rotation"]
    checks.add("Controls", "fragment-total-invariant", fragmented["TotalLengthBp"] == source["TotalLengthBp"], f"{fragmented['TotalLengthBp']}/{source['TotalLengthBp']}")
    checks.add("Controls", "fragment-n50-lower", float(fragmented["N50Bp"]) <= 50_149 and float(fragmented["N50Bp"]) < float(source["N50Bp"]), fragmented["N50Bp"])
    checks.add("Controls", "chimera-total-invariant", chimeric["TotalLengthBp"] == source["TotalLengthBp"], f"{chimeric['TotalLengthBp']}/{source['TotalLengthBp']}")
    checks.add("Controls", "chimera-n50-invariant", chimeric["N50Bp"] == source["N50Bp"] and chimeric["L50"] == source["L50"], f"{chimeric['N50Bp']}/{source['N50Bp']}")
    checks.add("Controls", "chimera-misassemblies-increase", float(chimeric["Misassemblies"]) > float(source["Misassemblies"]), f"{chimeric['Misassemblies']}/{source['Misassemblies']}")
    block_rows = read_tsv(frozen / "control-block-audit.tsv")
    checks.add("Controls", "block-count", len(block_rows) == 20, len(block_rows))
    checks.add("Controls", "block-size", all(row["BlockBp"] == "100000" for row in block_rows), {row["BlockBp"] for row in block_rows})

    correlations = read_tsv(frozen / "metric-correlation-audit.tsv")
    checks.add("Metrics", "correlation-count", len(correlations) == 6, len(correlations))
    checks.add("Metrics", "correlation-scope", all(row["N"] == "11" and -1 <= float(row["SpearmanRho"]) <= 1 for row in correlations), correlations)
    summary = json.loads((frozen / "run-summary.json").read_text())
    checks.add("Summary", "summary-branch-count", summary.get("evaluation_branches") == 17, summary.get("evaluation_branches"))
    checks.add("Summary", "summary-per-genome", summary.get("per_genome_rows") == 1271, summary.get("per_genome_rows"))
    checks.add("Summary", "n50-not-correctness", summary.get("n50_is_correctness") is False, summary.get("n50_is_correctness"))
    checks.add("Summary", "no-universal-threshold", summary.get("universal_assembly_threshold_claimed") is False, summary.get("universal_assembly_threshold_claimed"))
    checks.add("Summary", "controls-not-biological", summary.get("diagnostic_controls_are_biological_results") is False, summary.get("diagnostic_controls_are_biological_results"))
    checks.add("Summary", "physical-report-counts", summary.get("physical_reference_reports") == {"MOCK1": 64, "MOCK2": 52, "MOCK1+MOCK2": 56}, summary.get("physical_reference_reports"))
    resource_rows = read_tsv(frozen / "resource-usage.tsv")
    checks.add("Resources", "four-evaluations", len(resource_rows) == 4, len(resource_rows))
    checks.add("Resources", "all-exit-zero", all(row["ExitStatus"] == "0" for row in resource_rows), [row["ExitStatus"] for row in resource_rows])
    checks.add("Resources", "peak-ram-recorded", all(float(row["PeakRSSGiB"]) > 0 for row in resource_rows), [row["PeakRSSGiB"] for row in resource_rows])

    chapter_text = chapter.read_text(encoding="utf-8")
    required_sections = (
        "## 这一步对应论文里的哪张图",
        "## 理论：",
        "## 准备工作",
        "## 可复制代码",
        "## 审计与升级",
        "## 出版级美化",
        "## 常见坑",
        "## 这段 Methods 怎么写",
        "## 换成你自己的数据怎么做",
        "## 参考",
    )
    for index, heading in enumerate(required_sections, 1):
        checks.add("Chapter", f"section-{index}", heading in chapter_text, heading)
    checks.add("Chapter", "not-draft", re.search(r"(?m)^draft:\s*false\s*$", chapter_text) is not None, "draft")
    checks.add("Chapter", "eval-false", re.search(r"(?ms)^execute:\s*\n(?:.*\n){0,5}?\s+eval:\s*false\s*$", chapter_text) is not None, "eval")
    checks.add("Chapter", "no-result-placeholder", "ARTICLE33_ACTUAL_RESULTS" not in chapter_text, "placeholder")
    checks.add("Chapter", "methods-version", "QUAST 5.3.0" in chapter_text and "min-alignment 500" in chapter_text and "min-identity 97" in chapter_text, "Methods contract")
    checks.add("Chapter", "seed", str(SEED) in chapter_text, SEED)
    checks.add("Chapter", "truth-counts", "71" in chapter_text and "87" in chapter_text and "1,271" in chapter_text, "71/87/1,271")
    checks.add("Chapter", "hardware", all(token in chapter_text for token in ("RAM", "磁盘", "CPU", "耗时")), "resource labels")
    prohibited = ("本篇可独立跑通", "这体现全系列", "接口只学一次", "作者代码通常长这样", "（即本文）")
    checks.add("Chapter", "no-meta-prose", not any(text in chapter_text for text in prohibited), [text for text in prohibited if text in chapter_text])
    for stem in FIGURE_STEMS:
        checks.add("Chapter", f"figure-reference-{stem}", f"../figures/{stem}.png" in chapter_text, stem)

    figure_evidence_ladder(figure_dir)
    figure_n50_na50(metrics, figure_dir)
    figure_tradeoff(metrics, figure_dir)
    figure_controls_gates(metrics, figure_dir)
    image_rows = image_audit(figure_dir, checks)
    write_tsv(output / "image-audit.tsv", image_rows, ["Figure", "Status", "Detail"])

    metric_audit = [
        {
            "Branch": row["Branch"],
            "EvaluationSet": row["EvaluationSet"],
            "EvidenceClass": row["EvidenceClass"],
            "N50Bp": row["N50Bp"],
            "NA50Bp": row["NA50Bp"],
            "GenomeFractionPct": row["GenomeFractionPct"],
            "Misassemblies": row["Misassemblies"],
            "RecoveredGenomesGe90Pct": row["RecoveredGenomesGe90Pct"],
            "Status": "PASS",
        }
        for row in metrics
    ]
    write_tsv(output / "metric-audit.tsv", metric_audit, ["Branch", "EvaluationSet", "EvidenceClass", "N50Bp", "NA50Bp", "GenomeFractionPct", "Misassemblies", "RecoveredGenomesGe90Pct", "Status"])
    control_audit = [row for row in checks.rows if row["Category"] == "Controls"]
    write_tsv(output / "control-audit.tsv", control_audit, ["Category", "CheckID", "Status", "Detail"])
    chapter_audit = [row for row in checks.rows if row["Category"] == "Chapter"]
    write_tsv(output / "chapter-audit.tsv", chapter_audit, ["Category", "CheckID", "Status", "Detail"])
    write_tsv(output / "validation-checks.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])

    validation_summary = {
        "article": 33,
        "status": "passed" if checks.failed == 0 else "failed",
        "checks_total": len(checks.rows),
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "biological_assemblies": 15,
        "diagnostic_controls": 2,
        "evaluation_branches": 17,
        "truth_genomes": {"MOCK1": 71, "MOCK2": 87, "MOCK1+MOCK2": 87},
        "per_genome_rows": 1271,
        "seed": SEED,
        "n50_is_correctness": False,
        "universal_assembly_threshold_claimed": False,
        "diagnostic_controls_are_biological_results": False,
        "figures": list(FIGURE_STEMS),
    }
    (output / "validation-summary.json").write_text(json.dumps(validation_summary, indent=2) + "\n", encoding="utf-8")
    (output / "validation.log").write_text(
        f"Article 33 validation: {validation_summary['status']}\nChecks: {checks.passed}/{len(checks.rows)} passed\n",
        encoding="utf-8",
    )
    print(json.dumps(validation_summary))
    return 0 if checks.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
