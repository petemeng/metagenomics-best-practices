#!/usr/bin/env python3
"""Validate Article 34 frozen evidence and draw four publication-ready figures."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Iterator, TextIO

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article34-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image


SEED = 20260734
FIGURE_STEMS = (
    "34-gene-catalog-workflow",
    "34-gene-length-distributions",
    "34-strategy-truth-audit",
    "34-threshold-method-sensitivity",
)
COLORS = {
    "Individual": "#0072B2",
    "Co-assembly": "#E69F00",
    "Mix": "#009E73",
    "Complete": "#009E73",
    "Partial": "#E69F00",
    "Incomplete": "#CC79A7",
    "MMseqs2": "#0072B2",
    "CD-HIT": "#D55E00",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
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


def decompressed_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


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
    rows: list[dict[str, str]] = []
    expected_names: set[str] = set()
    manifest = frozen / "file-checksums.sha256"
    for number, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), start=1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if match is None:
            checks.add("Frozen input", f"checksum-line-{number}", False, line)
            continue
        expected, relative = match.groups()
        expected_names.add(relative)
        path = frozen / relative
        observed = sha256(path) if path.is_file() else "MISSING"
        status = observed == expected
        checks.add("Frozen input", f"sha256-{relative}", status, observed)
        rows.append(
            {
                "File": relative,
                "ExpectedSHA256": expected,
                "ObservedSHA256": observed,
                "Status": "PASS" if status else "FAIL",
            }
        )
    observed_names = {
        path.relative_to(frozen).as_posix()
        for path in frozen.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    }
    checks.add("Frozen input", "checksum-coverage", observed_names == expected_names, f"{len(observed_names)}/{len(expected_names)}")
    checks.add("Frozen input", "payload-count", len(rows) >= 35, len(rows))
    return rows


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


def box(ax: plt.Axes, x: float, y: float, w: float, h: float, text: str, color: str) -> None:
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.03,rounding_size=0.08", facecolor=color, edgecolor="white", linewidth=1.5
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color="white", fontsize=10, fontweight="bold")


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14, color="#455A64", linewidth=1.5))


def figure_workflow(figure_dir: Path) -> None:
    set_pub_style()
    fig, ax = plt.subplots(figsize=(15.5, 6.5))
    ax.set_xlim(0, 16.5)
    ax.set_ylim(0, 7)
    ax.axis("off")
    ax.text(0.3, 6.65, "A", fontsize=14, fontweight="bold")
    ax.text(0.8, 6.65, "Paper-compatible catalog construction", fontsize=14, fontweight="bold")
    box(ax, 0.6, 4.9, 2.0, 0.9, "M1 assembly\n2M read pairs", COLORS["Individual"])
    box(ax, 0.6, 3.5, 2.0, 0.9, "M2 assembly\n2M read pairs", COLORS["Individual"])
    box(ax, 0.6, 1.5, 2.0, 0.9, "Co-assembly\n4M read pairs", COLORS["Co-assembly"])
    box(ax, 3.3, 3.6, 2.0, 1.1, "Prodigal 2.6.3\nmetagenomic mode", "#6C757D")
    for y in (5.35, 3.95, 1.95):
        arrow(ax, (2.6, y), (3.3, 4.15))
    box(ax, 6.0, 4.4, 2.2, 1.0, "Individual pool\nMMseqs2 95/95", COLORS["Individual"])
    box(ax, 6.0, 1.4, 2.2, 1.0, "Co catalog\nno extra clustering", COLORS["Co-assembly"])
    arrow(ax, (5.3, 4.15), (6.0, 4.9))
    arrow(ax, (5.3, 4.05), (6.0, 1.9))
    box(ax, 9.2, 3.0, 2.2, 1.2, "Mix catalog\nsecond 95/95 cluster", COLORS["Mix"])
    arrow(ax, (8.2, 4.9), (9.2, 3.8))
    arrow(ax, (8.2, 1.9), (9.2, 3.4))
    version_box = FancyBboxPatch(
        (12.2, 3.0), 3.6, 1.2,
        boxstyle="round,pad=0.03,rounding_size=0.08",
        facecolor="#ECEFF1", edgecolor="#90A4AE", linewidth=1.2,
    )
    ax.add_patch(version_box)
    ax.text(
        14.0, 3.6,
        "Catalog version = input lineage +\ncaller + identity + coverage +\nrepresentative policy",
        ha="center", va="center", fontsize=9, color="#37474F",
    )
    arrow(ax, (11.4, 3.6), (12.2, 3.6))
    ax.text(0.3, 0.55, "B", fontsize=14, fontweight="bold")
    box(ax, 2.1, 0.15, 2.5, 0.75, "87 exact genomes", "#5E35B1")
    box(ax, 5.25, 0.15, 2.5, 0.75, "270,679 callable ORFs", "#5E35B1")
    box(ax, 8.4, 0.15, 3.1, 0.75, "Recovery: truth coverage\nSupport: catalog coverage", "#5E35B1")
    arrow(ax, (4.6, 0.52), (5.25, 0.52))
    arrow(ax, (7.75, 0.52), (8.4, 0.52))
    save_figure(fig, figure_dir, FIGURE_STEMS[0])


def figure_lengths(rows: list[dict[str, str]], figure_dir: Path) -> None:
    set_pub_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharex=True, sharey=True)
    order = ("Co-assembly", "Individual", "Mix")
    stack_order = ("Complete", "Partial", "Incomplete")
    for ax, strategy in zip(axes, order):
        subset = [row for row in rows if row["Assembler"] == "MEGAHIT" and row["Strategy"] == strategy and int(row["LengthBinStart"]) <= 2500]
        starts = sorted({int(row["LengthBinStart"]) for row in subset})
        bottom = np.zeros(len(starts))
        for completeness in stack_order:
            values = {int(row["LengthBinStart"]): int(row["Genes"]) / 1000 for row in subset if row["Completeness"] == completeness}
            heights = np.array([values.get(start, 0.0) for start in starts])
            ax.bar(starts, heights, width=48, bottom=bottom, color=COLORS[completeness], label=completeness, linewidth=0)
            bottom += heights
        ax.set_title(strategy)
        ax.set_xlabel("Predicted gene length (bp)")
        ax.set_xlim(0, 2500)
        ax.grid(axis="y", color="#ECEFF1", linewidth=0.8)
    axes[0].set_ylabel("Genes per 50-bp bin (thousands)")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("MEGAHIT catalog length distributions at 95% identity / 95% member coverage", y=1.10, fontweight="bold")
    save_figure(fig, figure_dir, FIGURE_STEMS[1])


def figure_strategy(
    catalog_rows: list[dict[str, str]], truth_rows: list[dict[str, str]], origin_rows: list[dict[str, str]], figure_dir: Path
) -> None:
    set_pub_style()
    primary = [row for row in catalog_rows if row["CatalogID"] in {
        "megahit-individual-primary", "megahit-co-primary", "megahit-mix-primary",
        "metaspades-individual-primary", "metaspades-co-primary", "metaspades-mix-primary",
    }]
    truth = {row["CatalogID"]: row for row in truth_rows}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

    strategies = ("Individual", "Co-assembly", "Mix")
    x = np.arange(len(strategies))
    width = 0.36
    for offset, assembler in ((-width / 2, "MEGAHIT"), (width / 2, "metaSPAdes")):
        values = [int(next(row for row in primary if row["Assembler"] == assembler and row["Strategy"] == strategy)["CatalogGenes"]) / 1000 for strategy in strategies]
        axes[0].bar(x + offset, values, width, label=assembler, color="#0072B2" if assembler == "MEGAHIT" else "#CC79A7")
    axes[0].set_xticks(x, strategies, rotation=15)
    axes[0].set_ylabel("Representative genes (thousands)")
    axes[0].set_title("A  Catalog size")
    axes[0].legend()
    axes[0].grid(axis="y", color="#ECEFF1")

    markers = {"MEGAHIT": "o", "metaSPAdes": "s"}
    for row in primary:
        audit = truth[row["CatalogID"]]
        axes[1].scatter(
            float(audit["TruthNRRecoveryPct"]), float(audit["CatalogTruthSupportPct"]),
            s=int(row["CatalogGenes"]) / 350, color=COLORS[row["Strategy"]], marker=markers[row["Assembler"]],
            edgecolor="white", linewidth=0.8, alpha=0.9,
        )
        axes[1].annotate(f"{row['Assembler']}\n{row['Strategy']}", (float(audit["TruthNRRecoveryPct"]), float(audit["CatalogTruthSupportPct"])), xytext=(4, 4), textcoords="offset points", fontsize=7)
    axes[1].set_xlabel("Callable truth clusters recovered (%)")
    axes[1].set_ylabel("Catalog representatives truth-supported (%)")
    axes[1].set_title("B  Recovery is not support")
    axes[1].grid(color="#ECEFF1")

    origin_counts: dict[tuple[str, str], int] = Counter()
    for row in origin_rows:
        origin_counts[(row["Assembler"], row["RepresentativeOrigin"])] += int(row["Genes"])
    assemblers = ("MEGAHIT", "metaSPAdes")
    bottoms = np.zeros(2)
    for origin in ("Individual", "Co-assembly"):
        values = np.array([origin_counts[(assembler, origin)] for assembler in assemblers])
        totals = np.array([sum(origin_counts[(assembler, item)] for item in ("Individual", "Co-assembly")) for assembler in assemblers])
        shares = 100 * values / totals
        axes[2].bar(assemblers, shares, bottom=bottoms, color=COLORS[origin], label=origin)
        bottoms += shares
    axes[2].set_ylim(0, 100)
    axes[2].set_ylabel("Mix representatives (%)")
    axes[2].set_title("C  Representative origin")
    axes[2].legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    fig.subplots_adjust(bottom=0.20, wspace=0.30)
    fig.suptitle("Catalog strategy changes size, evidence, and representative provenance", y=1.02, fontweight="bold")
    save_figure(fig, figure_dir, FIGURE_STEMS[2])


def figure_sensitivity(
    catalog_rows: list[dict[str, str]], truth_rows: list[dict[str, str]], genome_rows: list[dict[str, str]], figure_dir: Path
) -> None:
    set_pub_style()
    truth = {row["CatalogID"]: row for row in truth_rows}
    catalog = {row["CatalogID"]: row for row in catalog_rows}
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.1))

    for strategy in ("Individual", "Co-assembly", "Mix"):
        subset = [row for row in genome_rows if row["Assembler"] == "MEGAHIT" and row["Strategy"] == strategy and row["CatalogID"].endswith("primary")]
        x = np.array([float(row["ExpectedAbundanceMeanPct"]) for row in subset])
        y = np.array([float(row["GeneRecoveryPct"]) for row in subset])
        axes[0].scatter(x, y, s=14, alpha=0.38, color=COLORS[strategy], label=strategy)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Expected mean DNA abundance (%)")
    axes[0].set_ylabel("Per-genome callable gene recovery (%)")
    axes[0].set_title("A  Recovery remains abundance-dependent")
    axes[0].grid(color="#ECEFF1", which="both")
    axes[0].legend()

    ids = ("megahit-mix-id90-cov80", "megahit-mix-primary", "megahit-mix-id99-cov95", "megahit-mix-cdhit")
    labels = ("MMseqs2\n90/80", "MMseqs2 95/95", "MMseqs2\n99/95", "CD-HIT 95/95")
    annotation_offsets = {
        "megahit-mix-id90-cov80": (8, 8),
        "megahit-mix-primary": (20, 6),
        "megahit-mix-id99-cov95": (8, 8),
        "megahit-mix-cdhit": (20, -20),
    }
    for catalog_id, label in zip(ids, labels):
        row = catalog[catalog_id]
        audit = truth[catalog_id]
        method = "CD-HIT" if "cdhit" in catalog_id else "MMseqs2"
        marker = "D" if method == "CD-HIT" else "o"
        size = 70 if method == "CD-HIT" else 145
        point = (int(row["CatalogGenes"]) / 1000, float(audit["TruthNRRecoveryPct"]))
        axes[1].scatter(
            *point, s=size, marker=marker, color=COLORS[method], edgecolor="white",
            linewidth=1, zorder=3 if method == "CD-HIT" else 2,
        )
        axes[1].annotate(
            label, point, xytext=annotation_offsets[catalog_id], textcoords="offset points", fontsize=8,
            arrowprops={"arrowstyle": "-", "color": "#78909C", "linewidth": 0.7}
            if catalog_id in {"megahit-mix-primary", "megahit-mix-cdhit"} else None,
        )
    axes[1].set_xlabel("MEGAHIT mix representatives (thousands)")
    axes[1].set_ylabel("Callable truth clusters recovered (%)")
    axes[1].set_title("B  Threshold and software are model choices")
    axes[1].grid(color="#ECEFF1")
    fig.suptitle("Truth-aware sensitivity analysis prevents a single catalog definition from becoming invisible", y=1.02, fontweight="bold")
    save_figure(fig, figure_dir, FIGURE_STEMS[3])


def image_audit(figure_dir: Path, checks: Checks) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for stem in FIGURE_STEMS:
        png = figure_dir / f"{stem}.png"
        pdf = figure_dir / f"{stem}.pdf"
        tiff = figure_dir / f"{stem}.tiff"
        status = png.is_file() and pdf.is_file() and tiff.is_file()
        detail = "missing"
        if status:
            with Image.open(png) as image:
                width, height = image.size
                colors = image.convert("RGB").getcolors(maxcolors=10_000_000)
                color_count = len(colors) if colors else 10_000_001
            status = width >= 1600 and height >= 900 and color_count > 20 and pdf.stat().st_size > 5_000 and tiff.stat().st_size > 10_000
            detail = f"{width}x{height}; colors={color_count}; pdf={pdf.stat().st_size}; tiff={tiff.stat().st_size}"
        checks.add("Images", stem, status, detail)
        rows.append({"Figure": stem, "Status": "PASS" if status else "FAIL", "Detail": detail})
    return rows


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = (args.frozen_dir or root / "data/small/34-nonredundant-gene-catalog-frozen").resolve()
    output = (args.output_dir or root / "results/34-nonredundant-gene-catalog").resolve()
    figure_dir = (args.figure_dir or root / "figures").resolve()
    chapter = (args.chapter or root / "chapters/34-nonredundant-gene-catalog.qmd").resolve()
    output.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    checksum_rows = verify_checksum_manifest(frozen, checks)
    write_tsv(output / "checksum-audit.tsv", checksum_rows, ["File", "ExpectedSHA256", "ObservedSHA256", "Status"])

    contract = json.loads((frozen / "frozen-contract.json").read_text(encoding="utf-8"))
    checks.add("Contract", "article", contract.get("article") == 34, contract.get("article"))
    checks.add("Contract", "seed", contract.get("seed") == SEED, contract.get("seed"))
    checks.add("Contract", "truth-genomes", contract.get("truth_genomes") == 87, contract.get("truth_genomes"))
    checks.add("Contract", "benchmark-commit", contract.get("benchmark_commit") == "a429a3724d4593f35b8d7323b20252a6be90e1cd", contract.get("benchmark_commit"))

    lineage = read_tsv(frozen / "input-lineage.tsv")
    checks.add("Inputs", "six-assemblies", len(lineage) == 6, len(lineage))
    checks.add("Inputs", "two-assemblers", {row["Assembler"] for row in lineage} == {"MEGAHIT", "metaSPAdes"}, {row["Assembler"] for row in lineage})
    checks.add("Inputs", "strategy-counts", Counter(row["Strategy"] for row in lineage) == {"Individual": 4, "Co-assembly": 2}, Counter(row["Strategy"] for row in lineage))
    checks.add("Inputs", "contig-total", sum(int(row["Contigs"]) for row in lineage) == 96_275, sum(int(row["Contigs"]) for row in lineage))
    checks.add("Inputs", "assembly-bases", sum(int(row["AssemblyBases"]) for row in lineage) == 398_103_372, sum(int(row["AssemblyBases"]) for row in lineage))
    source_rows: list[dict[str, str]] = []
    for row in lineage:
        path = root / row["SourcePath"]
        observed = sha256(path) if path.is_file() else "MISSING"
        status = observed == row["SourceCompressedSHA256"]
        checks.add("Inputs", f"upstream-{row['Branch']}", status, observed)
        source_rows.append({"Branch": row["Branch"], "Path": row["SourcePath"], "ExpectedSHA256": row["SourceCompressedSHA256"], "ObservedSHA256": observed, "Status": "PASS" if status else "FAIL"})
    write_tsv(output / "source-audit.tsv", source_rows, ["Branch", "Path", "ExpectedSHA256", "ObservedSHA256", "Status"])

    truth_genomes = read_tsv(frozen / "truth-genomes.tsv")
    checks.add("Truth", "genome-count", len(truth_genomes) == 87, len(truth_genomes))
    checks.add("Truth", "mock1-subset-count", sum(row["InMOCK1"] == "yes" for row in truth_genomes) == 71, sum(row["InMOCK1"] == "yes" for row in truth_genomes))
    checks.add("Truth", "reference-bases", sum(int(row["ReferenceBases"]) for row in truth_genomes) == 292_802_018, sum(int(row["ReferenceBases"]) for row in truth_genomes))

    prediction = read_tsv(frozen / "gene-prediction-summary.tsv")
    by_branch = {row["Branch"]: row for row in prediction}
    checks.add("Prediction", "branch-count", len(prediction) == 7, len(prediction))
    expected_counts = {
        "megahit-m1": 62_571, "megahit-m2": 61_635, "megahit-co": 91_985,
        "metaspades-m1": 65_109, "metaspades-m2": 64_135, "metaspades-co": 95_972,
        "truth": 270_679,
    }
    for branch, expected in expected_counts.items():
        checks.add("Prediction", f"genes-{branch}", int(by_branch[branch]["Genes"]) == expected, by_branch[branch]["Genes"])
    checks.add("Prediction", "assembly-gene-total", sum(int(row["Genes"]) for row in prediction if row["Branch"] != "truth") == 441_407, sum(int(row["Genes"]) for row in prediction if row["Branch"] != "truth"))
    checks.add("Prediction", "class-conservation", all(int(row["Genes"]) == int(row["CompleteGenes"]) + int(row["PartialGenes"]) + int(row["IncompleteGenes"]) for row in prediction), "all branches")

    catalog_rows = read_tsv(frozen / "catalog-summary.tsv")
    catalog = {row["CatalogID"]: row for row in catalog_rows}
    checks.add("Catalog", "catalog-count", len(catalog_rows) == 9 and len(catalog) == 9, len(catalog_rows))
    primary_ids = {
        "megahit-individual-primary", "megahit-co-primary", "megahit-mix-primary",
        "metaspades-individual-primary", "metaspades-co-primary", "metaspades-mix-primary",
    }
    checks.add("Catalog", "primary-six", primary_ids <= set(catalog), sorted(primary_ids - set(catalog)))
    checks.add("Catalog", "primary-parameters", all(catalog[item]["MinimumAminoAcidIdentity"] == "0.95" and catalog[item]["MinimumMemberCoverage"] == "0.95" for item in primary_ids if catalog[item]["Strategy"] != "Co-assembly"), "95/95")
    checks.add("Catalog", "megahit-mix-raw", int(catalog["megahit-mix-primary"]["RawInputGenes"]) == 216_191, catalog["megahit-mix-primary"]["RawInputGenes"])
    checks.add("Catalog", "megahit-mix-representatives", int(catalog["megahit-mix-primary"]["CatalogGenes"]) == 93_782, catalog["megahit-mix-primary"]["CatalogGenes"])
    checks.add("Catalog", "compression-conservation", all(int(row["RawInputGenes"]) == int(row["CatalogGenes"]) + int(row["RemovedRedundantGenes"]) for row in catalog_rows), "all catalogs")

    truth_rows = read_tsv(frozen / "truth-audit-summary.tsv")
    truth = {row["CatalogID"]: row for row in truth_rows}
    checks.add("Truth audit", "catalog-count", len(truth_rows) == 9, len(truth_rows))
    checks.add("Truth audit", "truth-denominator", {row["TruthNRClusters"] for row in truth_rows} == {"260868"}, {row["TruthNRClusters"] for row in truth_rows})
    checks.add("Truth audit", "megahit-mix-recovery", abs(float(truth["megahit-mix-primary"]["TruthNRRecoveryPct"]) - 24.740865) < 1e-5, truth["megahit-mix-primary"]["TruthNRRecoveryPct"])
    checks.add("Truth audit", "megahit-mix-support", abs(float(truth["megahit-mix-primary"]["CatalogTruthSupportPct"]) - 92.632915) < 1e-5, truth["megahit-mix-primary"]["CatalogTruthSupportPct"])
    checks.add("Truth audit", "mix-gains-recovery", float(truth["megahit-mix-primary"]["TruthNRRecoveryPct"]) > float(truth["megahit-co-primary"]["TruthNRRecoveryPct"]) > float(truth["megahit-individual-primary"]["TruthNRRecoveryPct"]), [truth[item]["TruthNRRecoveryPct"] for item in ("megahit-individual-primary", "megahit-co-primary", "megahit-mix-primary")])

    per_genome = read_tsv(frozen / "per-genome-recovery.tsv")
    checks.add("Per genome", "row-count", len(per_genome) == 783, len(per_genome))
    checks.add("Per genome", "unique-keys", len({(row["CatalogID"], row["GenBankAssembly"]) for row in per_genome}) == 783, len({(row["CatalogID"], row["GenBankAssembly"]) for row in per_genome}))
    checks.add("Per genome", "bounded-recovery", all(0 <= float(row["GeneRecoveryPct"]) <= 100 for row in per_genome), "0..100")

    agreement = read_tsv(frozen / "method-agreement.tsv")
    checks.add("Sensitivity", "method-row", len(agreement) == 1, len(agreement))
    checks.add("Sensitivity", "same-raw-universe", int(agreement[0]["RawGenes"]) == 216_191, agreement[0]["RawGenes"])
    checks.add("Sensitivity", "bounded-agreement", all(0 <= float(agreement[0][field]) <= 1 for field in ("RepresentativeJaccard", "CoClusterPairPrecision", "CoClusterPairRecall", "CoClusterPairF1")), agreement[0])
    checks.add("Sensitivity", "cdhit-not-identical", int(agreement[0]["ClustersA"]) != int(agreement[0]["ClustersB"]), f"{agreement[0]['ClustersA']}/{agreement[0]['ClustersB']}")

    reps_meta = read_tsv(frozen / "primary-catalog-representatives.tsv.gz")
    membership = read_tsv(frozen / "primary-catalog-membership.tsv.gz")
    checks.add("Primary payload", "representative-metadata", len(reps_meta) == 93_782, len(reps_meta))
    checks.add("Primary payload", "membership-rows", len(membership) == 216_191, len(membership))
    checks.add("Primary payload", "unique-members", len({row["MemberID"] for row in membership}) == 216_191, len({row["MemberID"] for row in membership}))
    counts = Counter(row["RepresentativeID"] for row in membership)
    checks.add("Primary payload", "cluster-size-ledger", all(counts[row["RepresentativeID"]] == int(row["ClusterSize"]) for row in reps_meta), "all representatives")

    faa = frozen / "catalog/megahit-mix-primary.faa.gz"
    fna = frozen / "catalog/megahit-mix-primary.fna.gz"
    faa_ids = set()
    fna_ids = set()
    valid_faa = valid_fna = True
    divisible = True
    for header, sequence in fasta_records(faa):
        faa_ids.add(header.split()[0])
        valid_faa &= set(sequence.upper()) <= set("ABCDEFGHIKLMNPQRSTVWXYZ*UOJ")
    for header, sequence in fasta_records(fna):
        fna_ids.add(header.split()[0])
        valid_fna &= set(sequence.upper()) <= set("ACGTNRYKMSWBDHV")
        divisible &= len(sequence) % 3 == 0
    checks.add("Primary payload", "faa-count", len(faa_ids) == 93_782, len(faa_ids))
    checks.add("Primary payload", "fna-count", len(fna_ids) == 93_782, len(fna_ids))
    checks.add("Primary payload", "paired-identifiers", faa_ids == fna_ids == {row["RepresentativeID"] for row in reps_meta}, len(faa_ids & fna_ids))
    checks.add("Primary payload", "valid-alphabets", valid_faa and valid_fna, f"faa={valid_faa}; fna={valid_fna}")
    checks.add("Primary payload", "codon-length", divisible, divisible)
    checks.add("Primary payload", "uncompressed-fna-sha", decompressed_sha256(fna) == catalog["megahit-mix-primary"]["RepresentativeFNA_SHA256"], decompressed_sha256(fna))

    resources = read_tsv(frozen / "resource-usage.tsv")
    checks.add("Resources", "many-recorded-steps", len(resources) >= 140, len(resources))
    checks.add("Resources", "all-exit-zero", all(row["ExitStatus"] == "0" for row in resources), sorted({row["ExitStatus"] for row in resources}))
    checks.add("Resources", "rss-recorded", all(float(row["MaximumRSSKiB"]) > 0 for row in resources), min(float(row["MaximumRSSKiB"]) for row in resources))

    chapter_text = chapter.read_text(encoding="utf-8")
    required_sections = (
        "## 这一步对应论文里的哪张图", "## 理论：", "## 准备工作", "## 可复制代码",
        "## 审计与升级", "## 出版级美化", "## 常见坑", "## 这段 Methods 怎么写",
        "## 换成你自己的数据怎么做", "## 参考",
    )
    for index, heading in enumerate(required_sections, start=1):
        checks.add("Chapter", f"section-{index}", heading in chapter_text, heading)
    checks.add("Chapter", "not-draft", re.search(r"(?m)^draft:\s*false\s*$", chapter_text) is not None, "draft")
    checks.add("Chapter", "eval-false", re.search(r"(?ms)^execute:\s*\n(?:.*\n){0,5}?\s+eval:\s*false\s*$", chapter_text) is not None, "eval")
    checks.add("Chapter", "versions", all(token in chapter_text for token in ("Prodigal 2.6.3", "MMseqs2 9.d36de", "CD-HIT 4.8.1")), "tool versions")
    checks.add("Chapter", "parameters", all(token in chapter_text for token in ("--min-seq-id 0.95", "--cov-mode 1", "--cluster-mode 2", "-aS 0.95")), "cluster parameters")
    checks.add("Chapter", "seed", str(SEED) in chapter_text, SEED)
    checks.add("Chapter", "locked-counts", all(token in chapter_text for token in ("441,407", "270,679", "260,868", "93,782", "216,191")), "headline counts")
    checks.add("Chapter", "hardware", all(token in chapter_text for token in ("RAM", "磁盘", "CPU", "耗时")), "resource labels")
    checks.add("Chapter", "inline-theme", all(token in chapter_text for token in ("install.packages", "pal_pub", "scale_color_pub", "scale_fill_pub", "theme_pub", "save_pub")), "inline preparation")
    prohibited = ("本篇可独立跑通", "这体现全系列", "接口只学一次", "作者代码通常长这样", "（即本文）")
    checks.add("Chapter", "no-meta-prose", not any(text in chapter_text for text in prohibited), [text for text in prohibited if text in chapter_text])
    for stem in FIGURE_STEMS:
        checks.add("Chapter", f"figure-reference-{stem}", f"../figures/{stem}.png" in chapter_text, stem)

    length_rows = read_tsv(frozen / "gene-length-histogram.tsv")
    origin_rows = read_tsv(frozen / "mix-origin-summary.tsv")
    figure_workflow(figure_dir)
    figure_lengths(length_rows, figure_dir)
    figure_strategy(catalog_rows, truth_rows, origin_rows, figure_dir)
    figure_sensitivity(catalog_rows, truth_rows, per_genome, figure_dir)
    image_rows = image_audit(figure_dir, checks)
    write_tsv(output / "image-audit.tsv", image_rows, ["Figure", "Status", "Detail"])

    write_tsv(output / "catalog-audit.tsv", catalog_rows, list(catalog_rows[0]))
    write_tsv(output / "truth-audit.tsv", truth_rows, list(truth_rows[0]))
    write_tsv(output / "chapter-audit.tsv", [row for row in checks.rows if row["Category"] == "Chapter"], ["Category", "CheckID", "Status", "Detail"])
    write_tsv(output / "validation-checks.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])
    summary = {
        "article": 34,
        "status": "passed" if checks.failed == 0 else "failed",
        "checks_total": len(checks.rows),
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "assembly_branches": 6,
        "predicted_assembly_genes": 441_407,
        "truth_genomes": 87,
        "truth_genes": 270_679,
        "truth_nr_clusters": 260_868,
        "primary_catalog_genes": 93_782,
        "primary_truth_recovery_pct": float(truth["megahit-mix-primary"]["TruthNRRecoveryPct"]),
        "primary_catalog_support_pct": float(truth["megahit-mix-primary"]["CatalogTruthSupportPct"]),
        "cdhit_catalog_genes": int(catalog["megahit-mix-cdhit"]["CatalogGenes"]),
        "representative_jaccard": float(agreement[0]["RepresentativeJaccard"]),
        "catalog_size_is_gene_richness": False,
        "universal_mix_winner_claimed": False,
        "truth_absence_is_false_positive": False,
        "seed": SEED,
        "figures": list(FIGURE_STEMS),
    }
    (output / "validation-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "validation.log").write_text(
        f"Article 34 validation: {summary['status']}\nChecks: {checks.passed}/{len(checks.rows)} passed\n",
        encoding="utf-8",
    )
    print(json.dumps(summary))
    return 0 if checks.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
