#!/usr/bin/env python3
"""Shared checksum, chapter, figure, and scientific-contract validation for Articles 37–40."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article37-40-matplotlib")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#4D4D4D", "#BDBDBD"]
SAMPLE_COLORS = {
    "Catalog genes": "#4D4D4D", "BGC regions": "#4D4D4D", "BGC genes": "#56B4E9",
    "MOCK1 reads": "#0072B2", "MOCK2 reads": "#D55E00",
}
FIGURES = {
    37: ["37-tool-consensus", "37-cazyme-class-profile", "37-family-abundance", "37-cgc-substrate"],
    38: ["38-evidence-tiers", "38-primary-hit-quality", "38-drug-class-profile", "38-positive-controls"],
    39: ["39-threshold-database-sensitivity", "39-hit-quality", "39-vfc-category-profile", "39-context-controls"],
    40: ["40-tool-bgc-yield", "40-fragmentation-sensitivity", "40-bgc-type-profile", "40-mibig-similarity-abundance"],
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        names = fields or list(rows[0])
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass
class Checks:
    rows: list[dict[str, object]] = field(default_factory=list)

    def add(self, category: str, check: str, passed: bool, detail: object) -> None:
        self.rows.append({"Category": category, "CheckID": check, "Status": "PASS" if passed else "FAIL", "Detail": json.dumps(detail, ensure_ascii=False) if not isinstance(detail, str) else detail})

    @property
    def failed(self) -> int:
        return sum(row["Status"] == "FAIL" for row in self.rows)

    @property
    def passed(self) -> int:
        return sum(row["Status"] == "PASS" for row in self.rows)


def checksum_audit(frozen: Path, checks: Checks) -> list[dict[str, str]]:
    rows = []
    manifest = frozen / "file-checksums.sha256"
    for line in manifest.read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        path = frozen / relative
        observed = sha256(path) if path.is_file() else "MISSING"
        status = expected == observed
        checks.add("Checksum", relative, status, observed)
        rows.append({"File": relative, "ExpectedSHA256": expected, "ObservedSHA256": observed, "Status": "PASS" if status else "FAIL"})
    listed = {row["File"] for row in rows}
    actual = {p.relative_to(frozen).as_posix() for p in frozen.rglob("*") if p.is_file() and p.name != "file-checksums.sha256"}
    checks.add("Checksum", "manifest-complete", listed == actual, {"listed": len(listed), "actual": len(actual)})
    return rows


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.labelsize": 10, "axes.titlesize": 11,
        "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8.5,
        "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def save(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=350, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def grouped_bar(ax: plt.Axes, labels: list[str], series: list[tuple[str, list[float]]], xlabel: str, horizontal: bool = True) -> None:
    positions = np.arange(len(labels)); width = 0.8 / len(series)
    for index, (name, values) in enumerate(series):
        offset = (index - (len(series) - 1) / 2) * width
        color = SAMPLE_COLORS.get(name, COLORS[index])
        if horizontal:
            ax.barh(positions + offset, values, height=width, label=name, color=color)
        else:
            ax.bar(positions + offset, values, width=width, label=name, color=color)
    if horizontal:
        ax.set_yticks(positions, labels); ax.set_xlabel(xlabel); ax.grid(axis="x", color="#E6E6E6")
    else:
        ax.set_xticks(positions, labels, rotation=20, ha="right"); ax.set_ylabel(xlabel); ax.grid(axis="y", color="#E6E6E6")
    ax.legend(frameon=False)


def draw37(frozen: Path, figure_dir: Path) -> None:
    evidence = read_tsv(frozen / "evidence-tier-summary.tsv")
    labels = [row["EvidenceTier"] for row in evidence][::-1]
    fig, ax = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
    grouped_bar(ax, labels, [
        ("Catalog genes", [float(row["GenePercent"]) for row in evidence][::-1]),
        ("MOCK1 reads", [float(row["MOCK1ReadPercent"]) for row in evidence][::-1]),
        ("MOCK2 reads", [float(row["MOCK2ReadPercent"]) for row in evidence][::-1]),
    ], "Genes or assigned read mass (%)")
    save(fig, figure_dir / "37-tool-consensus")

    classes = read_tsv(frozen / "cazyme-class-summary.tsv")
    gene_total = sum(float(r["FractionalGeneEquivalent"]) for r in classes)
    m1_total = sum(float(r["MOCK1FractionalRawReads"]) for r in classes)
    m2_total = sum(float(r["MOCK2FractionalRawReads"]) for r in classes)
    fig, ax = plt.subplots(figsize=(8.2, 4.7), constrained_layout=True)
    grouped_bar(ax, [r["CAZyClass"] for r in classes][::-1], [
        ("Catalog genes", [100 * float(r["FractionalGeneEquivalent"]) / gene_total for r in classes][::-1]),
        ("MOCK1 reads", [100 * float(r["MOCK1FractionalRawReads"]) / m1_total for r in classes][::-1]),
        ("MOCK2 reads", [100 * float(r["MOCK2FractionalRawReads"]) / m2_total for r in classes][::-1]),
    ], "Fractional share within consensus CAZyme mass (%)")
    save(fig, figure_dir / "37-cazyme-class-profile")

    all_families = read_tsv(frozen / "cazyme-family-summary.tsv")
    families = all_families[:15][::-1]
    family_totals = {
        "Catalog genes": sum(float(r["FractionalGeneEquivalent"]) for r in all_families),
        "MOCK1 reads": sum(float(r["MOCK1FractionalRawReads"]) for r in all_families),
        "MOCK2 reads": sum(float(r["MOCK2FractionalRawReads"]) for r in all_families),
    }
    fig, ax = plt.subplots(figsize=(8.5, 6.2), constrained_layout=True)
    grouped_bar(ax, [r["CAZyFamily"] for r in families], [
        ("Catalog genes", [100 * float(r["FractionalGeneEquivalent"]) / family_totals["Catalog genes"] for r in families]),
        ("MOCK1 reads", [100 * float(r["MOCK1FractionalRawReads"]) / family_totals["MOCK1 reads"] for r in families]),
        ("MOCK2 reads", [100 * float(r["MOCK2FractionalRawReads"]) / family_totals["MOCK2 reads"] for r in families]),
    ], "Fractional share within family-assigned mass (%)")
    save(fig, figure_dir / "37-family-abundance")

    cgc = read_tsv(frozen / "btheta-cgc-summary.tsv")[:7]
    substrates = read_tsv(frozen / "btheta-substrate-summary.tsv")[:10][::-1]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), constrained_layout=True)
    axes[0].barh(np.arange(len(cgc)), [int(r["Genes"]) for r in cgc], color="#0072B2")
    axes[0].set_yticks(np.arange(len(cgc)), [r["GeneType"] for r in cgc]); axes[0].invert_yaxis(); axes[0].set_xlabel("Genes in predicted CGCs"); axes[0].set_title("CGC gene types")
    axes[1].barh(np.arange(len(substrates)), [int(r["CGCs"]) for r in substrates], color="#D55E00")
    axes[1].set_yticks(np.arange(len(substrates)), [r["PredictedSubstrate"] for r in substrates]); axes[1].set_xlabel("CGCs with substrate prediction"); axes[1].set_title("Top predicted substrates")
    for ax in axes: ax.grid(axis="x", color="#E6E6E6")
    save(fig, figure_dir / "37-cgc-substrate")


def draw38(frozen: Path, figure_dir: Path) -> None:
    tiers = read_tsv(frozen / "evidence-tier-summary.tsv")
    datasets = ["catalog", "coassembly", "pseudomonas", "staphylococcus"]
    fig, ax = plt.subplots(figsize=(8.8, 4.8), constrained_layout=True)
    x = np.arange(len(datasets)); width = 0.24
    for index, tier in enumerate(("Perfect", "Strict", "Loose")):
        values = [int(next(r for r in tiers if r["Dataset"] == ds and r["EvidenceTier"] == tier)["UniqueORFs"]) for ds in datasets]
        ax.bar(x + (index - 1) * width, values, width=width, label=tier, color=COLORS[index])
    ax.set_xticks(x, ["Catalog proteins", "Co-assembly", "P. aeruginosa", "S. aureus"])
    ax.set_yscale("symlog", linthresh=10); ax.set_ylabel("Unique ORFs (symlog scale)"); ax.legend(frameon=False); ax.grid(axis="y", color="#E6E6E6")
    save(fig, figure_dir / "38-evidence-tiers")

    quality = read_tsv(frozen / "primary-hit-quality.tsv")
    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    for tier, color in (("Perfect", "#009E73"), ("Strict", "#0072B2")):
        selected = [r for r in quality if r["EvidenceTier"] == tier]
        if selected:
            ax.scatter([float(r["IdentityPercent"]) for r in selected], [float(r["ReferenceLengthPercent"]) for r in selected], s=35, alpha=0.8, color=color, label=tier)
    ax.set_xlabel("Best-hit amino-acid identity (%)"); ax.set_ylabel("Length relative to CARD reference (%)"); ax.legend(frameon=False); ax.grid(color="#E6E6E6")
    save(fig, figure_dir / "38-primary-hit-quality")

    all_classes = read_tsv(frozen / "drug-class-summary.tsv")
    classes = all_classes[:12][::-1]
    class_totals = {
        "Catalog genes": sum(float(r["FractionalGeneEquivalent"]) for r in all_classes),
        "MOCK1 reads": sum(float(r["MOCK1FractionalRawReads"]) for r in all_classes),
        "MOCK2 reads": sum(float(r["MOCK2FractionalRawReads"]) for r in all_classes),
    }
    fig, ax = plt.subplots(figsize=(10.2, 6.2), constrained_layout=True)
    grouped_bar(ax, [r["DrugClass"] for r in classes], [
        ("Catalog genes", [100 * float(r["FractionalGeneEquivalent"]) / class_totals["Catalog genes"] for r in classes]),
        ("MOCK1 reads", [100 * float(r["MOCK1FractionalRawReads"]) / class_totals["MOCK1 reads"] for r in classes]),
        ("MOCK2 reads", [100 * float(r["MOCK2FractionalRawReads"]) / class_totals["MOCK2 reads"] for r in classes]),
    ], "Fractional share within drug-class-assigned mass (%)")
    save(fig, figure_dir / "38-drug-class-profile")

    controls = read_tsv(frozen / "positive-control-hits.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), constrained_layout=True)
    for ax, control in zip(axes, ("pseudomonas", "staphylococcus")):
        selected = [r for r in controls if r["Control"] == control]
        counts = [sum(r["EvidenceTier"] == tier for r in selected) for tier in ("Perfect", "Strict")]
        bars = ax.bar(["Perfect", "Strict"], counts, color=["#009E73", "#0072B2"])
        ax.set_title("P. aeruginosa ATCC 9027" if control == "pseudomonas" else "S. aureus USA300"); ax.set_ylabel("Primary CARD hits"); ax.grid(axis="y", color="#E6E6E6")
        for bar, value in zip(bars, counts): ax.text(bar.get_x() + bar.get_width()/2, value, str(value), ha="center", va="bottom")
    save(fig, figure_dir / "38-positive-controls")


def draw39(frozen: Path, figure_dir: Path) -> None:
    sensitivity = read_tsv(frozen / "sensitivity-summary.tsv")
    fig, ax = plt.subplots(figsize=(7.7, 4.5), constrained_layout=True)
    labels = [r["Branch"] for r in sensitivity]
    shared = [int(r["SharedWithCorePrimary"]) for r in sensitivity]
    additional = [int(r["AdditionalVsCorePrimary"]) for r in sensitivity]
    x = np.arange(len(labels)); ax.bar(x, shared, color="#0072B2", label="Shared with core 90/80"); ax.bar(x, additional, bottom=shared, color="#E69F00", label="Additional")
    ax.set_xticks(x, labels, rotation=15, ha="right"); ax.set_ylabel("Catalog genes"); ax.legend(frameon=False); ax.grid(axis="y", color="#E6E6E6")
    save(fig, figure_dir / "39-threshold-database-sensitivity")

    hits = [r for r in read_tsv(frozen / "all-hit-quality.tsv") if r["Branch"] == "Core 90/80"]
    categories = sorted({r["VFCCategory"] for r in hits})
    fig, ax = plt.subplots(figsize=(7.3, 5.2), constrained_layout=True)
    for index, label in enumerate(categories):
        selected = [r for r in hits if r["VFCCategory"] == label]
        ax.scatter([float(r["IdentityPercent"]) for r in selected], [float(r["ReferenceCoveragePercent"]) for r in selected], s=30, alpha=0.75, color=COLORS[index % len(COLORS)], label=label)
    ax.axvline(90, color="#666666", linestyle="--", linewidth=0.8); ax.axhline(80, color="#666666", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Nucleotide identity (%)"); ax.set_ylabel("Reference coverage (%)"); ax.set_xlim(89, 101); ax.set_ylim(79, 101); ax.grid(color="#E6E6E6"); ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save(fig, figure_dir / "39-hit-quality")

    cats = read_tsv(frozen / "vfc-category-summary.tsv")[::-1]
    vfc_totals = {
        "Catalog genes": sum(float(r["Genes"]) for r in cats),
        "MOCK1 reads": sum(float(r["MOCK1RawReads"]) for r in cats),
        "MOCK2 reads": sum(float(r["MOCK2RawReads"]) for r in cats),
    }
    fig, ax = plt.subplots(figsize=(9.5, 5.6), constrained_layout=True)
    grouped_bar(ax, [r["VFCCategory"] for r in cats], [
        ("Catalog genes", [100 * float(r["Genes"]) / vfc_totals["Catalog genes"] for r in cats]),
        ("MOCK1 reads", [100 * float(r["MOCK1RawReads"]) / vfc_totals["MOCK1 reads"] for r in cats]),
        ("MOCK2 reads", [100 * float(r["MOCK2RawReads"]) / vfc_totals["MOCK2 reads"] for r in cats]),
    ], "Share within primary VFDB-matched mass (%)")
    save(fig, figure_dir / "39-vfc-category-profile")

    context = read_tsv(frozen / "context-summary.tsv")
    labels = [r["Branch"].replace(" core 90/80", "") for r in context]
    fig, ax = plt.subplots(figsize=(8.4, 4.6), constrained_layout=True)
    grouped_bar(ax, labels, [
        ("Hits", [int(r["Hits"]) for r in context]),
        ("Sequences with hits", [int(r["SequencesWithHits"]) for r in context]),
        ("VFC categories", [int(r["VFCategories"]) for r in context]),
    ], "Count", horizontal=False)
    save(fig, figure_dir / "39-context-controls")


def draw40(frozen: Path, figure_dir: Path) -> None:
    yields = read_tsv(frozen / "tool-yield-summary.tsv")
    display = {
        "salinispora-full": "S. tropica, complete",
        "salinispora-fragmented": "S. tropica, 20-kb fragments",
        "nostoc": "Nostoc PCC 7120",
        "coassembly-ge20kb": "Metagenome contigs >=20 kb",
    }
    labels = [display[row["Dataset"]] for row in yields]
    fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    grouped_bar(ax, labels, [
        ("antiSMASH", [int(row["antiSMASHRegions"]) for row in yields]),
        ("GECCO", [int(row["GECCORegions"]) for row in yields]),
    ], "Predicted BGC regions", horizontal=False)
    save(fig, figure_dir / "40-tool-bgc-yield")

    fragmentation = read_tsv(frozen / "fragmentation-sensitivity.tsv")
    recovery = ("Recovered >=80%", "Partial 20-<80%", "Missed <20%")
    recovery_display = ("Recovered ≥80%", "Partially recovered (20–<80%)", "Missed (<20%)")
    colors = ("#009E73", "#E69F00", "#D55E00")
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    bottom = np.zeros(2)
    for label, display_label, color in zip(recovery, recovery_display, colors):
        values = [sum(row["Tool"] == tool and row["RecoveryClass"] == label for row in fragmentation) for tool in ("antiSMASH", "GECCO")]
        bars = ax.bar(["antiSMASH", "GECCO"], values, bottom=bottom, color=color, label=display_label)
        for bar, value, base in zip(bars, values, bottom):
            if value:
                ax.text(bar.get_x() + bar.get_width() / 2, base + value / 2, str(value), ha="center", va="center", color="white", fontweight="bold")
        bottom += np.asarray(values)
    ax.set_ylabel("Complete-genome BGC regions"); ax.legend(frameon=False); ax.grid(axis="y", color="#E6E6E6")
    save(fig, figure_dir / "40-fragmentation-sensitivity")

    types = read_tsv(frozen / "bgc-type-summary.tsv")
    totals = {
        "Regions": sum(float(row["FractionalRegions"]) for row in types),
        "Genes": sum(float(row["FractionalGenes"]) for row in types),
        "MOCK1": sum(float(row["MOCK1FractionalRawReads"]) for row in types),
        "MOCK2": sum(float(row["MOCK2FractionalRawReads"]) for row in types),
    }
    labels = [row["Category"] for row in types][::-1]
    fig, ax = plt.subplots(figsize=(8.7, 5.0), constrained_layout=True)
    grouped_bar(ax, labels, [
        ("BGC regions", [100 * float(row["FractionalRegions"]) / totals["Regions"] for row in types][::-1]),
        ("BGC genes", [100 * float(row["FractionalGenes"]) / totals["Genes"] for row in types][::-1]),
        ("MOCK1 reads", [100 * float(row["MOCK1FractionalRawReads"]) / totals["MOCK1"] for row in types][::-1]),
        ("MOCK2 reads", [100 * float(row["MOCK2FractionalRawReads"]) / totals["MOCK2"] for row in types][::-1]),
    ], "Fractional share within antiSMASH BGC mass (%)")
    save(fig, figure_dir / "40-bgc-type-profile")

    regions = read_tsv(frozen / "coassembly-bgc-abundance.tsv")
    category_colors = dict(zip(CATEGORY_ORDER := ("PKS", "NRPS", "RiPP", "Terpene", "Saccharide", "Other"), COLORS))
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.5), sharex=True, sharey=True, constrained_layout=True)
    for ax, sample in zip(axes, ("MOCK1", "MOCK2")):
        for category in CATEGORY_ORDER:
            selected = [row for row in regions if row["Categories"].split(";")[0] == category and row["TopMIBiGRiQ"] != ""]
            if not selected:
                continue
            x = [float(row["TopMIBiGRiQ"]) for row in selected]
            y = [float(row[f"{sample}FractionalRawReads"]) + 1 for row in selected]
            size = [18 + 42 * min(float(row["LengthBp"]) / 100_000, 1) for row in selected]
            ax.scatter(x, y, s=size, alpha=0.75, color=category_colors[category], label=category)
        ax.set_title(sample); ax.set_xlabel("Best MIBiG 4.0 region-to-region score"); ax.set_yscale("log"); ax.grid(color="#E6E6E6")
    axes[0].set_ylabel("Assigned raw reads + 1 (log scale)")
    handles, legend_labels = axes[1].get_legend_handles_labels()
    axes[1].legend(handles, legend_labels, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save(fig, figure_dir / "40-mibig-similarity-abundance")


def audit_chapter(article: int, chapter: Path, figures: list[str], checks: Checks) -> list[dict[str, str]]:
    text = chapter.read_text(encoding="utf-8")
    expected = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "target": "## 这一步对应论文里的哪张图" in text,
        "theory": "## 理论" in text, "setup": "## 准备工作" in text,
        "code": "## 可复制代码" in text, "audit": "## 审计与升级" in text,
        "publication": "## 出版级美化" in text, "pitfalls": "## 常见坑" in text,
        "methods": "## 这段 Methods 怎么写" in text, "own-data": "## 换成你自己的数据怎么做" in text,
        "references": "## 参考" in text, "resource-contract": all(x in text for x in ("RAM", "磁盘", "CPU", "耗时")),
        "seed": f"202607{article}" in text, "methods-template": "Methods template" in text,
        "four-images": all(f"../figures/{stem}.png" in text for stem in figures),
        "no-placeholders": re.search(r"__[A-Z0-9_]+__", text) is None,
        "no-meta-prose": not any(x in text for x in ("本篇可独立", "本文可独立", "全系列约定", "接口只学一次", "作者代码通常长这样", "（即本文）", "无头服务器")),
    }
    rows = []
    for key, status in expected.items():
        checks.add("Chapter", key, status, key); rows.append({"Check": key, "Status": "PASS" if status else "FAIL"})
    return rows


def audit_images(article: int, figure_dir: Path, checks: Checks) -> list[dict[str, str]]:
    rows = []
    for stem in FIGURES[article]:
        for suffix in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            ok = path.is_file() and path.stat().st_size > 0; detail = "exists" if ok else "missing"
            if ok and suffix in {"png", "tiff"}:
                with Image.open(path) as image:
                    dpi = image.info.get("dpi", (0, 0)); width, height = image.size
                    ok = width >= 1400 and height >= 900 and min(float(dpi[0]), float(dpi[1])) >= 340
                    detail = f"{width}x{height}; dpi={dpi}"
            checks.add("Figure", f"{stem}-{suffix}", ok, detail); rows.append({"Figure": path.name, "Status": "PASS" if ok else "FAIL", "Detail": detail})
    return rows


def scientific_checks(article: int, frozen: Path, checks: Checks) -> dict[str, object]:
    summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    checks.add("Identity", "article", summary["article"] == article, summary["article"])
    if article == 40:
        calls = read_tsv(frozen / "bgc-region-calls.tsv.gz")
        yields = read_tsv(frozen / "tool-yield-summary.tsv")
        fragmentation = read_tsv(frozen / "fragmentation-sensitivity.tsv")
        types = read_tsv(frozen / "bgc-type-summary.tsv")
        regions = read_tsv(frozen / "coassembly-bgc-abundance.tsv")
        for dataset in ("salinispora-full", "salinispora-fragmented", "nostoc", "coassembly-ge20kb"):
            for tool, key in (("antiSMASH", "antismash_regions"), ("GECCO", "gecco_regions")):
                observed = sum(row["Dataset"] == dataset and row["Tool"] == tool for row in calls)
                checks.add("Ledger", f"{dataset}-{tool}-count", observed == summary[key][dataset], observed)
        expected_antismash = {"salinispora-full": 19, "salinispora-fragmented": 29, "nostoc": 16, "coassembly-ge20kb": 71}
        expected_gecco = {"salinispora-full": 15, "salinispora-fragmented": 27, "nostoc": 9, "coassembly-ge20kb": 21}
        checks.add("Scientific", "antismash-counts", summary["antismash_regions"] == expected_antismash, summary["antismash_regions"])
        checks.add("Scientific", "gecco-counts", summary["gecco_regions"] == expected_gecco, summary["gecco_regions"])
        checks.add("Scientific", "coassembly-bgc-genes", summary["coassembly_antismash_genes"] == 1499, summary["coassembly_antismash_genes"])
        checks.add("Scientific", "coassembly-mock1-reads", abs(float(summary["coassembly_primary_raw_reads"]["MOCK1"]) - 69790) < 1e-8, summary["coassembly_primary_raw_reads"])
        checks.add("Scientific", "coassembly-mock2-reads", abs(float(summary["coassembly_primary_raw_reads"]["MOCK2"]) - 68816) < 1e-8, summary["coassembly_primary_raw_reads"])
        expected_fragmentation = {
            "antiSMASH": {"Recovered >=80%": 5, "Partial 20-<80%": 14, "Missed <20%": 0},
            "GECCO": {"Recovered >=80%": 3, "Partial 20-<80%": 10, "Missed <20%": 2},
        }
        checks.add("Scientific", "fragmentation-recovery", summary["salinispora_fragmentation"] == expected_fragmentation, summary["salinispora_fragmentation"])
        checks.add("Scientific", "salinosporamide-fragmentation", summary["salinosporamide_knowncluster_similarity_percent"] == {"complete_genome": 90, "twenty_kb_fragment": 35}, summary["salinosporamide_knowncluster_similarity_percent"])
        checks.add("Ledger", "yield-four-datasets", len(yields) == 4, len(yields))
        checks.add("Ledger", "fragmentation-rows-close", len(fragmentation) == summary["antismash_regions"]["salinispora-full"] + summary["gecco_regions"]["salinispora-full"], len(fragmentation))
        checks.add("Ledger", "type-regions-conserve", abs(sum(float(row["FractionalRegions"]) for row in types) - len(regions)) < 1e-5, sum(float(row["FractionalRegions"]) for row in types))
        checks.add("Ledger", "coassembly-region-count", len(regions) == summary["coassembly_antismash_regions"], len(regions))
        checks.add("Ledger", "catalog-membership-complete", summary["coassembly_genes_missing_catalog_membership"] == 0, summary["coassembly_genes_missing_catalog_membership"])
        checks.add("Boundary", "mibig-not-novelty", summary["mibig_similarity_is_not_compound_novelty"] is True, summary["mibig_similarity_is_not_compound_novelty"])
        checks.add("Boundary", "presence-not-production", summary["bgc_presence_is_not_expression_or_metabolite_production"] is True, summary["bgc_presence_is_not_expression_or_metabolite_production"])
        checks.add("Boundary", "fragmentation-aware", summary["fragmentation_changes_bgc_boundaries"] is True, summary["fragmentation_changes_bgc_boundaries"])
        checks.add("Ledger", "mock1-denominator", summary["sample_assigned_reads"]["MOCK1"] == 2_784_234, summary["sample_assigned_reads"])
        checks.add("Ledger", "mock2-denominator", summary["sample_assigned_reads"]["MOCK2"] == 2_777_443, summary["sample_assigned_reads"])
        return summary
    checks.add("Identity", "catalog-genes", summary["catalog_genes"] == 93_782, summary["catalog_genes"])
    expected = {
        37: {"primary_cazyme_genes": 2050, "single_tool_sensitivity_genes": 1849, "three_tool_consensus_genes": 1605, "btheta_cgcs": 117, "btheta_cgcs_with_substrate": 25},
        38: {"catalog_primary_genes": 36, "catalog_loose_only_genes": 6763, "coassembly_primary_orfs": 34, "pseudomonas_primary_orfs": 32, "staphylococcus_primary_orfs": 21},
        39: {"catalog_core_primary_genes": 93, "catalog_core_sensitive_genes": 123, "catalog_full_primary_genes": 184, "coassembly_core_primary_hits": 92, "pseudomonas_core_primary_hits": 199, "staphylococcus_core_primary_hits": 89},
    }[article]
    for key, value in expected.items(): checks.add("Scientific", key, summary[key] == value, summary[key])
    calls_name = {37: "cazyme-gene-calls.tsv.gz", 38: "resistome-gene-calls.tsv.gz", 39: "virulome-gene-calls.tsv.gz"}[article]
    calls = read_tsv(frozen / calls_name)
    checks.add("Ledger", "all-catalog-genes", len(calls) == 93_782 and len({r["GeneID"] for r in calls}) == 93_782, len(calls))
    if article == 37:
        evidence = read_tsv(frozen / "evidence-tier-summary.tsv")
        checks.add("Ledger", "evidence-closes", sum(int(r["Genes"]) for r in evidence) == 93_782, sum(int(r["Genes"]) for r in evidence))
        classes = read_tsv(frozen / "cazyme-class-summary.tsv")
        checks.add("Ledger", "class-fractional-conservation", abs(sum(float(r["FractionalGeneEquivalent"]) for r in classes) - 2050) < 1e-5, sum(float(r["FractionalGeneEquivalent"]) for r in classes))
        checks.add("Boundary", "substrate-not-activity", summary["substrate_predictions_are_not_activity"] is True, summary["substrate_predictions_are_not_activity"])
    elif article == 38:
        checks.add("Boundary", "loose-sensitivity-only", summary["loose_is_sensitivity_only"] is True, summary["loose_is_sensitivity_only"])
        checks.add("Boundary", "not-phenotype", summary["arg_presence_is_not_phenotypic_resistance"] is True, summary["arg_presence_is_not_phenotypic_resistance"])
        checks.add("Boundary", "no-nudge", summary["include_nudge"] is False, summary["include_nudge"])
    else:
        cats = read_tsv(frozen / "vfc-category-summary.tsv")
        checks.add("Ledger", "vfc-count-closes", sum(int(float(r["Genes"])) for r in cats) == 93, sum(int(float(r["Genes"])) for r in cats))
        checks.add("Boundary", "not-phenotype", summary["vfdb_match_is_not_virulence_phenotype"] is True, summary["vfdb_match_is_not_virulence_phenotype"])
    checks.add("Ledger", "mock1-denominator", summary["sample_assigned_reads"]["MOCK1"] == 2_784_234, summary["sample_assigned_reads"])
    checks.add("Ledger", "mock2-denominator", summary["sample_assigned_reads"]["MOCK2"] == 2_777_443, summary["sample_assigned_reads"])
    return summary


def main(article: int) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    args = parser.parse_args()
    frozen, output, figure_dir, chapter = args.frozen_dir.resolve(), args.output_dir.resolve(), args.figure_dir.resolve(), args.chapter.resolve()
    output.mkdir(parents=True, exist_ok=True); figure_dir.mkdir(parents=True, exist_ok=True); style(); checks = Checks()
    checksum_rows = checksum_audit(frozen, checks)
    summary = scientific_checks(article, frozen, checks)
    chapter_rows = audit_chapter(article, chapter, FIGURES[article], checks)
    {37: draw37, 38: draw38, 39: draw39, 40: draw40}[article](frozen, figure_dir)
    image_rows = audit_images(article, figure_dir, checks)
    write_tsv(output / "checksum-audit.tsv", checksum_rows)
    write_tsv(output / "chapter-audit.tsv", chapter_rows)
    write_tsv(output / "image-audit.tsv", image_rows)
    write_tsv(output / "validation-checks.tsv", checks.rows)
    result = {"article": article, "status": "passed" if checks.failed == 0 else "failed", "checks_passed": checks.passed, "checks_failed": checks.failed, **summary}
    (output / "validation-summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "validation.log").write_text("\n".join([f"Article {article} validation: {result['status'].upper()}", f"PASS={checks.passed}", f"FAIL={checks.failed}"] + [f"{r['Status']}\t{r['Category']}\t{r['CheckID']}\t{r['Detail']}" for r in checks.rows]) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if checks.failed == 0 else 1
