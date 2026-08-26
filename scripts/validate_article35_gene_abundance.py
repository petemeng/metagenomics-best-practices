#!/usr/bin/env python3
"""Validate Article 35 frozen evidence and draw four publication-ready figures."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import runpy
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, TextIO

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article35-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


SEED = 20260735
CATALOG_GENES = 93_782
EXPECTED_INPUT_READS = {"MOCK1": 3_999_706, "MOCK2": 3_999_776}
EXPECTED_PRIMARY_MAPPED = {"MOCK1": 3_346_575, "MOCK2": 3_335_721}
EXPECTED_MAIN_ASSIGNED = {"MOCK1": 2_784_234, "MOCK2": 2_777_443}
EXPECTED_STRICT_ASSIGNED = {"MOCK1": 2_586_966, "MOCK2": 2_572_767}
EXPECTED_LEGACY_ASSIGNED = {"MOCK1": 8_378, "MOCK2": 8_345}
EXPECTED_COMPLETENESS = {"Complete": 64_849, "Partial": 27_124, "Incomplete": 1_809}
EXPECTED_UNIREF90_HIT_GENES = 89_339
EXPECTED_REACTION_GENES = 15_392
EXPECTED_MULTI_REACTION_GENES = 4_488
EXPECTED_REACTION_READS = {"MOCK1": 596_658, "MOCK2": 594_503}
EXPECTED_REACTION_ROWS = 5_804
EXPECTED_DIAMOND_CANDIDATE_ROWS = 342_801
UNIREF90_SHA256 = "67a00a99ead2a00c737b4b9cb7e64ecc9085c2539bbd21a2d0c92913936995a8"
REACTION_MAP_SHA256 = "8419ce78a62ca9130914f2c347a9708111cedc7de52ba274659ce51ec7de7752"
ENV_YAML_SHA256 = "7a0280495b6712a6d590397e4876f13dd103112bfe14dba3f0687549f853874f"
ENV_LOCK_SHA256 = "91e495c6618400a829703bb2bea903d364de55f0a431f42f7cb5e134c2ef4c2d"
EXPECTED_LINEAGE_SHA256 = {
    "primary-gene-catalog-fna": "56f0be1fa7230517318dd745deba55be204da473ee3b6abbc24bd56ccaf3ceb6",
    "primary-gene-catalog-faa": "3db88ff78a548dddfc48caa8a17f04bcbb58dcafe2345d9dd31bc4e12f2a3569",
    "primary-gene-catalog-metadata": "677479da11ef41a6f27f11798c38dfd9b5830b5c564d9a8f436951413acd7c09",
    "MOCK1-R1": "d917a7241a29e2151bd1e2928994acd34b309d29ede5d6bd90933c8f31717148",
    "MOCK1-R2": "d224c825299ebdcc903b3e4665fbe13efada5e2855b18750a9cbc2ef661735ae",
    "MOCK2-R1": "a1908f80858cc44f19ba50a786c6c81fe85106de1073d02f4eaa1ed961379cd1",
    "MOCK2-R2": "6e945df7b952a849bf6b23e49afe175db8bca9bb07995bf596b37344abbc062d",
    "uniref90-diamond-db": UNIREF90_SHA256,
    "uniref90-metacyc-reaction-map": REACTION_MAP_SHA256,
    "delgado-paper-xml": "bb784d69a8ab260b90abe67d31a192a7f7bef951a98f2bed6ce3dae9092ce9ba",
}
FIGURE_STEMS = (
    "35-read-to-gene-ledger",
    "35-mapping-policy-sensitivity",
    "35-unit-normalization-audit",
    "35-functional-aggregation",
)
FROZEN_SCRIPT_NAMES = (
    "download_article35_gene_abundance_sources.sh",
    "prepare_article35_gene_abundance_inputs.py",
    "parse_article35_sam.py",
    "run_article35_gene_abundance.py",
    "summarize_article35_gene_abundance.py",
    "freeze_article35_gene_abundance.py",
    "validate_article35_gene_abundance.py",
)
SAMPLE_COLORS = {"MOCK1": "#0072B2", "MOCK2": "#D55E00"}
LEDGER_COLORS = {
    "Assigned": "#009E73",
    "Filtered mapped": "#E69F00",
    "Unmapped": "#BDBDBD",
}
LENGTH_COLORS = {
    "<300 bp": "#56B4E9",
    "300-599 bp": "#0072B2",
    "600-899 bp": "#009E73",
    "900-1,499 bp": "#E69F00",
    ">=1,500 bp": "#CC79A7",
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
    for number, line in enumerate((frozen / "file-checksums.sha256").read_text(encoding="utf-8").splitlines(), start=1):
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
            "axes.titlesize": 11,
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


def figure_ledger(mapping: list[dict[str, str]], legacy: list[dict[str, str]], figure_dir: Path) -> None:
    set_pub_style()
    main = {row["Sample"]: row for row in mapping if row["Policy"] == "Main"}
    legacy_by_sample = {row["Sample"]: row for row in legacy}
    samples = ("MOCK1", "MOCK2")
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    y = np.arange(len(samples))
    left = np.zeros(len(samples))
    for component in ("Assigned", "Filtered mapped", "Unmapped"):
        values = []
        for sample in samples:
            row = main[sample]
            total = int(row["InputReads"])
            if component == "Assigned":
                value = int(row["AssignedReads"])
            elif component == "Filtered mapped":
                value = int(row["FilteredMappedReads"])
            else:
                value = int(row["UnmappedReads"])
            values.append(100 * value / total)
        ax.barh(y, values, left=left, color=LEDGER_COLORS[component], label=component, height=0.55)
        left += np.asarray(values)
    for index, sample in enumerate(samples):
        row = main[sample]
        raw_rate = 100 * int(row["PrimaryMappedReads"]) / int(row["InputReads"])
        assigned_rate = float(row["AssignedPctInput"])
        legacy_rate = float(legacy_by_sample[sample]["AssignedPct"])
        ax.vlines(raw_rate, index - 0.28, index + 0.28, color="#37474F", linewidth=1.2)
        ax.text(assigned_rate / 2, index, f"{assigned_rate:.1f}%", ha="center", va="center", color="white", fontweight="bold")
        ax.scatter([legacy_rate], [index + 0.33], marker="D", s=45, color="#332288", zorder=4)
    ax.plot([], [], color="#37474F", linewidth=1.2, label="Raw primary mapping")
    ax.scatter([], [], marker="D", s=45, color="#332288", label="10k-read historical branch")
    ax.set_yticks(y, samples)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Input reads (%)")
    ax.set_title("Raw mapping and assignable abundance are different denominators", fontweight="bold")
    ax.grid(axis="x", color="#ECEFF1", linewidth=0.8)
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[0])


def figure_policy(mapping: list[dict[str, str]], figure_dir: Path) -> None:
    set_pub_style()
    order = ("AllPrimary", "IdentityQcov", "Main", "Strict")
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    for sample in ("MOCK1", "MOCK2"):
        rows = {row["Policy"]: row for row in mapping if row["Sample"] == sample}
        x = [float(rows[policy]["AssignedPctInput"]) for policy in order]
        y = [int(rows[policy]["DetectedGenes"]) / 1000 for policy in order]
        ax.plot(x, y, marker="o", linewidth=2, markersize=7, color=SAMPLE_COLORS[sample], label=sample)
    label_offsets = {
        "AllPrimary": (8, 4),
        "IdentityQcov": (8, 10),
        "Main": (-8, -19),
        "Strict": (8, 6),
    }
    for policy in order:
        rows = [row for row in mapping if row["Policy"] == policy]
        x_value = sum(float(row["AssignedPctInput"]) for row in rows) / len(rows)
        y_value = sum(int(row["DetectedGenes"]) / 1000 for row in rows) / len(rows)
        offset = label_offsets[policy]
        ax.annotate(
            policy.replace("IdentityQcov", "Identity + qcov"),
            (x_value, y_value),
            xytext=offset,
            textcoords="offset points",
            fontsize=8,
            color="#37474F",
            ha="right" if offset[0] < 0 else "left",
        )
    ax.set_xlabel("Assigned reads / input reads (%)")
    ax.set_ylabel("Detected catalog genes (thousands)")
    ax.set_title("Alignment policy changes both denominator and detected feature space", fontweight="bold")
    ax.grid(color="#ECEFF1", linewidth=0.8)
    ax.legend()
    fig.tight_layout()
    save_figure(fig, figure_dir, FIGURE_STEMS[1])


def figure_units(normalization: list[dict[str, str]], lengths: list[dict[str, str]], figure_dir: Path) -> None:
    set_pub_style()
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.2))
    units = ("CPM", "RPKM", "TPM")
    x = np.arange(len(units))
    width = 0.34
    by_sample = {row["Sample"]: row for row in normalization}
    for offset, sample in ((-width / 2, "MOCK1"), (width / 2, "MOCK2")):
        values = [float(by_sample[sample][f"{unit}Sum"]) / 1_000_000 for unit in units]
        axes[0].bar(x + offset, values, width, color=SAMPLE_COLORS[sample], label=sample)
    axes[0].axhline(1, color="#455A64", linestyle="--", linewidth=1)
    axes[0].set_xticks(x, units)
    axes[0].set_ylabel("Column sum (millions)")
    axes[0].set_title("A  Only CPM and TPM have a 1M closure")
    axes[0].grid(axis="y", color="#ECEFF1")
    axes[0].legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.13))

    bins = ("<300 bp", "300-599 bp", "600-899 bp", "900-1,499 bp", ">=1,500 bp")
    keys = {(row["Sample"], row["LengthBin"]): row for row in lengths}
    labels = ("M1 Count", "M1 TPM", "M2 Count", "M2 TPM")
    bottoms = np.zeros(len(labels))
    for bin_name in bins:
        values: list[float] = []
        for sample in ("MOCK1", "MOCK2"):
            row = keys[(sample, bin_name)]
            sample_total_count = sum(float(keys[(sample, item)]["AssignedReads"]) for item in bins)
            sample_total_tpm = sum(float(keys[(sample, item)]["TPM"]) for item in bins)
            values.extend((100 * float(row["AssignedReads"]) / sample_total_count, 100 * float(row["TPM"]) / sample_total_tpm))
        axes[1].bar(np.arange(4), values, bottom=bottoms, color=LENGTH_COLORS[bin_name], label=bin_name)
        bottoms += np.asarray(values)
    axes[1].set_xticks(np.arange(4), labels, rotation=15)
    axes[1].set_ylim(0, 100)
    axes[1].set_ylabel("Abundance share (%)")
    axes[1].set_title("B  Length correction changes composition")
    axes[1].legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.19), fontsize=8)
    axes[1].grid(axis="y", color="#ECEFF1")
    fig.suptitle("Normalization units answer different quantitative questions", y=1.02, fontweight="bold")
    fig.subplots_adjust(bottom=0.25, wspace=0.30)
    save_figure(fig, figure_dir, FIGURE_STEMS[2])


def figure_functions(audit: list[dict[str, str]], figure_dir: Path) -> None:
    set_pub_style()
    by_sample = {row["Sample"]: row for row in audit}
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.0))
    categories = ("Genes\nUniRef90 hit", "Genes\nwith reaction", "Assigned reads\nUniRef90 hit", "Assigned reads\nwith reaction")
    x = np.arange(len(categories))
    width = 0.34
    for offset, sample in ((-width / 2, "MOCK1"), (width / 2, "MOCK2")):
        row = by_sample[sample]
        values = (
            100 * int(row["GenesWithUniRef90Hit"]) / int(row["CatalogGenes"]),
            100 * int(row["GenesWithReaction"]) / int(row["CatalogGenes"]),
            100 * int(row["AssignedReadsWithUniRef90Hit"]) / int(row["AssignedReads"]),
            float(row["AssignedReadPctWithReaction"]),
        )
        axes[0].bar(x + offset, values, width, color=SAMPLE_COLORS[sample], label=sample)
    axes[0].set_xticks(x, categories)
    axes[0].set_ylabel("Coverage (%)")
    axes[0].set_ylim(0, 100)
    axes[0].set_title("A  Annotation coverage needs two denominators")
    axes[0].legend()
    axes[0].grid(axis="y", color="#ECEFF1")

    units = ("CPM", "TPM")
    x2 = np.arange(2)
    for offset, sample in ((-width / 2, "MOCK1"), (width / 2, "MOCK2")):
        row = by_sample[sample]
        values = (float(row["CopyCPMInflation"]), float(row["CopyTPMInflation"]))
        bars = axes[1].bar(x2 + offset, values, width, color=SAMPLE_COLORS[sample], label=sample)
        axes[1].bar_label(bars, labels=[f"{value:.3f}×" for value in values], padding=3, fontsize=8)
    axes[1].axhline(1, color="#455A64", linestyle="--", linewidth=1, label="Mass-conserving split")
    axes[1].set_xticks(x2, units)
    axes[1].set_ylabel("Total after copy / total after split")
    max_inflation = max(
        max(float(row["CopyCPMInflation"]), float(row["CopyTPMInflation"]))
        for row in audit
    )
    axes[1].set_ylim(0, max_inflation * 1.15)
    axes[1].set_title("B  Copying multi-mapped functions inflates mass")
    axes[1].grid(axis="y", color="#ECEFF1")
    axes[1].legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=8)
    fig.suptitle("Functional aggregation requires an explicit one-to-many policy", y=1.02, fontweight="bold")
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.22, wspace=0.25)
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
                preview = image.convert("RGB")
                preview.thumbnail((256, 256))
                colors = preview.getcolors(maxcolors=65_536)
                color_count = len(colors) if colors else 65_537
            with Image.open(tiff) as image:
                tiff_dpi = image.info.get("dpi", (0, 0))
                tiff_compression = image.info.get("compression", "")
            status = (
                width >= 1600
                and height >= 900
                and color_count > 20
                and pdf.stat().st_size > 5_000
                and tiff.stat().st_size > 10_000
                and all(abs(float(value) - 350) < 1 for value in tiff_dpi)
                and tiff_compression == "tiff_lzw"
            )
            detail = (
                f"{width}x{height}; colors={color_count}; pdf={pdf.stat().st_size}; "
                f"tiff={tiff.stat().st_size}; dpi={tiff_dpi}; compression={tiff_compression}"
            )
        checks.add("Images", stem, status, detail)
        rows.append({"Figure": stem, "Status": "PASS" if status else "FAIL", "Detail": detail})
    return rows


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = (args.frozen_dir or root / "data/small/35-gene-abundance-frozen").resolve()
    output = (args.output_dir or root / "results/35-gene-abundance").resolve()
    figure_dir = (args.figure_dir or root / "figures").resolve()
    chapter = (args.chapter or root / "chapters/35-gene-abundance.qmd").resolve()
    output.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    checksum_rows = verify_checksum_manifest(frozen, checks)
    write_tsv(output / "checksum-audit.tsv", checksum_rows, ["File", "ExpectedSHA256", "ObservedSHA256", "Status"])

    contract = json.loads((frozen / "frozen-contract.json").read_text(encoding="utf-8"))
    summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    checks.add("Contract", "article", contract.get("article") == 35 == summary.get("article"), contract.get("article"))
    checks.add("Contract", "seed", contract.get("seed") == SEED == summary.get("seed"), contract.get("seed"))
    checks.add("Contract", "catalog-size", contract.get("catalog_genes") == CATALOG_GENES == summary.get("catalog_genes"), contract.get("catalog_genes"))
    checks.add("Contract", "uniref-release", contract.get("uniref90_release") == "v201901b", contract.get("uniref90_release"))
    checks.add("Contract", "uniref-sha", contract.get("uniref90_db_sha256") == UNIREF90_SHA256, contract.get("uniref90_db_sha256"))
    checks.add("Contract", "reaction-sha", contract.get("reaction_map_sha256") == REACTION_MAP_SHA256, contract.get("reaction_map_sha256"))
    checks.add("Contract", "read-unit", contract.get("read_unit") == "R1 and R2 mapped independently as reads", contract.get("read_unit"))
    checks.add("Contract", "large-assets-excluded", set(contract.get("large_assets_excluded", [])) == {"FASTQ", "Bowtie2 index", "SAM/BAM", "UniRef90 DIAMOND database"}, contract.get("large_assets_excluded"))
    primary_policy = contract.get("primary_mapping_policy", {})
    checks.add("Contract", "identity-formula", primary_policy.get("identity_formula") == "1 - NM / sum(CIGAR M,I,D,=,X alignment columns)", primary_policy.get("identity_formula"))
    checks.add("Contract", "qcov-formula", primary_policy.get("query_coverage_formula") == "sum(CIGAR M,I,=,X query-aligned bases) / full query length", primary_policy.get("query_coverage_formula"))
    checks.add("Contract", "main-filters", primary_policy.get("filters") == {"minimum_mapq": 10, "minimum_identity": 0.95, "minimum_query_coverage": 0.8}, primary_policy.get("filters"))
    functional_policy = contract.get("functional_policy", {})
    checks.add("Contract", "diamond-iterate", functional_policy.get("iterate") == ["faster", "sensitive"], functional_policy.get("iterate"))
    checks.add("Contract", "diamond-top5", functional_policy.get("max_target_seqs") == 5, functional_policy.get("max_target_seqs"))
    checks.add("Contract", "diamond-masking", functional_policy.get("masking") == 1, functional_policy.get("masking"))
    checks.add("Contract", "diamond-thresholds", functional_policy.get("thresholds") == {"minimum_identity_pct": 50, "minimum_query_coverage_pct": 80, "maximum_evalue": 1e-5}, functional_policy.get("thresholds"))
    checks.add("Contract", "best-hit-order", functional_policy.get("best_hit_order") == "bitscore desc, evalue asc, identity desc, UniRef90 ID asc", functional_policy.get("best_hit_order"))

    lineage = read_tsv(frozen / "input-lineage.tsv")
    checks.add("Inputs", "lineage-rows", len(lineage) == 10, len(lineage))
    checks.add("Inputs", "lineage-unique", len({row["Asset"] for row in lineage}) == 10, len({row["Asset"] for row in lineage}))
    checks.add("Inputs", "two-mocks-four-fastqs", sum(row["Asset"].startswith(("MOCK1", "MOCK2")) for row in lineage) == 4, [row["Asset"] for row in lineage])
    checks.add("Inputs", "all-sha256", all(re.fullmatch(r"[0-9a-f]{64}", row["SHA256"]) for row in lineage), "10 source digests")
    observed_lineage = {row["Asset"]: row["SHA256"] for row in lineage}
    checks.add("Inputs", "locked-lineage-digests", observed_lineage == EXPECTED_LINEAGE_SHA256, sorted(set(observed_lineage) ^ set(EXPECTED_LINEAGE_SHA256)))
    checks.add("Inputs", "env-yaml-sha", sha256(frozen / "env/gene-abundance.yml") == ENV_YAML_SHA256, sha256(frozen / "env/gene-abundance.yml"))
    checks.add("Inputs", "env-lock-sha", sha256(frozen / "env/gene-abundance-linux-64.lock") == ENV_LOCK_SHA256, sha256(frozen / "env/gene-abundance-linux-64.lock"))
    lock_lines = (frozen / "env/gene-abundance-linux-64.lock").read_text(encoding="utf-8").splitlines()
    checks.add("Inputs", "env-lock-explicit-header", "@EXPLICIT" in lock_lines[:6], lock_lines[:6])
    checks.add("Inputs", "env-lock-package-urls", sum(line.startswith("https://") for line in lock_lines) >= 100, sum(line.startswith("https://") for line in lock_lines))
    checks.add("Inputs", "env-yaml-current", sha256(frozen / "env/gene-abundance.yml") == sha256(root / "env/gene-abundance.yml"), sha256(root / "env/gene-abundance.yml"))
    checks.add("Inputs", "env-lock-current", sha256(frozen / "env/gene-abundance-linux-64.lock") == sha256(root / "env/gene-abundance-linux-64.lock"), sha256(root / "env/gene-abundance-linux-64.lock"))
    for name in FROZEN_SCRIPT_NAMES:
        checks.add("Inputs", f"script-current-{name}", sha256(frozen / "scripts" / name) == sha256(root / "scripts" / name), name)

    tools = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    expected_tools = {"Bowtie2": "2.5.5", "SAMtools": "1.23.1", "HTSeq": "2.1.2", "DIAMOND": "2.2.4", "Python": "3.10.20"}
    for tool, version in expected_tools.items():
        checks.add("Tools", tool.lower(), tools.get(tool) == version, tools.get(tool))
    checks.add("Tools", "seqtk", tools.get("seqtk", "").startswith("1.5"), tools.get("seqtk"))

    command_log = read_tsv(frozen / "command-log.tsv")
    expected_steps = {"build-bowtie2-index", "annotate-uniref90"}
    for sample in ("MOCK1", "MOCK2"):
        expected_steps.update(
            {
                f"legacy-select-{sample}",
                f"legacy-map-{sample}",
                f"legacy-bam-{sample}",
                f"legacy-htseq-{sample}",
                f"audited-map-{sample}",
            }
        )
    commands = {row["Step"]: row["Command"] for row in command_log}
    checks.add("Commands", "exact-step-set", set(commands) == expected_steps, sorted(set(commands) ^ expected_steps))
    annotation_command = commands.get("annotate-uniref90", "")
    checks.add(
        "Commands",
        "diamond-contract",
        all(token in annotation_command for token in ("--iterate faster", "--sensitive", "--id 50", "--query-cover 80", "--evalue 1e-5", "--max-target-seqs 5", "--masking 1")),
        annotation_command,
    )
    annotation_log = (frozen / "logs/annotate-uniref90.stderr.log").read_text(encoding="utf-8", errors="replace")
    checks.add("Commands", "diamond-log-iterated", "Running iterated search mode with sensitivity steps: faster, sensitive" in annotation_log, "faster -> sensitive")
    checks.add("Commands", "diamond-log-database-size", "sequences: 87296736, letters: 29247941583" in annotation_log, "87,296,736 sequences")
    checks.add("Commands", "diamond-log-faster-pass", "Aligned 88370/93782 queries in this iteration" in annotation_log, "88,370/93,782")
    checks.add("Commands", "diamond-log-sensitive-pass", "Aligned 969/93782 queries in this iteration, 89339/93782 total." in annotation_log, "969 additional; 89,339 total")
    checks.add("Commands", "diamond-log-no-error", "Error:" not in annotation_log and "terminate called" not in annotation_log, "no DIAMOND error markers")
    for sample in ("MOCK1", "MOCK2"):
        command = commands.get(f"audited-map-{sample}", "")
        checks.add(
            "Commands",
            f"audited-map-{sample}",
            all(token in command for token in ("--very-sensitive-local", "-k 2", "--seed 20260735", "parse_article35_sam.py")),
            command,
        )

    parser_namespace = runpy.run_path(str(frozen / "scripts/parse_article35_sam.py"))
    alignment_metrics = parser_namespace["alignment_metrics"]
    _, _, deletion_identity, deletion_qcov = alignment_metrics("50M2D48M", "A" * 98, ["NM:i:3"])
    checks.add("Mapping", "identity-includes-deletion-columns", abs(deletion_identity - 0.97) < 1e-12, deletion_identity)
    checks.add("Mapping", "qcov-excludes-deletions", abs(deletion_qcov - 1.0) < 1e-12, deletion_qcov)
    _, _, clipped_identity, clipped_qcov = alignment_metrics("10S90M", "A" * 100, ["NM:i:0"])
    checks.add("Mapping", "soft-clipping-lowers-qcov", abs(clipped_identity - 1.0) < 1e-12 and abs(clipped_qcov - 0.9) < 1e-12, (clipped_identity, clipped_qcov))

    legacy = read_tsv(frozen / "legacy-mapping-audit.tsv")
    checks.add("Historical branch", "two-samples", len(legacy) == 2 and {row["Sample"] for row in legacy} == {"MOCK1", "MOCK2"}, len(legacy))
    for row in legacy:
        total = int(row["AssignedReads"]) + sum(int(row[field]) for field in ("NotAligned", "NoFeature", "AmbiguousFeature", "TooLowMAPQ", "AlignmentNotUnique"))
        checks.add("Historical branch", f"ledger-{row['Sample']}", total == 10_000, total)
        checks.add("Historical branch", f"locked-assigned-{row['Sample']}", int(row["AssignedReads"]) == EXPECTED_LEGACY_ASSIGNED[row["Sample"]], row["AssignedReads"])
        checks.add("Historical branch", f"mapq-zero-{row['Sample']}", row["MAPQThreshold"] == "0", row["MAPQThreshold"])

    mapping = read_tsv(frozen / "mapping-policy-summary.tsv")
    checks.add("Mapping", "eight-policy-rows", len(mapping) == 8, len(mapping))
    checks.add("Mapping", "four-policies", Counter(row["Policy"] for row in mapping) == {"AllPrimary": 2, "IdentityQcov": 2, "Main": 2, "Strict": 2}, Counter(row["Policy"] for row in mapping))
    mapping_by_key = {(row["Sample"], row["Policy"]): row for row in mapping}
    for sample in ("MOCK1", "MOCK2"):
        rows = [mapping_by_key[(sample, policy)] for policy in ("AllPrimary", "IdentityQcov", "Main", "Strict")]
        input_counts = {int(row["InputReads"]) for row in rows}
        primary_counts = {int(row["PrimaryMappedReads"]) for row in rows}
        unmapped_counts = {int(row["UnmappedReads"]) for row in rows}
        checks.add("Mapping", f"stable-denominators-{sample}", len(input_counts) == len(primary_counts) == len(unmapped_counts) == 1, (input_counts, primary_counts, unmapped_counts))
        main = mapping_by_key[(sample, "Main")]
        strict = mapping_by_key[(sample, "Strict")]
        checks.add("Mapping", f"ledger-{sample}", int(main["AssignedReads"]) + int(main["FilteredMappedReads"]) + int(main["UnmappedReads"]) == int(main["InputReads"]), main)
        checks.add("Mapping", f"locked-input-{sample}", int(main["InputReads"]) == EXPECTED_INPUT_READS[sample], main["InputReads"])
        checks.add("Mapping", f"locked-primary-{sample}", int(main["PrimaryMappedReads"]) == EXPECTED_PRIMARY_MAPPED[sample], main["PrimaryMappedReads"])
        checks.add("Mapping", f"locked-main-{sample}", int(main["AssignedReads"]) == EXPECTED_MAIN_ASSIGNED[sample], main["AssignedReads"])
        checks.add("Mapping", f"locked-strict-{sample}", int(strict["AssignedReads"]) == EXPECTED_STRICT_ASSIGNED[sample], strict["AssignedReads"])
        assigned = [int(row["AssignedReads"]) for row in rows]
        detected = [int(row["DetectedGenes"]) for row in rows]
        checks.add("Mapping", f"policy-monotone-reads-{sample}", assigned[0] >= assigned[1] >= assigned[2] >= assigned[3], assigned)
        checks.add("Mapping", f"policy-monotone-genes-{sample}", detected[0] >= detected[1] >= detected[2] >= detected[3], detected)
        checks.add("Mapping", f"main-rate-{sample}", 65 < float(main["AssignedPctInput"]) < 75, main["AssignedPctInput"])
        checks.add("Mapping", f"not-raw-map-rate-{sample}", int(main["AssignedReads"]) < int(main["PrimaryMappedReads"]), f"{main['AssignedReads']}/{main['PrimaryMappedReads']}")

    normalization = read_tsv(frozen / "normalization-audit.tsv")
    checks.add("Units", "two-samples", len(normalization) == 2, len(normalization))
    for row in normalization:
        sample = row["Sample"]
        main = mapping_by_key[(sample, "Main")]
        checks.add("Units", f"raw-count-closure-{sample}", int(row["RawCountSum"]) == int(row["AssignedReads"]) == int(main["AssignedReads"]), row["RawCountSum"])
        checks.add("Units", f"cpm-closure-{sample}", float(row["CPMClosureError"]) < 1e-6, row["CPMClosureError"])
        checks.add("Units", f"tpm-closure-{sample}", float(row["TPMClosureError"]) < 1e-6, row["TPMClosureError"])
        checks.add("Units", f"rpkm-not-closed-{sample}", abs(float(row["RPKMSum"]) - 1_000_000) > 1_000, row["RPKMSum"])

    abundance = read_tsv(frozen / "gene-abundance-long.tsv.gz")
    checks.add("Gene table", "row-count", len(abundance) == 2 * CATALOG_GENES, len(abundance))
    checks.add("Gene table", "sample-rows", Counter(row["Sample"] for row in abundance) == {"MOCK1": CATALOG_GENES, "MOCK2": CATALOG_GENES}, Counter(row["Sample"] for row in abundance))
    checks.add("Gene table", "unique-keys", len({(row["Sample"], row["GeneID"]) for row in abundance}) == len(abundance), len(abundance))
    checks.add("Gene table", "positive-lengths", all(int(row["NtLength"]) > 0 for row in abundance), min(int(row["NtLength"]) for row in abundance))
    for sample in ("MOCK1", "MOCK2"):
        completeness = Counter(row["Completeness"] for row in abundance if row["Sample"] == sample)
        checks.add("Gene table", f"completeness-lineage-{sample}", completeness == EXPECTED_COMPLETENESS, completeness)

    completeness_summary = read_tsv(frozen / "gene-completeness-summary.tsv")
    checks.add("Gene table", "completeness-summary-rows", len(completeness_summary) == 6, len(completeness_summary))
    for sample in ("MOCK1", "MOCK2"):
        rows = [row for row in completeness_summary if row["Sample"] == sample]
        checks.add("Gene table", f"completeness-summary-categories-{sample}", {row["Completeness"] for row in rows} == set(EXPECTED_COMPLETENESS), [row["Completeness"] for row in rows])
        checks.add("Gene table", f"completeness-cpm-closure-{sample}", abs(sum(float(row["CPM"]) for row in rows) - 1_000_000) < 1e-6, sum(float(row["CPM"]) for row in rows))
        checks.add("Gene table", f"completeness-tpm-closure-{sample}", abs(sum(float(row["TPM"]) for row in rows) - 1_000_000) < 1e-6, sum(float(row["TPM"]) for row in rows))

    annotations = read_tsv(frozen / "gene-functional-annotation.tsv.gz")
    checks.add("Functions", "annotation-row-count", len(annotations) == CATALOG_GENES, len(annotations))
    checks.add("Functions", "annotation-unique-genes", len({row["GeneID"] for row in annotations}) == CATALOG_GENES, len({row["GeneID"] for row in annotations}))
    hit_genes = sum(bool(row["UniRef90"]) for row in annotations)
    checks.add("Functions", "hit-count-summary", hit_genes == int(summary["uniref90_hit_genes"]), hit_genes)
    checks.add("Functions", "locked-uniref90-hit-count", hit_genes == EXPECTED_UNIREF90_HIT_GENES, hit_genes)
    checks.add("Functions", "uniref-id-format", all(not row["UniRef90"] or re.fullmatch(r"UniRef90_[A-Za-z0-9]+", row["UniRef90"]) for row in annotations), "UniRef90 IDs")
    checks.add(
        "Functions",
        "reaction-count-consistency",
        all(int(row["ReactionCount"]) == (len(row["Reactions"].split(";")) if row["Reactions"] else 0) for row in annotations),
        "ReactionCount equals parsed IDs",
    )
    checks.add("Functions", "hit-thresholds", all(not row["UniRef90"] or (float(row["Pident"]) >= 50 and float(row["QueryCoveragePct"]) >= 80 and float(row["Evalue"]) <= 1e-5) for row in annotations), "50% identity; 80% qcov; E<=1e-5")

    search_audit = read_tsv(frozen / "annotation-search-audit.tsv")
    checks.add("Functions", "annotation-search-audit-row", len(search_audit) == 1, len(search_audit))
    search_row = search_audit[0]
    checks.add("Functions", "locked-candidate-row-count", int(search_row["CandidateRows"]) == EXPECTED_DIAMOND_CANDIDATE_ROWS, search_row["CandidateRows"])
    checks.add("Functions", "candidate-query-count", int(search_row["QueriesWithCandidate"]) == EXPECTED_UNIREF90_HIT_GENES, search_row["QueriesWithCandidate"])
    checks.add("Functions", "candidate-cap-observed", int(search_row["MaxCandidatesPerQuery"]) == 5, search_row["MaxCandidatesPerQuery"])
    checks.add("Functions", "search-best-hit-count", int(search_row["BestHitGenes"]) == hit_genes, search_row["BestHitGenes"])
    checks.add("Functions", "search-threshold-audit", float(search_row["MinimumIdentityPct"]) == 50 and float(search_row["MinimumQueryCoveragePct"]) == 80 and float(search_row["MaximumEvalue"]) == 1e-5, search_row)

    functional = read_tsv(frozen / "functional-aggregation-audit.tsv")
    checks.add("Functions", "two-functional-rows", len(functional) == 2, len(functional))
    checks.add("Functions", "stable-reaction-gene-space", len({row["GenesWithReaction"] for row in functional}) == 1, [row["GenesWithReaction"] for row in functional])
    for row in functional:
        sample = row["Sample"]
        checks.add("Functions", f"locked-reaction-genes-{sample}", int(row["GenesWithReaction"]) == EXPECTED_REACTION_GENES, row["GenesWithReaction"])
        checks.add("Functions", f"locked-multi-reaction-genes-{sample}", int(row["MultiReactionGenes"]) == EXPECTED_MULTI_REACTION_GENES, row["MultiReactionGenes"])
        checks.add("Functions", f"locked-reaction-reads-{sample}", int(row["AssignedReadsWithReaction"]) == EXPECTED_REACTION_READS[sample], row["AssignedReadsWithReaction"])
        checks.add("Functions", f"split-reads-conserves-{sample}", abs(float(row["SplitReadEquivalentSum"]) - float(row["AssignedReads"])) < 1e-6, row["SplitReadEquivalentSum"])
        checks.add("Functions", f"copy-reads-inflate-{sample}", float(row["CopyReadInflation"]) > 1, row["CopyReadInflation"])
        checks.add("Functions", f"copy-read-cpm-agree-{sample}", abs(float(row["CopyReadInflation"]) - float(row["CopyCPMInflation"])) < 1e-12, (row["CopyReadInflation"], row["CopyCPMInflation"]))
        checks.add("Functions", f"split-cpm-conserves-{sample}", abs(float(row["SplitCPMSum"]) - float(row["GeneCPMSum"])) < 1e-6, row["SplitCPMSum"])
        checks.add("Functions", f"split-tpm-conserves-{sample}", abs(float(row["SplitTPMSum"]) - float(row["GeneTPMSum"])) < 1e-6, row["SplitTPMSum"])
        checks.add("Functions", f"copy-cpm-inflates-{sample}", float(row["CopyCPMInflation"]) > 1, row["CopyCPMInflation"])
        checks.add("Functions", f"copy-tpm-inflates-{sample}", float(row["CopyTPMInflation"]) > 1, row["CopyTPMInflation"])
        checks.add("Functions", f"reaction-read-bound-{sample}", 0 < float(row["AssignedReadPctWithReaction"]) < 100, row["AssignedReadPctWithReaction"])

    reactions = read_tsv(frozen / "reaction-abundance-long.tsv")
    checks.add("Functions", "locked-reaction-row-count", len(reactions) == EXPECTED_REACTION_ROWS, len(reactions))
    checks.add("Functions", "reaction-table-unique-keys", len({(row["Sample"], row["Reaction"]) for row in reactions}) == len(reactions), len(reactions))
    reaction_names = {row["Reaction"] for row in reactions}
    checks.add("Functions", "unannotated-buckets", {"NO_UNIREF90_HIT", "UNIREF90_NO_REACTION"} <= reaction_names, sorted({"NO_UNIREF90_HIT", "UNIREF90_NO_REACTION"} - reaction_names))

    resources = read_tsv(frozen / "resource-usage.tsv")
    checks.add("Resources", "recorded-steps", len(resources) == 12, len(resources))
    checks.add("Resources", "all-exit-zero", all(row["ExitStatus"] == "0" for row in resources), sorted({row["ExitStatus"] for row in resources}))
    checks.add("Resources", "rss-recorded", all(float(row["PeakRSSKiB"]) > 0 for row in resources), min(float(row["PeakRSSKiB"]) for row in resources))

    boundaries = summary.get("boundaries", {})
    checks.add("Boundaries", "seven-false-claims", len(boundaries) == 7 and all(value is False for value in boundaries.values()), boundaries)

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
    checks.add("Chapter", "no-placeholders", "ARTICLE35_" not in chapter_text, sorted(set(re.findall(r"ARTICLE35_[A-Z0-9_]+", chapter_text))))
    checks.add("Chapter", "versions", all(token in chapter_text for token in ("Bowtie2 2.5.5", "SAMtools 1.23.1", "HTSeq 2.1.2", "seqtk 1.5", "DIAMOND 2.2.4")), "tool versions")
    checks.add("Chapter", "upstream-input-versions", all(token in chapter_text for token in ("fastp 1.3.6", "Prodigal 2.6.3", "MMseqs2 9.d36de")), "input-generating tool versions")
    checks.add("Chapter", "mapping-parameters", all(token in chapter_text for token in ("--very-sensitive-local", "-k 2", "MAPQ >= 10", "identity >= 95%", "query coverage >= 80%")), "mapping policy")
    checks.add("Chapter", "annotation-parameters", all(token in chapter_text for token in ("--iterate faster", "--sensitive", "--id 50", "--query-cover 80", "--evalue 1e-5", "--max-target-seqs 5")), "DIAMOND policy")
    checks.add("Chapter", "seed", str(SEED) in chapter_text and "-s100" in chapter_text, SEED)
    checks.add("Chapter", "locked-counts", all(f"{value:,}" in chapter_text for value in summary["main_assigned_reads"].values()), summary["main_assigned_reads"])
    checks.add("Chapter", "hardware", all(token in chapter_text for token in ("RAM", "磁盘", "CPU", "耗时")), "resource labels")
    checks.add("Chapter", "inline-theme", all(token in chapter_text for token in ("install.packages", "pal_pub", "scale_color_pub", "scale_fill_pub", "theme_pub", "save_pub")), "inline preparation")
    checks.add("Chapter", "unit-boundaries", all(token in chapter_text for token in ("count-based", "absolute abundance", "NO_UNIREF90_HIT", "UNIREF90_NO_REACTION")), "unit and missingness boundaries")
    prohibited = ("本篇可独立跑通", "这体现全系列", "接口只学一次", "作者代码通常长这样", "（即本文）")
    checks.add("Chapter", "no-meta-prose", not any(text in chapter_text for text in prohibited), [text for text in prohibited if text in chapter_text])
    for stem in FIGURE_STEMS:
        checks.add("Chapter", f"figure-reference-{stem}", f"../figures/{stem}.png" in chapter_text, stem)

    lengths = read_tsv(frozen / "gene-length-bin-summary.tsv")
    figure_ledger(mapping, legacy, figure_dir)
    figure_policy(mapping, figure_dir)
    figure_units(normalization, lengths, figure_dir)
    figure_functions(functional, figure_dir)
    image_rows = image_audit(figure_dir, checks)
    write_tsv(output / "image-audit.tsv", image_rows, ["Figure", "Status", "Detail"])

    write_tsv(output / "mapping-audit.tsv", mapping, list(mapping[0]))
    write_tsv(output / "functional-audit.tsv", functional, list(functional[0]))
    write_tsv(output / "chapter-audit.tsv", [row for row in checks.rows if row["Category"] == "Chapter"], ["Category", "CheckID", "Status", "Detail"])
    write_tsv(output / "validation-checks.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])
    validation = {
        "article": 35,
        "status": "passed" if checks.failed == 0 else "failed",
        "checks_total": len(checks.rows),
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "catalog_genes": CATALOG_GENES,
        "input_reads": summary["input_reads"],
        "raw_mapping_rate_pct": summary["raw_mapping_rate_pct"],
        "main_assigned_reads": summary["main_assigned_reads"],
        "main_assigned_pct": summary["main_assigned_pct"],
        "main_detected_genes": summary["main_detected_genes"],
        "uniref90_hit_genes": summary["uniref90_hit_genes"],
        "genes_with_reaction": summary["genes_with_reaction"],
        "assigned_read_pct_with_reaction": summary["assigned_read_pct_with_reaction"],
        "copy_cpm_inflation": summary["copy_cpm_inflation"],
        "count_unit": "reads",
        "mapping_rate_is_assigned_fraction": False,
        "rpkm_is_raw_count_matrix": False,
        "tpm_is_absolute_abundance": False,
        "copy_aggregation_conserves_mass": False,
        "independent_prediction_claimed": False,
        "partial_orf_hit_is_complete_function": False,
        "seed": SEED,
        "figures": list(FIGURE_STEMS),
    }
    (output / "validation-summary.json").write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    (output / "validation.log").write_text(
        f"Article 35 validation: {validation['status']}\nChecks: {checks.passed}/{len(checks.rows)} passed\n",
        encoding="utf-8",
    )
    print(json.dumps(validation))
    return 0 if checks.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
