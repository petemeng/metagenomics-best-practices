#!/usr/bin/env python3
"""Validate Article 36 frozen evidence and draw four publication-ready figures."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, TextIO

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metagenome-article36-matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


CATALOG_GENES = 93_782
EXPECTED_SAMPLE_READS = {"MOCK1": 2_784_234, "MOCK2": 2_777_443}
EXPECTED_SEED_ORTHOLOG_GENES = 84_511
EXPECTED_STATES: dict[str, int] = {
    "No seed ortholog": 9_271,
    "Orthology only": 6_243,
    "Broad/family only": 24_965,
    "Specific identifier": 53_303,
}
EXPECTED_MAIN_GO_GENES = 8_821
EXPECTED_ALL_GO_GENES = 20_471
EXPECTED_UNIQUE_KOS = 5_860
ENV_YAML_SHA256 = "7931295cd72876a5cf61d6001301f38f05db3bb54dcebb7cd55ae029d1a8cf7d"
ENV_LOCK_SHA256 = "573e3f29d9e9fc033621346a4ee02812daaa89afd6918cff0c0caaeb2bd250b1"
FIGURE_STEMS = (
    "36-field-coverage",
    "36-functional-dark-matter",
    "36-cog-category-profile",
    "36-go-evidence-policy",
)
STATE_ORDER = (
    "No seed ortholog",
    "Orthology only",
    "Broad/family only",
    "Specific identifier",
)
STATE_COLORS = {
    "No seed ortholog": "#BDBDBD",
    "Orthology only": "#E69F00",
    "Broad/family only": "#56B4E9",
    "Specific identifier": "#009E73",
}
SAMPLE_COLORS = {"Catalog genes": "#4D4D4D", "MOCK1 reads": "#0072B2", "MOCK2 reads": "#D55E00"}
VISIBLE_LABELS = (
    "Genes or assigned read mass (%)",
    "Catalog genes",
    "MOCK1 reads",
    "MOCK2 reads",
    "Operational annotation state (%)",
    "Fractional share within COG-assigned mass (%)",
    "Non-electronic",
    "All evidence",
    "Genes with GO",
    "Unique GO terms",
    "Gene-GO links",
)


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
    checks.add("Frozen input", "payload-count", len(rows) >= 35, len(rows))
    return rows


def set_pub_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def draw_figures(frozen: Path, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    set_pub_style()
    coverage = read_tsv(frozen / "field-coverage-summary.tsv")
    states = read_tsv(frozen / "annotation-state-summary.tsv")
    cog = read_tsv(frozen / "cog-category-summary.tsv")
    go_audit = read_tsv(frozen / "go-evidence-audit.tsv")

    fields = [row["Field"] for row in coverage][::-1]
    y = np.arange(len(fields))
    fig, ax = plt.subplots(figsize=(8.4, 5.4), constrained_layout=True)
    series = (
        ("Catalog genes", [float(row["GenePercent"]) for row in coverage][::-1]),
        ("MOCK1 reads", [float(row["MOCK1ReadPercent"]) for row in coverage][::-1]),
        ("MOCK2 reads", [float(row["MOCK2ReadPercent"]) for row in coverage][::-1]),
    )
    offsets = (-0.22, 0.0, 0.22)
    for (label, values), offset in zip(series, offsets):
        ax.scatter(values, y + offset, s=38, color=SAMPLE_COLORS[label], label=label, zorder=3)
    ax.set_yticks(y, fields)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Genes or assigned read mass (%)")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.01))
    save_figure(fig, figure_dir / "36-field-coverage")

    labels = ["Catalog genes", "MOCK1 reads", "MOCK2 reads"]
    matrix = np.array(
        [
            [float(next(row for row in states if row["AnnotationState"] == state)["GenePercent"]) for state in STATE_ORDER],
            [float(next(row for row in states if row["AnnotationState"] == state)["MOCK1ReadPercent"]) for state in STATE_ORDER],
            [float(next(row for row in states if row["AnnotationState"] == state)["MOCK2ReadPercent"]) for state in STATE_ORDER],
        ]
    )
    fig, ax = plt.subplots(figsize=(7.4, 5.1), constrained_layout=True)
    bottom = np.zeros(3)
    x = np.arange(3)
    for index, state in enumerate(STATE_ORDER):
        values = matrix[:, index]
        bars = ax.bar(x, values, bottom=bottom, width=0.68, color=STATE_COLORS[state], label=state)
        for bar, value, base in zip(bars, values, bottom):
            if value >= 5:
                ax.text(bar.get_x() + bar.get_width() / 2, base + value / 2, f"{value:.1f}%", ha="center", va="center", fontsize=8)
        bottom += values
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Operational annotation state (%)")
    ax.legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.01))
    save_figure(fig, figure_dir / "36-functional-dark-matter")

    observed_cog = [row for row in cog if float(row["FractionalGeneEquivalent"]) > 0]
    observed_cog.sort(key=lambda row: float(row["FractionalGeneEquivalent"]), reverse=True)
    selected = observed_cog[:15][::-1]
    gene_total = sum(float(row["FractionalGeneEquivalent"]) for row in observed_cog)
    mock1_total = sum(float(row["MOCK1FractionalRawReads"]) for row in observed_cog)
    mock2_total = sum(float(row["MOCK2FractionalRawReads"]) for row in observed_cog)
    y = np.arange(len(selected))
    fig, ax = plt.subplots(figsize=(9.8, 6.8), constrained_layout=True)
    width = 0.24
    cog_labels = [f"{row['COGCategory']}  {row['Description']}" for row in selected]
    cog_series = (
        ("Catalog genes", [100 * float(row["FractionalGeneEquivalent"]) / gene_total for row in selected]),
        ("MOCK1 reads", [100 * float(row["MOCK1FractionalRawReads"]) / mock1_total for row in selected]),
        ("MOCK2 reads", [100 * float(row["MOCK2FractionalRawReads"]) / mock2_total for row in selected]),
    )
    for index, (label, values) in enumerate(cog_series):
        ax.barh(y + (index - 1) * width, values, height=width, color=SAMPLE_COLORS[label], label=label)
    ax.set_yticks(y, cog_labels)
    ax.set_xlabel("Fractional share within COG-assigned mass (%)")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.01))
    save_figure(fig, figure_dir / "36-cog-category-profile")

    selected_metrics = ["Genes with >=1 GO term", "Unique GO terms", "Gene-GO links"]
    metric_labels = ["Genes with GO", "Unique GO terms", "Gene-GO links"]
    fig, axes = plt.subplots(1, 3, figsize=(9.4, 3.8), constrained_layout=True)
    for ax, metric, label in zip(axes, selected_metrics, metric_labels):
        row = next(item for item in go_audit if item["Metric"] == metric)
        values = [int(row["NonElectronic"]), int(row["AllEvidence"])]
        bars = ax.bar([0, 1], values, color=["#0072B2", "#CC79A7"], width=0.65)
        ax.set_xticks([0, 1], ["Non-electronic", "All evidence"], rotation=20, ha="right")
        ax.set_title(label)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 4))
        ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, value, f"{value:,}", ha="center", va="bottom", fontsize=8)
    save_figure(fig, figure_dir / "36-go-evidence-policy")


def audit_images(figure_dir: Path, checks: Checks) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for stem in FIGURE_STEMS:
        for suffix in ("pdf", "png", "tiff"):
            path = figure_dir / f"{stem}.{suffix}"
            exists = path.is_file() and path.stat().st_size > 0
            checks.add("Figures", f"exists-{stem}-{suffix}", exists, path)
            detail = "missing"
            status = exists
            if exists and suffix in {"png", "tiff"}:
                with Image.open(path) as image:
                    width, height = image.size
                    dpi = image.info.get("dpi", (0, 0))
                    min_dpi = min(float(dpi[0]), float(dpi[1])) if dpi else 0
                    detail = f"{width}x{height}; dpi={dpi}"
                    status = width >= 1500 and height >= 900 and min_dpi >= 340
                    checks.add("Figures", f"raster-quality-{stem}-{suffix}", status, detail)
            rows.append({"Figure": path.name, "Status": "PASS" if status else "FAIL", "Detail": detail})
    visible_text = "\n".join(VISIBLE_LABELS)
    checks.add("Figures", "visible-labels-english", re.search(r"[\u3400-\u9fff]", visible_text) is None, visible_text)
    return rows


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = (args.frozen_dir or root / "data/small/36-eggnog-functional-annotation-frozen").resolve()
    output = (args.output_dir or root / "results/36-eggnog-functional-annotation").resolve()
    figure_dir = (args.figure_dir or root / "figures").resolve()
    chapter = (args.chapter or root / "chapters/36-eggnog-functional-annotation.qmd").resolve()
    output.mkdir(parents=True, exist_ok=True)
    checks = Checks()

    checksum_rows = verify_checksum_manifest(frozen, checks)
    contract = json.loads((frozen / "frozen-contract.json").read_text(encoding="utf-8"))
    summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    checks.add("Identity", "catalog-genes", summary["catalog_genes"] == CATALOG_GENES, summary["catalog_genes"])
    checks.add("Identity", "seed-ortholog-genes", summary["seed_ortholog_genes"] == EXPECTED_SEED_ORTHOLOG_GENES, summary["seed_ortholog_genes"])
    checks.add("Identity", "environment-yaml", sha256(frozen / "env/eggnog-annotation.yml") == ENV_YAML_SHA256, sha256(frozen / "env/eggnog-annotation.yml"))
    checks.add("Identity", "environment-lock", sha256(frozen / "env/eggnog-annotation-linux-64.lock") == ENV_LOCK_SHA256, sha256(frozen / "env/eggnog-annotation-linux-64.lock"))
    checks.add("Identity", "mapper-version", contract["eggnog_mapper"] == "2.1.15", contract["eggnog_mapper"])
    checks.add("Identity", "database-release", contract["eggnog_database"] == "5.0.2", contract["eggnog_database"])
    checks.add("Identity", "diamond-version", contract["diamond"] == "2.0.15", contract["diamond"])

    run_contract = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    search_contract = run_contract["search"]
    annotation_contract = run_contract["annotation"]
    checks.add(
        "Method contract",
        "search-policy",
        search_contract == {
            "mode": "diamond",
            "diamond": "2.0.15",
            "sensitivity": "sensitive",
            "iterate": "yes",
            "evalue": 0.001,
            "seed_ortholog_evalue": 0.001,
            "identity_filter": None,
            "query_coverage_filter": None,
            "subject_coverage_filter": None,
            "outfmt_short": True,
        },
        search_contract,
    )
    checks.add(
        "Method contract",
        "annotation-policy",
        annotation_contract == {
            "tax_scope": "auto",
            "tax_scope_mode": "inner_narrowest",
            "target_orthologs": "all",
            "go_evidence_primary": "non-electronic",
            "go_evidence_sensitivity": "all",
            "pfam_realign": "none",
        },
        annotation_contract,
    )

    database_rows = read_tsv(frozen / "database-manifest.tsv")
    checks.add("Database", "manifest-assets", len(database_rows) == 4, len(database_rows))
    checks.add("Database", "releases-locked", {row["Release"] for row in database_rows} == {"5.0.2"}, {row["Release"] for row in database_rows})
    checks.add("Database", "no-latest-url", all("latest" not in row["SourceURL"] for row in database_rows), "immutable 5.0.2 URLs")
    checks.add("Database", "installed-hashes", all(re.fullmatch(r"[0-9a-f]{64}", row["InstalledSHA256"]) for row in database_rows), "four SHA-256 values")
    write_tsv(output / "database-audit.tsv", database_rows, list(database_rows[0]))

    gene_rows = read_tsv(frozen / "gene-functional-annotation.tsv.gz")
    gene_ids = [row["GeneID"] for row in gene_rows]
    state_counts = Counter(row["AnnotationState"] for row in gene_rows)
    checks.add("Annotation", "gene-table-rows", len(gene_rows) == CATALOG_GENES, len(gene_rows))
    checks.add("Annotation", "unique-gene-ids", len(set(gene_ids)) == CATALOG_GENES, len(set(gene_ids)))
    checks.add("Annotation", "state-vocabulary", set(state_counts) == set(STATE_ORDER), dict(state_counts))
    checks.add("Annotation", "state-partition", sum(state_counts.values()) == CATALOG_GENES, sum(state_counts.values()))
    checks.add("Annotation", "exact-state-counts", dict(state_counts) == EXPECTED_STATES, dict(state_counts))
    checks.add("Annotation", "exact-main-go", summary["main_go_genes"] == EXPECTED_MAIN_GO_GENES, summary["main_go_genes"])
    checks.add("Annotation", "exact-all-go", summary["all_evidence_go_genes"] == EXPECTED_ALL_GO_GENES, summary["all_evidence_go_genes"])
    checks.add("Annotation", "exact-unique-kos", summary["unique_kos"] == EXPECTED_UNIQUE_KOS, summary["unique_kos"])
    checks.add("Annotation", "go-policy-primary", contract["go_evidence_primary"] == "non-electronic", contract["go_evidence_primary"])
    checks.add("Annotation", "go-all-sensitivity-only", contract["go_evidence_sensitivity"] == "all" and summary["go_all_is_primary"] is False, summary["go_all_is_primary"])
    checks.add("Annotation", "go-term-subset", summary["go_subset_violations"] == 0, summary["go_subset_violations"])
    checks.add("Annotation", "prediction-boundary", summary["functional_annotations_are_predictions"] is True, summary["functional_annotations_are_predictions"])
    checks.add("Annotation", "absence-boundary", summary["absence_is_not_gene_absence"] is True, summary["absence_is_not_gene_absence"])
    checks.add("Annotation", "fractional-cog", summary["multi_label_cog_uses_fractional_allocation"] is True, summary["multi_label_cog_uses_fractional_allocation"])

    state_rows = read_tsv(frozen / "annotation-state-summary.tsv")
    checks.add("Ledger", "state-summary-rows", len(state_rows) == 4, len(state_rows))
    checks.add("Ledger", "gene-percent-closes", abs(sum(float(row["GenePercent"]) for row in state_rows) - 100) < 1e-8, sum(float(row["GenePercent"]) for row in state_rows))
    for sample in ("MOCK1", "MOCK2"):
        read_sum = sum(int(row[f"{sample}RawReads"]) for row in state_rows)
        percent_sum = sum(float(row[f"{sample}ReadPercent"]) for row in state_rows)
        checks.add("Ledger", f"{sample.lower()}-read-total", read_sum == EXPECTED_SAMPLE_READS[sample], read_sum)
        checks.add("Ledger", f"{sample.lower()}-percent-closes", abs(percent_sum - 100) < 1e-8, percent_sum)

    overlap_rows = read_tsv(frozen / "annotation-evidence-overlap.tsv")
    allowed_combinations = {
        f"COG={cog} | KO={ko} | GO={go}"
        for cog in ("yes", "no")
        for ko in ("yes", "no")
        for go in ("yes", "no")
    }
    observed_combinations = {row["EvidenceCombination"] for row in overlap_rows}
    checks.add(
        "Ledger",
        "evidence-combination-vocabulary",
        bool(observed_combinations) and observed_combinations.issubset(allowed_combinations),
        sorted(observed_combinations),
    )
    checks.add("Ledger", "evidence-genes-close", sum(int(row["Genes"]) for row in overlap_rows) == CATALOG_GENES, sum(int(row["Genes"]) for row in overlap_rows))
    go_audit = read_tsv(frozen / "go-evidence-audit.tsv")
    go_genes = next(row for row in go_audit if row["Metric"] == "Genes with >=1 GO term")
    checks.add("GO evidence", "all-not-smaller", int(go_genes["AllEvidence"]) >= int(go_genes["NonElectronic"]), go_genes)
    checks.add("GO evidence", "main-summary-match", int(go_genes["NonElectronic"]) == summary["main_go_genes"], go_genes["NonElectronic"])
    checks.add("GO evidence", "all-summary-match", int(go_genes["AllEvidence"]) == summary["all_evidence_go_genes"], go_genes["AllEvidence"])

    fractional_rows = read_tsv(frozen / "fractional-allocation-audit.tsv")
    checks.add("Fractional allocation", "three-layers", {row["Layer"] for row in fractional_rows} == {"COG", "KO", "GO"}, [row["Layer"] for row in fractional_rows])
    for row in fractional_rows:
        layer = row["Layer"].lower()
        gene_tolerance = max(1e-6, 1e-9 * int(row["GenesWithAnyTerm"]))
        checks.add("Fractional allocation", f"{layer}-gene-conservation", abs(float(row["GeneDifference"])) <= gene_tolerance, row["GeneDifference"])
        for sample in ("MOCK1", "MOCK2"):
            read_tolerance = max(1e-6, 1e-9 * int(row[f"{sample}ReadsFromGenesWithAnyTerm"]))
            checks.add("Fractional allocation", f"{layer}-{sample.lower()}-read-conservation", abs(float(row[f"{sample}Difference"])) <= read_tolerance, row[f"{sample}Difference"])

    annotation_audit = [
        {"Metric": "Catalog genes", "Observed": len(gene_rows), "Expected": CATALOG_GENES, "Status": "PASS" if len(gene_rows) == CATALOG_GENES else "FAIL"},
        {"Metric": "Seed ortholog genes", "Observed": summary["seed_ortholog_genes"], "Expected": EXPECTED_SEED_ORTHOLOG_GENES, "Status": "PASS" if summary["seed_ortholog_genes"] == EXPECTED_SEED_ORTHOLOG_GENES else "FAIL"},
        {"Metric": "Main GO genes", "Observed": summary["main_go_genes"], "Expected": EXPECTED_MAIN_GO_GENES, "Status": "PASS" if summary["main_go_genes"] == EXPECTED_MAIN_GO_GENES else "FAIL"},
        {"Metric": "All-evidence GO genes", "Observed": summary["all_evidence_go_genes"], "Expected": EXPECTED_ALL_GO_GENES, "Status": "PASS" if summary["all_evidence_go_genes"] == EXPECTED_ALL_GO_GENES else "FAIL"},
    ]
    write_tsv(output / "annotation-audit.tsv", annotation_audit, ["Metric", "Observed", "Expected", "Status"])

    text = chapter.read_text(encoding="utf-8")
    chapter_checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "nine-target": "## 这一步对应论文里的哪张图" in text,
        "nine-theory": "## 理论" in text,
        "nine-setup": "## 准备工作" in text,
        "nine-code": "## 可复制代码" in text,
        "nine-audit": "## 审计与升级" in text,
        "nine-publication": "## 出版级美化" in text,
        "nine-pitfalls": "## 常见坑" in text,
        "nine-methods": "## 这段 Methods 怎么写" in text,
        "nine-own-data": "## 换成你自己的数据怎么做" in text,
        "resource-ram": "RAM" in text and "磁盘" in text and "CPU" in text and "耗时" in text,
        "versions": all(token in text for token in ("eggNOG-mapper 2.1.15", "eggNOG 5.0.2", "DIAMOND 2.0.15")),
        "go-policy": "non-electronic" in text and "all-evidence" in text,
        "seed": "20260736" in text,
        "methods-template": "Methods template" in text,
        "four-images": all(f"../figures/{stem}.png" in text for stem in FIGURE_STEMS),
        "references": "## 参考" in text,
        "no-result-placeholders": re.search(r"__[A-Z0-9_]+__", text) is None,
        "no-meta-prose": not any(phrase in text for phrase in ("本篇可独立", "本文可独立", "全系列约定", "接口只学一次", "作者代码通常长这样", "（即本文）", "无头服务器")),
    }
    chapter_rows = []
    for key, status in chapter_checks.items():
        checks.add("Chapter", key, status, key)
        chapter_rows.append({"Check": key, "Status": "PASS" if status else "FAIL"})
    write_tsv(output / "chapter-audit.tsv", chapter_rows, ["Check", "Status"])

    draw_figures(frozen, figure_dir)
    image_rows = audit_images(figure_dir, checks)
    write_tsv(output / "image-audit.tsv", image_rows, ["Figure", "Status", "Detail"])
    write_tsv(output / "checksum-audit.tsv", checksum_rows, ["File", "ExpectedSHA256", "ObservedSHA256", "Status"])
    write_tsv(output / "validation-checks.tsv", checks.rows, ["Category", "CheckID", "Status", "Detail"])

    validation_summary = {
        "article": 36,
        "status": "passed" if checks.failed == 0 else "failed",
        "checks_passed": checks.passed,
        "checks_failed": checks.failed,
        "catalog_genes": CATALOG_GENES,
        "seed_ortholog_genes": summary["seed_ortholog_genes"],
        "annotation_states": dict(state_counts),
        "main_go_genes": summary["main_go_genes"],
        "all_evidence_go_genes": summary["all_evidence_go_genes"],
        "unique_kos": summary["unique_kos"],
        "eggnog_mapper": contract["eggnog_mapper"],
        "eggnog_database": contract["eggnog_database"],
        "diamond": contract["diamond"],
        "functional_annotations_are_predictions": True,
        "absence_is_not_gene_absence": True,
        "go_all_is_primary": False,
        "multi_label_cog_uses_fractional_allocation": True,
        "multi_label_terms_use_fractional_allocation": True,
    }
    (output / "validation-summary.json").write_text(json.dumps(validation_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    log_lines = [f"Article 36 validation: {validation_summary['status'].upper()}", f"PASS={checks.passed}", f"FAIL={checks.failed}"]
    log_lines.extend(f"{row['Status']}\t{row['Category']}\t{row['CheckID']}\t{row['Detail']}" for row in checks.rows)
    (output / "validation.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(json.dumps(validation_summary, indent=2, sort_keys=True))
    return 0 if checks.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
