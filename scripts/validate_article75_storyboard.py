#!/usr/bin/env python3
"""Validate Article 75's frozen evidence, figures, chapter, and rendered HTML."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
import yaml
from PIL import Image


FIGURE_STEMS = (
    "75-wirbel-figure-arc",
    "75-five-figure-storyboard",
    "75-panel-budget",
    "75-main-supplement",
    "75-claim-evidence",
    "75-traceability-audit",
    "75-sensitivity-matrix",
    "75-style-contract",
    "75-reviewer-attack-map",
)
ANCHOR = "75-wirbel-figure1-original.jpg"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    parser.add_argument("--rendered-html", type=Path, required=True)
    parser.add_argument("--qa-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map({"true": True, "false": False})


class Audit:
    def __init__(self) -> None:
        self.rows: list[dict[str, object]] = []

    def add(self, category: str, check: str, passed: object, observed: object = "") -> None:
        self.rows.append(
            {
                "Category": category,
                "Check": check,
                "Passed": bool(passed),
                "Observed": str(observed),
            }
        )

    def write(self, output: Path) -> None:
        output.mkdir(parents=True, exist_ok=True)
        frame = pd.DataFrame(self.rows)
        frame.to_csv(output / "checks.tsv", sep="\t", index=False, lineterminator="\n")
        failed = frame.loc[~frame["Passed"]]
        report = {
            "article": 75,
            "status": "passed" if failed.empty else "failed",
            "checks": len(frame),
            "passed": int(frame["Passed"].sum()),
            "failed": len(failed),
            "failures": failed[["Category", "Check", "Observed"]].to_dict("records"),
        }
        (output / "qa_report.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(report, indent=2))
        if not failed.empty:
            raise SystemExit(1)


def audit_bundle(frozen: Path, audit: Audit) -> None:
    manifest = json.loads((frozen / "bundle-manifest.json").read_text())
    audit.add("Bundle", "article", manifest.get("article") == 75, manifest.get("article"))
    audit.add("Bundle", "payload-files", manifest.get("payload_files") == 19, manifest.get("payload_files"))
    audit.add("Bundle", "script-files", manifest.get("script_files") == 5, manifest.get("script_files"))
    audit.add("Bundle", "environment-files", manifest.get("environment_files") == 2, manifest.get("environment_files"))
    audit.add("Bundle", "cross-dataset-boundary", "cannot be combined" in manifest.get("source_boundary", ""), manifest.get("source_boundary"))

    checksum_path = frozen / "file-checksums.sha256"
    lines = [line for line in checksum_path.read_text().splitlines() if line.strip()]
    audit.add("Checksum", "record-count", len(lines) == 27, len(lines))
    recorded = set()
    for line in lines:
        digest, relative = line.split("  ", 1)
        path = frozen / relative
        recorded.add(relative)
        audit.add("Checksum file", relative, path.is_file() and sha256(path) == digest, path)
    payload = {
        path.relative_to(frozen).as_posix()
        for path in frozen.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    }
    audit.add("Checksum", "coverage", payload == recorded, f"payload={len(payload)} recorded={len(recorded)}")


def audit_sources(root: Path, frozen: Path, audit: Audit) -> None:
    sources = pd.read_csv(frozen / "source-artifact-manifest.tsv", sep="\t")
    audit.add("Sources", "artifact-count", len(sources) == 15, len(sources))
    audit.add("Sources", "unique-ids", sources["ArtifactID"].is_unique, sources["ArtifactID"].tolist())
    for row in sources.itertuples(index=False):
        path = root / row.RelativePath
        audit.add("Source artifact", row.ArtifactID, path.is_file() and path.stat().st_size == row.Bytes and sha256(path) == row.SHA256, row.RelativePath)

    paper = json.loads((frozen / "paper-source-manifest.json").read_text())
    audit.add("Paper", "article", paper.get("article") == 75, paper.get("article"))
    audit.add("Paper", "doi", paper["paper"]["doi"] == "10.1038/s41591-019-0406-6", paper["paper"]["doi"])
    audit.add("Paper", "pmcid", paper.get("pmcid") == "PMC7984229", paper.get("pmcid"))
    audit.add("Paper", "figure-members", len(paper["selected_figure_members"]) == 5, len(paper["selected_figure_members"]))
    anchor = frozen / "wirbel-figure1-original.jpg"
    audit.add("Paper", "anchor-bytes", anchor.stat().st_size == paper["anchor"]["bytes"], anchor.stat().st_size)
    audit.add("Paper", "anchor-sha", sha256(anchor) == paper["anchor"]["sha256"], sha256(anchor))
    with Image.open(anchor) as image:
        audit.add("Paper", "anchor-dimensions", image.size == (800, 941), image.size)
    audit.add("Paper", "rights-boundary", "excluded" in paper["rights_boundary"], paper["rights_boundary"])


def audit_tables(frozen: Path, audit: Audit) -> None:
    paper = pd.read_csv(frozen / "wirbel-main-figure-ledger.tsv", sep="\t")
    audit.add("Paper ledger", "figures", paper["Figure"].tolist() == [1, 2, 3, 4, 5], paper["Figure"].tolist())
    audit.add("Paper ledger", "panels", paper["Panels"].tolist() == [4, 4, 5, 6, 3], paper["Panels"].tolist())
    audit.add("Paper ledger", "panel-total", int(paper["Panels"].sum()) == 22, paper["Panels"].sum())
    for row in paper.itertuples(index=False):
        audit.add("Paper ledger hash", row.SourceLabel, bool(re.fullmatch(r"[0-9a-f]{64}", row.CaptionSHA256)) and bool(re.fullmatch(r"[0-9a-f]{64}", row.ImageSHA256)), row.SourceImage)

    storyboard = pd.read_csv(frozen / "main-figure-storyboard.tsv", sep="\t")
    audit.add("Storyboard", "figures", storyboard["Figure"].tolist() == [1, 2, 3, 4, 5], storyboard["Figure"].tolist())
    audit.add("Storyboard", "panel-budgets", storyboard["Panels"].tolist() == [4, 5, 6, 5, 5], storyboard["Panels"].tolist())
    audit.add("Storyboard", "validation-gates", storyboard["ValidationGate"].notna().all(), storyboard["ValidationGate"].tolist())
    audit.add("Storyboard", "stop-rules", storyboard["StopRule"].str.len().gt(20).all(), storyboard["StopRule"].tolist())

    panels = pd.read_csv(frozen / "panel-register.tsv", sep="\t")
    audit.add("Panels", "rows", len(panels) == 25, len(panels))
    audit.add("Panels", "ids-unique", panels["ResultID"].is_unique, panels["ResultID"].tolist())
    audit.add("Panels", "figures", set(panels["Figure"]) == {1, 2, 3, 4, 5}, sorted(set(panels["Figure"])))
    counts = panels.groupby("Figure").size().tolist()
    audit.add("Panels", "counts-match-storyboard", counts == storyboard["Panels"].tolist(), counts)
    audit.add("Panels", "units-complete", panels["AnalysisUnit"].str.len().gt(1).all(), panels["AnalysisUnit"].tolist())
    audit.add("Panels", "gates-complete", panels["Gate"].str.len().gt(5).all(), panels["Gate"].tolist())

    supplement = pd.read_csv(frozen / "main-supplement-map.tsv", sep="\t")
    audit.add("Supplement", "blocks", len(supplement) == 14, len(supplement))
    audit.add("Supplement", "full-detail", supplement["SupplementDetail"].eq(2).all(), supplement["SupplementDetail"].tolist())
    audit.add("Supplement", "main-levels", set(supplement["MainSpace"]) == {0, 1, 2}, sorted(set(supplement["MainSpace"])))
    audit.add("Supplement", "nulls-retained", supplement["EvidenceBlock"].str.contains("failed", case=False).any(), supplement["EvidenceBlock"].tolist())

    claims = pd.read_csv(frozen / "claim-evidence-matrix.tsv", sep="\t")
    audit.add("Claims", "cells", len(claims) == 30, len(claims))
    audit.add("Claims", "six-layers-each", claims.groupby("ClaimID").size().eq(6).all(), claims.groupby("ClaimID").size().to_dict())
    clinical = claims[claims["ClaimID"].eq("C5") & claims["EvidenceLayer"].eq("Perturbation or intervention")]
    audit.add("Claims", "causal-needs-intervention", len(clinical) == 1 and clinical.iloc[0]["Requirement"] == "Required", clinical.to_dict("records"))
    mechanism = claims[claims["ClaimID"].eq("C4") & claims["EvidenceLayer"].eq("Model or orthogonal assay")]
    audit.add("Claims", "mechanism-needs-assay", len(mechanism) == 1 and mechanism.iloc[0]["Requirement"] == "Required", mechanism.to_dict("records"))

    ladder = pd.read_csv(frozen / "evidence-language-ladder.tsv", sep="\t")
    audit.add("Evidence ladder", "rungs", ladder["Rung"].tolist() == [1, 2, 3, 4, 5, 6], ladder["Rung"].tolist())
    audit.add("Evidence ladder", "forbidden-leaps", ladder["ForbiddenLeap"].str.len().gt(10).all(), ladder["ForbiddenLeap"].tolist())

    sensitivity = pd.read_csv(frozen / "sensitivity-matrix.tsv", sep="\t")
    audit.add("Sensitivity", "cells", len(sensitivity) == 45, len(sensitivity))
    audit.add("Sensitivity", "nine-axes", sensitivity["SensitivityAxis"].nunique() == 9, sensitivity["SensitivityAxis"].nunique())
    audit.add("Sensitivity", "five-figures-each", sensitivity.groupby("SensitivityAxis").size().eq(5).all(), sensitivity.groupby("SensitivityAxis").size().to_dict())
    mag = sensitivity[sensitivity["SensitivityAxis"].eq("MAG quality threshold")].set_index("Figure")
    audit.add("Sensitivity", "mag-main-branch", mag.loc[3, "Requirement"] == "Required" and mag.loc[4, "Requirement"] == "Required", mag["Requirement"].to_dict())

    attacks = pd.read_csv(frozen / "reviewer-attack-map.tsv", sep="\t")
    audit.add("Reviewer map", "rows", len(attacks) == 10, len(attacks))
    audit.add("Reviewer map", "risk-range", attacks["Likelihood"].between(3, 5).all() and attacks["Impact"].between(4, 5).all(), attacks[["Likelihood", "Impact"]].to_dict("records"))
    for required in ("Outcome leakage", "Cohort confounding", "Missing external validation", "Incomplete MAG QC", "Unsupported host link", "Version drift"):
        audit.add("Reviewer attack", required, required in set(attacks["Attack"]), attacks["Attack"].tolist())

    style = pd.read_csv(frozen / "figure-style-contract.tsv", sep="\t")
    audit.add("Style", "rows", len(style) == 8, len(style))
    audit.add("Style", "hex-colors", style["Hex"].str.fullmatch(r"#[0-9A-Fa-f]{6}").all(), style["Hex"].tolist())
    audit.add("Style", "control-case-distinct", style.set_index("Semantic").loc["Control", "Hex"] != style.set_index("Semantic").loc["Case", "Hex"], style[["Semantic", "Hex"]].to_dict("records"))

    trace = pd.read_csv(frozen / "result-traceability-ledger.tsv", sep="\t")
    audit.add("Traceability", "rows", len(trace) == 25, len(trace))
    for column in ("UnitRecorded", "StatisticRecorded", "ValidationRecorded", "ChecksumRequired"):
        values = parse_bool(trace[column])
        audit.add("Traceability", column, values.notna().all() and values.all(), values.tolist())
    audit.add("Traceability", "code-targets", trace["CodeTarget"].str.fullmatch(r"analysis/f[1-5][a-f]\.R").all(), trace["CodeTarget"].tolist())
    audit.add("Traceability", "output-targets", trace["OutputTarget"].str.fullmatch(r"results/f[1-5][a-f]\.tsv").all(), trace["OutputTarget"].tolist())

    metrics = pd.read_csv(frozen / "series-evidence-metrics.tsv", sep="\t")
    audit.add("Evidence metrics", "rows", len(metrics) == 13, len(metrics))
    indexed = metrics.set_index("EvidenceID")
    expected = {"E01": 8, "E02": 771, "E05": 23, "E06": 11, "E08": 46, "E09": 12, "E10": 8128, "E11": 48, "E12": 6, "E13": 0}
    for key, value in expected.items():
        audit.add("Evidence metric", key, float(indexed.loc[key, "Value"]) == value, indexed.loc[key].to_dict())
    audit.add("Evidence metric", "macro-auroc", abs(float(indexed.loc["E03", "Value"]) - 0.7786393435) < 1e-9, indexed.loc["E03", "Value"])
    audit.add("Evidence metric", "heterogeneity", abs(float(indexed.loc["E04", "Value"]) - 52.09417378) < 1e-8, indexed.loc["E04", "Value"])
    audit.add("Evidence metric", "boundaries-complete", metrics["EvidenceBoundary"].str.len().gt(20).all(), metrics["EvidenceBoundary"].tolist())

    versions = pd.read_csv(frozen / "version-ledger-example.tsv", sep="\t")
    audit.add("Versions", "rows", len(versions) == 4, len(versions))
    audit.add("Versions", "named-releases", versions["Release"].str.len().gt(3).all(), versions["Release"].tolist())
    audit.add("Versions", "checksums", versions["Checksum"].str.len().gt(10).all(), versions["Checksum"].tolist())

    analysis = json.loads((frozen / "analysis-metrics.json").read_text())
    expected_analysis = {
        "article": 75,
        "analysis_seed": 75001,
        "plot_seed": 20260775,
        "anchor_paper_figures": 5,
        "anchor_paper_panels": 22,
        "storyboard_figures": 5,
        "storyboard_panels": 25,
        "supplement_blocks": 14,
        "claim_matrix_cells": 30,
        "sensitivity_cells": 45,
        "traceability_results": 25,
    }
    for key, value in expected_analysis.items():
        audit.add("Analysis metrics", key, analysis.get(key) == value, analysis.get(key))


def audit_figures(root: Path, frozen: Path, audit: Audit) -> None:
    manifest = json.loads((frozen / "figure-manifest.json").read_text())
    audit.add("Figure", "article", manifest.get("article") == 75, manifest.get("article"))
    audit.add("Figure", "plot-seed", manifest.get("plot_seed") == 20260775, manifest.get("plot_seed"))
    records = manifest.get("files", {})
    audit.add("Figure", "file-count", len(records) == 19, len(records))
    for stem in FIGURE_STEMS:
        for suffix in (".png", ".svg"):
            name = stem + suffix
            path = root / "figures/article75" / name
            record = records.get(name, {})
            audit.add("Figure file", name, path.is_file() and path.stat().st_size == record.get("bytes") and sha256(path) == record.get("sha256"), path)
            if suffix == ".png" and path.is_file():
                with Image.open(path) as image:
                    audit.add("Figure resolution", name, image.width >= 2800 and image.height >= 1600, image.size)
            if suffix == ".svg" and path.is_file():
                vector = path.read_text(errors="ignore")
                audit.add("Figure language", name, not bool(re.search(r"[\u4e00-\u9fff]", vector)), "English-only vector text")
    anchor = root / "figures/article75" / ANCHOR
    record = records.get(ANCHOR, {})
    audit.add("Figure file", ANCHOR, anchor.is_file() and anchor.stat().st_size == record.get("bytes") and sha256(anchor) == record.get("sha256"), anchor)
    with Image.open(anchor) as image:
        audit.add("Figure anchor", "dimensions", image.size == (800, 941), image.size)


def frontmatter(text: str) -> dict[str, object]:
    match = re.match(r"^---\n(.*?)\n---\n", text, flags=re.S)
    if not match:
        raise ValueError("Missing YAML frontmatter")
    return yaml.safe_load(match.group(1))


def audit_chapter(root: Path, chapter: Path, rendered: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    meta = frontmatter(text)
    audit.add("Chapter", "title", "第 75 篇" in str(meta.get("title")), meta.get("title"))
    audit.add("Chapter", "draft-false", meta.get("draft") is False, meta.get("draft"))
    audit.add("Chapter", "eval-true", meta.get("execute", {}).get("eval") is True, meta.get("execute"))
    audit.add("Chapter", "freeze-auto", meta.get("execute", {}).get("freeze") == "auto", meta.get("execute"))
    audit.add("Chapter", "bibliography", meta.get("bibliography") == "../references.bib", meta.get("bibliography"))
    audit.add("Chapter", "wechat-images", meta.get("wechat", {}).get("expected_images") == 10, meta.get("wechat", {}).get("expected_images"))
    headings = [
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
    ]
    for heading in headings:
        audit.add("Chapter section", heading, heading in text, heading)
    image_names = [ANCHOR] + [stem + ".png" for stem in FIGURE_STEMS]
    for name in image_names:
        audit.add("Chapter image", name, name in text, name)
    markdown_images = re.findall(r"!\[[^\]]*\]\([^)]*(75-[^/)]+\.(?:png|jpg))[^)]*\)", text)
    audit.add("Chapter", "image-count", len(markdown_images) == 10 and set(markdown_images) == set(image_names), markdown_images)
    audit.add("Chapter", "inline-helpers", all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")), "pal_pub/theme_pub/save_pub")
    audit.add("Chapter", "no-shared-source", 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text, "source() absent")
    audit.add("Chapter", "seed", "set.seed(20260775)" in text, "set.seed(20260775)")
    audit.add("Chapter", "anchor-rights", "权利归原权利人" in text and "不属于本仓库" in text, "rights boundary")
    audit.add("Chapter", "cross-dataset-boundary", "不能把这些数值拼成同一个生物学故事" in text, "tutorial examples remain separate")
    audit.add("Chapter", "methods-versions", all(token in text for token in ("R v4.4.1", "ggplot2 v3.5.1", "20260775")), "R/ggplot2/seed")
    for phrase in ("本篇可独立跑通", "这体现全系列", "接口只学一次", "已建好", "作者代码通常长这样", "（即本文）"):
        audit.add("Chapter prose", f"forbid:{phrase}", phrase not in text, phrase)
    audit.add("References", "wirbel-key", "wirbel2019crc" in (root / "references.bib").read_text(), "wirbel2019crc")
    audit.add("References", "citation-used", "[@wirbel2019crc]" in text, "[@wirbel2019crc]")

    audit.add("Render", "html-exists", rendered.is_file(), rendered)
    if rendered.is_file():
        html = rendered.read_text(encoding="utf-8", errors="ignore")
        audit.add("Render", "html-size", rendered.stat().st_size > 150_000, rendered.stat().st_size)
        audit.add("Render", "title-present", "一篇宏基因组论文的主图和补图如何组织" in html, "title")
        for name in image_names:
            audit.add("Render image", name, name in html, name)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    chapter = args.chapter.resolve()
    rendered = args.rendered_html.resolve()
    qa = args.qa_dir.resolve()
    audit = Audit()
    audit_bundle(frozen, audit)
    audit_sources(root, frozen, audit)
    audit_tables(frozen, audit)
    audit_figures(root, frozen, audit)
    audit_chapter(root, chapter, rendered, audit)
    audit.write(qa)


if __name__ == "__main__":
    main()
