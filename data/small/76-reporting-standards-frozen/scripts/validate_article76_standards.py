#!/usr/bin/env python3
"""Validate Article 76's frozen bundle, figures, chapter and rendered HTML."""

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
    "76-standard-selection",
    "76-reporting-layer-stack",
    "76-checklist-sections",
    "76-standards-crosswalk",
    "76-mimag-compliance",
    "76-miuvig-readiness",
    "76-owner-timeline",
    "76-na-ledger",
    "76-submission-readiness",
)
ANCHOR = "76-streams-figure1-original.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    parser.add_argument("--rendered-html", type=Path, required=True)
    parser.add_argument("--qa-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map({"true": True, "false": False})


class Audit:
    def __init__(self) -> None:
        self.rows: list[dict[str, object]] = []

    def add(self, category: str, check: str, passed: object, observed: object = "") -> None:
        self.rows.append(
            {"Category": category, "Check": check, "Passed": bool(passed), "Observed": str(observed)}
        )

    def write(self, output: Path) -> None:
        output.mkdir(parents=True, exist_ok=True)
        frame = pd.DataFrame(self.rows)
        frame.to_csv(output / "checks.tsv", sep="\t", index=False, lineterminator="\n")
        failed = frame.loc[~frame["Passed"]]
        report = {
            "article": 76,
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
    manifest = json.loads((frozen / "bundle-manifest.json").read_text(encoding="utf-8"))
    audit.add("Bundle", "article", manifest.get("article") == 76, manifest.get("article"))
    audit.add("Bundle", "payload-files", manifest.get("payload_files") == 20, manifest.get("payload_files"))
    audit.add("Bundle", "script-files", manifest.get("script_files") == 5, manifest.get("script_files"))
    audit.add("Bundle", "environment-files", manifest.get("environment_files") == 2, manifest.get("environment_files"))
    audit.add("Bundle", "independent-example-boundary", "cannot be combined" in manifest.get("source_boundary", ""), manifest.get("source_boundary"))

    lines = [line for line in (frozen / "file-checksums.sha256").read_text().splitlines() if line.strip()]
    audit.add("Checksum", "record-count", len(lines) == 28, len(lines))
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
    sources = json.loads((frozen / "source-manifest.json").read_text(encoding="utf-8"))
    audit.add("Sources", "article", sources.get("article") == 76, sources.get("article"))
    audit.add("Sources", "record-count", len(sources.get("records", [])) == 5, len(sources.get("records", [])))
    audit.add("Sources", "streams-sheet", sources.get("streams_sheet") == "STREAMS_final", sources.get("streams_sheet"))
    audit.add("Sources", "draft-sheets-excluded", "STREAMS v2.0" in sources.get("excluded_streams_sheets", []), sources.get("excluded_streams_sheets"))
    expected = {
        "STREAMS_Guidelines_Zenodo.xlsx": "c3905dbdd28f256968a06157682277a50712ad045b9624a0f11cce47414ff788",
        "STORMS_Excel_1.03.xlsx": "4d763f0d62ee27aa43217a0ebe47e75b9beb70c8534f867c5c2fc3b6f1714b60",
        "PMC6436528-fulltext.xml": "7fdec05ef70b15e852d21225fbcae8b0ecc24e4accdaf5075d1e976ab9295ba4",
        "PMC6871006-fulltext.xml": "289c221ece18f17a0a8930a354cfc3e01b712751abfdd52ff7b054c6f31c36fa",
        "streams-figure1-original.png": "993fb497cc09e2e3446b22233ca70ab9bacff9042babce39955251288a34279d",
    }
    records = {record["File"]: record for record in sources["records"]}
    audit.add("Sources", "expected-files", set(records) == set(expected), sorted(records))
    for filename, digest in expected.items():
        record = records.get(filename, {})
        audit.add("Source lock", filename, record.get("SHA256") == digest, record.get("SHA256"))
    audit.add("Sources", "rights-boundary", "excluded" in sources["anchor"]["rights_boundary"], sources["anchor"]["rights_boundary"])
    anchor = frozen / "streams-figure1-original.png"
    audit.add("Sources", "anchor-sha", anchor.is_file() and sha256(anchor) == expected[anchor.name], sha256(anchor) if anchor.is_file() else "missing")
    with Image.open(anchor) as image:
        audit.add("Sources", "anchor-dimensions", image.size == (2174, 1098), image.size)

    artifacts = pd.read_csv(frozen / "source-artifact-manifest.tsv", sep="\t")
    audit.add("Artifacts", "count", len(artifacts) == 8, len(artifacts))
    audit.add("Artifacts", "unique-id", artifacts["ArtifactID"].is_unique, artifacts["ArtifactID"].tolist())
    for row in artifacts.itertuples(index=False):
        path = root / row.RelativePath
        audit.add("Artifact lock", row.ArtifactID, path.is_file() and path.stat().st_size == row.Bytes and sha256(path) == row.SHA256, row.RelativePath)


def audit_tables(frozen: Path, audit: Audit) -> None:
    items = pd.read_csv(frozen / "checklist-items.tsv", sep="\t", dtype={"ItemNumber": str})
    audit.add("Checklist", "rows", len(items) == 136, len(items))
    storms = items[items["Standard"].eq("STORMS")]
    streams = items[items["Standard"].eq("STREAMS")]
    audit.add("STORMS", "expanded-rows", len(storms) == 69, len(storms))
    audit.add("STORMS", "top-level-items", storms["TopLevelItem"].nunique() == 17, storms["TopLevelItem"].nunique())
    audit.add("STORMS", "preserve-4.10", "4.10" in set(storms["ItemNumber"]), storms.loc[storms["Item"].eq("Replication"), "ItemNumber"].tolist())
    audit.add("STREAMS", "recommendations", len(streams) == 67, len(streams))
    audit.add("STREAMS", "top-level-items", streams["TopLevelItem"].nunique() == 18, streams["TopLevelItem"].nunique())
    audit.add("STREAMS", "sheet", streams["WorkbookSheet"].eq("STREAMS_final").all(), streams["WorkbookSheet"].unique())
    audit.add("Checklist", "six-sections-each", items.groupby("Standard")["ManuscriptSection"].nunique().eq(6).all(), items.groupby("Standard")["ManuscriptSection"].nunique().to_dict())
    audit.add("Checklist", "recommendations-present", items["Recommendation"].str.len().gt(10).all(), items["Recommendation"].str.len().min())

    sections = pd.read_csv(frozen / "checklist-section-counts.tsv", sep="\t")
    expected_storms = {"Abstract": 4, "Introduction": 2, "Methods": 46, "Results": 6, "Discussion": 6, "Other information": 5}
    expected_streams = {"Abstract": 6, "Introduction": 2, "Methods": 42, "Results": 5, "Discussion": 6, "Other information": 6}
    for standard, expected in (("STORMS", expected_storms), ("STREAMS", expected_streams)):
        observed = sections[sections["Standard"].eq(standard)].set_index("ManuscriptSection")["ExpandedRecommendations"].to_dict()
        audit.add("Checklist sections", standard, observed == expected, observed)

    selection = pd.read_csv(frozen / "standard-selection-matrix.tsv", sep="\t")
    audit.add("Selection", "rows", len(selection) == 28, len(selection))
    audit.add("Selection", "one-study-core", selection[selection["Standard"].isin(["STORMS", "STREAMS"])].groupby("Scenario")["RoleCode"].sum().eq(2).all(), selection.groupby("Scenario")["RoleCode"].sum().to_dict())
    human = selection[selection["Scenario"].str.startswith("Human") & selection["Standard"].eq("STORMS")]
    environment = selection[selection["Scenario"].str.startswith("Environmental") & selection["Standard"].eq("STREAMS")]
    audit.add("Selection", "human-storms", human["Role"].eq("Core").all(), human["Role"].tolist())
    audit.add("Selection", "environment-streams", environment["Role"].eq("Core").all(), environment["Role"].tolist())
    audit.add("Selection", "mimag-is-add-on", selection[selection["Standard"].eq("MIMAG")]["Role"].isin(["Add", "Not selected"]).all(), selection[selection["Standard"].eq("MIMAG")]["Role"].tolist())
    audit.add("Selection", "miuvig-is-add-on", selection[selection["Standard"].eq("MIUViG")]["Role"].isin(["Add", "Not selected"]).all(), selection[selection["Standard"].eq("MIUViG")]["Role"].tolist())

    layers = pd.read_csv(frozen / "reporting-layer-map.tsv", sep="\t")
    audit.add("Layers", "rows", len(layers) == 16, len(layers))
    audit.add("Layers", "four-by-four", layers.groupby("Layer").size().eq(4).all() and layers.groupby("Standard").size().eq(4).all(), layers.groupby("Layer").size().to_dict())

    cross = pd.read_csv(frozen / "standards-crosswalk.tsv", sep="\t")
    audit.add("Crosswalk", "rows", len(cross) == 40, len(cross))
    audit.add("Crosswalk", "ten-domains", cross["ReportingDomain"].nunique() == 10, cross["ReportingDomain"].nunique())
    audit.add("Crosswalk", "four-standards-each", cross.groupby("ReportingDomain").size().eq(4).all(), cross.groupby("ReportingDomain").size().to_dict())
    mag_domain = cross[cross["ReportingDomain"].eq("Bacterial/archaeal genome quality")].set_index("Standard")
    virus_domain = cross[cross["ReportingDomain"].eq("Virus genome quality")].set_index("Standard")
    audit.add("Crosswalk", "mimag-core-domain", mag_domain.loc["MIMAG", "Coverage"] == "Core", mag_domain["Coverage"].to_dict())
    audit.add("Crosswalk", "miuvig-core-domain", virus_domain.loc["MIUViG", "Coverage"] == "Core", virus_domain["Coverage"].to_dict())

    mimag_def = pd.read_csv(frozen / "mimag-quality-criteria.tsv", sep="\t")
    audit.add("MIMAG definition", "rows", len(mimag_def) == 10, len(mimag_def))
    audit.add("MIMAG definition", "four-levels", mimag_def["QualityLevel"].nunique() == 4, mimag_def["QualityLevel"].unique())
    hq = mimag_def[mimag_def["QualityLevel"].eq("High-quality draft")]
    audit.add("MIMAG definition", "hq-completeness", hq["Requirement"].str.contains(">90%", regex=False).any(), hq.to_dict("records"))
    audit.add("MIMAG definition", "hq-contamination", hq["Requirement"].str.contains("<5%", regex=False).any(), hq.to_dict("records"))
    audit.add("MIMAG definition", "hq-rna", hq["Requirement"].str.contains("23S, 16S, and 5S", regex=False).any(), hq.to_dict("records"))

    miuvig_fields = pd.read_csv(frozen / "miuvig-mandatory-metadata.tsv", sep="\t")
    categories = pd.read_csv(frozen / "miuvig-quality-categories.tsv", sep="\t")
    audit.add("MIUViG definition", "mandatory-fields", len(miuvig_fields) == 8, len(miuvig_fields))
    audit.add("MIUViG definition", "quality-cells", len(categories) == 9, len(categories))
    audit.add("MIUViG definition", "three-categories", categories["Category"].nunique() == 3, categories["Category"].unique())
    audit.add("MIUViG definition", "finished-review", categories[categories["Category"].eq("Finished genome")]["Requirement"].str.contains("manual review").any(), categories[categories["Category"].eq("Finished genome")].to_dict("records"))

    mag = pd.read_csv(frozen / "article44-mimag-compliance.tsv", sep="\t")
    audit.add("MAG audit", "rows", len(mag) == 23, len(mag))
    audit.add("MAG audit", "unique", mag["MAG"].is_unique, mag["MAG"].tolist())
    audit.add("MAG audit", "hq-count", mag["Article44ExtendedTier"].eq("High quality").sum() == 4, mag["Article44ExtendedTier"].value_counts().to_dict())
    audit.add("MAG audit", "mq-count", mag["Article44ExtendedTier"].eq("Medium quality").sum() == 19, mag["Article44ExtendedTier"].value_counts().to_dict())
    audit.add("MAG audit", "gunc-all-pass", bool_series(mag["GUNCPass"]).all(), mag["GUNCPass"].tolist())
    audit.add("MAG audit", "core-extension-agreement", bool_series(mag["CoreVsExtendedAgreement"]).all(), mag["CoreVsExtendedAgreement"].tolist())
    audit.add("MAG audit", "rRNA-complete-count", bool_series(mag["CompleteRRNASet"]).sum() == 5, bool_series(mag["CompleteRRNASet"]).sum())
    missing_gate = mag["HQMissingGate"].ne("No missing HQ gate")
    audit.add("MAG audit", "missing-gate-recorded", mag.loc[missing_gate, "HQMissingGate"].notna().all() and missing_gate.sum() == 13, mag["HQMissingGate"].value_counts().to_dict())

    virus = pd.read_csv(frozen / "article54-miuvig-compliance.tsv", sep="\t")
    audit.add("UViG audit", "eight-fields", len(virus) == 8, len(virus))
    audit.add("UViG audit", "field-match-source", set(virus["MandatoryMetadata"]) == set(miuvig_fields["MandatoryMetadata"]), virus["MandatoryMetadata"].tolist())
    audit.add("UViG audit", "complete-four", virus["Status"].eq("Complete").sum() == 4, virus["Status"].value_counts().to_dict())
    audit.add("UViG audit", "partial-two", virus["Status"].eq("Partial").sum() == 2, virus["Status"].value_counts().to_dict())
    audit.add("UViG audit", "missing-two", virus["Status"].eq("Missing").sum() == 2, virus["Status"].value_counts().to_dict())
    audit.add("UViG audit", "gaps-explicit", virus.loc[virus["Status"].ne("Complete"), "EvidenceOrGap"].str.len().gt(25).all(), virus.to_dict("records"))

    owners = pd.read_csv(frozen / "field-responsibility-matrix.tsv", sep="\t")
    audit.add("Ownership", "rows", len(owners) == 12, len(owners))
    audit.add("Ownership", "ordered", owners["Order"].tolist() == list(range(1, 13)), owners["Order"].tolist())
    audit.add("Ownership", "all-owned", owners["Owner"].str.len().gt(3).all(), owners["Owner"].tolist())
    audit.add("Ownership", "all-artifacts", owners["RequiredArtifact"].str.len().gt(15).all(), owners["RequiredArtifact"].tolist())

    nas = pd.read_csv(frozen / "not-applicable-ledger.tsv", sep="\t")
    audit.add("N/A", "rows", len(nas) == 8, len(nas))
    audit.add("N/A", "reasons", bool_series(nas["ReasonRecorded"]).all(), nas["ReasonRecorded"].tolist())
    audit.add("N/A", "approvers", nas["Approver"].str.len().gt(3).all(), nas["Approver"].tolist())
    audit.add("N/A", "finished-claims-declined", set(nas.loc[nas["Field"].str.startswith("Finished"), "Disposition"]) == {"Not applicable"}, nas.to_dict("records"))

    ready = pd.read_csv(frozen / "submission-readiness.tsv", sep="\t")
    audit.add("Readiness", "rows", len(ready) == 21, len(ready))
    audit.add("Readiness", "has-missing", ready["Status"].eq("Missing").any(), ready["Status"].value_counts().to_dict())
    audit.add("Readiness", "no-score", not any("score" in value.lower() for value in ready["ReadinessField"]), ready["ReadinessField"].tolist())
    audit.add("Readiness", "next-actions", ready["EvidenceOrNextAction"].str.len().gt(20).all(), ready["EvidenceOrNextAction"].tolist())

    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    expected = {
        "article": 76,
        "analysis_seed": 76001,
        "plot_seed": 20260776,
        "storms_top_level_items": 17,
        "storms_expanded_recommendations": 69,
        "streams_recommendations": 67,
        "streams_manuscript_sections": 6,
        "mimag_quality_levels": 4,
        "miuvig_mandatory_fields": 8,
        "article44_mags": 23,
        "article44_high_quality": 4,
        "article44_medium_quality": 19,
        "article44_core_extended_agreement": 23,
        "article54_uviqs": 46,
        "article54_miuvig_complete_fields": 4,
        "article54_miuvig_partial_fields": 2,
        "article54_miuvig_missing_fields": 2,
        "source_records": 5,
    }
    for key, value in expected.items():
        audit.add("Metric", key, metrics.get(key) == value, metrics.get(key))

    methods = json.loads((frozen / "methods-contract.json").read_text(encoding="utf-8"))
    audit.add("Methods", "four-standards", set(methods["standards"]) == {"STORMS", "STREAMS", "MIMAG", "MIUViG"}, methods["standards"])
    audit.add("Methods", "not-quality-score", "not a study-quality score" in methods["boundary"], methods["boundary"])
    audit.add("Methods", "streams-sheet", "STREAMS_final" in methods["streams_sheet_policy"], methods["streams_sheet_policy"])
    audit.add("Methods", "gunc-extension", "GUNC" in methods["extended_mag_gate"], methods["extended_mag_gate"])


def audit_figures(root: Path, frozen: Path, audit: Audit) -> None:
    figures = root / "figures/article76"
    manifest = json.loads((frozen / "figure-manifest.json").read_text(encoding="utf-8"))
    records = {record["file"]: record for record in manifest["figures"]}
    expected = {f"{stem}.{extension}" for stem in FIGURE_STEMS for extension in ("png", "svg")} | {ANCHOR}
    audit.add("Figures", "manifest-files", set(records) == expected, sorted(records))
    audit.add("Figures", "file-count", len(records) == 19, len(records))
    for filename in sorted(expected):
        path = figures / filename
        record = records.get(filename, {})
        audit.add("Figure file", filename, path.is_file() and path.stat().st_size == record.get("bytes") and sha256(path) == record.get("sha256"), path)
        if filename.endswith(".png"):
            with Image.open(path) as image:
                audit.add("Figure raster", filename, image.width >= 1800 and image.height >= 900, image.size)
        if filename.endswith(".svg"):
            text = path.read_text(encoding="utf-8")
            audit.add("Figure language", filename, re.search(r"[\u3400-\u9fff]", text) is None, "English-only" if re.search(r"[\u3400-\u9fff]", text) is None else "Han text found")


def audit_chapter(chapter: Path, rendered: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    match = re.match(r"^---\n(.*?)\n---\n", text, flags=re.S)
    audit.add("Chapter", "frontmatter", match is not None, chapter)
    metadata = yaml.safe_load(match.group(1)) if match else {}
    audit.add("Chapter", "draft-false", metadata.get("draft") is False, metadata.get("draft"))
    audit.add("Chapter", "eval-true", metadata.get("execute", {}).get("eval") is True, metadata.get("execute"))
    audit.add("Chapter", "freeze-auto", metadata.get("execute", {}).get("freeze") == "auto", metadata.get("execute"))
    audit.add("Chapter", "expected-images", metadata.get("wechat", {}).get("expected_images") == 10, metadata.get("wechat", {}).get("expected_images"))
    for heading in (
        "## 1. 对应论文里的哪张图",
        "## 2. 理论：为什么这么做",
        "## 3. 准备工作",
        "## 4. 可复制代码",
        "## 5. 审计与升级",
        "## 6. 出版级美化",
        "## 7. 常见坑",
        "## 8. 这段 Methods 怎么写",
        "## 9. 换成你自己的数据怎么做",
        "## 参考",
    ):
        audit.add("Chapter heading", heading, heading in text, heading)
    audit.add("Chapter", "inline-theme", all(token in text for token in ("pal_pub <-", "scale_color_pub <-", "scale_fill_pub <-", "theme_pub <-", "save_pub <-")), "five helpers")
    audit.add("Chapter", "no-source-theme", 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text, "inline only")
    audit.add("Chapter", "deterministic", "set.seed(20260776)" in text, "set.seed(20260776)")
    audit.add("Chapter", "methods-versions", all(value in text for value in ("STORMS v1.03", "STREAMS v1.0", "10.1038/nbt.3893", "10.1038/nbt.4306")), "version locks")
    audit.add("Chapter", "not-a-score", "不是研究质量评分" in text, "boundary stated")
    audit.add("Chapter", "rights-boundary", "不属于本仓库 CC BY/MIT" in text, "anchor rights")
    audit.add("Chapter", "four-citations", all(value in text for value in ("@mirzayi2021storms", "@kelliher2025streams", "@bowers2017mimag", "@roux2019miuvig")), "citations")
    for stem in FIGURE_STEMS:
        audit.add("Chapter figure", stem, f"../figures/article76/{stem}.png" in text, stem)
    audit.add("Chapter figure", ANCHOR, f"../figures/article76/{ANCHOR}" in text, ANCHOR)
    audit.add("Chapter", "no-placeholder", not any(term in text for term in ("Planned chapter", "TODO", "待补", "draft: true")), "placeholder scan")

    audit.add("Render", "exists", rendered.is_file(), rendered)
    if rendered.is_file():
        html = rendered.read_text(encoding="utf-8")
        audit.add("Render", "size", rendered.stat().st_size > 180_000, rendered.stat().st_size)
        audit.add("Render", "title", "报告标准" in html and "STORMS" in html and "STREAMS" in html, "title tokens")
        audit.add("Render", "ten-image-references", sum(name in html for name in [ANCHOR, *[f"{stem}.png" for stem in FIGURE_STEMS]]) == 10, "image references")
        audit.add("Render", "no-tofu", "□□" not in html, "tofu scan")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_bundle(frozen, audit)
    audit_sources(root, frozen, audit)
    audit_tables(frozen, audit)
    audit_figures(root, frozen, audit)
    audit_chapter(args.chapter.resolve(), args.rendered_html.resolve(), audit)
    audit.write(args.qa_dir.resolve())


if __name__ == "__main__":
    main()
