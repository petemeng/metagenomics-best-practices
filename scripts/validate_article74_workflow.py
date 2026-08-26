#!/usr/bin/env python3
"""Offline acceptance tests for Article 74's workflow reproducibility packet."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd
from PIL import Image


FIGURE_STEMS = (
    "74-release-lock",
    "74-parameter-precedence",
    "74-execution-unit-matrix",
    "74-profile-separation",
    "74-resume-audit",
    "74-resource-envelope",
    "74-database-lock",
    "74-provenance-bundle",
)
ANCHORS = ("74-mag-metromap-original.png", "74-funcscan-metromap-original.png")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--qa-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass
class Check:
    category: str
    check: str
    status: bool
    detail: str


class Audit:
    def __init__(self) -> None:
        self.rows: list[Check] = []

    def add(self, category: str, check: str, status: bool, detail: object = "") -> None:
        self.rows.append(Check(category, check, bool(status), str(detail)))

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {"Category": row.category, "Check": row.check, "Status": "PASS" if row.status else "FAIL", "Detail": row.detail}
                for row in self.rows
            ]
        )


def audit_bundle(frozen: Path, audit: Audit) -> None:
    checksums = frozen / "file-checksums.sha256"
    audit.add("Bundle", "checksum-file", checksums.is_file(), checksums)
    if checksums.is_file():
        lines = [line for line in checksums.read_text().splitlines() if line]
        audit.add("Bundle", "checksum-count", len(lines) == 53, len(lines))
        for line in lines:
            parts = line.split("  ", 1)
            valid = len(parts) == 2 and bool(re.fullmatch(r"[0-9a-f]{64}", parts[0]))
            audit.add("Checksum", f"format-{len(audit.rows)}", valid, line)
            if valid:
                digest, relative = parts
                path = frozen / relative
                audit.add("Checksum", relative, path.is_file() and sha256(path) == digest, digest)

    manifest_path = frozen / "bundle-manifest.json"
    audit.add("Bundle", "manifest", manifest_path.is_file(), manifest_path)
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        expected = {
            "article": 74,
            "payload_files": 26,
            "selected_source_files": 19,
            "script_files": 5,
            "environment_files": 2,
        }
        for field, value in expected.items():
            audit.add("Bundle", field, manifest.get(field) == value, manifest.get(field))
        contract = manifest.get("contract", "")
        for token in ("nf-core/mag 5.5.0", "nf-core/funcscan 4.0.0", "Nextflow 26.04.0", "resume trace"):
            audit.add("Bundle contract", token, token in contract, contract)


def audit_sources(frozen: Path, audit: Audit) -> None:
    manifest = json.loads((frozen / "source-manifest.json").read_text())
    audit.add("Source", "article", manifest.get("article") == 74, manifest.get("article"))
    audit.add("Source", "snapshot-date", manifest.get("snapshot_date") == "2026-08-24", manifest.get("snapshot_date"))
    downloads = manifest.get("immutable_downloads", {})
    audit.add("Source", "download-count", len(downloads) == 13, len(downloads))
    selected = manifest.get("selected_source_files", {})
    audit.add("Source", "selected-count", len(selected) == 19, len(selected))
    for name, record in selected.items():
        path = frozen / "source" / name
        audit.add("Source file", name, path.is_file() and path.stat().st_size == record["bytes"] and sha256(path) == record["sha256"], record["sha256"])

    releases = json.loads((frozen / "source/release-records.stable.json").read_text())
    expected = {
        "nf-core/mag": ("5.5.0", "56abab5b023ce953c9c43fe21090d156ad0e18af"),
        "nf-core/funcscan": ("4.0.0", "aee3dc965eb0c77267435544dda30da858763913"),
    }
    for name, (release, commit) in expected.items():
        audit.add("Release source", f"{name}-release", releases[name]["release"] == release, releases[name]["release"])
        audit.add("Release source", f"{name}-commit", releases[name]["commit"] == commit, releases[name]["commit"])
    audit.add("Release source", "nextflow-release", releases["Nextflow"]["release"] == "26.04.0", releases["Nextflow"]["release"])

    mag_config = (frozen / "source/mag-nextflow.config").read_text()
    for token in (
        "nextflowVersion = '!>=26.04.0'",
        "version         = '5.5.0'",
        "id 'nf-schema@2.7.2'",
        "gtdbtk_r232_data.tar.gz",
        "checkm2_db_version                   = 14897628",
        "run_checkm2                          = false",
        "run_gunc                             = false",
        "metabat_rng_seed                     = 1",
        "semibin_rng_seed                     = 1",
    ):
        audit.add("Pipeline source", token, token in mag_config, token)

    func_config = (frozen / "source/funcscan-nextflow.config").read_text()
    for token in ("run_amp_screening", "run_arg_screening", "run_bgc_screening", "run_cazyme_screening"):
        audit.add("Funcscan source", token, token in func_config, token)


def audit_release_and_databases(frozen: Path, audit: Audit) -> None:
    release = pd.read_csv(frozen / "release-lock.tsv", sep="\t", dtype=str)
    audit.add("Release", "row-count", len(release) == 3, len(release))
    audit.add("Release", "components", set(release["Component"]) == {"nf-core/mag", "nf-core/funcscan", "Nextflow"}, release["Component"].tolist())
    indexed = release.set_index("Component")
    locks = {
        "nf-core/mag": ("5.5.0", "871528217457f88677e1dcd5a3822bed9d1234c3f07366652cb5f49b7a2700b3"),
        "nf-core/funcscan": ("4.0.0", "1d387cdd7265529ac6cebceda890334578a35d71755378ccdbe32fd9b5c8f224"),
        "Nextflow": ("26.04.0", "5e2b4a354b4d7634d7211b71417d61606878fb49e9b224b50ded6e2c69114870"),
    }
    for name, (version, digest) in locks.items():
        audit.add("Release", f"{name}-version", indexed.loc[name, "Release"] == version, indexed.loc[name, "Release"])
        audit.add("Release", f"{name}-sha256", indexed.loc[name, "SHA256"] == digest, indexed.loc[name, "SHA256"])

    db = pd.read_csv(frozen / "database-lock.tsv", sep="\t", dtype=str)
    audit.add("Database", "row-count", len(db) == 3, len(db))
    expected = {
        "GTDB-Tk reference package": ("R11-RS232", "60806405195", "25a59e0352b1fd150c589f56559767d4"),
        "CheckM2 reference database": ("Zenodo 14897628; dataset version 3", "1735095710", "07c10655620843b517d0df0c160d911f"),
        "GUNC ProGenomes": ("ProGenomes 2.1; GUNC 1.1.0 download endpoint", "7185710760", "bc93a855e0760aad5c4e5f2d0e26da46"),
    }
    indexed_db = db.set_index("Database")
    for name, (release_name, size, checksum) in expected.items():
        audit.add("Database", f"{name}-release", indexed_db.loc[name, "Release"] == release_name, indexed_db.loc[name, "Release"])
        audit.add("Database", f"{name}-bytes", indexed_db.loc[name, "Bytes"] == size, indexed_db.loc[name, "Bytes"])
        audit.add("Database", f"{name}-checksum", indexed_db.loc[name, "Checksum"] == checksum, indexed_db.loc[name, "Checksum"])
        audit.add("Database", f"{name}-https", indexed_db.loc[name, "SourceURL"].startswith("https://"), indexed_db.loc[name, "SourceURL"])
        audit.add("Database", f"{name}-local-path", indexed_db.loc[name, "LocalPath"].startswith("/shared/db/"), indexed_db.loc[name, "LocalPath"])


def audit_input_and_parameters(root: Path, frozen: Path, audit: Audit) -> None:
    inputs = pd.read_csv(frozen / "input-manifest.tsv", sep="\t", dtype=str)
    audit.add("Input", "rows", len(inputs) == 2, len(inputs))
    audit.add("Input", "mates", inputs["Mate"].tolist() == ["R1", "R2"], inputs["Mate"].tolist())
    audit.add("Input", "records", inputs["Records"].tolist() == ["99991", "99991"], inputs["Records"].tolist())
    audit.add("Input", "run", inputs["Run"].eq("ERR9765746").all(), inputs["Run"].unique())
    expected_hashes = {
        "R1": "ce438de2814a6e605cc24b00e53b697d018f325bc2a9694cece5be158ff32101",
        "R2": "207d080cd5dc98bd10fd4af8c492fb3f2000f73e1458a9d1d77e57347f4fe459",
    }
    for row in inputs.itertuples(index=False):
        path = root / row.RelativePath
        audit.add("Input file", row.Mate, path.is_file() and path.stat().st_size == int(row.Bytes) and sha256(path) == row.SHA256 == expected_hashes[row.Mate], path)

    sheet = pd.read_csv(frozen / "samplesheet.csv")
    audit.add("Samplesheet", "header", sheet.columns.tolist() == ["sample", "run", "group", "short_reads_1", "short_reads_2", "short_reads_platform"], sheet.columns.tolist())
    audit.add("Samplesheet", "row", len(sheet) == 1, len(sheet))
    audit.add("Samplesheet", "platform", sheet.iloc[0]["short_reads_platform"] == "ILLUMINA", sheet.iloc[0]["short_reads_platform"])

    params = (frozen / "params.publication.yml").read_text()
    required = (
        "coassemble_group: false",
        "binning_map_mode: own",
        "skip_megahit: false",
        "skip_metabat2: false",
        "metabat_rng_seed: 1",
        "semibin_rng_seed: 1",
        "run_checkm2: true",
        "run_gunc: true",
        "gtdbtk_r232_data.tar.gz",
        "zenodo-14897628",
        "progenomes2.1",
    )
    for token in required:
        audit.add("Parameters", token, token in params, token)

    hpc = (frozen / "hpc.slurm.config").read_text()
    for token in ("executor = 'slurm'", "queue = 'compute'", "memory: 220.GB", "apptainer", "cacheDir = '/shared/containers/nf-core-mag-5.5.0'", "docker.enabled = false", "singularity.enabled = false", "conda.enabled = false"):
        audit.add("HPC config", token, token in hpc, token)

    command = (frozen / "production-command.sh").read_text()
    for token in ('PIPELINE_RELEASE="5.5.0"', 'NEXTFLOW_VERSION="26.04.0"', 'export NXF_VER="${NEXTFLOW_VERSION}"', '-r "${PIPELINE_RELEASE}"', "-profile apptainer", "-params-file params.publication.yml", "-c hpc.slurm.config", "-resume"):
        audit.add("Production command", token, token in command, token)

    parsed = pd.read_csv(frozen / "profile-parse-evidence.tsv", sep="\t", dtype=str)
    audit.add("Config parse", "rows", len(parsed) == 12, len(parsed))
    indexed = parsed.set_index("EffectiveConfigKey")["ParsedValue"]
    for key, value in {"manifest.version": "5.5.0", "process.executor": "slurm", "apptainer.enabled": "true", "docker.enabled": "false", "singularity.enabled": "false", "conda.enabled": "false"}.items():
        audit.add("Config parse", key, indexed.get(key) == value, indexed.get(key))


def audit_runtime(frozen: Path, audit: Audit) -> None:
    trace = pd.read_csv(frozen / "stub-runtime-trace.tsv", sep="\t", dtype={"Hash": str})
    audit.add("Runtime", "rows", len(trace) == 10, len(trace))
    first = trace[trace["Run"].eq("first-success")]
    resumed = trace[trace["Run"].eq("resume")]
    audit.add("Runtime", "first-tasks", len(first) == 5, len(first))
    audit.add("Runtime", "first-completed", first["Status"].eq("COMPLETED").all(), first["Status"].value_counts().to_dict())
    audit.add("Runtime", "resume-tasks", len(resumed) == 5, len(resumed))
    audit.add("Runtime", "resume-cached", int(resumed["Status"].eq("CACHED").sum()) == 4, resumed["Status"].value_counts().to_dict())
    rebuilt = resumed[~resumed["Status"].eq("CACHED")]
    audit.add("Runtime", "resume-multiqc-rebuilt", len(rebuilt) == 1 and rebuilt.iloc[0]["Process"].endswith("MULTIQC (mag)"), rebuilt["Process"].tolist())
    first_hashes = set(first["Hash"])
    for row in resumed[resumed["Status"].eq("CACHED")].itertuples(index=False):
        audit.add("Runtime cache", row.Process, row.Hash in first_hashes, row.Hash)
    audit.add("Runtime", "exit-zero", trace["Exit"].eq(0).all(), trace["Exit"].value_counts().to_dict())

    summary = pd.read_csv(frozen / "runtime-summary.tsv", sep="\t", dtype=str)
    audit.add("Runtime summary", "attempts", summary["Attempt"].tolist() == ["offline-schema", "first-success", "resume"], summary["Attempt"].tolist())
    audit.add("Runtime summary", "offline-boundary", summary.iloc[0]["Outcome"] == "FAILED_AS_EXPECTED" and "NXF_OFFLINE" in summary.iloc[0]["Interpretation"], summary.iloc[0].to_dict())

    scope = pd.read_csv(frozen / "stub-scope.tsv", sep="\t")
    audit.add("Boundary", "rows", len(scope) == 7, len(scope))
    supported = scope["SupportedByLocalRun"].astype(bool)
    audit.add("Boundary", "supported-count", int(supported.sum()) == 4, supported.value_counts().to_dict())
    biological = scope[scope["Claim"].str.contains("biologically valid", regex=False)]
    audit.add("Boundary", "no-biological-claim", len(biological) == 1 and not bool(biological.iloc[0]["SupportedByLocalRun"]), biological.to_dict("records"))


def audit_tables(frozen: Path, audit: Audit) -> None:
    defaults = pd.read_csv(frozen / "pipeline-defaults-audit.tsv", sep="\t", dtype=str)
    audit.add("Defaults", "rows", len(defaults) == 13, len(defaults))
    indexed = defaults.set_index("Parameter")
    audit.add("Defaults", "checkm2-off-to-on", indexed.loc["run_checkm2", "PipelineDefault5.5.0"] == "false" and indexed.loc["run_checkm2", "PublicationOverride"] == "true", indexed.loc["run_checkm2"].to_dict())
    audit.add("Defaults", "gunc-off-to-on", indexed.loc["run_gunc", "PipelineDefault5.5.0"] == "false" and indexed.loc["run_gunc", "PublicationOverride"] == "true", indexed.loc["run_gunc"].to_dict())

    precedence = pd.read_csv(frozen / "parameter-precedence.tsv", sep="\t")
    audit.add("Precedence", "ranks", precedence["Rank"].tolist() == [1, 2, 3, 4], precedence["Rank"].tolist())
    audit.add("Precedence", "params-file", precedence.iloc[2]["Layer"] == "Parameter file", precedence.iloc[2].to_dict())

    units = pd.read_csv(frozen / "execution-unit-matrix.tsv", sep="\t")
    audit.add("Units", "rows", len(units) == 4, len(units))
    primary = units.iloc[0]
    audit.add("Units", "primary-own", not bool(primary["coassemble_group"]) and primary["binning_map_mode"] == "own" and primary["AssemblyUnit"] == "sample", primary.to_dict())

    funcscan = pd.read_csv(frozen / "funcscan-branch-contract.tsv", sep="\t")
    audit.add("Funcscan", "branches", set(funcscan["Branch"]) == {"ARG", "AMP", "BGC", "CAZyme/CGC"}, funcscan["Branch"].tolist())
    audit.add("Funcscan", "flags", set(funcscan["ActivationFlag"]) == {"run_arg_screening", "run_amp_screening", "run_bgc_screening", "run_cazyme_screening"}, funcscan["ActivationFlag"].tolist())

    provenance = pd.read_csv(frozen / "provenance-bundle.tsv", sep="\t")
    present = provenance["PresentInLocalPacket"].astype(bool)
    audit.add("Provenance", "rows", len(provenance) == 9, len(provenance))
    audit.add("Provenance", "pending-count", int((~present).sum()) == 2, provenance.loc[~present, "Layer"].tolist())
    audit.add("Provenance", "pending-items", set(provenance.loc[~present, "Layer"]) == {"Container identity", "Scientific acceptance"}, provenance.loc[~present, "Layer"].tolist())


def audit_figures(root: Path, frozen: Path, audit: Audit) -> None:
    manifest = json.loads((frozen / "figure-manifest.json").read_text())
    audit.add("Figure", "article", manifest.get("article") == 74, manifest.get("article"))
    audit.add("Figure", "plot-seed", manifest.get("plot_seed") == 20260774, manifest.get("plot_seed"))
    records = manifest.get("files", {})
    audit.add("Figure", "file-count", len(records) == 18, len(records))
    for stem in FIGURE_STEMS:
        for suffix in (".png", ".svg"):
            name = stem + suffix
            path = root / "figures/article74" / name
            record = records.get(name, {})
            audit.add("Figure file", name, path.is_file() and path.stat().st_size == record.get("bytes") and sha256(path) == record.get("sha256"), path)
            if suffix == ".png" and path.is_file():
                with Image.open(path) as image:
                    audit.add("Figure resolution", name, image.width >= 2400 and image.height >= 1400, image.size)
            if suffix == ".svg" and path.is_file():
                text = path.read_text(errors="ignore")
                audit.add("Figure language", name, not bool(re.search(r"[\u4e00-\u9fff]", text)), "English-only vector text")
    for name in ANCHORS:
        path = root / "figures/article74" / name
        source_name = name.replace("74-", "", 1)
        top = frozen / source_name
        source = frozen / "source" / source_name
        audit.add("Anchor", name, path.is_file() and top.is_file() and source.is_file() and sha256(path) == sha256(top) == sha256(source), path)


def audit_chapter(root: Path, audit: Audit) -> None:
    chapter = root / "chapters/74-nfcore-mag-workflows.qmd"
    audit.add("Chapter", "file", chapter.is_file(), chapter)
    if not chapter.is_file():
        return
    text = chapter.read_text()
    front = text.split("---", 2)[1] if text.startswith("---") else ""
    audit.add("Chapter", "not-draft", "draft: true" not in front, front)
    audit.add("Chapter", "eval-false", "eval: false" in front, front)
    headings = (
        "## 这一步对应论文里的哪张图",
        "## 理论：工作流锁住的是分析契约",
        "## 准备工作",
        "## 可复制代码",
        "## 审计与升级",
        "## 出版级美化",
        "## 常见坑",
        "## 这段 Methods 怎么写",
        "## 换成你自己的数据怎么做",
        "## 参考",
    )
    positions = []
    for heading in headings:
        index = text.find(heading)
        positions.append(index)
        audit.add("Chapter structure", heading, index >= 0, index)
    audit.add("Chapter structure", "order", positions == sorted(positions) and all(index >= 0 for index in positions), positions)
    for stem in FIGURE_STEMS:
        audit.add("Chapter figure", stem, f"../figures/article74/{stem}.png" in text, stem)
    for anchor in ANCHORS:
        audit.add("Chapter figure", anchor, f"../figures/article74/{anchor}" in text, anchor)
    for token in ("-r 5.5.0", "-profile apptainer", "-params-file", "hpc.slurm.config", "GTDB R232", "14897628", "ProGenomes 2.1", "nf-core/funcscan", "4.0.0", "ARG", "AMP", "BGC", "-resume", "stub-run"):
        audit.add("Chapter content", token, token in text, token)
    audit.add("Chapter", "no-source-theme", 'source("R/theme_pub.R")' not in text, "inline prep required")
    for phrase in ("本篇可独立跑通", "这体现全系列", "作者代码通常长这样", "即本文"):
        audit.add("Chapter prose", f"forbidden-{phrase}", phrase not in text, phrase)
    audit.add("Chapter", "citations", text.count("@") >= 8, text.count("@"))
    audit.add("Chapter", "code-blocks", text.count("```") >= 24, text.count("```"))


class RenderParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.images: list[str] = []
        self.citations: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if tag == "img":
            self.images.append(values.get("src", "") or "")
        if tag == "a" and values.get("role") == "doc-biblioref":
            self.citations.append(values.get("href", "") or "")


def audit_render(root: Path, audit: Audit) -> None:
    html = root / "_site/chapters/74-nfcore-mag-workflows.html"
    audit.add("Render", "html", html.is_file() and html.stat().st_size > 250_000, html)
    if not html.is_file():
        return
    text = html.read_text()
    parser = RenderParser()
    parser.feed(text)
    audit.add("Render", "image-count", len(parser.images) == 10 and len(set(parser.images)) == 10, parser.images)
    for src in parser.images:
        if src.startswith(("http:", "https:", "data:")):
            continue
        path = (html.parent / src).resolve()
        audit.add("Render image", src, path.is_file(), path)
    audit.add("Render", "citation-count", len(parser.citations) >= 20 and len(set(parser.citations)) >= 10, f"{len(parser.citations)}/{len(set(parser.citations))}")
    audit.add("Render", "no-unresolved-crossref", "???" not in text, "???")
    audit.add("Render", "title", "第 74 篇 · 用 nf-core/mag" in text, "title")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    qa = args.qa_dir.resolve()
    qa.mkdir(parents=True, exist_ok=True)
    audit = Audit()
    audit_bundle(frozen, audit)
    audit_sources(frozen, audit)
    audit_release_and_databases(frozen, audit)
    audit_input_and_parameters(root, frozen, audit)
    audit_runtime(frozen, audit)
    audit_tables(frozen, audit)
    audit_figures(root, frozen, audit)
    audit_chapter(root, audit)
    audit_render(root, audit)
    frame = audit.frame()
    frame.to_csv(qa / "checks.tsv", sep="\t", index=False, lineterminator="\n")
    failures = frame[frame["Status"].eq("FAIL")]
    report = {
        "article": 74,
        "status": "passed" if failures.empty else "failed",
        "checks": len(frame),
        "passed": int(frame["Status"].eq("PASS").sum()),
        "failed": len(failures),
        "failed_checks": failures[["Category", "Check", "Detail"]].to_dict("records"),
    }
    (qa / "qa_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    raise SystemExit(0 if failures.empty else 1)


if __name__ == "__main__":
    main()
