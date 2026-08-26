#!/usr/bin/env python3
"""Validate Article 77's frozen packet, figures, chapter and rendered HTML."""

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
    "77-repository-routing",
    "77-accession-graph",
    "77-release-layers",
    "77-artifact-readiness",
    "77-identifier-state",
    "77-container-digest",
    "77-database-manifest",
    "77-release-gates",
    "77-availability-statement",
)
ANCHOR = "77-tenhoopen-figure2-original.jpg"
EXPECTED_SOURCE_SHA = {
    "PMC5737865-fulltext.xml": "8952fa2f6bc3c8960e67a71be152a16d6700aa81405e8a44c984c6cb4b7f0d90",
    "gix047fig2.jpg": "1a7ebe7ca72de90d9de9a6f7d5daec8b641a73fece7b3458d40a46fc45237477",
}


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


def booleans(series: pd.Series) -> pd.Series:
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
            "article": 77,
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
    audit.add("Bundle", "article", manifest.get("article") == 77, manifest.get("article"))
    audit.add("Bundle", "payload-files", manifest.get("payload_files") == 21, manifest.get("payload_files"))
    audit.add("Bundle", "script-files", manifest.get("script_files") == 5, manifest.get("script_files"))
    audit.add("Bundle", "environment-files", manifest.get("environment_files") == 2, manifest.get("environment_files"))
    audit.add("Bundle", "not-publication", "not an SRA/ENA" in manifest.get("publication_boundary", ""), manifest.get("publication_boundary"))
    audit.add("Bundle", "pending-not-identifier", "No pending token" in manifest.get("publication_boundary", ""), manifest.get("publication_boundary"))
    audit.add("Bundle", "ownership-boundary", "third-party" in manifest.get("ownership_boundary", ""), manifest.get("ownership_boundary"))

    lines = [line for line in (frozen / "file-checksums.sha256").read_text().splitlines() if line.strip()]
    audit.add("Checksum", "record-count", len(lines) == 29, len(lines))
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
    audit.add("Sources", "article", sources.get("article") == 77, sources.get("article"))
    audit.add("Sources", "snapshot", sources.get("snapshot_date") == "2026-08-23", sources.get("snapshot_date"))
    records = {row["File"]: row for row in sources.get("records", [])}
    audit.add("Sources", "two-byte-locks", set(records) == set(EXPECTED_SOURCE_SHA), sorted(records))
    for filename, digest in EXPECTED_SOURCE_SHA.items():
        audit.add("Source lock", filename, records.get(filename, {}).get("SHA256") == digest, records.get(filename, {}).get("SHA256"))
    anchor = sources["anchor"]
    audit.add("Anchor", "license", anchor.get("license") == "CC BY 4.0", anchor.get("license"))
    audit.add("Anchor", "attribution", "10.1093/gigascience/gix047" in anchor.get("attribution", ""), anchor.get("attribution"))
    audit.add("Anchor", "dimensions-recorded", anchor.get("dimensions") == [767, 508], anchor.get("dimensions"))
    frozen_anchor = frozen / "tenhoopen-figure2-original.jpg"
    audit.add("Anchor", "frozen-sha", frozen_anchor.is_file() and sha256(frozen_anchor) == EXPECTED_SOURCE_SHA["gix047fig2.jpg"], frozen_anchor)
    with Image.open(frozen_anchor) as image:
        audit.add("Anchor", "image-dimensions", image.size == (767, 508), image.size)
    audit.add("Policy snapshot", "live-pages-mutable", sources["policy_snapshot"]["mutable_pages"] == 17, sources["policy_snapshot"])
    audit.add("Policy snapshot", "recheck-rule", "Recheck" in sources["policy_snapshot"]["rule"], sources["policy_snapshot"]["rule"])
    audit.add("Download policy", "zip-wrapper-not-locked", "ZIP wrapper" in sources["download_manifest"]["archive_policy"], sources["download_manifest"]["archive_policy"])

    artifacts = pd.read_csv(frozen / "source-artifact-manifest.tsv", sep="\t")
    audit.add("Source artifacts", "count", len(artifacts) == 15, len(artifacts))
    audit.add("Source artifacts", "unique", artifacts["ArtifactID"].is_unique, artifacts["ArtifactID"].tolist())
    audit.add("Source artifacts", "relative-paths", ~artifacts["RelativePath"].str.startswith("/").any(), artifacts["RelativePath"].tolist())
    for row in artifacts.itertuples(index=False):
        path = root / row.RelativePath
        audit.add("Artifact lock", row.ArtifactID, path.is_file() and path.stat().st_size == row.Bytes and sha256(path) == row.SHA256, row.RelativePath)


def audit_policy_and_routing(frozen: Path, audit: Audit) -> None:
    sources = pd.read_csv(frozen / "policy-source-registry.tsv", sep="\t")
    audit.add("Policy", "sources", len(sources) == 18, len(sources))
    audit.add("Policy", "unique-id", sources["SourceID"].is_unique, sources["SourceID"].tolist())
    audit.add("Policy", "snapshot-date", sources["RetrievedDate"].eq("2026-08-23").all(), sources["RetrievedDate"].unique())
    audit.add("Policy", "https-or-doi", sources["URL"].str.startswith("https://").all(), sources["URL"].tolist())
    audit.add("Policy", "core-providers", {"NCBI", "ENA", "Zenodo", "GitHub", "Docker", "Springer Nature"}.issubset(set(sources["Provider"])), sorted(sources["Provider"].unique()))
    live = sources[~sources["SourceID"].eq("LIFECYCLE")]
    audit.add("Policy", "live-recheck", live["RecheckAtSubmission"].eq("Yes").all(), live["RecheckAtSubmission"].value_counts().to_dict())

    assertions = pd.read_csv(frozen / "policy-assertions.tsv", sep="\t")
    audit.add("Assertions", "rows", len(assertions) == 22, len(assertions))
    audit.add("Assertions", "unique-id", assertions["AssertionID"].is_unique, assertions["AssertionID"].tolist())
    audit.add("Assertions", "source-links", assertions["SourceIDs"].str.len().gt(3).all(), assertions["SourceIDs"].tolist())
    for token in ("raw or minimally processed", "not the accession", "validation and submission", "tag is mutable", "generated and reused", "not resubmitted"):
        audit.add("Assertion content", token, assertions["Assertion"].str.contains(token, case=False, regex=False).any(), token)

    routes = pd.read_csv(frozen / "repository-routing-matrix.tsv", sep="\t")
    audit.add("Routing", "rows", len(routes) == 15, len(routes))
    audit.add("Routing", "ordered", routes["RouteOrder"].tolist() == list(range(1, 16)), routes["RouteOrder"].tolist())
    audit.add("Routing", "unique-assets", routes["AssetClass"].is_unique, routes["AssetClass"].tolist())
    expected = {"Raw reads", "MAG", "UViG / vOTU", "Gene catalog", "Workflow source", "Container image", "Database manifest"}
    audit.add("Routing", "required-assets", expected.issubset(set(routes["AssetClass"])), sorted(routes["AssetClass"]))
    raw = routes.loc[routes["AssetClass"].eq("Raw reads")].iloc[0]
    mag = routes.loc[routes["AssetClass"].eq("MAG")].iloc[0]
    code = routes.loc[routes["AssetClass"].eq("Workflow source")].iloc[0]
    container = routes.loc[routes["AssetClass"].eq("Container image")].iloc[0]
    audit.add("Routing", "raw-to-archive", "SRA / ENA" in raw["Destination"], raw["Destination"])
    audit.add("Routing", "mag-to-insdc", "GenBank / ENA" in mag["Destination"], mag["Destination"])
    audit.add("Routing", "code-dual-preservation", "GitHub + Zenodo" in code["Destination"], code["Destination"])
    audit.add("Routing", "container-oci", container["Destination"] == "OCI registry", container["Destination"])


def audit_identifiers_and_objects(frozen: Path, audit: Audit) -> None:
    objects = pd.read_csv(frozen / "object-relationship.tsv", sep="\t")
    audit.add("Objects", "rows", len(objects) == 14, len(objects))
    audit.add("Objects", "unique-id", objects["ObjectID"].is_unique, objects["ObjectID"].tolist())
    audit.add("Objects", "study", ((objects["ObjectType"] == "Study") & (objects["Identifier"] == "PRJEB52977")).sum() == 1, objects.to_dict("records"))
    audit.add("Objects", "two-samples", objects["ObjectType"].eq("Sample").sum() == 2, objects["ObjectType"].value_counts().to_dict())
    audit.add("Objects", "two-runs", set(objects.loc[objects["ObjectType"].eq("Run"), "Identifier"]) == {"ERR9765746", "ERR9765747"}, objects.loc[objects["ObjectType"].eq("Run"), "Identifier"].tolist())
    experiments = objects[objects["ObjectType"].eq("Experiment")]
    audit.add("Objects", "experiments-not-invented", experiments["Identifier"].eq("EXISTING_NOT_CAPTURED").all(), experiments["Identifier"].tolist())
    audit.add("Objects", "derived-parent-links", objects.loc[objects["ObjectType"].isin(["Primary assembly", "Gene catalog", "MAG collection"]), "Parent"].str.len().gt(5).all(), objects.to_dict("records"))
    audit.add("Objects", "uvig-independent", "INDEPENDENT_CHECKV_FIXTURE" in set(objects["Parent"]), objects["Parent"].tolist())

    ids = pd.read_csv(frozen / "identifier-registry.tsv", sep="\t")
    audit.add("Identifiers", "rows", len(ids) == 15, len(ids))
    audit.add("Identifiers", "existing-five", ids["State"].eq("EXISTING_THIRD_PARTY").sum() == 5, ids["State"].value_counts().to_dict())
    audit.add("Identifiers", "one-local", ids["State"].eq("LOCAL_ONLY").sum() == 1, ids["State"].value_counts().to_dict())
    audit.add("Identifiers", "one-resolve", ids["State"].eq("RESOLVE_BEFORE_RELEASE").sum() == 1, ids["State"].value_counts().to_dict())
    audit.add("Identifiers", "eight-blocked", ids["State"].eq("BLOCKED").sum() == 8, ids["State"].value_counts().to_dict())
    audit.add("Identifiers", "pending-token", ids.loc[ids["State"].eq("BLOCKED"), "Value"].eq("PENDING_NOT_INVENTED").all(), ids.loc[ids["State"].eq("BLOCKED"), "Value"].tolist())
    audit.add("Identifiers", "no-fake-doi", ~ids["Value"].str.match(r"10\.\d{4,9}/").any(), ids["Value"].tolist())
    audit.add("Identifiers", "no-fake-ncbi-project", ~ids.loc[ids["State"].ne("EXISTING_THIRD_PARTY"), "Value"].str.match(r"PRJNA\d+").any(), ids["Value"].tolist())
    audit.add("Identifiers", "state-codes", ids["StateCode"].notna().all(), ids[["State", "StateCode"]].drop_duplicates().to_dict("records"))


def audit_release_packet(frozen: Path, audit: Audit) -> None:
    readiness = pd.read_csv(frozen / "release-readiness.tsv", sep="\t")
    audit.add("Readiness", "rows", len(readiness) == 15, len(readiness))
    audit.add("Readiness", "all-not-external", (~booleans(readiness["ExternalReady"])).all(), readiness["ExternalReady"].tolist())
    audit.add("Readiness", "external-blocked", readiness.loc[readiness["ReleaseObject"].eq("External release"), "Status"].eq("Blocked").all(), readiness.to_dict("records"))
    audit.add("Readiness", "two-existing", readiness["Status"].eq("Existing third-party").sum() == 2, readiness["Status"].value_counts().to_dict())
    audit.add("Readiness", "container-missing", readiness.loc[readiness["ReleaseObject"].eq("Container image"), "Status"].eq("Missing").all(), readiness.to_dict("records"))
    audit.add("Readiness", "next-actions", readiness["EvidenceOrNextAction"].str.len().gt(20).all(), readiness["EvidenceOrNextAction"].tolist())

    gates = pd.read_csv(frozen / "release-gate-ledger.tsv", sep="\t")
    audit.add("Gates", "rows", len(gates) == 18, len(gates))
    audit.add("Gates", "ordered", gates["GateOrder"].tolist() == list(range(1, 19)), gates["GateOrder"].tolist())
    audit.add("Gates", "blocked-many", gates["Status"].eq("Blocked").sum() >= 10, gates["Status"].value_counts().to_dict())
    audit.add("Gates", "three-pass", gates["Status"].eq("Pass").sum() == 3, gates["Status"].value_counts().to_dict())
    audit.add("Gates", "owners", gates["Owner"].str.len().gt(5).all(), gates["Owner"].tolist())
    for gate in ("Ownership and non-duplicate", "MAG manual review", "UViG MIUViG record", "Accession resolvability", "Tagged code release", "Code archive DOI", "Container digest"):
        audit.add("Mandatory gate", gate, gates.loc[gates["Gate"].eq(gate), "Status"].eq("Blocked").all(), gates.loc[gates["Gate"].eq(gate)].to_dict("records"))

    containers = pd.read_csv(frozen / "container-ledger.tsv", sep="\t")
    audit.add("Container", "rows", len(containers) == 8, len(containers))
    audit.add("Container", "recipe-local", containers.loc[containers["ContainerField"].eq("Build recipe"), "State"].eq("LOCAL_ONLY").all(), containers.to_dict("records"))
    audit.add("Container", "tag-blocked", containers.loc[containers["ContainerField"].eq("Image tag"), "State"].eq("BLOCKED").all(), containers.to_dict("records"))
    audit.add("Container", "two-digests", {"Multi-platform index digest", "Platform manifest digest"}.issubset(set(containers["ContainerField"])), containers["ContainerField"].tolist())
    digest_rows = containers[containers["ContainerField"].str.contains("digest", case=False)]
    audit.add("Container", "no-fake-digest", digest_rows["Value"].eq("PENDING_NOT_INVENTED").all(), digest_rows["Value"].tolist())

    databases = pd.read_csv(frozen / "database-manifest-public.tsv", sep="\t")
    audit.add("Database", "rows", len(databases) == 3, len(databases))
    audit.add("Database", "no-local-path-column", "LocalPath" not in databases.columns, databases.columns.tolist())
    audit.add("Database", "no-shared-path", not databases.astype(str).apply(lambda column: column.str.contains("/shared/").any()).any(), "public fields")
    audit.add("Database", "required-fields", {"Release", "Artifact", "ChecksumType", "Checksum", "SourceURL"}.issubset(databases.columns), databases.columns.tolist())
    audit.add("Database", "all-checksums", databases["Checksum"].str.len().ge(32).all(), databases["Checksum"].tolist())
    audit.add("Database", "all-sources", databases["SourceURL"].str.startswith("https://").all(), databases["SourceURL"].tolist())
    audit.add("Database", "version-locks", all(token in " ".join(databases["Release"].astype(str)) for token in ("R11-RS232", "version 3", "ProGenomes 2.1")), databases["Release"].tolist())

    release = pd.read_csv(frozen / "release-artifact-manifest.tsv", sep="\t")
    audit.add("Release artifacts", "rows", len(release) == 11, len(release))
    audit.add("Release artifacts", "four-fastqs", release["AssetClass"].eq("Raw read").sum() == 4, release["AssetClass"].value_counts().to_dict())
    audit.add("Release artifacts", "raw-existing", release.loc[release["AssetClass"].eq("Raw read"), "State"].eq("EXISTING_THIRD_PARTY").all(), release.to_dict("records"))
    audit.add("Release artifacts", "gene-hashes", release.loc[release["AssetClass"].str.startswith("Gene catalog"), "Checksum"].str.len().eq(64).all(), release.to_dict("records"))
    audit.add("Release artifacts", "not-staged-explicit", release.loc[release["Artifact"].str.contains("megahit-mix"), "Bytes"].astype(str).eq("NOT_STAGED").all(), release.to_dict("records"))
    audit.add("Release artifacts", "external-blockers", set(release.loc[release["Artifact"].isin(["MAG-collection", "UViG-collection", "Workflow-release", "Container-image"]), "State"]) == {"BLOCKED"}, release.to_dict("records"))

    packet = pd.read_csv(frozen / "package-index.tsv", sep="\t")
    audit.add("Package", "rows", len(packet) == 18, len(packet))
    audit.add("Package", "ordered", packet["Order"].tolist() == list(range(1, 19)), packet["Order"].tolist())
    audit.add("Package", "core-files", {"README", "Checksums", "CITATION", "License", "Availability statements"}.issubset(set(packet["PackageObject"])), packet["PackageObject"].tolist())


def audit_availability(frozen: Path, audit: Audit) -> None:
    data = pd.read_csv(frozen / "data-availability-components.tsv", sep="\t")
    code = pd.read_csv(frozen / "code-availability-components.tsv", sep="\t")
    audit.add("Data Availability", "rows", len(data) == 9, len(data))
    audit.add("Code Availability", "rows", len(code) == 8, len(code))
    audit.add("Data Availability", "reused-identifiers", all(token in ";".join(data["IdentifierOrEvidence"].astype(str)) for token in ("PRJEB52977", "SAMEA14435832", "ERR9765746")), data.to_dict("records"))
    audit.add("Data Availability", "generated-pending", data.loc[data["Component"].isin(["Primary assembly", "MAGs", "UViGs", "Gene catalog"]), "IdentifierOrEvidence"].eq("PENDING_NOT_INVENTED").all(), data.to_dict("records"))
    audit.add("Data Availability", "restrictions-clause", "Restrictions" in set(data["Component"]), data["Component"].tolist())
    audit.add("Code Availability", "repo-tag-doi", {"Repository URL", "Release tag", "Version DOI"}.issubset(set(code["Component"])), code["Component"].tolist())
    audit.add("Code Availability", "digest", "Container digest" in set(code["Component"]), code["Component"].tolist())
    audit.add("Code Availability", "citation-missing", code.loc[code["Component"].eq("CITATION.cff"), "Status"].eq("Missing").all(), code.to_dict("records"))

    statement = (frozen / "availability-statements.md").read_text(encoding="utf-8")
    for token in ("PRJEB52977", "SAMEA14435832", "SAMEA14435833", "ERR9765746", "ERR9765747", "have not been deposited", "must not be described", "Production Data Availability template", "Production Code Availability template", "OCI_INDEX_DIGEST", "OCI_PLATFORM_DIGEST"):
        audit.add("Statement", token, token in statement, token)
    audit.add("Statement", "separate-headings", "## Worked-example Data Availability" in statement and "## Worked-example Code Availability" in statement, "two statements")


def audit_metrics_and_methods(frozen: Path, audit: Audit) -> None:
    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    expected = {
        "article": 77,
        "analysis_seed": 77001,
        "plot_seed": 20260777,
        "snapshot_date": "2026-08-23",
        "policy_sources": 18,
        "policy_assertions": 22,
        "routing_assets": 15,
        "object_records": 14,
        "identifier_records": 15,
        "existing_third_party_identifiers": 5,
        "pending_or_blocked_identifiers": 9,
        "external_ready_records": 0,
        "release_gates": 18,
        "passed_release_gates": 3,
        "raw_runs": 2,
        "raw_fastq_files": 4,
        "primary_catalog_genes": 93782,
        "article49_blocked_gates": 5,
        "database_records": 3,
        "source_artifacts": 15,
        "release_artifacts": 11,
        "new_accessions_or_dois_claimed": 0,
    }
    for key, value in expected.items():
        audit.add("Metric", key, metrics.get(key) == value, metrics.get(key))
    audit.add("Metric", "blocked-gates", metrics.get("blocked_release_gates", 0) >= 10, metrics.get("blocked_release_gates"))
    audit.add("Metric", "raw-bytes-positive", metrics.get("raw_bytes", 0) > 8_000_000_000, metrics.get("raw_bytes"))
    audit.add("Metric", "worktree-not-clean", metrics.get("code_worktree_clean") is False, metrics.get("code_worktree_clean"))

    methods = json.loads((frozen / "methods-contract.json").read_text(encoding="utf-8"))
    audit.add("Methods", "article", methods.get("article") == 77, methods.get("article"))
    audit.add("Methods", "no-external-action", methods.get("external_action", "").startswith("None"), methods.get("external_action"))
    audit.add("Methods", "pending-policy", "never a citeable identifier" in methods.get("identifier_policy", ""), methods.get("identifier_policy"))
    audit.add("Methods", "ownership", "not resubmitted" in methods.get("ownership_policy", ""), methods.get("ownership_policy"))
    audit.add("Methods", "database-four-fields", all(token in methods.get("database_policy", "") for token in ("Release", "artifact", "checksum", "source URL")), methods.get("database_policy"))
    audit.add("Methods", "separate-statements", "separate" in methods.get("availability_policy", ""), methods.get("availability_policy"))


def audit_figures(root: Path, frozen: Path, audit: Audit) -> None:
    figures = root / "figures/article77"
    manifest = json.loads((frozen / "figure-manifest.json").read_text(encoding="utf-8"))
    records = {row["file"]: row for row in manifest["figures"]}
    expected = {f"{stem}.{extension}" for stem in FIGURE_STEMS for extension in ("png", "svg")} | {ANCHOR}
    audit.add("Figures", "manifest-article", manifest.get("article") == 77, manifest.get("article"))
    audit.add("Figures", "seed", manifest.get("plot_seed") == 20260777, manifest.get("plot_seed"))
    audit.add("Figures", "files", set(records) == expected, sorted(records))
    audit.add("Figures", "count", len(records) == 19, len(records))
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
    with Image.open(figures / ANCHOR) as image:
        audit.add("Figure anchor", "dimensions", image.size == (767, 508), image.size)


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
        "## 这一步对应论文里的哪张图", "## 理论：为什么这么做", "## 准备工作",
        "## 可复制代码", "## 审计与升级", "## 出版级美化",
        "## 常见坑", "## 这段 Methods 怎么写", "## 换成你自己的数据怎么做", "## 参考",
    ):
        audit.add("Chapter heading", heading, heading in text, heading)
    audit.add("Chapter", "inline-theme", all(token in text for token in ("pal_pub <-", "scale_color_pub <-", "scale_fill_pub <-", "theme_pub <-", "save_pub <-")), "five helpers")
    audit.add("Chapter", "no-source-theme", 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text, "inline only")
    audit.add("Chapter", "deterministic", "set.seed(20260777)" in text, "set.seed(20260777)")
    audit.add("Chapter", "official-citations", all(token in text for token in ("@tenhoopen2017metagenomic", "@ncbi2026srasubmit", "@ena2026readsubmission", "@zenodo2026records", "@github2026zenodo", "@docker2026digests", "@nature2026dataavailability")), "official sources")
    for token in ("PRJEB52977", "SAMEA14435832", "SAMEA14435833", "ERR9765746", "ERR9765747", "PENDING_NOT_INVENTED", "93,782", "GTDB R232", "CheckM2", "ProGenomes 2.1"):
        audit.add("Chapter evidence", token, token in text, token)
    for token in ("SUB#", "ERZ", "-validate", "-submit", "CITATION.cff", "index digest", "platform digest", "Data Availability", "Code Availability"):
        audit.add("Chapter concept", token, token in text, token)
    audit.add("Chapter", "third-party-boundary", "不能重新提交" in text and "第三方" in text, "ownership wording")
    audit.add("Chapter", "not-a-publication", "不是一次公开发布" in text, "publication boundary")
    audit.add("Chapter", "no-fake-new-doi", "10.5281/zenodo.PENDING" not in text, "fake DOI scan")
    audit.add("Chapter", "no-local-shared-path", "/shared/" not in text, "path scan")
    audit.add("Chapter", "anchor-attribution", "CC BY 4.0" in text and "10.1093/gigascience/gix047" in text, "rights")
    for stem in FIGURE_STEMS:
        audit.add("Chapter figure", stem, f"../figures/article77/{stem}.png" in text, stem)
    audit.add("Chapter figure", ANCHOR, f"../figures/article77/{ANCHOR}" in text, ANCHOR)
    audit.add("Chapter", "no-placeholder", not any(term in text for term in ("Planned chapter", "TODO", "待补", "draft: true")), "placeholder scan")

    audit.add("Render", "exists", rendered.is_file(), rendered)
    if rendered.is_file():
        html = rendered.read_text(encoding="utf-8")
        audit.add("Render", "size", rendered.stat().st_size > 180_000, rendered.stat().st_size)
        audit.add("Render", "title", "数据、代码与基因组提交" in html, "title")
        audit.add("Render", "ten-images", sum(name in html for name in [ANCHOR, *[f"{stem}.png" for stem in FIGURE_STEMS]]) == 10, "image references")
        audit.add("Render", "no-tofu", "□□" not in html, "tofu scan")
        audit.add("Render", "availability", "Data Availability" in html and "Code Availability" in html, "availability sections")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_bundle(frozen, audit)
    audit_sources(root, frozen, audit)
    audit_policy_and_routing(frozen, audit)
    audit_identifiers_and_objects(frozen, audit)
    audit_release_packet(frozen, audit)
    audit_availability(frozen, audit)
    audit_metrics_and_methods(frozen, audit)
    audit_figures(root, frozen, audit)
    audit_chapter(args.chapter.resolve(), args.rendered_html.resolve(), audit)
    audit.write(args.qa_dir.resolve())


if __name__ == "__main__":
    main()
