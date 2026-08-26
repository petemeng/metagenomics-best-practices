#!/usr/bin/env python3
"""Fail-closed validation for Article 45 dRep dereplication evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from article42_44_validation_utils import (
    Audit,
    as_bool,
    audit_checksums,
    audit_figures,
    finish,
    read_tsv,
)


FIGURES = (
    "45-dereplication-yield",
    "45-representative-quality",
    "45-ani-alignment-audit",
    "45-source-retention",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def gunzip_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    script = (frozen / "scripts/plot_article45_drep.R").read_text(encoding="utf-8")
    expected = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260745" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "耗时")),
        "real-study": all(token in text for token in ("PRJEB52977", "ERR9765746", "ERR9765747")),
        "frozen-evidence": "data/small/45-drep-dereplication-frozen" in text,
        "methods-results": all(token in text for token in ("Methods template", "Results template")),
        "versions": all(token in text for token in ("dRep 3.6.2", "Mash 2.3", "fastANI 1.34")),
        "thresholds": all(token in text for token in ("-pa 0.90", "-sa 0.95", "-sa 0.999", "-nc 0.30")),
        "truth-blind": "不参与聚类、阈值或 representative selection" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", script) is None,
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-meta-prose": not any(token in text for token in (
            "本篇可独立", "本文可独立", "全系列约定", "接口只学一次",
            "作者代码通常长这样", "（即本文）", "无头服务器",
        )),
    }
    required_sections = (
        "这一步对应论文里的哪张图",
        "理论：",
        "准备工作",
        "可复制代码",
        "审计与升级",
        "出版级美化",
        "常见坑",
        "这段 Methods 怎么写",
        "换成你自己的数据怎么做",
        "参考",
    )
    for section in required_sections:
        expected[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in expected.items():
        audit.add("Chapter", check, status, check)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    args = parser.parse_args()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_checksums(frozen, audit)

    summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    contract = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    inputs = read_tsv(frozen / "input-genomes.tsv")
    membership = read_tsv(frozen / "cluster-membership.tsv.gz")
    representatives = read_tsv(frozen / "representative-genomes.tsv")
    thresholds = read_tsv(frozen / "threshold-summary.tsv")
    stability = read_tsv(frozen / "representative-stability.tsv")
    pairs = read_tsv(frozen / "pairwise-ani.tsv.gz")
    resources = read_tsv(frozen / "resource-summary.tsv")
    commands = read_tsv(frozen / "command-log.tsv")
    tools = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}

    audit.add("Identity", "article", summary.get("article") == 45 and contract.get("article") == 45, summary)
    audit.add("Identity", "truth-blind", summary.get("truth_used_for_clustering_or_selection") is False and contract.get("truth_used_for_clustering_or_selection") is False, summary)
    audit.add("Identity", "deterministic-contract", contract.get("random_process") is False and contract.get("input_order") == "lexicographic genome basename", contract)
    audit.add("Contract", "thresholds", contract.get("main_secondary_ani") == 0.95 and contract.get("sensitivity_secondary_ani") == 0.999 and contract.get("minimum_alignment_fraction") == 0.30, contract)
    audit.add("Input", "124-unique", len(inputs) == 124 and len({row["Genome"] for row in inputs}) == 124, len(inputs))
    audit.add("Input", "source-stages", Counter(row["SourceStage"] for row in inputs) == {"Article42-QC-pass": 101, "Article44-selected": 23}, Counter(row["SourceStage"] for row in inputs))
    audit.add("Input", "quality-gate", all(float(row["Completeness"]) >= 50 and float(row["Contamination"]) < 10 for row in inputs), len(inputs))
    audit.add("Input", "sha256-formatted", all(re.fullmatch(r"[0-9a-f]{64}", row["SHA256"]) for row in inputs), len(inputs))
    audit.add("Input", "exact-sequence-identities-retained", len({row["SHA256"] for row in inputs}) == 108, len({row["SHA256"] for row in inputs}))

    by_branch: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in membership:
        by_branch[row["Branch"]].append(row)
    expected_clusters = {"Species 95% ANI": 24, "Near-clone 99.9% ANI": 25}
    input_names = {row["Genome"] for row in inputs}
    for branch, cluster_count in expected_clusters.items():
        rows = by_branch[branch]
        observed_clusters = {row["Cluster"] for row in rows}
        observed_reps = {row["Genome"] for row in rows if as_bool(row["IsRepresentative"])}
        listed_reps = {row["Representative"] for row in rows}
        audit.add("Cluster", f"{branch}-coordinate", len(rows) == 124 and {row["Genome"] for row in rows} == input_names, len(rows))
        audit.add("Cluster", f"{branch}-counts", len(observed_clusters) == cluster_count and len(observed_reps) == cluster_count, {"clusters": len(observed_clusters), "representatives": len(observed_reps)})
        audit.add("Cluster", f"{branch}-representatives", observed_reps == listed_reps, len(observed_reps))
        audit.add("Cluster", f"{branch}-cluster-size", all(int(row["ClusterSize"]) == sum(other["Cluster"] == row["Cluster"] for other in rows) for row in rows), branch)

    threshold_map = {row["Branch"]: row for row in thresholds}
    audit.add("Result", "species-summary", summary == {
        "article": 45,
        "input_genomes": 124,
        "near_clone_clusters": 25,
        "species_clusters": 24,
        "species_clusters_split_at_99_9": 1,
        "species_representatives_from_article44": 16,
        "truth_used_for_clustering_or_selection": False,
    }, summary)
    audit.add("Result", "species-threshold-row", int(threshold_map["Species 95% ANI"]["Representatives"]) == 24 and int(threshold_map["Species 95% ANI"]["GenomesRemoved"]) == 100 and float(threshold_map["Species 95% ANI"]["MinimumAlignmentFractionPct"]) == 30, threshold_map["Species 95% ANI"])
    audit.add("Result", "nearclone-threshold-row", int(threshold_map["Near-clone 99.9% ANI"]["Representatives"]) == 25 and int(threshold_map["Near-clone 99.9% ANI"]["GenomesRemoved"]) == 99, threshold_map["Near-clone 99.9% ANI"])
    audit.add("Result", "one-species-split", len(stability) == 24 and sum(as_bool(row["SplitsAt99_9Pct"]) for row in stability) == 1, len(stability))
    near_pairs = [row for row in pairs if row["Branch"] == "Near-clone 99.9% ANI"]
    audit.add("Result", "nearclone-pair-audit", len(near_pairs) == 270 and sum(not as_bool(row["SameCluster"]) for row in near_pairs) == 4, len(near_pairs))
    audit.add("Result", "representative-ledger", len(representatives) == 49 and sum(row["Branch"] == "Species 95% ANI" for row in representatives) == 24, len(representatives))
    expected_rep_hashes = {
        row["Genome"]: row["SHA256"]
        for row in membership
        if row["Branch"] == "Species 95% ANI" and as_bool(row["IsRepresentative"])
    }
    frozen_reps = sorted((frozen / "representative-genomes").glob("*.fna.gz"))
    observed_rep_hashes = {path.name.removesuffix(".gz"): gunzip_sha256(path) for path in frozen_reps}
    audit.add("Result", "representative-fasta-coordinate", len(frozen_reps) == 24 and observed_rep_hashes == expected_rep_hashes, len(frozen_reps))

    for branch, native, cutoff in (("Species 95% ANI", "species95", 0.05), ("Near-clone 99.9% ANI", "nearclone999", 0.001)):
        cdb = read_csv(frozen / f"raw/{native}/Cdb.csv")
        wdb = read_csv(frozen / f"raw/{native}/Wdb.csv")
        ndb = read_csv(frozen / f"raw/{native}/Ndb.csv")
        genome_info = read_csv(frozen / f"raw/{native}/genomeInformation.csv")
        summary_map = {row["Genome"]: row["Cluster"] for row in by_branch[branch]}
        audit.add("Native", f"{native}-Cdb", len(cdb) == 124 and {row["genome"]: row["secondary_cluster"] for row in cdb} == summary_map, len(cdb))
        audit.add("Native", f"{native}-Cdb-threshold", all(abs(float(row["threshold"]) - cutoff) < 1e-9 and row["comparison_algorithm"] == "fastANI" for row in cdb), cutoff)
        audit.add("Native", f"{native}-Wdb", len(wdb) == expected_clusters[branch] and {row["genome"] for row in wdb} == {row["Genome"] for row in by_branch[branch] if as_bool(row["IsRepresentative"])}, len(wdb))
        audit.add("Native", f"{native}-Ndb", len(ndb) == 664 and all(0 <= float(row["ani"]) <= 1 and 0 <= float(row["alignment_coverage"]) <= 1 for row in ndb), len(ndb))
        audit.add("Native", f"{native}-genome-info", len(genome_info) == 124 and {row["genome"] for row in genome_info} == input_names, len(genome_info))

    expected_versions = {
        "dRep": "3.6.2",
        "Mash": "2.3",
        "fastANI": "version 1.34",
        "Prodigal": "Prodigal V2.6.3: February, 2016",
    }
    audit.add("Version", "tools", tools == expected_versions, tools)
    audit.add("Execution", "commands", len(commands) == 2 and all(int(row["ExitStatus"]) == 0 for row in commands) and any("-sa 0.95" in row["Command"] for row in commands) and any("-sa 0.999" in row["Command"] for row in commands), len(commands))
    audit.add("Execution", "resources", len(resources) == 2 and all(int(row["ExitStatus"]) == 0 and float(row["WallSeconds"]) > 0 and float(row["PeakRAMGiB"]) > 0 for row in resources), resources)
    log_names = {path.name for path in (frozen / "logs").iterdir() if path.is_file()}
    audit.add("Execution", "logs-current-only", log_names == {
        "drep-species95.stdout.log", "drep-species95.stderr.log", "drep-species95.time.txt",
        "drep-nearclone999.stdout.log", "drep-nearclone999.stderr.log", "drep-nearclone999.time.txt",
    }, sorted(log_names))

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=45,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={"input_genomes": 124, "species_clusters": 24, "near_clone_clusters": 25},
    )


if __name__ == "__main__":
    raise SystemExit(main())
