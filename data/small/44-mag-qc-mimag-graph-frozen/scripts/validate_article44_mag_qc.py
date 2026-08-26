#!/usr/bin/env python3
"""Fail-closed validation for Article 44 complete MIMAG and graph evidence."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from article42_44_validation_utils import Audit, as_bool, audit_chapter, audit_checksums, audit_figures, finish, read_tsv


FIGURES = ("44-quality-landscape", "44-mimag-requirements", "44-assembly-graph-audit", "44-checkm-audit")


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
    quality = read_tsv(frozen / "mag-quality-summary.tsv")
    graph = read_tsv(frozen / "assembly-graph-audit.tsv")
    graph_node_map = read_tsv(frozen / "assembly-graph-node-map.tsv.gz")
    truth = read_tsv(frozen / "mag-quality-truth-audit.tsv")
    rrna = read_tsv(frozen / "rrna-features.tsv")
    trna = read_tsv(frozen / "trna-features.tsv")
    counts = read_tsv(frozen / "mimag-tier-counts.tsv")
    tools = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    raw_checkm2 = read_tsv(frozen / "raw/checkm2-quality-report.tsv")
    raw_gunc = read_tsv(frozen / "raw/gunc-maxcss.tsv")
    raw_checkm1 = read_tsv(frozen / "raw/checkm1-qa.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    commands = read_tsv(frozen / "command-log.tsv")
    domain_audit = read_tsv(frozen / "qc-domain-audit.tsv")
    reconstruction = read_tsv(frozen / "selected-mag-reconstruction-audit.tsv")

    audit.add("Identity", "article", summary.get("article") == 44 and contract.get("article") == 44, summary.get("article"))
    audit.add("Identity", "seed", contract.get("seed") == 20260744 and summary.get("seed") == 20260744, contract.get("seed"))
    audit.add("Identity", "truth-blinded", summary.get("truth_used_for_quality_classification") is False, summary)
    expected_tools = {"CheckM2": "1.1.0", "GUNC": "1.1.0", "CheckM1": "1.2.5", "barrnap": "1.10.5", "tRNAscan-SE": "2.0.13", "Prodigal": "2.6.3", "MEGAHIT toolkit": "1.2.9"}
    audit.add("Version", "tools", tools == expected_tools, tools)
    audit.add("Input", "reconstruction-pass", len(reconstruction) == len(quality) and all(row["Status"] == "PASS" for row in reconstruction), len(reconstruction))
    audit.add("QC", "raw-coordinate", len(raw_checkm2) == len(raw_gunc) == len(raw_checkm1) == len(quality), len(quality))
    audit.add("QC", "graph-coordinate", {row["MAG"] for row in graph} == {row["MAG"] for row in quality}, len(graph))
    audit.add("QC", "truth-coordinate", {row["MAG"] for row in truth} == {row["MAG"] for row in quality}, len(truth))
    domains = {row["MAG"]: row["Domain"] for row in domain_audit}
    audit.add("Domain", "coordinate", set(domains) == {row["MAG"] for row in quality}, len(domains))
    audit.add("Domain", "values", set(domains.values()) == {"Archaea", "Bacteria"}, Counter(domains.values()))
    command_by_label = {row["Label"]: row["Command"] for row in commands}
    for row in quality:
        high = float(row["CheckM2Completeness"]) > 90 and float(row["CheckM2Contamination"]) < 5 and as_bool(row["GUNCPass"]) and as_bool(row["Complete5S"]) and as_bool(row["Complete16S"]) and as_bool(row["Complete23S"]) and int(row["TRNAIsotypes"]) >= 18
        medium = float(row["CheckM2Completeness"]) >= 50 and float(row["CheckM2Contamination"]) < 10 and as_bool(row["GUNCPass"])
        expected = "High quality" if high else "Medium quality" if medium else "Low/failed"
        audit.add("MIMAGFormula", row["MAG"], row["MIMAGQuality"] == expected, {"observed": row["MIMAGQuality"], "expected": expected})
        audit.add("MIMAGEvidence", f"{row['MAG']}-coding", int(row["ProteinCodingGenes"]) > 0 and 0 < float(row["CodingDensityPct"]) <= 100, row["CodingDensityPct"])
        audit.add("MIMAGEvidence", f"{row['MAG']}-taxonomy-pending", row["Taxonomy"] == "Pending independent GTDB-Tk classification", row["Taxonomy"])
        audit.add("Domain", f"{row['MAG']}-summary", row["MarkerDomain"] == domains.get(row["MAG"]), row["MarkerDomain"])
        kingdom = "arc" if row["MarkerDomain"] == "Archaea" else "bac"
        trna_mode = "-A" if row["MarkerDomain"] == "Archaea" else "-B"
        barrnap_command = command_by_label.get(f"barrnap-{row['MAG']}", "")
        trna_command = command_by_label.get(f"trnascan-{row['MAG']}", "")
        audit.add("DomainMode", f"{row['MAG']}-barrnap", f"--kingdom {kingdom}" in barrnap_command, barrnap_command)
        audit.add("DomainMode", f"{row['MAG']}-trnascan", f" {trna_mode} " in f" {trna_command} ", trna_command)
    observed_counts = Counter(row["MIMAGQuality"] for row in quality)
    audit.add("MIMAG", "tier-counts", {row["MIMAGQuality"]: int(row["MAGs"]) for row in counts} == {tier: observed_counts[tier] for tier in ("High quality", "Medium quality", "Low/failed")}, observed_counts)
    audit.add("MIMAG", "summary-counts", summary["selected_bins"] == len(quality) and summary["mimag_counts"] == dict(observed_counts), summary)
    audit.add("MIMAG", "rrna-features-real", all(row["Marker"] in {"5S", "16S", "23S", "Other"} for row in rrna), len(rrna))
    audit.add("MIMAG", "trna-evidence", all(float(row["Score"]) >= 0 for row in trna), len(trna))
    audit.add("Graph", "no-circular-claim", all(row["CircularGenomeClaim"] == "No—graph evidence alone is insufficient" for row in graph) and summary["circular_genomes_claimed"] == 0, len(graph))
    audit.add("Graph", "components-positive", all(int(row["K141Components"]) >= 1 for row in graph), len(graph))
    audit.add("Graph", "nodes-bounded", all(0 <= int(row["K141NodesPresent"]) <= int(row["MAGContigs"]) for row in graph), len(graph))
    audit.add(
        "Graph", "node-map-exact",
        len(graph_node_map) == summary["fastg_nodes_exactly_mapped"]
        and len({row["FastgNode"] for row in graph_node_map}) == len(graph_node_map)
        and len({row["Contig"] for row in graph_node_map}) == len(graph_node_map)
        and all(row["MappingRule"] == "Exact forward-sequence SHA-256" for row in graph_node_map),
        len(graph_node_map),
    )
    audit.add(
        "Graph", "selected-nodes-present",
        summary["mags_with_all_k141_nodes_present"] == len(graph)
        and all(int(row["K141NodesPresent"]) == int(row["MAGContigs"]) for row in graph),
        summary["mags_with_all_k141_nodes_present"],
    )
    audit.add(
        "Graph", "real-adjacency-evidence",
        summary["fastg_undirected_edges"] > 0
        and summary["mags_with_k141_boundary_edges"] == len(graph)
        and all(int(row["K141BoundaryEdges"]) > 0 for row in graph),
        {"edges": summary["fastg_undirected_edges"], "mags_with_boundary": summary["mags_with_k141_boundary_edges"]},
    )
    audit.add("Graph", "boundary-not-tier", summary["graph_boundaries_reported_not_used_for_mimag"] is True and contract["graph_role"].startswith("triage evidence"), summary)
    audit.add("Truth", "posthoc-only", all(row["TruthUsedForMIMAGTier"] == "No" for row in truth), len(truth))
    audit.add("Execution", "commands-success", len(commands) >= 7 and all(int(row["ReturnCode"]) == 0 for row in commands), len(commands))
    audit.add("Execution", "resources-success", len(resources) == len(commands) and all(int(row["ExitStatus"]) == 0 for row in resources), len(resources))
    audit_chapter(args.chapter.resolve(), audit, article=44, figure_stems=FIGURES, tokens=("CheckM2 1.1.0", "GUNC 1.1.0", "CheckM1 1.2.5", "barrnap 1.10.5", "tRNAscan-SE 2.0.13", "MIMAG", "assembly graph"))
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(article=44, audit=audit, output=args.output_dir.resolve(), payload={"selected_bins": len(quality), "mimag_counts": dict(observed_counts)})


if __name__ == "__main__":
    raise SystemExit(main())
