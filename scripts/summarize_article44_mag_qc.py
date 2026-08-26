#!/usr/bin/env python3
"""Build the complete MIMAG and assembly-graph audit tables for Article 44."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from article41_44_utils import dump_json, fasta_records, fasta_summary, parse_time, read_tsv, sha256, write_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--article43-frozen", type=Path, required=True)
    return parser.parse_args()


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "pass"}


def qc_tables(work: Path):
    checkm2, gunc, checkm1 = {}, {}, {}
    for row in csv.DictReader((work / "qc/checkm2/quality_report.tsv").open(encoding="utf-8"), delimiter="\t"):
        name = row.get("Name") or row.get("name") or row.get("Genome")
        checkm2[name.removesuffix(".fna")] = row
    maxcss = sorted((work / "qc/gunc").glob("GUNC.*.maxCSS_level.tsv"))
    if len(maxcss) != 1:
        raise RuntimeError(f"Expected one GUNC maxCSS table: {maxcss}")
    for row in csv.DictReader(maxcss[0].open(encoding="utf-8"), delimiter="\t"):
        name = row.get("genome") or row.get("Genome")
        gunc[name.removesuffix(".fna")] = row
    for row in csv.DictReader((work / "qc/checkm1-qa.tsv").open(encoding="utf-8"), delimiter="\t"):
        name = row.get("Bin Id") or row.get("Bin ID") or row.get("Name")
        if not name:
            raise RuntimeError(f"Cannot identify CheckM1 bin ID column: {list(row)}")
        checkm1[name.removesuffix(".fna")] = row
    return checkm2, gunc, checkm1


def parse_barrnap(path: Path, mag: str) -> list[dict[str, object]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) != 9:
                raise RuntimeError(f"Unexpected barrnap GFF row: {raw[:200]}")
            attributes = parts[8]
            marker = next((name for name in ("5S", "16S", "23S") if name in attributes), "Other")
            incomplete = "partial" in attributes.lower() or "aligned only" in attributes.lower()
            rows.append({
                "MAG": mag, "Contig": parts[0], "Start": int(parts[3]), "End": int(parts[4]),
                "Strand": parts[6], "Marker": marker, "CompleteHit": not incomplete,
                "Attributes": attributes,
            })
    return rows


def parse_trnascan(path: Path, mag: str) -> list[dict[str, object]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("Sequence") or raw.startswith("Name") or set(raw.strip()) == {"-"}:
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9:
                parts = raw.split()
            if len(parts) < 9 or not parts[1].isdigit():
                continue
            note = " ".join(parts[9:]) if len(parts) > 9 else ""
            rows.append({
                "MAG": mag, "Contig": parts[0], "TRNANumber": int(parts[1]),
                "Start": int(parts[2]), "End": int(parts[3]), "Isotype": parts[4],
                "Anticodon": parts[5], "Score": float(parts[8]), "Note": note,
                "CountForMIMAG": parts[4] not in {"Pseudo", "Undet", "Sup"} and "pseudo" not in note.lower(),
            })
    return rows


def coding_summary(path: Path) -> tuple[int, int]:
    genes = 0
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) == 9 and parts[2] == "CDS":
                genes += 1
                intervals[parts[0]].append((int(parts[3]), int(parts[4])))
    coding_bp = 0
    for contig_intervals in intervals.values():
        ordered = sorted((min(left, right), max(left, right)) for left, right in contig_intervals)
        start, end = ordered[0]
        for left, right in ordered[1:]:
            if left > end + 1:
                coding_bp += end - start + 1
                start, end = left, right
            else:
                end = max(end, right)
        coding_bp += end - start + 1
    return genes, coding_bp


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def map_fastg_nodes(path: Path, source_path: Path) -> tuple[dict[str, str], list[dict[str, object]]]:
    """Map toolkit-renamed FASTG nodes back to k141 contigs by exact sequence hash."""
    source_by_hash: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    for order, (contig, sequence) in enumerate(fasta_records(source_path), start=1):
        source_by_hash[sequence_sha256(sequence)].append((contig, order, len(sequence)))
    node_map: dict[str, str] = {}
    audit_rows: list[dict[str, object]] = []
    pattern = re.compile(r"NODE_(\d+)_")
    for header, sequence in fasta_records(path):
        source_token = header.split(":", 1)[0].rstrip(";")
        if "'" in source_token:
            continue
        match = pattern.search(source_token)
        if not match:
            raise RuntimeError(f"Cannot parse FASTG node header: {header[:200]}")
        node = f"NODE_{int(match.group(1))}"
        digest = sequence_sha256(sequence)
        candidates = source_by_hash.get(digest, [])
        if len(candidates) != 1:
            raise RuntimeError(
                f"FASTG node {node} has {len(candidates)} exact k141 sequence matches"
            )
        contig, source_order, length_bp = candidates[0]
        previous = node_map.setdefault(node, contig)
        if previous != contig:
            raise RuntimeError(f"FASTG node {node} maps to multiple k141 contigs")
        audit_rows.append({
            "FastgNode": node, "Contig": contig, "SourceOrder": source_order,
            "LengthBp": length_bp, "SequenceSHA256": digest,
            "MappingRule": "Exact forward-sequence SHA-256",
        })
    source_count = sum(len(values) for values in source_by_hash.values())
    if len(node_map) != source_count or len({row["Contig"] for row in audit_rows}) != source_count:
        raise RuntimeError(
            f"FASTG/k141 coordinate mismatch: nodes={len(node_map)}, source={source_count}"
        )
    return node_map, sorted(audit_rows, key=lambda row: int(str(row["FastgNode"]).split("_")[1]))


def parse_fastg_edges(path: Path, node_map: dict[str, str]) -> tuple[set[str], set[tuple[str, str]]]:
    nodes: set[str] = set()
    edges: set[tuple[str, str]] = set()
    pattern = re.compile(r"NODE_(\d+)_")
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.startswith(">"):
                continue
            node_numbers = pattern.findall(raw)
            if not node_numbers:
                raise RuntimeError(f"Cannot parse FASTG adjacency header: {raw[:200]}")
            fastg_nodes = [f"NODE_{int(number)}" for number in node_numbers]
            missing = [node for node in fastg_nodes if node not in node_map]
            if missing:
                raise RuntimeError(f"FASTG adjacency references unmapped nodes: {missing[:5]}")
            names = [node_map[node] for node in fastg_nodes]
            source = names[0]
            nodes.add(source)
            for target in names[1:]:
                nodes.add(target)
                edges.add(tuple(sorted((source, target))))
    return nodes, edges


def components(nodes: set[str], edges: set[tuple[str, str]]) -> int:
    if not nodes:
        return 0
    adjacency: dict[str, set[str]] = {node: set() for node in nodes}
    for left, right in edges:
        if left in nodes and right in nodes and left != right:
            adjacency[left].add(right)
            adjacency[right].add(left)
    seen = set()
    count = 0
    for start in sorted(nodes):
        if start in seen:
            continue
        count += 1
        stack = [start]
        seen.add(start)
        while stack:
            node = stack.pop()
            for neighbour in adjacency[node]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
    return count


def first_float(row: dict[str, str], *keys: str) -> float:
    for key in keys:
        value = row.get(key)
        if value not in (None, "", "NA", "N/A"):
            return float(value)
    return math.nan


def main() -> int:
    args = parse_args()
    root, work, frozen43 = args.project_root.resolve(), args.work_dir.resolve(), args.article43_frozen.resolve()
    if not (work / ".article44-run-complete").is_file():
        raise FileNotFoundError("Run run_article44_mag_qc.py first")
    summary = work / "summary"
    summary.mkdir(exist_ok=True)
    checkm2, gunc, checkm1 = qc_tables(work)
    mags = sorted(path.stem for path in (work / "bins").glob("*.fna"))
    if set(mags) != set(checkm2) or set(mags) != set(gunc) or set(mags) != set(checkm1):
        raise RuntimeError("Article 44 MAG/QC coordinate mismatch")
    domain_rows = read_tsv(work / "qc-domain-audit.tsv")
    domain_by_mag = {row["MAG"]: row["Domain"] for row in domain_rows}
    if set(mags) != set(domain_by_mag) or any(domain not in {"Archaea", "Bacteria"} for domain in domain_by_mag.values()):
        raise RuntimeError("Article 44 marker-domain audit mismatch")

    rrna_rows, trna_rows = [], []
    basic = {}
    for mag in mags:
        fasta = work / "bins" / f"{mag}.fna"
        stats, _ = fasta_summary(fasta)
        genes, coding_bp = coding_summary(work / "features" / mag / "prodigal.gff")
        rrna_rows.extend(parse_barrnap(work / "features" / mag / "barrnap.gff", mag))
        trna_rows.extend(parse_trnascan(work / "features" / mag / "trnascan.tsv", mag))
        basic[mag] = {
            "Contigs": stats["Contigs"], "BinBp": stats["TotalBp"], "N50Bp": stats["N50Bp"],
            "GCPct": stats["GCPct"], "ProteinCodingGenes": genes,
            "CodingDensityPct": 100 * coding_bp / int(stats["TotalBp"]), "MAGFASTA_SHA256": sha256(fasta),
        }
    write_tsv(summary / "rrna-features.tsv", rrna_rows, fieldnames=["MAG", "Contig", "Start", "End", "Strand", "Marker", "CompleteHit", "Attributes"])
    write_tsv(summary / "trna-features.tsv", trna_rows, fieldnames=["MAG", "Contig", "TRNANumber", "Start", "End", "Isotype", "Anticodon", "Score", "Note", "CountForMIMAG"])

    k141_source = root / "data/raw/article30/work/assemblies/megahit-coassembly/intermediate_contigs/k141.contigs.fa"
    fastg_node_map, node_map_rows = map_fastg_nodes(work / "graph/megahit-k141.fastg", k141_source)
    write_tsv(summary / "assembly-graph-node-map.tsv.gz", node_map_rows)
    graph_nodes, graph_edges = parse_fastg_edges(work / "graph/megahit-k141.fastg", fastg_node_map)
    membership = read_tsv(frozen43 / "selected-refinement-membership.tsv.gz")
    mag_nodes: dict[str, set[str]] = defaultdict(set)
    for row in membership:
        mag_nodes[row["RefinedID"]].add(row["Contig"])
    full_source = root / "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz"
    full_order = [name for name, _ in fasta_records(full_source)]
    pair_rows = read_tsv(root / "data/small/41-read-mapping-depth-frozen/raw/paired-contigs.tsv")
    pair_edges = []
    for row in pair_rows:
        left_i, right_i = int(row["contigIdx"]), int(row["contigIdxMate"])
        if left_i >= len(full_order) or right_i >= len(full_order):
            raise RuntimeError("Article 41 paired-contig index exceeds Article 30 coordinate set")
        pair_edges.append((full_order[left_i], full_order[right_i], float(row["AvgCoverage"])))
    graph_audit = []
    graph_link_rows = []
    for mag in mags:
        members = mag_nodes[mag]
        present = members & graph_nodes
        internal = {edge for edge in graph_edges if edge[0] in members and edge[1] in members}
        boundary = {edge for edge in graph_edges if (edge[0] in members) ^ (edge[1] in members)}
        self_loops = {edge for edge in internal if edge[0] == edge[1]}
        nonself_internal = {edge for edge in internal if edge[0] != edge[1]}
        pair_internal = sum(weight for left, right, weight in pair_edges if left != right and left in members and right in members)
        pair_boundary = sum(weight for left, right, weight in pair_edges if (left in members) ^ (right in members))
        pair_fraction = 100 * pair_boundary / (pair_internal + pair_boundary) if pair_internal + pair_boundary else math.nan
        component_count = components(members, nonself_internal)
        if len(members) == 1:
            continuity = "Single contig"
        elif component_count == 1:
            continuity = "Connected in k141 graph"
        else:
            continuity = "Fragmented in k141 graph"
        for left, right in sorted(internal):
            graph_link_rows.append({"MAG": mag, "Evidence": "MEGAHIT k141", "LinkClass": "Internal", "Node1": left, "Node2": right, "Weight": 1})
        for left, right in sorted(boundary):
            graph_link_rows.append({"MAG": mag, "Evidence": "MEGAHIT k141", "LinkClass": "Boundary", "Node1": left, "Node2": right, "Weight": 1})
        graph_audit.append({
            "MAG": mag, "MAGContigs": len(members), "K141NodesPresent": len(present),
            "K141InternalEdges": len(nonself_internal), "K141BoundaryEdges": len(boundary),
            "K141SelfLoops": len(self_loops), "K141Components": component_count,
            "PairedInternalWeight": pair_internal, "PairedBoundaryWeight": pair_boundary,
            "PairedBoundaryPct": pair_fraction, "GraphContinuity": continuity,
            "CircularGenomeClaim": "No—graph evidence alone is insufficient",
        })
    write_tsv(summary / "assembly-graph-audit.tsv", graph_audit)
    write_tsv(summary / "assembly-graph-links.tsv.gz", graph_link_rows, fieldnames=["MAG", "Evidence", "LinkClass", "Node1", "Node2", "Weight"])

    graph_by_mag = {row["MAG"]: row for row in graph_audit}
    quality_rows = []
    for mag in mags:
        cm2, gu, cm1 = checkm2[mag], gunc[mag], checkm1[mag]
        completeness = float(cm2["Completeness"])
        contamination = float(cm2["Contamination"])
        gunc_pass = as_bool(gu.get("pass.GUNC", gu.get("pass_gunc", "")))
        mag_rrna = [row for row in rrna_rows if row["MAG"] == mag and row["CompleteHit"]]
        rrna_markers = {row["Marker"] for row in mag_rrna}
        mag_trna = [row for row in trna_rows if row["MAG"] == mag and row["CountForMIMAG"]]
        trna_isotypes = {row["Isotype"] for row in mag_trna}
        pre_marker_hq = completeness > 90 and contamination < 5 and gunc_pass
        if pre_marker_hq and {"5S", "16S", "23S"} <= rrna_markers and len(trna_isotypes) >= 18:
            tier = "High quality"
        elif completeness >= 50 and contamination < 10 and gunc_pass:
            tier = "Medium quality"
        else:
            tier = "Low/failed"
        lineage = cm1.get("Marker lineage") or cm1.get("Marker Lineage") or "Unknown"
        strain = first_float(cm1, "Strain heterogeneity", "Strain Heterogeneity")
        css = first_float(gu, "clade_separation_score", "CSS")
        quality_rows.append({
            "MAG": mag, **basic[mag], "CheckM2Completeness": completeness,
            "CheckM2Contamination": contamination, "CheckM2Model": cm2.get("Completeness_Model_Used", ""),
            "GUNCPass": gunc_pass, "GUNCCladeSeparationScore": css,
            "CheckM1MarkerLineage": lineage, "MarkerDomain": domain_by_mag[mag],
            "CheckM1Completeness": first_float(cm1, "Completeness"),
            "CheckM1Contamination": first_float(cm1, "Contamination"), "CheckM1StrainHeterogeneityPct": strain,
            "Complete5S": "5S" in rrna_markers, "Complete16S": "16S" in rrna_markers,
            "Complete23S": "23S" in rrna_markers, "CompleteRRNASet": {"5S", "16S", "23S"} <= rrna_markers,
            "TRNAGenes": len(mag_trna), "TRNAIsotypes": len(trna_isotypes),
            "PreMarkerHQ": pre_marker_hq, "MIMAGQuality": tier,
            "GraphContinuity": graph_by_mag[mag]["GraphContinuity"],
            "K141Components": graph_by_mag[mag]["K141Components"],
            "PairedBoundaryPct": graph_by_mag[mag]["PairedBoundaryPct"],
            "Taxonomy": "Pending independent GTDB-Tk classification",
        })
    write_tsv(summary / "mag-quality-summary.tsv", quality_rows)

    truth = {row["RefinedID"]: row for row in read_tsv(frozen43 / "selected-mag-candidates.tsv")}
    truth_rows = []
    for row in quality_rows:
        audit = truth[row["MAG"]]
        truth_rows.append({
            "MAG": row["MAG"], "MIMAGQuality": row["MIMAGQuality"],
            "DominantTruthAccession": audit["DominantTruthAccession"], "DominantTruthReference": audit["DominantTruthReference"],
            "DominantGenomeRecoveryPct": audit["DominantGenomeRecoveryPct"], "AlignedPurityPct": audit["AlignedPurityPct"],
            "TruthProxyTier": audit["TruthProxyTier"], "TruthUsedForMIMAGTier": "No",
        })
    write_tsv(summary / "mag-quality-truth-audit.tsv", truth_rows)
    counts = Counter(row["MIMAGQuality"] for row in quality_rows)
    write_tsv(summary / "mimag-tier-counts.tsv", [
        {"MIMAGQuality": tier, "MAGs": counts[tier]} for tier in ("High quality", "Medium quality", "Low/failed")
    ])
    resource_rows = [parse_time(path) for path in sorted((work / "logs").glob("*.time.txt"))]
    if any(int(row["ExitStatus"]) != 0 for row in resource_rows):
        raise RuntimeError("An Article 44 command has a non-zero resource exit status")
    write_tsv(summary / "resource-summary.tsv", resource_rows)
    payload = {
        "article": 44, "seed": 20260744, "selected_bins": len(mags),
        "mimag_counts": dict(counts), "pre_marker_hq": sum(bool(row["PreMarkerHQ"]) for row in quality_rows),
        "complete_rrna_sets": sum(bool(row["CompleteRRNASet"]) for row in quality_rows),
        "minimum_18_trna_isotypes": sum(int(row["TRNAIsotypes"]) >= 18 for row in quality_rows),
        "gunc_pass": sum(bool(row["GUNCPass"]) for row in quality_rows),
        "marker_domains": dict(Counter(row["MarkerDomain"] for row in quality_rows)),
        "fastg_nodes_exactly_mapped": len(fastg_node_map),
        "fastg_undirected_edges": len(graph_edges),
        "mags_with_all_k141_nodes_present": sum(int(row["K141NodesPresent"]) == int(row["MAGContigs"]) for row in graph_audit),
        "mags_with_k141_boundary_edges": sum(int(row["K141BoundaryEdges"]) > 0 for row in graph_audit),
        "graph_boundaries_reported_not_used_for_mimag": True,
        "circular_genomes_claimed": 0,
        "truth_used_for_quality_classification": False,
    }
    dump_json(summary / "run-summary.json", payload)
    (work / ".article44-summary-complete").write_text("PASS\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
