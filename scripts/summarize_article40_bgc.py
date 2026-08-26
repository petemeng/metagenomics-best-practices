#!/usr/bin/env python3
"""Summarize Article 40 antiSMASH/GECCO runs into compact audit tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from article37_40_utils import dump_json, parse_time, read_abundance, read_tsv, sha256, write_tsv


SEED = 20260740
DATASETS = ("salinispora-full", "salinispora-fragmented", "nostoc", "coassembly-ge20kb")
CATEGORY_ORDER = ("PKS", "NRPS", "RiPP", "Terpene", "Saccharide", "Other")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def normalize_category(value: str) -> str:
    label = value.strip().lower()
    if "pks" in label or "polyketide" in label:
        return "PKS"
    if "nrp" in label:
        return "NRPS"
    if "ripp" in label or "ribosomal" in label:
        return "RiPP"
    if "terpene" in label:
        return "Terpene"
    if "saccharide" in label:
        return "Saccharide"
    return "Other"


def ordered_categories(values: list[str]) -> list[str]:
    observed = {normalize_category(value) for value in values}
    return [category for category in CATEGORY_ORDER if category in observed] or ["Other"]


def best_mibig(record: dict[str, object], region_number: int) -> tuple[str, str, float | None]:
    modules = record.get("modules", {})
    compare = modules.get("antismash.modules.cluster_compare", {}) if isinstance(modules, dict) else {}
    databases = compare.get("db_results", {}) if isinstance(compare, dict) else {}
    mibig = databases.get("MIBiG", {}) if isinstance(databases, dict) else {}
    by_region = mibig.get("by_region", {}) if isinstance(mibig, dict) else {}
    region = by_region.get(str(region_number), {}) if isinstance(by_region, dict) else {}
    method = region.get("RegionToRegion_RiQ", {}) if isinstance(region, dict) else {}
    scores = method.get("scores_by_region", {}) if isinstance(method, dict) else {}
    references = method.get("reference_regions", {}) if isinstance(method, dict) else {}
    if not scores:
        return "", "", None
    reference, score = max(scores.items(), key=lambda item: float(item[1]))
    accession = reference.split(":", 1)[0]
    description = ""
    if isinstance(references, dict) and isinstance(references.get(reference), dict):
        description = str(references[reference].get("description", ""))
    return accession, description, float(score)


def best_knowncluster(record: dict[str, object], region_number: int) -> tuple[str, str, int | None]:
    modules = record.get("modules", {})
    clusterblast = modules.get("antismash.modules.clusterblast", {}) if isinstance(modules, dict) else {}
    known = clusterblast.get("knowncluster", {}) if isinstance(clusterblast, dict) else {}
    results = known.get("results", []) if isinstance(known, dict) else []
    hit = next((row for row in results if int(row.get("region_number", -1)) == region_number), None)
    if not hit or not hit.get("ranking"):
        return "", "", None
    reference, score = hit["ranking"][0]
    return str(reference.get("accession", "")), str(reference.get("description", "")), int(score.get("similarity", 0))


def parse_location(location: str) -> tuple[int, int]:
    numbers = [int(value) for value in re.findall(r"\d+", location)]
    if len(numbers) < 2:
        raise ValueError(f"Cannot parse antiSMASH feature location: {location}")
    return numbers[0], numbers[-1]


def parse_antismash(work: Path, dataset: str) -> tuple[list[dict[str, object]], dict[str, list[dict[str, object]]]]:
    path = work / "antismash" / dataset / f"{dataset}.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    lengths = fasta_lengths(work / "inputs" / f"{dataset}.fna")
    calls: list[dict[str, object]] = []
    genes: dict[str, list[dict[str, object]]] = {}
    for record in payload["records"]:
        sequence_id = record["id"]
        length = lengths[sequence_id]
        cds = []
        for feature in record.get("features", []):
            if feature.get("type") != "CDS":
                continue
            start, end = parse_location(feature["location"])
            cds.append({"Start0": start, "End0": end})
        cds.sort(key=lambda row: (int(row["Start0"]), int(row["End0"])))
        for ordinal, gene in enumerate(cds, 1):
            gene["GeneID"] = f"{sequence_id}_{ordinal}"
        genes[sequence_id] = cds
        for number, area in enumerate(record.get("areas", []), 1):
            categories = ordered_categories([str(item.get("category", "other")) for item in area.get("protoclusters", {}).values()])
            accession, description, riq = best_mibig(record, number)
            kc_accession, kc_description, kc_similarity = best_knowncluster(record, number)
            start, end = int(area["start"]), int(area["end"])
            calls.append({
                "Dataset": dataset,
                "Tool": "antiSMASH",
                "RegionID": f"{dataset}::{sequence_id}::r{number}",
                "SequenceID": sequence_id,
                "Start0": start,
                "End0": end,
                "LengthBp": end - start,
                "Products": ";".join(area.get("products", [])) or "Unknown",
                "Categories": ";".join(categories),
                "AtSequenceEdge": "yes" if start <= 0 or end >= length else "no",
                "MeanProbability": "",
                "MaximumProbability": "",
                "TopMIBiGAccession": accession,
                "TopMIBiGDescription": description,
                "TopMIBiGRiQ": "" if riq is None else f"{riq:.8f}",
                "TopKnownClusterAccession": kc_accession,
                "TopKnownClusterDescription": kc_description,
                "TopKnownClusterSimilarityPercent": "" if kc_similarity is None else kc_similarity,
            })
    return calls, genes


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name = ""
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0]
                lengths[name] = 0
            elif name:
                lengths[name] += len(line.strip())
    return lengths


def parse_gecco(work: Path, dataset: str) -> list[dict[str, object]]:
    path = work / "gecco" / dataset / f"{dataset}.clusters.tsv"
    lengths = fasta_lengths(work / "inputs" / f"{dataset}.fna")
    calls = []
    for row in read_tsv(path):
        sequence_id = row["sequence_id"]
        start, end = int(row["start"]), int(row["end"])
        calls.append({
            "Dataset": dataset,
            "Tool": "GECCO",
            "RegionID": f"{dataset}::{sequence_id}::{row['cluster_id']}",
            "SequenceID": sequence_id,
            "Start0": start,
            "End0": end,
            "LengthBp": end - start,
            "Products": row["type"] or "Unknown",
            "Categories": ";".join(ordered_categories((row["type"] or "Unknown").split(";"))),
            "AtSequenceEdge": "yes" if start <= 1 or end >= lengths[sequence_id] - 1 else "no",
            "MeanProbability": row["average_p"],
            "MaximumProbability": row["max_p"],
            "TopMIBiGAccession": "",
            "TopMIBiGDescription": "",
            "TopMIBiGRiQ": "",
            "TopKnownClusterAccession": "",
            "TopKnownClusterDescription": "",
            "TopKnownClusterSimilarityPercent": "",
        })
    return calls


def overlap_bp(left: dict[str, object], right: dict[str, object]) -> int:
    if left["SequenceID"] != right["SequenceID"]:
        return 0
    return max(0, min(int(left["End0"]), int(right["End0"])) - max(int(left["Start0"]), int(right["Start0"])))


def reciprocal_match(left: dict[str, object], right: dict[str, object], threshold: float = 0.25) -> bool:
    overlap = overlap_bp(left, right)
    return overlap / max(1, int(left["LengthBp"])) >= threshold and overlap / max(1, int(right["LengthBp"])) >= threshold


def tool_yield(calls: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for dataset in DATASETS:
        anti = [row for row in calls if row["Dataset"] == dataset and row["Tool"] == "antiSMASH"]
        gecco = [row for row in calls if row["Dataset"] == dataset and row["Tool"] == "GECCO"]
        anti_supported = sum(any(reciprocal_match(row, other) for other in gecco) for row in anti)
        gecco_supported = sum(any(reciprocal_match(row, other) for other in anti) for row in gecco)
        edges_anti = sum(row["AtSequenceEdge"] == "yes" for row in anti)
        edges_gecco = sum(row["AtSequenceEdge"] == "yes" for row in gecco)
        rows.append({
            "Dataset": dataset,
            "antiSMASHRegions": len(anti),
            "GECCORegions": len(gecco),
            "antiSMASHSupportedByGECCO": anti_supported,
            "GECCOSupportedByAntiSMASH": gecco_supported,
            "antiSMASHEdgeRegions": edges_anti,
            "GECCOEdgeRegions": edges_gecco,
            "ReciprocalOverlapThreshold": 0.25,
        })
    return rows


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    start, end = intervals[0]
    for next_start, next_end in intervals[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end
    return total + end - start


def fragmentation_summary(calls: list[dict[str, object]], fragment_map: list[dict[str, str]]) -> list[dict[str, object]]:
    offsets = {row["FragmentID"]: (int(row["GenomeStart1"]) - 1, int(row["LengthBp"])) for row in fragment_map}
    rows = []
    for tool in ("antiSMASH", "GECCO"):
        complete = [row for row in calls if row["Dataset"] == "salinispora-full" and row["Tool"] == tool]
        fragmented = [row for row in calls if row["Dataset"] == "salinispora-fragmented" and row["Tool"] == tool]
        mapped = []
        for row in fragmented:
            offset, fragment_length = offsets[str(row["SequenceID"])]
            mapped.append({**row, "GenomeStart0": offset + int(row["Start0"]), "GenomeEnd0": offset + int(row["End0"]), "FragmentLength": fragment_length})
        for region in complete:
            start, end = int(region["Start0"]), int(region["End0"])
            overlaps = []
            edge_predictions = 0
            for prediction in mapped:
                left = max(start, int(prediction["GenomeStart0"]))
                right = min(end, int(prediction["GenomeEnd0"]))
                if right > left:
                    overlaps.append((left, right))
                    edge_predictions += prediction["AtSequenceEdge"] == "yes"
            covered = union_length(overlaps)
            fraction = covered / max(1, end - start)
            recovery = "Recovered >=80%" if fraction >= 0.8 else "Partial 20-<80%" if fraction >= 0.2 else "Missed <20%"
            crossings = math.floor((end - 1) / 20_000) - math.floor(start / 20_000)
            rows.append({
                "Tool": tool,
                "FullRegionID": region["RegionID"],
                "FullStart0": start,
                "FullEnd0": end,
                "FullLengthBp": end - start,
                "FragmentPredictionsOverlapping": len(overlaps),
                "FragmentPredictionsAtEdge": edge_predictions,
                "FragmentUnionCoverageBp": covered,
                "FragmentUnionCoverageFraction": f"{fraction:.8f}",
                "RecoveryClass": recovery,
                "FragmentBoundaryCrossings": crossings,
            })
    return rows


def member_to_representative(root: Path) -> dict[str, str]:
    path = root / "data/small/34-nonredundant-gene-catalog-frozen/primary-catalog-membership.tsv.gz"
    mapping = {}
    for row in read_tsv(path):
        if row["MemberBranch"] == "megahit-co":
            mapping[row["MemberID"]] = row["RepresentativeID"]
    return mapping


def annotate_coassembly_abundance(
    root: Path,
    calls: list[dict[str, object]],
    anti_genes: dict[str, list[dict[str, object]]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    regions = [row for row in calls if row["Dataset"] == "coassembly-ge20kb" and row["Tool"] == "antiSMASH"]
    member_map = member_to_representative(root)
    abundance, denominators = read_abundance(root)
    memberships: dict[str, list[dict[str, object]]] = defaultdict(list)
    for region in regions:
        for gene in anti_genes.get(str(region["SequenceID"]), []):
            midpoint = (int(gene["Start0"]) + int(gene["End0"])) / 2
            if int(region["Start0"]) <= midpoint < int(region["End0"]):
                memberships[str(gene["GeneID"])].append(region)
    region_stats = {str(row["RegionID"]): {"Genes": 0.0, "Mapped": 0.0, "MOCK1": 0.0, "MOCK2": 0.0} for row in regions}
    gene_rows = []
    missing = set()
    for gene_id, assigned_regions in memberships.items():
        fraction = 1 / len(assigned_regions)
        representative = member_map.get(gene_id, "")
        if not representative:
            missing.add(gene_id)
        counts = abundance.get(representative, {}) if representative else {}
        for region in assigned_regions:
            stats = region_stats[str(region["RegionID"])]
            stats["Genes"] += fraction
            stats["Mapped"] += fraction if representative else 0
            stats["MOCK1"] += fraction * counts.get("MOCK1", 0)
            stats["MOCK2"] += fraction * counts.get("MOCK2", 0)
            gene_rows.append({
                "RegionID": region["RegionID"], "MemberGeneID": gene_id,
                "RepresentativeID": representative, "AssignmentFraction": f"{fraction:.8f}",
                "MOCK1FractionalRawReads": f"{fraction * counts.get('MOCK1', 0):.8f}",
                "MOCK2FractionalRawReads": f"{fraction * counts.get('MOCK2', 0):.8f}",
            })
    output = []
    for region in regions:
        stats = region_stats[str(region["RegionID"])]
        output.append({
            **region,
            "FractionalGenes": f"{stats['Genes']:.8f}",
            "MappedFractionalGenes": f"{stats['Mapped']:.8f}",
            "MOCK1FractionalRawReads": f"{stats['MOCK1']:.8f}",
            "MOCK2FractionalRawReads": f"{stats['MOCK2']:.8f}",
            "MOCK1AssignedReadPercent": f"{100 * stats['MOCK1'] / denominators['MOCK1']:.10f}",
            "MOCK2AssignedReadPercent": f"{100 * stats['MOCK2'] / denominators['MOCK2']:.10f}",
        })
    audit = {
        "genes_in_antismash_regions": len(memberships),
        "genes_missing_catalog_membership": len(missing),
        "sample_assigned_reads": denominators,
        "region_fractional_genes": sum(float(row["FractionalGenes"]) for row in output),
        "region_mock1_raw_reads": sum(float(row["MOCK1FractionalRawReads"]) for row in output),
        "region_mock2_raw_reads": sum(float(row["MOCK2FractionalRawReads"]) for row in output),
    }
    return output, gene_rows, audit


def type_summary(regions: list[dict[str, object]]) -> list[dict[str, object]]:
    totals: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for region in regions:
        categories = str(region["Categories"]).split(";")
        weight = 1 / len(categories)
        for category in categories:
            totals[category]["Regions"] += weight
            totals[category]["Genes"] += weight * float(region["FractionalGenes"])
            totals[category]["MOCK1"] += weight * float(region["MOCK1FractionalRawReads"])
            totals[category]["MOCK2"] += weight * float(region["MOCK2FractionalRawReads"])
    rows = []
    for category in CATEGORY_ORDER:
        if category not in totals:
            continue
        values = totals[category]
        rows.append({
            "Category": category,
            "FractionalRegions": f"{values['Regions']:.8f}",
            "FractionalGenes": f"{values['Genes']:.8f}",
            "MOCK1FractionalRawReads": f"{values['MOCK1']:.8f}",
            "MOCK2FractionalRawReads": f"{values['MOCK2']:.8f}",
        })
    return rows


def similarity_summary(regions: list[dict[str, object]]) -> list[dict[str, object]]:
    labels = ("High (>=0.75)", "Intermediate (0.50-<0.75)", "Low (<0.50)", "No score")
    totals: dict[str, dict[str, float]] = {label: defaultdict(float) for label in labels}
    for region in regions:
        value = region["TopMIBiGRiQ"]
        score = float(value) if value != "" else None
        label = "No score" if score is None else "High (>=0.75)" if score >= 0.75 else "Intermediate (0.50-<0.75)" if score >= 0.5 else "Low (<0.50)"
        totals[label]["Regions"] += 1
        totals[label]["MOCK1"] += float(region["MOCK1FractionalRawReads"])
        totals[label]["MOCK2"] += float(region["MOCK2FractionalRawReads"])
    return [{
        "SimilarityClass": label,
        "Regions": int(totals[label]["Regions"]),
        "MOCK1FractionalRawReads": f"{totals[label]['MOCK1']:.8f}",
        "MOCK2FractionalRawReads": f"{totals[label]['MOCK2']:.8f}",
    } for label in labels]


def main() -> int:
    args = parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    if not (work / ".article40-run-complete").is_file():
        raise FileNotFoundError("Run Article 40 antiSMASH/GECCO workflow first")
    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    all_calls: list[dict[str, object]] = []
    genes_by_dataset: dict[str, dict[str, list[dict[str, object]]]] = {}
    output_audit = []
    for dataset in DATASETS:
        anti_calls, genes = parse_antismash(work, dataset)
        gecco_calls = parse_gecco(work, dataset)
        all_calls.extend(anti_calls)
        all_calls.extend(gecco_calls)
        genes_by_dataset[dataset] = genes
        for tool, path in (
            ("antiSMASH", work / "antismash" / dataset / f"{dataset}.json"),
            ("GECCO", work / "gecco" / dataset / f"{dataset}.clusters.tsv"),
        ):
            output_audit.append({"Dataset": dataset, "Tool": tool, "File": str(path.relative_to(work)), "Bytes": path.stat().st_size, "SHA256": sha256(path)})

    yield_rows = tool_yield(all_calls)
    fragment_rows = fragmentation_summary(all_calls, read_tsv(work / "inputs/salinispora-fragment-map.tsv"))
    coassembly_regions, gene_rows, mapping_audit = annotate_coassembly_abundance(root, all_calls, genes_by_dataset["coassembly-ge20kb"])
    type_rows = type_summary(coassembly_regions)
    similarity_rows = similarity_summary(coassembly_regions)
    positive_rows = sorted(
        [
            row for row in all_calls
            if row["Dataset"] == "salinispora-full"
            and row["Tool"] == "antiSMASH"
            and row["TopKnownClusterSimilarityPercent"] != ""
        ],
        key=lambda row: (-int(row["TopKnownClusterSimilarityPercent"]), str(row["RegionID"])),
    )
    resource_rows = [parse_time(path) for path in sorted((work / "logs").glob("*.time.txt"))]

    write_tsv(summary_dir / "bgc-region-calls.tsv.gz", all_calls)
    write_tsv(summary_dir / "tool-yield-summary.tsv", yield_rows)
    write_tsv(summary_dir / "fragmentation-sensitivity.tsv", fragment_rows)
    write_tsv(summary_dir / "coassembly-bgc-abundance.tsv", coassembly_regions)
    write_tsv(summary_dir / "gene-to-bgc.tsv.gz", gene_rows)
    write_tsv(summary_dir / "bgc-type-summary.tsv", type_rows)
    write_tsv(summary_dir / "mibig-similarity-summary.tsv", similarity_rows)
    write_tsv(summary_dir / "salinispora-positive-control.tsv", positive_rows)
    write_tsv(summary_dir / "tool-output-audit.tsv", output_audit)
    write_tsv(summary_dir / "resource-usage.tsv", resource_rows)

    counts = {(row["Dataset"], row["Tool"]): 0 for row in all_calls}
    for row in all_calls:
        counts[(str(row["Dataset"]), str(row["Tool"]))] += 1
    fragment_counts = Counter((row["Tool"], row["RecoveryClass"]) for row in fragment_rows)
    co_anti = counts[("coassembly-ge20kb", "antiSMASH")]
    co_gecco = counts[("coassembly-ge20kb", "GECCO")]
    salinosporamide_full = next(
        row for row in positive_rows if "salinosporamide A" in str(row["TopKnownClusterDescription"])
    )
    salinosporamide_fragment = max(
        (
            row for row in all_calls
            if row["Dataset"] == "salinispora-fragmented"
            and row["Tool"] == "antiSMASH"
            and "salinosporamide A" in str(row["TopKnownClusterDescription"])
        ),
        key=lambda row: int(row["TopKnownClusterSimilarityPercent"]),
    )
    summary = {
        "article": 40,
        "seed": SEED,
        "antismash_regions": {dataset: counts[(dataset, "antiSMASH")] for dataset in DATASETS},
        "gecco_regions": {dataset: counts[(dataset, "GECCO")] for dataset in DATASETS},
        "salinispora_fragmentation": {
            tool: {label: fragment_counts[(tool, label)] for label in ("Recovered >=80%", "Partial 20-<80%", "Missed <20%")}
            for tool in ("antiSMASH", "GECCO")
        },
        "salinosporamide_knowncluster_similarity_percent": {
            "complete_genome": int(salinosporamide_full["TopKnownClusterSimilarityPercent"]),
            "twenty_kb_fragment": int(salinosporamide_fragment["TopKnownClusterSimilarityPercent"]),
        },
        "coassembly_antismash_regions": co_anti,
        "coassembly_gecco_regions": co_gecco,
        "coassembly_antismash_genes": mapping_audit["genes_in_antismash_regions"],
        "coassembly_genes_missing_catalog_membership": mapping_audit["genes_missing_catalog_membership"],
        "coassembly_primary_raw_reads": {"MOCK1": round(mapping_audit["region_mock1_raw_reads"], 8), "MOCK2": round(mapping_audit["region_mock2_raw_reads"], 8)},
        "sample_assigned_reads": mapping_audit["sample_assigned_reads"],
        "database_versions": {"PFAM": "35.0", "MIBiG": "4.0", "MITE": "1.3"},
        "reciprocal_overlap_threshold": 0.25,
        "mibig_similarity_is_not_compound_novelty": True,
        "bgc_presence_is_not_expression_or_metabolite_production": True,
        "fragmentation_changes_bgc_boundaries": True,
    }
    dump_json(summary_dir / "run-summary.json", summary)
    (work / ".article40-summary-complete").write_text(json.dumps(summary, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
