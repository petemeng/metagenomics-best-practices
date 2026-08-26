#!/usr/bin/env python3
"""Summarize Article 34 gene prediction, clustering, and exact-truth evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator, TextIO


SEED = 20260734
PRIMARY_CATALOGS = {
    "megahit-individual-primary",
    "megahit-co-primary",
    "megahit-mix-primary",
    "metaspades-individual-primary",
    "metaspades-co-primary",
    "metaspades-mix-primary",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def raw_pairs(path: Path) -> Iterator[tuple[str, str]]:
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Expected at least two columns in {path}:{line_number}")
            yield fields[0], fields[1]


def deterministic_gzip_text(path: Path) -> TextIO:
    raw = path.open("wb")
    gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(gz, encoding="utf-8", newline="")


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = deterministic_gzip_text(path) if path.suffix == ".gz" else path.open("w", encoding="utf-8", newline="")
    with handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def pct(numerator: int | float, denominator: int | float) -> str:
    return "NA" if denominator == 0 else f"{100 * numerator / denominator:.6f}"


def median(values: list[int]) -> str:
    return f"{statistics.median(values):.3f}" if values else "NA"


def load_membership(path: Path) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = defaultdict(list)
    seen: set[str] = set()
    for representative, member in raw_pairs(path):
        if member in seen:
            raise ValueError(f"Member {member} occurs more than once in {path}")
        clusters[representative].append(member)
        seen.add(member)
    return dict(clusters)


def expand_membership(work: Path, contract: dict[str, str]) -> dict[str, list[str]]:
    final = load_membership(work / contract["MembershipTSV"])
    stage1_path = contract.get("Stage1MembershipTSV", "")
    if not stage1_path:
        return final
    stage1 = load_membership(work / stage1_path)
    expanded: dict[str, list[str]] = {}
    for representative, members in final.items():
        raw_members: list[str] = []
        for member in members:
            raw_members.extend(stage1.get(member, [member]))
        if len(raw_members) != len(set(raw_members)):
            raise ValueError(f"Expanded cluster {representative} contains duplicate raw genes")
        expanded[representative] = raw_members
    return expanded


def branch_prediction_summary(metadata: list[dict[str, str]]) -> list[dict]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in metadata:
        grouped[row["Branch"]].append(row)
    rows: list[dict] = []
    for branch in sorted(grouped):
        genes = grouped[branch]
        classes = Counter(row["Completeness"] for row in genes)
        lengths = [int(row["NtLength"]) for row in genes]
        rows.append(
            {
                "Branch": branch,
                "Assembler": genes[0]["Assembler"],
                "OriginStrategy": genes[0]["OriginStrategy"],
                "Mock": genes[0]["Mock"],
                "Genes": len(genes),
                "CompleteGenes": classes["Complete"],
                "PartialGenes": classes["Partial"],
                "IncompleteGenes": classes["Incomplete"],
                "CompletePct": pct(classes["Complete"], len(genes)),
                "MedianGeneBp": median(lengths),
                "MeanGeneBp": f"{statistics.fmean(lengths):.3f}",
                "TotalGeneBp": sum(lengths),
            }
        )
    return rows


def size_bin(size: int) -> str:
    if size == 1:
        return "1"
    if size == 2:
        return "2"
    if size <= 5:
        return "3-5"
    if size <= 10:
        return "6-10"
    if size <= 50:
        return "11-50"
    return ">50"


def length_bin(length: int) -> tuple[str, int]:
    if length > 2_500:
        return ">2500", 2_501
    lower = (max(length, 1) - 1) // 50 * 50 + 1
    upper = lower + 49
    return f"{lower}-{upper}", lower


def catalog_summaries(
    work: Path,
    contracts: list[dict[str, str]],
    metadata: dict[str, dict[str, str]],
) -> tuple[list[dict], list[dict], list[dict], dict[str, dict[str, list[str]]]]:
    summary_rows: list[dict] = []
    size_rows: list[dict] = []
    length_rows: list[dict] = []
    memberships: dict[str, dict[str, list[str]]] = {}
    for contract in contracts:
        catalog_id = contract["CatalogID"]
        clusters = expand_membership(work, contract)
        memberships[catalog_id] = clusters
        representatives = set(clusters)
        missing = representatives - set(metadata)
        if missing:
            raise ValueError(f"Missing metadata for {len(missing)} representatives in {catalog_id}")
        raw_members = {member for members in clusters.values() for member in members}
        sizes = [len(members) for members in clusters.values()]
        rep_rows = [metadata[identifier] for identifier in representatives]
        classes = Counter(row["Completeness"] for row in rep_rows)
        origins = Counter(row["OriginStrategy"] for row in rep_rows)
        mixed_clusters = sum(
            len({metadata[member]["OriginStrategy"] for member in members}) > 1 for members in clusters.values()
        )
        lengths = [int(row["NtLength"]) for row in rep_rows]
        summary_rows.append(
            {
                **contract,
                "RawInputGenes": len(raw_members),
                "CatalogGenes": len(representatives),
                "RemovedRedundantGenes": len(raw_members) - len(representatives),
                "CompressionPct": pct(len(raw_members) - len(representatives), len(raw_members)),
                "SingletonClusters": sum(size == 1 for size in sizes),
                "SingletonPct": pct(sum(size == 1 for size in sizes), len(sizes)),
                "LargestCluster": max(sizes),
                "MedianClusterSize": median(sizes),
                "CompleteGenes": classes["Complete"],
                "PartialGenes": classes["Partial"],
                "IncompleteGenes": classes["Incomplete"],
                "CompletePct": pct(classes["Complete"], len(representatives)),
                "MedianGeneBp": median(lengths),
                "TotalGeneBp": sum(lengths),
                "RepresentativeFromIndividual": origins["Individual"],
                "RepresentativeFromCoassembly": origins["Co-assembly"],
                "MixedOriginClusters": mixed_clusters,
                "RepresentativeFAA_SHA256": file_sha256(work / contract["RepresentativesFAA"]),
                "RepresentativeFNA_SHA256": file_sha256(work / contract["RepresentativesFNA"]),
            }
        )
        counts = Counter(size_bin(size) for size in sizes)
        for label in ("1", "2", "3-5", "6-10", "11-50", ">50"):
            size_rows.append(
                {
                    "CatalogID": catalog_id,
                    "Assembler": contract["Assembler"],
                    "Strategy": contract["Strategy"],
                    "ClusterSizeBin": label,
                    "Clusters": counts[label],
                    "ClusterPct": pct(counts[label], len(sizes)),
                }
            )
        if catalog_id in PRIMARY_CATALOGS:
            hist = Counter()
            for row in rep_rows:
                label, start = length_bin(int(row["NtLength"]))
                hist[(row["Completeness"], label, start)] += 1
            for (completeness, label, start), count in sorted(hist.items(), key=lambda item: (item[0][2], item[0][0])):
                length_rows.append(
                    {
                        "CatalogID": catalog_id,
                        "Assembler": contract["Assembler"],
                        "Strategy": contract["Strategy"],
                        "Completeness": completeness,
                        "LengthBin": label,
                        "LengthBinStart": start,
                        "Genes": count,
                    }
                )
    return summary_rows, size_rows, length_rows, memberships


def file_sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_hits(path: Path) -> list[tuple[str, str]]:
    return [(query, target) for query, target in raw_pairs(path)]


def truth_summaries(
    work: Path,
    contracts: list[dict[str, str]],
    memberships: dict[str, dict[str, list[str]]],
    truth_meta: dict[str, dict[str, str]],
    truth_genomes: list[dict[str, str]],
) -> tuple[list[dict], list[dict]]:
    truth_clusters = load_membership(work / "truth/nonredundant-primary/membership.tsv")
    truth_member_to_rep = {member: rep for rep, members in truth_clusters.items() for member in members}
    truth_cluster_classes = Counter(truth_meta[rep]["Completeness"] for rep in truth_clusters)
    genes_by_genome = Counter(row["GenBankAssembly"] for row in truth_meta.values())
    genome_info = {row["GenBankAssembly"]: row for row in truth_genomes}
    summary_rows: list[dict] = []
    genome_rows: list[dict] = []

    for contract in contracts:
        catalog_id = contract["CatalogID"]
        forward = read_hits(work / f"truth-audit/{catalog_id}/catalog-to-truth-recovery.tsv")
        support = read_hits(work / f"truth-audit/{catalog_id}/catalog-to-truth-support.tsv")
        hit_truth_genes = {target for _, target in forward}
        hit_truth_clusters = {truth_member_to_rep[target] for target in hit_truth_genes}
        hit_catalog_reps = {query for query, _ in support}
        catalog_reps = set(memberships[catalog_id])
        if not hit_catalog_reps <= catalog_reps:
            raise ValueError(f"Truth support search returned unknown representatives for {catalog_id}")
        recovered_classes = Counter(truth_meta[rep]["Completeness"] for rep in hit_truth_clusters)
        summary_rows.append(
            {
                "CatalogID": catalog_id,
                "Assembler": contract["Assembler"],
                "Strategy": contract["Strategy"],
                "Method": contract["Method"],
                "TruthGenes": len(truth_meta),
                "RecoveredTruthGenes": len(hit_truth_genes),
                "TruthGeneRecoveryPct": pct(len(hit_truth_genes), len(truth_meta)),
                "TruthNRClusters": len(truth_clusters),
                "RecoveredTruthNRClusters": len(hit_truth_clusters),
                "TruthNRRecoveryPct": pct(len(hit_truth_clusters), len(truth_clusters)),
                "TruthCompleteClusters": truth_cluster_classes["Complete"],
                "RecoveredTruthCompleteClusters": recovered_classes["Complete"],
                "TruthCompleteClusterRecoveryPct": pct(recovered_classes["Complete"], truth_cluster_classes["Complete"]),
                "CatalogGenes": len(catalog_reps),
                "TruthSupportedCatalogGenes": len(hit_catalog_reps),
                "CatalogTruthSupportPct": pct(len(hit_catalog_reps), len(catalog_reps)),
            }
        )
        recovered_by_genome = Counter(truth_meta[gene]["GenBankAssembly"] for gene in hit_truth_genes)
        for accession in sorted(genome_info):
            denominator = genes_by_genome[accession]
            recovered = recovered_by_genome[accession]
            info = genome_info[accession]
            genome_rows.append(
                {
                    "CatalogID": catalog_id,
                    "Assembler": contract["Assembler"],
                    "Strategy": contract["Strategy"],
                    "Method": contract["Method"],
                    "GenBankAssembly": accession,
                    "Reference": info["Reference"],
                    "InMOCK1": info["InMOCK1"],
                    "ExpectedAbundanceMOCK1Pct": info["ExpectedAbundanceMOCK1Pct"],
                    "ExpectedAbundanceMOCK2Pct": info["ExpectedAbundanceMOCK2Pct"],
                    "ExpectedAbundanceMeanPct": info["ExpectedAbundanceMeanPct"],
                    "TruthGenes": denominator,
                    "RecoveredTruthGenes": recovered,
                    "GeneRecoveryPct": pct(recovered, denominator),
                }
            )
    return summary_rows, genome_rows


def mix_origin_summary(
    contracts: list[dict[str, str]],
    memberships: dict[str, dict[str, list[str]]],
    metadata: dict[str, dict[str, str]],
) -> list[dict]:
    rows: list[dict] = []
    for contract in contracts:
        catalog_id = contract["CatalogID"]
        if catalog_id not in {"megahit-mix-primary", "metaspades-mix-primary"}:
            continue
        clusters = memberships[catalog_id]
        grouped: Counter[tuple[str, str]] = Counter()
        for representative in clusters:
            item = metadata[representative]
            grouped[(item["OriginStrategy"], item["Completeness"])] += 1
        total = len(clusters)
        for origin in ("Individual", "Co-assembly"):
            for completeness in ("Complete", "Partial", "Incomplete"):
                count = grouped[(origin, completeness)]
                rows.append(
                    {
                        "CatalogID": catalog_id,
                        "Assembler": contract["Assembler"],
                        "RepresentativeOrigin": origin,
                        "Completeness": completeness,
                        "Genes": count,
                        "CatalogPct": pct(count, total),
                    }
                )
    return rows


def pair_count(size: int) -> int:
    return size * (size - 1) // 2


def method_agreement(memberships: dict[str, dict[str, list[str]]]) -> list[dict]:
    left_id = "megahit-mix-primary"
    right_id = "megahit-mix-cdhit"
    left = memberships[left_id]
    right = memberships[right_id]
    left_label = {member: rep for rep, members in left.items() for member in members}
    right_label = {member: rep for rep, members in right.items() for member in members}
    if set(left_label) != set(right_label):
        raise ValueError("MMseqs2 and CD-HIT do not cover the same expanded raw gene set")
    contingency = Counter((left_label[member], right_label[member]) for member in left_label)
    same_both = sum(pair_count(count) for count in contingency.values())
    same_left = sum(pair_count(len(members)) for members in left.values())
    same_right = sum(pair_count(len(members)) for members in right.values())
    precision = same_both / same_right if same_right else 1.0
    recall = same_both / same_left if same_left else 1.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    reps_left = set(left)
    reps_right = set(right)
    return [
        {
            "MethodA": "MMseqs2 9.d36de",
            "MethodB": "CD-HIT 4.8.1 local identity",
            "RawGenes": len(left_label),
            "ClustersA": len(left),
            "ClustersB": len(right),
            "RepresentativeIntersection": len(reps_left & reps_right),
            "RepresentativeJaccard": f"{len(reps_left & reps_right) / len(reps_left | reps_right):.9f}",
            "CoClusterPairPrecision": f"{precision:.9f}",
            "CoClusterPairRecall": f"{recall:.9f}",
            "CoClusterPairF1": f"{f1:.9f}",
        }
    ]


def parse_elapsed(value: str) -> float:
    parts = value.split(":")
    if len(parts) == 2:
        minutes, seconds = parts
        return int(minutes) * 60 + float(seconds)
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    raise ValueError(f"Unexpected elapsed time: {value}")


def resource_usage(log_dir: Path, valid_steps: set[str]) -> list[dict]:
    rows: list[dict] = []
    for path in sorted(log_dir.glob("*.time.txt")):
        if path.name.removesuffix(".time.txt") not in valid_steps:
            continue
        if path.stat().st_size == 0:
            continue
        values: dict[str, str] = {}
        for raw in path.read_text(encoding="utf-8").splitlines():
            if ": " in raw:
                key, value = raw.strip().split(": ", 1)
                values[key] = value
        elapsed_key = "Elapsed (wall clock) time (h:mm:ss or m:ss)"
        rows.append(
            {
                "Step": path.name.removesuffix(".time.txt"),
                "UserSeconds": values.get("User time (seconds)", "NA"),
                "SystemSeconds": values.get("System time (seconds)", "NA"),
                "ElapsedSeconds": f"{parse_elapsed(values[elapsed_key]):.3f}" if elapsed_key in values else "NA",
                "MaximumRSSKiB": values.get("Maximum resident set size (kbytes)", "NA"),
                "FileSystemInputs": values.get("File system inputs", "NA"),
                "FileSystemOutputs": values.get("File system outputs", "NA"),
                "ExitStatus": values.get("Exit status", "NA"),
            }
        )
    return rows


def primary_membership_rows(
    membership: dict[str, list[str]], metadata: dict[str, dict[str, str]]
) -> tuple[list[dict], list[dict]]:
    members: list[dict] = []
    reps: list[dict] = []
    for representative in sorted(membership):
        raw_members = membership[representative]
        rep = metadata[representative]
        reps.append(
            {
                "RepresentativeID": representative,
                "RepresentativeOrigin": rep["OriginStrategy"],
                "Completeness": rep["Completeness"],
                "PartialCode": rep["PartialCode"],
                "NtLength": rep["NtLength"],
                "AaLength": rep["AaLength"],
                "ClusterSize": len(raw_members),
                "IndividualMembers": sum(metadata[member]["OriginStrategy"] == "Individual" for member in raw_members),
                "CoassemblyMembers": sum(metadata[member]["OriginStrategy"] == "Co-assembly" for member in raw_members),
            }
        )
        for member in sorted(raw_members):
            members.append(
                {
                    "RepresentativeID": representative,
                    "MemberID": member,
                    "MemberOrigin": metadata[member]["OriginStrategy"],
                    "MemberBranch": metadata[member]["Branch"],
                }
            )
    return members, reps


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article34-run-complete").is_file():
        raise FileNotFoundError("Article 34 run is incomplete")
    summary = work / "summary"
    assembly_rows = read_tsv(summary / "gene-metadata.tsv.gz")
    truth_rows = read_tsv(summary / "truth-gene-metadata.tsv.gz")
    assembly_meta = {row["GeneID"]: row for row in assembly_rows}
    truth_meta = {row["GeneID"]: row for row in truth_rows}
    contracts = read_tsv(summary / "catalog-contracts.tsv")

    prediction = branch_prediction_summary(assembly_rows + truth_rows)
    catalog, cluster_bins, length_hist, memberships = catalog_summaries(work, contracts, assembly_meta)
    truth, per_genome = truth_summaries(
        work,
        contracts,
        memberships,
        truth_meta,
        read_tsv(summary / "truth-genomes.tsv"),
    )
    origins = mix_origin_summary(contracts, memberships, assembly_meta)
    agreement = method_agreement(memberships)
    command_rows = read_tsv(summary / "command-log.tsv")
    valid_steps = {re.sub(r"[^A-Za-z0-9_.-]+", "-", row["Step"]) for row in command_rows}
    resources = resource_usage(work / "logs", valid_steps)
    primary_members, primary_reps = primary_membership_rows(memberships["megahit-mix-primary"], assembly_meta)

    write_tsv(summary / "gene-prediction-summary.tsv", prediction)
    write_tsv(summary / "catalog-summary.tsv", catalog)
    write_tsv(summary / "truth-audit-summary.tsv", truth)
    write_tsv(summary / "per-genome-recovery.tsv", per_genome)
    write_tsv(summary / "mix-origin-summary.tsv", origins)
    write_tsv(summary / "cluster-size-bins.tsv", cluster_bins)
    write_tsv(summary / "gene-length-histogram.tsv", length_hist)
    write_tsv(summary / "method-agreement.tsv", agreement)
    write_tsv(summary / "resource-usage.tsv", resources)
    write_tsv(summary / "primary-catalog-membership.tsv.gz", primary_members)
    write_tsv(summary / "primary-catalog-representatives.tsv.gz", primary_reps)

    primary_catalog = next(row for row in catalog if row["CatalogID"] == "megahit-mix-primary")
    primary_truth = next(row for row in truth if row["CatalogID"] == "megahit-mix-primary")
    run_summary = {
        "seed": SEED,
        "predicted_assembly_genes": len(assembly_rows),
        "predicted_truth_genes": len(truth_rows),
        "truth_nr_clusters": int(primary_truth["TruthNRClusters"]),
        "catalogs_audited": len(contracts),
        "primary_catalog": {
            "catalog_id": "megahit-mix-primary",
            "raw_input_genes": int(primary_catalog["RawInputGenes"]),
            "catalog_genes": int(primary_catalog["CatalogGenes"]),
            "compression_pct": float(primary_catalog["CompressionPct"]),
            "complete_genes": int(primary_catalog["CompleteGenes"]),
            "truth_nr_recovery_pct": float(primary_truth["TruthNRRecoveryPct"]),
            "truth_support_pct": float(primary_truth["CatalogTruthSupportPct"]),
            "representatives_fna_sha256": primary_catalog["RepresentativeFNA_SHA256"],
        },
        "method_agreement": agreement[0],
        "resource_steps": len(resources),
    }
    (summary / "run-summary.json").write_text(json.dumps(run_summary, indent=2) + "\n", encoding="utf-8")
    (work / ".article34-summary-complete").write_text(json.dumps(run_summary, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(run_summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
