#!/usr/bin/env python3
"""Summarize Article 45 dRep outputs without using mock truth."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

from article41_44_utils import dump_json, parse_time, read_tsv, write_tsv


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    inputs = {row["Genome"]: row for row in read_tsv(work / "input-genomes.tsv")}
    if len(inputs) != 124:
        raise ValueError("Article 45 input ledger must contain 124 genomes")

    branch_specs = [("Species 95% ANI", "species95", 0.95), ("Near-clone 99.9% ANI", "nearclone999", 0.999)]
    all_membership: list[dict[str, object]] = []
    all_representatives: list[dict[str, object]] = []
    threshold_rows: list[dict[str, object]] = []
    source_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []
    branch_maps: dict[str, dict[str, str]] = {}

    for branch_label, branch_dir, threshold in branch_specs:
        tables = work / "drep" / branch_dir / "data_tables"
        cdb = read_csv(tables / "Cdb.csv")
        wdb = read_csv(tables / "Wdb.csv")
        sdb = {row["genome"]: float(row["score"]) for row in read_csv(tables / "Sdb.csv")}
        ndb = read_csv(tables / "Ndb.csv")
        if len(cdb) != len(inputs):
            raise ValueError(f"{branch_dir}: expected 124 Cdb rows")
        representative = {row["cluster"]: row["genome"] for row in wdb}
        cluster_map = {row["genome"]: row["secondary_cluster"] for row in cdb}
        branch_maps[branch_dir] = cluster_map
        counts = Counter(cluster_map.values())
        if set(representative) != set(counts):
            raise ValueError(f"{branch_dir}: Wdb/Cdb cluster mismatch")

        for genome in sorted(inputs):
            source = inputs[genome]
            cluster = cluster_map[genome]
            rep = representative[cluster]
            all_membership.append(
                {
                    "Branch": branch_label,
                    "ANIThresholdPct": threshold * 100,
                    "Genome": genome,
                    "Cluster": cluster,
                    "Representative": rep,
                    "IsRepresentative": genome == rep,
                    "ClusterSize": counts[cluster],
                    "SourceStage": source["SourceStage"],
                    "SourceBranch": source["SourceBranch"],
                    "Completeness": float(source["Completeness"]),
                    "Contamination": float(source["Contamination"]),
                    "N50Bp": int(source["N50Bp"]),
                    "GenomeBp": int(source["GenomeBp"]),
                    "SHA256": source["SHA256"],
                    "MIMAGQuality": source["MIMAGQuality"],
                    "dRepScore": sdb[genome],
                }
            )

        for cluster in sorted(counts, key=lambda value: tuple(int(x) for x in value.split("_"))):
            members = sorted(genome for genome, assigned in cluster_map.items() if assigned == cluster)
            rep = representative[cluster]
            source = inputs[rep]
            all_representatives.append(
                {
                    "Branch": branch_label,
                    "ANIThresholdPct": threshold * 100,
                    "Cluster": cluster,
                    "Representative": rep,
                    "Members": counts[cluster],
                    "MemberGenomes": ";".join(members),
                    "MemberSourceStages": ";".join(sorted({inputs[g]["SourceStage"] for g in members})),
                    "RepresentativeSourceStage": source["SourceStage"],
                    "RepresentativeSourceBranch": source["SourceBranch"],
                    "Completeness": float(source["Completeness"]),
                    "Contamination": float(source["Contamination"]),
                    "N50Bp": int(source["N50Bp"]),
                    "GenomeBp": int(source["GenomeBp"]),
                    "MIMAGQuality": source["MIMAGQuality"],
                    "dRepScore": sdb[rep],
                }
            )

        threshold_rows.append(
            {
                "Branch": branch_label,
                "ANIThresholdPct": threshold * 100,
                "MinimumAlignmentFractionPct": 30,
                "InputGenomes": len(inputs),
                "Clusters": len(counts),
                "Representatives": len(representative),
                "GenomesRemoved": len(inputs) - len(representative),
                "SingletonClusters": sum(size == 1 for size in counts.values()),
                "MultiGenomeClusters": sum(size > 1 for size in counts.values()),
                "LargestCluster": max(counts.values()),
                "RepresentativesFromArticle44": sum(inputs[g]["SourceStage"] == "Article44-selected" for g in representative.values()),
            }
        )

        grouped_input = Counter(row["SourceBranch"] for row in inputs.values())
        grouped_reps = Counter(inputs[g]["SourceBranch"] for g in representative.values())
        for source_branch in sorted(grouped_input):
            source_rows.append(
                {
                    "Branch": branch_label,
                    "SourceBranch": source_branch,
                    "InputGenomes": grouped_input[source_branch],
                    "Representatives": grouped_reps[source_branch],
                    "RetentionPct": 100 * grouped_reps[source_branch] / grouped_input[source_branch],
                }
            )

        seen: set[tuple[str, str]] = set()
        for row in ndb:
            first, second = sorted((row["reference"], row["querry"]))
            if first == second or (first, second) in seen:
                continue
            seen.add((first, second))
            pair_rows.append(
                {
                    "Branch": branch_label,
                    "Genome1": first,
                    "Genome2": second,
                    "ANIPct": 100 * float(row["ani"]),
                    "AlignmentFractionPct": 100 * float(row["alignment_coverage"]),
                    "SameCluster": cluster_map[first] == cluster_map[second],
                    "Cluster1": cluster_map[first],
                    "Cluster2": cluster_map[second],
                }
            )

    stability_rows: list[dict[str, object]] = []
    species_groups: dict[str, list[str]] = defaultdict(list)
    for genome, cluster in branch_maps["species95"].items():
        species_groups[cluster].append(genome)
    near_map = branch_maps["nearclone999"]
    for cluster, members in sorted(species_groups.items()):
        near_clusters = sorted({near_map[g] for g in members})
        stability_rows.append(
            {
                "SpeciesCluster": cluster,
                "SpeciesMembers": len(members),
                "NearCloneClusters": len(near_clusters),
                "SplitsAt99_9Pct": len(near_clusters) > 1,
                "NearCloneClusterIDs": ";".join(near_clusters),
            }
        )

    write_tsv(summary_dir / "cluster-membership.tsv.gz", all_membership)
    write_tsv(summary_dir / "representative-genomes.tsv", all_representatives)
    write_tsv(summary_dir / "threshold-summary.tsv", threshold_rows)
    write_tsv(summary_dir / "source-retention.tsv", source_rows)
    write_tsv(summary_dir / "pairwise-ani.tsv.gz", pair_rows)
    write_tsv(summary_dir / "representative-stability.tsv", stability_rows)
    resources = [
        parse_time(work / "logs" / "drep-species95.time.txt"),
        parse_time(work / "logs" / "drep-nearclone999.time.txt"),
    ]
    write_tsv(summary_dir / "resource-summary.tsv", resources)

    result = {
        "article": 45,
        "input_genomes": len(inputs),
        "species_clusters": threshold_rows[0]["Clusters"],
        "near_clone_clusters": threshold_rows[1]["Clusters"],
        "species_clusters_split_at_99_9": sum(bool(row["SplitsAt99_9Pct"]) for row in stability_rows),
        "species_representatives_from_article44": threshold_rows[0]["RepresentativesFromArticle44"],
        "truth_used_for_clustering_or_selection": False,
    }
    dump_json(summary_dir / "run-summary.json", result)
    (work / ".article45-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
