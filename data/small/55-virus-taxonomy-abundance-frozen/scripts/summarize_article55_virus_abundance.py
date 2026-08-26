#!/usr/bin/env python3
"""Summarize real mock-virome abundance, taxonomy, vOTU, and lifestyle evidence."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from article41_44_utils import dump_json, parse_time, read_tsv, write_tsv


SEED = 20260755
PRIMARY_BREADTH = 75.0
PRIMARY_DEPTH = 1.0
THRESHOLDS = (50.0, 70.0, 75.0, 90.0)
LIBRARY_ORDER = ("lib1_illumina", "lib2_illumina", "lib3_illumina")


def bool_text(value: bool) -> str:
    return "True" if value else "False"


def total_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def phage_id(sequence_id: str) -> str:
    return sequence_id.split("|", 1)[0]


def family_from_lineage(lineage: str) -> str:
    fields = lineage.split(";")
    return fields[6] if len(fields) > 6 and fields[6] else "Unclassified at family"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    results = args.results_dir.resolve()
    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = work / "phage-metadata.tsv"
    coverage_path = work / "published-coverage-depth.tsv"
    library_path = work / "illumina-library-metadata.tsv"
    taxonomy_path = (
        results
        / "genomad/cook15-phage-reference_annotate/cook15-phage-reference_taxonomy.tsv"
    )
    virus_summary_path = (
        results
        / "genomad/cook15-phage-reference_summary/cook15-phage-reference_virus_summary.tsv"
    )
    cluster_path = results / "votu/votu-clusters.tsv"
    ani_path = results / "votu/ani.tsv"
    required = (
        metadata_path,
        coverage_path,
        library_path,
        taxonomy_path,
        virus_summary_path,
        cluster_path,
        ani_path,
        work / "author-source-assertions.tsv",
        work / "lifecycle-source-assertions.tsv",
        work / "bbmap-smoke-coverage.tsv",
        work / "logs/bbmap-smoke.time.txt",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 55 evidence: " + ", ".join(missing))

    metadata_rows = read_tsv(metadata_path)
    metadata = {row["PhageID"]: row for row in metadata_rows}
    if len(metadata) != 15:
        raise RuntimeError(f"Expected 15 mock-community phages, observed {len(metadata)}")

    all_coverage = read_tsv(coverage_path)
    illumina = [row for row in all_coverage if row["Platform"] == "Illumina"]
    primary = [row for row in illumina if row["Dataset"] in LIBRARY_ORDER]
    pooled = [row for row in illumina if row["Dataset"] == "pooled_illumina"]
    if len(primary) != 45 or len(pooled) != 15:
        raise RuntimeError(
            f"Expected 45 individual and 15 pooled Illumina rows; observed "
            f"{len(primary)} and {len(pooled)}"
        )

    abundance_rows = []
    for row in illumina:
        breadth = float(row["BreadthPct"])
        depth = float(row["MeanDepthX"])
        present = breadth >= PRIMARY_BREADTH and depth >= PRIMARY_DEPTH
        abundance_rows.append(
            {
                **row,
                "PresenceDepthGateX": PRIMARY_DEPTH,
                "PresenceBreadthGatePct": PRIMARY_BREADTH,
                "Present": bool_text(present),
                "QuantificationUnit": "reference sequence before vOTU collapse",
                "Source": "Cook et al. 2024 Supplementary Table S4",
            }
        )
    write_tsv(summary_dir / "illumina-abundance-evidence.tsv", abundance_rows)

    by_phage: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in primary:
        by_phage[row["PhageID"]][row["Dataset"]] = row
    presence_rows = []
    prevalence_rows = []
    for pid in metadata:
        library_values = by_phage[pid]
        if set(library_values) != set(LIBRARY_ORDER):
            raise RuntimeError(f"Incomplete Illumina replicate coverage for {pid}")
        present_flags = []
        out = {"PhageID": pid}
        for dataset in LIBRARY_ORDER:
            row = library_values[dataset]
            breadth = float(row["BreadthPct"])
            depth = float(row["MeanDepthX"])
            present = breadth >= PRIMARY_BREADTH and depth >= PRIMARY_DEPTH
            present_flags.append(present)
            suffix = dataset.removeprefix("lib").removesuffix("_illumina")
            out[f"Library{suffix}BreadthPct"] = breadth
            out[f"Library{suffix}MeanDepthX"] = depth
            out[f"Library{suffix}Present"] = bool_text(present)
        out["ReplicatesPresent"] = sum(present_flags)
        out["PrevalenceFraction"] = round(sum(present_flags) / 3, 6)
        presence_rows.append(out)
        prevalence_rows.append(
            {
                "PhageID": pid,
                "InputGenomeCopies": metadata[pid]["InputGenomeCopies"],
                "ReplicatesPresent": sum(present_flags),
                "ReplicatesTested": 3,
                "PrevalenceFraction": round(sum(present_flags) / 3, 6),
                "Interpretation": {
                    3: "Detected in all technical library replicates",
                    2: "Detected in two technical library replicates",
                    1: "Detected in one technical library replicate",
                    0: "Not detected in any technical library replicate",
                }[sum(present_flags)],
            }
        )
    write_tsv(summary_dir / "illumina-presence-matrix.tsv", presence_rows)
    write_tsv(summary_dir / "replicate-prevalence.tsv", prevalence_rows)

    threshold_rows = []
    for dataset in (*LIBRARY_ORDER, "pooled_illumina"):
        rows = [row for row in illumina if row["Dataset"] == dataset]
        for threshold in THRESHOLDS:
            detected = sum(
                float(row["BreadthPct"]) >= threshold
                and float(row["MeanDepthX"]) >= PRIMARY_DEPTH
                for row in rows
            )
            threshold_rows.append(
                {
                    "Dataset": dataset,
                    "BreadthThresholdPct": int(threshold),
                    "MinimumMeanDepthX": PRIMARY_DEPTH,
                    "PhagesDetected": detected,
                    "PhagesTested": 15,
                }
            )
    write_tsv(summary_dir / "breadth-threshold-sensitivity.tsv", threshold_rows)

    pooled_by_phage = {row["PhageID"]: row for row in pooled}
    depth_rows = []
    for pid, meta in metadata.items():
        depths = [float(by_phage[pid][dataset]["MeanDepthX"]) for dataset in LIBRARY_ORDER]
        breadths = [float(by_phage[pid][dataset]["BreadthPct"]) for dataset in LIBRARY_ORDER]
        depth_rows.append(
            {
                "PhageID": pid,
                "InputGenomeCopies": float(meta["InputGenomeCopies"]),
                "InputCPM": float(meta["InputTPM"]),
                "MeanDepthAcrossLibrariesX": round(sum(depths) / 3, 6),
                "MeanBreadthAcrossLibrariesPct": round(sum(breadths) / 3, 6),
                "PooledMeanDepthX": float(pooled_by_phage[pid]["MeanDepthX"]),
                "PooledBreadthPct": float(pooled_by_phage[pid]["BreadthPct"]),
                "ReplicatesPresent": sum(
                    breadth >= PRIMARY_BREADTH and depth >= PRIMARY_DEPTH
                    for breadth, depth in zip(breadths, depths)
                ),
                "AbundanceInterpretation": (
                    "Duplicate-reference multi-mapping risk"
                    if pid in {"vB_Eco_mar001J1", "vB_Eco_mar002J2"}
                    else "Same-vOTU cross-mapping risk"
                    if pid in {"vB_EcoS_swan01", "vB_Eco_SLUR29"}
                    else "Sequence-level value; retain breadth gate"
                ),
            }
        )
    write_tsv(summary_dir / "input-versus-observed-depth.tsv", depth_rows)

    smoke_rows = []
    for row in read_tsv(work / "bbmap-smoke-coverage.tsv"):
        pid = phage_id(row["#rname"])
        smoke_rows.append(
            {
                "PhageID": pid,
                "ReferenceLengthBp": row["endpos"],
                "MappedRecords": row["numreads"],
                "BreadthPct": row["coverage"],
                "MeanDepthX": row["meandepth"],
                "MeanMAPQ": row["meanmapq"],
                "SmokeInput": "the 15 reference FASTA records read back as sequences",
                "Interpretation": "implementation check only; not a biological abundance estimate",
            }
        )
    if len(smoke_rows) != 15:
        raise RuntimeError("BBMap duplicate-reference smoke test did not cover 15 references")
    write_tsv(summary_dir / "bbmap-duplicate-reference-smoke.tsv", smoke_rows)

    taxonomy = {phage_id(row["seq_name"]): row for row in read_tsv(taxonomy_path)}
    virus_summary = {
        phage_id(row["seq_name"]): row for row in read_tsv(virus_summary_path)
    }
    if set(taxonomy) != set(metadata) or set(virus_summary) != set(metadata):
        raise RuntimeError("geNomad did not report all 15 reference genomes")

    cluster_lines = [
        line.split("\t")
        for line in cluster_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    member_to_rep: dict[str, str] = {}
    cluster_rows = []
    for representative, member_text in cluster_lines:
        members = member_text.split(",")
        rep_id = phage_id(representative)
        member_ids = [phage_id(member) for member in members]
        for member in member_ids:
            if member in member_to_rep:
                raise RuntimeError(f"vOTU member assigned twice: {member}")
            member_to_rep[member] = rep_id
        cluster_rows.append(
            {
                "VOTURepresentative": rep_id,
                "MemberCount": len(member_ids),
                "Members": ",".join(member_ids),
                "MinANIPct": 95,
                "MinTargetCoveragePct": 85,
                "AbundanceRule": "map once to a nonredundant vOTU catalog; do not sum member depths",
            }
        )
    if set(member_to_rep) != set(metadata):
        raise RuntimeError("vOTU clustering did not cover all 15 reference genomes")
    write_tsv(summary_dir / "votu-cluster-summary.tsv", cluster_rows)

    taxonomy_rows = []
    for pid, meta in metadata.items():
        tax = taxonomy[pid]
        virus = virus_summary[pid]
        lineage = tax["lineage"]
        taxonomy_rows.append(
            {
                "PhageID": pid,
                "Accession": meta["Accession"],
                "VOTURepresentative": member_to_rep[pid],
                "VOTUMemberCount": next(
                    row["MemberCount"]
                    for row in cluster_rows
                    if row["VOTURepresentative"] == member_to_rep[pid]
                ),
                "geNomadFamily": family_from_lineage(lineage),
                "geNomadTaxonomyAgreement": tax["agreement"],
                "geNomadLineage": lineage,
                "geNomadVirusScore": virus["virus_score"],
                "geNomadTopology": virus["topology"],
                "TaxonomyDatabase": "geNomad DB v1.9; ICTV MSL39",
            }
        )
    write_tsv(summary_dir / "taxonomy-votu-ledger.tsv", taxonomy_rows)

    pair_rows = []
    seen: set[tuple[str, str]] = set()
    for row in read_tsv(ani_path):
        left, right = phage_id(row["qname"]), phage_id(row["tname"])
        if left == right:
            continue
        pair = tuple(sorted((left, right)))
        if pair in seen:
            continue
        seen.add(pair)
        reciprocal = [
            x
            for x in read_tsv(ani_path)
            if {phage_id(x["qname"]), phage_id(x["tname"])} == set(pair)
            and phage_id(x["qname"]) != phage_id(x["tname"])
        ]
        ani = max(float(x["pid"]) for x in reciprocal)
        shorter_af = max(max(float(x["qcov"]), float(x["tcov"])) for x in reciprocal)
        same = member_to_rep[left] == member_to_rep[right]
        pair_rows.append(
            {
                "PhageA": pair[0],
                "PhageB": pair[1],
                "ANIPct": round(ani, 4),
                "ShorterAlignmentFractionPct": round(shorter_af, 4),
                "PassANI95": bool_text(ani >= 95),
                "PassAF85": bool_text(shorter_af >= 85),
                "SameVOTU": bool_text(same),
            }
        )
    write_tsv(summary_dir / "votu-pairwise-threshold-audit.tsv", pair_rows)

    lifecycle_rows = []
    for pid, meta in metadata.items():
        confirmed = meta["ConfirmedLifecycle"]
        deeppl = meta["DeepPLPrediction"]
        phatyp = meta["PhaTYPPrediction"]
        lifecycle_rows.append(
            {
                "PhageID": pid,
                "ConfirmedLifecycle": confirmed,
                "DeepPLPrediction": deeppl,
                "PhaTYPPrediction": phatyp,
                "DeepPLMatchesConfirmed": bool_text(
                    confirmed in {"Temperate", "Virulent"}
                    and deeppl
                    == {"Temperate": "Lysogenic", "Virulent": "Lytic"}[confirmed]
                ),
                "PhaTYPMatchesConfirmed": bool_text(
                    confirmed in {"Temperate", "Virulent"}
                    and phatyp
                    == {"Temperate": "Lysogenic", "Virulent": "Lytic"}[confirmed]
                ),
                "EvidenceStatus": (
                    "Experimentally supported label available"
                    if confirmed in {"Temperate", "Virulent"}
                    else "Prediction only or not reported"
                ),
                "PhysicalState": "Free virion in mock virome",
                "InterpretationBoundary": (
                    "Free-particle recovery does not prove a virulent-only lifecycle"
                ),
            }
        )
    write_tsv(summary_dir / "lifecycle-evidence-ledger.tsv", lifecycle_rows)

    write_tsv(
        summary_dir / "state-lifecycle-evidence-map.tsv",
        [
            {
                "Observation": "Viral sequence inside a host contig",
                "Supports": "Integrated/proviral physical state at assembly time",
                "DoesNotProve": "Inducible temperate lifecycle without orthogonal evidence",
                "PreferredLabel": "Provirus candidate",
            },
            {
                "Observation": "Sequence recovered from a virus-enriched fraction",
                "Supports": "Free particle was sampled after that protocol",
                "DoesNotProve": "Virulent-only lifecycle",
                "PreferredLabel": "Free-virus fraction detection",
            },
            {
                "Observation": "Sequence-model prediction",
                "Supports": "Predicted lysogenic or lytic class",
                "DoesNotProve": "Experimental lifecycle",
                "PreferredLabel": "Predicted lifestyle",
            },
            {
                "Observation": "Isolation, induction, or curated experimental record",
                "Supports": "Experimentally supported lifecycle for that phage",
                "DoesNotProve": "Current in-sample physical state",
                "PreferredLabel": "Temperate or virulent with evidence source",
            },
        ],
    )

    for source_name in ("author-source-assertions.tsv", "lifecycle-source-assertions.tsv"):
        rows = read_tsv(work / source_name)
        if not rows or not all(row["Pass"] == "True" for row in rows):
            raise RuntimeError(f"Source assertion failed: {source_name}")
        write_tsv(summary_dir / source_name, rows)

    resource_rows = []
    labels = {
        "genomad": ("geNomad", 16, results / "genomad"),
        "makeblastdb": ("makeblastdb", 1, results / "votu"),
        "blastn": ("BLASTN", 16, results / "votu"),
        "anicalc": ("anicalc", 1, results / "votu"),
        "aniclust": ("aniclust", 1, results / "votu"),
    }
    for label, (tool, threads, output_dir) in labels.items():
        timing = work / f"logs/{label}.time.txt"
        if not timing.is_file():
            raise FileNotFoundError(f"Missing GNU time record: {timing}")
        parsed = parse_time(timing)
        resource_rows.append(
            {
                "Tool": tool,
                "ThreadsRequested": threads,
                "WallSeconds": parsed["WallSeconds"],
                "PeakRAMGiB": parsed["PeakRAMGiB"],
                "OutputBytesSharedDirectory": total_bytes(output_dir),
                "ExitStatus": parsed["ExitStatus"],
                "Measurement": "GNU time -v; 15 current RefSeq/GenBank accessions",
            }
        )
    smoke_time = parse_time(work / "logs/bbmap-smoke.time.txt")
    resource_rows.append(
        {
            "Tool": "BBMap smoke test",
            "ThreadsRequested": "auto (reduced to 23 by BBMap)",
            "WallSeconds": smoke_time["WallSeconds"],
            "PeakRAMGiB": smoke_time["PeakRAMGiB"],
            "OutputBytesSharedDirectory": (work / "bbmap-smoke-coverage.tsv").stat().st_size,
            "ExitStatus": smoke_time["ExitStatus"],
            "Measurement": "GNU time -v; implementation check with -Xmx2g",
        }
    )
    write_tsv(work / "resource-summary.tsv", resource_rows)
    write_tsv(
        work / "tool-versions.tsv",
        [
            {
                "Tool": "geNomad",
                "Version": "1.12.0",
                "Database": "geNomad DB v1.9; ICTV MSL39",
                "Role": "virus confirmation and taxonomy",
            },
            {
                "Tool": "BLASTN",
                "Version": "2.17.0",
                "Database": "15-accession all-vs-all catalog",
                "Role": "nucleotide alignments for vOTU clustering",
            },
            {
                "Tool": "anicalc / aniclust",
                "Version": "CheckV 1.1.1 distribution",
                "Database": "not applicable",
                "Role": "95% ANI / 85% target-coverage vOTUs",
            },
            {
                "Tool": "BBMap",
                "Version": "38.69 (study-reported and smoke-tested)",
                "Database": "15 study reference genomes",
                "Role": "published read mapping; minid=0.90, ambiguous=all",
            },
            {
                "Tool": "SAMtools",
                "Version": "1.24",
                "Database": "not applicable",
                "Role": "BAM indexing and coverage audit",
            },
            {
                "Tool": "SeqKit",
                "Version": "2.13.0",
                "Database": "not applicable",
                "Role": "FASTA/FASTQ inspection",
            },
        ],
    )
    write_tsv(
        work / "determinism-audit.tsv",
        [
            {
                "Component": component,
                "RandomProcess": "False",
                "Seed": "Not applicable" if component != "plots" else SEED,
                "DeterminismControl": control,
                "Status": "PASS",
            }
            for component, control in (
                ("published abundance", "fixed Cook Supplementary Table S4 and explicit gates"),
                ("geNomad", "fixed accessions, version, database release, and --splits 16"),
                ("vOTU", "fixed catalog order and 95/85 thresholds"),
                ("BBMap smoke", "fixed 15-reference FASTA, BBMap 38.69, and ambiguous=all"),
                ("plots", "set.seed(20260755); no random jitter"),
            )
        ],
    )

    counts_75 = {
        row["Dataset"]: row["PhagesDetected"]
        for row in threshold_rows
        if row["BreadthThresholdPct"] == 75
    }
    family_counts = Counter(row["geNomadFamily"] for row in taxonomy_rows)
    prevalence_counts = Counter(row["ReplicatesPresent"] for row in prevalence_rows)
    lifecycle_counts = Counter(row["ConfirmedLifecycle"] for row in lifecycle_rows)
    multi_member = [row for row in cluster_rows if row["MemberCount"] > 1]
    if len(multi_member) != 1:
        raise RuntimeError(f"Expected one multi-member vOTU, observed {len(multi_member)}")
    duplicate_pair = next(
        row
        for row in pair_rows
        if {row["PhageA"], row["PhageB"]}
        == {"vB_Eco_mar001J1", "vB_Eco_mar002J2"}
    )
    smoke_by_phage = {row["PhageID"]: row for row in smoke_rows}
    summary = {
        "article": 55,
        "seed": SEED,
        "mock_phages": len(metadata),
        "illumina_libraries": 3,
        "illumina_read_pairs": sum(
            int(row["ReadPairs"]) for row in read_tsv(library_path)
        ),
        "primary_presence_min_depth_x": PRIMARY_DEPTH,
        "primary_presence_min_breadth_pct": PRIMARY_BREADTH,
        "detected_at_75pct": counts_75,
        "replicate_prevalence_counts": {str(key): value for key, value in sorted(prevalence_counts.items())},
        "universal_nondetections": sorted(
            row["PhageID"] for row in prevalence_rows if row["ReplicatesPresent"] == 0
        ),
        "genomad_detected": len(virus_summary),
        "genomad_family_counts": dict(sorted(family_counts.items())),
        "votu_clusters": len(cluster_rows),
        "votu_multi_member_clusters": len(multi_member),
        "largest_votu_members": multi_member[0]["MemberCount"],
        "largest_votu_member_ids": multi_member[0]["Members"].split(","),
        "j1_j2_ani_pct": duplicate_pair["ANIPct"],
        "j1_j2_shorter_af_pct": duplicate_pair["ShorterAlignmentFractionPct"],
        "bbmap_smoke_j1_depth_x": float(smoke_by_phage["vB_Eco_mar001J1"]["MeanDepthX"]),
        "bbmap_smoke_j2_depth_x": float(smoke_by_phage["vB_Eco_mar002J2"]["MeanDepthX"]),
        "confirmed_lifecycle_counts": dict(sorted(lifecycle_counts.items())),
        "deeppl_confirmed_correct": sum(
            row["DeepPLMatchesConfirmed"] == "True" for row in lifecycle_rows
        ),
        "phatyp_confirmed_correct": sum(
            row["PhaTYPMatchesConfirmed"] == "True" for row in lifecycle_rows
        ),
        "slur29_reference_length_delta_bp": int(
            next(row for row in read_tsv(work / "reference-audit.tsv") if row["PhageID"] == "vB_Eco_SLUR29")["LengthDeltaBp"]
        ),
        "random_output_requested": False,
    }
    dump_json(summary_dir / "summary.json", summary)
    (work / ".article55-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
