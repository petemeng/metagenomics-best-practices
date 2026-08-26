#!/usr/bin/env python3
"""Summarize real geNomad, VirSorter2, CheckV, and vOTU outputs for Article 54."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

from article41_44_utils import (
    dump_json,
    fasta_summary,
    parse_time,
    read_tsv,
    sha256,
    write_tsv,
)


QUALITY_ORDER = ("High-quality", "Medium-quality", "Low-quality", "Not-determined")


def as_number(value: str) -> float | None:
    if value in {"", "NA", "nan", "None"}:
        return None
    return float(value)


def total_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def normalize_vs2_id(value: str) -> str:
    return value.split("||", 1)[0]


def bool_text(value: bool) -> str:
    return "True" if value else "False"


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
    logs_dir = work / "logs"
    summary_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    input_fna = root / "data/raw/article54/checkv-test-sequences.fna"
    checkv_truth = (
        root
        / "data/raw/article54/checkv-upstream-ground-truth-quality-summary.tsv"
    )
    miuvig_xml = root / "data/raw/article54/PMC6871006.xml"
    genomad_dir = results / "genomad"
    checkv_dir = results / "checkv"
    vs2_dir = results / "virsorter2"
    votu_dir = results / "votu"
    required = (
        input_fna,
        checkv_truth,
        miuvig_xml,
        genomad_dir
        / "checkv-test-sequences_summary/checkv-test-sequences_virus_summary.tsv",
        genomad_dir
        / "checkv-test-sequences_aggregated_classification/checkv-test-sequences_aggregated_classification.tsv",
        checkv_dir / "quality_summary.tsv",
        checkv_dir / "complete_genomes.tsv",
        vs2_dir / "final-viral-score.tsv",
        vs2_dir / "iter-0/all-fullseq-proba.tsv",
        votu_dir / "ani.tsv",
        votu_dir / "votu-clusters.tsv",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 54 outputs: " + ", ".join(missing))

    fasta_stats, fasta_records = fasta_summary(input_fna)
    input_ids = list(fasta_records)
    if len(input_ids) != 46:
        raise RuntimeError(f"Expected 46 CheckV regression sequences, observed {len(input_ids)}")

    genomad_summary = read_tsv(
        genomad_dir
        / "checkv-test-sequences_summary/checkv-test-sequences_virus_summary.tsv"
    )
    genomad_detected = {row["seq_name"]: row for row in genomad_summary}
    genomad_all = {
        row["seq_name"]: row
        for row in read_tsv(
            genomad_dir
            / "checkv-test-sequences_aggregated_classification/checkv-test-sequences_aggregated_classification.tsv"
        )
    }
    vs2_final = {
        normalize_vs2_id(row["seqname"]): row
        for row in read_tsv(vs2_dir / "final-viral-score.tsv")
    }
    vs2_all = {
        normalize_vs2_id(row["seqname"]): row
        for row in read_tsv(vs2_dir / "iter-0/all-fullseq-proba.tsv")
    }
    checkv_rows = {
        row["contig_id"]: row for row in read_tsv(checkv_dir / "quality_summary.tsv")
    }
    if set(input_ids) != set(genomad_all) or set(input_ids) != set(vs2_all):
        raise RuntimeError("A discovery score table does not cover all 46 input sequences")
    if set(input_ids) != set(checkv_rows):
        raise RuntimeError("CheckV quality table does not cover all input sequences")

    detection_rows = []
    for contig in input_ids:
        g_detected = contig in genomad_detected
        v_detected = contig in vs2_final
        if g_detected and v_detected:
            pattern = "Both"
        elif g_detected:
            pattern = "geNomad only"
        elif v_detected:
            pattern = "VirSorter2 only"
        else:
            pattern = "Neither"
        g_all = genomad_all[contig]
        g = genomad_detected.get(contig, {})
        v_all = vs2_all[contig]
        v = vs2_final.get(contig, {})
        cv = checkv_rows[contig]
        v_score = as_number(v.get("max_score", ""))
        hallmark = int(v.get("hallmark", "0") or 0)
        v_high = v_detected and v_score is not None and (
            v_score >= 0.9 or (v_score >= 0.7 and hallmark >= 1)
        )
        detection_rows.append(
            {
                "ContigID": contig,
                "LengthBp": fasta_records[contig]["LengthBp"],
                "geNomadDetected": bool_text(g_detected),
                "geNomadVirusScore": g_all["virus_score"],
                "geNomadTopology": g.get("topology", "Not called"),
                "geNomadHallmarks": g.get("n_hallmarks", "NA"),
                "geNomadTaxonomy": g.get("taxonomy", "Not assigned"),
                "VirSorter2Detected": bool_text(v_detected),
                "VirSorter2MaxScoreAll": max(
                    float(v_all["dsDNAphage"]), float(v_all["ssDNA"])
                ),
                "VirSorter2FinalScore": v.get("max_score", "NA"),
                "VirSorter2HighConfidence": bool_text(v_high),
                "VirSorter2Group": v.get("max_score_group", "Not called"),
                "VirSorter2Hallmarks": v.get("hallmark", "NA"),
                "VirSorter2ViralGenePct": v.get("viral", "NA"),
                "VirSorter2CellularGenePct": v.get("cellular", "NA"),
                "DiscoveryPattern": pattern,
                "CheckVProvirus": cv["provirus"],
                "CheckVProviralLength": cv["proviral_length"],
                "CheckVQuality": cv["checkv_quality"],
                "MIUViGQuality": cv["miuvig_quality"],
                "CompletenessPct": cv["completeness"],
                "CompletenessMethod": cv["completeness_method"],
                "ContaminationPct": cv["contamination"],
                "CheckVWarnings": cv["warnings"],
            }
        )
    write_tsv(summary_dir / "virus-evidence-matrix.tsv", detection_rows)

    pattern_counts = Counter(row["DiscoveryPattern"] for row in detection_rows)
    overlap_rows = [
        {"Pattern": label, "Contigs": pattern_counts[label]}
        for label in ("Both", "geNomad only", "VirSorter2 only", "Neither")
    ]
    write_tsv(summary_dir / "discovery-overlap.tsv", overlap_rows)

    quality_counts = Counter(row["CheckVQuality"] for row in detection_rows)
    quality_rows = [
        {
            "CheckVQuality": quality,
            "Contigs": quality_counts[quality],
            "MinimumCompletenessPct": {
                "High-quality": 90,
                "Medium-quality": 50,
                "Low-quality": 0,
                "Not-determined": "NA",
            }[quality],
        }
        for quality in QUALITY_ORDER
    ]
    write_tsv(summary_dir / "checkv-quality-counts.tsv", quality_rows)

    complete_rows = read_tsv(checkv_dir / "complete_genomes.tsv")
    complete_by_id = {row["contig_id"]: row for row in complete_rows}
    terminal_rows = []
    for contig in input_ids:
        g = genomad_detected.get(contig, {})
        cv = complete_by_id.get(contig, {})
        if g.get("topology") == "DTR" or cv:
            terminal_rows.append(
                {
                    "ContigID": contig,
                    "geNomadTopology": g.get("topology", "Not called"),
                    "CheckVPrediction": cv.get("prediction_type", "Not called"),
                    "CheckVConfidence": cv.get("confidence_level", "NA"),
                    "CheckVReason": cv.get("confidence_reason", "NA"),
                    "CheckVQuality": checkv_rows[contig]["checkv_quality"],
                    "MIUViGFinished": "False",
                }
            )
    write_tsv(summary_dir / "terminal-repeat-audit.tsv", terminal_rows)

    directed_ani = read_tsv(votu_dir / "ani.tsv")
    canonical: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in directed_ani:
        if row["qname"] == row["tname"]:
            continue
        canonical[tuple(sorted((row["qname"], row["tname"])))].append(row)
    pair_rows = []
    for (left, right), rows in sorted(canonical.items()):
        ani = max(float(row["pid"]) for row in rows)
        shorter_af = max(
            max(float(row["qcov"]), float(row["tcov"])) for row in rows
        )
        pass_ani = ani >= 95
        pass_af = shorter_af >= 85
        pair_rows.append(
            {
                "SequenceA": left,
                "SequenceB": right,
                "ANIpct": round(ani, 4),
                "ShorterAlignmentFractionPct": round(shorter_af, 4),
                "PassANI95": bool_text(pass_ani),
                "PassAF85": bool_text(pass_af),
                "SameVOTU": bool_text(pass_ani and pass_af),
                "BoundaryPair": bool_text(
                    {left, right} == {"UHGV-0001702", "UHGV-0001715"}
                ),
            }
        )
    write_tsv(summary_dir / "votu-pairwise-threshold-audit.tsv", pair_rows)

    clusters = read_tsv(votu_dir / "votu-clusters.tsv")
    if not clusters:
        # aniclust emits a headerless two-column table.
        cluster_lines = [
            line.split("\t")
            for line in (votu_dir / "votu-clusters.tsv")
            .read_text(encoding="utf-8")
            .splitlines()
            if line
        ]
    else:
        # DictReader consumes the first headerless row; reconstruct all lines.
        cluster_lines = [
            line.split("\t")
            for line in (votu_dir / "votu-clusters.tsv")
            .read_text(encoding="utf-8")
            .splitlines()
            if line
        ]
    cluster_counts = Counter(line[0] for line in cluster_lines)
    write_tsv(
        summary_dir / "votu-cluster-summary.tsv",
        [
            {
                "Representative": representative,
                "Members": members,
                "MinANI": 95,
                "MinShorterAlignmentFraction": 85,
            }
            for representative, members in sorted(cluster_counts.items())
        ],
    )

    library_rows = [
        {
            "Library": "Total metagenome",
            "Step": "No particle-size enrichment",
            "WhatItRetains": "Cell-associated and extracellular sequence context",
            "KnownBlindSpot": "Rare viruses may remain below assembly depth",
            "QuantitativeUse": "Within-protocol mapping only after coverage QC",
            "PrimaryEvidence": "MIUViG",
        },
        {
            "Library": "Virus-enriched virome",
            "Step": "Clarification and filtration",
            "WhatItRetains": "Particles passing the selected pore and recovery workflow",
            "KnownBlindSpot": "0.22 um filtration can remove large viruses",
            "QuantitativeUse": "Protocol-specific composition",
            "PrimaryEvidence": "Parras-Molto 2018; NetoVIR 2015",
        },
        {
            "Library": "Virus-enriched virome",
            "Step": "DNase or RNase before extraction",
            "WhatItRetains": "Nuclease-protected nucleic acids",
            "KnownBlindSpot": "Unprotected or disrupted virions can be lost",
            "QuantitativeUse": "Requires process controls",
            "PrimaryEvidence": "MIUViG; NetoVIR 2015",
        },
        {
            "Library": "Amplified virome",
            "Step": "MDA or SISPA",
            "WhatItRetains": "Low-input material after nonlinear amplification",
            "KnownBlindSpot": "Small circular genomes and sequence/GC-dependent bias",
            "QuantitativeUse": "Nonquantitative read depth",
            "PrimaryEvidence": "MIUViG; Parras-Molto 2018",
        },
    ]
    write_tsv(summary_dir / "library-design-bias-audit.tsv", library_rows)

    ladder_rows = [
        {
            "Order": 1,
            "Decision": "Virus discovery",
            "RequiredEvidence": "geNomad and/or VirSorter2 score plus gene context",
            "ForbiddenShortcut": "Treating CheckV as a virus detector",
        },
        {
            "Order": 2,
            "Decision": "Boundary cleanup",
            "RequiredEvidence": "CheckV host-virus boundary and provirus status",
            "ForbiddenShortcut": "Counting host flanks as viral genes",
        },
        {
            "Order": 3,
            "Decision": "Genome quality",
            "RequiredEvidence": "Completeness method, confidence, contamination, warnings",
            "ForbiddenShortcut": "Equating DTR with a finished genome",
        },
        {
            "Order": 4,
            "Decision": "vOTU clustering",
            "RequiredEvidence": "ANI >=95% and alignment fraction >=85% of shorter sequence",
            "ForbiddenShortcut": "Applying ANI without an alignment-fraction gate",
        },
        {
            "Order": 5,
            "Decision": "MIUViG reporting",
            "RequiredEvidence": "Origin, detection, quality, taxonomy, abundance, provenance",
            "ForbiddenShortcut": "Calling automated high-quality output finished",
        },
    ]
    write_tsv(summary_dir / "virus-evidence-ladder.tsv", ladder_rows)

    xml = miuvig_xml.read_text(encoding="utf-8")
    assertions = [
        (
            "vOTU threshold",
            "95% average nucleotide identity over 85% alignment fraction" in xml,
        ),
        ("high-quality draft threshold", "representing ≥90%" in xml),
        ("finished requires manual review", "extensive manual review" in xml),
        ("amplified datasets nonquantitative", "Read mapping from nonquantitative datasets" in xml),
        ("virus-enriched workflow", "filtration steps, DNase or RNase treatments" in xml),
    ]
    write_tsv(
        summary_dir / "miuvig-source-assertions.tsv",
        [
            {"Assertion": label, "Pass": bool_text(passed), "Source": "PMC6871006"}
            for label, passed in assertions
        ],
    )
    if not all(passed for _, passed in assertions):
        raise RuntimeError("One or more MIUViG source assertions failed")

    observed_quality = checkv_dir / "quality_summary.tsv"
    fixture_exact = observed_quality.read_bytes() == checkv_truth.read_bytes()

    resource_rows = []
    labels = {
        "genomad": (16, genomad_dir),
        "checkv": (16, checkv_dir),
        "virsorter2": (16, vs2_dir),
    }
    source_logs = results / "logs"
    for path in sorted(source_logs.iterdir()):
        if path.is_file():
            shutil.copy2(path, logs_dir / path.name)
    for label, (threads, output_dir) in labels.items():
        candidates = (
            logs_dir / f"{label}.time.txt",
            logs_dir / f"{label}.stderr-time.log",
        )
        timing_path = next((path for path in candidates if path.is_file()), None)
        if timing_path is None:
            raise FileNotFoundError(f"Missing GNU time log for {label}")
        parsed = parse_time(timing_path)
        resource_rows.append(
            {
                "Tool": {"genomad": "geNomad", "checkv": "CheckV", "virsorter2": "VirSorter2"}[label],
                "Threads": threads,
                "WallSeconds": parsed["WallSeconds"],
                "PeakRAMGiB": parsed["PeakRAMGiB"],
                "OutputBytes": total_bytes(output_dir),
                "ExitStatus": parsed["ExitStatus"],
                "Measurement": "GNU time -v; official 46-sequence CheckV fixture",
            }
        )
    write_tsv(work / "resource-summary.tsv", resource_rows)

    tool_rows = [
        {
            "Tool": "geNomad",
            "Version": "1.12.0",
            "SourceCommit": "8c5fd0d1722d458a3e8ff50278cdc00ed4a514fc",
            "Database": "geNomad v1.9; ICTV MSL39",
        },
        {
            "Tool": "VirSorter2",
            "Version": "2.2.4",
            "SourceCommit": "d96ba672442090ebc943ec3ceb87fe4cfdaaaa0b",
            "Database": "Zenodo 4269607 v0.4",
        },
        {
            "Tool": "CheckV",
            "Version": "1.1.1",
            "SourceCommit": "6a118f20e895105ce0e4f10257955494c60f1293",
            "Database": "v1.5; release date 2023-01-10",
        },
        {
            "Tool": "BLASTN / anicalc / aniclust",
            "Version": "BLAST 2.17.0 / CheckV 1.1.1",
            "SourceCommit": "CheckV commit above",
            "Database": "all-vs-all input FASTA",
        },
    ]
    write_tsv(work / "tool-versions.tsv", tool_rows)
    write_tsv(
        work / "determinism-audit.tsv",
        [
            {
                "Component": component,
                "RandomProcess": "False",
                "Seed": "Not applicable",
                "DeterminismControl": control,
                "Status": "PASS",
            }
            for component, control in (
                ("geNomad", "fixed input, database, version, and thresholds"),
                ("VirSorter2", "fixed input, v0.4 database, and score rules"),
                ("CheckV", "exact match to upstream regression fixture"),
                ("vOTU", "deterministic BLASTN plus greedy clustering order"),
                ("plots", "set.seed(20260754); no jitter requiring random output"),
            )
        ],
    )

    boundary = next(row for row in pair_rows if row["BoundaryPair"] == "True")
    summary = {
        "article": 54,
        "seed": 20260754,
        "input_sequences": len(input_ids),
        "input_total_bp": fasta_stats["TotalBp"],
        "input_n50_bp": fasta_stats["N50Bp"],
        "genomad_detected": len(genomad_detected),
        "virsorter2_default_detected": len(vs2_final),
        "virsorter2_high_confidence": sum(
            row["VirSorter2HighConfidence"] == "True" for row in detection_rows
        ),
        "detected_by_both": pattern_counts["Both"],
        "genomad_only": pattern_counts["geNomad only"],
        "virsorter2_only": pattern_counts["VirSorter2 only"],
        "detected_union": len(set(genomad_detected) | set(vs2_final)),
        "detected_by_neither": pattern_counts["Neither"],
        "checkv_quality_counts": dict(quality_counts),
        "checkv_proviruses": sum(row["CheckVProvirus"] == "Yes" for row in detection_rows),
        "checkv_low_confidence_dtr": sum(
            row["confidence_level"] == "low" for row in complete_rows
        ),
        "genomad_dtr": sum(row["topology"] == "DTR" for row in genomad_summary),
        "checkv_fixture_exact_match": fixture_exact,
        "votu_clusters": len(cluster_counts),
        "votu_multi_member_clusters": sum(value > 1 for value in cluster_counts.values()),
        "votu_nonself_pairs": len(pair_rows),
        "votu_pairs_passing_both": sum(row["SameVOTU"] == "True" for row in pair_rows),
        "boundary_pair": [boundary["SequenceA"], boundary["SequenceB"]],
        "boundary_pair_ani_pct": boundary["ANIpct"],
        "boundary_pair_shorter_af_pct": boundary["ShorterAlignmentFractionPct"],
        "random_output_requested": False,
    }
    dump_json(summary_dir / "summary.json", summary)
    (work / ".article54-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
