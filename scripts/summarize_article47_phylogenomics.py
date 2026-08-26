#!/usr/bin/env python3
"""Summarize marker recovery, trees and conservative novelty decisions."""

from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path

from Bio import AlignIO, Phylo

from article41_44_utils import dump_json, parse_time, read_tsv, sha256, write_tsv


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1", "yes"}


def model_from_report(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    patterns = (
        r"Best-fit model according to BIC:\s*(\S+)",
        r"Best-fit model:\s*(\S+)",
    )
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return match.group(1)
    return "MODEL_NOT_PARSED"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article47-run-complete").is_file():
        raise FileNotFoundError("Run run_article47_phylogenomics.py first")

    genome_ledger = {row["SGB"]: row for row in read_tsv(work / "genome-ledger.tsv")}
    alignments = read_tsv(work / "alignment-paths.tsv")
    alignment_rows: list[dict[str, object]] = []
    domain_rows: list[dict[str, object]] = []
    tree_rows: list[dict[str, object]] = []
    query_sgbs_in_trees: set[str] = set()
    for row in alignments:
        domain = row["Domain"]
        alignment_path = Path(row["Alignment"])
        tree_path = Path(row["TreeFile"])
        report_path = Path(row["IQTreeReport"])
        if not alignment_path.is_file() or not tree_path.is_file() or not report_path.is_file():
            raise FileNotFoundError(f"Missing alignment/tree/report for {domain}")
        alignment = AlignIO.read(alignment_path, "fasta")
        tree = Phylo.read(tree_path, "newick")
        alignment_tips = {record.id for record in alignment}
        tree_tips = {tip.name for tip in tree.get_terminals()}
        if alignment_tips != tree_tips:
            raise ValueError(f"Alignment/tree tip mismatch for {domain}")
        query_ids = {
            tip.removesuffix("_QUERY_MAG")
            for tip in tree_tips
            if tip.startswith("SGB_")
        }
        query_sgbs_in_trees.update(query_ids)
        query_count = len(query_ids)
        reference_count = sum(tip.startswith("REF_") for tip in tree_tips)
        expected_query_ids = {
            sgb for sgb, genome in genome_ledger.items() if genome["Domain"] == domain
        }
        if not query_ids or not query_ids.issubset(expected_query_ids):
            raise ValueError(f"Unexpected query tips for {domain}: {sorted(query_ids)}")
        if reference_count < 1:
            raise ValueError(f"No reference tips retained for {domain}")
        for record in alignment:
            sequence = str(record.seq)
            observed = sum(base not in {"-", ".", "?", "X"} for base in sequence.upper())
            kind = "Query MAG" if record.id.startswith("SGB_") else "Reference genome"
            alignment_rows.append({
                "Domain": domain,
                "Tip": record.id,
                "Type": kind,
                "AlignmentSites": len(sequence),
                "ObservedAminoAcids": observed,
                "OccupancyPct": round(100 * observed / len(sequence), 6),
            })
        support_pairs: list[tuple[float, float]] = []
        for clade in tree.get_nonterminals():
            match = re.fullmatch(
                r"([0-9]+(?:\.[0-9]+)?)/([0-9]+(?:\.[0-9]+)?)",
                clade.name or "",
            )
            if match:
                support_pairs.append((float(match.group(1)), float(match.group(2))))
        if not support_pairs:
            raise ValueError(f"No SH-aLRT/UFBoot branch supports parsed for {domain}")
        sh_alrt = [pair[0] for pair in support_pairs]
        ufboot = [pair[1] for pair in support_pairs]
        domain_rows.append({
            "Domain": domain,
            "InputQueryMAGs": len(expected_query_ids),
            "QueryMAGs": query_count,
            "ExcludedQueryMAGs": len(expected_query_ids - query_ids),
            "ExcludedQueryIDs": ";".join(sorted(expected_query_ids - query_ids)),
            "UniqueReferenceTips": reference_count,
            "TreeTips": len(tree_tips),
            "AlignmentSites": alignment.get_alignment_length(),
            "MedianTipOccupancyPct": statistics.median(
                100 * sum(base not in {"-", ".", "?", "X"} for base in str(record.seq).upper()) / len(record.seq)
                for record in alignment
            ),
            "BestFitModel": model_from_report(report_path),
            "Seed": int(row["Seed"]),
            "AlignmentSHA256": sha256(alignment_path),
            "TreeSHA256": sha256(tree_path),
            "SupportedInternalBranches": len(support_pairs),
            "MedianSHaLRT": round(statistics.median(sh_alrt), 6),
            "PctSHaLRTGe80": round(100 * sum(value >= 80 for value in sh_alrt) / len(sh_alrt), 6),
            "MedianUFBoot": round(statistics.median(ufboot), 6),
            "PctUFBootGe95": round(100 * sum(value >= 95 for value in ufboot) / len(ufboot), 6),
        })
        for tip in sorted(tree_tips):
            tree_rows.append({
                "Domain": domain,
                "Tip": tip,
                "StableID": tip.removesuffix("_QUERY_MAG").removesuffix("_REFERENCE"),
                "Type": "Query MAG" if tip.startswith("SGB_") else "Reference genome",
                "InAlignment": tip in alignment_tips,
            })

    gtotree_rows: list[dict[str, object]] = []
    for domain in ("Bacteria", "Archaea"):
        source = work / f"gtotree/{domain.lower()}/Genomes_summary_info.tsv"
        for row in read_tsv(source):
            assembly_id = row["assembly_id"]
            in_tree = as_bool(row["in_final_tree"])
            if assembly_id.startswith("SGB_") and in_tree != (assembly_id in query_sgbs_in_trees):
                raise ValueError(f"GToTree inclusion/tree mismatch for {assembly_id}")
            gtotree_rows.append({
                "Domain": domain,
                "AssemblyID": assembly_id,
                "Label": row["label"],
                "Type": "Query MAG" if assembly_id.startswith("SGB_") else "Reference genome",
                "NumSCGHits": row["num_SCG_hits"],
                "UniqueSCGHits": row["uniq_SCG_hits"],
                "SCGCompletenessPct": row["perc_comp"],
                "SCGRedundancyPct": row["perc_redund"],
                "SCGsAfterLengthFilter": row["num_SCG_hits_after_len_filt"],
                "InFinalTree": in_tree,
            })

    pretree = read_tsv(work / "novelty-pretree-audit.tsv")
    final_novelty: list[dict[str, object]] = []
    for row in pretree:
        candidate = as_bool(row["PreTreeNovelSpeciesCandidate"])
        included = row["SGB"] in query_sgbs_in_trees
        final_novelty.append({
            **row,
            "PhylogenomicTreeIncluded": included,
            "TreeEligibilityNote": (
                "Retained by GToTree SCG filters"
                if included else "Excluded by GToTree SCG filters; ANI-based novelty audit retained"
            ),
            "FinalNovelSpeciesStatus": (
                "CANDIDATE_REQUIRES_BROADER_REFERENCE_AND_NOMENCLATURAL_REVIEW"
                if candidate else "NOVEL_SPECIES_NOT_SUPPORTED"
            ),
            "CandidatusNameProposed": False,
            "NovelSpeciesClaim": False,
        })

    references = read_tsv(work / "reference-request-ledger.tsv")
    reference_inventory: list[dict[str, object]] = []
    for ref in references:
        accession = ref["ReferenceAccession"]
        matches = [
            path for path in (work / "gtotree").rglob("*")
            if path.is_file() and accession in path.name
        ]
        reference_inventory.append({
            **ref,
            "DownloadedArtifactCount": len(matches),
            "DownloadedArtifacts": ";".join(path.relative_to(work).as_posix() for path in matches),
            "ArtifactSHA256": ";".join(sha256(path) for path in matches),
        })

    resources = [
        parse_time(work / f"logs/{label}.time.txt")
        for label in ("gtotree-bacteria", "iqtree-bacteria", "gtotree-archaea", "iqtree-archaea")
    ]
    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(summary_dir / "alignment-occupancy.tsv", alignment_rows)
    write_tsv(summary_dir / "domain-tree-summary.tsv", domain_rows)
    write_tsv(summary_dir / "tree-tip-ledger.tsv", tree_rows)
    write_tsv(summary_dir / "gtotree-genome-audit.tsv", gtotree_rows)
    write_tsv(summary_dir / "novelty-decision-audit.tsv", final_novelty)
    write_tsv(summary_dir / "reference-download-inventory.tsv", reference_inventory)
    write_tsv(summary_dir / "resource-summary.tsv", resources)
    result = {
        "article": 47,
        "input_sgbs": 24,
        "tree_included_sgbs": len(query_sgbs_in_trees),
        "tree_excluded_sgbs": len(genome_ledger) - len(query_sgbs_in_trees),
        "tree_excluded_ids": sorted(set(genome_ledger) - query_sgbs_in_trees),
        "domains": {
            row["Domain"]: {
                "input_query_mags": row["InputQueryMAGs"],
                "query_mags": row["QueryMAGs"],
                "excluded_query_mags": row["ExcludedQueryMAGs"],
                "reference_tips": row["UniqueReferenceTips"],
                "alignment_sites": row["AlignmentSites"],
                "best_fit_model": row["BestFitModel"],
                "median_ufboot": row["MedianUFBoot"],
                "pct_ufboot_ge95": row["PctUFBootGe95"],
            }
            for row in domain_rows
        },
        "species_fields_assigned": sum(as_bool(row["SpeciesAssigned"]) for row in final_novelty),
        "species_fields_unresolved": sum(not as_bool(row["SpeciesAssigned"]) for row in final_novelty),
        "pretree_novelty_candidates": sum(as_bool(row["PreTreeNovelSpeciesCandidate"]) for row in final_novelty),
        "novel_species_claims": 0,
        "candidatus_names_proposed": 0,
        "truth_used_for_novelty_or_reference_selection": False,
    }
    dump_json(summary_dir / "run-summary.json", result)
    (work / ".article47-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
