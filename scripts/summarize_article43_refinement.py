#!/usr/bin/env python3
"""Summarize Article 43 refinement, provenance, and post-hoc mock truth audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

from article41_44_utils import dump_json, fasta_records, fasta_summary, parse_time, read_tsv, sha256, write_tsv
from summarize_article42_binning import n50, parse_coords, truth_assignments, union_length


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--article42-frozen", type=Path, required=True)
    return parser.parse_args()


def read_qc(checkm_path: Path, gunc_path: Path):
    checkm, gunc = {}, {}
    for row in csv.DictReader(checkm_path.open(encoding="utf-8"), delimiter="\t"):
        name = row.get("Name") or row.get("name") or row.get("Genome")
        if name:
            checkm[name.removesuffix(".fna")] = row
    for row in csv.DictReader(gunc_path.open(encoding="utf-8"), delimiter="\t"):
        name = row.get("genome") or row.get("Genome")
        if name:
            gunc[name.removesuffix(".fna")] = row
    return checkm, gunc


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes", "pass"}


def main() -> int:
    args = parse_args()
    root, work, frozen42 = args.project_root.resolve(), args.work_dir.resolve(), args.article42_frozen.resolve()
    if not (work / ".article43-run-complete").is_file():
        raise FileNotFoundError("Run run_article43_refinement.py first")
    summary = work / "summary"
    summary.mkdir(exist_ok=True)
    common = work / "inputs/megahit-coassembly.ge1500.fna"
    assembly_summary, contig_stats = fasta_summary(common)
    sequences = dict(fasta_records(common))
    if set(sequences) != set(contig_stats):
        raise RuntimeError("Article 43 coordinate set is internally inconsistent")

    truth_coords = root / "data/raw/article33/work/metaquast/MOCK1_MOCK2/combined_reference/contigs_reports/minimap_output/sr-megahit-co.coords"
    manifest_path = root / "data/small/33-assembly-qc-frozen/truth-manifest.tsv"
    manifest_rows = read_tsv(manifest_path)
    truth_manifest = {
        row["GenBankAssembly"]: row for row in manifest_rows if row["EvaluationSet"] == "MOCK1+MOCK2"
    }
    if len(truth_manifest) != 87:
        raise RuntimeError(f"Expected 87 co-assembly truth genomes, observed {len(truth_manifest)}")
    query_intervals, ref_intervals, identities, alignment_rows = parse_coords(truth_coords, set(contig_stats))
    truth_rows = truth_assignments(contig_stats, query_intervals, identities, truth_manifest)
    truth_by_contig = {row["Contig"]: row for row in truth_rows}
    frozen_truth = read_tsv(frozen42 / "truth-contig-assignment.tsv.gz")
    frozen_truth_core = {
        row["Contig"]: (row["BestTruthAccession"], row["BestAlignedQueryBp"], row["AmbiguousTruthAssignment"])
        for row in frozen_truth
    }
    current_truth_core = {
        row["Contig"]: (str(row["BestTruthAccession"]), str(row["BestAlignedQueryBp"]), str(row["AmbiguousTruthAssignment"]))
        for row in truth_rows
    }
    if frozen_truth_core != current_truth_core:
        raise RuntimeError("Article 42 frozen truth coordinate audit has drifted")
    write_tsv(summary / "truth-input-audit.tsv", [
        {"Role": "Article42-frozen-truth-assignment", "Path": str((frozen42 / "truth-contig-assignment.tsv.gz").relative_to(root)), "Bytes": (frozen42 / "truth-contig-assignment.tsv.gz").stat().st_size, "SHA256": sha256(frozen42 / "truth-contig-assignment.tsv.gz"), "Status": "PASS"},
        {"Role": "MetaQUAST-coassembly-coordinates", "Path": str(truth_coords.relative_to(root)), "Bytes": truth_coords.stat().st_size, "SHA256": sha256(truth_coords), "Status": "PASS"},
        {"Role": "Article33-truth-manifest", "Path": str(manifest_path.relative_to(root)), "Bytes": manifest_path.stat().st_size, "SHA256": sha256(manifest_path), "Status": "PASS"},
    ])

    maxcss_candidates = sorted((work / "qc/gunc").glob("GUNC.*.maxCSS_level.tsv"))
    if len(maxcss_candidates) != 1:
        raise RuntimeError(f"Expected one GUNC maxCSS table: {maxcss_candidates}")
    checkm, gunc = read_qc(work / "qc/checkm2/quality_report.tsv", maxcss_candidates[0])
    source_rows = read_tsv(work / "refined-bin-source.tsv")
    source_by_id = {row["RefinedID"]: row for row in source_rows}
    flat = work / "qc/all-refined"
    refined_fastas = sorted(flat.glob("*.fna"))
    if len(refined_fastas) != len(source_rows):
        raise RuntimeError("Refined FASTA/source manifest count mismatch")

    membership_rows = []
    quality_rows = []
    for fasta in refined_fastas:
        refined_id = fasta.stem
        source = source_by_id[refined_id]
        records = list(fasta_records(fasta))
        names = [name for name, _ in records]
        lengths = [int(contig_stats[name]["LengthBp"]) for name in names]
        bin_bp = sum(lengths)
        truth_bp = Counter()
        for name in names:
            truth = truth_by_contig[name]
            membership_rows.append(
                {
                    "Method": source["Method"], "RefinedID": refined_id, "SourceBin": source["SourceBin"],
                    "Contig": name, "LengthBp": contig_stats[name]["LengthBp"], "GCPct": contig_stats[name]["GCPct"],
                    "BestTruthAccession": truth["BestTruthAccession"],
                }
            )
            if truth["BestTruthAccession"] != "Unaligned":
                truth_bp[str(truth["BestTruthAccession"])] += int(contig_stats[name]["LengthBp"])
        dominant, dominant_bp = truth_bp.most_common(1)[0] if truth_bp else ("Unaligned", 0)
        assigned_bp = sum(truth_bp.values())
        purity = 100 * dominant_bp / assigned_bp if assigned_bp else math.nan
        assigned_fraction = 100 * assigned_bp / bin_bp if bin_bp else math.nan
        recovery_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
        if dominant != "Unaligned":
            for name in names:
                for reference_record, intervals in ref_intervals.get(name, {}).get(dominant, {}).items():
                    recovery_intervals[reference_record].extend(intervals)
        recovered_bp = sum(union_length(intervals) for intervals in recovery_intervals.values())
        reference_bp = int(truth_manifest[dominant]["ReferenceBases"]) if dominant != "Unaligned" else 0
        recovery = 100 * recovered_bp / reference_bp if reference_bp else math.nan
        truth_contamination = 100 - purity if not math.isnan(purity) else math.nan
        if assigned_fraction >= 80 and recovery >= 90 and truth_contamination < 5:
            truth_tier = "HQ proxy"
        elif assigned_fraction >= 80 and recovery >= 50 and truth_contamination < 10:
            truth_tier = "MQ proxy"
        else:
            truth_tier = "Below MQ proxy"
        cm, gu = checkm.get(refined_id), gunc.get(refined_id)
        if cm is None or gu is None:
            raise RuntimeError(f"QC output missing refined bin: {refined_id}")
        completeness = float(cm["Completeness"])
        contamination = float(cm["Contamination"])
        gunc_pass = as_bool(gu.get("pass.GUNC", gu.get("pass_gunc", "")))
        minimum_pass = completeness >= 50 and contamination < 10 and gunc_pass
        quality_rows.append(
            {
                "Method": source["Method"], "RefinedID": refined_id, "SourceBin": source["SourceBin"],
                "Contigs": len(names), "BinBp": bin_bp, "N50Bp": n50(lengths),
                "CheckM2Completeness": completeness, "CheckM2Contamination": contamination,
                "GUNCPass": gunc_pass, "ReferenceFreeMinimumPass": minimum_pass,
                "ReferenceFreeScore": completeness - 5 * contamination,
                "TruthAssignedBinPct": assigned_fraction, "DominantTruthAccession": dominant,
                "DominantTruthReference": truth_manifest[dominant]["Reference"] if dominant != "Unaligned" else "Unaligned",
                "AlignedPurityPct": purity, "AlignedContaminationProxyPct": truth_contamination,
                "DominantGenomeRecoveryPct": recovery, "TruthProxyTier": truth_tier,
                "RefinedFASTA": f"${{ARTICLE43_WORK_DIR}}/qc/all-refined/{fasta.name}", "RefinedSHA256": sha256(fasta),
            }
        )
    write_tsv(summary / "refined-bin-membership.tsv.gz", membership_rows)
    write_tsv(summary / "refinement-quality-truth-audit.tsv", quality_rows)

    # Trace each refined genome to every overlapping input bin. These links are
    # descriptive only; neither mock truth nor this provenance chooses a method.
    input_membership = read_tsv(frozen42 / "bin-membership.tsv.gz")
    input_by_contig: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for row in input_membership:
        input_by_contig[row["Contig"]].append((row["Branch"], row["CandidateID"]))
    links = []
    provenance = []
    for quality in quality_rows:
        refined_id = quality["RefinedID"]
        rows = [row for row in membership_rows if row["RefinedID"] == refined_id]
        overlap_bp: Counter[tuple[str, str]] = Counter()
        for row in rows:
            for branch, candidate in input_by_contig[row["Contig"]]:
                overlap_bp[(branch, candidate)] += int(row["LengthBp"])
        bin_bp = int(quality["BinBp"])
        for (branch, candidate), bp in sorted(overlap_bp.items(), key=lambda item: (-item[1], item[0])):
            links.append(
                {
                    "Method": quality["Method"], "RefinedID": refined_id, "InputBranch": branch,
                    "InputCandidateID": candidate, "OverlapBp": bp, "RefinedBinPct": 100 * bp / bin_bp,
                }
            )
        (dominant_branch, dominant_candidate), dominant_overlap = overlap_bp.most_common(1)[0]
        provenance.append(
            {
                "Method": quality["Method"], "RefinedID": refined_id,
                "ContributingBranches": len({branch for branch, _ in overlap_bp}),
                "ContributingInputBins": len(overlap_bp), "DominantInputBranch": dominant_branch,
                "DominantInputCandidateID": dominant_candidate, "DominantInputOverlapBp": dominant_overlap,
                "DominantInputCoveragePct": 100 * dominant_overlap / bin_bp,
            }
        )
    write_tsv(summary / "refinement-provenance-links.tsv.gz", links)
    write_tsv(summary / "refinement-provenance.tsv", provenance)

    method_rows = []
    for method in ("DAS Tool", "Binette"):
        selected = [row for row in quality_rows if row["Method"] == method]
        passing = [row for row in selected if row["ReferenceFreeMinimumPass"]]
        tiers = Counter(row["TruthProxyTier"] for row in selected)
        method_rows.append(
            {
                "Method": method, "Bins": len(selected), "BinnedContigs": sum(int(row["Contigs"]) for row in selected),
                "BinnedBp": sum(int(row["BinBp"]) for row in selected), "MinimumPassBins": len(passing),
                "PassingScoreSum": sum(float(row["ReferenceFreeScore"]) for row in passing),
                "HQTruthProxyBins": tiers["HQ proxy"], "MQTruthProxyBins": tiers["MQ proxy"],
                "BelowMQTruthProxyBins": tiers["Below MQ proxy"],
                "DistinctDominantTruthGenomes": len({row["DominantTruthAccession"] for row in selected if row["DominantTruthAccession"] != "Unaligned"}),
            }
        )
    write_tsv(summary / "refinement-summary.tsv", method_rows)
    selected_method_row = sorted(
        method_rows,
        key=lambda row: (int(row["MinimumPassBins"]), float(row["PassingScoreSum"]), row["Method"] == "Binette"),
        reverse=True,
    )[0]
    selected_method = selected_method_row["Method"]
    selected_quality = [
        row for row in quality_rows if row["Method"] == selected_method and row["ReferenceFreeMinimumPass"]
    ]
    selected_ids = {row["RefinedID"] for row in selected_quality}
    selected_membership = [row for row in membership_rows if row["RefinedID"] in selected_ids]
    if not selected_quality or not selected_membership:
        raise RuntimeError("Predeclared Article 43 selection yielded no final MAG candidates")
    write_tsv(summary / "selected-mag-candidates.tsv", selected_quality)
    write_tsv(summary / "selected-refinement-membership.tsv.gz", selected_membership)
    write_tsv(
        summary / "final-method-selection.tsv",
        [{
            "SelectedMethod": selected_method,
            "PrimaryCriterion": "Maximum number of CheckM2 completeness >=50%, contamination <10%, GUNC-pass bins",
            "PrimaryValue": selected_method_row["MinimumPassBins"],
            "TieBreak1": "Maximum sum of CheckM2 completeness - 5*contamination among passing bins",
            "TieBreak1Value": selected_method_row["PassingScoreSum"],
            "TieBreak2": "Binette if still tied",
            "MockTruthUsed": "No",
            "SelectedBins": len(selected_quality),
        }],
    )

    input_quality = read_tsv(frozen42 / "bin-quality-truth-audit.tsv")
    comparison_rows = []
    for branch in dict.fromkeys(row["Branch"] for row in input_quality):
        rows = [row for row in input_quality if row["Branch"] == branch]
        comparison_rows.append(
            {
                "Stage": "Input binner", "Method": branch, "Bins": len(rows),
                "MinimumPassBins": sum(as_bool(row["QCMinimumPass"]) for row in rows),
                "MedianCompleteness": median(float(row["CheckM2Completeness"]) for row in rows),
                "MedianContamination": median(float(row["CheckM2Contamination"]) for row in rows),
            }
        )
    for method in ("DAS Tool", "Binette"):
        rows = [row for row in quality_rows if row["Method"] == method]
        comparison_rows.append(
            {
                "Stage": "Refinement", "Method": method, "Bins": len(rows),
                "MinimumPassBins": sum(bool(row["ReferenceFreeMinimumPass"]) for row in rows),
                "MedianCompleteness": median(float(row["CheckM2Completeness"]) for row in rows),
                "MedianContamination": median(float(row["CheckM2Contamination"]) for row in rows),
            }
        )
    write_tsv(summary / "input-vs-refinement-summary.tsv", comparison_rows)
    resource_rows = [parse_time(path) for path in sorted((work / "logs").glob("*.time.txt"))]
    if any(int(row["ExitStatus"]) != 0 for row in resource_rows):
        raise RuntimeError("An Article 43 command has a non-zero resource exit status")
    write_tsv(summary / "resource-summary.tsv", resource_rows)
    payload = {
        "article": 43,
        "seed": 20260743,
        "coordinate_set": assembly_summary,
        "input_branches": 5,
        "refinement_methods": method_rows,
        "selected_method": selected_method,
        "selected_bins": len(selected_quality),
        "selected_binned_bp": sum(int(row["BinBp"]) for row in selected_quality),
        "truth_coordinate_rows_used": alignment_rows,
        "truth_used_for_selection": False,
        "selection_rule": "max minimum-pass bins; then max passing completeness-5*contamination sum; then Binette",
    }
    dump_json(summary / "run-summary.json", payload)
    (work / ".article43-summary-complete").write_text("PASS\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
