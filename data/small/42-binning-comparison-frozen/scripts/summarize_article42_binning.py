#!/usr/bin/env python3
"""Summarize five Article 42 binner branches with a truth-blinded mock audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

from article41_44_utils import dump_json, fasta_records, fasta_summary, parse_time, sha256, write_tsv


BRANCH_SPECS = {
    "MetaBAT2-MOCK1-only": ("bins/metabat2-MOCK1-only", None),
    "MetaBAT2-multisample": ("bins/metabat2-multisample", None),
    "SemiBin2-self-supervised": ("bins/semibin2-self-supervised/output_bins", None),
    "VAMB-taxonomy-free": ("bins/vamb-taxonomy-free/bins", None),
    "TaxVAMB-Kraken2": ("bins/taxvamb-kraken2/bins", None),
}
BRANCH_SLUG = {
    "MetaBAT2-MOCK1-only": "metabat2-mock1",
    "MetaBAT2-multisample": "metabat2-multisample",
    "SemiBin2-self-supervised": "semibin2",
    "VAMB-taxonomy-free": "vamb",
    "TaxVAMB-Kraken2": "taxvamb",
}
FASTA_SUFFIXES = {".fa", ".fna", ".fasta"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--checkm2-quality", type=Path)
    parser.add_argument("--gunc-maxcss", type=Path)
    return parser.parse_args()


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    left, right = ordered[0]
    for start, end in ordered[1:]:
        if start > right + 1:
            total += right - left + 1
            left, right = start, end
        else:
            right = max(right, end)
    return total + right - left + 1


def n50(lengths: list[int]) -> int:
    target = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            return length
    return 0


def safe_id(text: str) -> str:
    value = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()
    return value or "bin"


def parse_coords(path: Path, allowed_contigs: set[str]):
    query_intervals: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    ref_intervals: dict[str, dict[str, dict[str, list[tuple[int, int]]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    identities: dict[str, dict[str, list[tuple[float, int]]]] = defaultdict(lambda: defaultdict(list))
    parsed = 0
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            parts = [part.strip() for part in raw.split("|")]
            if len(parts) < 5:
                continue
            ref_coords = [int(value) for value in parts[0].split()[:2]]
            query_coords = [int(value) for value in parts[1].split()[:2]]
            identity = float(parts[3].split()[0])
            names = parts[4].split()
            if len(names) != 2:
                raise ValueError(f"Unexpected MetaQUAST coordinate names: {parts[4]}")
            reference_record, contig = names
            if contig not in allowed_contigs:
                continue
            match = re.search(r"GCA_\d+\.\d+", reference_record)
            if not match:
                raise ValueError(f"Reference accession is missing: {reference_record}")
            accession = match.group(0)
            q_interval = (min(query_coords), max(query_coords))
            r_interval = (min(ref_coords), max(ref_coords))
            query_intervals[contig][accession].append(q_interval)
            ref_intervals[contig][accession][reference_record].append(r_interval)
            identities[contig][accession].append((identity, q_interval[1] - q_interval[0] + 1))
            parsed += 1
    return query_intervals, ref_intervals, identities, parsed


def truth_assignments(
    contigs: dict[str, dict[str, object]],
    query_intervals,
    identities,
    truth_manifest: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    rows = []
    for contig, stats in contigs.items():
        spans = {
            accession: union_length(intervals)
            for accession, intervals in query_intervals.get(contig, {}).items()
        }
        ordered = sorted(spans.items(), key=lambda item: (-item[1], item[0]))
        if ordered:
            accession, best_span = ordered[0]
            second_span = ordered[1][1] if len(ordered) > 1 else 0
            weighted = identities[contig][accession]
            mean_identity = sum(value * length for value, length in weighted) / sum(length for _, length in weighted)
            reference = truth_manifest[accession]["Reference"]
        else:
            accession, best_span, second_span, mean_identity, reference = "Unaligned", 0, 0, math.nan, "Unaligned"
        length = int(stats["LengthBp"])
        rows.append(
            {
                "Contig": contig,
                "LengthBp": length,
                "BestTruthAccession": accession,
                "BestTruthReference": reference,
                "BestAlignedQueryBp": best_span,
                "SecondAlignedQueryBp": second_span,
                "BestAlignedFractionPct": 100 * best_span / length,
                "SecondToBestSpanRatio": second_span / best_span if best_span else 0,
                "CandidateTruthGenomes": len(ordered),
                "LengthWeightedIdentityPct": mean_identity,
                "AmbiguousTruthAssignment": int(best_span > 0 and second_span / best_span >= 0.8),
            }
        )
    return rows


def locate_bins(work: Path) -> dict[str, list[Path]]:
    result = {}
    for branch, (relative, _) in BRANCH_SPECS.items():
        directory = work / relative
        if not directory.is_dir():
            raise FileNotFoundError(f"Missing {branch} output directory: {directory}")
        paths = sorted(
            path for path in directory.iterdir()
            if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES and path.stat().st_size > 0
        )
        if not paths:
            raise RuntimeError(f"No FASTA bins for {branch}: {directory}")
        result[branch] = paths
    return result


def load_optional_qc(checkm_path: Path | None, gunc_path: Path | None):
    checkm, gunc = {}, {}
    if checkm_path:
        rows = list(csv.DictReader(checkm_path.open(encoding="utf-8"), delimiter="\t"))
        for row in rows:
            name = row.get("Name") or row.get("name") or row.get("Genome")
            if name:
                checkm[name.removesuffix(".fna")] = row
    if gunc_path:
        rows = list(csv.DictReader(gunc_path.open(encoding="utf-8"), delimiter="\t"))
        for row in rows:
            name = row.get("genome") or row.get("Genome")
            if name:
                gunc[name.removesuffix(".fna")] = row
    return checkm, gunc


def comb2(value: int) -> int:
    return value * (value - 1) // 2


def adjusted_rand(assign_a: dict[str, str], assign_b: dict[str, str]) -> tuple[float, int]:
    shared = sorted(set(assign_a) & set(assign_b))
    n = len(shared)
    if n < 2:
        return math.nan, n
    contingency = Counter((assign_a[name], assign_b[name]) for name in shared)
    row_totals = Counter(assign_a[name] for name in shared)
    col_totals = Counter(assign_b[name] for name in shared)
    sum_comb = sum(comb2(value) for value in contingency.values())
    rows_comb = sum(comb2(value) for value in row_totals.values())
    cols_comb = sum(comb2(value) for value in col_totals.values())
    total_comb = comb2(n)
    expected = rows_comb * cols_comb / total_comb if total_comb else 0
    maximum = 0.5 * (rows_comb + cols_comb)
    return ((sum_comb - expected) / (maximum - expected) if maximum != expected else 1.0), n


def main() -> int:
    args = parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    if not (work / ".article42-run-complete").is_file():
        raise FileNotFoundError("Run run_article42_binning.py first")
    summary_dir = work / "summary"
    summary_dir.mkdir(exist_ok=True)
    common_fasta = work / "inputs/megahit-coassembly.ge1500.fna"
    assembly_summary, contigs = fasta_summary(common_fasta)
    sequences = dict(fasta_records(common_fasta))
    if set(sequences) != set(contigs):
        raise RuntimeError("Common FASTA sequence/statistic coordinate mismatch")
    truth_path = root / "data/raw/article33/work/metaquast/MOCK1_MOCK2/combined_reference/contigs_reports/minimap_output/sr-megahit-co.coords"
    manifest_path = root / "data/small/33-assembly-qc-frozen/truth-manifest.tsv"
    if not truth_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError("Article 33 truth coordinates/manifest are required for the post-hoc audit")
    manifest_rows = list(csv.DictReader(manifest_path.open(encoding="utf-8"), delimiter="\t"))
    truth_manifest = {
        row["GenBankAssembly"]: row
        for row in manifest_rows
        if row["EvaluationSet"] == "MOCK1+MOCK2"
    }
    if len(truth_manifest) != 87:
        raise RuntimeError(f"Expected 87 co-assembly truth genomes, observed {len(truth_manifest)}")
    query_intervals, ref_intervals, identities, alignment_rows = parse_coords(truth_path, set(contigs))
    truth_rows = truth_assignments(contigs, query_intervals, identities, truth_manifest)
    truth_by_contig = {row["Contig"]: row for row in truth_rows}
    write_tsv(summary_dir / "truth-contig-assignment.tsv.gz", truth_rows)
    write_tsv(
        summary_dir / "truth-input-audit.tsv",
        [
            {"Role": "MetaQUAST-coassembly-coordinates", "Path": str(truth_path.relative_to(root)), "Bytes": truth_path.stat().st_size, "SHA256": sha256(truth_path), "Status": "PASS"},
            {"Role": "Article33-truth-manifest", "Path": str(manifest_path.relative_to(root)), "Bytes": manifest_path.stat().st_size, "SHA256": sha256(manifest_path), "Status": "PASS"},
        ],
    )

    bins_by_branch = locate_bins(work)
    candidate_dir = work / "candidates"
    if candidate_dir.exists():
        shutil.rmtree(candidate_dir)
    candidate_dir.mkdir()
    checkm, gunc = load_optional_qc(args.checkm2_quality, args.gunc_maxcss)
    membership_rows = []
    bin_rows = []
    branch_assignments: dict[str, dict[str, str]] = {}
    total_coordinate_bp = int(assembly_summary["TotalBp"])
    for branch, paths in bins_by_branch.items():
        branch_seen: set[str] = set()
        assignments: dict[str, str] = {}
        branch_candidate_dir = candidate_dir / BRANCH_SLUG[branch]
        branch_candidate_dir.mkdir()
        for position, path in enumerate(paths, start=1):
            source_records = list(fasta_records(path))
            names = [name for name, _ in source_records]
            if len(names) != len(set(names)):
                raise RuntimeError(f"Duplicate contig in {path}")
            unknown = set(names) - set(contigs)
            overlap = set(names) & branch_seen
            if unknown or overlap:
                raise RuntimeError(f"Invalid bin membership in {path}: unknown={len(unknown)}, overlap={len(overlap)}")
            branch_seen.update(names)
            candidate_id = f"{BRANCH_SLUG[branch]}__{position:03d}"
            candidate_path = branch_candidate_dir / f"{candidate_id}.fna"
            with candidate_path.open("w", encoding="utf-8", newline="\n") as handle:
                for name in sorted(names):
                    sequence = sequences[name]
                    handle.write(f">{name}\n")
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80] + "\n")
            lengths = [int(contigs[name]["LengthBp"]) for name in names]
            bin_bp = sum(lengths)
            if bin_bp < 200_000:
                raise RuntimeError(f"Output bin below locked 200-kb threshold: {branch} {path} {bin_bp}")
            truth_bp = Counter()
            for name in names:
                assignments[name] = candidate_id
                truth = truth_by_contig[name]
                membership_rows.append(
                    {
                        "Branch": branch,
                        "CandidateID": candidate_id,
                        "SourceBin": path.name,
                        "Contig": name,
                        "LengthBp": contigs[name]["LengthBp"],
                        "GCPct": contigs[name]["GCPct"],
                        "BestTruthAccession": truth["BestTruthAccession"],
                    }
                )
                if truth["BestTruthAccession"] != "Unaligned":
                    truth_bp[truth["BestTruthAccession"]] += int(contigs[name]["LengthBp"])
            dominant, dominant_bp = truth_bp.most_common(1)[0] if truth_bp else ("Unaligned", 0)
            truth_assigned_bp = sum(truth_bp.values())
            aligned_purity = 100 * dominant_bp / truth_assigned_bp if truth_assigned_bp else math.nan
            assigned_fraction = 100 * truth_assigned_bp / bin_bp
            recovery_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
            if dominant != "Unaligned":
                for name in names:
                    for record, intervals in ref_intervals.get(name, {}).get(dominant, {}).items():
                        recovery_intervals[record].extend(intervals)
            recovered_reference_bp = sum(union_length(intervals) for intervals in recovery_intervals.values())
            reference_bp = int(truth_manifest[dominant]["ReferenceBases"]) if dominant != "Unaligned" else 0
            recovery = 100 * recovered_reference_bp / reference_bp if reference_bp else math.nan
            contamination = 100 - aligned_purity if not math.isnan(aligned_purity) else math.nan
            if assigned_fraction >= 80 and recovery >= 90 and contamination < 5:
                truth_tier = "HQ proxy"
            elif assigned_fraction >= 80 and recovery >= 50 and contamination < 10:
                truth_tier = "MQ proxy"
            else:
                truth_tier = "Below MQ proxy"
            cm = checkm.get(candidate_id, {})
            gu = gunc.get(candidate_id, {})
            completeness = float(cm.get("Completeness", "nan")) if cm else math.nan
            checkm_contamination = float(cm.get("Contamination", "nan")) if cm else math.nan
            pass_gunc_text = gu.get("pass.GUNC", gu.get("pass_gunc", "")) if gu else ""
            pass_gunc = pass_gunc_text.lower() == "true" if pass_gunc_text else None
            qc_pass = (
                completeness >= 50 and checkm_contamination < 10 and pass_gunc is True
                if cm and gu else None
            )
            bin_rows.append(
                {
                    "Branch": branch,
                    "CandidateID": candidate_id,
                    "SourceBin": path.name,
                    "Contigs": len(names),
                    "BinBp": bin_bp,
                    "N50Bp": n50(lengths),
                    "TruthAssignedBinPct": assigned_fraction,
                    "DominantTruthAccession": dominant,
                    "DominantTruthReference": truth_manifest[dominant]["Reference"] if dominant != "Unaligned" else "Unaligned",
                    "DominantTruthContigBp": dominant_bp,
                    "AlignedPurityPct": aligned_purity,
                    "AlignedContaminationProxyPct": contamination,
                    "RecoveredReferenceBp": recovered_reference_bp,
                    "ReferenceBp": reference_bp,
                    "DominantGenomeRecoveryPct": recovery,
                    "TruthProxyTier": truth_tier,
                    "CheckM2Completeness": completeness,
                    "CheckM2Contamination": checkm_contamination,
                    "GUNCPass": "" if pass_gunc is None else pass_gunc,
                    "QCMinimumPass": "" if qc_pass is None else qc_pass,
                    "CandidateFASTA": f"${{ARTICLE42_WORK_DIR}}/candidates/{BRANCH_SLUG[branch]}/{candidate_id}.fna",
                    "CandidateSHA256": sha256(candidate_path),
                }
            )
        branch_assignments[branch] = assignments
    write_tsv(summary_dir / "bin-membership.tsv.gz", membership_rows)
    write_tsv(summary_dir / "bin-quality-truth-audit.tsv", bin_rows)

    branch_rows = []
    for branch in BRANCH_SPECS:
        selected = [row for row in bin_rows if row["Branch"] == branch]
        binned_contigs = {row["Contig"] for row in membership_rows if row["Branch"] == branch}
        tiers = Counter(row["TruthProxyTier"] for row in selected)
        quality_available = any(row["QCMinimumPass"] != "" for row in selected)
        branch_rows.append(
            {
                "Branch": branch,
                "Bins": len(selected),
                "BinnedContigs": len(binned_contigs),
                "BinnedBp": sum(int(contigs[name]["LengthBp"]) for name in binned_contigs),
                "BinnedAssemblyPct": 100 * sum(int(contigs[name]["LengthBp"]) for name in binned_contigs) / total_coordinate_bp,
                "HQTruthProxyBins": tiers["HQ proxy"],
                "MQTruthProxyBins": tiers["MQ proxy"],
                "BelowMQTruthProxyBins": tiers["Below MQ proxy"],
                "DistinctDominantTruthGenomes": len({row["DominantTruthAccession"] for row in selected if row["DominantTruthAccession"] != "Unaligned"}),
                "CheckM2GUNCQCPassBins": sum(str(row["QCMinimumPass"]).lower() == "true" for row in selected) if quality_available else "",
            }
        )
    write_tsv(summary_dir / "binner-summary.tsv", branch_rows)

    ari, shared = adjusted_rand(
        branch_assignments["MetaBAT2-MOCK1-only"],
        branch_assignments["MetaBAT2-multisample"],
    )
    single = set(branch_assignments["MetaBAT2-MOCK1-only"])
    multi = set(branch_assignments["MetaBAT2-multisample"])
    write_tsv(
        summary_dir / "single-vs-multisample.tsv",
        [
            {
                "Comparison": "MetaBAT2 MOCK1-only vs multisample",
                "SingleBinnedContigs": len(single),
                "MultiBinnedContigs": len(multi),
                "SharedBinnedContigs": len(single & multi),
                "SingleOnlyContigs": len(single - multi),
                "MultiOnlyContigs": len(multi - single),
                "AdjustedRandOnSharedBinnedContigs": ari,
                "AdjustedRandContigs": shared,
            }
        ],
    )
    taxonomy_rows = list(csv.DictReader((work / "taxonomy/taxonomy-summary.tsv").open(encoding="utf-8"), delimiter="\t"))
    write_tsv(summary_dir / "taxonomy-summary.tsv", taxonomy_rows)
    resource_rows = [parse_time(path) for path in sorted((work / "logs").glob("*.time.txt"))]
    override_path = work / "resource-overrides.tsv"
    if override_path.is_file():
        overrides = {
            row["Label"]: row
            for row in csv.DictReader(override_path.open(encoding="utf-8"), delimiter="\t")
        }
        resource_rows = [
            {
                **row,
                **({
                    "WallSeconds": float(overrides[row["Label"]]["WallSeconds"]),
                    "PeakRAMGiB": overrides[row["Label"]]["PeakRAMGiB"],
                    "ExitStatus": int(overrides[row["Label"]]["ExitStatus"]),
                    "Command": overrides[row["Label"]]["Command"],
                    "MeasurementStatus": overrides[row["Label"]]["MeasurementStatus"],
                } if row["Label"] in overrides else {}),
            }
            for row in resource_rows
        ]
    if any(int(row["ExitStatus"]) != 0 for row in resource_rows):
        raise RuntimeError("A binner command has a non-zero resource exit status")
    write_tsv(summary_dir / "resource-summary.tsv", resource_rows)
    payload = {
        "article": 42,
        "seed": 20260742,
        "coordinate_set": assembly_summary,
        "truth_genomes": len(truth_manifest),
        "truth_coordinate_rows_used": alignment_rows,
        "truth_assigned_contigs": sum(row["BestTruthAccession"] != "Unaligned" for row in truth_rows),
        "ambiguous_truth_contigs": sum(int(row["AmbiguousTruthAssignment"]) for row in truth_rows),
        "branches": branch_rows,
        "single_vs_multi": {
            "shared_binned_contigs": len(single & multi),
            "adjusted_rand": ari,
            "adjusted_rand_contigs": shared,
        },
        "taxonomy": taxonomy_rows,
        "quality_databases_applied": bool(checkm and gunc),
        "boundaries": {
            "truth_never_used_by_binners": True,
            "truth_proxy_is_not_checkm2_or_mimag": True,
            "kraken_taxonomy_only_used_by_taxvamb": True,
            "checkm2_gunc_do_not_choose_taxonomic_correctness": True,
        },
    }
    dump_json(summary_dir / "run-summary.json", payload)
    (work / ".article42-summary-complete").write_text("PASS\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
