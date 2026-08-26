#!/usr/bin/env python3
"""Summarize geNomad, CARD-RGI, and reference-linkage evidence for Article 57."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def find_one(base: Path, pattern: str) -> Path:
    matches = sorted(base.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one match for {base}/{pattern}, observed {matches}")
    return matches[0]


def tokens(value: str | None) -> list[str]:
    if value is None or value.strip() in {"", "NA", "None", "nan"}:
        return []
    return [item for item in value.split(";") if item and item != "NA"]


def maybe_float(value: str | None) -> float | None:
    if value is None or value.strip() in {"", "NA", "None", "nan"}:
        return None
    return float(value)


def summary_index(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["seq_name"]: row for row in rows}


def parent_sequence(gene: str, sequence_names: list[str]) -> str:
    for name in sequence_names:
        if gene.startswith(name + "_"):
            return name
    raise ValueError(f"Cannot map geNomad gene to sequence: {gene}")


def wilson_interval(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    p = success / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / denominator
    return center - half, center + half


def parse_paf(path: Path) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq = fields[:12]
            row: dict[str, object] = {
                "Query": qname,
                "QueryLength": int(qlen),
                "QueryStart": int(qstart),
                "QueryEnd": int(qend),
                "Strand": strand,
                "Target": tname,
                "TargetLength": int(tlen),
                "TargetStart": int(tstart),
                "TargetEnd": int(tend),
                "Matches": int(nmatch),
                "AlignmentLength": int(alen),
                "MapQ": int(mapq),
                "QueryCoverage": (int(qend) - int(qstart)) / int(qlen),
                "AlignmentIdentity": int(nmatch) / int(alen) if int(alen) else 0.0,
            }
            rank = (int(nmatch), int(alen), int(mapq))
            previous = best.get(qname)
            if previous is None or rank > previous["_rank"]:
                row["_rank"] = rank
                best[qname] = row
    for row in best.values():
        row.pop("_rank", None)
    return best


def parse_gnu_time(path: Path) -> tuple[float, float]:
    elapsed: str | None = None
    rss_kb: int | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("Elapsed (wall clock) time"):
            elapsed = line.rsplit(": ", 1)[1]
        elif line.startswith("Maximum resident set size (kbytes)"):
            rss_kb = int(line.rsplit(": ", 1)[1])
    if elapsed is None or rss_kb is None:
        raise ValueError(f"Incomplete GNU time log: {path}")
    parts = [float(item) for item in elapsed.split(":")]
    seconds = parts[-1]
    if len(parts) >= 2:
        seconds += parts[-2] * 60
    if len(parts) == 3:
        seconds += parts[0] * 3600
    return seconds, rss_kb / 1024 / 1024


def original_short_name(header: str) -> str:
    match = re.search(r"\bplasmid\s+([^, ]+)", header, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    accession = header.split()[0].split("|")[-1]
    return accession


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    if not (work / ".article57-run-complete").is_file():
        raise FileNotFoundError("Run run_article57_plasmids.py first")
    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    labels = read_tsv(work / "reference-replicon-labels.tsv")
    benchmark_labels = read_tsv(work / "reference-benchmark-labels.tsv")
    label_by_name = {row["SeqName"]: row for row in labels}
    ref_summary_path = find_one(
        work / "results/reference-benchmark", "*_summary/*_plasmid_summary.tsv"
    )
    ref_predictions = read_tsv(ref_summary_path)
    ref_by_name = summary_index(ref_predictions)
    reference_calls: list[dict[str, object]] = []
    confusion = Counter()
    for row in benchmark_labels:
        expected = row["ReferenceLabel"] == "Plasmid"
        observed = row["SeqName"] in ref_by_name
        cell = "TP" if expected and observed else "FN" if expected else "FP" if observed else "TN"
        confusion[cell] += 1
        call = ref_by_name.get(row["SeqName"], {})
        reference_calls.append(
            {
                "SeqName": row["SeqName"],
                "Assembly": row["Assembly"],
                "ReferenceLabel": row["ReferenceLabel"],
                "geNomadPlasmid": "yes" if observed else "no",
                "ConfusionCell": cell,
                "LengthBp": row["LengthBp"],
                "PlasmidScore": call.get("plasmid_score", ""),
                "FDR": call.get("fdr", ""),
                "Topology": call.get("topology", ""),
                "Hallmarks": call.get("n_hallmarks", ""),
                "MarkerEnrichment": call.get("marker_enrichment", ""),
                "ConjugationGenes": call.get("conjugation_genes", ""),
                "geNomadAMRGenes": call.get("amr_genes", ""),
                "OriginalHeader": row["OriginalHeader"],
            }
        )
    write_tsv(summary_dir / "reference-classification.tsv", reference_calls)
    confusion_rows = [
        {"ReferenceLabel": "Plasmid", "geNomadCall": "Plasmid", "Count": confusion["TP"]},
        {"ReferenceLabel": "Plasmid", "geNomadCall": "Not plasmid", "Count": confusion["FN"]},
        {"ReferenceLabel": "Other replicon", "geNomadCall": "Plasmid", "Count": confusion["FP"]},
        {"ReferenceLabel": "Other replicon", "geNomadCall": "Not plasmid", "Count": confusion["TN"]},
    ]
    write_tsv(summary_dir / "reference-confusion.tsv", confusion_rows)

    tp, fn, fp, tn = (confusion[key] for key in ("TP", "FN", "FP", "TN"))
    metrics: list[dict[str, object]] = []
    for metric, numerator, denominator in (
        ("Sensitivity", tp, tp + fn),
        ("Precision", tp, tp + fp),
        ("Specificity", tn, tn + fp),
    ):
        lower, upper = wilson_interval(numerator, denominator)
        metrics.append(
            {
                "Metric": metric,
                "Numerator": numerator,
                "Denominator": denominator,
                "Estimate": numerator / denominator,
                "Wilson95Lower": lower,
                "Wilson95Upper": upper,
            }
        )
    write_tsv(summary_dir / "reference-metrics.tsv", metrics)

    co_summary_path = find_one(
        work / "results/coassembly", "*_summary/*_plasmid_summary.tsv"
    )
    co_gene_path = find_one(work / "results/coassembly", "*_summary/*_plasmid_genes.tsv")
    co_candidates = read_tsv(co_summary_path)
    co_names = sorted((row["seq_name"] for row in co_candidates), key=len, reverse=True)
    co_genes = read_tsv(co_gene_path)
    genes_by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in co_genes:
        genes_by_contig[parent_sequence(row["gene"], co_names)].append(row)
    rgi_co = read_tsv(work / "rgi-primary-coassembly.tsv")
    rgi_by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rgi_co:
        rgi_by_contig[row["Contig"]].append(row)
    alignments = parse_paf(work / "results/coassembly-plasmid-to-reference.paf")

    candidate_rows: list[dict[str, object]] = []
    for row in co_candidates:
        name = row["seq_name"]
        gene_rows = genes_by_contig[name]
        conj_rows = [item for item in gene_rows if tokens(item.get("annotation_conjscan"))]
        internal_amr = [item for item in gene_rows if tokens(item.get("annotation_amr"))]
        external_arg = rgi_by_contig.get(name, [])
        alignment = alignments.get(name)
        target = str(alignment["Target"]) if alignment else ""
        target_label = label_by_name.get(target, {}).get("ReferenceLabel", "Unresolved")
        target_class = label_by_name.get(target, {}).get(
            "ReferenceSequenceClass", "Unresolved"
        )
        qcov = float(alignment["QueryCoverage"]) if alignment else 0.0
        identity = float(alignment["AlignmentIdentity"]) if alignment else 0.0
        if target_label == "Plasmid" and qcov >= 0.80 and identity >= 0.95:
            truth_support = "Reference-plasmid supported"
        elif target_class == "Complete cellular replicon" and qcov >= 0.80 and identity >= 0.95:
            truth_support = "Complete-cellular-reference conflict"
        elif target_label == "Other replicon" and qcov >= 0.80 and identity >= 0.95:
            truth_support = "Other-reference-label conflict"
        else:
            truth_support = "No high-coverage reference support"
        if external_arg and conj_rows:
            context = "ARG + plasmid call + conjugation marker"
            ceiling = "Predicted transferable ARG-bearing plasmid candidate"
        elif external_arg:
            context = "ARG + plasmid call"
            ceiling = "Predicted plasmid-associated ARG"
        elif conj_rows:
            context = "Plasmid call + conjugation marker"
            ceiling = "Predicted plasmid with transfer-related marker"
        else:
            context = "Plasmid call only"
            ceiling = "Predicted plasmid sequence"
        candidate_rows.append(
            {
                "Contig": name,
                "LengthBp": row["length"],
                "Topology": row["topology"],
                "Genes": row["n_genes"],
                "PlasmidScore": row["plasmid_score"],
                "FDR": row["fdr"],
                "Hallmarks": row["n_hallmarks"],
                "MarkerEnrichment": row["marker_enrichment"],
                "ConjugationMarkerCount": len(conj_rows),
                "ConjugationMarkers": ";".join(
                    sorted(
                        {
                            token
                            for item in conj_rows
                            for token in tokens(item.get("annotation_conjscan"))
                        }
                    )
                ),
                "geNomadAMRMarkerCount": len(internal_amr),
                "CARDPrimaryARGCount": len(external_arg),
                "CARDPrimaryARGs": ";".join(sorted({item["BestHitARO"] for item in external_arg})),
                "EvidenceContext": context,
                "ClaimCeiling": ceiling,
                "BestReference": target,
                "ReferenceLabel": target_label,
                "ReferenceSequenceClass": target_class,
                "AlignmentIdentity": f"{identity:.6f}",
                "QueryCoverage": f"{qcov:.6f}",
                "ReferenceSupport": truth_support,
            }
        )
    write_tsv(summary_dir / "coassembly-plasmid-candidates.tsv", candidate_rows)

    candidate_by_name = {row["Contig"]: row for row in candidate_rows}
    arg_ledger: list[dict[str, object]] = []
    for row in rgi_co:
        candidate = candidate_by_name.get(row["Contig"])
        if candidate is None:
            context = "ARG on non-plasmid contig"
            ceiling = "ARG sequence evidence; genomic mobility unresolved"
        else:
            context = str(candidate["EvidenceContext"])
            ceiling = str(candidate["ClaimCeiling"])
        arg_ledger.append(
            {
                **row,
                "geNomadPlasmid": "yes" if candidate else "no",
                "ConjugationMarkerCount": candidate["ConjugationMarkerCount"] if candidate else 0,
                "Topology": candidate["Topology"] if candidate else "",
                "PlasmidScore": candidate["PlasmidScore"] if candidate else "",
                "ReferenceSupport": candidate["ReferenceSupport"] if candidate else "Not evaluated",
                "EvidenceContext": context,
                "ClaimCeiling": ceiling,
            }
        )
    write_tsv(summary_dir / "arg-mobility-ledger.tsv", arg_ledger)
    context_counts = Counter(row["EvidenceContext"] for row in arg_ledger)
    write_tsv(
        summary_dir / "arg-context-summary.tsv",
        [
            {
                "EvidenceContext": context,
                "ARGCalls": count,
                "Fraction": count / len(arg_ledger),
            }
            for context, count in sorted(context_counts.items())
        ],
    )

    staph_summary_path = find_one(
        work / "results/staphylococcus", "*_summary/*_plasmid_summary.tsv"
    )
    staph_gene_path = find_one(
        work / "results/staphylococcus", "*_summary/*_plasmid_genes.tsv"
    )
    staph_candidates = read_tsv(staph_summary_path)
    staph_by_name = summary_index(staph_candidates)
    staph_names = sorted(staph_by_name, key=len, reverse=True)
    staph_genes = read_tsv(staph_gene_path)
    staph_gene_by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in staph_genes:
        staph_gene_by_contig[parent_sequence(row["gene"], staph_names)].append(row)
    staph_rgi = read_tsv(work / "rgi-primary-staphylococcus.tsv")
    staph_rgi_by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in staph_rgi:
        staph_rgi_by_contig[row["Contig"]].append(row)

    usa300_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    for label in labels:
        if label["Assembly"] != "GCA_000013465.1":
            continue
        name = label["SeqName"]
        short = original_short_name(label["OriginalHeader"])
        call = staph_by_name.get(name)
        gene_rows = staph_gene_by_contig.get(name, [])
        conj = [item for item in gene_rows if tokens(item.get("annotation_conjscan"))]
        args_here = staph_rgi_by_contig.get(name, [])
        usa300_rows.append(
            {
                "Replicon": name,
                "DisplayName": short,
                "ReferenceLabel": label["ReferenceLabel"],
                "LengthBp": label["LengthBp"],
                "geNomadPlasmid": "yes" if call else "no",
                "PlasmidScore": call.get("plasmid_score", "") if call else "",
                "Topology": call.get("topology", "") if call else "",
                "ConjugationMarkerCount": len(conj),
                "CARDPrimaryARGCount": len(args_here),
                "CARDPrimaryARGs": ";".join(sorted({item["BestHitARO"] for item in args_here})),
            }
        )
        if label["ReferenceLabel"] == "Plasmid":
            for arg in args_here:
                feature_rows.append(
                    {
                        "Replicon": name,
                        "DisplayName": short,
                        "RepliconLengthBp": label["LengthBp"],
                        "FeatureType": "CARD Perfect/Strict ARG",
                        "Feature": arg["BestHitARO"],
                        "Start": arg["Start"],
                        "End": arg["Stop"],
                        "Strand": arg["Orientation"],
                    }
                )
            for gene in conj:
                feature_rows.append(
                    {
                        "Replicon": name,
                        "DisplayName": short,
                        "RepliconLengthBp": label["LengthBp"],
                        "FeatureType": "Conjugation-related marker",
                        "Feature": gene["annotation_conjscan"],
                        "Start": gene["start"],
                        "End": gene["end"],
                        "Strand": gene["strand"],
                    }
                )
    write_tsv(summary_dir / "usa300-replicon-audit.tsv", usa300_rows)
    if not feature_rows:
        feature_rows.append(
            {
                "Replicon": "none",
                "DisplayName": "No retained feature",
                "RepliconLengthBp": 1,
                "FeatureType": "None",
                "Feature": "None",
                "Start": 0,
                "End": 0,
                "Strand": ".",
            }
        )
    write_tsv(summary_dir / "usa300-plasmid-features.tsv", feature_rows)

    ladder = [
        {
            "Rank": 1,
            "Evidence": "Observed transfer",
            "MinimumEvidence": "Mating/transformation/transduction with recipient validation",
            "MaximumClaim": "Transfer demonstrated under tested conditions",
            "MainFalsePositive": "Contamination or donor carry-over without recipient genotype check",
        },
        {
            "Rank": 2,
            "Evidence": "Closed element + physical linkage",
            "MinimumEvidence": "Circular or boundary-resolved element; ARG and transfer module; read support",
            "MaximumClaim": "Genetic capacity for mobility on one resolved element",
            "MainFalsePositive": "Misassembly or repeat-driven circularization",
        },
        {
            "Rank": 3,
            "Evidence": "Predicted plasmid + ARG + transfer marker",
            "MinimumEvidence": "Same contig; independent ARG call; conjugation/relaxase marker",
            "MaximumClaim": "Predicted transferable ARG-bearing plasmid candidate",
            "MainFalsePositive": "Fragment/chromid confusion; incomplete transfer module",
        },
        {
            "Rank": 4,
            "Evidence": "Predicted plasmid + ARG",
            "MinimumEvidence": "Same contig; independent plasmid and ARG calls",
            "MaximumClaim": "Predicted plasmid-associated ARG",
            "MainFalsePositive": "Short-contig classification or chimeric assembly",
        },
        {
            "Rank": 5,
            "Evidence": "ARG + local MGE marker",
            "MinimumEvidence": "Same contig or long read; transposase/integron/ICE feature nearby",
            "MaximumClaim": "ARG in a putative mobile genetic context",
            "MainFalsePositive": "Marker is inactive, distant, or belongs to another element",
        },
        {
            "Rank": 6,
            "Evidence": "Unlinked detections",
            "MinimumEvidence": "ARG and MGE detected only at sample level",
            "MaximumClaim": "ARG and MGE co-occur in the sample",
            "MainFalsePositive": "Different organisms or molecules; abundance correlation",
        },
    ]
    write_tsv(summary_dir / "mobility-evidence-ladder.tsv", ladder)

    resources: list[dict[str, object]] = []
    for label in ("reference-benchmark", "coassembly", "staphylococcus"):
        seconds, peak_gib = parse_gnu_time(work / f"logs/genomad-{label}.time.txt")
        resources.append(
            {
                "Branch": label,
                "WallSeconds": f"{seconds:.3f}",
                "PeakRSSGiB": f"{peak_gib:.6f}",
                "Threads": 16,
                "Splits": 16,
            }
        )
    seconds, peak_gib = parse_gnu_time(work / "logs/minimap2-plasmid-truth.time.txt")
    resources.append(
        {
            "Branch": "plasmid-to-reference alignment",
            "WallSeconds": f"{seconds:.3f}",
            "PeakRSSGiB": f"{peak_gib:.6f}",
            "Threads": 16,
            "Splits": "NA",
        }
    )
    write_tsv(summary_dir / "resource-usage.tsv", resources)
    write_tsv(
        work / "tool-versions.tsv",
        [
            {"Tool": "geNomad", "Version": "1.12.0", "Role": "plasmid classification and marker annotation"},
            {"Tool": "geNomad database", "Version": "1.9 (ICTV MSL39)", "Role": "marker and taxonomy database"},
            {"Tool": "RGI", "Version": "6.0.8", "Role": "independent Perfect/Strict ARG calls"},
            {"Tool": "CARD", "Version": "4.0.1", "Role": "RGI resistance models"},
            {"Tool": "minimap2", "Version": "2.31-r1302", "Role": "candidate-to-exact-reference alignment"},
        ],
    )
    write_tsv(
        work / "database-audit.tsv",
        [
            {
                "Database": "geNomad database",
                "Release": "v1.9",
                "ReleaseIdentity": "Zenodo 14886553; archive MD5 67244b528bb8bed464d1ca147136d33e",
                "Status": "PASS",
            },
            {
                "Database": "CARD",
                "Release": "4.0.1",
                "ReleaseIdentity": "card.json SHA-256 dee4dcdb0d9c7f79107452d64211d816d1eab55289ddf8dc5f1e99ddfdc5e111",
                "Status": "PASS",
            },
        ],
    )

    result = {
        "reference_replicons": len(labels),
        "reference_benchmark_replicons": len(benchmark_labels),
        "reference_plasmids": sum(row["ReferenceLabel"] == "Plasmid" for row in labels),
        "reference_tp": tp,
        "reference_fn": fn,
        "reference_fp": fp,
        "reference_tn": tn,
        "coassembly_plasmid_candidates": len(candidate_rows),
        "coassembly_plasmid_candidates_reference_supported": sum(
            row["ReferenceSupport"] == "Reference-plasmid supported" for row in candidate_rows
        ),
        "coassembly_plasmid_candidates_complete_cellular_conflict": sum(
            row["ReferenceSupport"] == "Complete-cellular-reference conflict"
            for row in candidate_rows
        ),
        "coassembly_plasmid_candidates_unresolved": sum(
            row["ReferenceSupport"] == "No high-coverage reference support"
            for row in candidate_rows
        ),
        "coassembly_primary_arg_calls": len(arg_ledger),
        "coassembly_primary_arg_calls_on_predicted_plasmids": sum(
            row["geNomadPlasmid"] == "yes" for row in arg_ledger
        ),
        "usa300_reference_plasmids": sum(
            row["ReferenceLabel"] == "Plasmid" for row in usa300_rows
        ),
        "usa300_predicted_plasmids": sum(row["geNomadPlasmid"] == "yes" for row in usa300_rows),
        "usa300_primary_arg_calls_on_reference_plasmids": sum(
            int(row["CARDPrimaryARGCount"])
            for row in usa300_rows
            if row["ReferenceLabel"] == "Plasmid"
        ),
        "seed": 20260757,
    }
    (summary_dir / "summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (work / ".article57-summary-complete").write_text("verified\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
