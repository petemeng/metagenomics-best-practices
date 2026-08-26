#!/usr/bin/env python3
"""Summarize Article 58 EukDetect2 and EukRep evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator


REFERENCE_EUKARYOTES = {"Saccharomyces cerevisiae", "Cryptococcus neoformans"}
EXPECTED_GROUP_ALIASES = {
    "Saccharomyces cerevisiae": {"Saccharomyces cerevisiae"},
    "Cryptococcus neoformans complex": {
        "Cryptococcus neoformans",
        "Cryptococcus deneoformans",
    },
}
MODES = ("strict", "balanced", "lenient")
FRAGMENT_LENGTHS = (3000, 5000, 10000, 20000)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_float(value: str) -> float:
    text = value.strip().rstrip("%")
    if text in {"", "NA", "None"}:
        return math.nan
    return float(text)


def read_names(path: Path) -> set[str]:
    names: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            name = line.strip()
            if name:
                names.add(name.lstrip(">").split()[0])
    return names


def parse_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].split()[0]
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before header in {path}")
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def safe_div(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else math.nan


def mock_target(name: str) -> str | None:
    for target, aliases in EXPECTED_GROUP_ALIASES.items():
        if name in aliases:
            return target
    return None


def elapsed_seconds(value: str) -> float:
    fields = value.strip().split(":")
    if len(fields) == 2:
        minutes, seconds = fields
        return 60 * float(minutes) + float(seconds)
    if len(fields) == 3:
        hours, minutes, seconds = fields
        return 3600 * float(hours) + 60 * float(minutes) + float(seconds)
    raise ValueError(f"Unsupported GNU time elapsed value: {value}")


def parse_time_file(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8")
    fields: dict[str, str] = {}
    for line in text.splitlines():
        if ": " not in line:
            continue
        key, value = line.strip().rsplit(": ", 1)
        fields[key] = value
    elapsed_key = next(
        (key for key in fields if key.startswith("Elapsed (wall clock) time")), None
    )
    if elapsed_key is None:
        raise RuntimeError(f"Cannot find elapsed time in {path}")
    user = float(fields["User time (seconds)"])
    system = float(fields["System time (seconds)"])
    peak_kib = float(fields["Maximum resident set size (kbytes)"])
    return {
        "WallSeconds": elapsed_seconds(fields[elapsed_key]),
        "CPUSeconds": user + system,
        "PeakRSSGiB": peak_kib / 1024 / 1024,
    }


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    if not (work / ".article58-run-complete").is_file():
        raise FileNotFoundError("Run run_article58_eukaryotes.py first")
    results = work / "results"
    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    primary = read_tsv(results / "eukdetect" / "Zymo_D6300_filtered_hits_table.txt")
    all_hits = read_tsv(results / "eukdetect" / "filtering" / "Zymo_D6300_all_hits_table.txt")
    eukfrac = read_tsv(results / "eukdetect" / "Zymo_D6300_filtered_hits_eukfrac.txt")
    normalized = read_tsv(results / "eukdetect" / "Zymo_D6300.normalized.tsv")
    all_by_taxid = {row["Taxid"]: row for row in all_hits}
    normalized_by_taxid = {row["TaxID"]: row for row in normalized}
    eukfrac_by_taxid = {row["TaxID"]: row for row in eukfrac}
    evidence: list[dict[str, object]] = []
    for row in primary:
        if row["Rank"] != "species":
            continue
        taxid = row["Taxid"]
        raw = all_by_taxid.get(taxid, {})
        norm = normalized_by_taxid.get(taxid, {})
        frac = eukfrac_by_taxid.get(taxid, {})
        evidence.append(
            {
                "Name": row["Name"],
                "TaxID": taxid,
                "ExpectedInMock": str(mock_target(row["Name"]) is not None).lower(),
                "MockTarget": mock_target(row["Name"]) or "Unexpected eukaryote",
                "ObservedMarkers": int(float(raw.get("Observed_markers", "0"))),
                "TotalMarkerReads": float(row["Total_reads"]),
                "DirectMarkerReads": float(row["Reads_aligned"]),
                "ReassignedReads": float(row["Reads_reassigned"]),
                "PercentIdentity": parse_float(raw.get("Percent_identity", row["PID_aligned"])),
                "MarkerCoveragePercent": parse_float(raw.get("Total_marker_coverage", "NA")),
                "RPKS": parse_float(row["RPKS"]),
                "RPKSB": parse_float(norm.get("Zymo_D6300_RPKSB", "NA")),
                "RelEukPercent": parse_float(frac.get("Relative_abundance", "NA")),
                "Genomes": row["Genomes"],
            }
        )
    if not evidence:
        raise RuntimeError("EukDetect2 produced no species-level evidence")
    write_tsv(summary_dir / "eukdetect-species-evidence.tsv", evidence)

    observed_species = {str(row["Name"]) for row in evidence}
    detection_audit: list[dict[str, object]] = []
    for target, aliases in EXPECTED_GROUP_ALIASES.items():
        matched = sorted(observed_species & aliases)
        detection_audit.append(
            {
                "MockTarget": target,
                "ObservedName": "; ".join(matched) if matched else "Not detected",
                "Expected": "true",
                "DetectedAtSpeciesRank": str(bool(matched)).lower(),
                "EvidenceBoundary": ">=2 marker genes and >=4 reads after EukDetect2 filtering",
            }
        )
    for species in sorted(name for name in observed_species if mock_target(name) is None):
        detection_audit.append(
            {
                "MockTarget": "Unexpected eukaryote",
                "ObservedName": species,
                "Expected": "false",
                "DetectedAtSpeciesRank": "true",
                "EvidenceBoundary": "unexpected species-level call requiring orthogonal review",
            }
        )
    write_tsv(summary_dir / "eukdetect-detection-audit.tsv", detection_audit)

    composition: list[dict[str, object]] = [
        {"Denominator": "Expected total DNA", "Component": "Saccharomyces cerevisiae", "Percent": 2.0},
        {"Denominator": "Expected total DNA", "Component": "Cryptococcus neoformans complex", "Percent": 2.0},
        {"Denominator": "Expected total DNA", "Component": "Other organisms", "Percent": 96.0},
        {"Denominator": "Expected eukaryote-only", "Component": "Saccharomyces cerevisiae", "Percent": 50.0},
        {"Denominator": "Expected eukaryote-only", "Component": "Cryptococcus neoformans complex", "Percent": 50.0},
        {"Denominator": "Expected eukaryote-only", "Component": "Other organisms", "Percent": 0.0},
    ]
    observed_total = sum(float(row["RelEukPercent"]) for row in evidence if not math.isnan(float(row["RelEukPercent"])))
    for row in evidence:
        component = str(row["MockTarget"]) if row["ExpectedInMock"] == "true" else "Other eukaryotes"
        composition.append(
            {"Denominator": "EukDetect RelEuk", "Component": component, "Percent": row["RelEukPercent"]}
        )
    if observed_total < 99.999:
        composition.append(
            {"Denominator": "EukDetect RelEuk", "Component": "Other eukaryotes", "Percent": max(0.0, 100.0 - observed_total)}
        )
    write_tsv(summary_dir / "eukdetect-composition.tsv", composition)

    fragment_ledger = read_tsv(work / "eukrep-fragment-ledger.tsv")
    truth = {row["FragmentID"]: row for row in fragment_ledger}
    benchmark_calls: list[dict[str, object]] = []
    confusion_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    species_rows: list[dict[str, object]] = []
    for fragment_length in FRAGMENT_LENGTHS:
        expected_ids = {row["FragmentID"] for row in fragment_ledger if int(row["FragmentLength"]) == fragment_length}
        for mode in MODES:
            euk_names = read_names(results / "eukrep-reference" / f"length-{fragment_length}.{mode}.euk.names")
            prok_names = read_names(results / "eukrep-reference" / f"length-{fragment_length}.{mode}.prok.names")
            if euk_names & prok_names:
                raise RuntimeError(f"Overlapping EukRep calls for {fragment_length}, {mode}")
            if not (euk_names | prok_names) <= expected_ids:
                raise RuntimeError(f"Unknown EukRep fragment identifier for {fragment_length}, {mode}")
            counts: Counter[tuple[str, str]] = Counter()
            species_counts: Counter[tuple[str, str]] = Counter()
            for fragment_id in sorted(expected_ids):
                prediction = "Eukaryote" if fragment_id in euk_names else "Prokaryote" if fragment_id in prok_names else "Unclassified"
                row = truth[fragment_id]
                counts[(row["DomainTruth"], prediction)] += 1
                species_counts[(row["Species"], prediction)] += 1
                benchmark_calls.append(
                    {
                        "FragmentID": fragment_id,
                        "FragmentLength": fragment_length,
                        "Mode": mode.capitalize(),
                        "Species": row["Species"],
                        "Truth": row["DomainTruth"],
                        "Prediction": prediction,
                        "Correct": str(prediction == row["DomainTruth"]).lower(),
                    }
                )
            for truth_label in ("Eukaryote", "Prokaryote"):
                for prediction in ("Eukaryote", "Prokaryote", "Unclassified"):
                    confusion_rows.append(
                        {
                            "FragmentLength": fragment_length,
                            "Mode": mode.capitalize(),
                            "Truth": truth_label,
                            "Prediction": prediction,
                            "Count": counts[(truth_label, prediction)],
                        }
                    )
            tp = counts[("Eukaryote", "Eukaryote")]
            fn = counts[("Eukaryote", "Prokaryote")] + counts[("Eukaryote", "Unclassified")]
            fp = counts[("Prokaryote", "Eukaryote")]
            tn = counts[("Prokaryote", "Prokaryote")] + counts[("Prokaryote", "Unclassified")]
            classified = sum(counts[(truth_label, prediction)] for truth_label in ("Eukaryote", "Prokaryote") for prediction in ("Eukaryote", "Prokaryote"))
            total = len(expected_ids)
            metrics = {
                "Sensitivity": safe_div(tp, tp + fn),
                "Specificity": safe_div(tn, tn + fp),
                "Precision": safe_div(tp, tp + fp),
                "Classification rate": safe_div(classified, total),
            }
            for metric, value in metrics.items():
                metric_rows.append(
                    {
                        "FragmentLength": fragment_length,
                        "Mode": mode.capitalize(),
                        "Metric": metric,
                        "Estimate": value,
                        "Percent": 100 * value,
                        "Numerator": {"Sensitivity": tp, "Specificity": tn, "Precision": tp, "Classification rate": classified}[metric],
                        "Denominator": {"Sensitivity": tp + fn, "Specificity": tn + fp, "Precision": tp + fp, "Classification rate": total}[metric],
                    }
                )
            for species in sorted({row["Species"] for row in fragment_ledger}):
                domain = next(row["DomainTruth"] for row in fragment_ledger if row["Species"] == species)
                denom = sum(species_counts[(species, prediction)] for prediction in ("Eukaryote", "Prokaryote", "Unclassified"))
                species_rows.append(
                    {
                        "FragmentLength": fragment_length,
                        "Mode": mode.capitalize(),
                        "Species": species,
                        "Truth": domain,
                        "PredictedEukaryote": species_counts[(species, "Eukaryote")],
                        "PredictedProkaryote": species_counts[(species, "Prokaryote")],
                        "Unclassified": species_counts[(species, "Unclassified")],
                        "CorrectPercent": 100 * safe_div(species_counts[(species, domain)], denom),
                    }
                )
    write_tsv(summary_dir / "eukrep-reference-calls.tsv", benchmark_calls)
    write_tsv(summary_dir / "eukrep-reference-confusion.tsv", confusion_rows)
    write_tsv(summary_dir / "eukrep-reference-metrics.tsv", metric_rows)
    write_tsv(summary_dir / "eukrep-reference-species.tsv", species_rows)

    contig_lengths = {header: len(sequence) for header, sequence in parse_fasta(results / "megahit-contigs.oneline.fna")}
    lengths_desc = sorted(contig_lengths.values(), reverse=True)
    total_length = sum(lengths_desc)
    running = 0
    n50 = 0
    for length in lengths_desc:
        running += length
        if running >= total_length / 2:
            n50 = length
            break
    best_alignment: dict[str, dict[str, object]] = {}
    paf_path = results / "megahit-to-zymo-v2.paf"
    with paf_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query, qlen = fields[0], int(fields[1])
            qstart, qend = int(fields[2]), int(fields[3])
            target = fields[5]
            matches, block, mapq = int(fields[9]), int(fields[10]), int(fields[11])
            candidate = {
                "Target": target,
                "Species": target.split("|")[1].replace("_", " "),
                "Identity": matches / block if block else 0.0,
                # PAF alignment block length includes indels and can exceed the
                # query span. Coverage is the aligned query interval instead.
                "QueryCoverage": (qend - qstart) / qlen if qlen else 0.0,
                "Matches": matches,
                "MapQ": mapq,
            }
            previous = best_alignment.get(query)
            rank = (candidate["Matches"], candidate["Identity"], candidate["QueryCoverage"], candidate["MapQ"])
            prev_rank = (-1, -1.0, -1.0, -1) if previous is None else (previous["Matches"], previous["Identity"], previous["QueryCoverage"], previous["MapQ"])
            if rank > prev_rank:
                best_alignment[query] = candidate

    assembly_calls: list[dict[str, object]] = []
    assembly_metrics: list[dict[str, object]] = []
    eligible = {name for name, length in contig_lengths.items() if length >= 3000}
    truth_by_contig: dict[str, tuple[str, str]] = {}
    for contig in eligible:
        hit = best_alignment.get(contig)
        if hit and float(hit["Identity"]) >= 0.95 and float(hit["QueryCoverage"]) >= 0.80:
            truth_domain = "Eukaryote" if hit["Species"] in REFERENCE_EUKARYOTES else "Prokaryote"
            truth_by_contig[contig] = (truth_domain, str(hit["Species"]))
        else:
            truth_by_contig[contig] = ("Unresolved", "Unresolved")
    for mode in MODES:
        euk_names = read_names(results / "eukrep-assembly" / f"megahit.{mode}.euk.names")
        prok_names = read_names(results / "eukrep-assembly" / f"megahit.{mode}.prok.names")
        counts: Counter[tuple[str, str]] = Counter()
        for contig in sorted(eligible):
            prediction = "Eukaryote" if contig in euk_names else "Prokaryote" if contig in prok_names else "Unclassified"
            truth_domain, species = truth_by_contig[contig]
            hit = best_alignment.get(contig, {})
            counts[(truth_domain, prediction)] += 1
            assembly_calls.append(
                {
                    "Contig": contig,
                    "LengthBp": contig_lengths[contig],
                    "Mode": mode.capitalize(),
                    "Truth": truth_domain,
                    "ReferenceSpecies": species,
                    "Prediction": prediction,
                    "IdentityPercent": 100 * float(hit.get("Identity", 0.0)),
                    "QueryCoveragePercent": 100 * float(hit.get("QueryCoverage", 0.0)),
                    "MapQ": int(hit.get("MapQ", 0)),
                }
            )
        tp = counts[("Eukaryote", "Eukaryote")]
        fn = counts[("Eukaryote", "Prokaryote")] + counts[("Eukaryote", "Unclassified")]
        fp = counts[("Prokaryote", "Eukaryote")]
        tn = counts[("Prokaryote", "Prokaryote")] + counts[("Prokaryote", "Unclassified")]
        for metric, numerator, denominator in (
            ("Sensitivity", tp, tp + fn),
            ("Specificity", tn, tn + fp),
            ("Precision", tp, tp + fp),
        ):
            estimate = safe_div(numerator, denominator)
            assembly_metrics.append(
                {
                    "Mode": mode.capitalize(),
                    "Metric": metric,
                    "Estimate": estimate,
                    "Percent": 100 * estimate,
                    "Numerator": numerator,
                    "Denominator": denominator,
                    "ResolvedEukaryoteContigs": sum(counts[("Eukaryote", p)] for p in ("Eukaryote", "Prokaryote", "Unclassified")),
                    "ResolvedProkaryoteContigs": sum(counts[("Prokaryote", p)] for p in ("Eukaryote", "Prokaryote", "Unclassified")),
                    "UnresolvedContigs": sum(counts[("Unresolved", p)] for p in ("Eukaryote", "Prokaryote", "Unclassified")),
                }
            )
    write_tsv(summary_dir / "eukrep-assembly-calls.tsv", assembly_calls)
    write_tsv(summary_dir / "eukrep-assembly-metrics.tsv", assembly_metrics)

    fastp_euk = json.loads((results / "qc" / "fastp-eukdetect.json").read_text(encoding="utf-8"))
    fastp_asm = json.loads((results / "qc" / "fastp-assembly.json").read_text(encoding="utf-8"))
    qc_rows = []
    for branch, payload in (("EukDetect", fastp_euk), ("Assembly", fastp_asm)):
        before = payload["summary"]["before_filtering"]
        after = payload["summary"]["after_filtering"]
        qc_rows.append(
            {
                "Branch": branch,
                "ReadPairsBefore": int(before["total_reads"]) // 2,
                "ReadPairsAfter": int(after["total_reads"]) // 2,
                "BasesBefore": before["total_bases"],
                "BasesAfter": after["total_bases"],
                "ReadRetentionPercent": 100 * after["total_reads"] / before["total_reads"],
                "Q30BeforePercent": 100 * before["q30_rate"],
                "Q30AfterPercent": 100 * after["q30_rate"],
            }
        )
    write_tsv(summary_dir / "fastp-qc-summary.tsv", qc_rows)

    assembly_summary = [
        {
            "RunAccession": "SRR12324253",
            "ContigsGE1000": len(contig_lengths),
            "TotalBases": total_length,
            "N50Bp": n50,
            "ContigsGE3000": len(eligible),
            "TruthResolved": sum(1 for truth_domain, _ in truth_by_contig.values() if truth_domain != "Unresolved"),
            "TruthEukaryote": sum(1 for truth_domain, _ in truth_by_contig.values() if truth_domain == "Eukaryote"),
            "TruthProkaryote": sum(1 for truth_domain, _ in truth_by_contig.values() if truth_domain == "Prokaryote"),
            "TruthUnresolved": sum(1 for truth_domain, _ in truth_by_contig.values() if truth_domain == "Unresolved"),
            "TruthGate": ">=95% identity and >=80% query coverage to official Zymo v2 reference",
        }
    ]
    write_tsv(summary_dir / "assembly-summary.tsv", assembly_summary)

    evidence_ladder = [
        {"Rank": 1, "Evidence": "Marker reads", "ClaimCeiling": "Taxon detected"},
        {"Rank": 2, "Evidence": "Eukaryotic contigs", "ClaimCeiling": "Domain-classified sequence"},
        {"Rank": 3, "Evidence": "Genome bin", "ClaimCeiling": "Coherent genome reconstruction"},
        {"Rank": 4, "Evidence": "RNA, protein or growth", "ClaimCeiling": "Activity or viability supported"},
        {"Rank": 5, "Evidence": "Longitudinal or experimental evidence", "ClaimCeiling": "Persistence or colonization supported"},
    ]
    write_tsv(summary_dir / "evidence-ladder.tsv", evidence_ladder)

    resource_rows: list[dict[str, object]] = []
    command_rows = read_tsv(work / "command-log.tsv")
    for row in command_rows:
        timing = Path(row["Timing"])
        if not timing.is_absolute():
            timing = root / timing
        if not timing.is_file():
            raise FileNotFoundError(timing)
        parsed = parse_time_file(timing)
        resource_rows.append(
            {
                "Label": row["Label"],
                "WallSeconds": round(parsed["WallSeconds"], 3),
                "CPUSeconds": round(parsed["CPUSeconds"], 3),
                "PeakRSSGiB": round(parsed["PeakRSSGiB"], 4),
                "Command": row["Command"],
            }
        )
    write_tsv(summary_dir / "resource-usage.tsv", resource_rows)

    summary = {
        "run_accession": "SRR12324253",
        "eukdetect_species_calls": len(evidence),
        "eukdetect_expected_species_detected": sum(
            bool(observed_species & aliases) for aliases in EXPECTED_GROUP_ALIASES.values()
        ),
        "eukdetect_unexpected_species_calls": sum(
            mock_target(name) is None for name in observed_species
        ),
        "reference_benchmark_fragments": len(fragment_ledger),
        "assembly_contigs_ge1000": len(contig_lengths),
        "assembly_total_bases": total_length,
        "assembly_n50": n50,
        "assembly_contigs_ge3000": len(eligible),
        "assembly_truth_resolved": sum(1 for truth_domain, _ in truth_by_contig.values() if truth_domain != "Unresolved"),
        "assembly_truth_eukaryote": sum(1 for truth_domain, _ in truth_by_contig.values() if truth_domain == "Eukaryote"),
        "assembly_truth_prokaryote": sum(1 for truth_domain, _ in truth_by_contig.values() if truth_domain == "Prokaryote"),
        "assembly_truth_unresolved": sum(1 for truth_domain, _ in truth_by_contig.values() if truth_domain == "Unresolved"),
        "measured_commands": len(resource_rows),
    }
    (summary_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (work / ".article58-summary-complete").write_text("verified\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
