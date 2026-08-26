#!/usr/bin/env python3
"""Summarize GTDB-Tk R232 taxonomy without using mock reference truth."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

from article41_44_utils import dump_json, parse_time, read_tsv, sha256, write_tsv


RANKS = (
    ("Domain", "d__"), ("Phylum", "p__"), ("Class", "c__"),
    ("Order", "o__"), ("Family", "f__"), ("Genus", "g__"), ("Species", "s__"),
)


def value(row: dict[str, str], key: str, default: str = "") -> str:
    candidate = row.get(key, default).strip()
    return "" if candidate in {"N/A", "NA", "none", "None"} else candidate


def first_value(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        candidate = value(row, key)
        if candidate:
            return candidate
    return ""


def optional_float(row: dict[str, str], key: str, multiplier: float = 1.0) -> float | str:
    candidate = value(row, key)
    if not candidate:
        return ""
    return float(candidate) * multiplier


def split_taxonomy(taxonomy: str) -> dict[str, str]:
    pieces = taxonomy.split(";")
    mapping = {piece[:3]: piece for piece in pieces if len(piece) >= 3 and piece[1:3] == "__"}
    return {rank: mapping.get(prefix, prefix) for rank, prefix in RANKS}


def assigned(taxon: str) -> bool:
    return len(taxon) > 3 and bool(taxon[3:].strip())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article46-run-complete").is_file():
        raise FileNotFoundError("Run run_article46_gtdbtk.py first")
    output = work / "gtdbtk"
    summary_files = sorted(output.glob("classify/article46.*.summary.tsv"))
    if not summary_files:
        raise FileNotFoundError("No GTDB-Tk summary tables found")
    native: list[dict[str, str]] = []
    for path in summary_files:
        with path.open(encoding="utf-8", newline="") as handle:
            native.extend(csv.DictReader(handle, delimiter="\t"))
    ledger = {row["SGB"]: row for row in read_tsv(work / "genome-ledger.tsv")}
    if {row["user_genome"] for row in native} != set(ledger):
        raise ValueError("GTDB-Tk output does not cover all 24 SGBs")

    taxonomy_rows: list[dict[str, object]] = []
    for row in sorted(native, key=lambda item: item["user_genome"]):
        sgb = row["user_genome"]
        taxonomy = value(row, "classification")
        ranks = split_taxonomy(taxonomy)
        species_assigned = assigned(ranks["Species"])
        fastani_reference = first_value(row, "closest_genome_reference", "fastani_reference")
        fastani_radius = first_value(row, "closest_genome_reference_radius", "fastani_reference_radius")
        fastani_taxonomy = first_value(row, "closest_genome_taxonomy", "fastani_taxonomy")
        fastani_ani = first_value(row, "closest_genome_ani", "fastani_ani")
        fastani_af_raw = first_value(row, "closest_genome_af", "fastani_af")
        fastani_af = float(fastani_af_raw) * 100 if fastani_af_raw else ""
        closest_af = optional_float(row, "closest_placement_af", 100)
        result: dict[str, object] = {
            "SGB": sgb,
            "Representative": ledger[sgb]["Representative"],
            "RepresentativeSHA256": ledger[sgb]["RepresentativeSHA256"],
            "GTDBRelease": "R232",
            "GTDBTaxonomy": taxonomy,
            **ranks,
            "SpeciesAssigned": species_assigned,
            "GTDBSpeciesUnresolvedCandidate": not species_assigned,
            "ClassificationMethod": value(row, "classification_method"),
            "FastANIReference": fastani_reference,
            "FastANIReferenceRadiusPct": float(fastani_radius) if fastani_radius else "",
            "FastANIReferenceTaxonomy": fastani_taxonomy,
            "FastANIANI": float(fastani_ani) if fastani_ani else "",
            "FastANIAFPct": fastani_af,
            "ClosestPlacementReference": value(row, "closest_placement_reference"),
            "ClosestPlacementRadiusPct": optional_float(row, "closest_placement_radius"),
            "ClosestPlacementTaxonomy": value(row, "closest_placement_taxonomy"),
            "ClosestPlacementANI": optional_float(row, "closest_placement_ani"),
            "ClosestPlacementAFPct": closest_af,
            "PplacerTaxonomy": value(row, "pplacer_taxonomy"),
            "MSAPercent": optional_float(row, "msa_percent"),
            "TranslationTable": value(row, "translation_table"),
            "REDValue": value(row, "red_value"),
            "Note": value(row, "note"),
            "Warnings": value(row, "warnings"),
            "Completeness": float(ledger[sgb]["Completeness"]),
            "Contamination": float(ledger[sgb]["Contamination"]),
            "MIMAGQuality": ledger[sgb]["MIMAGQuality"],
        }
        taxonomy_rows.append(result)

    rank_rows: list[dict[str, object]] = []
    for rank, _ in RANKS:
        for domain in sorted({str(row["Domain"]) for row in taxonomy_rows}):
            subset = [row for row in taxonomy_rows if row["Domain"] == domain]
            resolved = sum(assigned(str(row[rank])) for row in subset)
            rank_rows.append({
                "Domain": domain,
                "Rank": rank,
                "ResolvedSGBs": resolved,
                "UnresolvedSGBs": len(subset) - resolved,
                "SGBs": len(subset),
            })
    phyla = Counter((str(row["Domain"]), str(row["Phylum"])) for row in taxonomy_rows)
    phylum_rows = [
        {"Domain": domain, "Phylum": phylum, "SGBs": count}
        for (domain, phylum), count in sorted(phyla.items())
    ]
    methods = Counter(str(row["ClassificationMethod"]) for row in taxonomy_rows)
    method_rows = [
        {"ClassificationMethod": method, "SGBs": count}
        for method, count in sorted(methods.items())
    ]
    ani_rows = [
        {
            "SGB": row["SGB"], "Domain": row["Domain"],
            "SpeciesAssigned": row["SpeciesAssigned"],
            "FastANIReference": row["FastANIReference"],
            "FastANIReferenceRadiusPct": row["FastANIReferenceRadiusPct"],
            "FastANIANI": row["FastANIANI"], "FastANIAFPct": row["FastANIAFPct"],
            "ANIMarginToReferenceRadiusPctPoints": (
                float(row["FastANIANI"]) - float(row["FastANIReferenceRadiusPct"])
                if row["FastANIANI"] != "" and row["FastANIReferenceRadiusPct"] != "" else ""
            ),
        }
        for row in taxonomy_rows if row["FastANIReference"]
    ]
    marker_files: list[dict[str, object]] = []
    for path in sorted(output.rglob("*")):
        relative = path.relative_to(output).as_posix()
        if path.is_file() and ("single_copy" in relative or relative.endswith(".msa.fasta.gz")):
            marker_files.append({
                "Path": relative,
                "Bytes": path.stat().st_size,
                "SHA256": sha256(path),
            })

    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(summary_dir / "taxonomy-summary.tsv", taxonomy_rows)
    write_tsv(summary_dir / "rank-resolution.tsv", rank_rows)
    write_tsv(summary_dir / "phylum-summary.tsv", phylum_rows)
    write_tsv(summary_dir / "classification-method-summary.tsv", method_rows)
    write_tsv(
        summary_dir / "fastani-reference-audit.tsv", ani_rows,
        fieldnames=["SGB", "Domain", "SpeciesAssigned", "FastANIReference", "FastANIReferenceRadiusPct", "FastANIANI", "FastANIAFPct", "ANIMarginToReferenceRadiusPctPoints"],
    )
    write_tsv(
        summary_dir / "marker-file-inventory.tsv", marker_files,
        fieldnames=["Path", "Bytes", "SHA256"],
    )
    resources = [
        parse_time(work / "logs/gtdbtk-check-install.time.txt"),
        parse_time(work / "logs/gtdbtk-classify-wf.time.txt"),
    ]
    write_tsv(summary_dir / "resource-summary.tsv", resources)
    domain_counts = Counter(str(row["Domain"]).removeprefix("d__") for row in taxonomy_rows)
    result = {
        "article": 46,
        "input_sgbs": 24,
        "classified_sgbs": len(taxonomy_rows),
        "domain_counts": dict(sorted(domain_counts.items())),
        "species_assigned": sum(bool(row["SpeciesAssigned"]) for row in taxonomy_rows),
        "species_unresolved_candidates": sum(bool(row["GTDBSpeciesUnresolvedCandidate"]) for row in taxonomy_rows),
        "classification_methods": dict(sorted(methods.items())),
        "warnings_nonempty": sum(bool(row["Warnings"]) for row in taxonomy_rows),
        "fastani_reference_hits": len(ani_rows),
        "single_copy_or_msa_files": len(marker_files),
        "gtdb_release": "R232",
        "gtdbtk_version": "2.7.2",
        "truth_used_for_taxonomy": False,
    }
    dump_json(summary_dir / "run-summary.json", result)
    (work / ".article46-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
