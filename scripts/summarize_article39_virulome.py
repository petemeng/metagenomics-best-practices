#!/usr/bin/env python3
"""Summarize ABRicate/VFDB thresholds, database sensitivity, categories, and controls."""

from __future__ import annotations

import argparse
import collections
import re
from pathlib import Path

from article37_40_utils import dump_json, parse_time, read_abundance, read_metadata, read_tsv, write_tsv


CATALOG_GENES = 93_782
SAMPLES = ("MOCK1", "MOCK2")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def category(product: str) -> tuple[str, str]:
    factor = re.search(r"\[([^\[]+?) \(VF\d+\) - ([^\[]+?) \(VFC\d+\)\]", product)
    if not factor:
        return "Unparsed factor", "Unparsed category"
    return factor.group(1).strip(), factor.group(2).strip()


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    if not (work / ".article39-run-complete").is_file():
        raise FileNotFoundError("Article 39 run is incomplete")
    out = work / "summary"
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "Core 90/80": work / "catalog-core-primary/hits.tsv",
        "Core 80/80": work / "catalog-core-sensitive/hits.tsv",
        "Full 90/80": work / "catalog-full-primary/hits.tsv",
        "Co-assembly core 90/80": work / "coassembly/hits.tsv",
        "Pseudomonas core 90/80": work / "pseudomonas/hits.tsv",
        "Staphylococcus core 90/80": work / "staphylococcus/hits.tsv",
    }
    datasets = {label: read_tsv(path) for label, path in paths.items()}
    metadata = read_metadata(root)
    abundance, sample_totals = read_abundance(root)
    sets = {label: {row["SEQUENCE"] for row in rows} for label, rows in datasets.items()}
    primary_hits = {row["SEQUENCE"]: row for row in datasets["Core 90/80"]}
    if len(metadata) != CATALOG_GENES or len(primary_hits) != len(datasets["Core 90/80"]):
        raise ValueError("Catalog identity or one-hit-per-gene contract failed")

    gene_rows: list[dict[str, object]] = []
    for meta in metadata:
        gene = meta["RepresentativeID"]
        hit = primary_hits.get(gene)
        factor, vfc = category(hit["PRODUCT"]) if hit else ("-", "-")
        gene_rows.append({
            "GeneID": gene, "Completeness": meta["Completeness"], "NtLength": meta["NtLength"],
            "CorePrimary": "yes" if gene in sets["Core 90/80"] else "no",
            "CoreSensitive": "yes" if gene in sets["Core 80/80"] else "no",
            "FullPrimary": "yes" if gene in sets["Full 90/80"] else "no",
            "Gene": hit["GENE"] if hit else "-", "Accession": hit["ACCESSION"] if hit else "-",
            "VirulenceFactor": factor, "VFCCategory": vfc,
            "ReferenceCoveragePercent": hit["%COVERAGE"] if hit else "", "IdentityPercent": hit["%IDENTITY"] if hit else "",
            "MOCK1RawReads": abundance[gene]["MOCK1"], "MOCK2RawReads": abundance[gene]["MOCK2"],
        })
    write_tsv(out / "virulome-gene-calls.tsv.gz", gene_rows)

    sensitivity_rows = []
    for label in ("Core 90/80", "Core 80/80", "Full 90/80"):
        current = sets[label]
        sensitivity_rows.append({
            "Branch": label, "Hits": len(datasets[label]), "UniqueGenes": len(current),
            "SharedWithCorePrimary": len(current & sets["Core 90/80"]), "AdditionalVsCorePrimary": len(current - sets["Core 90/80"]),
        })
    write_tsv(out / "sensitivity-summary.tsv", sensitivity_rows)
    hit_rows = []
    for label, rows in datasets.items():
        for row in rows:
            factor, vfc = category(row["PRODUCT"])
            hit_rows.append({
                "Branch": label, "Sequence": row["SEQUENCE"], "Gene": row["GENE"], "Accession": row["ACCESSION"],
                "VirulenceFactor": factor, "VFCCategory": vfc, "ReferenceCoveragePercent": row["%COVERAGE"],
                "IdentityPercent": row["%IDENTITY"], "AlignedBp": int(row["END"]) - int(row["START"]) + 1,
            })
    write_tsv(out / "all-hit-quality.tsv", hit_rows)

    primary = [row for row in gene_rows if row["CorePrimary"] == "yes"]
    accumulator: dict[str, dict[str, float]] = collections.defaultdict(lambda: collections.defaultdict(float))
    for row in primary:
        label = str(row["VFCCategory"])
        accumulator[label]["Genes"] += 1
        for sample in SAMPLES:
            accumulator[label][f"{sample}RawReads"] += int(row[f"{sample}RawReads"])
    category_rows = [{"VFCCategory": label, **{k: round(v, 8) for k, v in values.items()}} for label, values in accumulator.items()]
    category_rows.sort(key=lambda row: (-float(row["Genes"]), row["VFCCategory"]))
    write_tsv(out / "vfc-category-summary.tsv", category_rows)

    control_rows = [row for row in hit_rows if row["Branch"] in {"Pseudomonas core 90/80", "Staphylococcus core 90/80"}]
    write_tsv(out / "positive-control-hits.tsv", control_rows)
    context_rows = []
    for label in ("Core 90/80", "Co-assembly core 90/80", "Pseudomonas core 90/80", "Staphylococcus core 90/80"):
        rows = datasets[label]
        context_rows.append({"Branch": label, "Hits": len(rows), "SequencesWithHits": len({row["SEQUENCE"] for row in rows}), "VFCategories": len({category(row["PRODUCT"])[1] for row in rows})})
    write_tsv(out / "context-summary.tsv", context_rows)
    write_tsv(out / "resource-usage.tsv", [parse_time(path) for path in sorted((work / "logs").glob("*.time.txt"))])
    dump_json(out / "run-summary.json", {
        "article": 39, "catalog_genes": CATALOG_GENES, "catalog_core_primary_genes": len(primary),
        "catalog_core_sensitive_genes": len(sets["Core 80/80"]), "catalog_full_primary_genes": len(sets["Full 90/80"]),
        "coassembly_core_primary_hits": len(datasets["Co-assembly core 90/80"]),
        "pseudomonas_core_primary_hits": len(datasets["Pseudomonas core 90/80"]),
        "staphylococcus_core_primary_hits": len(datasets["Staphylococcus core 90/80"]),
        "catalog_primary_raw_reads": {sample: sum(int(row[f"{sample}RawReads"]) for row in primary) for sample in SAMPLES},
        "sample_assigned_reads": sample_totals, "primary_database": "VFDB core set A 2026-07-24",
        "primary_thresholds": {"identity_percent": 90, "reference_coverage_percent": 80},
        "vfdb_match_is_not_virulence_phenotype": True, "strain_context_and_regulation_required": True,
    })
    (work / ".article39-summary-complete").write_text("complete\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
