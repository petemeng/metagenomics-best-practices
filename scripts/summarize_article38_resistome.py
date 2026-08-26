#!/usr/bin/env python3
"""Summarize RGI evidence tiers, CARD models, drug classes, and read weighting."""

from __future__ import annotations

import argparse
import collections
from pathlib import Path

from article37_40_utils import dump_json, parse_time, read_abundance, read_metadata, read_tsv, write_tsv


CATALOG_GENES = 93_782
SAMPLES = ("MOCK1", "MOCK2")
PRIMARY = {"Perfect", "Strict"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def split_labels(value: str) -> list[str]:
    return list(dict.fromkeys(x.strip() for x in value.split(";") if x.strip() and x.strip().lower() != "n/a"))


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    if not (work / ".article38-run-complete").is_file():
        raise FileNotFoundError("Article 38 run is incomplete")
    out = work / "summary"
    out.mkdir(parents=True, exist_ok=True)
    tool_rows = read_tsv(work / "tool-versions.tsv")
    for row in tool_rows:
        if row["Tool"] == "RGI" and "6.0.8" in row["VersionEvidence"]:
            row["VersionEvidence"] = "6.0.8"
    write_tsv(work / "tool-versions.tsv", tool_rows)
    datasets = {name: read_tsv(work / name / f"{name}.txt") for name in ("catalog", "coassembly", "pseudomonas", "staphylococcus")}
    metadata = read_metadata(root)
    abundance, sample_totals = read_abundance(root)
    catalog_hits = {row["ORF_ID"].split()[0]: row for row in datasets["catalog"]}
    if len(metadata) != CATALOG_GENES or len(catalog_hits) != len(datasets["catalog"]):
        raise ValueError("Catalog identity or one-hit-per-protein contract failed")

    gene_rows: list[dict[str, object]] = []
    for meta in metadata:
        gene = meta["RepresentativeID"]
        hit = catalog_hits.get(gene)
        tier = hit["Cut_Off"] if hit else "No hit"
        drug_classes = split_labels(hit["Drug Class"]) if hit else []
        mechanisms = split_labels(hit["Resistance Mechanism"]) if hit else []
        families = split_labels(hit["AMR Gene Family"]) if hit else []
        gene_rows.append({
            "GeneID": gene, "Completeness": meta["Completeness"], "AaLength": meta["AaLength"], "EvidenceTier": tier,
            "PrimaryARG": "yes" if tier in PRIMARY else "no", "BestHitARO": hit["Best_Hit_ARO"] if hit else "-",
            "ARO": hit["ARO"] if hit else "-", "ModelType": hit["Model_type"] if hit else "-",
            "IdentityPercent": hit["Best_Identities"] if hit else "", "ReferenceLengthPercent": hit["Percentage Length of Reference Sequence"] if hit else "",
            "DrugClasses": ";".join(drug_classes) or "-", "ResistanceMechanisms": ";".join(mechanisms) or "-",
            "AMRGeneFamilies": ";".join(families) or "-", "Nudged": hit["Nudged"] if hit else "",
            "MOCK1RawReads": abundance[gene]["MOCK1"], "MOCK2RawReads": abundance[gene]["MOCK2"],
        })
    write_tsv(out / "resistome-gene-calls.tsv.gz", gene_rows)

    tier_rows: list[dict[str, object]] = []
    for dataset, rows in datasets.items():
        for tier in ("Perfect", "Strict", "Loose"):
            selected = [row for row in rows if row["Cut_Off"] == tier]
            tier_rows.append({"Dataset": dataset, "EvidenceTier": tier, "Hits": len(selected), "UniqueORFs": len({row["ORF_ID"].split()[0] for row in selected})})
    write_tsv(out / "evidence-tier-summary.tsv", tier_rows)

    primary = [row for row in gene_rows if row["PrimaryARG"] == "yes"]
    def allocate(field: str, label_name: str) -> list[dict[str, object]]:
        accumulator: dict[str, dict[str, float]] = collections.defaultdict(lambda: collections.defaultdict(float))
        for row in primary:
            labels = [] if row[field] == "-" else str(row[field]).split(";")
            for label in labels:
                accumulator[label]["GenesWithLabel"] += 1
                accumulator[label]["FractionalGeneEquivalent"] += 1 / len(labels)
                for sample in SAMPLES:
                    reads = int(row[f"{sample}RawReads"])
                    accumulator[label][f"{sample}ReadsFromGenesWithLabel"] += reads
                    accumulator[label][f"{sample}FractionalRawReads"] += reads / len(labels)
        result = [{label_name: label, **{k: round(v, 8) for k, v in values.items()}} for label, values in accumulator.items()]
        return sorted(result, key=lambda row: (-float(row["FractionalGeneEquivalent"]), str(row[label_name])))

    write_tsv(out / "drug-class-summary.tsv", allocate("DrugClasses", "DrugClass"))
    write_tsv(out / "mechanism-summary.tsv", allocate("ResistanceMechanisms", "ResistanceMechanism"))
    write_tsv(out / "amr-family-summary.tsv", allocate("AMRGeneFamilies", "AMRGeneFamily"))
    model_rows = []
    for model, count in collections.Counter(row["ModelType"] for row in primary).most_common():
        model_rows.append({"ModelType": model, "PrimaryGenes": count, "GenePercent": 100 * count / len(primary)})
    write_tsv(out / "model-type-summary.tsv", model_rows)

    quality_rows = []
    for row in primary:
        quality_rows.append({
            "GeneID": row["GeneID"], "EvidenceTier": row["EvidenceTier"], "BestHitARO": row["BestHitARO"],
            "IdentityPercent": row["IdentityPercent"], "ReferenceLengthPercent": row["ReferenceLengthPercent"], "ModelType": row["ModelType"],
        })
    write_tsv(out / "primary-hit-quality.tsv", quality_rows)
    control_rows = []
    for dataset in ("pseudomonas", "staphylococcus"):
        for row in datasets[dataset]:
            if row["Cut_Off"] in PRIMARY:
                control_rows.append({"Control": dataset, "EvidenceTier": row["Cut_Off"], "BestHitARO": row["Best_Hit_ARO"], "ARO": row["ARO"], "IdentityPercent": row["Best_Identities"], "ReferenceLengthPercent": row["Percentage Length of Reference Sequence"], "ModelType": row["Model_type"]})
    write_tsv(out / "positive-control-hits.tsv", control_rows)
    write_tsv(out / "resource-usage.tsv", [parse_time(path) for path in sorted((work / "logs").glob("*.time.txt"))])
    primary_read_totals = {sample: sum(int(row[f"{sample}RawReads"]) for row in primary) for sample in SAMPLES}
    dump_json(out / "run-summary.json", {
        "article": 38, "catalog_genes": CATALOG_GENES, "catalog_hits_all_tiers": len(catalog_hits),
        "catalog_primary_genes": len(primary), "catalog_loose_only_genes": sum(row["EvidenceTier"] == "Loose" for row in gene_rows),
        "coassembly_primary_orfs": len({row["ORF_ID"].split()[0] for row in datasets["coassembly"] if row["Cut_Off"] in PRIMARY}),
        "pseudomonas_primary_orfs": len({row["ORF_ID"].split()[0] for row in datasets["pseudomonas"] if row["Cut_Off"] in PRIMARY}),
        "staphylococcus_primary_orfs": len({row["ORF_ID"].split()[0] for row in datasets["staphylococcus"] if row["Cut_Off"] in PRIMARY}),
        "catalog_primary_raw_reads": primary_read_totals, "sample_assigned_reads": sample_totals,
        "primary_tiers": ["Perfect", "Strict"], "loose_is_sensitivity_only": True, "include_nudge": False,
        "arg_presence_is_not_phenotypic_resistance": True, "protein_catalog_captures_rrna_variants": False,
    })
    (work / ".article38-summary-complete").write_text("complete\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
