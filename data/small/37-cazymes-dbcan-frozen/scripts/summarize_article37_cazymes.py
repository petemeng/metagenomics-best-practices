#!/usr/bin/env python3
"""Summarize dbCAN consensus, CAZy classes, substrates, CGCs, and read weighting."""

from __future__ import annotations

import argparse
import collections
import re
from pathlib import Path

from article37_40_utils import dump_json, parse_time, read_abundance, read_metadata, read_tsv, write_tsv


CATALOG_GENES = 93_782
SAMPLES = ("MOCK1", "MOCK2")
CLASSES = ("GH", "GT", "PL", "CE", "AA", "CBM", "Other module")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def tokens(value: str) -> list[str]:
    if not value or value == "-":
        return []
    found: list[str] = []
    for token in re.split(r"[|;+]", value):
        match = re.search(r"((?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?)", token)
        if match:
            family = re.sub(r"_e\d+$", "", match.group(1))
            if family not in found:
                found.append(family)
        elif token.strip().lower() in {"cohesin", "dockerin", "slh"}:
            family = token.strip().capitalize()
            if family not in found:
                found.append(family)
    return found


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    if not (work / ".article37-run-complete").is_file():
        raise FileNotFoundError("Article 37 run is incomplete")
    summary = work / "summary"
    summary.mkdir(parents=True, exist_ok=True)
    overview = {row["Gene ID"]: row for row in read_tsv(work / "catalog/overview.tsv")}
    metadata = read_metadata(root)
    abundance, sample_totals = read_abundance(root)
    if len(metadata) != CATALOG_GENES or sample_totals != {"MOCK1": 2_784_234, "MOCK2": 2_777_443}:
        raise ValueError("Catalog or abundance identity changed")

    gene_rows: list[dict[str, object]] = []
    for meta in metadata:
        gene = meta["RepresentativeID"]
        call = overview.get(gene, {})
        number = int(call.get("#ofTools", "0") or 0)
        families = tokens(call.get("Recommend Results", "-")) if number >= 2 else []
        classes = list(dict.fromkeys(
            re.match(r"[A-Z]+", family).group(0) if re.match(r"(?:GH|GT|PL|CE|AA|CBM)\d", family) else "Other module"
            for family in families
        ))
        substrates = [x.strip() for x in call.get("Substrate", "-").split(";") if x.strip() and x.strip() != "-"] if number >= 2 else []
        gene_rows.append({
            "GeneID": gene, "Completeness": meta["Completeness"], "AaLength": meta["AaLength"],
            "EvidenceTools": number, "EvidenceTier": "Primary consensus" if number >= 2 else ("Single-tool sensitivity" if number == 1 else "No CAZyme evidence"),
            "dbCANHMM": call.get("dbCAN_hmm", "-"), "dbCANsub": call.get("dbCAN_sub", "-"), "DIAMOND": call.get("DIAMOND", "-"),
            "RecommendedFamilies": ";".join(families) or "-", "CAZyClasses": ";".join(classes) or "-",
            "PredictedSubstrates": ";".join(substrates) or "-", "PrimaryCAZyme": "yes" if number >= 2 else "no",
            "MOCK1RawReads": abundance[gene]["MOCK1"], "MOCK2RawReads": abundance[gene]["MOCK2"],
        })
    write_tsv(summary / "cazyme-gene-calls.tsv.gz", gene_rows)

    evidence_rows = []
    for number, label in ((0, "No CAZyme evidence"), (1, "Single-tool sensitivity"), (2, "Two-tool consensus"), (3, "Three-tool consensus")):
        selected = [row for row in gene_rows if int(row["EvidenceTools"]) == number]
        row: dict[str, object] = {"EvidenceTools": number, "EvidenceTier": label, "Genes": len(selected), "GenePercent": 100 * len(selected) / CATALOG_GENES}
        for sample in SAMPLES:
            reads = sum(int(x[f"{sample}RawReads"]) for x in selected)
            row[f"{sample}RawReads"] = reads
            row[f"{sample}ReadPercent"] = 100 * reads / sample_totals[sample]
        evidence_rows.append(row)
    write_tsv(summary / "evidence-tier-summary.tsv", evidence_rows)

    primary = [row for row in gene_rows if row["PrimaryCAZyme"] == "yes"]
    class_acc: dict[str, dict[str, float]] = collections.defaultdict(lambda: collections.defaultdict(float))
    family_acc: dict[str, dict[str, float]] = collections.defaultdict(lambda: collections.defaultdict(float))
    substrate_acc: dict[str, dict[str, float]] = collections.defaultdict(lambda: collections.defaultdict(float))
    for row in primary:
        families = row["RecommendedFamilies"].split(";")
        classes = row["CAZyClasses"].split(";")
        substrates = [] if row["PredictedSubstrates"] == "-" else row["PredictedSubstrates"].split(";")
        for collection, labels in ((family_acc, families), (class_acc, classes), (substrate_acc, substrates)):
            if not labels:
                continue
            for label in labels:
                collection[label]["GenesWithLabel"] += 1
                collection[label]["FractionalGeneEquivalent"] += 1 / len(labels)
                for sample in SAMPLES:
                    reads = int(row[f"{sample}RawReads"])
                    collection[label][f"{sample}ReadsFromGenesWithLabel"] += reads
                    collection[label][f"{sample}FractionalRawReads"] += reads / len(labels)

    def allocation_rows(accumulator: dict[str, dict[str, float]], label_field: str) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for label, values in accumulator.items():
            rows.append({label_field: label, **{key: round(value, 8) for key, value in values.items()}})
        return sorted(rows, key=lambda row: (-float(row["FractionalGeneEquivalent"]), str(row[label_field])))

    class_rows = allocation_rows(class_acc, "CAZyClass")
    if {row["CAZyClass"] for row in class_rows} - set(CLASSES):
        raise ValueError("Unexpected CAZy class")
    write_tsv(summary / "cazyme-class-summary.tsv", class_rows)
    write_tsv(summary / "cazyme-family-summary.tsv", allocation_rows(family_acc, "CAZyFamily"))
    write_tsv(summary / "substrate-summary.tsv", allocation_rows(substrate_acc, "PredictedSubstrate"))

    completeness_rows = []
    for label in sorted({row["Completeness"] for row in gene_rows}):
        selected = [row for row in gene_rows if row["Completeness"] == label]
        primary_count = sum(row["PrimaryCAZyme"] == "yes" for row in selected)
        completeness_rows.append({"Completeness": label, "CatalogGenes": len(selected), "PrimaryCAZymeGenes": primary_count, "PrimaryPercent": 100 * primary_count / len(selected)})
    write_tsv(summary / "completeness-cazyme-summary.tsv", completeness_rows)

    cgc_rows = read_tsv(work / "btheta/cgc_standard_out.tsv")
    substrate_rows = read_tsv(work / "btheta/substrate_prediction.tsv")
    cgcs = {row["CGC#"] for row in cgc_rows}
    predicted = [row for row in substrate_rows if (row.get("dbCAN-PUL substrate", "") or row.get("dbCAN-sub substrate", "")).strip()]
    cgc_summary = []
    for gene_type, count in sorted(collections.Counter(row["Gene Type"] for row in cgc_rows).items(), key=lambda x: (-x[1], x[0])):
        cgc_summary.append({"GeneType": gene_type, "Genes": count, "CGCs": len(cgcs), "CGCsWithPredictedSubstrate": len(predicted)})
    write_tsv(summary / "btheta-cgc-summary.tsv", cgc_summary)
    sub_counts: collections.Counter[str] = collections.Counter()
    for row in substrate_rows:
        substrate = (row.get("dbCAN-PUL substrate", "") or row.get("dbCAN-sub substrate", "")).strip()
        if substrate:
            for item in substrate.split(";"):
                if item.strip():
                    label = item.strip()
                    if label == "hostglycan":
                        label = "host glycan"
                    sub_counts[label] += 1
    write_tsv(summary / "btheta-substrate-summary.tsv", [{"PredictedSubstrate": label, "CGCs": count} for label, count in sub_counts.most_common()])

    resource_rows = [parse_time(path) for path in sorted((work / "logs").glob("*.time.txt"))]
    write_tsv(summary / "resource-usage.tsv", resource_rows)
    run_summary = {
        "article": 37, "catalog_genes": CATALOG_GENES, "candidate_genes": len(overview),
        "primary_cazyme_genes": len(primary), "single_tool_sensitivity_genes": sum(row["EvidenceTools"] == 1 for row in gene_rows),
        "three_tool_consensus_genes": sum(row["EvidenceTools"] == 3 for row in gene_rows),
        "primary_substrate_genes": sum(row["PredictedSubstrates"] != "-" for row in primary),
        "btheta_cgcs": len(cgcs), "btheta_cgcs_with_substrate": len(predicted),
        "sample_assigned_reads": sample_totals, "primary_call": "at least two of DIAMOND, dbCAN HMM, dbCAN-sub",
        "substrate_predictions_are_not_activity": True,
    }
    dump_json(summary / "run-summary.json", run_summary)
    (work / ".article37-summary-complete").write_text("complete\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
