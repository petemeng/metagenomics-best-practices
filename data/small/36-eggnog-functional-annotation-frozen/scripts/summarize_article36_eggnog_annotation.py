#!/usr/bin/env python3
"""Summarize Article 36 KO/COG/GO coverage and operational functional dark matter."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, TextIO


CATALOG_GENES = 93_782
STATE_ORDER = (
    "No seed ortholog",
    "Orthology only",
    "Broad/family only",
    "Specific identifier",
)
COG_NAMES = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control and chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation and ribosome biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall and membrane biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolite biosynthesis and transport",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking and secretion",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages and transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}
UNINFORMATIVE_DESCRIPTION = re.compile(
    r"^(?:unknown|hypothetical protein|uncharacteri[sz]ed protein|function unknown|predicted protein|protein of unknown function)$",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_emapper(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("#query\t"):
                header = line[1:].split("\t")
                continue
            if not line or line.startswith("#"):
                continue
            if header is None:
                raise ValueError(f"Data row before #query header in {path}")
            values = line.split("\t")
            if len(values) != len(header):
                raise ValueError(f"Column mismatch in {path}: {len(values)} != {len(header)}")
            rows.append(dict(zip(header, values)))
    if header is None:
        raise ValueError(f"Missing #query header in {path}")
    return header, rows


def present(value: str | None) -> bool:
    return value not in (None, "", "-")


def informative_description(value: str | None) -> bool:
    if not present(value):
        return False
    return UNINFORMATIVE_DESCRIPTION.fullmatch(value.strip()) is None


def split_terms(value: str | None) -> list[str]:
    if not present(value):
        return []
    return sorted({term.strip() for term in value.split(",") if term.strip() and term.strip() != "-"})


def cog_letters(value: str | None) -> list[str]:
    if not present(value):
        return []
    return sorted({letter for letter in value.upper() if letter in COG_NAMES})


def annotation_state(row: dict[str, str] | None) -> str:
    if row is None:
        return "No seed ortholog"
    specific = any(
        present(row.get(field))
        for field in ("Preferred_name", "GOs", "EC", "KEGG_ko")
    )
    if specific:
        return "Specific identifier"
    broad = bool(cog_letters(row.get("COG_category"))) or present(row.get("PFAMs")) or informative_description(row.get("Description"))
    return "Broad/family only" if broad else "Orthology only"


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str] | None = None) -> None:
    materialized = list(rows)
    if not materialized:
        raise ValueError(f"Refusing empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    names = fields or list(materialized[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(materialized)


def gzip_text_writer(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(compressed, encoding="utf-8", newline="")


def parse_time_file(path: Path, step: str) -> dict[str, object]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if ": " in stripped:
            key, value = stripped.rsplit(": ", 1)
            values[key] = value
    return {
        "Step": step,
        "UserSeconds": values.get("User time (seconds)", "NA"),
        "SystemSeconds": values.get("System time (seconds)", "NA"),
        "Elapsed": values.get("Elapsed (wall clock) time (h:mm:ss or m:ss)", "NA"),
        "MaxRSSKiB": values.get("Maximum resident set size (kbytes)", "NA"),
        "FileSystemInputs": values.get("File system inputs", "NA"),
        "FileSystemOutputs": values.get("File system outputs", "NA"),
        "ExitStatus": values.get("Exit status", "NA"),
    }


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    if not (work / ".article36-run-complete").is_file():
        raise FileNotFoundError("Article 36 annotation run is incomplete")
    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    metadata_rows = read_gzip_tsv(
        root / "data/small/34-nonredundant-gene-catalog-frozen/primary-catalog-representatives.tsv.gz"
    )
    metadata = {row["RepresentativeID"]: row for row in metadata_rows}
    if len(metadata) != CATALOG_GENES:
        raise ValueError(f"Expected {CATALOG_GENES} metadata rows, observed {len(metadata)}")

    raw_counts: dict[str, dict[str, int]] = {"MOCK1": {}, "MOCK2": {}}
    abundance_path = root / "data/small/35-gene-abundance-frozen/gene-abundance-long.tsv.gz"
    with gzip.open(abundance_path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            raw_counts[row["Sample"]][row["GeneID"]] = int(row["RawCount"])
    for sample in raw_counts:
        if set(raw_counts[sample]) != set(metadata):
            raise ValueError(f"{sample} raw-count IDs do not match catalog metadata")
    sample_totals = {sample: sum(counts.values()) for sample, counts in raw_counts.items()}

    main_header, main_rows = read_emapper(work / "annotation/main/catalog-main.emapper.annotations")
    all_header, all_rows = read_emapper(work / "annotation/go-all/catalog-go-all.emapper.annotations")
    main_by_gene = {row["query"]: row for row in main_rows}
    all_by_gene = {row["query"]: row for row in all_rows}
    if len(main_by_gene) != len(main_rows) or len(all_by_gene) != len(all_rows):
        raise ValueError("Duplicate query IDs in eggNOG-mapper output")
    if not set(main_by_gene).issubset(metadata) or not set(all_by_gene).issubset(metadata):
        raise ValueError("eggNOG-mapper returned query IDs absent from the catalog")
    if set(main_by_gene) != set(all_by_gene):
        raise ValueError("Primary and all-GO annotation branches do not share the same seed-hit genes")

    state_count: Counter[str] = Counter()
    state_reads: dict[str, Counter[str]] = {sample: Counter() for sample in raw_counts}
    field_count: Counter[str] = Counter()
    field_reads: dict[str, Counter[str]] = {sample: Counter() for sample in raw_counts}
    overlap_count: Counter[str] = Counter()
    overlap_reads: dict[str, Counter[str]] = {sample: Counter() for sample in raw_counts}
    completeness_count: Counter[tuple[str, str]] = Counter()
    length_count: Counter[tuple[str, str]] = Counter()
    length_totals: Counter[str] = Counter()
    cog_gene_count: Counter[str] = Counter()
    cog_fractional: Counter[str] = Counter()
    cog_reads_fractional: dict[str, Counter[str]] = {sample: Counter() for sample in raw_counts}
    ko_gene_count: Counter[str] = Counter()
    go_gene_count: Counter[str] = Counter()
    ko_fractional: Counter[str] = Counter()
    go_fractional: Counter[str] = Counter()
    ko_reads_fractional: dict[str, Counter[str]] = {sample: Counter() for sample in raw_counts}
    go_reads_fractional: dict[str, Counter[str]] = {sample: Counter() for sample in raw_counts}
    all_go_gene_count: Counter[str] = Counter()
    main_go_links = 0
    all_go_links = 0
    genes_gaining_go = 0
    go_subset_violations = 0

    compact_fields = [
        "GeneID", "Completeness", "AaLength", "RawCountMOCK1", "RawCountMOCK2",
        "SeedOrtholog", "Evalue", "Score", "MaxAnnotLevel", "COGCategory",
        "Description", "PreferredName", "GOs", "EC", "KEGGKO", "PFAMs",
        "AnnotationState", "HasCOG", "HasKO", "HasGO", "HasEC", "HasPFAM",
    ]
    compact_path = summary_dir / "gene-functional-annotation.tsv.gz"
    with gzip_text_writer(compact_path) as compact_handle:
        writer = csv.DictWriter(compact_handle, fieldnames=compact_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene_id in sorted(metadata):
            meta = metadata[gene_id]
            row = main_by_gene.get(gene_id)
            state = annotation_state(row)
            state_count[state] += 1
            for sample in raw_counts:
                state_reads[sample][state] += raw_counts[sample][gene_id]

            has_seed = row is not None
            has_description = has_seed and informative_description(row.get("Description"))
            has_preferred = has_seed and present(row.get("Preferred_name"))
            cogs = cog_letters(row.get("COG_category") if row else None)
            kos = split_terms(row.get("KEGG_ko") if row else None)
            gos = split_terms(row.get("GOs") if row else None)
            ecs = split_terms(row.get("EC") if row else None)
            pfams = split_terms(row.get("PFAMs") if row else None)
            field_flags = {
                "Seed ortholog": has_seed,
                "Informative description": has_description,
                "Preferred name": has_preferred,
                "COG category": bool(cogs),
                "KEGG Orthology": bool(kos),
                "GO (non-electronic)": bool(gos),
                "EC number": bool(ecs),
                "PFAM transfer": bool(pfams),
            }
            for field, flag in field_flags.items():
                if flag:
                    field_count[field] += 1
                    for sample in raw_counts:
                        field_reads[sample][field] += raw_counts[sample][gene_id]

            combination = f"COG={'yes' if cogs else 'no'} | KO={'yes' if kos else 'no'} | GO={'yes' if gos else 'no'}"
            overlap_count[combination] += 1
            for sample in raw_counts:
                overlap_reads[sample][combination] += raw_counts[sample][gene_id]

            completeness = meta["Completeness"]
            completeness_count[(completeness, state)] += 1
            aa_length = int(meta["AaLength"])
            if aa_length < 100:
                length_bin = "<100 aa"
            elif aa_length < 200:
                length_bin = "100-199 aa"
            elif aa_length < 400:
                length_bin = "200-399 aa"
            else:
                length_bin = ">=400 aa"
            length_count[(length_bin, state)] += 1
            length_totals[length_bin] += 1

            if cogs:
                weight = 1.0 / len(cogs)
                for category in cogs:
                    cog_gene_count[category] += 1
                    cog_fractional[category] += weight
                    for sample in raw_counts:
                        cog_reads_fractional[sample][category] += raw_counts[sample][gene_id] * weight
            if kos:
                weight = 1.0 / len(kos)
                for term in kos:
                    ko_gene_count[term] += 1
                    ko_fractional[term] += weight
                    for sample in raw_counts:
                        ko_reads_fractional[sample][term] += raw_counts[sample][gene_id] * weight
            if gos:
                weight = 1.0 / len(gos)
                for term in gos:
                    go_gene_count[term] += 1
                    go_fractional[term] += weight
                    for sample in raw_counts:
                        go_reads_fractional[sample][term] += raw_counts[sample][gene_id] * weight

            all_gos = split_terms(all_by_gene.get(gene_id, {}).get("GOs"))
            main_go_links += len(gos)
            all_go_links += len(all_gos)
            genes_gaining_go += int(bool(set(all_gos) - set(gos)))
            go_subset_violations += int(not set(gos).issubset(all_gos))
            for term in all_gos:
                all_go_gene_count[term] += 1

            writer.writerow(
                {
                    "GeneID": gene_id,
                    "Completeness": completeness,
                    "AaLength": aa_length,
                    "RawCountMOCK1": raw_counts["MOCK1"][gene_id],
                    "RawCountMOCK2": raw_counts["MOCK2"][gene_id],
                    "SeedOrtholog": row.get("seed_ortholog", "-") if row else "-",
                    "Evalue": row.get("evalue", "-") if row else "-",
                    "Score": row.get("score", "-") if row else "-",
                    "MaxAnnotLevel": row.get("max_annot_lvl", "-") if row else "-",
                    "COGCategory": row.get("COG_category", "-") if row else "-",
                    "Description": row.get("Description", "-") if row else "-",
                    "PreferredName": row.get("Preferred_name", "-") if row else "-",
                    "GOs": row.get("GOs", "-") if row else "-",
                    "EC": row.get("EC", "-") if row else "-",
                    "KEGGKO": row.get("KEGG_ko", "-") if row else "-",
                    "PFAMs": row.get("PFAMs", "-") if row else "-",
                    "AnnotationState": state,
                    "HasCOG": int(bool(cogs)),
                    "HasKO": int(bool(kos)),
                    "HasGO": int(bool(gos)),
                    "HasEC": int(bool(ecs)),
                    "HasPFAM": int(bool(pfams)),
                }
            )

    state_rows = []
    for state in STATE_ORDER:
        row: dict[str, object] = {
            "AnnotationState": state,
            "Genes": state_count[state],
            "GenePercent": 100 * state_count[state] / CATALOG_GENES,
        }
        for sample in ("MOCK1", "MOCK2"):
            row[f"{sample}RawReads"] = state_reads[sample][state]
            row[f"{sample}ReadPercent"] = 100 * state_reads[sample][state] / sample_totals[sample]
        state_rows.append(row)
    write_tsv(summary_dir / "annotation-state-summary.tsv", state_rows)

    field_order = (
        "Seed ortholog", "Informative description", "Preferred name", "COG category",
        "KEGG Orthology", "GO (non-electronic)", "EC number", "PFAM transfer",
    )
    coverage_rows = []
    for field in field_order:
        row = {"Field": field, "Genes": field_count[field], "GenePercent": 100 * field_count[field] / CATALOG_GENES}
        for sample in ("MOCK1", "MOCK2"):
            row[f"{sample}RawReads"] = field_reads[sample][field]
            row[f"{sample}ReadPercent"] = 100 * field_reads[sample][field] / sample_totals[sample]
        coverage_rows.append(row)
    write_tsv(summary_dir / "field-coverage-summary.tsv", coverage_rows)

    overlap_rows = []
    for combination, genes in sorted(overlap_count.items(), key=lambda item: (-item[1], item[0])):
        row = {"EvidenceCombination": combination, "Genes": genes, "GenePercent": 100 * genes / CATALOG_GENES}
        for sample in ("MOCK1", "MOCK2"):
            row[f"{sample}RawReads"] = overlap_reads[sample][combination]
            row[f"{sample}ReadPercent"] = 100 * overlap_reads[sample][combination] / sample_totals[sample]
        overlap_rows.append(row)
    write_tsv(summary_dir / "annotation-evidence-overlap.tsv", overlap_rows)

    completeness_totals = Counter(row["Completeness"] for row in metadata_rows)
    completeness_rows = []
    for completeness in ("Complete", "Partial", "Incomplete"):
        for state in STATE_ORDER:
            count = completeness_count[(completeness, state)]
            completeness_rows.append(
                {
                    "Completeness": completeness,
                    "AnnotationState": state,
                    "Genes": count,
                    "WithinCompletenessPercent": 100 * count / completeness_totals[completeness],
                }
            )
    write_tsv(summary_dir / "completeness-annotation-summary.tsv", completeness_rows)

    length_rows = []
    for length_bin in ("<100 aa", "100-199 aa", "200-399 aa", ">=400 aa"):
        for state in STATE_ORDER:
            count = length_count[(length_bin, state)]
            length_rows.append(
                {
                    "LengthBin": length_bin,
                    "AnnotationState": state,
                    "Genes": count,
                    "WithinLengthBinPercent": 100 * count / length_totals[length_bin],
                }
            )
    write_tsv(summary_dir / "length-annotation-summary.tsv", length_rows)

    cog_rows = []
    for category in sorted(COG_NAMES):
        row = {
            "COGCategory": category,
            "Description": COG_NAMES[category],
            "GenesWithCategory": cog_gene_count[category],
            "FractionalGeneEquivalent": cog_fractional[category],
        }
        for sample in ("MOCK1", "MOCK2"):
            row[f"{sample}FractionalRawReads"] = cog_reads_fractional[sample][category]
        cog_rows.append(row)
    write_tsv(summary_dir / "cog-category-summary.tsv", cog_rows)

    def term_rows(
        gene_counts: Counter[str],
        fractional: Counter[str],
        read_fractional: dict[str, Counter[str]],
    ) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for term, genes in sorted(gene_counts.items(), key=lambda item: (-item[1], item[0])):
            rows.append(
                {
                    "Term": term,
                    "GenesWithTerm": genes,
                    "FractionalGeneEquivalent": fractional[term],
                    "MOCK1FractionalRawReads": read_fractional["MOCK1"][term],
                    "MOCK2FractionalRawReads": read_fractional["MOCK2"][term],
                }
            )
        return rows

    write_tsv(summary_dir / "ko-term-summary.tsv", term_rows(ko_gene_count, ko_fractional, ko_reads_fractional))
    write_tsv(summary_dir / "go-term-summary.tsv", term_rows(go_gene_count, go_fractional, go_reads_fractional))

    fractional_rows = []
    for layer, field, fractional, read_fractional in (
        ("COG", "COG category", cog_fractional, cog_reads_fractional),
        ("KO", "KEGG Orthology", ko_fractional, ko_reads_fractional),
        ("GO", "GO (non-electronic)", go_fractional, go_reads_fractional),
    ):
        row: dict[str, object] = {
            "Layer": layer,
            "GenesWithAnyTerm": field_count[field],
            "FractionalGeneEquivalent": sum(fractional.values()),
            "GeneDifference": sum(fractional.values()) - field_count[field],
        }
        for sample in ("MOCK1", "MOCK2"):
            assigned_term_reads = field_reads[sample][field]
            fractional_term_reads = sum(read_fractional[sample].values())
            row[f"{sample}ReadsFromGenesWithAnyTerm"] = assigned_term_reads
            row[f"{sample}FractionalRawReads"] = fractional_term_reads
            row[f"{sample}Difference"] = fractional_term_reads - assigned_term_reads
        fractional_rows.append(row)
    write_tsv(summary_dir / "fractional-allocation-audit.tsv", fractional_rows)

    main_go_genes = sum(bool(split_terms(row.get("GOs"))) for row in main_rows)
    all_go_genes = sum(bool(split_terms(row.get("GOs"))) for row in all_rows)
    go_audit_rows = [
        {"Metric": "Genes with >=1 GO term", "NonElectronic": main_go_genes, "AllEvidence": all_go_genes, "Difference": all_go_genes - main_go_genes},
        {"Metric": "Unique GO terms", "NonElectronic": len(go_gene_count), "AllEvidence": len(all_go_gene_count), "Difference": len(all_go_gene_count) - len(go_gene_count)},
        {"Metric": "Gene-GO links", "NonElectronic": main_go_links, "AllEvidence": all_go_links, "Difference": all_go_links - main_go_links},
        {"Metric": "Genes gaining >=1 GO term", "NonElectronic": 0, "AllEvidence": genes_gaining_go, "Difference": genes_gaining_go},
    ]
    write_tsv(summary_dir / "go-evidence-audit.tsv", go_audit_rows)

    paper_rows = [
        {"Completeness": "Complete", "EggNOGAnnotatedGenes": 5_354_169},
        {"Completeness": "Partial", "EggNOGAnnotatedGenes": 8_374_034},
        {"Completeness": "Incomplete", "EggNOGAnnotatedGenes": 7_865_395},
        {"Completeness": "Total", "EggNOGAnnotatedGenes": 21_593_598},
    ]
    write_tsv(summary_dir / "delgado-table3-eggnog-anchor.tsv", paper_rows)

    resource_rows = [
        parse_time_file(work / "logs/emapper-main.time.txt", "Main non-electronic GO"),
        parse_time_file(work / "logs/emapper-go-all.time.txt", "All-evidence GO sensitivity"),
    ]
    write_tsv(summary_dir / "resource-usage.tsv", resource_rows)

    md5_field = "md5" if "md5" in main_header else None
    md5_nonempty = sum(present(row.get(md5_field)) for row in main_rows) if md5_field else 0
    run_summary = {
        "article": 36,
        "catalog_genes": CATALOG_GENES,
        "seed_ortholog_genes": len(main_rows),
        "annotation_states": dict(state_count),
        "field_gene_counts": dict(field_count),
        "sample_raw_count_totals": sample_totals,
        "main_go_genes": main_go_genes,
        "all_evidence_go_genes": all_go_genes,
        "genes_gaining_go": genes_gaining_go,
        "go_subset_violations": go_subset_violations,
        "main_go_links": main_go_links,
        "all_evidence_go_links": all_go_links,
        "unique_kos": len(ko_gene_count),
        "unique_go_terms_non_electronic": len(go_gene_count),
        "cog_categories_observed": sum(cog_gene_count[letter] > 0 for letter in COG_NAMES),
        "main_output_rows": len(main_rows),
        "go_all_output_rows": len(all_rows),
        "main_output_columns": main_header,
        "go_all_output_columns": all_header,
        "main_md5_nonempty": md5_nonempty,
        "state_partition_closes": sum(state_count.values()) == CATALOG_GENES,
        "catalog_percent_denominator": CATALOG_GENES,
        "paper_bags_total_genes": 67_566_251,
        "paper_bags_eggnog_annotated_genes": 21_593_598,
        "paper_bags_eggnog_percent": 100 * 21_593_598 / 67_566_251,
        "functional_annotations_are_predictions": True,
        "absence_is_not_gene_absence": True,
        "go_all_is_primary": False,
        "multi_label_cog_uses_fractional_allocation": True,
        "multi_label_terms_use_fractional_allocation": True,
    }
    (summary_dir / "run-summary.json").write_text(json.dumps(run_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (work / ".article36-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(run_summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
