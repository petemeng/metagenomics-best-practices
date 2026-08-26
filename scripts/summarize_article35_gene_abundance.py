#!/usr/bin/env python3
"""Summarize Article 35 counts, units, functional aggregation, and resource ledgers."""

from __future__ import annotations

import argparse
import bz2
import csv
import gzip
import hashlib
import io
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO


SEED = 20260735
CATALOG_GENES = 93_782
POLICIES = ("AllPrimary", "IdentityQcov", "Main", "Strict")
MIN_UNIREF_PIDENT = 50.0
MIN_UNIREF_QUERY_COVERAGE = 80.0
MAX_UNIREF_EVALUE = 1e-5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--reaction-map", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_gzip_text(path: Path) -> TextIO:
    raw = path.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0)
    return io.TextIOWrapper(zipped, encoding="utf-8", newline="")


def write_tsv(path: Path, rows: Iterable[dict], fieldnames: list[str] | None = None) -> int:
    iterator = iter(rows)
    try:
        first = next(iterator)
    except StopIteration as exc:
        raise ValueError(f"Refusing to write empty table: {path}") from exc
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = deterministic_gzip_text(path) if path.suffix == ".gz" else path.open("w", encoding="utf-8", newline="")
    count = 0
    with handle:
        columns = fieldnames or list(first)
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(first)
        count = 1
        for row in iterator:
            writer.writerow(row)
            count += 1
    return count


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_htseq(path: Path) -> tuple[dict[str, int], dict[str, int]]:
    genes: dict[str, int] = {}
    special: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            key, value = raw.rstrip("\n").split("\t")
            if key.startswith("__"):
                special[key] = int(value)
            else:
                genes[key] = int(value)
    if len(genes) != CATALOG_GENES:
        raise ValueError(f"HTSeq output {path} has {len(genes)} genes")
    return genes, special


def legacy_summary(work: Path) -> list[dict]:
    rows: list[dict] = []
    for sample in ("MOCK1", "MOCK2"):
        genes, special = parse_htseq(work / "legacy" / f"{sample}.htseq-counts.tsv")
        assigned = sum(genes.values())
        total = assigned + sum(special.values())
        if total != 10_000:
            raise ValueError(f"Legacy HTSeq ledger does not close for {sample}: {total}")
        rows.append(
            {
                "Sample": sample,
                "SelectedForwardReads": 10_000,
                "AssignedReads": assigned,
                "AssignedPct": assigned / 100,
                "DetectedGenes": sum(value > 0 for value in genes.values()),
                "NotAligned": special.get("__not_aligned", 0),
                "NoFeature": special.get("__no_feature", 0),
                "AmbiguousFeature": special.get("__ambiguous", 0),
                "TooLowMAPQ": special.get("__too_low_aQual", 0),
                "AlignmentNotUnique": special.get("__alignment_not_unique", 0),
                "BowtieReportsOneAlignment": "yes",
                "MAPQThreshold": 0,
            }
        )
    return rows


def load_catalog(work: Path) -> tuple[list[str], dict[str, dict]]:
    rows = read_tsv(work / "reference/catalog-audit.tsv")
    if len(rows) != CATALOG_GENES:
        raise ValueError(f"Catalog metadata has {len(rows)} rows")
    catalog = {row["GeneID"]: row for row in rows}
    return sorted(catalog), catalog


def load_policy_counts(work: Path, genes: list[str]) -> tuple[dict[str, dict[str, dict[str, int]]], list[dict], list[dict]]:
    all_counts: dict[str, dict[str, dict[str, int]]] = {}
    policy_rows: list[dict] = []
    quality_rows: list[dict] = []
    for sample in ("MOCK1", "MOCK2"):
        sparse = {row["GeneID"]: {policy: int(row[policy]) for policy in POLICIES} for row in read_tsv(work / "mapping" / f"{sample}.policy-counts.tsv")}
        all_counts[sample] = {gene: sparse.get(gene, {policy: 0 for policy in POLICIES}) for gene in genes}
        summary = json.loads((work / "mapping" / f"{sample}.mapping-summary.json").read_text())
        for policy in POLICIES:
            assigned = int(summary["policy_assigned_reads"][policy])
            primary = int(summary["primary_mapped_reads"])
            total = int(summary["total_input_reads"])
            policy_rows.append(
                {
                    "Sample": sample,
                    "Policy": policy,
                    "InputReads": total,
                    "PrimaryMappedReads": primary,
                    "AssignedReads": assigned,
                    "FilteredMappedReads": primary - assigned,
                    "UnmappedReads": int(summary["unmapped_reads"]),
                    "AssignedPctInput": 100 * assigned / total,
                    "DetectedGenes": int(summary["policy_detected_genes"][policy]),
                    "SecondaryAlignments": int(summary["secondary_alignments"]),
                }
            )
        quality_rows.extend(read_tsv(work / "mapping" / f"{sample}.quality-histogram.tsv"))
    return all_counts, policy_rows, quality_rows


def best_uniref_hits(path: Path) -> tuple[dict[str, dict], dict]:
    hits: dict[str, dict] = {}
    candidate_counts: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Unexpected DIAMOND field count: {len(fields)}")
            query, subject = fields[:2]
            candidate_counts[query] += 1
            if candidate_counts[query] > 5:
                raise ValueError(f"DIAMOND returned more than five candidates for {query}")
            match = re.search(r"UniRef90_[A-Za-z0-9]+", subject)
            if not match:
                raise ValueError(f"No UniRef90 identifier in DIAMOND subject: {subject}")
            row = {
                "GeneID": query,
                "UniRef90": match.group(0),
                "Pident": float(fields[2]),
                "AlignmentLength": int(fields[3]),
                "QueryLength": int(fields[4]),
                "SubjectLength": int(fields[5]),
                "Evalue": float(fields[6]),
                "Bitscore": float(fields[7]),
                "QueryCoveragePct": float(fields[8]),
            }
            numeric = (row["Pident"], row["Evalue"], row["Bitscore"], row["QueryCoveragePct"])
            if not all(math.isfinite(value) for value in numeric):
                raise ValueError(f"Non-finite DIAMOND metric for {query}: {numeric}")
            if (
                row["Pident"] < MIN_UNIREF_PIDENT
                or row["QueryCoveragePct"] < MIN_UNIREF_QUERY_COVERAGE
                or row["Evalue"] > MAX_UNIREF_EVALUE
            ):
                continue
            previous = hits.get(query)
            key = (-row["Bitscore"], row["Evalue"], -row["Pident"], row["UniRef90"])
            if previous is None:
                hits[query] = row
            else:
                previous_key = (-previous["Bitscore"], previous["Evalue"], -previous["Pident"], previous["UniRef90"])
                if key < previous_key:
                    hits[query] = row
    audit = {
        "CandidateRows": sum(candidate_counts.values()),
        "QueriesWithCandidate": len(candidate_counts),
        "MaxCandidatesPerQuery": max(candidate_counts.values(), default=0),
        "BestHitGenes": len(hits),
        "MinimumIdentityPct": MIN_UNIREF_PIDENT,
        "MinimumQueryCoveragePct": MIN_UNIREF_QUERY_COVERAGE,
        "MaximumEvalue": MAX_UNIREF_EVALUE,
    }
    return hits, audit


def reaction_crosswalk(path: Path, wanted: set[str]) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    reactions: dict[str, set[str]] = defaultdict(set)
    ecs: dict[str, set[str]] = defaultdict(set)
    with bz2.open(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            reaction, ec = fields[:2]
            for identifier in fields[2:]:
                if identifier in wanted:
                    reactions[identifier].add(reaction)
                    if ec:
                        ecs[reaction].add(ec)
    return reactions, ecs


def ranks(values: dict[str, float]) -> dict[str, int]:
    return {key: rank for rank, (key, _) in enumerate(sorted(values.items(), key=lambda item: (-item[1], item[0])), start=1)}


def length_bin(length: int) -> str:
    if length < 300:
        return "<300 bp"
    if length < 600:
        return "300-599 bp"
    if length < 900:
        return "600-899 bp"
    if length < 1_500:
        return "900-1,499 bp"
    return ">=1,500 bp"


def parse_time(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ": " in raw:
            key, value = raw.strip().split(": ", 1)
            values[key] = value
    return values


def resource_rows(work: Path) -> list[dict]:
    commands = {row["Step"]: row for row in read_tsv(work / "command-log.tsv")}
    rows: list[dict] = []
    for path in sorted((work / "logs").glob("*.time.txt")):
        step = path.name[: -len(".time.txt")]
        values = parse_time(path)
        rows.append(
            {
                "Step": step,
                "Elapsed": values.get("Elapsed (wall clock) time (h:mm:ss or m:ss)", "NA"),
                "UserSeconds": values.get("User time (seconds)", "NA"),
                "SystemSeconds": values.get("System time (seconds)", "NA"),
                "PeakRSSKiB": values.get("Maximum resident set size (kbytes)", "NA"),
                "ExitStatus": values.get("Exit status", "NA"),
                "Command": commands.get(step, {}).get("Command", ""),
            }
        )
    if not rows:
        raise ValueError("No GNU time resource logs found")
    return rows


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article35-run-complete").is_file():
        raise SystemExit("Run run_article35_gene_abundance.py first")
    genes, catalog = load_catalog(work)
    counts, policy_rows, quality_rows = load_policy_counts(work, genes)
    legacy_rows = legacy_summary(work)
    write_tsv(work / "summary/legacy-mapping-audit.tsv", legacy_rows)
    write_tsv(work / "summary/mapping-policy-summary.tsv", policy_rows)
    write_tsv(work / "summary/alignment-quality-histogram.tsv", quality_rows)

    hits, annotation_search_audit = best_uniref_hits(work / "annotation/uniref90-top5.tsv")
    unexpected_hits = sorted(set(hits) - set(genes))
    if unexpected_hits:
        raise ValueError(f"DIAMOND returned query IDs outside the catalog: {unexpected_hits[:5]}")
    if not hits:
        raise ValueError("No catalog proteins passed the locked UniRef90 thresholds")
    write_tsv(work / "summary/annotation-search-audit.tsv", [annotation_search_audit])
    reaction_by_uniref, ec_by_reaction = reaction_crosswalk(args.reaction_map.resolve(), {row["UniRef90"] for row in hits.values()})
    annotation_rows: list[dict] = []
    for gene in genes:
        hit = hits.get(gene)
        uniref = hit["UniRef90"] if hit else ""
        reactions = sorted(reaction_by_uniref.get(uniref, set()))
        annotation_rows.append(
            {
                "GeneID": gene,
                "UniRef90": uniref,
                "Pident": "" if hit is None else hit["Pident"],
                "QueryCoveragePct": "" if hit is None else hit["QueryCoveragePct"],
                "Evalue": "" if hit is None else hit["Evalue"],
                "Bitscore": "" if hit is None else hit["Bitscore"],
                "ReactionCount": len(reactions),
                "Reactions": ";".join(reactions),
            }
        )
    annotation = {row["GeneID"]: row for row in annotation_rows}
    write_tsv(work / "summary/gene-functional-annotation.tsv.gz", annotation_rows)

    abundance_rows: list[dict] = []
    normalization_rows: list[dict] = []
    rank_rows: list[dict] = []
    length_accumulator: dict[tuple[str, str], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    completeness_accumulator: dict[tuple[str, str], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    sample_metrics: dict[str, dict] = {}
    for sample in ("MOCK1", "MOCK2"):
        raw = {gene: counts[sample][gene]["Main"] for gene in genes}
        assigned = sum(raw.values())
        expected = next(int(row["AssignedReads"]) for row in policy_rows if row["Sample"] == sample and row["Policy"] == "Main")
        if assigned != expected:
            raise ValueError(f"Main count denominator mismatch for {sample}: {assigned} != {expected}")
        rpk = {gene: raw[gene] / (int(catalog[gene]["NtLength"]) / 1_000) for gene in genes}
        rpk_sum = sum(rpk.values())
        cpm = {gene: raw[gene] * 1_000_000 / assigned for gene in genes}
        rpkm = {gene: raw[gene] * 1_000_000_000 / (int(catalog[gene]["NtLength"]) * assigned) for gene in genes}
        tpm = {gene: rpk[gene] * 1_000_000 / rpk_sum for gene in genes}
        count_rank = ranks(raw)
        tpm_rank = ranks(tpm)
        selected_rank_genes = set(sorted(genes, key=lambda gene: (-raw[gene], gene))[:100]) | set(
            sorted(genes, key=lambda gene: (-tpm[gene], gene))[:100]
        )
        for gene in sorted(selected_rank_genes):
            rank_rows.append(
                {
                    "Sample": sample,
                    "GeneID": gene,
                    "NtLength": catalog[gene]["NtLength"],
                    "RawCount": raw[gene],
                    "CPM": cpm[gene],
                    "RPKM": rpkm[gene],
                    "TPM": tpm[gene],
                    "CountRank": count_rank[gene],
                    "TPMRank": tpm_rank[gene],
                    "RankShiftTPMMinusCount": tpm_rank[gene] - count_rank[gene],
                }
            )
        for gene in genes:
            meta = catalog[gene]
            ann = annotation[gene]
            abundance_rows.append(
                {
                    "Sample": sample,
                    "GeneID": gene,
                    "NtLength": meta["NtLength"],
                    "Completeness": meta["Completeness"],
                    "RepresentativeOrigin": meta["RepresentativeOrigin"],
                    "RawCount": raw[gene],
                    "CPM": cpm[gene],
                    "RPKM": rpkm[gene],
                    "TPM": tpm[gene],
                    "UniRef90": ann["UniRef90"],
                    "ReactionCount": ann["ReactionCount"],
                }
            )
            key = (sample, length_bin(int(meta["NtLength"])))
            length_accumulator[key]["CatalogGenes"] += 1
            length_accumulator[key]["DetectedGenes"] += raw[gene] > 0
            length_accumulator[key]["AssignedReads"] += raw[gene]
            length_accumulator[key]["CPM"] += cpm[gene]
            length_accumulator[key]["TPM"] += tpm[gene]
            completeness_key = (sample, meta["Completeness"])
            completeness_accumulator[completeness_key]["CatalogGenes"] += 1
            completeness_accumulator[completeness_key]["DetectedGenes"] += raw[gene] > 0
            completeness_accumulator[completeness_key]["AssignedReads"] += raw[gene]
            completeness_accumulator[completeness_key]["CPM"] += cpm[gene]
            completeness_accumulator[completeness_key]["TPM"] += tpm[gene]
        normalization_rows.append(
            {
                "Sample": sample,
                "AssignedReads": assigned,
                "DetectedGenes": sum(value > 0 for value in raw.values()),
                "RawCountSum": sum(raw.values()),
                "CPMSum": sum(cpm.values()),
                "RPKMSum": sum(rpkm.values()),
                "TPMSum": sum(tpm.values()),
                "RPKSum": rpk_sum,
                "CPMClosureError": abs(sum(cpm.values()) - 1_000_000),
                "TPMClosureError": abs(sum(tpm.values()) - 1_000_000),
            }
        )
        sample_metrics[sample] = {"raw": raw, "cpm": cpm, "rpkm": rpkm, "tpm": tpm, "assigned": assigned}

    write_tsv(work / "summary/gene-abundance-long.tsv.gz", abundance_rows)
    write_tsv(work / "summary/unit-rank-audit.tsv", rank_rows)
    write_tsv(work / "summary/normalization-audit.tsv", normalization_rows)
    length_rows = [
        {"Sample": sample, "LengthBin": bin_name, **{key: value for key, value in values.items()}}
        for (sample, bin_name), values in sorted(length_accumulator.items())
    ]
    write_tsv(work / "summary/gene-length-bin-summary.tsv", length_rows)
    completeness_rows = [
        {"Sample": sample, "Completeness": completeness, **{key: value for key, value in values.items()}}
        for (sample, completeness), values in sorted(completeness_accumulator.items())
    ]
    write_tsv(work / "summary/gene-completeness-summary.tsv", completeness_rows)

    reaction_rows: list[dict] = []
    functional_audit: list[dict] = []
    for sample in ("MOCK1", "MOCK2"):
        metrics = sample_metrics[sample]
        aggregate: dict[str, dict[str, object]] = defaultdict(
            lambda: {
                "genes": set(),
                "split_reads": 0.0,
                "split_cpm": 0.0,
                "split_tpm": 0.0,
                "copy_reads": 0.0,
                "copy_cpm": 0.0,
                "copy_tpm": 0.0,
            }
        )
        genes_with_hit = 0
        genes_with_reaction = 0
        assigned_reads_with_hit = 0
        assigned_reads_with_reaction = 0
        multi_reaction_genes = 0
        for gene in genes:
            ann = annotation[gene]
            uniref = ann["UniRef90"]
            reactions = ann["Reactions"].split(";") if ann["Reactions"] else []
            if uniref:
                genes_with_hit += 1
                assigned_reads_with_hit += metrics["raw"][gene]
            if reactions:
                genes_with_reaction += 1
                assigned_reads_with_reaction += metrics["raw"][gene]
                multi_reaction_genes += len(reactions) > 1
                destinations = reactions
            elif uniref:
                destinations = ["UNIREF90_NO_REACTION"]
            else:
                destinations = ["NO_UNIREF90_HIT"]
            divisor = len(destinations)
            for reaction in destinations:
                item = aggregate[reaction]
                item["genes"].add(gene)
                item["split_reads"] += metrics["raw"][gene] / divisor
                item["split_cpm"] += metrics["cpm"][gene] / divisor
                item["split_tpm"] += metrics["tpm"][gene] / divisor
                item["copy_reads"] += metrics["raw"][gene]
                item["copy_cpm"] += metrics["cpm"][gene]
                item["copy_tpm"] += metrics["tpm"][gene]
        for reaction in sorted(aggregate):
            item = aggregate[reaction]
            reaction_rows.append(
                {
                    "Sample": sample,
                    "Reaction": reaction,
                    "EC": ";".join(sorted(ec_by_reaction.get(reaction, set()))),
                    "GenesContributing": len(item["genes"]),
                    "SplitReadEquivalent": item["split_reads"],
                    "CopyReadEquivalent": item["copy_reads"],
                    "SplitCPM": item["split_cpm"],
                    "SplitTPM": item["split_tpm"],
                    "CopyCPM": item["copy_cpm"],
                    "CopyTPM": item["copy_tpm"],
                }
            )
        split_cpm_total = sum(item["split_cpm"] for item in aggregate.values())
        split_tpm_total = sum(item["split_tpm"] for item in aggregate.values())
        split_read_total = sum(item["split_reads"] for item in aggregate.values())
        copy_cpm_total = sum(item["copy_cpm"] for item in aggregate.values())
        copy_tpm_total = sum(item["copy_tpm"] for item in aggregate.values())
        copy_read_total = sum(item["copy_reads"] for item in aggregate.values())
        functional_audit.append(
            {
                "Sample": sample,
                "CatalogGenes": len(genes),
                "GenesWithUniRef90Hit": genes_with_hit,
                "GenesWithReaction": genes_with_reaction,
                "MultiReactionGenes": multi_reaction_genes,
                "AssignedReads": metrics["assigned"],
                "AssignedReadsWithUniRef90Hit": assigned_reads_with_hit,
                "AssignedReadsWithReaction": assigned_reads_with_reaction,
                "AssignedReadPctWithReaction": 100 * assigned_reads_with_reaction / metrics["assigned"],
                "SplitReadEquivalentSum": split_read_total,
                "CopyReadEquivalentSum": copy_read_total,
                "CopyReadInflation": copy_read_total / split_read_total,
                "GeneCPMSum": sum(metrics["cpm"].values()),
                "SplitCPMSum": split_cpm_total,
                "CopyCPMSum": copy_cpm_total,
                "CopyCPMInflation": copy_cpm_total / split_cpm_total,
                "GeneTPMSum": sum(metrics["tpm"].values()),
                "SplitTPMSum": split_tpm_total,
                "CopyTPMSum": copy_tpm_total,
                "CopyTPMInflation": copy_tpm_total / split_tpm_total,
            }
        )
    write_tsv(work / "summary/reaction-abundance-long.tsv", reaction_rows)
    write_tsv(work / "summary/functional-aggregation-audit.tsv", functional_audit)

    resources = resource_rows(work)
    write_tsv(work / "summary/resource-usage.tsv", resources)
    primary_rows = {row["Sample"]: row for row in policy_rows if row["Policy"] == "Main"}
    summary = {
        "article": 35,
        "status": "completed",
        "seed": SEED,
        "catalog_genes": len(genes),
        "samples": 2,
        "input_reads": {sample: int(primary_rows[sample]["InputReads"]) for sample in primary_rows},
        "raw_mapping_rate_pct": {sample: 100 * int(primary_rows[sample]["PrimaryMappedReads"]) / int(primary_rows[sample]["InputReads"]) for sample in primary_rows},
        "main_assigned_reads": {sample: int(primary_rows[sample]["AssignedReads"]) for sample in primary_rows},
        "main_assigned_pct": {sample: float(primary_rows[sample]["AssignedPctInput"]) for sample in primary_rows},
        "main_detected_genes": {sample: int(primary_rows[sample]["DetectedGenes"]) for sample in primary_rows},
        "legacy_assigned_pct": {row["Sample"]: row["AssignedPct"] for row in legacy_rows},
        "uniref90_hit_genes": len(hits),
        "diamond_candidate_rows": annotation_search_audit["CandidateRows"],
        "diamond_max_candidates_per_query": annotation_search_audit["MaxCandidatesPerQuery"],
        "genes_with_reaction": {row["Sample"]: row["GenesWithReaction"] for row in functional_audit},
        "assigned_read_pct_with_reaction": {row["Sample"]: row["AssignedReadPctWithReaction"] for row in functional_audit},
        "copy_cpm_inflation": {row["Sample"]: row["CopyCPMInflation"] for row in functional_audit},
        "resource_steps": len(resources),
        "boundaries": {
            "mapping_rate_is_assigned_fraction": False,
            "rpkm_is_raw_count_matrix": False,
            "tpm_is_absolute_abundance": False,
            "copy_aggregation_conserves_mass": False,
            "two_mocks_are_biological_replicates": False,
            "catalog_is_independent_of_profiled_reads": False,
            "partial_orf_hit_proves_complete_function": False,
        },
        "checksums": {
            name: sha256(work / "summary" / name)
            for name in (
                "mapping-policy-summary.tsv",
                "normalization-audit.tsv",
                "functional-aggregation-audit.tsv",
                "reaction-abundance-long.tsv",
                "gene-abundance-long.tsv.gz",
            )
        },
    }
    (work / "summary/run-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (work / ".article35-summary-complete").write_text("completed\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
