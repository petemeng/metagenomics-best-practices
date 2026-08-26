#!/usr/bin/env python3
"""Summarize the checksum-locked Article 59 DRAM/METABOLIC real-data run."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_GENOMES = 28
EXPECTED_PRIMARY = 24
MODULE_CUTOFF = 0.75
KO_PATTERN = re.compile(r"K\d{5}")
DRAM_HEATMAP_MODULES = (
    "M00001", "M00004", "M00008", "M00009", "M00012", "M00165",
    "M00173", "M00374", "M00375", "M00376", "M00377", "M00422",
    "M00567",
)
KEY_FUNCTIONS = (
    "Methane oxidation - Partculate methane monooxygenase",
    "Methane oxidation - Soluble methane monoxygenase",
    "Methane production",
    "CBB cycle - Rubisco (Form I)",
    "CBB cycle - Rubisco (Form II)",
    "Wood Ljungdahl pathway",
    "Reverse TCA cycle",
    "Ammonia oxidation",
    "N2 fixation",
    "Nitrite oxidation",
    "Nitrate reduction",
    "Nitrite reduction to ammonia",
    "Nitrite reduction",
    "Nitric oxide reduction",
    "Nitrous oxide reduction",
    "Anammox",
    "Sulfide oxidation",
    "Sulfite reduction",
    "Sulfur oxidation",
    "Sulfur reduction",
    "Thiosulfate oxidation",
    "Sulfate reduction",
    "Thiosulfate disproportionation",
    "Arsenite oxidation",
    "Arsenate reduction",
    "Metal (Iron/Manganese) reduction",
)


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
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: str | None) -> float:
    if value is None or value.strip() in {"", "NA", "N/A", "None", "nan"}:
        return math.nan
    return float(value)


def as_bool_text(value: bool) -> str:
    return "true" if value else "false"


def extract_kos(value: str | None) -> set[str]:
    return set(KO_PATTERN.findall(value or ""))


def parse_dram_kos(path: Path) -> tuple[dict[str, set[str]], dict[str, int]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "fasta" not in reader.fieldnames:
            raise RuntimeError("DRAM annotations.tsv lacks the fasta column")
        identifier_columns = [
            name
            for name in reader.fieldnames
            if name.lower() in {"ko_id", "kegg_id", "kofam_id"}
            or (
                ("kegg" in name.lower() or "kofam" in name.lower())
                and name.lower().endswith("_id")
            )
        ]
        if not identifier_columns:
            identifier_columns = [
                name for name in reader.fieldnames if "kegg" in name.lower() or "kofam" in name.lower()
            ]
        if not identifier_columns:
            raise RuntimeError(f"No KOfam/KEGG identifier column in {path}")
        ko_sets: dict[str, set[str]] = defaultdict(set)
        annotation_rows: Counter[str] = Counter()
        for row in reader:
            genome = row["fasta"]
            if not genome:
                raise RuntimeError("DRAM annotation row has an empty fasta identifier")
            annotation_rows[genome] += 1
            for column in identifier_columns:
                ko_sets[genome].update(extract_kos(row.get(column)))
    return dict(ko_sets), dict(annotation_rows)


def parse_metabolic_kos(result_dir: Path) -> tuple[dict[str, set[str]], dict[str, dict[str, int]]]:
    ko_sets: dict[str, set[str]] = {}
    ko_counts: dict[str, dict[str, int]] = {}
    files = sorted(result_dir.glob("*.result.txt"))
    if len(files) != EXPECTED_GENOMES:
        raise RuntimeError(f"Expected {EXPECTED_GENOMES} METABOLIC KO tables, observed {len(files)}")
    for path in files:
        genome = path.name.removesuffix(".result.txt")
        counts: dict[str, int] = {}
        with path.open(encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                fields = raw.rstrip("\n").split("\t")
                ko = fields[0].strip() if fields else ""
                if len(fields) < 2 or not KO_PATTERN.fullmatch(ko):
                    raise RuntimeError(f"Malformed METABOLIC KO row {path}:{line_number}")
                count = int(fields[1]) if fields[1].strip() else 0
                if ko in counts:
                    raise RuntimeError(f"Duplicate METABOLIC KO row after normalization: {path}:{ko}")
                counts[ko] = count
        ko_counts[genome] = counts
        ko_sets[genome] = {ko for ko, count in counts.items() if count > 0}
    return ko_sets, ko_counts


def read_dram_product(path: Path) -> tuple[dict[str, dict[str, float]], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError("DRAM product table has no header")
        genome_column = reader.fieldnames[0]
        features = reader.fieldnames[1:]
        result: dict[str, dict[str, float]] = {}
        for row in reader:
            genome = row[genome_column]
            values: dict[str, float] = {}
            for feature in features:
                raw = row.get(feature, "")
                if raw in {"", "NA", "NaN", "nan", None}:
                    values[feature] = math.nan
                elif str(raw).lower() in {"true", "false"}:
                    values[feature] = 1.0 if str(raw).lower() == "true" else 0.0
                else:
                    values[feature] = float(raw)
            result[genome] = values
    return result, features


def module_name_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for row in read_tsv(path):
        module = row["module"]
        name = row["module_name"]
        if module in mapping and mapping[module] != name:
            raise RuntimeError(f"Conflicting DRAM module name for {module}")
        mapping[module] = name
    return mapping


def parse_metabolic_modules(path: Path) -> tuple[dict[str, dict[str, bool]], dict[str, str]]:
    rows = read_tsv(path)
    if not rows:
        raise RuntimeError("METABOLIC worksheet3 is empty")
    genome_columns = [column for column in rows[0] if column.endswith(" Module presence")]
    if len(genome_columns) != EXPECTED_GENOMES:
        raise RuntimeError(f"Expected {EXPECTED_GENOMES} worksheet3 genome columns")
    modules: dict[str, dict[str, bool]] = {}
    names: dict[str, str] = {}
    for row in rows:
        module = row["Module ID"]
        names[module] = row["Module"]
        modules[module] = {
            column.removesuffix(" Module presence"): row[column] == "Present"
            for column in genome_columns
        }
    return modules, names


def parse_metabolic_functions(path: Path) -> list[dict[str, object]]:
    rows = read_tsv(path)
    if not rows:
        raise RuntimeError("METABOLIC worksheet2 is empty")
    genome_columns = [column for column in rows[0] if column.endswith(" Function presence")]
    if len(genome_columns) != EXPECTED_GENOMES:
        raise RuntimeError(f"Expected {EXPECTED_GENOMES} worksheet2 genome columns")
    output: list[dict[str, object]] = []
    observed_functions = {row["Function"] for row in rows}
    missing = sorted(set(KEY_FUNCTIONS) - observed_functions)
    if missing:
        raise RuntimeError(f"Expected curated METABOLIC functions are missing: {missing}")
    for row in rows:
        if row["Function"] not in KEY_FUNCTIONS:
            continue
        for column in genome_columns:
            output.append(
                {
                    "Genome": column.removesuffix(" Function presence"),
                    "Category": row["Category"],
                    "Function": row["Function"],
                    "GeneRule": row["Gene abbreviation"],
                    "Present": as_bool_text(row[column] == "Present"),
                    "EvidenceType": "METABOLIC curated HMM/motif rule",
                }
            )
    return output


def elapsed_seconds(value: str) -> float:
    fields = value.strip().split(":")
    if len(fields) == 2:
        return 60 * float(fields[0]) + float(fields[1])
    if len(fields) == 3:
        return 3600 * float(fields[0]) + 60 * float(fields[1]) + float(fields[2])
    raise ValueError(f"Unsupported elapsed time: {value}")


def parse_time_file(path: Path) -> dict[str, float]:
    fields: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ": " in line:
            key, value = line.strip().rsplit(": ", 1)
            fields[key] = value
    elapsed_key = next((key for key in fields if key.startswith("Elapsed (wall clock)")), None)
    if elapsed_key is None:
        raise RuntimeError(f"No GNU-time elapsed field in {path}")
    return {
        "WallSeconds": elapsed_seconds(fields[elapsed_key]),
        "CPUSeconds": float(fields["User time (seconds)"]) + float(fields["System time (seconds)"]),
        "PeakRSSGiB": float(fields["Maximum resident set size (kbytes)"]) / 1024 / 1024,
        "FileSystemInputs": float(fields.get("File system inputs", "0")),
        "FileSystemOutputs": float(fields.get("File system outputs", "0")),
    }


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    if not (work / ".article59-runs-complete").is_file():
        raise FileNotFoundError("Run scripts/run_article59_metabolism.py first")
    summary_dir = work / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    ledger = read_tsv(work / "input-mag-ledger.tsv")
    if len(ledger) != EXPECTED_GENOMES:
        raise RuntimeError(f"Expected {EXPECTED_GENOMES} ledger rows")
    metadata = {row["Genome"]: row for row in ledger}
    genomes = sorted(metadata)
    primary = [genome for genome in genomes if metadata[genome]["AnalysisSet"] == "Primary real MAG"]
    if len(primary) != EXPECTED_PRIMARY:
        raise RuntimeError(f"Expected {EXPECTED_PRIMARY} primary MAGs")

    dram_full_kos, dram_annotation_rows = parse_dram_kos(work / "dram-annotation/annotations.tsv")
    metabolic_kos, metabolic_ko_counts = parse_metabolic_kos(
        work / "metabolic-output/KEGG_identifier_result"
    )
    if set(dram_full_kos) != set(genomes) or set(metabolic_kos) != set(genomes):
        raise RuntimeError("DRAM/METABOLIC genome identifiers do not match the input ledger")
    comparison_universes = {frozenset(counts) for counts in metabolic_ko_counts.values()}
    if len(comparison_universes) != 1:
        raise RuntimeError("METABOLIC per-genome KO grids are not identical")
    comparison_universe = set(next(iter(comparison_universes)))
    if len(comparison_universe) != 2_678:
        raise RuntimeError(
            f"Expected a 2,678-KO comparison universe, observed {len(comparison_universe)}"
        )
    dram_comparable_kos = {
        genome: dram_full_kos[genome] & comparison_universe for genome in genomes
    }

    ko_counts_rows: list[dict[str, object]] = []
    ko_concordance_rows: list[dict[str, object]] = []
    for genome in genomes:
        dram_set = dram_comparable_kos[genome]
        metabolic_set = metabolic_kos[genome]
        shared = dram_set & metabolic_set
        union = dram_set | metabolic_set
        row_meta = metadata[genome]
        ko_counts_rows.append(
            {
                "Genome": genome,
                "AnalysisSet": row_meta["AnalysisSet"],
                "Species": row_meta["Species"],
                "Domain": row_meta["Domain"],
                "Phylum": row_meta["Phylum"],
                "CompletenessPct": row_meta["CompletenessPct"],
                "ContaminationPct": row_meta["ContaminationPct"],
                "GenomeBp": row_meta["GenomeBp"],
                "DRAMAnnotationRows": dram_annotation_rows[genome],
                "DRAMFullKOCount": len(dram_full_kos[genome]),
                "DRAMComparableKOCount": len(dram_set),
                "METABOLICKOCount": len(metabolic_set),
                "SharedKOCount": len(shared),
                "DRAMOnlyKOCount": len(dram_set - metabolic_set),
                "METABOLICOnlyKOCount": len(metabolic_set - dram_set),
                "UnionKOCount": len(union),
                "Jaccard": len(shared) / len(union) if union else math.nan,
                "ComparisonKOUniverse": len(comparison_universe),
            }
        )
        for ko in sorted(union):
            dram_present = ko in dram_set
            metabolic_present = ko in metabolic_set
            pattern = "Both" if dram_present and metabolic_present else "DRAM only" if dram_present else "METABOLIC only"
            ko_concordance_rows.append(
                {
                    "Genome": genome,
                    "KO": ko,
                    "DRAMPresent": as_bool_text(dram_present),
                    "METABOLICPresent": as_bool_text(metabolic_present),
                    "METABOLICHitCount": metabolic_ko_counts[genome].get(ko, 0),
                    "PresencePattern": pattern,
                }
            )
    write_tsv(summary_dir / "ko-counts-by-genome.tsv", ko_counts_rows)
    write_tsv(summary_dir / "ko-tool-concordance.tsv", ko_concordance_rows)

    dram_product, product_features = read_dram_product(work / "dram-distill/product.tsv")
    dram_module_names = module_name_map(work / "database/dram-kofam-2026-06-01/module_step_form.tsv")
    metabolic_modules, metabolic_module_names = parse_metabolic_modules(
        work / "metabolic-output/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet3.tsv"
    )
    missing_modules = [module for module in DRAM_HEATMAP_MODULES if module not in metabolic_modules]
    if missing_modules:
        raise RuntimeError(f"METABOLIC lacks shared modules: {missing_modules}")
    module_rows: list[dict[str, object]] = []
    for module in DRAM_HEATMAP_MODULES:
        dram_name = dram_module_names[module]
        if dram_name not in product_features:
            raise RuntimeError(f"DRAM product lacks expected module column: {module} {dram_name}")
        for genome in genomes:
            coverage = dram_product[genome][dram_name]
            if math.isnan(coverage):
                raise RuntimeError(f"Missing DRAM module coverage for {genome}/{module}")
            dram_present = coverage >= MODULE_CUTOFF
            metabolic_present = metabolic_modules[module][genome]
            row_meta = metadata[genome]
            module_rows.append(
                {
                    "Genome": genome,
                    "AnalysisSet": row_meta["AnalysisSet"],
                    "Species": row_meta["Species"],
                    "Phylum": row_meta["Phylum"],
                    "ModuleID": module,
                    "DRAMModuleName": dram_name,
                    "METABOLICModuleName": metabolic_module_names[module],
                    "DRAMCoverage": coverage,
                    "DRAMPresentAt075": as_bool_text(dram_present),
                    "METABOLICPresentAt075": as_bool_text(metabolic_present),
                    "Agreement": as_bool_text(dram_present == metabolic_present),
                    "ComparisonBoundary": "same KOfam snapshot; tool-specific module-step definitions",
                }
            )
    write_tsv(summary_dir / "module-tool-concordance.tsv", module_rows)
    write_tsv(
        summary_dir / "pathway-evidence-matrix.tsv",
        [row for row in module_rows if row["AnalysisSet"] == "Primary real MAG"],
    )
    module_summary_rows: list[dict[str, object]] = []
    for module in DRAM_HEATMAP_MODULES:
        selected = [row for row in module_rows if row["ModuleID"] == module and row["AnalysisSet"] == "Primary real MAG"]
        agreements = sum(row["Agreement"] == "true" for row in selected)
        module_summary_rows.append(
            {
                "ModuleID": module,
                "DRAMModuleName": selected[0]["DRAMModuleName"],
                "PrimaryGenomes": len(selected),
                "Agreements": agreements,
                "AgreementRate": agreements / len(selected),
                "DRAMPresent": sum(row["DRAMPresentAt075"] == "true" for row in selected),
                "METABOLICPresent": sum(row["METABOLICPresentAt075"] == "true" for row in selected),
            }
        )
    write_tsv(summary_dir / "module-agreement-summary.tsv", module_summary_rows)

    function_rows = parse_metabolic_functions(
        work / "metabolic-output/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet2.tsv"
    )
    for row in function_rows:
        row_meta = metadata[str(row["Genome"])]
        row.update(
            {
                "AnalysisSet": row_meta["AnalysisSet"],
                "Species": row_meta["Species"],
                "Phylum": row_meta["Phylum"],
            }
        )
    function_rows.sort(key=lambda row: (str(row["Category"]), str(row["Function"]), str(row["Genome"])))
    write_tsv(summary_dir / "key-process-evidence.tsv", function_rows)

    counts_by_genome = {str(row["Genome"]): row for row in ko_counts_rows}
    absence_rows: list[dict[str, object]] = []
    for genome in primary:
        row_meta = metadata[genome]
        completeness = as_float(row_meta["CompletenessPct"])
        contamination = as_float(row_meta["ContaminationPct"])
        if contamination >= 5:
            risk = "Composite-pathway risk"
            rule = "Review contig provenance before assigning a complete pathway to one genome"
        elif completeness < 70:
            risk = "High missingness risk"
            rule = "Do not interpret an unobserved gene or module as biological absence"
        elif completeness < 90:
            risk = "Moderate missingness risk"
            rule = "Treat absence as weak evidence and inspect assembly/gene-call neighborhoods"
        else:
            risk = "Lower missingness risk"
            rule = "Absence remains provisional until key-gene and neighborhood review"
        module_subset = [row for row in module_rows if row["Genome"] == genome]
        absence_rows.append(
            {
                "Genome": genome,
                "Species": row_meta["Species"],
                "Phylum": row_meta["Phylum"],
                "CompletenessPct": completeness,
                "ContaminationPct": contamination,
                "Contigs": row_meta["Contigs"],
                "DRAMFullKOCount": counts_by_genome[genome]["DRAMFullKOCount"],
                "DRAMComparableKOCount": counts_by_genome[genome]["DRAMComparableKOCount"],
                "METABOLICKOCount": counts_by_genome[genome]["METABOLICKOCount"],
                "KOJaccard": counts_by_genome[genome]["Jaccard"],
                "DRAMModulesAt075": sum(row["DRAMPresentAt075"] == "true" for row in module_subset),
                "METABOLICModulesAt075": sum(row["METABOLICPresentAt075"] == "true" for row in module_subset),
                "AbsenceRisk": risk,
                "InterpretationRule": rule,
            }
        )
    write_tsv(summary_dir / "completeness-absence-audit.tsv", absence_rows)

    parent = "SGB_002"
    truncation_rows: list[dict[str, object]] = []
    for genome in ("TRUNC_050", "TRUNC_070", "TRUNC_090", "TRUNC_100"):
        row_meta = metadata[genome]
        dram_parent = dram_full_kos[parent]
        dram_comparable_parent = dram_comparable_kos[parent]
        metabolic_parent = metabolic_kos[parent]
        dram_current = dram_full_kos[genome]
        dram_comparable_current = dram_comparable_kos[genome]
        metabolic_current = metabolic_kos[genome]
        parent_module_rows = {row["ModuleID"]: row for row in module_rows if row["Genome"] == parent}
        current_module_rows = {row["ModuleID"]: row for row in module_rows if row["Genome"] == genome}
        dram_parent_present = {module for module, row in parent_module_rows.items() if row["DRAMPresentAt075"] == "true"}
        dram_current_present = {module for module, row in current_module_rows.items() if row["DRAMPresentAt075"] == "true"}
        metabolic_parent_present = {module for module, row in parent_module_rows.items() if row["METABOLICPresentAt075"] == "true"}
        metabolic_current_present = {module for module, row in current_module_rows.items() if row["METABOLICPresentAt075"] == "true"}
        truncation_rows.append(
            {
                "Genome": genome,
                "ParentGenome": parent,
                "RetentionTargetPct": row_meta["RetentionTargetPct"],
                "RetentionObservedPct": row_meta["RetentionObservedPct"],
                "GenomeBp": row_meta["GenomeBp"],
                "Contigs": row_meta["Contigs"],
                "DRAMFullKOCount": len(dram_current),
                "DRAMFullKORetentionPct": 100 * len(dram_current & dram_parent) / len(dram_parent),
                "DRAMFullNewKOCount": len(dram_current - dram_parent),
                "DRAMComparableKOCount": len(dram_comparable_current),
                "DRAMComparableKORetentionPct": 100 * len(dram_comparable_current & dram_comparable_parent) / len(dram_comparable_parent),
                "DRAMComparableNewKOCount": len(dram_comparable_current - dram_comparable_parent),
                "METABOLICKOCount": len(metabolic_current),
                "METABOLICKORetentionPct": 100 * len(metabolic_current & metabolic_parent) / len(metabolic_parent),
                "METABOLICNewKOCount": len(metabolic_current - metabolic_parent),
                "DRAMParentModules": len(dram_parent_present),
                "DRAMRetainedParentModules": len(dram_current_present & dram_parent_present),
                "METABOLICParentModules": len(metabolic_parent_present),
                "METABOLICRetainedParentModules": len(metabolic_current_present & metabolic_parent_present),
                "ExactParentSequence": as_bool_text(
                    row_meta["SequenceSetSHA256"] == metadata[parent]["SequenceSetSHA256"]
                ),
                "ExactDRAMKOReproduction": as_bool_text(
                    dram_current == dram_parent
                    and dram_comparable_current == dram_comparable_parent
                ),
                "ExactMETABOLICKOReproduction": as_bool_text(metabolic_current == metabolic_parent),
                "DeterministicSeed": row_meta["DeterministicSeed"],
            }
        )
    write_tsv(summary_dir / "truncation-sensitivity.tsv", truncation_rows)

    evidence_ladder = [
        {"Rank": 1, "Evidence": "Adaptive-threshold profile-HMM hit", "ClaimCeiling": "A homologous gene family was detected", "RequiredUpgrade": "Inspect alignment, domains and competing annotations"},
        {"Rank": 2, "Evidence": "Curated multi-gene or motif rule", "ClaimCeiling": "A diagnostic metabolic trait is supported", "RequiredUpgrade": "Check copy context, gene order and contamination"},
        {"Rank": 3, "Evidence": "Module steps above a declared cutoff", "ClaimCeiling": "A pathway is genomically encoded in part or in full", "RequiredUpgrade": "Report missing steps and alternative routes"},
        {"Rank": 4, "Evidence": "Quality-audited MAG plus neighborhood evidence", "ClaimCeiling": "The genome has metabolic potential", "RequiredUpgrade": "Confirm taxon abundance and pathway coherence across samples"},
        {"Rank": 5, "Evidence": "RNA, protein, metabolite, isotope or flux evidence", "ClaimCeiling": "Activity or flux is supported under measured conditions", "RequiredUpgrade": "Use condition-matched experiments and causal validation"},
    ]
    write_tsv(summary_dir / "evidence-ladder.tsv", evidence_ladder)

    time_labels = (
        "dram-database",
        "metabolic-dbcan-hmmpress",
        "metabolic-merops-diamond",
        "metabolic-prodigal",
        "dram-annotate",
        "dram-merge",
        "dram-distill",
        "metabolic-g",
    )
    resources: list[dict[str, object]] = []
    for label in time_labels:
        measurements = parse_time_file(work / "logs" / f"{label}.time.txt")
        resources.append({"Label": label, **measurements})
    write_tsv(summary_dir / "resource-usage.tsv", resources)

    versions = [
        {"Tool": "DRAM", "Version": "1.5.0", "ReleaseEvidence": "Python package metadata and exact conda lock"},
        {"Tool": "DRAM source forms", "Version": "commit fe61d759303f30db058d5d505c448b28e41b03f1", "ReleaseEvidence": "five checksum-locked GPL-3.0 forms"},
        {"Tool": "METABOLIC-G", "Version": "v4.0 commit 97236332519180f1d76a242dedb0aaa8191fdbb3", "ReleaseEvidence": "three checksum-locked GPL-3.0 archives"},
        {"Tool": "KOfamKOALA", "Version": "archive 2026-06-01", "ReleaseEvidence": "profiles.tar.gz and ko_list.gz SHA-256 manifest"},
        {"Tool": "HMMER", "Version": "3.4", "ReleaseEvidence": "exact conda locks"},
        {"Tool": "Prodigal", "Version": "2.6.3", "ReleaseEvidence": "meta mode in both workflows"},
        {"Tool": "MMseqs2", "Version": "18.8cc5c", "ReleaseEvidence": "DRAM exact conda lock"},
        {"Tool": "DIAMOND", "Version": "2.2.2", "ReleaseEvidence": "METABOLIC exact conda lock"},
        {"Tool": "dbCAN HMM database", "Version": "5.2.9 release 2026-05-05", "ReleaseEvidence": "SHA-256 396c791f...baa7a5"},
        {"Tool": "MEROPS pepunit.lib", "Version": "current file dated 2023-02-22", "ReleaseEvidence": "SHA-256 5deb3c49...9dcd74"},
        {"Tool": "GTDB taxonomy context", "Version": "R232", "ReleaseEvidence": "Article 46 frozen genome ledger"},
        {"Tool": "Plotting stack", "Version": "pandas 2.2.3; matplotlib 3.11.1; seaborn 0.13.2; Pillow 12.3.0", "ReleaseEvidence": "env/drep-linux-64.lock reused as the exact figure environment"},
    ]
    write_tsv(summary_dir / "tool-versions.tsv", versions)

    primary_ko_rows = [row for row in ko_counts_rows if row["AnalysisSet"] == "Primary real MAG"]
    primary_module_rows = [row for row in module_rows if row["AnalysisSet"] == "Primary real MAG"]
    module_agreements = sum(row["Agreement"] == "true" for row in primary_module_rows)
    key_primary = [row for row in function_rows if row["AnalysisSet"] == "Primary real MAG"]
    summary = {
        "article": 59,
        "primary_real_mags": len(primary),
        "deterministic_truncation_genomes": EXPECTED_GENOMES - len(primary),
        "kofam_release": "2026-06-01",
        "metabolic_compatible_profiles": 2643,
        "metabolic_requested_profiles": 2644,
        "metabolic_missing_profile": "K18513",
        "metabolic_per_genome_ko_grid_rows": 2678,
        "metabolic_custom_hmm_to_ko_rules_beyond_kofam_subset": 35,
        "dram_full_kofam_profiles": 27754,
        "ko_comparison_universe": len(comparison_universe),
        "module_cutoff": MODULE_CUTOFF,
        "shared_modules": len(DRAM_HEATMAP_MODULES),
        "primary_module_comparisons": len(primary_module_rows),
        "primary_module_agreements": module_agreements,
        "primary_module_agreement_rate": module_agreements / len(primary_module_rows),
        "primary_median_ko_jaccard": statistics.median(float(row["Jaccard"]) for row in primary_ko_rows),
        "primary_median_dram_full_ko_count": statistics.median(int(row["DRAMFullKOCount"]) for row in primary_ko_rows),
        "primary_median_dram_comparable_ko_count": statistics.median(int(row["DRAMComparableKOCount"]) for row in primary_ko_rows),
        "primary_median_metabolic_ko_count": statistics.median(int(row["METABOLICKOCount"]) for row in primary_ko_rows),
        "key_curated_functions": len(KEY_FUNCTIONS),
        "primary_key_function_positive_cells": sum(row["Present"] == "true" for row in key_primary),
        "truncation_parent": parent,
        "truncation_seed": 59002,
        "truncation_100_exact_dram_ko_reproduction": truncation_rows[-1]["ExactDRAMKOReproduction"] == "true",
        "truncation_100_exact_metabolic_ko_reproduction": truncation_rows[-1]["ExactMETABOLICKOReproduction"] == "true",
        "gene_hit_is_activity": False,
        "tool_output_is_flux": False,
        "random_process_after_input_preparation": False,
    }
    (summary_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (work / ".article59-summary-complete").write_text(
        "Article 59 summaries completed with fail-closed checks.\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
