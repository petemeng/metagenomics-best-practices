#!/usr/bin/env python3
"""Fail-closed validation for Article 59 DRAM/METABOLIC evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from PIL import Image

from article42_44_validation_utils import (
    Audit,
    as_bool,
    audit_chapter,
    audit_checksums,
    audit_figures,
    finish,
    read_tsv,
    sha256,
)


FIGURES = (
    "59-pathway-module-heatmap",
    "59-key-process-evidence",
    "59-tool-concordance",
    "59-completeness-absence-audit",
    "59-truncation-sensitivity",
    "59-metabolism-evidence-ladder",
)
MODULES = {
    "M00001", "M00004", "M00008", "M00009", "M00012", "M00165",
    "M00173", "M00374", "M00375", "M00376", "M00377", "M00422",
    "M00567",
}
TIME_LABELS = {
    "dram-database",
    "metabolic-dbcan-hmmpress",
    "metabolic-merops-diamond",
    "metabolic-prodigal",
    "dram-annotate",
    "dram-merge",
    "dram-distill",
    "metabolic-g",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stage-dir", type=Path, required=True)
    parser.add_argument("--plot-python", type=Path, default=Path(sys.executable))
    return parser.parse_args()


def near(value: object, expected: float, tolerance: float = 1e-8) -> bool:
    try:
        return math.isclose(float(value), expected, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def image_pixel_sha(path: Path) -> str:
    with Image.open(path) as image:
        normalized = image.convert("RGBA")
        payload = normalized.size[0].to_bytes(8, "little") + normalized.size[1].to_bytes(8, "little") + normalized.tobytes()
    return hashlib.sha256(payload).hexdigest()


def stage_reanalysis(
    root: Path, frozen: Path, stage: Path, audit: Audit, plot_python: Path
) -> None:
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    shutil.copy2(frozen / "input-mag-ledger.tsv", stage / "input-mag-ledger.tsv")
    (stage / ".article59-runs-complete").write_text("staged frozen evidence\n", encoding="utf-8")

    annotation_target = stage / "dram-annotation/annotations.tsv"
    annotation_target.parent.mkdir(parents=True)
    with gzip.open(frozen / "dram/annotations.tsv.gz", "rb") as source, annotation_target.open("wb") as target:
        shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
    distill = stage / "dram-distill"
    distill.mkdir(parents=True)
    shutil.copy2(frozen / "dram/product.tsv", distill / "product.tsv")
    module_target = stage / "database/dram-kofam-2026-06-01/module_step_form.tsv"
    module_target.parent.mkdir(parents=True)
    shutil.copy2(frozen / "dram/module_step_form.tsv", module_target)

    worksheet_target = stage / "metabolic-output/METABOLIC_result_each_spreadsheet"
    worksheet_target.mkdir(parents=True)
    for path in sorted((frozen / "metabolic/worksheets").glob("*.tsv")):
        shutil.copy2(path, worksheet_target / path.name)
    ko_target = stage / "metabolic-output/KEGG_identifier_result"
    ko_target.mkdir(parents=True)
    for path in sorted((frozen / "metabolic/KEGG_identifier_result").glob("*.txt")):
        shutil.copy2(path, ko_target / path.name)
    log_target = stage / "logs"
    log_target.mkdir(parents=True)
    for label in TIME_LABELS:
        shutil.copy2(frozen / "logs" / f"{label}.time.txt", log_target / f"{label}.time.txt")

    environment = os.environ.copy()
    environment["MPLCONFIGDIR"] = str(stage / "matplotlib")
    summary_command = [
        "python3",
        str(root / "scripts/summarize_article59_metabolism.py"),
        "--project-root",
        str(root),
        "--work-dir",
        str(stage),
    ]
    completed = subprocess.run(summary_command, cwd=root, env=environment, capture_output=True, text=True, check=False)
    audit.add("Reanalysis", "summary-exit", completed.returncode == 0, completed.stderr[-2000:])
    if completed.returncode != 0:
        return
    frozen_summary_files = sorted(
        path.name for path in frozen.glob("*") if path.is_file() and (path.suffix in {".tsv", ".json"})
        and path.name not in {
            "frozen-contract.json", "preparation-contract.json", "run-contract.json", "dram-config.json"
        }
    )
    staged_summary_files = sorted(path.name for path in (stage / "summary").glob("*") if path.is_file())
    expected_summary_files = sorted(
        name for name in frozen_summary_files if (frozen / name).is_file() and (stage / "summary" / name).is_file()
    )
    audit.add("Reanalysis", "summary-file-set", set(staged_summary_files) == set(expected_summary_files), {"stage": staged_summary_files, "expected": expected_summary_files})
    for name in expected_summary_files:
        audit.add(
            "Reanalysis",
            f"summary-hash-{name}",
            sha256(stage / "summary" / name) == sha256(frozen / name),
            sha256(stage / "summary" / name),
        )

    staged_figures = stage / "figures"
    plot_command = [
        str(plot_python),
        str(root / "scripts/plot_article59_metabolism.py"),
        "--summary-dir",
        str(stage / "summary"),
        "--figure-dir",
        str(staged_figures),
    ]
    plotted = subprocess.run(plot_command, cwd=root, env=environment, capture_output=True, text=True, check=False)
    audit.add("Reanalysis", "plot-exit", plotted.returncode == 0, plotted.stderr[-2000:])
    if plotted.returncode != 0:
        return
    audit_figures(staged_figures, audit, FIGURES)
    for stem in FIGURES:
        staged_png = staged_figures / f"{stem}.png"
        published_png = root / "figures" / f"{stem}.png"
        if staged_png.is_file() and published_png.is_file():
            audit.add(
                "Reanalysis",
                f"pixel-match-{stem}",
                image_pixel_sha(staged_png) == image_pixel_sha(published_png),
                image_pixel_sha(staged_png),
            )


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    chapter = args.chapter.resolve()
    figures = args.figure_dir.resolve()
    output = args.output_dir.resolve()
    stage = args.stage_dir.resolve()
    audit = Audit()

    audit_checksums(frozen, audit)
    contract = json.loads((frozen / "frozen-contract.json").read_text(encoding="utf-8"))
    audit.add("Contract", "article", contract.get("article") == 59, contract)
    audit.add("Contract", "real-primary-count", contract.get("primary_real_mags") == 24, contract)
    audit.add("Contract", "sensitivity-count", contract.get("deterministic_sensitivity_genomes") == 4, contract)
    audit.add("Contract", "seed", contract.get("seed") == 59002, contract)
    audit.add("Contract", "large-db-excluded", contract.get("large_databases_included") is False, contract)
    audit.add("Contract", "mag-fastas-excluded", contract.get("full_mag_fastas_included") is False, contract)

    manifest = read_tsv(frozen / "database-manifest.tsv")
    audit.add("Database", "manifest-assets", len(manifest) == 12, len(manifest))
    audit.add("Database", "manifest-tools", {row["Tool"] for row in manifest} == {"KOfamKOALA", "METABOLIC", "run_dbcan", "MEROPS", "DRAM"}, sorted({row["Tool"] for row in manifest}))
    expected_assets = {
        "profiles.tar.gz": ("1557643213", "038cb13c41bda8f97ba0d57299aebe74042fcc90fbf450cba191448e8c19606f"),
        "ko_list.gz": ("909392", "d418b074049c6687e6f6d119af2981dc4e36c378d43027ecf5cc6778fd27041a"),
        "dbCAN.hmm": ("129842960", "396c791fb13defab152864d8046687ef03b492937ec6f9bbab008bc77cbaa7a5"),
        "pepunit.lib": ("457531072", "5deb3c49b2c8ec17b040d0afc26dcdcbc78481c400cdbdf77d2bd0c7ff9dcd74"),
    }
    manifest_by_asset = {row["Asset"]: row for row in manifest}
    for asset, (size, digest) in expected_assets.items():
        row = manifest_by_asset.get(asset, {})
        audit.add("Database", f"locked-{asset}", row.get("ExpectedBytes") == size and row.get("SHA256") == digest, row)
    database_audit = read_tsv(frozen / "database-audit.tsv")
    audit.add("Database", "database-audit-rows", len(database_audit) == 12, len(database_audit))
    for row in database_audit:
        audit.add("Database", f"checksum-{row['Asset']}", as_bool(row["ChecksumPass"]) and row["ExpectedSHA256"] == row["ObservedSHA256"] and row["ExpectedBytes"] == row["ObservedBytes"], row)

    compatibility = read_tsv(frozen / "kofam-compatibility-audit.tsv")
    present = [row for row in compatibility if as_bool(row["PresentInKOfam2026_06_01"])]
    missing = [row for row in compatibility if not as_bool(row["PresentInKOfam2026_06_01"])]
    audit.add("Database", "metabolic-requested-profiles", len(compatibility) == 2644, len(compatibility))
    audit.add("Database", "metabolic-compatible-profiles", len(present) == 2643, len(present))
    audit.add(
        "Database",
        "missing-profile-fail-closed",
        len(missing) == 1
        and missing[0]["Profile"] == "K18513.hmm"
        and missing[0]["CompatibilityAction"] == "excluded_missing_in_locked_snapshot",
        missing,
    )

    ledger = read_tsv(frozen / "input-mag-ledger.tsv")
    audit.add("Input", "ledger-count", len(ledger) == 28, len(ledger))
    audit.add("Input", "genome-ids-unique", len({row["Genome"] for row in ledger}) == 28, [row["Genome"] for row in ledger])
    analysis_counts = Counter(row["AnalysisSet"] for row in ledger)
    audit.add("Input", "analysis-sets", analysis_counts == Counter({"Primary real MAG": 24, "Deterministic truncation sensitivity": 4}), analysis_counts)
    primary = [row for row in ledger if row["AnalysisSet"] == "Primary real MAG"]
    truncations = [row for row in ledger if row["AnalysisSet"] == "Deterministic truncation sensitivity"]
    audit.add("Input", "gtdb-release", all(row["GTDBRelease"] == "R232" for row in ledger), sorted({row["GTDBRelease"] for row in ledger}))
    audit.add("Input", "primary-quality-complete", all(row["CompletenessPct"] not in {"", "NA"} and row["ContaminationPct"] not in {"", "NA"} for row in primary), primary)
    audit.add("Input", "primary-sha256", all(len(row["RepresentativeSHA256"]) == 64 for row in primary), primary)
    audit.add(
        "Input",
        "sequence-set-sha256",
        all(len(row.get("SequenceSetSHA256", "")) == 64 for row in ledger),
        [row.get("SequenceSetSHA256", "") for row in ledger],
    )
    audit.add("Input", "truncation-seed", all(row["DeterministicSeed"] == "59002" and row["ParentGenome"] == "SGB_002" for row in truncations), truncations)
    audit.add("Input", "truncation-targets", {row["RetentionTargetPct"] for row in truncations} == {"50", "70", "90", "100"}, truncations)
    parent_hash = next(row["SequenceSetSHA256"] for row in primary if row["Genome"] == "SGB_002")
    exact_hash = next(row["SequenceSetSHA256"] for row in truncations if row["Genome"] == "TRUNC_100")
    audit.add("Input", "truncation-100-sequence-set", exact_hash == parent_hash, {"parent": parent_hash, "control": exact_hash})

    with gzip.open(frozen / "dram/annotations.tsv.gz", "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        annotation_rows = sum(1 for _ in handle)
    audit.add("Raw", "dram-annotation-header", "fasta" in header and any("kegg" in column.lower() or "kofam" in column.lower() for column in header), header)
    audit.add("Raw", "dram-annotation-size", annotation_rows > 30_000, annotation_rows)
    product = read_tsv(frozen / "dram/product.tsv")
    audit.add("Raw", "dram-product-genomes", len(product) == 28, len(product))
    genome_stats = read_tsv(frozen / "dram/genome_stats.tsv")
    audit.add("Raw", "dram-genome-stats", len(genome_stats) == 28, len(genome_stats))
    worksheets = frozen / "metabolic/worksheets"
    for index in range(1, 7):
        rows = read_tsv(worksheets / f"METABOLIC_result_worksheet{index}.tsv")
        audit.add("Raw", f"worksheet-{index}-nonempty", len(rows) > 0, len(rows))
    worksheet2 = read_tsv(worksheets / "METABOLIC_result_worksheet2.tsv")
    worksheet3 = read_tsv(worksheets / "METABOLIC_result_worksheet3.tsv")
    audit.add("Raw", "worksheet2-genomes", len([column for column in worksheet2[0] if column.endswith(" Function presence")]) == 28, list(worksheet2[0]))
    audit.add("Raw", "worksheet3-genomes", len([column for column in worksheet3[0] if column.endswith(" Module presence")]) == 28, list(worksheet3[0]))
    ko_results = sorted((frozen / "metabolic/KEGG_identifier_result").glob("*.result.txt"))
    ko_hits = sorted((frozen / "metabolic/KEGG_identifier_result").glob("*.hits.txt"))
    audit.add("Raw", "metabolic-result-pairs", len(ko_results) == 28 and len(ko_hits) == 28, {"results": len(ko_results), "hits": len(ko_hits)})
    for path in ko_results:
        lines = path.read_text(encoding="utf-8").splitlines()
        audit.add("Raw", f"ko-grid-{path.stem}", len(lines) == 2678 and all(line.startswith("K") and "\t" in line for line in lines), len(lines))
    metabolic_log = (frozen / "metabolic/METABOLIC_log.log").read_text(encoding="utf-8", errors="replace")
    for message in ("The hmmsearch is finished", "dbCAN2 searching is done", "MEROPS peptidase searching is done", "METABOLIC-G was done"):
        audit.add("Raw", f"metabolic-log-{message}", message in metabolic_log, message)
    audit.add("Raw", "metabolic-no-bioperl-failure", "Can't locate Bio/SeqIO.pm" not in metabolic_log, "BioPerl error absent")
    protein_id_audit = read_tsv(frozen / "metabolic-protein-id-audit.tsv")
    audit.add(
        "Raw",
        "globally-unique-protein-ids",
        len(protein_id_audit) == 1
        and int(protein_id_audit[0]["ProteinFiles"]) == 28
        and int(protein_id_audit[0]["ProteinRecords"]) > 30_000
        and int(protein_id_audit[0]["OutputIDsGenomePrefixed"])
        == int(protein_id_audit[0]["UniqueOutputProteinIDs"])
        and int(protein_id_audit[0]["DuplicateOutputProteinIDs"]) == 0,
        protein_id_audit,
    )

    counts = read_tsv(frozen / "ko-counts-by-genome.tsv")
    audit.add("KO", "count-table", len(counts) == 28, len(counts))
    count_by_genome = {row["Genome"]: row for row in counts}
    for row in counts:
        union = int(row["UnionKOCount"])
        shared = int(row["SharedKOCount"])
        arithmetic = shared + int(row["DRAMOnlyKOCount"]) + int(row["METABOLICOnlyKOCount"])
        audit.add("KO", f"arithmetic-{row['Genome']}", union == arithmetic and near(row["Jaccard"], shared / union) and int(row["DRAMFullKOCount"]) >= int(row["DRAMComparableKOCount"]) > 0 and int(row["METABOLICKOCount"]) > 0 and int(row["ComparisonKOUniverse"]) == 2678, row)
    concordance = read_tsv(frozen / "ko-tool-concordance.tsv")
    concordance_by_genome: dict[str, Counter[str]] = defaultdict(Counter)
    for row in concordance:
        concordance_by_genome[row["Genome"]][row["PresencePattern"]] += 1
        audit.add("KO-row", f"{row['Genome']}-{row['KO']}", row["PresencePattern"] in {"Both", "DRAM only", "METABOLIC only"} and as_bool(row["DRAMPresent"]) == (row["PresencePattern"] in {"Both", "DRAM only"}) and as_bool(row["METABOLICPresent"]) == (row["PresencePattern"] in {"Both", "METABOLIC only"}), row)
    for genome, row in count_by_genome.items():
        counter = concordance_by_genome[genome]
        audit.add("KO", f"long-table-{genome}", counter["Both"] == int(row["SharedKOCount"]) and counter["DRAM only"] == int(row["DRAMOnlyKOCount"]) and counter["METABOLIC only"] == int(row["METABOLICOnlyKOCount"]), counter)

    modules = read_tsv(frozen / "module-tool-concordance.tsv")
    audit.add("Module", "row-grid", len(modules) == 28 * 13, len(modules))
    audit.add("Module", "module-set", {row["ModuleID"] for row in modules} == MODULES, sorted({row["ModuleID"] for row in modules}))
    for row in modules:
        coverage = float(row["DRAMCoverage"])
        expected_dram = coverage >= 0.75
        expected_agreement = expected_dram == as_bool(row["METABOLICPresentAt075"])
        audit.add("Module-row", f"{row['Genome']}-{row['ModuleID']}", 0 <= coverage <= 1 and as_bool(row["DRAMPresentAt075"]) == expected_dram and as_bool(row["Agreement"]) == expected_agreement, row)
    pathway = read_tsv(frozen / "pathway-evidence-matrix.tsv")
    audit.add("Module", "primary-grid", len(pathway) == 24 * 13 and {row["AnalysisSet"] for row in pathway} == {"Primary real MAG"}, len(pathway))

    key_processes = read_tsv(frozen / "key-process-evidence.tsv")
    functions = {row["Function"] for row in key_processes}
    audit.add("Trait", "key-process-grid", len(key_processes) == 28 * 26 and len(functions) == 26, {"rows": len(key_processes), "functions": len(functions)})
    audit.add("Trait", "critical-functions", {"Ammonia oxidation", "Methane production", "Sulfate reduction", "N2 fixation", "Anammox"} <= functions, sorted(functions))
    for genome in {row["Genome"] for row in key_processes}:
        rows = [row for row in key_processes if row["Genome"] == genome]
        audit.add("Trait", f"grid-{genome}", len(rows) == 26 and all(row["Present"] in {"true", "false"} for row in rows), rows)

    absence = read_tsv(frozen / "completeness-absence-audit.tsv")
    audit.add("Absence", "primary-count", len(absence) == 24, len(absence))
    audit.add("Absence", "risk-language", all("absence" in row["InterpretationRule"].lower() or "pathway" in row["InterpretationRule"].lower() for row in absence), absence)
    truncation = read_tsv(frozen / "truncation-sensitivity.tsv")
    audit.add("Sensitivity", "truncation-count", len(truncation) == 4, len(truncation))
    audit.add("Sensitivity", "retention-ranges", all(0 <= float(row["DRAMFullKORetentionPct"]) <= 100 and 0 <= float(row["DRAMComparableKORetentionPct"]) <= 100 and 0 <= float(row["METABOLICKORetentionPct"]) <= 100 for row in truncation), truncation)
    exact = next(row for row in truncation if row["Genome"] == "TRUNC_100")
    audit.add("Sensitivity", "exact-parent-reproduction", as_bool(exact["ExactParentSequence"]) and as_bool(exact["ExactDRAMKOReproduction"]) and as_bool(exact["ExactMETABOLICKOReproduction"]), exact)
    audit.add("Sensitivity", "no-new-kos", all(int(row["DRAMFullNewKOCount"]) == 0 and int(row["DRAMComparableNewKOCount"]) == 0 and int(row["METABOLICNewKOCount"]) == 0 for row in truncation), truncation)

    ladder = read_tsv(frozen / "evidence-ladder.tsv")
    audit.add("Evidence", "ladder", len(ladder) == 5 and [int(row["Rank"]) for row in ladder] == [1, 2, 3, 4, 5], ladder)
    audit.add("Evidence", "claim-ceilings", ladder[0]["ClaimCeiling"] == "A homologous gene family was detected" and "Activity or flux" in ladder[-1]["ClaimCeiling"], ladder)
    resources = read_tsv(frozen / "resource-usage.tsv")
    audit.add("Resource", "labels", {row["Label"] for row in resources} == TIME_LABELS, sorted(row["Label"] for row in resources))
    for row in resources:
        audit.add("Resource", f"positive-{row['Label']}", float(row["WallSeconds"]) > 0 and float(row["CPUSeconds"]) >= 0 and float(row["PeakRSSGiB"]) > 0, row)
    versions = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    exact_versions = {
        "DRAM": "1.5.0",
        "METABOLIC-G": "v4.0 commit 97236332519180f1d76a242dedb0aaa8191fdbb3",
        "KOfamKOALA": "archive 2026-06-01",
        "HMMER": "3.4",
        "Prodigal": "2.6.3",
        "MMseqs2": "18.8cc5c",
        "DIAMOND": "2.2.2",
        "GTDB taxonomy context": "R232",
    }
    audit.add("Environment", "exact-versions", all(versions.get(tool) == version for tool, version in exact_versions.items()), versions)
    dram_lock = (frozen / "env/dram-linux-64.lock").read_text(encoding="utf-8")
    metabolic_lock = (frozen / "env/metabolic-linux-64.lock").read_text(encoding="utf-8")
    plotting_lock = (frozen / "env/drep-linux-64.lock").read_text(encoding="utf-8")
    audit.add("Environment", "dram-lock", all(token in dram_lock for token in ("dram-1.5.0", "mmseqs2-18.8cc5c", "numpy-1.23.5", "setuptools-80.10.2")), "DRAM exact lock")
    audit.add("Environment", "metabolic-lock", all(token in metabolic_lock for token in ("diamond-2.2.2", "hmmer-3.4", "prodigal-2.6.3", "r-base-4.4.3")), "METABOLIC exact lock")
    audit.add("Environment", "plotting-lock", all(token in plotting_lock for token in ("pandas-2.2.3", "matplotlib-base-3.11.1", "seaborn-0.13.2", "pillow-12.3.0")), "Article 59 exact plotting lock")

    summary = json.loads((frozen / "summary.json").read_text(encoding="utf-8"))
    audit.add("Summary", "core-contract", summary.get("article") == 59 and summary.get("primary_real_mags") == 24 and summary.get("shared_modules") == 13 and summary.get("metabolic_compatible_profiles") == 2643 and summary.get("metabolic_missing_profile") == "K18513" and summary.get("dram_full_kofam_profiles") == 27754 and summary.get("ko_comparison_universe") == 2678, summary)
    audit.add("Summary", "claim-boundary", summary.get("gene_hit_is_activity") is False and summary.get("tool_output_is_flux") is False, summary)
    audit.add("Summary", "determinism", summary.get("truncation_seed") == 59002 and summary.get("random_process_after_input_preparation") is False, summary)
    run_contract = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    audit.add(
        "Contract",
        "globally-unique-protein-ids",
        run_contract.get("metabolic_protein_ids_globally_prefixed") is True,
        run_contract,
    )

    audit_figures(figures, audit, FIGURES)
    audit_chapter(
        chapter,
        audit,
        article=59,
        figure_stems=FIGURES,
        tokens=(
            "DRAM 1.5.0",
            "METABOLIC-G v4.0",
            "KOfam 2026-06-01",
            "K18513",
            "DRAM2",
            "beta",
            "PYTHONNOUSERSITE",
            "CheckM2",
            "GTDB-Tk R232",
            "0.75",
            "gene hit",
            "metabolic potential",
            "flux",
        ),
    )
    plot_python = args.plot_python.resolve()
    audit.add("Environment", "plot-python", plot_python.is_file(), str(plot_python))
    if plot_python.is_file():
        stage_reanalysis(root, frozen, stage, audit, plot_python)
    return finish(
        article=59,
        audit=audit,
        output=output,
        payload={
            "primary_real_mags": 24,
            "sensitivity_genomes": 4,
            "kofam_release": "2026-06-01",
            "shared_modules": 13,
            "figure_families": len(FIGURES),
            "frozen_payload_files": len(list(frozen.rglob("*"))),
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
