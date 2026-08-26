#!/usr/bin/env python3
"""Fail-closed validation for Article 60 gapseq/CarveMe evidence."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
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
    "60-model-size-gapfill",
    "60-gapfill-burden",
    "60-medium-feasibility",
    "60-truncation-sensitivity",
    "60-model-audit",
    "60-evidence-ladder",
)

SUMMARY_FILES = (
    "model-structure-summary.tsv",
    "medium-feasibility.tsv",
    "gapfill-burden.tsv",
    "determinism-control.tsv",
    "truncation-sensitivity.tsv",
    "resource-usage.tsv",
    "evidence-ladder.tsv",
    "summary.json",
)

EXPECTED_DATABASE_ASSETS = {
    "Bacteria.tar.gz": ("345520909", "ec669d568a2e5459bfe6619865a3cd2612fb9edc3a9688ad45fdf619d2e1afeb"),
    "Archaea.tar.gz": ("16004564", "62e3cf8ed9b37494f864f840bf3d5f2a6d47869d9daf56e512f3fce091df51c7"),
    "md5sums.txt": ("397", "6fd16815ec42806cadc5ac34a06e4009535d9c7af9e748a1b6f1967314bc0a9c"),
    "universe_bacteria.xml.gz": ("595770", "b983aec83c4f6c30ec58ded6684d15e90355e7d6c9d573747dde9d0bbae35d19"),
    "universe_archaea.xml.gz": ("616926", "e09da6a02ee95f2b1c1eff5fca0dc62c5b7e288ae8bf08078bcabd214fa38b2f"),
    "bigg_proteins.dmnd": ("12147703", "9f8f675dc43c1e18f040e76183be64440396760b180ff908ce39dc30e12440f7"),
    "media_db.tsv": ("8711", "fab9ebde307ec626473fa8f546d22227b28f32289ecd4e26ea79d8cf086b78a2"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stage-dir", type=Path, required=True)
    parser.add_argument("--gapseq-env", type=Path, required=True)
    parser.add_argument("--model-python", type=Path, required=True)
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
        payload = (
            normalized.size[0].to_bytes(8, "little")
            + normalized.size[1].to_bytes(8, "little")
            + normalized.tobytes()
        )
    return hashlib.sha256(payload).hexdigest()


def model_stage_path(stage: Path, genome: str, tool: str, model_stage: str) -> Path:
    if tool == "gapseq":
        base = stage / "gapseq" / genome
        if model_stage == "draft":
            return base / f"{genome}-draft.xml"
        return base / "filled-permissive" / f"{genome}.xml"
    base = stage / "carveme" / genome
    if model_stage == "draft":
        return base / f"{genome}-draft.xml"
    return base / f"{genome}-filled-LB.xml"


def gunzip_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rb") as source_handle, target.open("wb") as target_handle:
        shutil.copyfileobj(source_handle, target_handle, length=8 * 1024 * 1024)


def stage_reanalysis(
    root: Path,
    frozen: Path,
    stage: Path,
    gapseq_env: Path,
    model_python: Path,
    plot_python: Path,
    audit: Audit,
) -> None:
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    for name in (
        "input-mag-ledger.tsv",
        "protein-id-audit.tsv",
        "truncation-ledger.tsv",
    ):
        shutil.copy2(frozen / name, stage / name)
    for path in sorted((frozen / "logs").rglob("*.time.txt")):
        target = stage / "logs" / path.relative_to(frozen / "logs")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    ledger = read_tsv(frozen / "input-mag-ledger.tsv")
    for row in ledger:
        genome = row["Genome"]
        for tool in ("gapseq", "CarveMe"):
            for model_stage in ("draft", "gapfilled"):
                source = frozen / "models" / tool / genome / f"{model_stage}.xml.gz"
                gunzip_copy(source, model_stage_path(stage, genome, tool, model_stage))
    (stage / ".article60-models-complete").write_text(
        "staged frozen Article 60 models\n", encoding="utf-8"
    )

    environment = os.environ.copy()
    environment["MPLCONFIGDIR"] = str(stage / "matplotlib")
    summarized = subprocess.run(
        [
            str(model_python),
            str(root / "scripts/summarize_article60_gem.py"),
            "--project-root",
            str(root),
            "--work-dir",
            str(stage),
            "--gapseq-env",
            str(gapseq_env),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    audit.add(
        "Reanalysis",
        "summary-exit",
        summarized.returncode == 0,
        (summarized.stdout + summarized.stderr)[-3000:],
    )
    if summarized.returncode != 0:
        return
    for name in SUMMARY_FILES:
        audit.add(
            "Reanalysis",
            f"summary-hash-{name}",
            sha256(stage / name) == sha256(frozen / name),
            sha256(stage / name),
        )

    staged_figures = stage / "figures"
    plotted = subprocess.run(
        [
            str(plot_python),
            str(root / "scripts/plot_article60_gem.py"),
            "--summary-dir",
            str(stage),
            "--figure-dir",
            str(staged_figures),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    audit.add(
        "Reanalysis",
        "plot-exit",
        plotted.returncode == 0,
        (plotted.stdout + plotted.stderr)[-3000:],
    )
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
    gapseq_env = args.gapseq_env.resolve()
    model_python = args.model_python.resolve()
    plot_python = args.plot_python.resolve()
    audit = Audit()

    audit_checksums(frozen, audit)
    contract = json.loads((frozen / "frozen-contract.json").read_text(encoding="utf-8"))
    audit.add("Contract", "article", contract.get("article") == 60, contract)
    audit.add("Contract", "primary-real-mags", contract.get("primary_real_mags") == 8, contract)
    audit.add("Contract", "sensitivity-genomes", contract.get("deterministic_sensitivity_genomes") == 4, contract)
    audit.add("Contract", "protein-fastas", contract.get("shared_protein_fastas_included") == 12, contract)
    audit.add("Contract", "sbml-models", contract.get("sbml_models_included") == 48, contract)
    audit.add("Contract", "large-db-excluded", contract.get("large_databases_included") is False, contract)
    audit.add("Contract", "seed", contract.get("truncation_seed") == 59002, contract)

    manifest = read_tsv(frozen / "database-manifest.tsv")
    by_asset = {row["Asset"]: row for row in manifest}
    audit.add("Database", "manifest-assets", len(manifest) == 7 and len(by_asset) == 7, len(manifest))
    for asset, (size, digest) in EXPECTED_DATABASE_ASSETS.items():
        row = by_asset.get(asset, {})
        audit.add(
            "Database",
            f"locked-{asset}",
            row.get("ExpectedBytes") == size and row.get("SHA256") == digest,
            row,
        )
    database_audit = json.loads((frozen / "database-audit.json").read_text(encoding="utf-8"))
    audit.add("Database", "gapseq-release", database_audit.get("gapseq_release") == "2.1.0", database_audit)
    audit.add("Database", "sequence-db", database_audit.get("gapseq_sequence_db_version") == "1.5" and database_audit.get("gapseq_zenodo_record") == 20446806, database_audit)
    audit.add("Database", "nested-archives", database_audit.get("gapseq_nested_archives") == 6, database_audit)
    audit.add("Database", "extracted-fastas", database_audit.get("gapseq_extracted_fasta") == {"Archaea": 43672, "Bacteria": 38888}, database_audit)
    audit.add("Database", "carveme-release", database_audit.get("carveme_release") == "1.6.6" and database_audit.get("carveme_assets") == 4, database_audit)

    ledger = read_tsv(frozen / "input-mag-ledger.tsv")
    audit.add("Input", "ledger-count", len(ledger) == 12, len(ledger))
    audit.add("Input", "genome-ids-unique", len({row["Genome"] for row in ledger}) == 12, [row["Genome"] for row in ledger])
    audit.add("Input", "analysis-sets", Counter(row["AnalysisSet"] for row in ledger) == Counter({"Primary real MAG": 8, "Deterministic truncation sensitivity": 4}), Counter(row["AnalysisSet"] for row in ledger))
    audit.add("Input", "gtdb-release", all(row["GTDBRelease"] == "R232" for row in ledger), sorted({row["GTDBRelease"] for row in ledger}))
    primary = [row for row in ledger if row["AnalysisSet"] == "Primary real MAG"]
    truncations = [row for row in ledger if row["AnalysisSet"].startswith("Deterministic")]
    audit.add("Input", "primary-quality", all(row["CompletenessPct"] not in {"", "NA"} and row["ContaminationPct"] not in {"", "NA"} for row in primary), primary)
    audit.add("Input", "truncation-seed", all(row["DeterministicSeed"] == "59002" and row["ParentGenome"] == "SGB_002" for row in truncations), truncations)
    audit.add("Input", "truncation-targets", {row["RetentionTargetPct"] for row in truncations} == {"50", "70", "90", "100"}, truncations)
    parent_sequence = next(row["SequenceSetSHA256"] for row in primary if row["Genome"] == "SGB_002")
    control_sequence = next(row["SequenceSetSHA256"] for row in truncations if row["Genome"] == "TRUNC_100")
    audit.add("Input", "exact-sequence-control", parent_sequence == control_sequence, {"parent": parent_sequence, "control": control_sequence})

    proteins = read_tsv(frozen / "protein-id-audit.tsv")
    audit.add("Input", "protein-files", len(proteins) == 12, len(proteins))
    audit.add("Input", "protein-id-gates", all(as_bool(row["UniqueProteinIDs"]) and as_bool(row["GenomePrefixPass"]) for row in proteins), proteins)
    protein_by_genome = {row["Genome"]: row for row in proteins}
    audit.add("Input", "exact-protein-control", protein_by_genome["SGB_002"]["ProteinSequenceMultisetSHA256"] == protein_by_genome["TRUNC_100"]["ProteinSequenceMultisetSHA256"], protein_by_genome["TRUNC_100"])
    compressed_proteins = sorted((frozen / "inputs/proteins").glob("*.faa.gz"))
    audit.add("Input", "bundled-proteins", len(compressed_proteins) == 12 and all(path.stat().st_size > 0 for path in compressed_proteins), len(compressed_proteins))

    model_archives = sorted((frozen / "models").rglob("*.xml.gz"))
    audit.add("Model", "bundled-model-count", len(model_archives) == 48, len(model_archives))
    for path in model_archives:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            prefix = handle.read(2000).lower()
        audit.add("Model", f"sbml-{path.relative_to(frozen)}", "<sbml" in prefix, path.stat().st_size)

    structures = read_tsv(frozen / "model-structure-summary.tsv")
    triples = {(row["Genome"], row["Tool"], row["Stage"]) for row in structures}
    audit.add("Model", "structure-rows", len(structures) == 48 and len(triples) == 48, len(structures))
    audit.add("Model", "tool-stage-grid", {row["Tool"] for row in structures} == {"gapseq", "CarveMe"} and {row["Stage"] for row in structures} == {"Draft", "Gap-filled"}, sorted(triples))
    audit.add("Model", "positive-structure-counts", all(int(row["Reactions"]) > 0 and int(row["Metabolites"]) > 0 and int(row["Genes"]) >= 0 for row in structures), structures)
    audit.add("Model", "model-hashes", all(len(row["SBMLSHA256"]) == 64 and len(row["ReactionStructureSHA256"]) == 64 for row in structures), structures)

    gaps = read_tsv(frozen / "gapfill-burden.tsv")
    audit.add("Gapfill", "comparison-count", len(gaps) == 24, len(gaps))
    audit.add("Gapfill", "draft-preserved", all(int(row["RemovedReactions"]) == 0 for row in gaps), [row for row in gaps if int(row["RemovedReactions"]) != 0])
    added_by_tool = {
        tool: sum(int(row["AddedReactions"]) for row in gaps if row["Tool"] == tool)
        for tool in {row["Tool"] for row in gaps}
    }
    audit.add(
        "Gapfill",
        "repair-added",
        set(added_by_tool) == {"gapseq", "CarveMe"}
        and all(total > 0 for total in added_by_tool.values()),
        added_by_tool,
    )
    audit.add("Gapfill", "added-reactions-audited", all(int(row["AddedWithoutGPR"]) + int(row["AddedGPRBackedReactions"]) == int(row["AddedReactions"]) for row in gaps), gaps)

    feasibility = read_tsv(frozen / "medium-feasibility.tsv")
    audit.add("FBA", "audit-count", len(feasibility) == 168, len(feasibility))
    no_uptake = [row for row in feasibility if row["Medium"] == "No uptake audit"]
    audit.add("FBA", "no-uptake-count", len(no_uptake) == 48, len(no_uptake))
    audit.add("FBA", "no-uptake-no-growth", all(not as_bool(row["GrowthAbove1e-6"]) for row in no_uptake), [row for row in no_uptake if as_bool(row["GrowthAbove1e-6"])])
    reconstruction_media = {"gapseq": "Construction medium", "CarveMe": "Construction medium"}
    repaired_growth = [
        row for row in feasibility
        if row["Stage"] == "Gap-filled" and row["Medium"] == reconstruction_media[row["Tool"]]
    ]
    audit.add("FBA", "reconstruction-medium-grid", len(repaired_growth) == 24, len(repaired_growth))
    audit.add("FBA", "reconstruction-medium-feasible", all(as_bool(row["GrowthAbove1e-6"]) for row in repaired_growth), [row for row in repaired_growth if not as_bool(row["GrowthAbove1e-6"])])

    controls = read_tsv(frozen / "determinism-control.tsv")
    audit.add("Determinism", "control-grid", len(controls) == 4, controls)
    audit.add("Determinism", "reaction-jaccard", all(near(row["ReactionJaccard"], 1.0) for row in controls), controls)
    audit.add("Determinism", "structure-hash", all(as_bool(row["ReactionStructureHashEqual"]) for row in controls), controls)
    sensitivity = read_tsv(frozen / "truncation-sensitivity.tsv")
    audit.add("Determinism", "sensitivity-grid", len(sensitivity) == 8 and {row["Genome"] for row in sensitivity} == {"TRUNC_050", "TRUNC_070", "TRUNC_090", "TRUNC_100"}, sensitivity)

    summary = json.loads((frozen / "summary.json").read_text(encoding="utf-8"))
    audit.add("Summary", "article", summary.get("article") == 60, summary)
    audit.add("Summary", "counts", summary.get("input_genomes") == 12 and summary.get("models") == 48 and summary.get("fba_audits") == 168 and summary.get("gapfill_comparisons") == 24, summary)
    audit.add("Summary", "exact-controls", summary.get("exact_duplicate_hash_matches") == 4, summary)
    audit.add("Summary", "no-leaks", summary.get("no_uptake_growth_flags") == 0, summary)
    audit.add("Summary", "namespace-rule", summary.get("cross_tool_reaction_id_overlap_compared") is False, summary)

    run_contract = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    audit.add("Run", "no-production-failures", run_contract.get("failures") == [], run_contract)
    audit.add("Run", "solver-locks", run_contract.get("gapseq_solver") == "GLPK 5.0" and run_contract.get("carveme_solver") == "SCIP 10.0.3", run_contract)
    audit.add("Run", "fixed-draft-gapfill", run_contract.get("carveme_gapfill_mode") == "independent gapfill CLI applied to the fixed draft model", run_contract)
    resources = read_tsv(frozen / "resource-usage.tsv")
    expected_steps = {
        (row["Genome"], "gapseq", step)
        for row in ledger
        for step in (
            "find",
            "transport",
            "draft",
            "gapfill-construction-medium" if row["Domain"] == "Archaea" else "gapfill-permissive",
        )
    } | {
        (genome, "CarveMe", step)
        for genome in {row["Genome"] for row in ledger}
        for step in ("diamond", "draft", "gapfill-rich")
    }
    observed_success = {
        (row["Genome"], row["Tool"], row["Step"])
        for row in resources if row["ExitStatus"] == "0"
    }
    audit.add("Run", "resource-ledger-grid", expected_steps <= observed_success, sorted(expected_steps - observed_success))

    audit_figures(figures, audit, FIGURES)
    plot_source = (root / "scripts/plot_article60_gem.py").read_text(encoding="utf-8")
    audit.add("Figure", "english-only-source", re.search(r"[\u4e00-\u9fff]", plot_source) is None, "plot source contains no CJK text")
    audit_chapter(
        chapter,
        audit,
        article=60,
        figure_stems=FIGURES,
        tokens=(
            "gapseq 2.1.0",
            "CarveMe 1.6.6",
            "GPR",
            "ALLmed.csv",
            "SCIP 10.0.3",
            "gap filling",
        ),
    )
    stage_reanalysis(
        root,
        frozen,
        stage,
        gapseq_env,
        model_python,
        plot_python,
        audit,
    )
    return finish(
        article=60,
        audit=audit,
        output=output,
        payload={
            "inputs": len(ledger),
            "models": len(structures),
            "figures": len(FIGURES),
            "fba_audits": len(feasibility),
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
