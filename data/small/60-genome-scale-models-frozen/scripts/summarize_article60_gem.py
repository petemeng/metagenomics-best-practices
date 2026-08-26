#!/usr/bin/env python3
"""Summarize Article 60 GEM structure, gap filling, feasibility, and audits."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

from carveme.cli.carve import config, load_media_db, project_dir
from reframed import Environment, FBA, load_cbmodel


GROWTH_THRESHOLD = 1e-6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--gapseq-env", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def finite_number(value: float | None) -> str:
    if value is None:
        return "NA"
    if not math.isfinite(float(value)):
        return str(value)
    return f"{float(value):.9g}"


def reaction_signature(model) -> str:
    rows: list[object] = []
    for reaction in sorted(model.reactions.values(), key=lambda item: item.id):
        rows.append(
            [
                reaction.id,
                finite_number(reaction.lb),
                finite_number(reaction.ub),
                finite_number(reaction.objective),
                [[metabolite, finite_number(coefficient)] for metabolite, coefficient in sorted(reaction.stoichiometry.items())],
            ]
        )
    payload = json.dumps(rows, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(payload).hexdigest()


def reaction_has_gpr(reaction) -> bool:
    association = getattr(reaction, "gpr", None)
    return association is not None and bool(str(association).strip())


def reaction_class(reaction_id: str) -> str:
    if reaction_id.startswith("R_EX_"):
        return "Exchange"
    if reaction_id.startswith("R_DM_") or reaction_id.startswith("R_sink_"):
        return "Demand/sink"
    if "bio" in reaction_id.lower() or reaction_id == "Growth":
        return "Biomass"
    return "Internal"


def topological_dead_ends(model) -> int:
    producers: dict[str, bool] = defaultdict(bool)
    consumers: dict[str, bool] = defaultdict(bool)
    for reaction in model.reactions.values():
        forward = reaction.ub > 0
        reverse = reaction.lb < 0
        for metabolite, coefficient in reaction.stoichiometry.items():
            if (forward and coefficient > 0) or (reverse and coefficient < 0):
                producers[metabolite] = True
            if (forward and coefficient < 0) or (reverse and coefficient > 0):
                consumers[metabolite] = True
    count = 0
    for metabolite in model.metabolites.values():
        compartment = str(getattr(metabolite, "compartment", "")).lower()
        if compartment in {"e", "e0", "extracellular"}:
            continue
        if not producers[metabolite.id] or not consumers[metabolite.id]:
            count += 1
    return count


def fba_result(model, environment: Environment) -> tuple[str, float | None, bool]:
    constraints = environment.apply(model, exclusive=True, inplace=False, warning=False)
    solution = FBA(model, constraints=constraints)
    value = None if solution.fobj is None else float(solution.fobj)
    feasible_growth = value is not None and math.isfinite(value) and value > GROWTH_THRESHOLD
    return str(solution.status).removeprefix("Status."), value, feasible_growth


def gapseq_environment(path: Path) -> Environment:
    environment = Environment()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            reaction = f"R_EX_{row['compounds']}_e0"
            max_flux = float(row["maxFlux"])
            environment[reaction] = (-max_flux, math.inf)
    return environment


def empty_environment(model) -> Environment:
    return Environment.empty(model)


def model_paths(work: Path, genome: str, tool: str) -> dict[str, Path]:
    if tool == "gapseq":
        base = work / "gapseq" / genome
        return {
            "Draft": base / f"{genome}-draft.xml",
            "Gap-filled": base / "filled-permissive" / f"{genome}.xml",
        }
    base = work / "carveme" / genome
    return {
        "Draft": base / f"{genome}-draft.xml",
        "Gap-filled": base / f"{genome}-filled-LB.xml",
    }


def load_model(path: Path, tool: str):
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    # Both pipelines export SBML Level 3 with the FBC v2 objective extension.
    # Passing ``None`` makes reframed fall back to an invalid/legacy flavor and
    # silently drops the active objective from gapseq drafts, so the biomass
    # reaction must be parsed explicitly through the fbc2 reader for both tools.
    return load_cbmodel(str(path), flavor="fbc2")


def summarize_models(
    work: Path, ledger: list[dict[str, str]], gapseq_env: Path
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    dict[tuple[str, str, str], object],
]:
    carveme_media = load_media_db(project_dir + config.get("input", "media_library"))
    gapseq_media_dir = gapseq_env / "share/gapseq/dat/media"
    gapseq_profiles = {
        "SGB_008": ("autotrophic.csv", "highH2", "autotrophic hyperthermophile proxy"),
        "SGB_010": ("MM_anaerobic_CO2_H2.csv", "highH2", "official methanogen tutorial medium"),
        "SGB_018": ("meerwasser.csv", "none", "marine host-rich proxy for an obligate symbiont"),
    }
    common_allmed = gapseq_environment(gapseq_media_dir / "ALLmed.csv")
    common_minimal = gapseq_environment(gapseq_media_dir / "MM_glu.csv")
    models: dict[tuple[str, str, str], object] = {}
    summaries: list[dict[str, object]] = []
    growth_rows: list[dict[str, object]] = []
    gap_rows: list[dict[str, object]] = []
    for row in ledger:
        genome = row["Genome"]
        for tool in ("gapseq", "CarveMe"):
            for stage, path in model_paths(work, genome, tool).items():
                model = load_model(path, tool)
                models[(genome, tool, stage)] = model
                counts = defaultdict(int)
                gpr_count = 0
                for reaction in model.reactions.values():
                    counts[reaction_class(reaction.id)] += 1
                    gpr_count += reaction_has_gpr(reaction)
                summaries.append(
                    {
                        "Genome": genome,
                        "Species": row["Species"],
                        "Domain": row["Domain"],
                        "Tool": tool,
                        "Stage": stage,
                        "Reactions": len(model.reactions),
                        "InternalReactions": counts["Internal"],
                        "ExchangeReactions": counts["Exchange"],
                        "DemandSinkReactions": counts["Demand/sink"],
                        "Metabolites": len(model.metabolites),
                        "Genes": len(model.genes),
                        "GPRBackedReactions": gpr_count,
                        "NoGPRReactions": len(model.reactions) - gpr_count,
                        "TopologicalDeadEndMetabolites": topological_dead_ends(model),
                        "BiomassReaction": model.biomass_reaction or "NA",
                        "SBMLSHA256": digest(path),
                        "ReactionStructureSHA256": reaction_signature(model),
                    }
                )
                if tool == "CarveMe":
                    media = {
                        "Construction medium": (
                            "CarveMe built-in LB",
                            "none",
                            Environment.from_compounds(carveme_media["LB"]),
                        ),
                        "Minimal proxy": (
                            "CarveMe built-in M9",
                            "none",
                            Environment.from_compounds(carveme_media["M9"]),
                        ),
                    }
                else:
                    profile, condition, rationale = gapseq_profiles.get(
                        genome,
                        ("ALLmed.csv", "none", "common permissive bacterial reconstruction proxy"),
                    )
                    media = {
                        "Construction medium": (
                            f"gapseq {profile}; {rationale}",
                            condition,
                            gapseq_environment(gapseq_media_dir / profile),
                        ),
                        "Common rich proxy": (
                            "gapseq ALLmed.csv",
                            "none",
                            common_allmed,
                        ),
                        "Minimal proxy": (
                            "gapseq MM_glu.csv",
                            "none",
                            common_minimal,
                        ),
                    }
                for medium, (definition, condition, environment) in media.items():
                    status, objective, grows = fba_result(model, environment)
                    growth_rows.append(
                        {
                            "Genome": genome,
                            "Species": row["Species"],
                            "Tool": tool,
                            "Stage": stage,
                            "Medium": medium,
                            "MediumDefinition": definition,
                            "EnvironmentalConditionDuringGapfill": condition,
                            "SolverStatus": status,
                            "BiomassObjective": finite_number(objective),
                            "GrowthAbove1e-6": str(grows).lower(),
                            "Interpretation": "constraint feasibility; not measured growth",
                        }
                    )
                status, objective, grows = fba_result(model, empty_environment(model))
                growth_rows.append(
                    {
                        "Genome": genome,
                        "Species": row["Species"],
                        "Tool": tool,
                        "Stage": stage,
                        "Medium": "No uptake audit",
                        "MediumDefinition": "all exchange uptake disabled",
                        "EnvironmentalConditionDuringGapfill": "none",
                        "SolverStatus": status,
                        "BiomassObjective": finite_number(objective),
                        "GrowthAbove1e-6": str(grows).lower(),
                        "Interpretation": "positive biomass would flag a leak or energy-generating cycle",
                    }
                )

            draft = models[(genome, tool, "Draft")]
            filled = models[(genome, tool, "Gap-filled")]
            draft_ids = set(draft.reactions)
            filled_ids = set(filled.reactions)
            added = sorted(filled_ids - draft_ids)
            removed = sorted(draft_ids - filled_ids)
            added_gpr = sum(reaction_has_gpr(filled.reactions[reaction]) for reaction in added)
            internal_added = sum(reaction_class(reaction) == "Internal" for reaction in added)
            top_added = [
                f"{reaction}:{filled.reactions[reaction].name or 'NA'}" for reaction in added[:20]
            ]
            gap_rows.append(
                {
                    "Genome": genome,
                    "Species": row["Species"],
                    "Domain": row["Domain"],
                    "Tool": tool,
                    "DraftReactions": len(draft_ids),
                    "FilledReactions": len(filled_ids),
                    "AddedReactions": len(added),
                    "AddedInternalReactions": internal_added,
                    "AddedGPRBackedReactions": added_gpr,
                    "AddedWithoutGPR": len(added) - added_gpr,
                    "RemovedReactions": len(removed),
                    "AddedFractionOfFilledPct": f"{100 * len(added) / max(len(filled_ids), 1):.6f}",
                    "HighGapfillBurdenGt10Pct": str(len(added) / max(len(filled_ids), 1) > 0.10).lower(),
                    "TopAddedReactionIDsAndNames": ";".join(top_added),
                }
            )
    return summaries, growth_rows, gap_rows, models


def duplicate_control(models: dict[tuple[str, str, str], object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for tool in ("gapseq", "CarveMe"):
        for stage in ("Draft", "Gap-filled"):
            parent = models[("SGB_002", tool, stage)]
            control = models[("TRUNC_100", tool, stage)]
            parent_set = set(parent.reactions)
            control_set = set(control.reactions)
            union = parent_set | control_set
            rows.append(
                {
                    "Tool": tool,
                    "Stage": stage,
                    "Parent": "SGB_002",
                    "ExactControl": "TRUNC_100",
                    "ParentReactions": len(parent_set),
                    "ControlReactions": len(control_set),
                    "ReactionJaccard": f"{len(parent_set & control_set) / max(len(union), 1):.9f}",
                    "ReactionStructureHashEqual": str(
                        reaction_signature(parent) == reaction_signature(control)
                    ).lower(),
                    "Expected": "identical protein sequences; model structure should match",
                }
            )
    return rows


def truncation_table(
    input_rows: list[dict[str, str]],
    protein_rows: list[dict[str, str]],
    summaries: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
    growth_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    inputs = {row["Genome"]: row for row in input_rows}
    proteins = {row["Genome"]: row for row in protein_rows}
    summary_index = {
        (str(row["Genome"]), str(row["Tool"]), str(row["Stage"])): row
        for row in summaries
    }
    gap_index = {(str(row["Genome"]), str(row["Tool"])): row for row in gap_rows}
    growth_index = {
        (str(row["Genome"]), str(row["Tool"]), str(row["Stage"]), str(row["Medium"])): row
        for row in growth_rows
    }
    output: list[dict[str, object]] = []
    for genome in ("TRUNC_050", "TRUNC_070", "TRUNC_090", "TRUNC_100"):
        for tool in ("gapseq", "CarveMe"):
            draft = summary_index[(genome, tool, "Draft")]
            gap = gap_index[(genome, tool)]
            rich_name = "Construction medium"
            minimal_name = "Minimal proxy"
            rich = growth_index[(genome, tool, "Gap-filled", rich_name)]
            minimal = growth_index[(genome, tool, "Gap-filled", minimal_name)]
            output.append(
                {
                    "Genome": genome,
                    "Tool": tool,
                    "RetentionObservedPct": inputs[genome]["RetentionObservedPct"],
                    "Genes": proteins[genome]["Genes"],
                    "DraftReactions": draft["Reactions"],
                    "GapfillAddedReactions": gap["AddedReactions"],
                    "GapfillAddedFractionPct": gap["AddedFractionOfFilledPct"],
                    "FilledRichBiomass": rich["BiomassObjective"],
                    "FilledRichGrowth": rich["GrowthAbove1e-6"],
                    "FilledMinimalBiomass": minimal["BiomassObjective"],
                    "FilledMinimalGrowth": minimal["GrowthAbove1e-6"],
                }
            )
    return output


def parse_elapsed(value: str) -> float:
    parts = value.strip().split(":")
    try:
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        if len(parts) == 2:
            return int(parts[0]) * 60 + float(parts[1])
        return float(parts[0])
    except ValueError:
        return math.nan


def resource_usage(work: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    pattern = re.compile(r"^(.*?):\s*(.*)$")
    for tool in ("gapseq", "carveme"):
        for path in sorted((work / "logs" / tool).glob("*.time.txt")):
            values: dict[str, str] = {}
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                stripped = line.strip()
                # GNU time embeds colons in the elapsed-time field label
                # (``h:mm:ss or m:ss``); splitting at the first colon therefore
                # corrupts both the key and value.  Preserve that field under a
                # stable key and parse the value after the final ``: `` token.
                if stripped.startswith("Elapsed (wall clock) time"):
                    values["Elapsed (wall clock) time"] = stripped.rsplit(": ", 1)[-1]
                    continue
                match = pattern.match(stripped)
                if match:
                    values[match.group(1)] = match.group(2)
            match = re.match(r"(?:gapseq|CarveMe)-(.+)-((?:SGB|TRUNC)_.+)$", path.stem.removesuffix(".time"))
            if not match:
                continue
            elapsed_key = next(
                (key for key in values if key.startswith("Elapsed (wall clock) time")), ""
            )
            rows.append(
                {
                    "Genome": match.group(2),
                    "Tool": "CarveMe" if tool == "carveme" else "gapseq",
                    "Step": match.group(1),
                    "ElapsedSeconds": f"{parse_elapsed(values.get(elapsed_key, 'NA')):.3f}",
                    "MaximumResidentSetKB": values.get("Maximum resident set size (kbytes)", "NA"),
                    "UserSeconds": values.get("User time (seconds)", "NA"),
                    "SystemSeconds": values.get("System time (seconds)", "NA"),
                    "ExitStatus": values.get("Exit status", "NA"),
                }
            )
    if not rows:
        raise RuntimeError("No /usr/bin/time resource ledgers found")
    return rows


def evidence_ladder() -> list[dict[str, object]]:
    return [
        {
            "Level": 1,
            "Evidence": "Checksum-locked MAG and shared proteins",
            "SupportedClaim": "The exact reconstruction input is known",
            "UnsupportedClaim": "The MAG is complete or uncontaminated",
        },
        {
            "Level": 2,
            "Evidence": "Sequence-supported reaction with GPR",
            "SupportedClaim": "A homolog supports enzymatic potential",
            "UnsupportedClaim": "The reaction is active in the sampled community",
        },
        {
            "Level": 3,
            "Evidence": "Draft network topology",
            "SupportedClaim": "The encoded reaction set forms this model",
            "UnsupportedClaim": "Failure to grow proves auxotrophy",
        },
        {
            "Level": 4,
            "Evidence": "Gap-filled reaction without GPR",
            "SupportedClaim": "The optimization required a topology repair",
            "UnsupportedClaim": "The organism contains the missing enzyme",
        },
        {
            "Level": 5,
            "Evidence": "Positive FBA objective under a named medium",
            "SupportedClaim": "Growth is feasible under model assumptions",
            "UnsupportedClaim": "The organism grows at the predicted rate in vitro",
        },
        {
            "Level": 6,
            "Evidence": "Independent phenotype or flux validation",
            "SupportedClaim": "A specific prediction has external support",
            "UnsupportedClaim": "Every model flux is validated",
        },
    ]


def main() -> None:
    args = parse_args()
    args.project_root = args.project_root.resolve()
    args.work_dir = args.work_dir.resolve()
    args.gapseq_env = args.gapseq_env.resolve()
    if not (args.work_dir / ".article60-models-complete").is_file():
        raise FileNotFoundError("Run run_article60_gem.py successfully first")
    inputs = read_tsv(args.work_dir / "input-mag-ledger.tsv")
    proteins = read_tsv(args.work_dir / "protein-id-audit.tsv")
    summaries, growth, gaps, models = summarize_models(args.work_dir, inputs, args.gapseq_env)
    controls = duplicate_control(models)
    truncations = truncation_table(inputs, proteins, summaries, gaps, growth)
    resources = resource_usage(args.work_dir)
    write_tsv(args.work_dir / "model-structure-summary.tsv", summaries)
    write_tsv(args.work_dir / "medium-feasibility.tsv", growth)
    write_tsv(args.work_dir / "gapfill-burden.tsv", gaps)
    write_tsv(args.work_dir / "determinism-control.tsv", controls)
    write_tsv(args.work_dir / "truncation-sensitivity.tsv", truncations)
    write_tsv(args.work_dir / "resource-usage.tsv", resources)
    write_tsv(args.work_dir / "evidence-ladder.tsv", evidence_ladder())
    leak_flags = [
        row for row in growth if row["Medium"] == "No uptake audit" and row["GrowthAbove1e-6"] == "true"
    ]
    summary = {
        "article": 60,
        "input_genomes": len(inputs),
        "tools": 2,
        "models": len(summaries),
        "fba_audits": len(growth),
        "gapfill_comparisons": len(gaps),
        "exact_duplicate_controls": len(controls),
        "exact_duplicate_hash_matches": sum(
            row["ReactionStructureHashEqual"] == "true" for row in controls
        ),
        "no_uptake_growth_flags": len(leak_flags),
        "high_gapfill_burden_flags": sum(
            row["HighGapfillBurdenGt10Pct"] == "true" for row in gaps
        ),
        "growth_threshold": GROWTH_THRESHOLD,
        "cross_tool_reaction_id_overlap_compared": False,
        "reason": "gapseq uses ModelSEED identifiers while CarveMe uses BiGG identifiers",
    }
    (args.work_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
