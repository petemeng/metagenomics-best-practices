#!/usr/bin/env python3
"""Fail-closed compatibility audit for the Article 61 SMETANA branch."""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import resource
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from reframed import Environment, FBA, set_default_solver
from smetana.interface import build_cache
from smetana.legacy import Community


SEED = 61001
GLOBAL_STEP = "SMETANA global MIP/MRO"
DETAILED_STEP = "SMETANA detailed SCS/MUS/MPS"
POSITIVE = 1e-9


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def complete_constraints(model) -> dict[str, tuple[float, float]]:
    return Environment.complete(model, max_uptake=1000).apply(
        model, inplace=False, warning=False
    )


def objective_value(model, objective: dict[str, float] | None = None) -> tuple[str, float]:
    solution = FBA(
        model,
        objective=objective,
        constraints=complete_constraints(model),
    )
    value = float(solution.fobj) if solution.fobj is not None else math.nan
    return str(solution.status).split(".")[-1], value


def member_growth(community: Community) -> dict[str, tuple[str, float]]:
    constraints = complete_constraints(community.merged)
    results: dict[str, tuple[str, float]] = {}
    for model_id, biomass_id in community.organisms_biomass_reactions.items():
        solution = FBA(
            community.merged,
            objective={biomass_id: 1.0},
            constraints=constraints,
        )
        value = float(solution.fobj) if solution.fobj is not None else math.nan
        results[model_id] = (str(solution.status).split(".")[-1], value)
    return results


def positive_count(frame: pd.DataFrame, column: str) -> int:
    return int(pd.to_numeric(frame[column], errors="coerce").fillna(0).gt(0).sum())


def main() -> None:
    args = parse_args()
    random.seed(SEED)
    np.random.seed(SEED)
    for variable in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[variable] = "1"
    set_default_solver("scip")

    work = args.work_dir.resolve()
    subcommunity = pd.read_csv(work / "smetana-subcommunity.tsv", sep="\t")
    models = subcommunity.sort_values("Rank").SBML.tolist()
    if len(models) != 6 or len(set(models)) != 6:
        raise RuntimeError("Expected six uniquely named SMETANA SBML models")

    medium = pd.read_csv(work / "smetana-media.tsv", sep="\t")
    if len(medium) != 159 or medium.compound.duplicated().any():
        raise RuntimeError("Expected 159 declared anoxic medium compounds")

    source_global = work / "smetana/western_global.tsv"
    source_detailed = work / "smetana/western_detailed.tsv"
    source_log = work / "smetana/run.log"
    for path in (source_global, source_detailed, source_log):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    global_frame = pd.read_csv(source_global, sep="\t", keep_default_na=False)
    detailed = pd.read_csv(source_detailed, sep="\t")
    if len(global_frame) != 1 or len(detailed) != 7853:
        raise RuntimeError(
            f"Locked SMETANA output dimensions changed: global={len(global_frame)}, "
            f"detailed={len(detailed)}"
        )
    finite_global = global_frame[["mip", "mro"]].apply(
        pd.to_numeric, errors="coerce"
    ).notna().any(axis=None)
    if finite_global:
        raise RuntimeError("Expected the audited global MIP/MRO output to be non-estimable")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        cache = build_cache(models, flavor="fbc2")
        model_objects = [
            cache.get_model(model_id, reset_id=True) for model_id in cache.get_ids()
        ]
    if {model.id for model in model_objects} != set(subcommunity.ModelID):
        raise RuntimeError("Loaded SMETANA model identifiers changed")

    interacting = Community(
        "all", model_objects, copy_models=False, interacting=True
    )
    pooled_reactions = set(interacting.merged.reactions)
    compatibility = medium.copy()
    compatibility["PoolExchange"] = compatibility.compound.map(
        lambda compound: f"R_EX_M_{compound}_e_pool"
    )
    compatibility["PresentInSubcommunity"] = compatibility.PoolExchange.isin(
        pooled_reactions
    )
    matched_medium = int(compatibility.PresentInSubcommunity.sum())
    if matched_medium != 120:
        raise RuntimeError(
            f"Locked pooled-exchange match changed: {matched_medium}/159"
        )
    write_tsv(work / "smetana-medium-compatibility.tsv", compatibility)

    standalone: dict[str, tuple[str, float]] = {}
    for model in model_objects:
        standalone[model.id] = objective_value(model)
    noninteracting = Community(
        "all", model_objects, copy_models=False, interacting=False
    )
    interacting_growth = member_growth(interacting)
    noninteracting_growth = member_growth(noninteracting)

    rows: list[dict[str, object]] = []
    for model_id in subcommunity.sort_values("Rank").ModelID:
        standalone_status, standalone_value = standalone[model_id]
        interacting_status, interacting_value = interacting_growth[model_id]
        noninteracting_status, noninteracting_value = noninteracting_growth[model_id]
        rows.append(
            {
                "ModelID": model_id,
                "StandaloneStatus": standalone_status,
                "StandaloneCompleteGrowth": standalone_value,
                "InteractingStatus": interacting_status,
                "InteractingCompleteGrowth": interacting_value,
                "NoninteractingStatus": noninteracting_status,
                "NoninteractingCompleteGrowth": noninteracting_value,
                "StandalonePositive": standalone_value > POSITIVE,
                "InteractingPositive": interacting_value > POSITIVE,
                "NoninteractingPositive": noninteracting_value > POSITIVE,
            }
        )
    audit = pd.DataFrame(rows)
    if not audit.StandalonePositive.all() or not audit.InteractingPositive.all():
        raise RuntimeError("A model failed the complete-environment positive-growth control")
    if audit.NoninteractingPositive.any():
        raise RuntimeError("Locked legacy non-interacting diagnostic unexpectedly changed")

    write_tsv(work / "smetana-compatibility-audit.tsv", audit)
    shutil.copy2(source_global, work / "smetana-global.tsv")
    shutil.copy2(source_detailed, work / "smetana-detailed.tsv")
    shutil.copy2(source_log, work / "smetana-run.log")

    ledger_path = work / "run-ledger.tsv"
    ledger = pd.read_csv(ledger_path, sep="\t")
    previous = ledger.set_index("Step").to_dict("index")
    ledger = ledger.loc[~ledger.Step.isin({GLOBAL_STEP, DETAILED_STEP})].copy()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    global_elapsed = float(previous.get(GLOBAL_STEP, {}).get("ElapsedSeconds", 0.0))
    detailed_elapsed = float(previous.get(DETAILED_STEP, {}).get("ElapsedSeconds", 0.0))
    global_rss = int(previous.get(GLOBAL_STEP, {}).get("CumulativeMaxRSSKB", rss))
    detailed_rss = int(previous.get(DETAILED_STEP, {}).get("CumulativeMaxRSSKB", rss))
    limitation_rows = pd.DataFrame(
        [
            {
                "Step": GLOBAL_STEP,
                "Status": "not_estimable",
                "ElapsedSeconds": global_elapsed,
                "CumulativeMaxRSSKB": global_rss,
                "Detail": (
                    "MIP: SMETANA 1.2.1 legacy non-interacting merge has zero "
                    "complete-medium growth for all 6 exported AGORA2 SBML members; "
                    "MRO: member minimal-medium solve failed for Bacteroides_vulgatus."
                ),
            },
            {
                "Step": DETAILED_STEP,
                "Status": "passed_with_limitation",
                "ElapsedSeconds": detailed_elapsed,
                "CumulativeMaxRSSKB": detailed_rss,
                "Detail": (
                    "7,853 rows returned; SCS and composite SMETANA were zero for "
                    "all rows, while MUS/MPS component outputs were non-zero. Do not "
                    "interpret zero composite scores as evidence that exchange is absent."
                ),
            },
        ]
    )
    ledger = pd.concat([ledger, limitation_rows], ignore_index=True)
    if len(ledger) != 10 or int(ledger.Status.eq("passed").sum()) != 8:
        raise RuntimeError("Expected eight passing MICOM steps plus two audited SMETANA rows")
    write_tsv(ledger_path, ledger)

    score_counts = {
        "positive_scs_rows": positive_count(detailed, "scs"),
        "positive_mus_rows": positive_count(detailed, "mus"),
        "positive_mps_rows": positive_count(detailed, "mps"),
        "positive_smetana_rows": positive_count(detailed, "smetana"),
    }
    expected_counts = {
        "positive_scs_rows": 0,
        "positive_mus_rows": 1330,
        "positive_mps_rows": 395,
        "positive_smetana_rows": 0,
    }
    if score_counts != expected_counts:
        raise RuntimeError(f"Locked detailed-score counts changed: {score_counts}")

    result = {
        "article": 61,
        "seed": SEED,
        "declared_medium_compounds": len(compatibility),
        "matched_pooled_exchanges": matched_medium,
        "global_rows": len(global_frame),
        "global_mip_status": "not_estimable_legacy_noninteracting_merge",
        "global_mro_status": "not_estimable_member_minimal_medium",
        "detailed_rows": len(detailed),
        **score_counts,
        "standalone_positive_models": int(audit.StandalonePositive.sum()),
        "interacting_positive_models": int(audit.InteractingPositive.sum()),
        "legacy_noninteracting_positive_models": int(audit.NoninteractingPositive.sum()),
        "interpretation": (
            "software_model_interface_limitation_not_biological_absence"
        ),
    }
    (work / "smetana-compatibility-summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
