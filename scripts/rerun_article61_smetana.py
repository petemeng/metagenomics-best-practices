#!/usr/bin/env python3
"""Re-run only the SMETANA branch after validating SBML/medium compatibility."""

from __future__ import annotations

import argparse
import contextlib
import json
import os
from pathlib import Path
import random
import resource
import time

import numpy as np
import pandas as pd
from reframed import set_default_solver
from smetana.interface import build_cache, main as run_smetana
from smetana.legacy import Community


SEED = 61001
GLOBAL_STEP = "SMETANA global MIP/MRO"
DETAILED_STEP = "SMETANA detailed SCS/MUS/MPS"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def timed_step(name: str, function) -> dict[str, object]:
    started = time.perf_counter()
    status = "passed"
    detail = ""
    try:
        function()
    except Exception as error:
        status = "failed"
        detail = f"{type(error).__name__}: {error}"
        raise
    finally:
        row = {
            "Step": name,
            "Status": status,
            "ElapsedSeconds": round(time.perf_counter() - started, 6),
            "CumulativeMaxRSSKB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "Detail": detail,
        }
    return row


def main() -> None:
    args = parse_args()
    random.seed(SEED)
    np.random.seed(SEED)
    for variable in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[variable] = "1"

    work = args.work_dir.resolve()
    subcommunity = pd.read_csv(work / "smetana-subcommunity.tsv", sep="\t")
    models = subcommunity.sort_values("Rank").SBML.tolist()
    if len(models) != 6 or len(set(models)) != 6:
        raise RuntimeError("SMETANA recovery requires six unique SBML models")
    medium = pd.read_csv(work / "smetana-media.tsv", sep="\t")
    if len(medium) != 159 or medium.compound.duplicated().any():
        raise RuntimeError("Unexpected SMETANA medium membership")

    cache = build_cache(models, flavor="fbc2")
    model_objects = [cache.get_model(model_id, reset_id=True) for model_id in cache.get_ids()]
    community = Community("all", model_objects, copy_models=False)
    pooled_reactions = set(community.merged.reactions)
    compatibility = medium.copy()
    compatibility["PoolExchange"] = compatibility.compound.map(
        lambda compound: f"R_EX_M_{compound}_e_pool"
    )
    compatibility["PresentInSubcommunity"] = compatibility.PoolExchange.isin(pooled_reactions)
    matched = int(compatibility.PresentInSubcommunity.sum())
    if matched < 100:
        raise RuntimeError(
            f"Only {matched}/159 declared medium compounds match pooled exchanges"
        )
    write_tsv(work / "smetana-medium-compatibility.tsv", compatibility)

    smetana_dir = work / "smetana"
    smetana_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = str(smetana_dir / "western")
    for suffix in ("global.tsv", "detailed.tsv"):
        path = Path(output_prefix + "_" + suffix)
        if path.exists():
            path.unlink()
    set_default_solver("scip")
    common = {
        "output": output_prefix,
        "flavor": "fbc2",
        "media": "western_diet_gut_anoxic_membership",
        "mediadb": str(work / "smetana-media.tsv"),
        "aerobic": False,
        "zeros": True,
        "verbose": True,
        "min_mol_weight": True,
        "use_lp": False,
    }
    rows: list[dict[str, object]] = []
    log_path = smetana_dir / "run.log"
    with log_path.open("w", encoding="utf-8") as log:
        def run_global() -> None:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                run_smetana(models, mode="global", **common)
                log.flush()

        def run_detailed() -> None:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                run_smetana(models, mode="detailed", ignore_coupling=False, **common)
                log.flush()

        rows.append(timed_step(GLOBAL_STEP, run_global))
        rows.append(timed_step(DETAILED_STEP, run_detailed))

    global_output = smetana_dir / "western_global.tsv"
    detailed_output = smetana_dir / "western_detailed.tsv"
    if not global_output.is_file() or not detailed_output.is_file():
        raise RuntimeError("SMETANA did not produce both expected tables")
    global_frame = pd.read_csv(global_output, sep="\t", keep_default_na=False)
    detailed_frame = pd.read_csv(detailed_output, sep="\t")
    finite_global = global_frame[["mip", "mro"]].apply(
        pd.to_numeric, errors="coerce"
    ).notna().any(axis=None)
    if detailed_frame.empty:
        raise RuntimeError("SMETANA detailed output remained empty")
    positive_scs = int(
        pd.to_numeric(detailed_frame.scs, errors="coerce").fillna(0).gt(0).sum()
    )
    positive_composite = int(
        pd.to_numeric(detailed_frame.smetana, errors="coerce").fillna(0).gt(0).sum()
    )
    if not finite_global:
        rows[0]["Status"] = "not_estimable"
        rows[0]["Detail"] = (
            "Global MIP/MRO returned no finite value; run "
            "audit_article61_smetana_compatibility.py before interpretation."
        )
    if positive_scs == 0 and positive_composite == 0:
        rows[1]["Status"] = "passed_with_limitation"
        rows[1]["Detail"] = (
            f"{len(detailed_frame):,} rows returned, but all SCS and composite "
            "SMETANA scores are zero; this is not evidence of biological absence."
        )
    write_tsv(work / "smetana-global.tsv", global_frame)
    write_tsv(work / "smetana-detailed.tsv", detailed_frame)

    ledger_path = work / "run-ledger.tsv"
    ledger = pd.read_csv(ledger_path, sep="\t")
    ledger = ledger.loc[~ledger.Step.isin({GLOBAL_STEP, DETAILED_STEP})].copy()
    ledger = pd.concat([ledger, pd.DataFrame(rows)], ignore_index=True)
    expected_status = {
        GLOBAL_STEP: "passed" if finite_global else "not_estimable",
        DETAILED_STEP: (
            "passed_with_limitation"
            if positive_scs == 0 and positive_composite == 0
            else "passed"
        ),
    }
    observed_status = ledger.set_index("Step").Status.to_dict()
    if len(ledger) != 10 or any(
        observed_status.get(step) != status for step, status in expected_status.items()
    ):
        raise RuntimeError("Recovered run ledger did not preserve audited SMETANA statuses")
    write_tsv(ledger_path, ledger)

    result = {
        "article": 61,
        "seed": SEED,
        "declared_medium_compounds": len(compatibility),
        "matched_pooled_exchanges": matched,
        "global_rows": len(global_frame),
        "global_estimable": bool(finite_global),
        "detailed_rows": len(detailed_frame),
        "positive_scs_rows": positive_scs,
        "positive_detailed_scores": positive_composite,
        "steps": rows,
    }
    (work / "smetana-recovery-audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
