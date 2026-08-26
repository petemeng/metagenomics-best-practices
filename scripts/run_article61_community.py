#!/usr/bin/env python3
"""Run the locked MICOM and SMETANA analyses for Article 61."""

from __future__ import annotations

import argparse
import contextlib
import importlib.metadata
import json
import os
import random
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd
from micom import load_pickle
from micom.workflows import build, grow, tradeoff
from micom.workflows.core import workflow
from micom.workflows.media import process_medium
from micom.workflows.results import GrowthResults, combine_results
from micom.logger import micom_console
from reframed import set_default_solver
from smetana.interface import main as run_smetana


SEED = 61001
CUTOFF = 0.001
TRADEOFF_GRID = tuple(round(x, 1) for x in np.arange(0.1, 1.0, 0.1))
TRADEOFF_BATCHES = ((0.1, 0.2), (0.3, 0.4), (0.5, 0.6), (0.7, 0.8, 0.9))
PRIMARY_TRADEOFF = 0.5
ATOL = 1e-6
RTOL = 1e-6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    return parser.parse_args()


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def package_version(name: str) -> str:
    return importlib.metadata.version(name)


class Ledger:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.rows: list[dict[str, object]] = []

    @contextlib.contextmanager
    def step(self, name: str):
        started = time.perf_counter()
        status = "passed"
        detail = ""
        try:
            yield
        except Exception as error:
            status = "failed"
            detail = f"{type(error).__name__}: {error}"
            raise
        finally:
            elapsed = time.perf_counter() - started
            self.rows.append(
                {
                    "Step": name,
                    "Status": status,
                    "ElapsedSeconds": round(elapsed, 6),
                    "CumulativeMaxRSSKB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    "Detail": detail,
                }
            )
            write_tsv(self.path, pd.DataFrame(self.rows))


def stable_sort(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    if any(column in frame.index.names for column in columns):
        frame = frame.reset_index(drop=True)
    return frame.sort_values(columns, kind="mergesort").reset_index(drop=True)


def _single_tradeoff_worker(args: list[object]) -> pd.DataFrame:
    """Version-safe one-value MICOM trade-off worker."""
    model_path, value, medium, atol, rtol, presolve = args
    community = load_pickle(model_path)
    exchange_ids = [reaction.id for reaction in community.exchanges]
    community.medium = medium[medium.index.isin(exchange_ids)]
    community.solver.configuration.presolve = presolve
    community.optimize(rtol=rtol, atol=atol)
    solution = community.cooperative_tradeoff(
        fraction=value,
        fluxes=False,
        pfba=False,
        atol=atol,
        rtol=rtol,
    )
    rates = solution.members.copy()
    rates["taxon"] = rates.index
    rates["tradeoff"] = value
    rates["sample_id"] = community.id
    return rates.loc[rates.taxon.ne("medium")]


def single_tradeoff(
    manifest: pd.DataFrame,
    model_folder: Path,
    medium: pd.DataFrame,
    value: float,
    threads: int,
) -> pd.DataFrame:
    """Run one trade-off value without MICOM 0.39.0's singleton wrapper bug."""
    samples = manifest.sample_id.astype(str).unique()
    paths = {
        sample: model_folder
        / str(manifest.loc[manifest.sample_id.eq(sample), "file"].iloc[0])
        for sample in samples
    }
    processed = process_medium(medium, samples)
    arguments = [
        [
            str(paths[sample]),
            value,
            processed.flux.loc[processed.sample_id.eq(sample)],
            ATOL,
            RTOL,
            True,
        ]
        for sample in samples
    ]
    results = workflow(
        _single_tradeoff_worker,
        arguments,
        threads=threads,
        description=f"Trade-off {value:.1f}",
        progress=False,
    )
    return pd.concat(results, ignore_index=True)


def validate_tradeoff_frame(
    frame: pd.DataFrame,
    manifest: pd.DataFrame,
    values: tuple[float, ...],
    path: Path,
) -> pd.DataFrame:
    expected_samples = set(manifest.sample_id.astype(str))
    required = {"sample_id", "taxon", "tradeoff", "growth_rate", "abundance"}
    if not required.issubset(frame.columns):
        raise RuntimeError(f"Incomplete MICOM checkpoint: {path}")
    if set(frame.sample_id.astype(str)) != expected_samples:
        raise RuntimeError(f"MICOM checkpoint has the wrong sample set: {path}")
    observed_values = set(np.round(frame.tradeoff.to_numpy(float), 6))
    expected_values = set(np.round(np.asarray(values, dtype=float), 6))
    if frame.empty or observed_values != expected_values:
        raise RuntimeError(f"MICOM checkpoint has the wrong trade-offs: {path}")
    return frame


def tradeoff_checkpoint(
    manifest: pd.DataFrame,
    model_folder: Path,
    medium: pd.DataFrame,
    value: float,
    threads: int,
    path: Path,
) -> pd.DataFrame:
    """Run one trade-off value and retain a restart-safe checkpoint."""
    if path.is_file():
        frame = pd.read_csv(path, sep="\t")
    else:
        frame = single_tradeoff(manifest, model_folder, medium, value, threads)
        frame = stable_sort(frame, ["sample_id", "taxon"])
        write_tsv(path, frame)
    return validate_tradeoff_frame(frame, manifest, (value,), path)


def tradeoff_batch_checkpoint(
    manifest: pd.DataFrame,
    model_folder: Path,
    medium: pd.DataFrame,
    values: tuple[float, ...],
    threads: int,
    path: Path,
) -> pd.DataFrame:
    """Run two or more trade-off values in one restart-safe batch."""
    if len(values) < 2:
        raise ValueError("A MICOM trade-off batch requires at least two values")
    if path.is_file():
        frame = pd.read_csv(path, sep="\t")
    else:
        frame = tradeoff(
            manifest,
            str(model_folder),
            medium,
            tradeoffs=values,
            threads=threads,
            atol=ATOL,
            rtol=RTOL,
            presolve=True,
        )
        frame = frame.loc[frame.tradeoff.notna()].copy()
        frame = stable_sort(frame, ["sample_id", "tradeoff", "taxon"])
        write_tsv(path, frame)
    return validate_tradeoff_frame(frame, manifest, values, path)


def primary_checkpoint(
    manifest_row: pd.DataFrame,
    model_folder: Path,
    medium: pd.DataFrame,
    path: Path,
) -> GrowthResults:
    """Run pFBA for one sample and retain a restart-safe ZIP checkpoint."""
    sample_id = str(manifest_row.sample_id.iloc[0])
    if path.is_file():
        result = GrowthResults.load(str(path))
    else:
        result = grow(
            manifest_row,
            str(model_folder),
            medium,
            tradeoff=PRIMARY_TRADEOFF,
            threads=1,
            strategy="pFBA",
            atol=ATOL,
            rtol=RTOL,
            presolve=True,
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(f".{path.name}.tmp")
        result.save(str(temporary))
        temporary.replace(path)
    if set(result.growth_rates.sample_id.astype(str)) != {sample_id}:
        raise RuntimeError(f"pFBA checkpoint has the wrong sample: {path}")
    if result.growth_rates.empty or result.exchanges.empty:
        raise RuntimeError(f"pFBA checkpoint is empty: {path}")
    return result


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    random.seed(SEED)
    np.random.seed(SEED)
    if os.environ.get("ARTICLE61_QUIET") == "1":
        micom_console.quiet = True
    for variable in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[variable] = "1"

    work = args.work_dir.resolve()
    database = Path((work / "database-path.txt").read_text(encoding="utf-8").strip())
    if not (database / "manifest.csv").is_file():
        raise FileNotFoundError(database / "manifest.csv")
    taxonomy = pd.read_csv(work / "micom-taxonomy.tsv", sep="\t")
    equal_taxonomy = pd.read_csv(work / "micom-taxonomy-equal.tsv", sep="\t")
    medium = pd.read_csv(work / "medium.tsv", sep="\t")[["reaction", "flux"]]
    ledger = Ledger(work / "run-ledger.tsv")

    versions = {
        package: package_version(package)
        for package in (
            "micom", "smetana", "cobra", "reframed", "optlang", "osqp",
            "highspy", "pyscipopt", "pandas", "numpy", "scipy",
        )
    }
    versions.update(
        {
            "seed": SEED,
            "micom_solver": "hybrid (OSQP QP + HiGHS LP)",
            "smetana_solver": "SCIP via PySCIPOpt",
            "threads": args.threads,
            "atol": ATOL,
            "rtol": RTOL,
            "primary_tradeoff": PRIMARY_TRADEOFF,
            "flux_strategy": "pFBA",
        }
    )
    (work / "software-versions.json").write_text(
        json.dumps(versions, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    model_folder = work / "micom/models-observed"
    with ledger.step("MICOM build: observed abundances"):
        observed_manifest = build(
            taxonomy,
            str(database),
            str(model_folder),
            cutoff=CUTOFF,
            threads=args.threads,
            solver="hybrid",
        )
        write_tsv(work / "micom-observed-manifest.tsv", observed_manifest)

    with ledger.step("MICOM cooperative-tradeoff grid"):
        tradeoff_dir = work / "micom/checkpoints/tradeoff"
        tradeoff_rates = pd.concat(
            [
                tradeoff_batch_checkpoint(
                    observed_manifest,
                    model_folder,
                    medium,
                    values,
                    args.threads,
                    tradeoff_dir
                    / f"tradeoff-{values[0]:.1f}-{values[-1]:.1f}.tsv",
                )
                for values in TRADEOFF_BATCHES
            ],
            ignore_index=True,
        )
        tradeoff_rates = stable_sort(tradeoff_rates, ["sample_id", "tradeoff", "taxon"])
        write_tsv(work / "micom-tradeoff-growth.tsv", tradeoff_rates)

    with ledger.step("MICOM primary pFBA"):
        primary_dir = work / "micom/checkpoints/primary-pfba"
        primary = combine_results(
            primary_checkpoint(
                observed_manifest.loc[observed_manifest.sample_id.eq(sample_id)].copy(),
                model_folder,
                medium,
                primary_dir / f"{sample_id}.zip",
            )
            for sample_id in observed_manifest.sample_id.astype(str)
        )
        primary.growth_rates = stable_sort(primary.growth_rates, ["sample_id", "taxon"])
        primary.exchanges = stable_sort(primary.exchanges, ["sample_id", "taxon", "reaction"])
        primary.annotations = stable_sort(primary.annotations, ["reaction"])
        primary.save(str(work / "micom-primary-results.zip"))
        write_tsv(work / "micom-primary-growth.tsv", primary.growth_rates)
        write_tsv(work / "micom-primary-exchanges.tsv", primary.exchanges)
        write_tsv(work / "micom-exchange-annotations.tsv", primary.annotations)

    medium_sensitivity: list[pd.DataFrame] = []
    for scale in (0.5, 1.0, 2.0):
        with ledger.step(f"MICOM medium sensitivity: {scale:g}x"):
            scaled = medium.copy()
            scaled["flux"] *= scale
            rates = tradeoff_checkpoint(
                observed_manifest,
                model_folder,
                scaled,
                PRIMARY_TRADEOFF,
                args.threads,
                work / "micom/checkpoints/medium" / f"medium-{scale:.1f}.tsv",
            )
            rates["MediumScale"] = scale
            medium_sensitivity.append(rates)
    write_tsv(
        work / "micom-medium-sensitivity.tsv",
        stable_sort(pd.concat(medium_sensitivity, ignore_index=True), ["sample_id", "MediumScale", "taxon"]),
    )

    equal_folder = work / "micom/models-equal"
    with ledger.step("MICOM build: equal-abundance sensitivity"):
        equal_manifest = build(
            equal_taxonomy,
            str(database),
            str(equal_folder),
            cutoff=CUTOFF,
            threads=args.threads,
            solver="hybrid",
        )
        write_tsv(work / "micom-equal-manifest.tsv", equal_manifest)
    with ledger.step("MICOM equal-abundance tradeoff"):
        equal_rates = tradeoff_checkpoint(
            equal_manifest,
            equal_folder,
            medium,
            PRIMARY_TRADEOFF,
            args.threads,
            work / "micom/checkpoints/equal" / "tradeoff-0.5.tsv",
        )
        write_tsv(
            work / "micom-equal-abundance-growth.tsv",
            stable_sort(equal_rates, ["sample_id", "taxon"]),
        )

    subcommunity = pd.read_csv(work / "smetana-subcommunity.tsv", sep="\t")
    models = subcommunity.sort_values("Rank").SBML.tolist()
    if len(models) != len(set(models)) or len(models) != 6:
        raise RuntimeError("SMETANA requires six uniquely named model files")
    media_db = str(work / "smetana-media.tsv")
    smetana_dir = work / "smetana"
    smetana_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = str(smetana_dir / "western")
    set_default_solver("scip")
    with (smetana_dir / "run.log").open("w", encoding="utf-8") as log:
        with ledger.step("SMETANA global MIP/MRO"):
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                run_smetana(
                    models,
                    mode="global",
                    output=output_prefix,
                    flavor="fbc2",
                    media="western_diet_gut_anoxic_membership",
                    mediadb=media_db,
                    aerobic=False,
                    zeros=True,
                    verbose=True,
                    min_mol_weight=True,
                    use_lp=False,
                )
                log.flush()
        with ledger.step("SMETANA detailed SCS/MUS/MPS"):
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                run_smetana(
                    models,
                    mode="detailed",
                    output=output_prefix,
                    flavor="fbc2",
                    media="western_diet_gut_anoxic_membership",
                    mediadb=media_db,
                    aerobic=False,
                    zeros=True,
                    verbose=True,
                    min_mol_weight=True,
                    use_lp=False,
                    ignore_coupling=False,
                )
                log.flush()

    global_output = smetana_dir / "western_global.tsv"
    detailed_output = smetana_dir / "western_detailed.tsv"
    if not global_output.is_file() or not detailed_output.is_file():
        raise RuntimeError("SMETANA did not produce both expected output tables")
    global_frame = pd.read_csv(global_output, sep="\t", keep_default_na=False)
    detailed_frame = pd.read_csv(detailed_output, sep="\t")
    if detailed_frame.empty:
        raise RuntimeError("SMETANA detailed output is empty")
    finite_global = global_frame[["mip", "mro"]].apply(
        pd.to_numeric, errors="coerce"
    ).notna().any(axis=None)
    positive_scs = int(
        pd.to_numeric(detailed_frame.scs, errors="coerce").fillna(0).gt(0).sum()
    )
    positive_composite = int(
        pd.to_numeric(detailed_frame.smetana, errors="coerce").fillna(0).gt(0).sum()
    )
    for row in ledger.rows:
        if row["Step"] == "SMETANA global MIP/MRO" and not finite_global:
            row["Status"] = "not_estimable"
            row["Detail"] = (
                "Global MIP/MRO returned no finite value; run the locked "
                "SMETANA compatibility audit before interpretation."
            )
        if (
            row["Step"] == "SMETANA detailed SCS/MUS/MPS"
            and positive_scs == 0
            and positive_composite == 0
        ):
            row["Status"] = "passed_with_limitation"
            row["Detail"] = (
                f"{len(detailed_frame):,} rows returned, but all SCS and composite "
                "SMETANA scores are zero; this is not evidence of biological absence."
            )
    write_tsv(work / "smetana-global.tsv", global_frame)
    write_tsv(work / "smetana-detailed.tsv", detailed_frame)
    write_tsv(ledger.path, pd.DataFrame(ledger.rows))
    print(json.dumps({"versions": versions, "steps": ledger.rows}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
