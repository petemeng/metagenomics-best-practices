#!/usr/bin/env python3
"""Run the locked community- and MAG-level C/N/S audit for Article 62."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import json
import os
import random
import re
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


SEED = 62001
PERMUTATIONS = 9999
BOOTSTRAPS = 2000
MIN_NONZERO_SPRINGS = 12
KO_PATTERN = re.compile(r"K\d{5}")
FORMULA_PATTERN = re.compile(r"^[K0-9+()*/.\-\s]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, default=Path("db/element-cycle-cache"))
    parser.add_argument(
        "--rules", type=Path, default=Path("data/small/62-element-cycle-marker-rules.tsv")
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def digest(path: Path, algorithm: str = "sha256") -> str:
    checksum = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def write_tsv(path: Path, frame: pd.DataFrame, compress: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    target = path.with_suffix(path.suffix + ".gz") if compress else path
    frame.to_csv(
        target,
        sep="\t",
        index=False,
        lineterminator="\n",
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0}
        if compress
        else None,
    )


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
            self.rows.append(
                {
                    "Step": name,
                    "Status": status,
                    "ElapsedSeconds": round(time.perf_counter() - started, 6),
                    "CumulativeMaxRSSKB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    "Detail": detail,
                }
            )
            write_tsv(self.path, pd.DataFrame(self.rows))


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def rng_for(label: str) -> np.random.Generator:
    """Create a stable random stream that does not depend on loop order."""
    payload = hashlib.sha256(f"{SEED}|{label}".encode("utf-8")).digest()
    return np.random.default_rng(int.from_bytes(payload[:8], "little"))


def read_detected_table(path: Path, id_column: str, required_kos: list[str]) -> tuple[pd.DataFrame, list[str]]:
    header = pd.read_csv(path, sep="\t", comment="#", nrows=0)
    available = set(header.columns)
    if id_column not in available:
        raise RuntimeError(f"{path.name} lacks identifier column {id_column}")
    selected = [id_column] + sorted(set(required_kos) & available)
    frame = pd.read_csv(path, sep="\t", comment="#", usecols=selected)
    missing = sorted(set(required_kos) - available)
    for ko in missing:
        frame[ko] = 0.0
    return frame[[id_column] + sorted(required_kos)], missing


def evaluate_formula(frame: pd.DataFrame, formula: str) -> np.ndarray:
    if FORMULA_PATTERN.fullmatch(formula) is None:
        raise RuntimeError(f"Unsafe formula: {formula}")
    names = sorted(set(KO_PATTERN.findall(formula)))
    local = {name: frame[name].to_numpy(dtype=float) for name in names}
    value = eval(formula, {"__builtins__": {}}, local)  # noqa: S307 - grammar is locked above
    result = np.asarray(value, dtype=float)
    if result.shape != (len(frame),) or not np.isfinite(result).all() or (result < 0).any():
        raise RuntimeError(f"Invalid values produced by formula: {formula}")
    return result


def parse_dnf(text: str) -> list[list[str]]:
    groups = [group.split(",") for group in text.split(";")]
    if not groups or any(not group for group in groups):
        raise RuntimeError(f"Invalid DNF rule: {text}")
    if any(KO_PATTERN.fullmatch(ko) is None for group in groups for ko in group):
        raise RuntimeError(f"Invalid KO in DNF rule: {text}")
    return groups


def carrier_evidence(
    frame: pd.DataFrame, rule: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    alternatives = parse_dnf(rule)
    fractions = []
    observed = []
    totals = []
    strict = []
    for group in alternatives:
        presence = frame[group].to_numpy(dtype=float) > 0
        fractions.append(presence.mean(axis=1))
        observed.append(presence.sum(axis=1))
        totals.append(np.full(len(frame), len(group), dtype=int))
        strict.append(presence.all(axis=1))
    fraction_matrix = np.column_stack(fractions)
    best = fraction_matrix.argmax(axis=1)
    rows = np.arange(len(frame))
    return (
        np.column_stack(strict).any(axis=1),
        fraction_matrix[rows, best],
        np.column_stack(observed)[rows, best].astype(int),
        np.column_stack(totals)[rows, best].astype(int),
    )


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    result = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    p = values[finite]
    if len(p) == 0:
        return result
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    result[finite] = adjusted
    return result


def pearson_rows(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x_centered = x - x.mean(axis=1, keepdims=True)
    y_centered = y - y.mean(axis=1, keepdims=True)
    numerator = np.sum(x_centered * y_centered, axis=1)
    denominator = np.sqrt(
        np.sum(x_centered**2, axis=1) * np.sum(y_centered**2, axis=1)
    )
    return np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan, dtype=float),
        where=denominator > 0,
    )


def association(
    x: np.ndarray,
    y: np.ndarray,
    rng: np.random.Generator,
    permutations: int = PERMUTATIONS,
    bootstraps: int = BOOTSTRAPS,
) -> tuple[float, float, float, float, int]:
    keep = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x[keep], dtype=float)
    y = np.asarray(y[keep], dtype=float)
    n = len(x)
    if n < 8 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan, np.nan, np.nan, np.nan, n
    xr = rankdata(x)
    yr = rankdata(y)
    observed = pearson_rows(xr[None, :], yr[None, :])[0]
    permuted = np.vstack([rng.permutation(yr) for _ in range(permutations)])
    null = pearson_rows(np.broadcast_to(xr, permuted.shape), permuted)
    pvalue = (1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1)
    indices = rng.integers(0, n, size=(bootstraps, n))
    bootstrap_xr = rankdata(x[indices], axis=1)
    bootstrap_yr = rankdata(y[indices], axis=1)
    boot = pearson_rows(bootstrap_xr, bootstrap_yr)
    boot = boot[np.isfinite(boot)]
    low, high = np.quantile(boot, [0.025, 0.975]) if len(boot) else (np.nan, np.nan)
    return float(observed), float(pvalue), float(low), float(high), n


def parse_biom(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("matrix_type") != "sparse" or payload.get("shape") != [780, 500]:
        raise RuntimeError(f"Unexpected BIOM contract: {payload.get('shape')}")
    row_ids = [row["id"] for row in payload["rows"]]
    column_ids = [column["id"] for column in payload["columns"]]
    matrix = np.zeros(payload["shape"], dtype=float)
    for row, column, value in payload["data"]:
        matrix[int(row), int(column)] = float(value)
    frame = pd.DataFrame(matrix, index=row_ids, columns=column_ids)
    audit = {
        "shape": payload["shape"],
        "matrix_type": payload["matrix_type"],
        "nonzero_entries": int(np.count_nonzero(matrix)),
        "minimum_column_sum": float(matrix.sum(axis=0).min()),
        "maximum_column_sum": float(matrix.sum(axis=0).max()),
    }
    return frame, audit


def phylum_from_taxonomy(value: object) -> str:
    parts = [part.strip() for part in str(value).split(";") if part.strip()]
    return parts[1] if len(parts) > 1 else "Unclassified"


def main() -> None:
    args = parse_args()
    random.seed(SEED)
    np.random.seed(SEED)
    for variable in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[variable] = "1"

    cache = args.cache_dir.resolve()
    rules_path = args.rules.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    ledger = Ledger(out / "run-ledger.tsv")
    rules = pd.read_csv(rules_path, sep="\t", dtype=str)
    expected_elements = {"Carbon", "Carbon/Nitrogen", "Nitrogen", "Sulfur"}
    if len(rules) != 20 or set(rules.Element) != expected_elements or rules.ProcessID.duplicated().any():
        raise RuntimeError("The Article 62 rule table does not satisfy its 20-process contract")
    all_kos = sorted(
        set().union(
            *(set(KO_PATTERN.findall(row.CommunityFormula)) for row in rules.itertuples()),
            *(set(KO_PATTERN.findall(row.StrictCarrierDNF)) for row in rules.itertuples()),
        )
    )

    with ledger.step("Verify locked resources"):
        manifest = json.loads((cache / "download-manifest.json").read_text(encoding="utf-8"))
        if manifest["dataset_doi"] != "10.6084/m9.figshare.30284068.v2":
            raise RuntimeError("Unexpected dataset release")
        for record in manifest["resources"]:
            path = cache / record["local_file"]
            if not path.is_file() or path.stat().st_size != record["bytes"]:
                raise RuntimeError(f"Resource contract failed for {path}")
            if digest(path) != record["sha256"]:
                raise RuntimeError(f"SHA-256 failed for {path}")
        if digest(rules_path) != "f0e2761ac2220c401be16b6d3d2607ce30fc2eadb8bf34c5c0a1bf9536b28dca":
            raise RuntimeError("Marker rule table changed without contract update")

    with ledger.step("Read published KO and metadata tables"):
        community_source_kos = len(
            pd.read_csv(
                cache / "ko-proportions-in-metagenomes.tsv.gz",
                sep="\t",
                comment="#",
                nrows=0,
            ).columns
        ) - 1
        mag_source_kos = len(
            pd.read_csv(
                cache / "kos-in-mags.tsv.gz",
                sep="\t",
                comment="#",
                nrows=0,
            ).columns
        ) - 1
        if community_source_kos != 10_011 or mag_source_kos != 6_841:
            raise RuntimeError(
                "Published KO table widths changed: "
                f"community={community_source_kos}, MAG={mag_source_kos}"
            )
        community_kos, missing_community = read_detected_table(
            cache / "ko-proportions-in-metagenomes.tsv.gz", "sample", all_kos
        )
        mag_kos, missing_mag = read_detected_table(
            cache / "kos-in-mags.tsv.gz", "MAG", all_kos
        )
        metadata = pd.read_csv(
            cache / "sample-metadata.tsv",
            sep="\t",
            usecols=[
                "sample", "hotspring", "hotspring_common_name", "year", "pH",
                "temperature", "temperature_regime", "broad_region_short",
                "is_hotspring", "filtered_Nsequences", "metagenome2MAG_read_recruitment_rate",
            ],
        )
        mag_metadata = pd.read_csv(
            cache / "mag-metadata.tsv",
            sep="\t",
            usecols=[
                "MAG", "taxonomy", "Completeness", "Contamination", "Ncontigs",
                "total_contig_length", "GenBank_accession", "BioSample",
            ],
        )
        if len(community_kos) != 500 or len(metadata) != 500 or len(mag_kos) != 780:
            raise RuntimeError("Published table dimensions changed")
        if set(community_kos["sample"]) != set(metadata["sample"]):
            raise RuntimeError("Community KO samples and metadata do not match exactly")
        if set(mag_kos["MAG"]) != set(mag_metadata["MAG"]):
            raise RuntimeError("MAG KO and MAG metadata identifiers do not match exactly")
        metadata["pH"] = pd.to_numeric(metadata.pH, errors="coerce")
        metadata["temperature"] = pd.to_numeric(metadata.temperature, errors="coerce")
        metadata["metagenome2MAG_read_recruitment_rate"] = pd.to_numeric(
            metadata.metagenome2MAG_read_recruitment_rate, errors="coerce"
        )
        write_tsv(
            out / "missing-ko-columns.tsv",
            pd.concat(
                [
                    pd.DataFrame({"Table": "Community", "KO": missing_community}),
                    pd.DataFrame({"Table": "MAG", "KO": missing_mag}),
                ],
                ignore_index=True,
            ),
        )

    with ledger.step("Calculate community process indices"):
        community_rows: list[pd.DataFrame] = []
        for row in rules.itertuples(index=False):
            raw_index = evaluate_formula(community_kos, row.CommunityFormula)
            complete_gate = carrier_evidence(
                community_kos, row.StrictCarrierDNF
            )[0]
            community_rows.append(
                pd.DataFrame(
                    {
                        "SampleID": community_kos["sample"],
                        "ProcessID": row.ProcessID,
                        "RawCommunityIndex": raw_index,
                        "CommunityCompleteGate": complete_gate,
                        "CommunityIndex": np.where(complete_gate, raw_index, 0.0),
                    }
                )
            )
        community_long = pd.concat(community_rows, ignore_index=True).merge(
            rules[["ProcessID", "Element", "Process", "DirectionCaveat"]],
            on="ProcessID",
            how="left",
            validate="many_to_one",
        )
        community_long = community_long.merge(
            metadata.rename(columns={"sample": "SampleID"}),
            on="SampleID",
            how="left",
            validate="many_to_one",
        )
        write_tsv(out / "sample-process-index.tsv", community_long, compress=True)
        spring_process = (
            community_long.groupby(
                ["hotspring", "ProcessID", "Element", "Process"],
                observed=True,
                sort=True,
            )
            .agg(
                RawCommunityIndex=("RawCommunityIndex", "median"),
                CommunityIndex=("CommunityIndex", "median"),
                CompleteGateFraction=("CommunityCompleteGate", "mean"),
                Samples=("SampleID", "nunique"),
                Temperature=("temperature", "median"),
                pH=("pH", "median"),
            )
            .reset_index()
        )
        write_tsv(out / "spring-process-index.tsv", spring_process)
        regime = (
            community_long.groupby(
                ["temperature_regime", "ProcessID", "Element", "Process"],
                observed=True,
                sort=True,
            )
            .agg(
                RawMedianIndex=("RawCommunityIndex", "median"),
                MedianIndex=("CommunityIndex", "median"),
                CompleteGateFraction=("CommunityCompleteGate", "mean"),
                Q25=("CommunityIndex", lambda x: x.quantile(0.25)),
                Q75=("CommunityIndex", lambda x: x.quantile(0.75)),
                Samples=("SampleID", "nunique"),
                Springs=("hotspring", "nunique"),
            )
            .reset_index()
        )
        write_tsv(out / "temperature-regime-summary.tsv", regime)

    with ledger.step("Test spring-level environment associations"):
        association_rows: list[dict[str, object]] = []
        for variable in ("Temperature", "pH"):
            for row in rules.itertuples(index=False):
                frame = spring_process.loc[spring_process.ProcessID.eq(row.ProcessID)]
                nonzero_springs = int(frame.CommunityIndex.gt(0).sum())
                if nonzero_springs < MIN_NONZERO_SPRINGS:
                    rho, pvalue, low, high, n = (np.nan, np.nan, np.nan, np.nan, len(frame))
                else:
                    rho, pvalue, low, high, n = association(
                        frame.CommunityIndex.to_numpy(float),
                        frame[variable].to_numpy(float),
                        rng_for(f"environment|{variable}|{row.ProcessID}"),
                    )
                association_rows.append(
                    {
                        "ProcessID": row.ProcessID,
                        "Element": row.Element,
                        "Process": row.Process,
                        "Variable": variable,
                        "SpearmanRho": rho,
                        "PermutationP": pvalue,
                        "BootstrapLow95": low,
                        "BootstrapHigh95": high,
                        "IndependentSprings": n,
                        "NonzeroCommunitySprings": nonzero_springs,
                        "MinimumNonzeroSprings": MIN_NONZERO_SPRINGS,
                        "Permutations": PERMUTATIONS,
                        "Bootstraps": BOOTSTRAPS,
                    }
                )
        associations = pd.DataFrame(association_rows)
        associations["FDR"] = associations.groupby("Variable", sort=False)[
            "PermutationP"
        ].transform(lambda x: bh_adjust(x.to_numpy()))
        write_tsv(out / "environment-associations.tsv", associations)

    with ledger.step("Call strict and relaxed MAG carriers"):
        evidence_rows: list[pd.DataFrame] = []
        carrier_matrix = pd.DataFrame(index=mag_kos.MAG)
        for row in rules.itertuples(index=False):
            strict, completeness, observed, best_total = carrier_evidence(
                mag_kos, row.StrictCarrierDNF
            )
            max_markers = max(len(group) for group in parse_dnf(row.StrictCarrierDNF))
            relaxed = completeness >= 0.8 if max_markers >= 3 else strict.copy()
            carrier_matrix[row.ProcessID] = strict
            evidence_rows.append(
                pd.DataFrame(
                    {
                        "MAG": mag_kos.MAG,
                        "ProcessID": row.ProcessID,
                        "Element": row.Element,
                        "Process": row.Process,
                        "StrictCarrier": strict,
                        "Relaxed80Carrier": relaxed,
                        "BestMarkerCompleteness": completeness,
                        "ObservedMarkersInBestAlternative": observed,
                        "BestAlternativeMarkerCount": best_total,
                    }
                )
            )
        evidence = pd.concat(evidence_rows, ignore_index=True).merge(
            mag_metadata,
            on="MAG",
            how="left",
            validate="many_to_one",
        )
        evidence["Phylum"] = evidence.taxonomy.map(phylum_from_taxonomy)
        write_tsv(out / "mag-process-evidence.tsv", evidence, compress=True)
        carrier_summary = (
            evidence.groupby(["ProcessID", "Element", "Process"], sort=True)
            .agg(
                StrictCarrierMAGs=("StrictCarrier", "sum"),
                Relaxed80CarrierMAGs=("Relaxed80Carrier", "sum"),
                MedianCompleteness=(
                    "Completeness",
                    lambda x: float(np.nanmedian(x[evidence.loc[x.index, "StrictCarrier"]]))
                    if evidence.loc[x.index, "StrictCarrier"].any()
                    else np.nan,
                ),
                MedianContamination=(
                    "Contamination",
                    lambda x: float(np.nanmedian(x[evidence.loc[x.index, "StrictCarrier"]]))
                    if evidence.loc[x.index, "StrictCarrier"].any()
                    else np.nan,
                ),
                CarrierPhyla=(
                    "Phylum",
                    lambda x: x[evidence.loc[x.index, "StrictCarrier"]].nunique(),
                ),
            )
            .reset_index()
        )
        write_tsv(out / "mag-carrier-summary.tsv", carrier_summary)

    with ledger.step("Integrate MAG abundance and recovery ceiling"):
        abundance, biom_audit = parse_biom(cache / "mag-abundances-per-sample.biom")
        abundance = abundance.reindex(index=carrier_matrix.index, columns=metadata["sample"])
        if abundance.isna().any().any():
            raise RuntimeError("BIOM identifiers do not align to the KO and metadata tables")
        denominator = abundance.sum(axis=0).replace(0, np.nan)
        carrier_fraction = pd.DataFrame(index=abundance.columns)
        for process_id in rules.ProcessID:
            carrier_fraction[process_id] = (
                abundance.loc[carrier_matrix[process_id].to_numpy()].sum(axis=0) / denominator
            )
        carrier_long = (
            carrier_fraction.rename_axis("SampleID")
            .reset_index()
            .melt(id_vars="SampleID", var_name="ProcessID", value_name="RecoveredMAGCarrierFraction")
            .merge(
                rules[["ProcessID", "Element", "Process", "DirectionCaveat"]],
                on="ProcessID",
                how="left",
                validate="many_to_one",
            )
            .merge(
                metadata.rename(columns={"sample": "SampleID"}),
                on="SampleID",
                how="left",
                validate="many_to_one",
            )
        )
        write_tsv(out / "sample-carrier-fraction.tsv", carrier_long, compress=True)
        spring_carrier = (
            carrier_long.groupby(
                ["hotspring", "ProcessID", "Element", "Process"],
                observed=True,
                sort=True,
            )
            .agg(
                RecoveredMAGCarrierFraction=("RecoveredMAGCarrierFraction", "median"),
                RecruitmentRate=("metagenome2MAG_read_recruitment_rate", "median"),
                Samples=("SampleID", "nunique"),
            )
            .reset_index()
        )
        write_tsv(out / "spring-carrier-fraction.tsv", spring_carrier)
        recovery = (
            metadata.groupby("hotspring", sort=True)
            .agg(
                Samples=("sample", "nunique"),
                MedianRecruitment=("metagenome2MAG_read_recruitment_rate", "median"),
                Q25Recruitment=(
                    "metagenome2MAG_read_recruitment_rate", lambda x: x.quantile(0.25)
                ),
                Q75Recruitment=(
                    "metagenome2MAG_read_recruitment_rate", lambda x: x.quantile(0.75)
                ),
                MedianTemperature=("temperature", "median"),
                MedianpH=("pH", "median"),
            )
            .reset_index()
        )
        write_tsv(out / "mag-recovery-ceiling.tsv", recovery)

        phyla = mag_metadata.assign(Phylum=mag_metadata.taxonomy.map(phylum_from_taxonomy))[
            ["MAG", "Phylum"]
        ].set_index("MAG")
        phylum_rows: list[dict[str, object]] = []
        mean_abundance = abundance.mean(axis=1)
        for process_id in rules.ProcessID:
            carriers = carrier_matrix.index[carrier_matrix[process_id]]
            if len(carriers) == 0:
                continue
            temp = pd.DataFrame(
                {
                    "Phylum": phyla.loc[carriers, "Phylum"],
                    "MeanConditionalAbundance": mean_abundance.loc[carriers],
                }
            )
            for phylum, group in temp.groupby("Phylum", sort=True):
                phylum_rows.append(
                    {
                        "ProcessID": process_id,
                        "Phylum": phylum,
                        "CarrierMAGs": len(group),
                        "MeanRecoveredMAGAbundance": group.MeanConditionalAbundance.sum(),
                    }
                )
        carrier_phylum = pd.DataFrame(phylum_rows).merge(
            rules[["ProcessID", "Element", "Process"]],
            on="ProcessID",
            how="left",
            validate="many_to_one",
        )
        write_tsv(out / "carrier-phylum-summary.tsv", carrier_phylum)

    with ledger.step("Audit community-to-MAG concordance"):
        joined = spring_process.merge(
            spring_carrier,
            on=["hotspring", "ProcessID", "Element", "Process"],
            how="inner",
            validate="one_to_one",
        )
        concordance_rows: list[dict[str, object]] = []
        for row in rules.itertuples(index=False):
            frame = joined.loc[joined.ProcessID.eq(row.ProcessID)]
            nonzero_community = int(frame.CommunityIndex.gt(0).sum())
            nonzero_carrier = int(frame.RecoveredMAGCarrierFraction.gt(0).sum())
            if min(nonzero_community, nonzero_carrier) < MIN_NONZERO_SPRINGS:
                rho, pvalue, low, high, n = (np.nan, np.nan, np.nan, np.nan, len(frame))
            else:
                rho, pvalue, low, high, n = association(
                    frame.CommunityIndex.to_numpy(float),
                    frame.RecoveredMAGCarrierFraction.to_numpy(float),
                    rng_for(f"concordance|{row.ProcessID}"),
                )
            concordance_rows.append(
                {
                    "ProcessID": row.ProcessID,
                    "Element": row.Element,
                    "Process": row.Process,
                    "SpearmanRho": rho,
                    "PermutationP": pvalue,
                    "BootstrapLow95": low,
                    "BootstrapHigh95": high,
                    "IndependentSprings": n,
                    "NonzeroCommunitySprings": nonzero_community,
                    "NonzeroCarrierSprings": nonzero_carrier,
                    "MinimumNonzeroSprings": MIN_NONZERO_SPRINGS,
                }
            )
        concordance = pd.DataFrame(concordance_rows)
        concordance["FDR"] = bh_adjust(concordance.PermutationP.to_numpy())
        write_tsv(out / "community-mag-concordance.tsv", concordance)

    with ledger.step("Summarize nitrogen-step co-occurrence"):
        nitrogen_steps = ["nitrate_to_nitrite", "nitrite_to_no", "no_to_n2o", "n2o_to_n2"]
        complete = carrier_matrix[nitrogen_steps].all(axis=1)
        nitrogen_rows = []
        for process_id in nitrogen_steps:
            carriers = carrier_matrix[process_id]
            nitrogen_rows.append(
                {
                    "ProcessID": process_id,
                    "StrictCarrierMAGs": int(carriers.sum()),
                    "CarrierPhyla": int(
                        mag_metadata.set_index("MAG").loc[carrier_matrix.index[carriers], "taxonomy"]
                        .map(phylum_from_taxonomy)
                        .nunique()
                    ),
                    "MedianConditionalFraction": float(
                        carrier_fraction[process_id].median()
                    ),
                    "PositiveSamples": int(carrier_fraction[process_id].gt(0).sum()),
                    "CompleteChainMAGs": int(complete.sum()),
                }
            )
        nitrogen_summary = pd.DataFrame(nitrogen_rows).merge(
            rules[["ProcessID", "Process"]], on="ProcessID", how="left", validate="one_to_one"
        )
        write_tsv(out / "nitrogen-step-summary.tsv", nitrogen_summary)

    with ledger.step("Write contracts and versions"):
        write_tsv(out / "process-rules.tsv", rules)
        write_tsv(
            out / "selected-sample-metadata.tsv",
            metadata.sort_values("sample", kind="mergesort"),
        )
        versions = {
            package: package_version(package)
            for package in ("python", "pandas", "numpy", "scipy")
        }
        versions.update(
            {
                "python_runtime": os.sys.version.split()[0],
                "seed": SEED,
                "permutations": PERMUTATIONS,
                "bootstraps": BOOTSTRAPS,
                "minimum_nonzero_springs": MIN_NONZERO_SPRINGS,
                "analysis_unit": "56 hot-spring medians",
                "community_input": "published KOfam r111 KO proportions",
                "mag_input": "published 780-MAG KO count table",
            }
        )
        (out / "software-versions.json").write_text(
            json.dumps(versions, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        metrics = {
            "article": 62,
            "seed": SEED,
            "samples": int(metadata["sample"].nunique()),
            "hot_springs": int(metadata.hotspring.nunique()),
            "mags": int(mag_metadata.MAG.nunique()),
            "processes": int(len(rules)),
            "community_source_kos": community_source_kos,
            "mag_source_kos": mag_source_kos,
            "required_marker_kos": len(all_kos),
            "community_marker_columns_loaded": int(community_kos.shape[1] - 1),
            "mag_marker_columns_loaded": int(mag_kos.shape[1] - 1),
            "missing_community_marker_columns": len(missing_community),
            "missing_mag_marker_columns": len(missing_mag),
            "median_mag_recruitment_rate": float(
                metadata.metagenome2MAG_read_recruitment_rate.median()
            ),
            "minimum_mag_recruitment_rate": float(
                metadata.metagenome2MAG_read_recruitment_rate.min()
            ),
            "maximum_mag_recruitment_rate": float(
                metadata.metagenome2MAG_read_recruitment_rate.max()
            ),
            "strict_carrier_calls": int(evidence.StrictCarrier.sum()),
            "relaxed_carrier_calls": int(evidence.Relaxed80Carrier.sum()),
            "community_process_sample_pairs_passing_complete_gate": int(
                community_long.CommunityCompleteGate.sum()
            ),
            "significant_temperature_associations_fdr05": int(
                associations.Variable.eq("Temperature").mul(associations.FDR.lt(0.05)).sum()
            ),
            "significant_ph_associations_fdr05": int(
                associations.Variable.eq("pH").mul(associations.FDR.lt(0.05)).sum()
            ),
            "significant_community_mag_concordance_fdr05": int(concordance.FDR.lt(0.05).sum()),
            "complete_denitrification_chain_mags": int(complete.sum()),
            "biom": biom_audit,
            "rules_sha256": digest(rules_path),
            "diting_formula_sha256": digest(cache / "diting-pathway-formulas-v0.3.txt"),
            "completed_steps": len(ledger.rows) + 1,
        }
        (out / "analysis-metrics.json").write_text(
            json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        contract = {
            "dataset_doi": manifest["dataset_doi"],
            "paper_doi": manifest["paper_doi"],
            "figshare_version": manifest["figshare_version"],
            "kofam_release_in_source_study": "r111",
            "kofam_evalue_threshold_in_source_study": "1e-10",
            "diting_commit": "53e1d3edb84be08b7aacb79ac588be671250b477",
            "primary_inference_unit": "hot spring median",
            "seed": SEED,
            "permutations": PERMUTATIONS,
            "bootstraps": BOOTSTRAPS,
            "minimum_nonzero_springs": MIN_NONZERO_SPRINGS,
            "interpretation": "genetic potential; no activity, direction, or rate claim",
        }
        (out / "analysis-contract.json").write_text(
            json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
