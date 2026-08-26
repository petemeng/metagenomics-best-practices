#!/usr/bin/env python3
"""Summarize PanPhlAn presence/absence, annotations, distances and threshold sensitivity."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist, squareform

from article41_44_utils import dump_json, read_tsv, sha256, write_tsv


ANNOTATION_FIELDS = ("NR90", "NR50", "GO", "KO", "KEGG", "Pfam", "EC", "eggNOG")
CATEGORY_ORDER = ("Core >=95%", "Accessory 5-<95%", "Rare >0-<5%", "Undetected")


def category(count: int, total: int) -> str:
    prevalence = count / total
    if prevalence >= 0.95:
        return "Core >=95%"
    if prevalence >= 0.05:
        return "Accessory 5-<95%"
    if count > 0:
        return "Rare >0-<5%"
    return "Undetected"


def parse_annotation(path: Path) -> tuple[dict[str, dict[str, str]], list[dict[str, object]]]:
    values: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {field: set() for field in ANNOTATION_FIELDS[1:]}
    )
    anomalies: list[dict[str, object]] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = raw.split("\t")
        if tuple(fields[:8]) == ANNOTATION_FIELDS:
            continue
        if not fields or not fields[0].startswith("UniRef90_"):
            anomalies.append({
                "Line": line_number, "ObservedFields": len(fields),
                "Action": "excluded non-UniRef record", "Prefix": raw[:100],
            })
            continue
        action = "none"
        if len(fields) != 8:
            action = "retained first 8 fields; removed embedded repeated header"
        fields = (fields + [""] * 8)[:8]
        if action != "none" and fields[7].endswith("NR90"):
            fields[7] = fields[7][:-4]
            action += "; stripped appended NR90 token from eggNOG"
        if action != "none":
            anomalies.append({
                "Line": line_number, "ObservedFields": len(raw.split("\t")),
                "Action": action, "Prefix": raw[:100],
            })
        gene = fields[0]
        for field, cell in zip(ANNOTATION_FIELDS[1:], fields[1:]):
            for value in cell.split(","):
                value = value.strip()
                if value:
                    values[gene][field].add(value)
    collapsed = {
        gene: {field: ";".join(sorted(field_values)) for field, field_values in annotations.items()}
        for gene, annotations in values.items()
    }
    return collapsed, anomalies


def parse_plateau_log(path: Path) -> dict[str, dict[str, object]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    rows: dict[str, dict[str, object]] = {}
    pattern = re.compile(
        r"\[I\]\s+(\S+) median coverage: ([0-9.]+); left-side cov: ([0-9.]+); "
        r"right-side cov: ([0-9.]+); out-plateau cov: ([0-9.]+)"
    )
    for sample, median, left, right, out in pattern.findall(text):
        rows[sample] = {
            "MedianCoverage": float(median), "LeftCoverage": float(left),
            "RightCoverage": float(right), "OutPlateauCoverage": float(out),
            "MultiStrainWarning": f"{sample} WARNING: sample may contain multiple strains" in text,
        }
    return rows


def orient_axis(values: np.ndarray) -> np.ndarray:
    index = int(np.argmax(np.abs(values)))
    return values if values[index] >= 0 else -values


def pcoa(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = distance.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ (distance ** 2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    positive = eigenvalues > max(eigenvalues[0] * 1e-10, 1e-12)
    if positive.sum() < 2:
        raise ValueError("Jaccard distance has fewer than two positive PCoA axes")
    eigenvalues, eigenvectors = eigenvalues[positive], eigenvectors[:, positive]
    coordinates = eigenvectors[:, :2] * np.sqrt(eigenvalues[:2])
    coordinates[:, 0] = orient_axis(coordinates[:, 0])
    coordinates[:, 1] = orient_axis(coordinates[:, 1])
    explained = eigenvalues[:2] / eigenvalues.sum() * 100
    return coordinates, explained


def deterministic_gzip_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as input_handle, target.open("wb") as raw_output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, compresslevel=9, mtime=0) as output_handle:
            shutil.copyfileobj(input_handle, output_handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article52-run-complete").is_file():
        raise FileNotFoundError("Run run_article52_panphlan.py first")
    summary_dir = work / "summary"
    if summary_dir.exists():
        shutil.rmtree(summary_dir)
    summary_dir.mkdir(parents=True)

    metadata = pd.read_csv(work / "sample-metadata.tsv", sep="\t", dtype=str)
    metadata = metadata[["Sample", "Study", "Country"]].drop_duplicates().set_index("Sample")
    primary = pd.read_csv(work / "output/primary-matrix.tsv", sep="\t", index_col=0)
    sensitive = pd.read_csv(work / "output/sensitive-matrix.tsv", sep="\t", index_col=0)
    if primary.index.duplicated().any() or sensitive.index.duplicated().any():
        raise ValueError("Duplicate gene-family identifiers")
    if not primary.index.equals(sensitive.index):
        raise ValueError("Primary and sensitivity gene coordinates differ")
    if not set(np.unique(primary.to_numpy())).issubset({0, 1}):
        raise ValueError("Primary matrix is not binary")
    if not set(np.unique(sensitive.to_numpy())).issubset({0, 1}):
        raise ValueError("Sensitivity matrix is not binary")

    primary_samples = [column for column in primary if not column.startswith("REF_")]
    sensitive_samples = [column for column in sensitive if not column.startswith("REF_")]
    references = [column for column in primary if column.startswith("REF_")]
    if (len(primary_samples), len(sensitive_samples), len(references)) != (22, 25, 15):
        raise ValueError("Unexpected retained sample/reference counts")
    if set(sensitive_samples) != set(metadata.index):
        raise ValueError("Sensitivity sample set does not match the 25-sample metadata")
    if not np.array_equal(primary[primary_samples].to_numpy(), sensitive[primary_samples].to_numpy()):
        raise ValueError("Gene calls changed for shared samples despite identical gene thresholds")

    annotation_path = work / "pangenome/Eubacterium_rectale/panphlan_Eubacterium_rectale_annot.tsv"
    annotations, annotation_anomalies = parse_annotation(annotation_path)
    write_tsv(summary_dir / "annotation-format-audit.tsv", annotation_anomalies)
    missing_annotation_rows = sorted(set(primary.index) - set(annotations))
    if len(annotation_anomalies) != 5 or len(missing_annotation_rows) != 8:
        raise ValueError(
            f"Unexpected annotation snapshot structure: {len(annotation_anomalies)} anomalies, "
            f"{len(missing_annotation_rows)} missing gene families"
        )

    primary_counts = primary[primary_samples].sum(axis=1).astype(int)
    sensitive_counts = sensitive[sensitive_samples].sum(axis=1).astype(int)
    reference_counts = primary[references].sum(axis=1).astype(int)
    prevalence_rows: list[dict[str, object]] = []
    for gene in primary.index:
        annot = annotations.get(gene, {field: "" for field in ANNOTATION_FIELDS[1:]})
        row = {
            "GeneFamily": gene,
            "PrimaryCount": int(primary_counts[gene]),
            "PrimaryN": len(primary_samples),
            "PrimaryPrevalence": round(primary_counts[gene] / len(primary_samples), 8),
            "PrimaryCategory": category(int(primary_counts[gene]), len(primary_samples)),
            "SensitivityCount": int(sensitive_counts[gene]),
            "SensitivityN": len(sensitive_samples),
            "SensitivityPrevalence": round(sensitive_counts[gene] / len(sensitive_samples), 8),
            "SensitivityCategory": category(int(sensitive_counts[gene]), len(sensitive_samples)),
            "ReferenceCount": int(reference_counts[gene]),
            "ReferenceN": len(references),
            "ReferencePrevalence": round(reference_counts[gene] / len(references), 8),
            **annot,
        }
        prevalence_rows.append(row)
    write_tsv(summary_dir / "gene-family-prevalence.tsv.gz", prevalence_rows)
    write_tsv(summary_dir / "gene-family-prevalence-head.tsv", prevalence_rows[:8])

    branch_rows: list[dict[str, object]] = []
    for branch, frame, samples in (
        ("Primary", primary, primary_samples), ("Sensitivity", sensitive, sensitive_samples)
    ):
        counts = frame[samples].sum(axis=1).astype(int)
        sample_gene_counts = frame[samples].sum(axis=0).astype(int)
        categories = [category(int(value), len(samples)) for value in counts]
        branch_rows.append({
            "Branch": branch, "InputSamples": 25, "RetainedSamples": len(samples),
            "ReferenceGenomes": len(references), "PangenomeGeneFamilies": len(frame),
            "MedianGeneFamiliesPerSample": float(sample_gene_counts.median()),
            "MinimumGeneFamiliesPerSample": int(sample_gene_counts.min()),
            "MaximumGeneFamiliesPerSample": int(sample_gene_counts.max()),
            "StrictCore100Pct": int((counts == len(samples)).sum()),
            "OperationalCore95Pct": categories.count("Core >=95%"),
            "Accessory5To95Pct": categories.count("Accessory 5-<95%"),
            "RareBelow5Pct": categories.count("Rare >0-<5%"),
            "Undetected": categories.count("Undetected"),
        })
    reference_gene_counts = primary[references].sum(axis=0).astype(int)
    branch_rows.append({
        "Branch": "Reference genomes", "InputSamples": 15, "RetainedSamples": 15,
        "ReferenceGenomes": 15, "PangenomeGeneFamilies": len(primary),
        "MedianGeneFamiliesPerSample": float(reference_gene_counts.median()),
        "MinimumGeneFamiliesPerSample": int(reference_gene_counts.min()),
        "MaximumGeneFamiliesPerSample": int(reference_gene_counts.max()),
        "StrictCore100Pct": int((reference_counts == len(references)).sum()),
        "OperationalCore95Pct": int((reference_counts / len(references) >= 0.95).sum()),
        "Accessory5To95Pct": int(((reference_counts / len(references) >= 0.05) & (reference_counts / len(references) < 0.95)).sum()),
        "RareBelow5Pct": 0, "Undetected": int((reference_counts == 0).sum()),
    })
    write_tsv(summary_dir / "pangenome-summary.tsv", branch_rows)

    primary_log = parse_plateau_log(work / "logs/panphlan-primary.stdout.log")
    sensitive_log = parse_plateau_log(work / "logs/panphlan-sensitivity.stdout.log")
    if set(primary_log) != set(metadata.index) or set(sensitive_log) != set(metadata.index):
        raise ValueError("Plateau audit log is missing samples")
    sample_rows = []
    for sample, meta in metadata.iterrows():
        metrics = primary_log[sample]
        primary_retained = sample in primary_samples
        sensitivity_retained = sample in sensitive_samples
        reason = "Passed all primary plateau gates" if primary_retained else "Failed primary left-side plateau gate"
        sample_rows.append({
            "Sample": sample, "Study": meta["Study"], "Country": meta["Country"],
            "MedianCoverage": metrics["MedianCoverage"],
            "LeftCoverage": metrics["LeftCoverage"],
            "RightCoverage": metrics["RightCoverage"],
            "OutPlateauCoverage": metrics["OutPlateauCoverage"],
            "PrimaryRetained": primary_retained,
            "SensitivityRetained": sensitivity_retained,
            "PrimaryGeneFamilies": int(primary[sample].sum()) if primary_retained else "NA",
            "SensitivityGeneFamilies": int(sensitive[sample].sum()),
            "MultiStrainWarning": metrics["MultiStrainWarning"],
            "PrimaryDecision": reason,
        })
    write_tsv(summary_dir / "sample-filter-audit.tsv", sample_rows)

    primary_array = primary[primary_samples].T.to_numpy(dtype=np.uint8)
    distance = squareform(pdist(primary_array, metric="jaccard"))
    pair_rows: list[dict[str, object]] = []
    for left in range(len(primary_samples)):
        for right in range(left + 1, len(primary_samples)):
            sample1, sample2 = primary_samples[left], primary_samples[right]
            meta1, meta2 = metadata.loc[sample1], metadata.loc[sample2]
            if meta1["Study"] == meta2["Study"]:
                stratum = "Same study"
            elif meta1["Country"] == meta2["Country"]:
                stratum = "Same country, different study"
            else:
                stratum = "Different country"
            pair_rows.append({
                "Sample1": sample1, "Sample2": sample2,
                "Study1": meta1["Study"], "Study2": meta2["Study"],
                "Country1": meta1["Country"], "Country2": meta2["Country"],
                "PairStratum": stratum, "JaccardDistance": round(float(distance[left, right]), 8),
            })
    write_tsv(summary_dir / "pairwise-jaccard.tsv", pair_rows)
    pair_frame = pd.DataFrame(pair_rows)
    pair_summary = []
    for stratum in ("Same study", "Same country, different study", "Different country"):
        values = pair_frame.loc[pair_frame.PairStratum == stratum, "JaccardDistance"].to_numpy(float)
        pair_summary.append({
            "PairStratum": stratum, "Pairs": len(values),
            "MedianJaccard": round(float(np.median(values)), 8),
            "Q1Jaccard": round(float(np.quantile(values, 0.25)), 8),
            "Q3Jaccard": round(float(np.quantile(values, 0.75)), 8),
            "MeanJaccard": round(float(np.mean(values)), 8),
        })
    write_tsv(summary_dir / "pairwise-jaccard-summary.tsv", pair_summary)

    distance_rows = []
    for index, sample in enumerate(primary_samples):
        distance_rows.append({
            "Sample": sample,
            **{other: round(float(distance[index, j]), 8) for j, other in enumerate(primary_samples)},
        })
    write_tsv(summary_dir / "sample-jaccard-matrix.tsv", distance_rows)
    coordinates, explained = pcoa(distance)
    pcoa_rows = []
    for index, sample in enumerate(primary_samples):
        pcoa_rows.append({
            "Sample": sample, "Study": metadata.loc[sample, "Study"],
            "Country": metadata.loc[sample, "Country"],
            "PCoA1": round(float(coordinates[index, 0]), 8),
            "PCoA2": round(float(coordinates[index, 1]), 8),
            "PCoA1Pct": round(float(explained[0]), 4),
            "PCoA2Pct": round(float(explained[1]), 4),
        })
    write_tsv(summary_dir / "pcoa-jaccard.tsv", pcoa_rows)

    primary_categories = pd.Series(
        [category(int(value), len(primary_samples)) for value in primary_counts], index=primary.index
    )
    sensitive_categories = pd.Series(
        [category(int(value), len(sensitive_samples)) for value in sensitive_counts], index=sensitive.index
    )
    transition_rows = []
    for source in CATEGORY_ORDER:
        for target in CATEGORY_ORDER:
            transition_rows.append({
                "PrimaryCategory": source, "SensitivityCategory": target,
                "GeneFamilies": int(((primary_categories == source) & (sensitive_categories == target)).sum()),
            })
    write_tsv(summary_dir / "category-transition.tsv", transition_rows)

    annotation_rows = []
    for group in CATEGORY_ORDER:
        genes = primary.index[primary_categories == group]
        for database in ANNOTATION_FIELDS[1:]:
            annotated = sum(bool(annotations.get(gene, {}).get(database, "")) for gene in genes)
            annotation_rows.append({
                "Category": group, "Database": database, "GeneFamilies": len(genes),
                "AnnotatedGeneFamilies": annotated,
                "AnnotatedPct": round(100 * annotated / len(genes), 4) if len(genes) else math.nan,
            })
    write_tsv(summary_dir / "annotation-coverage.tsv", annotation_rows)

    accessory_genes = [gene for gene in primary.index if primary_categories[gene] == "Accessory 5-<95%"]
    ranked = []
    for gene in accessory_genes:
        annot = annotations.get(gene, {})
        nonempty = sum(bool(annot.get(field, "")) for field in ANNOTATION_FIELDS[1:])
        if nonempty:
            ranked.append((abs(primary_counts[gene] / len(primary_samples) - 0.5), -nonempty, gene))
    selected = [gene for _, _, gene in sorted(ranked)[:20]]
    selected_array = primary.loc[selected, primary_samples].to_numpy(dtype=np.uint8)
    sample_order = leaves_list(linkage(squareform(distance, checks=False), method="average"))
    feature_order = leaves_list(linkage(pdist(selected_array, metric="jaccard"), method="average"))
    sample_rank = {primary_samples[index]: rank + 1 for rank, index in enumerate(sample_order)}
    feature_rank = {selected[index]: rank + 1 for rank, index in enumerate(feature_order)}
    heatmap_rows = []
    for gene in selected:
        annot = annotations[gene]
        label_value = next(
            (annot[field].split(";")[0] for field in ("KO", "Pfam", "GO", "eggNOG") if annot.get(field)),
            "annotated",
        )
        label = f"{gene.removeprefix('UniRef90_')} | {label_value}"
        for sample in primary_samples:
            heatmap_rows.append({
                "GeneFamily": gene, "FeatureLabel": label,
                "FeatureOrder": feature_rank[gene], "Sample": sample,
                "SampleOrder": sample_rank[sample], "Country": metadata.loc[sample, "Country"],
                "Present": int(primary.loc[gene, sample]),
                "PrevalencePct": round(100 * primary_counts[gene] / len(primary_samples), 2),
            })
    write_tsv(summary_dir / "accessory-feature-heatmap.tsv", heatmap_rows)

    for source_name, target_name in (
        ("primary-matrix.tsv", "primary-presence-absence.tsv.gz"),
        ("sensitive-matrix.tsv", "sensitive-presence-absence.tsv.gz"),
        ("coverage-matrix.tsv", "coverage-matrix.tsv.gz"),
    ):
        deterministic_gzip_copy(work / "output" / source_name, summary_dir / target_name)
    native_rows = []
    for path in sorted(summary_dir.glob("*.tsv.gz")):
        native_rows.append({"File": path.name, "Bytes": path.stat().st_size, "SHA256": sha256(path)})
    write_tsv(summary_dir / "native-output-manifest.tsv", native_rows)

    nearest = min(pair_rows, key=lambda row: float(row["JaccardDistance"]))
    result = {
        "article": 52,
        "species_label_in_official_tutorial": "Eubacterium rectale",
        "official_samples": 25,
        "studies": int(metadata.Study.nunique()),
        "countries": int(metadata.Country.nunique()),
        "reference_genomes": len(references),
        "pangenome_gene_families": len(primary),
        "coverage_matrix_gene_families": int(pd.read_csv(work / "output/coverage-matrix.tsv", sep="\t", usecols=[0]).shape[0]),
        "primary_retained_samples": len(primary_samples),
        "sensitivity_retained_samples": len(sensitive_samples),
        "primary_excluded_samples": sorted(set(sensitive_samples) - set(primary_samples)),
        "primary_multistrain_warnings": sum(bool(row["MultiStrainWarning"]) for row in sample_rows if row["PrimaryRetained"]),
        "primary_operational_core": int((primary_categories == "Core >=95%").sum()),
        "primary_accessory": int((primary_categories == "Accessory 5-<95%").sum()),
        "primary_rare": int((primary_categories == "Rare >0-<5%").sum()),
        "primary_undetected": int((primary_categories == "Undetected").sum()),
        "primary_median_genes_per_sample": float(primary[primary_samples].sum(axis=0).median()),
        "sensitivity_operational_core": int((sensitive_categories == "Core >=95%").sum()),
        "reference_strict_core": int((reference_counts == len(references)).sum()),
        "common_sample_calls_byte_equivalent": True,
        "annotation_unique_gene_families": len(annotations),
        "annotation_missing_gene_families": len(missing_annotation_rows),
        "annotation_format_anomalies": len(annotation_anomalies),
        "nearest_pair": [nearest["Sample1"], nearest["Sample2"]],
        "nearest_pair_jaccard": nearest["JaccardDistance"],
        "pcoa_axis1_pct": pcoa_rows[0]["PCoA1Pct"],
        "pcoa_axis2_pct": pcoa_rows[0]["PCoA2Pct"],
    }
    for row in pair_summary:
        key = row["PairStratum"].lower().replace(",", "").replace(" ", "_")
        result[f"{key}_pairs"] = row["Pairs"]
        result[f"{key}_median_jaccard"] = row["MedianJaccard"]
    dump_json(summary_dir / "summary.json", result)
    (work / ".article52-summary-complete").write_text("complete\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
