#!/usr/bin/env python3
"""Prepare a checksum-identified redundant MAG pool for Article 45."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path

from article41_44_utils import dump_json, fasta_summary, read_tsv, sha256, write_tsv


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()

    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    genomes_dir = work / "inputs" / "genomes"
    genomes_dir.mkdir(parents=True, exist_ok=True)

    article42_work = root / "work" / "article42"
    article44_work = root / "work" / "article44"
    article42_table = article42_work / "summary" / "bin-quality-truth-audit.tsv"
    article44_table = root / "data" / "small" / "44-mag-qc-mimag-graph-frozen" / "mag-quality-summary.tsv"
    if not article42_table.is_file() or not article44_table.is_file():
        raise FileNotFoundError("Articles 42 and 44 evidence must exist before Article 45")

    records: list[dict[str, object]] = []
    for row in read_tsv(article42_table):
        if not is_true(row["QCMinimumPass"]):
            continue
        candidate = row["CandidateFASTA"].replace("${ARTICLE42_WORK_DIR}", str(article42_work))
        source = Path(candidate)
        if not source.is_file():
            raise FileNotFoundError(source)
        target = genomes_dir / f"{row['CandidateID']}.fna"
        shutil.copy2(source, target)
        summary, _ = fasta_summary(target)
        records.append(
            {
                "Genome": target.name,
                "GenomeID": row["CandidateID"],
                "SourceStage": "Article42-QC-pass",
                "SourceBranch": row["Branch"],
                "Completeness": float(row["CheckM2Completeness"]),
                "Contamination": float(row["CheckM2Contamination"]),
                "MIMAGQuality": "Not re-audited for complete markers",
                "Contigs": summary["Contigs"],
                "GenomeBp": summary["TotalBp"],
                "N50Bp": summary["N50Bp"],
                "SHA256": sha256(target),
                "Path": str(target),
            }
        )

    for row in read_tsv(article44_table):
        source = article44_work / "bins" / f"{row['MAG']}.fna"
        if not source.is_file():
            raise FileNotFoundError(source)
        target = genomes_dir / source.name
        shutil.copy2(source, target)
        summary, _ = fasta_summary(target)
        if sha256(target) != row["MAGFASTA_SHA256"]:
            raise ValueError(f"Article 44 checksum mismatch: {source.name}")
        records.append(
            {
                "Genome": target.name,
                "GenomeID": row["MAG"],
                "SourceStage": "Article44-selected",
                "SourceBranch": "Binette selected set",
                "Completeness": float(row["CheckM2Completeness"]),
                "Contamination": float(row["CheckM2Contamination"]),
                "MIMAGQuality": row["MIMAGQuality"],
                "Contigs": summary["Contigs"],
                "GenomeBp": summary["TotalBp"],
                "N50Bp": summary["N50Bp"],
                "SHA256": sha256(target),
                "Path": str(target),
            }
        )

    records.sort(key=lambda row: str(row["Genome"]))
    if len(records) != 124 or len({row["Genome"] for row in records}) != 124:
        raise ValueError(f"Expected 124 unique genomes, observed {len(records)}")
    if any(float(row["Completeness"]) < 50 or float(row["Contamination"]) >= 10 for row in records):
        raise ValueError("Input pool violates the pre-registered 50%/<10% gate")

    write_tsv(work / "input-genomes.tsv", records)
    with (work / "genome-info.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome", "completeness", "contamination"])
        writer.writeheader()
        for row in records:
            writer.writerow(
                {
                    "genome": row["Genome"],
                    "completeness": row["Completeness"],
                    "contamination": row["Contamination"],
                }
            )
    (work / "genomes.txt").write_text(
        "\n".join(str(row["Path"]) for row in records) + "\n", encoding="utf-8"
    )
    write_tsv(
        work / "input-lineage.tsv",
        [
            {
                "Output": "Article 45 redundant MAG pool",
                "ImmediateInput": "101 Article 42 CheckM2+GUNC-pass bins + 23 Article 44 selected MAGs",
                "Transformation": "checksum-exact copy with collision-safe genome names",
                "TruthUsed": "No",
                "Evidence": "input-genomes.tsv",
            },
            {
                "Output": "dRep quality input",
                "ImmediateInput": "release-matched Article 42/44 CheckM2 estimates",
                "Transformation": "three-column genomeInfo.csv; no CheckM rerun inside dRep",
                "TruthUsed": "No",
                "Evidence": "genome-info.csv",
            },
        ],
    )
    dump_json(
        work / "run-contract.json",
        {
            "article": 45,
            "input_genomes": 124,
            "quality_gate": {"completeness_min_pct": 50, "contamination_max_exclusive_pct": 10},
            "main_secondary_ani": 0.95,
            "sensitivity_secondary_ani": 0.999,
            "minimum_alignment_fraction": 0.30,
            "secondary_algorithm": "fastANI",
            "input_order": "lexicographic genome basename",
            "truth_used_for_clustering_or_selection": False,
            "random_process": False,
        },
    )
    print(f"Prepared {len(records)} checksum-identified genomes in {genomes_dir}")


if __name__ == "__main__":
    main()
