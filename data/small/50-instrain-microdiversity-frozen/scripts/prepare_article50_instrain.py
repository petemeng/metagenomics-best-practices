#!/usr/bin/env python3
"""Prepare checksum-tracked catalog/BAM inputs for Article 50 inStrain analysis."""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

from article41_44_utils import dump_json, fasta_records, read_tsv, sha256, write_tsv


def cluster_key(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in value.split("_"))


def verify_manifest(root: Path) -> None:
    manifest = root / "file-checksums.sha256"
    if not manifest.is_file():
        raise FileNotFoundError(manifest)
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        path = root / relative
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Frozen checksum mismatch: {root.name}/{relative}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    if work.exists():
        shutil.rmtree(work)
    inputs = work / "inputs"
    inputs.mkdir(parents=True)
    (work / "logs").mkdir()

    frozen45 = root / "data/small/45-drep-dereplication-frozen"
    frozen46 = root / "data/small/46-gtdbtk-taxonomy-frozen"
    frozen48 = root / "data/small/48-mag-abundance-coverm-frozen"
    for frozen in (frozen45, frozen46, frozen48):
        verify_manifest(frozen)

    representatives = [
        row
        for row in read_tsv(frozen45 / "cluster-membership.tsv.gz")
        if row["Branch"] == "Species 95% ANI"
        and row["IsRepresentative"].lower() == "true"
    ]
    representatives.sort(key=lambda row: cluster_key(row["Cluster"]))
    if len(representatives) != 24:
        raise ValueError(f"Expected 24 SGB representatives, observed {len(representatives)}")
    taxonomy = {
        row["SGB"]: row for row in read_tsv(frozen46 / "taxonomy-summary.tsv")
    }
    coverage = {
        (row["SGB"], row["Sample"]): row
        for row in read_tsv(frozen48 / "coverm-long.tsv.gz")
        if row["Branch"] == "Primary 95% identity"
    }

    catalog_path = inputs / "article50-sgb-catalog.fna"
    stb_path = inputs / "article50-scaffold-to-bin.tsv"
    genome_rows: list[dict[str, object]] = []
    stb_rows: list[tuple[str, str]] = []
    seen_scaffolds: set[str] = set()
    with catalog_path.open("w", encoding="ascii", newline="\n") as catalog:
        for index, row in enumerate(representatives, start=1):
            sgb = f"SGB_{index:03d}"
            source = frozen45 / "representative-genomes" / f"{row['Genome']}.gz"
            if not source.is_file():
                raise FileNotFoundError(source)
            total_bp = 0
            contigs = 0
            for original, sequence in fasta_records(source):
                scaffold = f"{sgb}~{original}"
                if scaffold in seen_scaffolds:
                    raise ValueError(f"Duplicate scaffold identifier: {scaffold}")
                seen_scaffolds.add(scaffold)
                stb_rows.append((scaffold, sgb))
                total_bp += len(sequence)
                contigs += 1
                catalog.write(f">{scaffold}\n")
                for start in range(0, len(sequence), 80):
                    catalog.write(sequence[start : start + 80] + "\n")
            if total_bp != int(row["GenomeBp"]):
                raise ValueError(f"Genome length mismatch: {sgb}")
            tax = taxonomy[sgb]
            genome_rows.append(
                {
                    "SGB": sgb,
                    "dRepCluster": row["Cluster"],
                    "Representative": row["Genome"],
                    "RepresentativeSequenceSHA256": row["SHA256"],
                    "GenomeBp": total_bp,
                    "Contigs": contigs,
                    "GTDBRelease": "R232",
                    "GTDBTaxonomy": tax["GTDBTaxonomy"],
                    "MOCK1CoverMBreadthPct": float(
                        coverage[(sgb, "MOCK1")]["CoveredFractionPct"]
                    ),
                    "MOCK2CoverMBreadthPct": float(
                        coverage[(sgb, "MOCK2")]["CoveredFractionPct"]
                    ),
                }
            )
    with stb_path.open("w", encoding="utf-8", newline="\n") as handle:
        for scaffold, sgb in stb_rows:
            handle.write(f"{scaffold}\t{sgb}\n")
    write_tsv(work / "genome-ledger.tsv", genome_rows)

    bam_specs = (
        (
            "MOCK1",
            "ERR9765746",
            root
            / "work/article48/bam/identity95/coverm-genome.ERR9765746_clean_R1.fastq.gz.bam",
        ),
        (
            "MOCK2",
            "ERR9765747",
            root
            / "work/article48/bam/identity95/coverm-genome.ERR9765747_clean_R1.fastq.gz.bam",
        ),
    )
    sample_rows: list[dict[str, object]] = []
    for sample, accession, source in bam_specs:
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(source)
        target = inputs / f"{sample}.identity95.unfiltered.bam"
        target.symlink_to(os.path.relpath(source, target.parent))
        sample_rows.append(
            {
                "Sample": sample,
                "RunAccession": accession,
                "BAM": str(target),
                "BAMSource": str(source),
                "BAMSHA256": sha256(source),
                "BAMBytes": source.stat().st_size,
                "Mapper": "Strobealign 0.17.0 via CoverM 0.8.0",
                "BAMState": "coordinate-sorted unfiltered cache; inStrain applies pair filters",
            }
        )
    write_tsv(work / "sample-ledger.tsv", sample_rows)
    write_tsv(
        work / "input-lineage.tsv",
        [
            {
                "Output": "Article 50 concatenated reference and scaffold-to-bin map",
                "ImmediateInput": "24 checksum-covered Article 45 species representatives",
                "Transformation": "prefix each contig with stable SGB ID; concatenate without sequence changes",
                "TruthUsed": "No",
                "Evidence": "genome-ledger.tsv; article50-scaffold-to-bin.tsv",
            },
            {
                "Output": "Article 50 inStrain BAM inputs",
                "ImmediateInput": "Article 48 CoverM/Strobealign unfiltered coordinate-sorted BAM cache",
                "Transformation": "symlink and SHA-256 audit; no alignment editing",
                "TruthUsed": "No",
                "Evidence": "sample-ledger.tsv",
            },
        ],
    )
    dump_json(
        work / "run-contract.json",
        {
            "article": 50,
            "seed": 20260750,
            "random_process": False,
            "catalog_sgbs": 24,
            "samples": ["MOCK1", "MOCK2"],
            "run_accessions": ["ERR9765746", "ERR9765747"],
            "min_read_ani": 0.95,
            "min_mapq": 2,
            "pairing_filter": "paired_only",
            "min_cov": 5,
            "min_freq": 0.05,
            "fdr": 1e-6,
            "rarefied_coverage": 1000000000,
            "rarefied_metric_enabled": False,
            "rarefaction_audit": "disabled because inStrain 1.10.0 uses unseeded NumPy sampling and exposes no CLI seed",
            "compare_database_mode": True,
            "compare_presence_breadth": 0.5,
            "same_strain_reporting_rule": {
                "popANI_min_pct": 99.999,
                "percent_genome_compared_min_pct": 50,
            },
            "truth_used_for_profiling_or_comparison": False,
            "cached_bam_filtering_note": "Article 48 BAM is unfiltered; inStrain applies its own explicit read-pair filters",
        },
    )
    (work / ".article50-inputs-complete").write_text("complete\n", encoding="utf-8")
    print(
        f"Prepared {len(genome_rows)} SGBs, {len(stb_rows)} scaffolds, "
        f"and {len(sample_rows)} BAM inputs"
    )


if __name__ == "__main__":
    main()
