#!/usr/bin/env python3
"""Prepare checksum-identified species representatives and reads for Article 48."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

from article41_44_utils import dump_json, fasta_summary, read_tsv, sha256, write_tsv


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
            raise ValueError(f"Frozen checksum mismatch: {relative}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    frozen45 = root / "data/small/45-drep-dereplication-frozen"
    frozen41 = root / "data/small/41-read-mapping-depth-frozen"
    verify_manifest(frozen45)
    verify_manifest(frozen41)

    genomes_dir = work / "inputs/genomes"
    if genomes_dir.exists():
        shutil.rmtree(genomes_dir)
    genomes_dir.mkdir(parents=True)

    membership = [
        row for row in read_tsv(frozen45 / "cluster-membership.tsv.gz")
        if row["Branch"] == "Species 95% ANI" and row["IsRepresentative"].lower() == "true"
    ]
    if len(membership) != 24:
        raise ValueError(f"Expected 24 Article 45 species representatives, observed {len(membership)}")
    genomes: list[dict[str, object]] = []
    for index, row in enumerate(sorted(membership, key=lambda item: cluster_key(item["Cluster"])), start=1):
        sgb = f"SGB_{index:03d}"
        source = frozen45 / "representative-genomes" / f"{row['Genome']}.gz"
        target = genomes_dir / f"{sgb}.fna"
        if not source.is_file():
            raise FileNotFoundError(source)
        with gzip.open(source, "rb") as input_handle, target.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle)
        observed_hash = sha256(target)
        if observed_hash != row["SHA256"]:
            raise ValueError(f"Representative checksum mismatch: {row['Genome']}")
        fasta, _ = fasta_summary(target)
        genomes.append({
            "SGB": sgb,
            "dRepCluster": row["Cluster"],
            "Representative": row["Genome"],
            "RepresentativeSHA256": observed_hash,
            "GenomeBp": fasta["TotalBp"],
            "Contigs": fasta["Contigs"],
            "N50Bp": fasta["N50Bp"],
            "Completeness": float(row["Completeness"]),
            "Contamination": float(row["Contamination"]),
            "MIMAGQuality": row["MIMAGQuality"],
            "FrozenSource": f"data/small/45-drep-dereplication-frozen/representative-genomes/{row['Genome']}.gz",
            "Path": str(target),
        })
    write_tsv(work / "genome-ledger.tsv", genomes)

    input_audit = {row["Role"]: row for row in read_tsv(frozen41 / "input-audit.tsv")}
    sample_specs = (
        ("MOCK1", "ERR9765746", "MOCK1-R1", "MOCK1-R2", 1_999_853),
        ("MOCK2", "ERR9765747", "MOCK2-R1", "MOCK2-R2", 1_999_888),
    )
    samples: list[dict[str, object]] = []
    for sample, accession, role1, role2, pairs in sample_specs:
        first, second = input_audit[role1], input_audit[role2]
        r1, r2 = root / first["Path"], root / second["Path"]
        for path, expected in ((r1, first["SHA256"]), (r2, second["SHA256"])):
            if not path.is_file() or sha256(path) != expected:
                raise ValueError(f"Read checksum mismatch: {path}")
        samples.append({
            "Sample": sample,
            "RunAccession": accession,
            "ReadPairs": pairs,
            "R1": str(r1),
            "R2": str(r2),
            "R1SHA256": first["SHA256"],
            "R2SHA256": second["SHA256"],
            "Source": "PRJEB52977 clean fixed subset from Article 30",
        })
    write_tsv(work / "samples.tsv", samples)
    write_tsv(work / "input-lineage.tsv", [
        {
            "Output": "Article 48 species catalog",
            "ImmediateInput": "24 Article 45 95%-ANI representatives",
            "Transformation": "decompress checksum-covered FASTA and assign SGB_001..SGB_024",
            "TruthUsed": "No",
            "Evidence": "genome-ledger.tsv",
        },
        {
            "Output": "Article 48 read inputs",
            "ImmediateInput": "ERR9765746 and ERR9765747 clean fixed subsets",
            "Transformation": "SHA-256 verification only",
            "TruthUsed": "No",
            "Evidence": "samples.tsv",
        },
    ])
    dump_json(work / "run-contract.json", {
        "article": 48,
        "seed": 20260748,
        "random_process": False,
        "catalog_genomes": 24,
        "catalog_source": "Article 45 species95 representatives",
        "mapper": "strobealign",
        "main_read_identity_pct": 95,
        "strict_read_identity_pct": 97,
        "minimum_read_aligned_pct": 75,
        "proper_pairs_only": True,
        "exclude_supplementary": True,
        "minimum_covered_fraction_pct": 0,
        "reporting_detection_rule": {"covered_fraction_min_pct": 50, "mean_depth_min": 1},
        "reporting_high_breadth_pct": 90,
        "contig_end_exclusion_bp": 75,
        "trim_percentiles": [5, 95],
        "methods": ["mean", "trimmed_mean", "covered_fraction", "relative_abundance", "count", "anir", "length"],
        "sample_order": ["MOCK1", "MOCK2"],
        "genome_order": [row["SGB"] for row in genomes],
        "truth_used_for_mapping_or_detection": False,
    })
    (work / ".article48-inputs-complete").write_text("complete\n", encoding="utf-8")
    print(f"Prepared {len(genomes)} representatives and {len(samples)} samples")


if __name__ == "__main__":
    main()
