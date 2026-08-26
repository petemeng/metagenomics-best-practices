#!/usr/bin/env python3
"""Prepare GTDB-audited MAGs and closest references for domain-specific trees."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

from article41_44_utils import dump_json, read_tsv, sha256, write_tsv


def cluster_key(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in value.split("_"))


def verify_manifest(root: Path) -> None:
    for line in (root / "file-checksums.sha256").read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        path = root / relative
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Frozen checksum mismatch: {root.name}/{relative}")


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1", "yes"}


def optional_float(value: str) -> float | None:
    return float(value) if value.strip() else None


def ncbi_accession(value: str) -> str:
    value = value.strip()
    for prefix in ("RS_", "GB_"):
        if value.startswith(prefix):
            value = value.removeprefix(prefix)
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    if work.exists():
        shutil.rmtree(work)
    (work / "inputs/genomes").mkdir(parents=True)
    (work / "inputs/trees").mkdir(parents=True)
    (work / "logs").mkdir()

    frozen45 = root / "data/small/45-drep-dereplication-frozen"
    frozen46 = root / "data/small/46-gtdbtk-taxonomy-frozen"
    verify_manifest(frozen45)
    verify_manifest(frozen46)
    taxonomy = {row["SGB"]: row for row in read_tsv(frozen46 / "taxonomy-summary.tsv")}
    representatives = [
        row for row in read_tsv(frozen45 / "cluster-membership.tsv.gz")
        if row["Branch"] == "Species 95% ANI" and as_bool(row["IsRepresentative"])
    ]
    representatives.sort(key=lambda row: cluster_key(row["Cluster"]))
    if len(representatives) != 24 or len(taxonomy) != 24:
        raise ValueError("Article 47 requires the complete 24-SGB Article 45/46 coordinate")

    genomes: list[dict[str, object]] = []
    references: dict[str, dict[str, object]] = {}
    novelty: list[dict[str, object]] = []
    for index, row in enumerate(representatives, start=1):
        sgb = f"SGB_{index:03d}"
        tax = taxonomy[sgb]
        domain = tax["Domain"].removeprefix("d__")
        if domain not in {"Bacteria", "Archaea"}:
            raise ValueError(f"Unsupported or unresolved domain for {sgb}: {domain}")
        source = frozen45 / "representative-genomes" / f"{row['Genome']}.gz"
        target = work / "inputs/genomes" / f"{sgb}.fna"
        with gzip.open(source, "rb") as input_handle, target.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle)
        if sha256(target) != row["SHA256"]:
            raise ValueError(f"Representative checksum mismatch: {row['Genome']}")

        reference_raw = tax["FastANIReference"] or tax["ClosestPlacementReference"]
        reference = ncbi_accession(reference_raw)
        if reference and not reference.startswith(("GCA_", "GCF_")):
            raise ValueError(f"Non-NCBI closest reference identifier for {sgb}: {reference_raw}")
        if reference:
            entry = references.setdefault(reference, {
                "ReferenceAccession": reference,
                "Domain": domain,
                "LinkedSGBs": [],
                "GTDBReferenceID": reference_raw,
            })
            if entry["Domain"] != domain:
                raise ValueError(f"Reference linked across domains: {reference}")
            entry["LinkedSGBs"].append(sgb)

        ani = optional_float(tax["FastANIANI"])
        radius = optional_float(tax["FastANIReferenceRadiusPct"])
        af = optional_float(tax["FastANIAFPct"])
        radius_cleared = bool(
            ani is not None and radius is not None and af is not None
            and ani >= radius and af >= 50
        )
        quality_eligible = float(row["Completeness"]) >= 90 and float(row["Contamination"]) < 5
        pre_tree_candidate = not as_bool(tax["SpeciesAssigned"]) and not radius_cleared and quality_eligible
        novelty.append({
            "SGB": sgb,
            "Domain": domain,
            "Representative": row["Genome"],
            "Completeness": float(row["Completeness"]),
            "Contamination": float(row["Contamination"]),
            "GTDBTaxonomy": tax["GTDBTaxonomy"],
            "SpeciesAssigned": as_bool(tax["SpeciesAssigned"]),
            "ClosestReference": reference,
            "FastANIANI": ani if ani is not None else "",
            "SpeciesRadiusPct": radius if radius is not None else "",
            "AlignmentFractionPct": af if af is not None else "",
            "SpeciesRadiusAndAFCleared": radius_cleared,
            "QualityEligibleForNoveltyReview": quality_eligible,
            "PreTreeNovelSpeciesCandidate": pre_tree_candidate,
            "NovelSpeciesClaim": False,
            "Interpretation": (
                "GTDB species assignment supported by reference-specific ANI/AF"
                if as_bool(tax["SpeciesAssigned"]) and radius_cleared
                else "Taxonomic boundary requires phylogenomic and broader reference review"
                if pre_tree_candidate
                else "Not eligible for a novel-species claim from this evidence"
            ),
        })
        genomes.append({
            "SGB": sgb,
            "Domain": domain,
            "Representative": row["Genome"],
            "RepresentativeSHA256": row["SHA256"],
            "Completeness": float(row["Completeness"]),
            "Contamination": float(row["Contamination"]),
            "Path": str(target),
            "ClosestReference": reference,
        })

    reference_rows = []
    for accession, row in sorted(references.items()):
        reference_rows.append({
            **row,
            "LinkedSGBs": ",".join(row["LinkedSGBs"]),
            "ImmutableAccessionVersion": "." in accession,
        })
    if not reference_rows or any(not row["ImmutableAccessionVersion"] for row in reference_rows):
        raise ValueError("Each selected NCBI reference must have an immutable accession version")

    for domain in ("Bacteria", "Archaea"):
        domain_key = domain.lower()
        domain_genomes = [row for row in genomes if row["Domain"] == domain]
        domain_refs = [row for row in reference_rows if row["Domain"] == domain]
        if len(domain_genomes) < 2 or len(domain_refs) < 1:
            raise ValueError(f"Insufficient query/reference sampling for {domain}")
        paths = work / "inputs/trees" / f"{domain_key}-query-fastas.txt"
        paths.write_text("".join(f"{row['Path']}\n" for row in domain_genomes), encoding="utf-8")
        accessions = work / "inputs/trees" / f"{domain_key}-reference-accessions.txt"
        accessions.write_text("".join(f"{row['ReferenceAccession']}\n" for row in domain_refs), encoding="utf-8")
        mappings = []
        for row in domain_genomes:
            mappings.append(f"{Path(row['Path']).name}\t{row['SGB']}\tQUERY_MAG\n")
        for row in domain_refs:
            mappings.append(f"{row['ReferenceAccession']}\tREF_{row['ReferenceAccession']}\tREFERENCE\n")
        (work / "inputs/trees" / f"{domain_key}-labels.tsv").write_text("".join(mappings), encoding="utf-8")

    write_tsv(work / "genome-ledger.tsv", genomes)
    write_tsv(work / "reference-request-ledger.tsv", reference_rows)
    write_tsv(work / "novelty-pretree-audit.tsv", novelty)
    write_tsv(work / "input-lineage.tsv", [{
        "Output": "Article 47 domain-specific phylogenomic inputs",
        "ImmediateInput": "Article 45 representatives plus Article 46 closest GTDB references",
        "Transformation": "checksum-exact decompression; immutable NCBI accession selection",
        "TruthUsed": "No—known mock identities do not select references or novelty status",
        "Evidence": "genome-ledger.tsv; reference-request-ledger.tsv; novelty-pretree-audit.tsv",
    }])
    dump_json(work / "run-contract.json", {
        "article": 47,
        "seed": 20260747,
        "input_sgbs": 24,
        "domains": ["Bacteria", "Archaea"],
        "reference_selection": "GTDB-Tk R232 fastANI reference, else closest placement reference",
        "reference_accessions_versioned": True,
        "gtt_hmm_sets": {"Bacteria": "Bacteria", "Archaea": "Archaea"},
        "gtotree_version": "1.8.17",
        "iqtree_version": "3.1.3",
        "iqtree_model": "MFP",
        "ultrafast_bootstrap_replicates": 1000,
        "sh_alrt_replicates": 1000,
        "truth_used_for_novelty_or_reference_selection": False,
        "novel_species_claim_requires_manual_taxonomic_review": True,
    })
    (work / ".article47-inputs-complete").write_text("complete\n", encoding="utf-8")
    print(f"Prepared 24 SGBs and {len(reference_rows)} unique closest references")


if __name__ == "__main__":
    main()
