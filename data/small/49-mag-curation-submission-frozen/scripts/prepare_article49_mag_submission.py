#!/usr/bin/env python3
"""Build a fail-closed MAG curation and NCBI pre-submission review package."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import math
import shutil
import statistics
from collections import defaultdict
from pathlib import Path

from article41_44_utils import (
    dump_json,
    fasta_records,
    read_tsv,
    sha256,
    write_tsv,
)


RUNS = ("ERR9765746", "ERR9765747")
MISSING_SOURCE = "NOT_AVAILABLE_FROM_TUTORIAL_SOURCE—INVESTIGATOR_MUST_SUPPLY"


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


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1", "yes"}


def robust_z(values: list[float]) -> list[float]:
    if not values:
        return []
    center = statistics.median(values)
    mad = statistics.median(abs(value - center) for value in values)
    if mad == 0:
        return [0.0 for _ in values]
    return [(value - center) / (1.4826 * mad) for value in values]


def sequence_digest(records: list[tuple[str, str]]) -> str:
    digest = hashlib.sha256()
    for _, sequence in records:
        digest.update(sequence.encode("ascii"))
    return digest.hexdigest()


def gunzip_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_review_fasta(source: Path, target: Path, isolate: str) -> dict[str, object]:
    records = list(fasta_records(source))
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as zipped:
            for index, (_, sequence) in enumerate(records, start=1):
                header = f">{isolate}_contig{index:05d} [SRA={','.join(RUNS)}]\n"
                zipped.write(header.encode("ascii"))
                for start in range(0, len(sequence), 80):
                    zipped.write(sequence[start:start + 80].encode("ascii") + b"\n")
    return {
        "Contigs": len(records),
        "TotalBp": sum(len(sequence) for _, sequence in records),
        "SequenceOnlySHA256": sequence_digest(records),
        "ReviewFASTA_SHA256": sha256(target),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    (work / "logs").mkdir()

    frozen = {
        number: root / f"data/small/{slug}"
        for number, slug in {
            41: "41-read-mapping-depth-frozen",
            43: "43-bin-refinement-frozen",
            44: "44-mag-qc-mimag-graph-frozen",
            45: "45-drep-dereplication-frozen",
            46: "46-gtdbtk-taxonomy-frozen",
            48: "48-mag-abundance-coverm-frozen",
        }.items()
    }
    for directory in frozen.values():
        verify_manifest(directory)

    membership45 = [
        row for row in read_tsv(frozen[45] / "cluster-membership.tsv.gz")
        if row["Branch"] == "Species 95% ANI" and as_bool(row["IsRepresentative"])
    ]
    if len(membership45) != 24:
        raise ValueError("Expected 24 Article 45 species representatives")
    membership45.sort(key=lambda row: cluster_key(row["Cluster"]))
    sgb_rows = []
    for index, row in enumerate(membership45, start=1):
        sgb_rows.append({"SGB": f"SGB_{index:03d}", **row})

    quality44 = {row["MAG"]: row for row in read_tsv(frozen[44] / "mag-quality-summary.tsv")}
    graph44 = {row["MAG"]: row for row in read_tsv(frozen[44] / "assembly-graph-audit.tsv")}
    taxonomy46 = {row["SGB"]: row for row in read_tsv(frozen[46] / "taxonomy-summary.tsv")}
    if set(taxonomy46) != {row["SGB"] for row in sgb_rows}:
        raise ValueError("Article 46 taxonomy coordinate does not match Article 45 catalog")

    coverage48 = [
        row for row in read_tsv(frozen[48] / "coverm-long.tsv.gz")
        if row["Branch"] == "Primary 95% identity"
    ]
    coverage = {(row["SGB"], row["Sample"]): row for row in coverage48}
    if len(coverage) != 48:
        raise ValueError("Expected two primary CoverM rows for each of 24 SGBs")

    depth41 = read_tsv(frozen[41] / "contig-depth-long.tsv.gz")
    depth = {(row["Contig"], row["Sample"]): row for row in depth41}
    membership43 = read_tsv(frozen[43] / "selected-refinement-membership.tsv.gz")
    contigs_by_mag: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in membership43:
        contigs_by_mag[row["RefinedID"]].append(row)

    dispositions: list[dict[str, object]] = []
    review_candidates: list[dict[str, object]] = []
    for row in sgb_rows:
        mag = row["Genome"].removesuffix(".fna")
        article44_audited = row["SourceStage"] == "Article44-selected" and mag in quality44
        numeric_gate = float(row["Completeness"]) >= 90 and int(row["GenomeBp"]) >= 100_000
        if article44_audited and numeric_gate:
            disposition = "TECHNICAL_REVIEW_SET"
            reason = "Article 44 MIMAG/GUNC record present; NCBI numeric genome gates pass"
        elif article44_audited:
            disposition = "HOLD_NUMERIC_GENOME_GATE"
            reason = "CheckM2 completeness is below the current NCBI 90% minimum"
        else:
            disposition = "HOLD_FULL_MAG_AUDIT_MISSING"
            reason = "Representative was not re-audited through the complete Article 44 MIMAG/GUNC/manual-review contract"
        item = {
            "SGB": row["SGB"],
            "Representative": row["Genome"],
            "SourceStage": row["SourceStage"],
            "Completeness": float(row["Completeness"]),
            "Contamination": float(row["Contamination"]),
            "GenomeBp": int(row["GenomeBp"]),
            "Article44CompleteAudit": article44_audited,
            "NCBINumericGenomeGate": numeric_gate,
            "Disposition": disposition,
            "Reason": reason,
        }
        dispositions.append(item)
        if disposition == "TECHNICAL_REVIEW_SET":
            review_candidates.append({**row, "MAG": mag})
    if len(review_candidates) != 12:
        raise ValueError(f"Expected 12 technical-review candidates, observed {len(review_candidates)}")

    anomaly_rows: list[dict[str, object]] = []
    curation_rows: list[dict[str, object]] = []
    fasta_rows: list[dict[str, object]] = []
    mimag_rows: list[dict[str, object]] = []
    biosample_rows: list[dict[str, object]] = []
    assembly_rows: list[dict[str, object]] = []
    taxonomy_request: list[dict[str, object]] = []

    review_fasta_dir = work / "review-fasta"
    for number, catalog in enumerate(review_candidates, start=1):
        isolate = f"MBPMAG{number:03d}"
        sgb, mag = catalog["SGB"], catalog["MAG"]
        quality = quality44[mag]
        graph = graph44[mag]
        taxonomy = taxonomy46[sgb]
        source = frozen[45] / "representative-genomes" / f"{catalog['Genome']}.gz"
        if not source.is_file():
            raise FileNotFoundError(source)
        if gunzip_sha256(source) != catalog["SHA256"]:
            raise ValueError(f"Article 45 representative checksum mismatch: {mag}")
        if catalog["SHA256"] != quality["MAGFASTA_SHA256"]:
            raise ValueError(f"Article 44/45 representative identity mismatch: {mag}")

        contigs = sorted(contigs_by_mag[mag], key=lambda row: row["Contig"])
        if len(contigs) != int(quality["Contigs"]):
            raise ValueError(f"Contig membership mismatch: {mag}")
        gc = [float(row["GCPct"]) for row in contigs]
        log_depth_1 = [math.log2(float(depth[(row["Contig"], "MOCK1")]["JgiMeanDepth"]) + 0.1) for row in contigs]
        log_depth_2 = [math.log2(float(depth[(row["Contig"], "MOCK2")]["JgiMeanDepth"]) + 0.1) for row in contigs]
        z_gc, z_d1, z_d2 = robust_z(gc), robust_z(log_depth_1), robust_z(log_depth_2)
        outlier_bp = 0
        outlier_count = 0
        for index, row in enumerate(contigs):
            flags = []
            if abs(z_gc[index]) > 3.5:
                flags.append("GC")
            if abs(z_d1[index]) > 3.5:
                flags.append("MOCK1 depth")
            if abs(z_d2[index]) > 3.5:
                flags.append("MOCK2 depth")
            flagged = bool(flags)
            if flagged:
                outlier_count += 1
                outlier_bp += int(row["LengthBp"])
            anomaly_rows.append({
                "Isolate": isolate,
                "SGB": sgb,
                "MAG": mag,
                "Contig": row["Contig"],
                "LengthBp": int(row["LengthBp"]),
                "GCPct": float(row["GCPct"]),
                "MOCK1MeanDepth": float(depth[(row["Contig"], "MOCK1")]["JgiMeanDepth"]),
                "MOCK2MeanDepth": float(depth[(row["Contig"], "MOCK2")]["JgiMeanDepth"]),
                "RobustZGC": round(z_gc[index], 6),
                "RobustZLog2DepthMOCK1": round(z_d1[index], 6),
                "RobustZLog2DepthMOCK2": round(z_d2[index], 6),
                "AutomatedOutlier": flagged,
                "OutlierSignals": "; ".join(flags),
                "AutomatedAction": "REVIEW_ONLY—DO_NOT_REMOVE_AUTOMATICALLY" if flagged else "NONE",
            })

        target = review_fasta_dir / f"{isolate}.fsa.gz"
        fasta_audit = write_review_fasta(source, target, isolate)
        fasta_rows.append({
            "Isolate": isolate,
            "SGB": sgb,
            "SourceRepresentative": catalog["Genome"],
            "SourceSequenceSHA256": catalog["SHA256"],
            "ReviewFile": target.name,
            **fasta_audit,
            "HeaderQualifier": f"[SRA={','.join(RUNS)}]",
            "SubmissionStatus": "DEMONSTRATION_REVIEW_FILE—DO_NOT_SUBMIT",
        })

        curation_rows.append({
            "Isolate": isolate,
            "SGB": sgb,
            "MAG": mag,
            "Contigs": int(quality["Contigs"]),
            "GenomeBp": int(quality["BinBp"]),
            "AutomatedOutlierContigs": outlier_count,
            "AutomatedOutlierBpPct": round(100 * outlier_bp / int(quality["BinBp"]), 6),
            "GUNCPass": as_bool(quality["GUNCPass"]),
            "GraphContinuity": graph["GraphContinuity"],
            "GraphComponents": int(graph["K141Components"]),
            "GraphBoundaryPct": float(graph["PairedBoundaryPct"]),
            "CoverageGCAutomatedScreen": "PASS_WITH_REVIEW_FLAGS" if outlier_count else "PASS_NO_ROBUST_OUTLIERS",
            "AnvioManualSignoff": "PENDING_INVESTIGATOR_REVIEW",
            "BandageManualSignoff": "PENDING_INVESTIGATOR_REVIEW",
            "ContigRemovalPerformed": False,
            "ExternalSubmissionReady": False,
        })

        mean_depths = {
            sample: float(coverage[(sgb, sample)]["MeanDepth"])
            for sample in ("MOCK1", "MOCK2")
        }
        breadth = {
            sample: float(coverage[(sgb, sample)]["CoveredFractionPct"])
            for sample in ("MOCK1", "MOCK2")
        }
        mimag_rows.append({
            "Isolate": isolate,
            "SGB": sgb,
            "GTDBRelease": "R232",
            "GTDBTaxonomy": taxonomy["GTDBTaxonomy"],
            "GTDBSpeciesAssigned": as_bool(taxonomy["SpeciesAssigned"]),
            "CheckM2Completeness": float(quality["CheckM2Completeness"]),
            "CheckM2Contamination": float(quality["CheckM2Contamination"]),
            "GUNCPass": as_bool(quality["GUNCPass"]),
            "Contigs": int(quality["Contigs"]),
            "GenomeBp": int(quality["BinBp"]),
            "N50Bp": int(quality["N50Bp"]),
            "GCPct": float(quality["GCPct"]),
            "Complete5S": as_bool(quality["Complete5S"]),
            "Complete16S": as_bool(quality["Complete16S"]),
            "Complete23S": as_bool(quality["Complete23S"]),
            "TRNAGenes": int(quality["TRNAGenes"]),
            "TRNAIsotypes": int(quality["TRNAIsotypes"]),
            "CodingDensityPct": float(quality["CodingDensityPct"]),
            "MIMAGQuality": quality["MIMAGQuality"],
            "MOCK1MeanDepth": mean_depths["MOCK1"],
            "MOCK2MeanDepth": mean_depths["MOCK2"],
            "MOCK1BreadthPct": breadth["MOCK1"],
            "MOCK2BreadthPct": breadth["MOCK2"],
            "DerivedFrom": ",".join(RUNS),
            "ManualCurationStatus": "PENDING_INVESTIGATOR_REVIEW",
        })
        taxonomy_request.append({
            "isolate": isolate,
            "unmodified_GTDB_R232_lineage": taxonomy["GTDBTaxonomy"],
        })
        biosample_rows.append({
            "sample_name": isolate,
            "organism": "AWAITING_NCBI_TAXONOMY_COORDINATION",
            "isolate": isolate,
            "collection_date": MISSING_SOURCE,
            "env_broad_scale": MISSING_SOURCE,
            "env_local_scale": MISSING_SOURCE,
            "env_medium": MISSING_SOURCE,
            "geo_loc_name": "not applicable—synthetic mock community",
            "isolation_source": "ZymoBIOMICS Microbial Community Standard II (Log Distribution), catalog D6311",
            "lat_lon": "not applicable—synthetic mock community",
            "sample_type": "metagenomic assembly",
            "derived_from": ",".join(RUNS),
            "metagenome_source": "synthetic microbial community",
            "package": "MIMAG 6.0 review draft",
            "review_status": "INCOMPLETE—DO_NOT_UPLOAD",
        })
        assembly_rows.append({
            "file_name": target.name,
            "organism": "AWAITING_NCBI_TAXONOMY_COORDINATION",
            "isolate": isolate,
            "BioProject": "UNREGISTERED",
            "BioSample": "UNREGISTERED",
            "assembly_method": "MEGAHIT 1.2.9; binette 1.0.5",
            "genome_coverage": f"MOCK1 {mean_depths['MOCK1']:.3f}x; MOCK2 {mean_depths['MOCK2']:.3f}x",
            "sequencing_technology": "Illumina paired-end shotgun metagenomics",
            "review_status": "INCOMPLETE—DO_NOT_UPLOAD",
        })

    checklist = [
        {"Gate": "95%-ANI representative coordinate fixed", "Status": "PASS", "Evidence": "catalog-disposition.tsv; Article 45 checksums"},
        {"Gate": "CheckM or CheckM2 completeness >=90%", "Status": "PASS", "Evidence": "12/12 technical-review MAGs"},
        {"Gate": "Total sequence >=100,000 nt", "Status": "PASS", "Evidence": "12/12 technical-review MAGs"},
        {"Gate": "GUNC and complete MIMAG audit", "Status": "PASS", "Evidence": "mimag-quality-supplement.tsv; Article 44"},
        {"Gate": "Coverage/GC anomaly screen", "Status": "PASS_REVIEW_FLAGS_RETAINED", "Evidence": "contig-anomaly-audit.tsv.gz"},
        {"Gate": "anvi'o and Bandage manual signoff", "Status": "BLOCKED", "Evidence": "manual-review-sheet.tsv"},
        {"Gate": "Investigator confirms ownership and non-duplicate submission", "Status": "BLOCKED", "Evidence": "public tutorial mock data are not a new submitter-owned genome study"},
        {"Gate": "NCBI organism names coordinated from GTDB lineages", "Status": "BLOCKED", "Evidence": "taxonomy-name-request.tsv is a draft; no NCBI reply recorded"},
        {"Gate": "BioProject and per-MAG BioSample accessions", "Status": "BLOCKED", "Evidence": "no accessions invented"},
        {"Gate": "Mandatory source metadata investigator-confirmed", "Status": "BLOCKED", "Evidence": "biosample-review-draft.tsv explicitly marks missing fields"},
    ]

    summary_dir = work / "summary"
    summary_dir.mkdir()
    write_tsv(summary_dir / "catalog-disposition.tsv", dispositions)
    write_tsv(summary_dir / "manual-review-sheet.tsv", curation_rows)
    write_tsv(summary_dir / "contig-anomaly-audit.tsv.gz", anomaly_rows)
    write_tsv(summary_dir / "review-fasta-manifest.tsv", fasta_rows)
    write_tsv(summary_dir / "mimag-quality-supplement.tsv", mimag_rows)
    write_tsv(summary_dir / "taxonomy-name-request.tsv", taxonomy_request)
    write_tsv(summary_dir / "biosample-review-draft.tsv", biosample_rows)
    write_tsv(summary_dir / "genome-batch-review-draft.tsv", assembly_rows)
    write_tsv(summary_dir / "submission-readiness-checklist.tsv", checklist)
    write_tsv(work / "input-lineage.tsv", [{
        "Output": "Article 49 pre-submission review package",
        "ImmediateInput": "Articles 41/43/44/45/46/48 checksum-covered evidence",
        "Transformation": "catalog gate, robust GC/depth screen, deterministic header rewrite, metadata completeness audit",
        "TruthUsed": "No—mock reference identities do not drive curation or submission gates",
        "Evidence": "catalog-disposition.tsv; contig-anomaly-audit.tsv.gz; review-fasta-manifest.tsv",
    }])
    dump_json(work / "run-contract.json", {
        "article": 49,
        "seed": 20260749,
        "random_process": False,
        "catalog_sgbs": 24,
        "technical_review_set_rule": {
            "article44_complete_audit": True,
            "checkm2_completeness_min_pct": 90,
            "total_sequence_min_nt": 100_000,
        },
        "contig_outlier_rule": "absolute robust z > 3.5 for GC or log2(mean depth + 0.1)",
        "outlier_action": "review only; never automatic deletion",
        "mimag_package": "NCBI MIMAG 6.0",
        "derived_from_runs": list(RUNS),
        "external_submission_performed": False,
        "accessions_invented": False,
        "truth_used_for_curation_or_submission": False,
    })
    dump_json(summary_dir / "run-summary.json", {
        "article": 49,
        "catalog_sgbs": 24,
        "article44_complete_audit_representatives": sum(row["Article44CompleteAudit"] for row in dispositions),
        "technical_review_set": len(review_candidates),
        "held_numeric_gate": sum(row["Disposition"] == "HOLD_NUMERIC_GENOME_GATE" for row in dispositions),
        "held_full_audit_missing": sum(row["Disposition"] == "HOLD_FULL_MAG_AUDIT_MISSING" for row in dispositions),
        "automated_outlier_contigs": sum(row["AutomatedOutlier"] for row in anomaly_rows),
        "manual_signoffs_complete": 0,
        "external_submission_ready": 0,
        "external_submission_performed": False,
        "accessions_invented": False,
    })
    (work / ".article49-summary-complete").write_text("complete\n", encoding="utf-8")
    print(f"Prepared {len(review_candidates)} MAG review files; external-ready MAGs: 0")


if __name__ == "__main__":
    main()
