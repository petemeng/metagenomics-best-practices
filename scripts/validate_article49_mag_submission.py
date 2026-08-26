#!/usr/bin/env python3
"""Fail-closed validation for Article 49 MAG curation/submission evidence."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

from article42_44_validation_utils import (
    Audit,
    as_bool,
    audit_checksums,
    audit_figures,
    finish,
    read_tsv,
)


FIGURES = (
    "49-submission-funnel",
    "49-curation-flags",
    "49-readiness-matrix",
    "49-mimag-quality",
)
RUNS = "ERR9765746,ERR9765747"
MISSING_SOURCE = "NOT_AVAILABLE_FROM_TUTORIAL_SOURCE—INVESTIGATOR_MUST_SUPPLY"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_review_fasta(path: Path, isolate: str) -> dict[str, object]:
    headers: list[str] = []
    sequences: list[str] = []
    current: list[str] = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
                headers.append(line)
            else:
                current.append(line)
        if current:
            sequences.append("".join(current))
    sequence_digest = hashlib.sha256()
    for sequence in sequences:
        sequence_digest.update(sequence.encode("ascii"))
    header_pattern = re.compile(
        rf"^{re.escape('>' + isolate)}_contig\d{{5}} \[SRA={RUNS}\]$"
    )
    return {
        "headers": headers,
        "sequences": sequences,
        "header_contract": len(headers) == len(sequences)
        and all(header_pattern.fullmatch(header) for header in headers),
        "unique_headers": len(headers) == len(set(headers)),
        "contigs": len(sequences),
        "total_bp": sum(map(len, sequences)),
        "sequence_sha256": sequence_digest.hexdigest(),
        "file_sha256": sha256(path),
    }


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot_script = (frozen / "scripts/plot_article49_mag_submission.R").read_text(
        encoding="utf-8"
    )
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260749" in text,
        "inline-theme": all(
            token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")
        ),
        "resource-contract": all(
            token in text for token in ("CPU", "RAM", "磁盘", "耗时", "没有服务器")
        ),
        "real-study": all(
            token in text for token in ("PRJEB52977", "ERR9765746", "ERR9765747")
        ),
        "frozen-input": "data/small/49-mag-curation-submission-frozen" in text,
        "lineage-inputs": all(
            token in text
            for token in (
                "41-read-mapping-depth-frozen",
                "43-bin-refinement-frozen",
                "44-mag-qc-mimag-graph-frozen",
                "45-drep-dereplication-frozen",
                "46-gtdbtk-taxonomy-frozen",
                "48-mag-abundance-coverm-frozen",
            )
        ),
        "methods-results": all(
            token in text for token in ("Methods template", "Results template")
        ),
        "curation-contract": all(
            token in text
            for token in (
                "anvi-interactive",
                "Bandage image",
                "REVIEW_ONLY",
                "没有自动删除任何 contig",
            )
        ),
        "submission-contract": all(
            token in text
            for token in (
                "MIMAG 6.0",
                "genomes@ncbi.nlm.nih.gov",
                "one MIMAG BioSample",
                "Webin-CLI",
                "-validate",
                "external_submission_performed=false",
                "accessions_invented=false",
            )
        ),
        "official-citations": all(
            token in text
            for token in (
                "@ncbi2026magfaq",
                "@ncbi2026mimag",
                "@ncbi2026biosampleorganism",
                "@ena2026magsubmission",
            )
        ),
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(
            token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')
        ),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot_script) is None,
        "no-source-theme": 'source("R/theme_pub.R")' not in text
        and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-meta-prose": not any(
            token in text
            for token in (
                "本篇可独立",
                "本文可独立",
                "全系列约定",
                "接口只学一次",
                "作者代码通常长这样",
                "（即本文）",
                "无头服务器",
            )
        ),
    }
    sections = (
        "这一步对应论文里的哪张图",
        "理论：",
        "准备工作",
        "可复制代码",
        "审计与升级",
        "出版级美化",
        "常见坑",
        "这段 Methods 怎么写",
        "换成你自己的数据怎么做",
        "参考",
    )
    for section in sections:
        checks[f"section-{section}"] = (
            re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
        )
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    args = parser.parse_args()
    frozen = args.frozen_dir.resolve()
    audit = Audit()
    audit_checksums(frozen, audit)

    summary = json.loads((frozen / "run-summary.json").read_text(encoding="utf-8"))
    contract = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    catalog = read_tsv(frozen / "catalog-disposition.tsv")
    curation = read_tsv(frozen / "manual-review-sheet.tsv")
    anomalies = read_tsv(frozen / "contig-anomaly-audit.tsv.gz")
    manifest = read_tsv(frozen / "review-fasta-manifest.tsv")
    mimag = read_tsv(frozen / "mimag-quality-supplement.tsv")
    taxonomy = read_tsv(frozen / "taxonomy-name-request.tsv")
    biosample = read_tsv(frozen / "biosample-review-draft.tsv")
    genome_batch = read_tsv(frozen / "genome-batch-review-draft.tsv")
    checklist = read_tsv(frozen / "submission-readiness-checklist.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    commands = read_tsv(frozen / "command-log.tsv")
    tools = read_tsv(frozen / "tool-versions.tsv")

    expected_isolates = [f"MBPMAG{i:03d}" for i in range(1, 13)]
    expected_review_sgbs = [
        "SGB_002",
        "SGB_003",
        "SGB_005",
        "SGB_007",
        "SGB_012",
        "SGB_016",
        "SGB_017",
        "SGB_019",
        "SGB_021",
        "SGB_022",
        "SGB_023",
        "SGB_024",
    ]

    audit.add(
        "Identity",
        "article",
        summary.get("article") == 49 and contract.get("article") == 49,
        summary,
    )
    audit.add("Identity", "seed", contract.get("seed") == 20260749, contract)
    audit.add(
        "Identity",
        "truth-blind",
        contract.get("truth_used_for_curation_or_submission") is False,
        contract,
    )
    audit.add(
        "Contract",
        "review-only",
        contract.get("random_process") is False
        and contract.get("outlier_action") == "review only; never automatic deletion"
        and "absolute robust z > 3.5" in contract.get("contig_outlier_rule", ""),
        contract,
    )
    audit.add(
        "Contract",
        "submission-false",
        contract.get("external_submission_performed") is False
        and contract.get("accessions_invented") is False,
        contract,
    )
    audit.add(
        "Contract",
        "technical-gate",
        contract.get("technical_review_set_rule")
        == {
            "article44_complete_audit": True,
            "checkm2_completeness_min_pct": 90,
            "total_sequence_min_nt": 100000,
        },
        contract,
    )

    dispositions = Counter(row["Disposition"] for row in catalog)
    audit.add(
        "Catalog",
        "coordinate",
        len(catalog) == 24
        and [row["SGB"] for row in catalog] == [f"SGB_{i:03d}" for i in range(1, 25)],
        len(catalog),
    )
    audit.add(
        "Catalog",
        "dispositions",
        dispositions
        == {
            "TECHNICAL_REVIEW_SET": 12,
            "HOLD_NUMERIC_GENOME_GATE": 4,
            "HOLD_FULL_MAG_AUDIT_MISSING": 8,
        },
        dispositions,
    )
    audit.add(
        "Catalog",
        "gate-formula",
        all(
            as_bool(row["NCBINumericGenomeGate"])
            == (float(row["Completeness"]) >= 90 and int(row["GenomeBp"]) >= 100000)
            for row in catalog
        )
        and all(
            (row["Disposition"] == "TECHNICAL_REVIEW_SET")
            == (
                as_bool(row["Article44CompleteAudit"])
                and as_bool(row["NCBINumericGenomeGate"])
            )
            for row in catalog
        ),
        len(catalog),
    )

    audit.add(
        "Curation",
        "coordinate",
        len(curation) == 12
        and [row["Isolate"] for row in curation] == expected_isolates
        and [row["SGB"] for row in curation] == expected_review_sgbs,
        len(curation),
    )
    audit.add(
        "Curation",
        "manual-pending",
        all(
            row["AnvioManualSignoff"] == "PENDING_INVESTIGATOR_REVIEW"
            and row["BandageManualSignoff"] == "PENDING_INVESTIGATOR_REVIEW"
            and not as_bool(row["ContigRemovalPerformed"])
            and not as_bool(row["ExternalSubmissionReady"])
            for row in curation
        ),
        len(curation),
    )
    audit.add(
        "Curation",
        "graph-fragmented",
        all(row["GraphContinuity"] == "Fragmented in k141 graph" for row in curation),
        len(curation),
    )
    audit.add(
        "Curation",
        "contig-total",
        sum(int(row["Contigs"]) for row in curation) == 1829,
        sum(int(row["Contigs"]) for row in curation),
    )
    audit.add(
        "Curation",
        "flag-summary",
        sum(int(row["AutomatedOutlierContigs"]) for row in curation) == 135
        and sum(int(row["AutomatedOutlierContigs"]) == 0 for row in curation) == 2,
        curation,
    )

    flagged = [row for row in anomalies if as_bool(row["AutomatedOutlier"])]
    audit.add(
        "Anomaly",
        "rows",
        len(anomalies) == 1829 and len(flagged) == 135,
        {"rows": len(anomalies), "flagged": len(flagged)},
    )
    audit.add(
        "Anomaly",
        "review-only",
        all(
            row["AutomatedAction"]
            == (
                "REVIEW_ONLY—DO_NOT_REMOVE_AUTOMATICALLY"
                if as_bool(row["AutomatedOutlier"])
                else "NONE"
            )
            for row in anomalies
        ),
        len(anomalies),
    )
    audit.add(
        "Anomaly",
        "flag-formula",
        all(
            as_bool(row["AutomatedOutlier"])
            == any(
                abs(float(row[key])) > 3.5
                for key in (
                    "RobustZGC",
                    "RobustZLog2DepthMOCK1",
                    "RobustZLog2DepthMOCK2",
                )
            )
            for row in anomalies
        ),
        len(anomalies),
    )
    audit.add(
        "Anomaly",
        "flagged-bp",
        sum(int(row["LengthBp"]) for row in flagged) == 706042,
        sum(int(row["LengthBp"]) for row in flagged),
    )

    audit.add(
        "FASTA",
        "manifest-coordinate",
        len(manifest) == 12
        and [row["Isolate"] for row in manifest] == expected_isolates
        and [row["SGB"] for row in manifest] == expected_review_sgbs,
        len(manifest),
    )
    for row in manifest:
        isolate = row["Isolate"]
        path = frozen / "review-fasta" / row["ReviewFile"]
        observed = inspect_review_fasta(path, isolate) if path.is_file() else {}
        audit.add("FASTA", f"{isolate}-exists", path.is_file(), str(path))
        audit.add(
            "FASTA",
            f"{isolate}-content",
            bool(observed)
            and observed["header_contract"]
            and observed["unique_headers"]
            and observed["contigs"] == int(row["Contigs"])
            and observed["total_bp"] == int(row["TotalBp"])
            and observed["sequence_sha256"] == row["SequenceOnlySHA256"]
            and observed["file_sha256"] == row["ReviewFASTA_SHA256"]
            and row["HeaderQualifier"] == f"[SRA={RUNS}]"
            and row["SubmissionStatus"]
            == "DEMONSTRATION_REVIEW_FILE—DO_NOT_SUBMIT",
            observed or "missing",
        )

    quality_counts = Counter(row["MIMAGQuality"] for row in mimag)
    audit.add(
        "MIMAG",
        "coordinate",
        len(mimag) == 12
        and [row["Isolate"] for row in mimag] == expected_isolates
        and [row["SGB"] for row in mimag] == expected_review_sgbs,
        len(mimag),
    )
    audit.add(
        "MIMAG",
        "quality-tiers",
        quality_counts == {"High quality": 4, "Medium quality": 8},
        quality_counts,
    )
    audit.add(
        "MIMAG",
        "numeric-gates",
        all(
            float(row["CheckM2Completeness"]) >= 90
            and int(row["GenomeBp"]) >= 100000
            and as_bool(row["GUNCPass"])
            and row["GTDBRelease"] == "R232"
            and as_bool(row["GTDBSpeciesAssigned"])
            and row["DerivedFrom"] == RUNS
            and row["ManualCurationStatus"] == "PENDING_INVESTIGATOR_REVIEW"
            for row in mimag
        ),
        len(mimag),
    )

    audit.add(
        "Metadata",
        "taxonomy-request",
        len(taxonomy) == 12
        and [row["isolate"] for row in taxonomy] == expected_isolates
        and all(
            row["unmodified_GTDB_R232_lineage"].startswith(("d__Bacteria;", "d__Archaea;"))
            and row["unmodified_GTDB_R232_lineage"].count(";") == 6
            for row in taxonomy
        ),
        len(taxonomy),
    )
    source_fields = (
        "collection_date",
        "env_broad_scale",
        "env_local_scale",
        "env_medium",
    )
    audit.add(
        "Metadata",
        "biosample-draft",
        len(biosample) == 12
        and [row["sample_name"] for row in biosample] == expected_isolates
        and all(row[field] == MISSING_SOURCE for row in biosample for field in source_fields)
        and all(
            row["organism"] == "AWAITING_NCBI_TAXONOMY_COORDINATION"
            and row["sample_type"] == "metagenomic assembly"
            and row["derived_from"] == RUNS
            and row["package"] == "MIMAG 6.0 review draft"
            and row["review_status"] == "INCOMPLETE—DO_NOT_UPLOAD"
            for row in biosample
        ),
        len(biosample),
    )
    audit.add(
        "Metadata",
        "genome-batch-draft",
        len(genome_batch) == 12
        and [row["isolate"] for row in genome_batch] == expected_isolates
        and all(
            row["BioProject"] == "UNREGISTERED"
            and row["BioSample"] == "UNREGISTERED"
            and row["review_status"] == "INCOMPLETE—DO_NOT_UPLOAD"
            for row in genome_batch
        ),
        len(genome_batch),
    )
    accession_pattern = re.compile(
        r"\b(?:PRJNA\d+|SAMN\d+|PRJEB\d+|SAMEA\d+|GCA_\d+(?:\.\d+)?)\b"
    )
    draft_text = "\n".join(
        "\t".join(row.values()) for table in (biosample, genome_batch) for row in table
    )
    audit.add(
        "Metadata",
        "no-invented-accessions",
        accession_pattern.search(draft_text) is None,
        accession_pattern.findall(draft_text),
    )

    checklist_counts = Counter(row["Status"] for row in checklist)
    audit.add(
        "Readiness",
        "ten-gates",
        len(checklist) == 10
        and checklist_counts
        == {"PASS": 4, "PASS_REVIEW_FLAGS_RETAINED": 1, "BLOCKED": 5},
        checklist_counts,
    )
    audit.add(
        "Readiness",
        "manual-and-accessions-blocked",
        any(
            row["Gate"] == "anvi'o and Bandage manual signoff"
            and row["Status"] == "BLOCKED"
            for row in checklist
        )
        and any(
            row["Gate"] == "BioProject and per-MAG BioSample accessions"
            and row["Status"] == "BLOCKED"
            for row in checklist
        ),
        checklist,
    )

    expected_summary = {
        "article": 49,
        "catalog_sgbs": 24,
        "article44_complete_audit_representatives": 16,
        "technical_review_set": 12,
        "held_numeric_gate": 4,
        "held_full_audit_missing": 8,
        "automated_outlier_contigs": 135,
        "manual_signoffs_complete": 0,
        "external_submission_ready": 0,
        "external_submission_performed": False,
        "accessions_invented": False,
    }
    audit.add("Result", "summary-exact", summary == expected_summary, summary)
    audit.add(
        "Execution",
        "command",
        len(commands) == 1
        and commands[0]["Label"] == "article49-prep"
        and int(commands[0]["ExitStatus"]) == 0
        and "prepare_article49_mag_submission.py" in commands[0]["Command"],
        commands,
    )
    audit.add(
        "Execution",
        "resource",
        len(resources) == 1
        and int(resources[0]["ExitStatus"]) == 0
        and float(resources[0]["WallSeconds"]) > 0
        and 0 < float(resources[0]["PeakRAMGiB"]) < 1,
        resources,
    )
    audit.add(
        "Execution",
        "tool-version",
        len(tools) == 1
        and tools[0]["Tool"] == "Python"
        and re.fullmatch(r"\d+\.\d+\.\d+", tools[0]["Version"]) is not None,
        tools,
    )
    upstream = sorted((frozen / "upstream").glob("article*-file-checksums.sha256"))
    audit.add(
        "Lineage",
        "six-upstream-manifests",
        [path.name for path in upstream]
        == [
            f"article{number}-file-checksums.sha256"
            for number in (41, 43, 44, 45, 46, 48)
        ]
        and all(path.stat().st_size > 0 for path in upstream),
        [path.name for path in upstream],
    )

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=49,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "catalog_sgbs": 24,
            "technical_review_set": 12,
            "automated_outlier_contigs": 135,
            "external_submission_ready": 0,
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
