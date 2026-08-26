#!/usr/bin/env python3
"""Summarize seven truth-aware Article 32 assembly branches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import statistics
from pathlib import Path


BRANCHES = [
    "spades-short-only",
    "spades-illumina-ont",
    "spades-illumina-hifi",
    "flye-ont",
    "flye-ont-polypolish-default",
    "flye-ont-polypolish-careful",
    "flye-hifi",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metaquast-dir", type=Path, required=True)
    parser.add_argument("--normalized-dir", type=Path, required=True)
    parser.add_argument("--truth-manifest", type=Path, required=True)
    parser.add_argument("--resource-dir", type=Path, required=True)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    fields = fields or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value: str | None):
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text in {"", "-", "NA", "None"}:
        return None
    parsed = float(text)
    return int(parsed) if parsed.is_integer() else parsed


def first_unaligned(value: str | None):
    if value is None:
        return None
    match = re.match(r"\s*([0-9,]+)", value)
    return int(match.group(1).replace(",", "")) if match else None


def parse_elapsed(value: str) -> float:
    parts = [float(part) for part in value.split(":")]
    if len(parts) == 3:
        return parts[0] * 3600 + parts[1] * 60 + parts[2]
    if len(parts) == 2:
        return parts[0] * 60 + parts[1]
    return parts[0]


def parse_time_file(path: Path) -> dict:
    values = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if ": " in line:
                key, value = line.split(": ", 1)
                values[key] = value
    elapsed_key = next((key for key in values if key.startswith("Elapsed (wall clock) time")), None)
    elapsed = parse_elapsed(values[elapsed_key]) if elapsed_key else None
    rss_kb = numeric(values.get("Maximum resident set size (kbytes)"))
    return {
        "Step": path.stem,
        "ElapsedSeconds": f"{elapsed:.3f}" if elapsed is not None else "NA",
        "PeakRSSGiB": f"{rss_kb / 1024 / 1024:.3f}" if rss_kb is not None else "NA",
        "CPUPercent": values.get("Percent of CPU this job got", "NA").rstrip("%"),
        "FileSystemInputs": numeric(values.get("File system inputs")) or 0,
        "FileSystemOutputs": numeric(values.get("File system outputs")) or 0,
        "ExitStatus": numeric(values.get("Exit status")),
    }


def parse_polypolish_log(path: Path, mode: str) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    input_lengths = [
        int(value.replace(",", ""))
        for value in re.findall(
            r"^Polishing .+ \(([0-9,]+) bp\):$", text, flags=re.MULTILINE
        )
    ]
    output_lengths = [
        int(value.replace(",", ""))
        for value in re.findall(
            r"^\s{2}.+_polypolish \(([0-9,]+) bp\)$",
            text,
            flags=re.MULTILINE,
        )
    ]
    changed = [
        int(value.replace(",", ""))
        for value in re.findall(
            r"^\s{2}([0-9,]+) positions changed", text, flags=re.MULTILINE
        )
    ]
    alignment_rows = re.findall(
        r":\s*([0-9,]+) alignments from ([0-9,]+) reads",
        text,
    )
    kept = re.findall(r"^\s{2}([0-9,]+) alignments kept$", text, flags=re.MULTILINE)
    discarded = re.findall(
        r"^\s{2}([0-9,]+) alignments discarded$", text, flags=re.MULTILINE
    )
    row = {
        "Mode": mode,
        "LogBytes": path.stat().st_size,
        "AlignmentFiles": len(alignment_rows),
        "InputAlignments": sum(
            int(alignments.replace(",", ""))
            for alignments, _ in alignment_rows
        ),
        "InputMappedReads": sum(
            int(reads.replace(",", "")) for _, reads in alignment_rows
        ),
        "HighQualityAlignmentsKept": int(kept[-1].replace(",", ""))
        if kept
        else "NA",
        "HighQualityAlignmentsDiscarded": int(discarded[-1].replace(",", ""))
        if discarded
        else "NA",
        "PolishedSequences": len(input_lengths),
        "InputBases": sum(input_lengths),
        "OutputBases": sum(output_lengths) if output_lengths else "NA",
        "ChangedPositions": sum(changed),
    }
    row["CarefulFlagVisible"] = "yes" if mode == "careful" else "no"
    return row


def fasta_sequence_inventory(path: Path) -> list[tuple[str, str, int, str]]:
    inventory: list[tuple[str, str, int, str]] = []
    identifier: str | None = None
    header: str | None = None
    sequence = bytearray()

    def finish_record() -> None:
        nonlocal identifier, header, sequence
        if identifier is None or header is None:
            return
        inventory.append(
            (identifier, header, len(sequence), hashlib.sha256(sequence).hexdigest())
        )
        identifier = None
        header = None
        sequence = bytearray()

    with path.open("rb") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                finish_record()
                header = line[1:].decode("utf-8")
                identifier = header.split(None, 1)[0]
            else:
                if identifier is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                sequence.extend(line.upper())
    finish_record()
    identifiers = [row[0] for row in inventory]
    if not inventory or len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Empty FASTA or duplicate sequence IDs: {path}")
    return inventory


def main() -> int:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    combined_path = args.metaquast_dir / "combined_reference" / "transposed_report.tsv"
    combined = read_tsv(combined_path)
    by_branch = {row["Assembly"]: row for row in combined if row["Assembly"] in BRANCHES}
    if set(by_branch) != set(BRANCHES):
        raise ValueError(f"Unexpected combined-reference branches: {sorted(by_branch)}")
    auxiliary = [row for row in combined if row["Assembly"] not in BRANCHES]
    allowed_auxiliary = {f"{branch}_broken" for branch in BRANCHES}
    unexpected_auxiliary = {row["Assembly"] for row in auxiliary} - allowed_auxiliary
    if unexpected_auxiliary:
        raise ValueError(
            f"Unexpected combined-reference auxiliary branches: {sorted(unexpected_auxiliary)}"
        )
    if auxiliary:
        split_rows = []
        for row in auxiliary:
            parent = row["Assembly"].removesuffix("_broken")
            split_rows.append(
                {
                    "Assembly": row["Assembly"],
                    "ParentBranch": parent,
                    "SequencesGe1kb": numeric(row.get("# contigs")),
                    "TotalLengthBp": numeric(row.get("Total length")),
                    "N50Bp": numeric(row.get("N50")),
                    "Misassemblies": numeric(row.get("# misassemblies")),
                    "GenomeFractionPct": numeric(row.get("Genome fraction (%)")),
                    "MismatchesPer100Kbp": numeric(row.get("# mismatches per 100 kbp")),
                    "IndelsPer100Kbp": numeric(row.get("# indels per 100 kbp")),
                }
            )
        write_tsv(output / "split-scaffold-sensitivity.tsv", split_rows)

    structure = {
        row["Branch"]: row
        for row in read_tsv(args.normalized_dir / "assembly-structure.tsv")
    }
    if set(structure) != set(BRANCHES):
        raise ValueError("Normalized assembly branches do not match the seven-branch contract")

    metrics = []
    for branch in BRANCHES:
        row = by_branch[branch]
        local = structure[branch]
        metrics.append(
            {
                "Branch": branch,
                "SequencesGe1kb": numeric(local["Sequences"]),
                "TotalLengthBp": numeric(local["TotalLengthBp"]),
                "LargestBp": numeric(local["LargestBp"]),
                "N50Bp": numeric(local["N50Bp"]),
                "L50": numeric(local["L50"]),
                "NBases": numeric(local["NBases"]),
                "Misassemblies": numeric(row.get("# misassemblies")),
                "MisassembledContigs": numeric(row.get("# misassembled contigs")),
                "FullyUnalignedContigs": first_unaligned(row.get("# unaligned contigs")),
                "UnalignedLengthBp": numeric(row.get("Unaligned length")),
                "GenomeFractionPct": numeric(row.get("Genome fraction (%)")),
                "DuplicationRatio": numeric(row.get("Duplication ratio")),
                "NsPer100Kbp": numeric(row.get("# N's per 100 kbp")),
                "MismatchesPer100Kbp": numeric(row.get("# mismatches per 100 kbp")),
                "IndelsPer100Kbp": numeric(row.get("# indels per 100 kbp")),
                "LargestAlignmentBp": numeric(row.get("Largest alignment")),
                "TotalAlignedLengthBp": numeric(row.get("Total aligned length")),
                "NA50Bp": numeric(row.get("NA50")),
                "LA50": numeric(row.get("LA50")),
                "AssemblySHA256": local["SHA256"],
            }
        )

    truth = {row["CurrentGenomeLabel"]: row for row in read_tsv(args.truth_manifest)}
    if len(truth) != 71:
        raise ValueError(f"Expected 71 truth genomes, found {len(truth)}")
    truth_by_quast_label = {
        re.sub(r"[^A-Za-z0-9_.-]", "_", label): (label, row)
        for label, row in truth.items()
    }
    per_genome = []
    report_dirs: dict[str, Path] = {}
    for reference_dir in sorted((args.metaquast_dir / "runs_per_reference").iterdir()):
        report = reference_dir / "transposed_report.tsv"
        if not reference_dir.is_dir() or not report.is_file():
            continue
        quast_label = reference_dir.name
        if quast_label in truth:
            label = quast_label
        elif quast_label in truth_by_quast_label:
            label = truth_by_quast_label[quast_label][0]
        else:
            raise ValueError(f"MetaQUAST reference is absent from truth manifest: {quast_label}")
        if label in report_dirs:
            raise ValueError(f"Duplicate MetaQUAST report for reference {label}")
        report_dirs[label] = reference_dir

    for label, truth_row in truth.items():
        reference_dir = report_dirs.get(label)
        report_rows = (
            {row["Assembly"]: row for row in read_tsv(reference_dir / "transposed_report.tsv")}
            if reference_dir is not None
            else {}
        )
        unexpected = set(report_rows) - set(BRANCHES) - {
            f"{branch}_broken" for branch in BRANCHES
        }
        if unexpected:
            raise ValueError(f"Unexpected branches for reference {label}: {sorted(unexpected)}")
        abundance = float(truth_row["ExpectedAbundancePct"])
        abundance_bin = "<0.1%" if abundance < 0.1 else ("0.1-<1%" if abundance < 1 else ">=1%")
        for branch in BRANCHES:
            row = report_rows.get(branch)
            fraction = numeric(row.get("Genome fraction (%)")) if row else None
            mismatch = (
                numeric(row.get("# mismatches per 100 kbp")) if row else None
            )
            indel = numeric(row.get("# indels per 100 kbp")) if row else None
            per_genome.append(
                {
                    "Reference": label,
                    "ExpectedAbundancePct": f"{abundance:.12g}",
                    "AbundanceBin": abundance_bin,
                    "Branch": branch,
                    "GenomeFractionPct": fraction if fraction is not None else 0,
                    "RecoveredGe90Pct": "yes" if fraction is not None and fraction >= 90 else "no",
                    "FullGenomeGe99Pct": "yes" if fraction is not None and fraction >= 99 else "no",
                    "MismatchesPer100Kbp": mismatch if mismatch is not None else "NA",
                    "IndelsPer100Kbp": indel if indel is not None else "NA",
                    "ReferenceReportPresent": "yes" if reference_dir is not None else "no",
                    "BranchReportPresent": "yes" if row is not None else "no",
                }
            )
    if len(per_genome) != 71 * len(BRANCHES):
        raise ValueError(f"Expected {71 * len(BRANCHES)} per-genome rows, found {len(per_genome)}")
    write_tsv(output / "per-genome-metaquast.tsv", per_genome)

    for metric in metrics:
        rows = [row for row in per_genome if row["Branch"] == metric["Branch"]]
        metric["RecoveredGenomesGe90Pct"] = sum(row["RecoveredGe90Pct"] == "yes" for row in rows)
        metric["FullGenomesGe99Pct"] = sum(row["FullGenomeGe99Pct"] == "yes" for row in rows)
    write_tsv(output / "branch-metrics.tsv", metrics)

    draft_inventory = fasta_sequence_inventory(
        args.normalized_dir / "flye-ont.ge1000.fasta"
    )
    draft_order = [row[0] for row in draft_inventory]
    draft_map = {row[0]: row[2:] for row in draft_inventory}
    polishing_sequence_rows = []
    for mode in ("default", "careful"):
        branch = f"flye-ont-polypolish-{mode}"
        output_inventory = fasta_sequence_inventory(
            args.normalized_dir / f"{branch}.ge1000.fasta"
        )
        output_order = [row[0] for row in output_inventory]
        if len(output_order) != len(set(output_order)):
            raise ValueError(
                f"Duplicate canonical Polypolish sequence IDs for {mode}"
            )
        output_map = {row[0]: row[2:] for row in output_inventory}
        expected_header_annotations = sum(
            header == f"{identifier} polypolish"
            for identifier, header, _, _ in output_inventory
        )
        shared = set(draft_map) & set(output_map)
        identical = sum(draft_map[name] == output_map[name] for name in shared)
        draft_total = sum(row[2] for row in draft_inventory)
        output_total = sum(row[2] for row in output_inventory)
        polishing_sequence_rows.append(
            {
                "Mode": mode,
                "Branch": branch,
                "DraftSequences": len(draft_inventory),
                "OutputSequences": len(output_inventory),
                "CanonicalSharedSequenceIDs": len(shared),
                "LostCanonicalSequenceIDs": len(set(draft_map) - set(output_map)),
                "NewCanonicalSequenceIDs": len(set(output_map) - set(draft_map)),
                "CanonicalSequenceIDSetEqual": "yes"
                if set(draft_map) == set(output_map)
                else "no",
                "CanonicalSequenceOrderEqual": "yes"
                if draft_order == output_order
                else "no",
                "ExpectedPolypolishHeaderAnnotations": expected_header_annotations,
                "UnexpectedOutputHeaders": len(output_inventory)
                - expected_header_annotations,
                "IdenticalSequences": identical,
                "ChangedSequences": len(shared) - identical,
                "DraftTotalBp": draft_total,
                "OutputTotalBp": output_total,
                "LengthDeltaBp": output_total - draft_total,
                "DraftAssemblySHA256": structure["flye-ont"]["SHA256"],
                "OutputAssemblySHA256": structure[branch]["SHA256"],
            }
        )
    write_tsv(
        output / "polishing-sequence-audit.tsv", polishing_sequence_rows
    )

    bins = []
    for branch in BRANCHES:
        for abundance_bin in ("<0.1%", "0.1-<1%", ">=1%"):
            rows = [
                row for row in per_genome
                if row["Branch"] == branch and row["AbundanceBin"] == abundance_bin
            ]
            bins.append(
                {
                    "Branch": branch,
                    "AbundanceBin": abundance_bin,
                    "TruthGenomes": len(rows),
                    "MedianGenomeFractionPct": f"{statistics.median(float(row['GenomeFractionPct']) for row in rows):.6f}",
                    "RecoveredGenomesGe90Pct": sum(row["RecoveredGe90Pct"] == "yes" for row in rows),
                    "FullGenomesGe99Pct": sum(row["FullGenomeGe99Pct"] == "yes" for row in rows),
                }
            )
    write_tsv(output / "abundance-bin-recovery.tsv", bins)

    resources = [parse_time_file(path) for path in sorted(args.resource_dir.glob("*.txt"))]
    write_tsv(output / "resource-usage.tsv", resources)
    polish_rows = []
    for mode in ("default", "careful"):
        log = args.log_dir / f"polypolish-{mode}.log"
        if log.is_file():
            polish_rows.append(parse_polypolish_log(log, mode))
    if polish_rows:
        write_tsv(output / "polypolish-log-audit.tsv", polish_rows)

    summary = {
        "status": "computed",
        "branches": len(BRANCHES),
        "truth_genomes": len(truth),
        "per_genome_rows": len(per_genome),
        "physical_reference_reports": len(report_dirs),
        "missing_reference_reports": len(truth) - len(report_dirs),
        "missing_branch_reports_zero_recovery": sum(
            row["BranchReportPresent"] == "no" for row in per_genome
        ),
        "polishing_sequence_audit_modes": len(polishing_sequence_rows),
        "split_scaffold_sensitivity_branches": len(auxiliary),
        "minimum_contig_bp": 1000,
        "minimum_alignment_bp": 500,
        "minimum_identity_pct": 97,
        "full_genome_threshold_pct": 99,
        "seed": 20260732,
    }
    (output / "run-summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
