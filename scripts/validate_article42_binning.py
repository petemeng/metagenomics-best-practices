#!/usr/bin/env python3
"""Fail-closed validation for Article 42 binner comparison evidence."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from article42_44_validation_utils import Audit, as_bool, audit_chapter, audit_checksums, audit_figures, finish, read_tsv


FIGURES = ("42-binner-quality-yield", "42-recovery-purity", "42-single-vs-multisample", "42-taxonomy-coverage")
BRANCHES = {"MetaBAT2-MOCK1-only", "MetaBAT2-multisample", "SemiBin2-self-supervised", "VAMB-taxonomy-free", "TaxVAMB-Kraken2"}


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
    bins = read_tsv(frozen / "bin-quality-truth-audit.tsv")
    membership = read_tsv(frozen / "bin-membership.tsv.gz")
    branches = read_tsv(frozen / "binner-summary.tsv")
    truth = read_tsv(frozen / "truth-contig-assignment.tsv.gz")
    taxonomy = read_tsv(frozen / "taxonomy-summary.tsv")
    tools = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    commands = read_tsv(frozen / "command-log.tsv")
    qc_commands = read_tsv(frozen / "qc-command-log.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    raw_checkm = read_tsv(frozen / "raw/checkm2-quality-report.tsv")
    raw_gunc = read_tsv(frozen / "raw/gunc-maxcss.tsv")

    audit.add("Identity", "article", summary.get("article") == 42 and contract.get("article") == 42, summary.get("article"))
    audit.add("Identity", "seed", contract.get("seed") == 20260742 and summary.get("seed") == 20260742, contract.get("seed"))
    audit.add("Identity", "truth-blinded", contract.get("truth_blinding") is True and summary["boundaries"]["truth_never_used_by_binners"] is True, summary["boundaries"])
    coordinate = contract["coordinate_set"]
    for key, expected in {"Contigs": 10203, "TotalBp": 74932939, "MinimumBp": 1500, "MaximumBp": 1064594, "N50Bp": 23167}.items():
        audit.add("Coordinate", key, int(coordinate[key]) == expected, coordinate[key])
    audit.add("Branch", "five-branches", {row["Branch"] for row in branches} == BRANCHES, [row["Branch"] for row in branches])
    audit.add("Version", "tools", tools == {"MetaBAT2": "2.18", "SemiBin2": "2.3.0", "Vamb": "5.0.4", "Kraken2": "2.17.1"}, tools)
    audit.add("Truth", "87-genomes", summary["truth_genomes"] == 87, summary["truth_genomes"])
    audit.add("Truth", "coordinate-rows", len(truth) == 10203 and len({row["Contig"] for row in truth}) == 10203, len(truth))
    audit.add("Taxonomy", "contig-conservation", sum(int(row["Contigs"]) for row in taxonomy) == 10203, taxonomy)
    audit.add("Taxonomy", "unclassified-retained", any(row.get("DeepestRank") == "unclassified" or row.get("KrakenStatus") == "U" for row in taxonomy), taxonomy)

    by_branch_contigs: dict[str, list[str]] = defaultdict(list)
    for row in membership:
        by_branch_contigs[row["Branch"]].append(row["Contig"])
    for branch in sorted(BRANCHES):
        rows = [row for row in bins if row["Branch"] == branch]
        contigs = by_branch_contigs[branch]
        audit.add("Branch", f"{branch}-bins-positive", len(rows) > 0, len(rows))
        audit.add("Branch", f"{branch}-membership-unique", len(contigs) == len(set(contigs)), len(contigs))
        reported = next(row for row in branches if row["Branch"] == branch)
        audit.add("Branch", f"{branch}-bin-count", int(reported["Bins"]) == len(rows), reported["Bins"])
        audit.add("Branch", f"{branch}-contig-count", int(reported["BinnedContigs"]) == len(contigs), reported["BinnedContigs"])
        audit.add("Branch", f"{branch}-minimum-bin", all(int(row["BinBp"]) >= 200000 for row in rows), min(int(row["BinBp"]) for row in rows))
    candidate_ids = [row["CandidateID"] for row in bins]
    audit.add("Candidate", "global-identifiers", len(candidate_ids) == len(set(candidate_ids)), len(candidate_ids))
    audit.add("Candidate", "membership-coverage", set(candidate_ids) == {row["CandidateID"] for row in membership}, len(membership))
    for row in bins:
        expected = float(row["CheckM2Completeness"]) >= 50 and float(row["CheckM2Contamination"]) < 10 and as_bool(row["GUNCPass"])
        audit.add("QCFormula", row["CandidateID"], as_bool(row["QCMinimumPass"]) == expected, row["QCMinimumPass"])
    audit.add("QC", "checkm-row-count", len(raw_checkm) == len(bins), len(raw_checkm))
    audit.add("QC", "gunc-row-count", len(raw_gunc) == len(bins), len(raw_gunc))
    audit.add("QC", "all-values-finite", all(row["CheckM2Completeness"] not in {"", "nan"} and row["CheckM2Contamination"] not in {"", "nan"} for row in bins), len(bins))
    audit.add("Execution", "binner-commands", len(commands) == 7 and all(int(row["ReturnCode"]) == 0 for row in commands), len(commands))
    audit.add("Execution", "qc-commands", len(qc_commands) == 2 and all(int(row["ReturnCode"]) == 0 for row in qc_commands), len(qc_commands))
    audit.add("Execution", "resources-success", all(int(row["ExitStatus"]) == 0 for row in resources), len(resources))
    semibin_resource = next(row for row in resources if row["Label"] == "semibin2-self-supervised")
    audit.add("Execution", "semibin-resource-disclosure", semibin_resource["MeasurementStatus"].startswith("wall from timestamped"), semibin_resource)
    audit.add("Execution", "seeded-binners", all("20260742" in row["Command"] for row in commands if row["Label"] in {"metabat2-MOCK1-only", "metabat2-multisample", "semibin2-self-supervised", "vamb-taxonomy-free", "taxvamb-kraken2"}), "commands")
    single = read_tsv(frozen / "single-vs-multisample.tsv")[0]
    audit.add("Design", "ari-bounded", -1 <= float(single["AdjustedRandOnSharedBinnedContigs"]) <= 1, single)
    audit.add("Design", "shared-positive", int(single["SharedBinnedContigs"]) > 0, single["SharedBinnedContigs"])
    audit_chapter(args.chapter.resolve(), audit, article=42, figure_stems=FIGURES, tokens=("MetaBAT2 v2.18", "SemiBin2 v2.3.0", "VAMB v5.0.4", "TaxVAMB", "Kraken2 Standard-8", "truth-blind"))
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(article=42, audit=audit, output=args.output_dir.resolve(), payload={"candidate_bins": len(bins), "branches": len(branches)})


if __name__ == "__main__":
    raise SystemExit(main())
