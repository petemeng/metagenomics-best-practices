#!/usr/bin/env python3
"""Fail-closed validation for Article 43 bin-refinement evidence."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from article42_44_validation_utils import Audit, as_bool, audit_chapter, audit_checksums, audit_figures, finish, read_tsv


FIGURES = ("43-refinement-yield", "43-quality-landscape", "43-refinement-provenance", "43-method-selection")


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
    quality = read_tsv(frozen / "refinement-quality-truth-audit.tsv")
    membership = read_tsv(frozen / "refined-bin-membership.tsv.gz")
    methods = read_tsv(frozen / "refinement-summary.tsv")
    selected = read_tsv(frozen / "selected-mag-candidates.tsv")
    selected_membership = read_tsv(frozen / "selected-refinement-membership.tsv.gz")
    selection = read_tsv(frozen / "final-method-selection.tsv")[0]
    provenance = read_tsv(frozen / "refinement-provenance.tsv")
    links = read_tsv(frozen / "refinement-provenance-links.tsv.gz")
    tools = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    reconstruction = read_tsv(frozen / "candidate-reconstruction-audit.tsv")
    input_sets = read_tsv(frozen / "input-binsets.tsv")
    raw_checkm = read_tsv(frozen / "raw/checkm2-quality-report.tsv")
    raw_gunc = read_tsv(frozen / "raw/gunc-maxcss.tsv")
    resources = read_tsv(frozen / "resource-summary.tsv")
    commands = read_tsv(frozen / "command-log.tsv")

    audit.add("Identity", "article", summary.get("article") == 43 and contract.get("article") == 43, summary.get("article"))
    audit.add("Identity", "seed", summary.get("seed") == 20260743 and contract.get("seed") == 20260743, contract.get("seed"))
    audit.add("Identity", "truth-blinded", contract.get("truth_blinding") is True and summary.get("truth_used_for_selection") is False, summary)
    audit.add("Version", "tools", tools == {"DAS Tool": "1.1.7", "Binette": "1.2.1", "Prodigal": "2.6.3", "CheckM2": "1.1.0"}, tools)
    audit.add("Input", "five-binsets", len(input_sets) == 5 and sum(int(row["Bins"]) for row in input_sets) == len(reconstruction), input_sets)
    audit.add("Input", "reconstruction-pass", all(row["Status"] == "PASS" for row in reconstruction), len(reconstruction))
    audit.add("Output", "two-methods", {row["Method"] for row in methods} == {"DAS Tool", "Binette"}, methods)
    audit.add("Output", "qc-row-count", len(raw_checkm) == len(quality) == len(raw_gunc), len(quality))
    ids = [row["RefinedID"] for row in quality]
    audit.add("Output", "global-identifiers", len(ids) == len(set(ids)), len(ids))
    by_method: dict[str, list[str]] = defaultdict(list)
    for row in membership:
        by_method[row["Method"]].append(row["Contig"])
    for method in ("DAS Tool", "Binette"):
        rows = [row for row in quality if row["Method"] == method]
        report = next(row for row in methods if row["Method"] == method)
        audit.add("Method", f"{method}-bins", int(report["Bins"]) == len(rows) > 0, report["Bins"])
        audit.add("Method", f"{method}-partition", len(by_method[method]) == len(set(by_method[method])), len(by_method[method]))
        audit.add("Method", f"{method}-pass-count", int(report["MinimumPassBins"]) == sum(as_bool(row["ReferenceFreeMinimumPass"]) for row in rows), report["MinimumPassBins"])
    for row in quality:
        expected = float(row["CheckM2Completeness"]) >= 50 and float(row["CheckM2Contamination"]) < 10 and as_bool(row["GUNCPass"])
        audit.add("QCFormula", row["RefinedID"], as_bool(row["ReferenceFreeMinimumPass"]) == expected, row["ReferenceFreeMinimumPass"])
    ranked = sorted(methods, key=lambda row: (int(row["MinimumPassBins"]), float(row["PassingScoreSum"]), row["Method"] == "Binette"), reverse=True)
    audit.add("Selection", "rule-recomputed", selection["SelectedMethod"] == ranked[0]["Method"] == summary["selected_method"], ranked)
    audit.add("Selection", "truth-not-used", selection["MockTruthUsed"] == "No", selection)
    audit.add("Selection", "selected-count", int(selection["SelectedBins"]) == len(selected) == summary["selected_bins"], len(selected))
    selected_ids = {row["RefinedID"] for row in selected}
    audit.add("Selection", "selected-pass-only", all(as_bool(row["ReferenceFreeMinimumPass"]) for row in selected), len(selected))
    audit.add("Selection", "selected-membership", {row["RefinedID"] for row in selected_membership} == selected_ids, len(selected_membership))
    audit.add("Selection", "selected-partition", len(selected_membership) == len({row["Contig"] for row in selected_membership}), len(selected_membership))
    audit.add("Provenance", "one-summary-per-bin", {row["RefinedID"] for row in provenance} == set(ids), len(provenance))
    audit.add("Provenance", "links-cover-bins", {row["RefinedID"] for row in links} == set(ids), len(links))
    audit.add("Provenance", "dominant-bounded", all(0 < float(row["DominantInputCoveragePct"]) <= 100 for row in provenance), len(provenance))
    audit.add("Execution", "five-commands", len(commands) == 5 and all(int(row["ReturnCode"]) == 0 for row in commands), len(commands))
    audit.add("Execution", "resources-success", all(int(row["ExitStatus"]) == 0 for row in resources), len(resources))
    audit_chapter(args.chapter.resolve(), audit, article=43, figure_stems=FIGURES, tokens=("DAS Tool 1.1.7", "Binette 1.2.1", "CheckM2 1.1.0", "GUNC 1.1.0", "truth-blind"))
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(article=43, audit=audit, output=args.output_dir.resolve(), payload={"refined_bins": len(quality), "selected_method": selection["SelectedMethod"], "selected_bins": len(selected)})


if __name__ == "__main__":
    raise SystemExit(main())
