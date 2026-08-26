#!/usr/bin/env python3
"""Fail-closed validation for Article 57 plasmids and ARG mobility."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
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
    "57-reference-plasmid-audit",
    "57-coassembly-plasmid-evidence",
    "57-arg-mobility-context",
    "57-usa300-positive-control",
    "57-mobility-evidence-ladder",
)


def near(value: object, expected: float, tolerance: float = 1e-6) -> bool:
    try:
        return math.isclose(float(value), expected, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def valid_tilde_fences(text: str) -> bool:
    opened = False
    for line in text.splitlines():
        if not line.startswith("~~~"):
            continue
        if opened:
            if line.strip() != "~~~":
                return False
            opened = False
        else:
            if line.strip() not in {"~~~bash", "~~~r", "~~~text"}:
                return False
            opened = True
    return not opened


def audit_chapter(chapter: Path, frozen: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    plot = (frozen / "scripts/plot_article57_plasmids.R").read_text(encoding="utf-8")
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "seed": "20260757" in text,
        "inline-theme": all(
            token in text
            for token in ("pal_pub <-", "scale_color_pub <-", "scale_fill_pub", "theme_pub <-", "save_pub <-")
        ),
        "resource-contract": all(token in text for token in ("CPU", "RAM", "磁盘", "耗时")),
        "real-data": all(
            token in text
            for token in (
                "PRJEB52977",
                "GCA_000013465.1",
                "ERR9765746",
                "ERR9765747",
                "a429a3724d4593f35b8d7323b20252a6be90e1cd",
            )
        ),
        "frozen-input": "data/small/57-plasmids-mobile-elements-frozen" in text,
        "version-contract": all(
            token in text
            for token in ("geNomad 1.12.0", "database v1.9", "RGI 6.0.8", "CARD 4.0.1", "minimap2 2.31-r1302")
        ),
        "mobility-boundary": all(
            token in text
            for token in (
                "同一 contig",
                "不等于已经发生转移",
                "conjugation_genes",
                "CONJscan",
                "MOB-suite",
                "IntegronFinder",
                "实验",
            )
        ),
        "calibration-boundary": all(
            token in text for token in ("score ≥ 0.7", "score calibration", "FDR", "NA")
        ),
        "methods-results": all(token in text for token in ("Methods template", "Results template")),
        "citations": all(
            token in text
            for token in (
                "@camargo2024genomad",
                "@meslier2022platforms",
                "@alcock2023card",
                "@pike2022broadhostmge",
                "@partridge2018mgeamr",
                "@robertson2018mobsuite",
                "@cury2016integronfinder",
            )
        ),
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "plot-labels-english": re.search(r"[\u3400-\u9fff]", plot) is None,
        "code-fences": valid_tilde_fences(text),
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z][A-Z0-9_]*__|\bTODO\b|\bTBD\b|\bNNN\b", text) is None,
        "no-meta-prose": not any(
            token in text
            for token in ("本篇可独立", "本文可独立", "全系列约定", "接口只学一次", "作者代码通常长这样", "（即本文）", "无头服务器")
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
        checks[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
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

    run = json.loads((frozen / "run-contract.json").read_text(encoding="utf-8"))
    exact_run = {
        "article": 57,
        "seed": 20260757,
        "benchmark_commit": "a429a3724d4593f35b8d7323b20252a6be90e1cd",
        "reference_replicons": 399,
        "reference_plasmid_labels": 43,
        "reference_benchmark_replicons": 86,
        "coassembly_contigs": 18354,
        "coassembly_bases": 84811518,
        "rgi_primary_coassembly_calls": 34,
        "rgi_primary_staphylococcus_calls": 21,
        "genomad_min_score": 0.7,
        "genomad_score_calibration": False,
        "random_output_requested": False,
    }
    for key, expected in exact_run.items():
        observed = run.get(key)
        status = near(observed, expected) if isinstance(expected, float) else observed == expected
        audit.add("Run contract", key, status, {"expected": expected, "observed": observed})

    assets = read_tsv(frozen / "asset-check-audit.tsv")
    audit.add("Input", "asset-count", len(assets) == 4, len(assets))
    audit.add("Input", "all-checksums-pass", all(as_bool(row["ChecksumPass"]) for row in assets), assets)
    labels = read_tsv(frozen / "reference-replicon-labels.tsv")
    audit.add("Input", "reference-sequences", len(labels) == 399, len(labels))
    audit.add("Input", "reference-plasmid-labels", Counter(row["ReferenceLabel"] for row in labels) == Counter({"Other replicon": 356, "Plasmid": 43}), Counter(row["ReferenceLabel"] for row in labels))
    audit.add(
        "Input",
        "reference-sequence-class-present",
        all(row.get("ReferenceSequenceClass") for row in labels),
        Counter(row.get("ReferenceSequenceClass") for row in labels),
    )

    benchmark = read_tsv(frozen / "reference-benchmark-labels.tsv")
    pairs: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in benchmark:
        pairs[row["PairID"]].append(row)
    audit.add("Input", "matched-benchmark-size", len(benchmark) == 86 and len(pairs) == 43, {"rows": len(benchmark), "pairs": len(pairs)})
    audit.add(
        "Input",
        "matched-benchmark-balance",
        all(Counter(row["ReferenceLabel"] for row in rows) == Counter({"Plasmid": 1, "Other replicon": 1}) for rows in pairs.values()),
        "one labelled plasmid and one other sequence per pair",
    )

    summary = json.loads((frozen / "summary.json").read_text(encoding="utf-8"))
    audit.add("Summary", "reference-inventory", summary["reference_replicons"] == 399 and summary["reference_benchmark_replicons"] == 86 and summary["reference_plasmids"] == 43, summary)
    audit.add("Summary", "confusion-denominators", summary["reference_tp"] + summary["reference_fn"] == 43 and summary["reference_fp"] + summary["reference_tn"] == 43, summary)
    audit.add("Summary", "coassembly-arg-count", summary["coassembly_primary_arg_calls"] == 34, summary)
    audit.add("Summary", "usa300-reference-plasmids", summary["usa300_reference_plasmids"] == 3, summary)
    audit.add("Summary", "usa300-plasmid-arg-positive-control", summary["usa300_primary_arg_calls_on_reference_plasmids"] == 3, summary)

    reference_calls = read_tsv(frozen / "reference-classification.tsv")
    confusion = read_tsv(frozen / "reference-confusion.tsv")
    metrics = read_tsv(frozen / "reference-metrics.tsv")
    audit.add("Table", "reference-call-rows", len(reference_calls) == 86, len(reference_calls))
    audit.add("Table", "confusion-total", sum(int(row["Count"]) for row in confusion) == 86, confusion)
    metric_names = {row["Metric"] for row in metrics}
    audit.add(
        "Table",
        "reference-metric-set",
        metric_names == {"Sensitivity", "Precision", "Specificity"},
        sorted(metric_names),
    )

    candidates = read_tsv(frozen / "coassembly-plasmid-candidates.tsv")
    audit.add("Table", "candidate-count", len(candidates) == summary["coassembly_plasmid_candidates"], len(candidates))
    audit.add("Table", "candidate-score-gate", all(float(row["PlasmidScore"]) >= 0.7 for row in candidates), [row["PlasmidScore"] for row in candidates])
    audit.add(
        "Table",
        "uncalibrated-fdr",
        all(row["FDR"] in {"", "NA"} for row in candidates),
        sorted({row["FDR"] for row in candidates}),
    )
    audit.add("Table", "claim-ceiling-populated", all(row["ClaimCeiling"] for row in candidates), len(candidates))
    support_counts = Counter(row["ReferenceSupport"] for row in candidates)
    audit.add(
        "Table",
        "candidate-reference-audit",
        support_counts == Counter({
            "Reference-plasmid supported": 70,
            "Complete-cellular-reference conflict": 120,
            "No high-coverage reference support": 28,
        }),
        support_counts,
    )

    args_table = read_tsv(frozen / "arg-mobility-ledger.tsv")
    audit.add("Table", "arg-ledger-count", len(args_table) == 34, len(args_table))
    audit.add("Table", "arg-tier-domain", {row["EvidenceTier"] for row in args_table} <= {"Perfect", "Strict"}, sorted({row["EvidenceTier"] for row in args_table}))
    audit.add("Table", "arg-context-populated", all(row["EvidenceContext"] and row["ClaimCeiling"] for row in args_table), len(args_table))

    usa = read_tsv(frozen / "usa300-replicon-audit.tsv")
    audit.add("Table", "usa300-four-replicons", len(usa) == 4, usa)
    audit.add("Table", "usa300-three-labelled-plasmids", sum(row["ReferenceLabel"] == "Plasmid" for row in usa) == 3, usa)
    audit.add("Table", "usa300-three-plasmid-arg-calls", sum(int(row["CARDPrimaryARGCount"]) for row in usa if row["ReferenceLabel"] == "Plasmid") == 3, usa)

    ladder = read_tsv(frozen / "mobility-evidence-ladder.tsv")
    audit.add("Table", "evidence-ladder", len(ladder) == 6 and [int(row["Rank"]) for row in ladder] == list(range(1, 7)), ladder)
    audit.add("Table", "experimental-top-rank", ladder[0]["Evidence"] == "Observed transfer", ladder[0])

    resources = read_tsv(frozen / "resource-usage.tsv")
    audit.add("Resource", "four-measured-branches", len(resources) == 4, resources)
    audit.add("Resource", "positive-runtime-and-memory", all(float(row["WallSeconds"]) > 0 and float(row["PeakRSSGiB"]) > 0 for row in resources), resources)

    versions = {row["Tool"]: row["Version"] for row in read_tsv(frozen / "tool-versions.tsv")}
    audit.add("Environment", "tool-versions", versions == {"geNomad": "1.12.0", "geNomad database": "1.9 (ICTV MSL39)", "RGI": "6.0.8", "CARD": "4.0.1", "minimap2": "2.31-r1302"}, versions)
    virus_lock = (frozen / "env/virus-discovery-linux-64.lock").read_text(encoding="utf-8")
    resistome_lock = (frozen / "env/resistome-linux-64.lock").read_text(encoding="utf-8")
    assembly_lock = (frozen / "env/assembly-linux-64.lock").read_text(encoding="utf-8")
    audit.add("Environment", "explicit-locks", all("@EXPLICIT" in lock for lock in (virus_lock, resistome_lock, assembly_lock)), "three Linux exact locks")
    audit.add("Environment", "genomad-lock", "genomad-1.12.0" in virus_lock, "geNomad 1.12.0")
    audit.add("Environment", "rgi-lock", "rgi-6.0.8" in resistome_lock, "RGI 6.0.8")
    audit.add("Environment", "minimap2-lock", "minimap2-2.31" in assembly_lock, "minimap2 2.31")

    audit_chapter(args.chapter.resolve(), frozen, audit)
    audit_figures(args.figure_dir.resolve(), audit, FIGURES)
    return finish(
        article=57,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "frozen_dir": str(frozen),
            "chapter": str(args.chapter.resolve()),
            "figures": list(FIGURES),
            "reference_benchmark": "43 labelled plasmids + 43 deterministic length-matched other sequences",
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
