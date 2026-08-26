#!/usr/bin/env python3
"""Fail-closed validation for Article 62 element-cycle evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image

from article42_44_validation_utils import (
    Audit,
    audit_checksums,
    audit_figures,
    finish,
    read_tsv,
    sha256,
)


FIGURES = (
    "62-process-landscape",
    "62-environment-associations",
    "62-mag-carrier-map",
    "62-community-mag-concordance",
    "62-nitrogen-step-cooccurrence",
    "62-mag-recovery-ceiling",
)

DETERMINISTIC_ANALYSIS = (
    "analysis-contract.json",
    "analysis-metrics.json",
    "carrier-phylum-summary.tsv",
    "community-mag-concordance.tsv",
    "environment-associations.tsv",
    "mag-carrier-summary.tsv",
    "mag-process-evidence.tsv.gz",
    "mag-recovery-ceiling.tsv",
    "missing-ko-columns.tsv",
    "nitrogen-step-summary.tsv",
    "process-rules.tsv",
    "sample-carrier-fraction.tsv.gz",
    "sample-process-index.tsv.gz",
    "selected-sample-metadata.tsv",
    "software-versions.json",
    "spring-carrier-fraction.tsv",
    "spring-process-index.tsv",
    "temperature-regime-summary.tsv",
)

EXPECTED_RESOURCES = {
    "sample-metadata.tsv": (
        133160,
        "8dc11ab08a4c53e800038f3ab87958b085267417aa7ced06ce0b78a9985f1c37",
    ),
    "ko-proportions-in-metagenomes.tsv.gz": (
        15712433,
        "7a25f3bd22eef266061cbc8819eb87a84606c5c7e5a7eb625f301689e8e30c9d",
    ),
    "mag-metadata.tsv": (
        3750601,
        "fd667ac8afe505a894dfe0e621b64c63f44cb7e6a40ce695f86f54196d60825f",
    ),
    "mag-abundances-per-sample.biom": (
        2017689,
        "d3188968399136dfe595bacf9b2fda01568284a03945550f7ae0153eb4a7b131",
    ),
    "kos-in-mags.tsv.gz": (
        684507,
        "4733be481e754fcc7718ae4784dea507b8c7419ccda2744d3fe8138ff583197c",
    ),
    "diting-pathway-formulas-v0.3.txt": (
        70799,
        "96e32f16e1a260930334408dc00e0f7584adc4c3dc4eac678ce389f6eb74fecd",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--chapter", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stage-dir", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    return parser.parse_args()


def near(value: object, expected: float, tolerance: float = 1e-10) -> bool:
    try:
        return math.isclose(float(value), expected, abs_tol=tolerance, rel_tol=tolerance)
    except (TypeError, ValueError):
        return False


def image_pixel_sha(path: Path) -> str:
    with Image.open(path) as image:
        normalized = image.convert("RGBA")
        payload = (
            normalized.size[0].to_bytes(8, "little")
            + normalized.size[1].to_bytes(8, "little")
            + normalized.tobytes()
        )
    return hashlib.sha256(payload).hexdigest()


def audit_chapter_62(chapter: Path, audit: Audit) -> None:
    text = chapter.read_text(encoding="utf-8")
    sections = (
        "这一步对应论文里的哪张图",
        "理论",
        "准备工作",
        "可复制代码",
        "审计与升级",
        "出版级美化",
        "常见坑",
        "这段 Methods 怎么写",
        "换成你自己的数据怎么做",
        "参考",
    )
    checks = {
        "draft-false": "draft: false" in text,
        "eval-false": re.search(r"execute:\s*\n\s+eval:\s+false", text) is not None,
        "freeze-auto": "freeze: auto" in text,
        "bibliography": "bibliography: ../references.bib" in text,
        "analysis-seed": "62001" in text,
        "plot-seed": "20260762" in text,
        "inline-theme": all(token in text for token in ("pal_pub <-", "theme_pub <-", "save_pub <-")),
        "resource-contract": all(token in text for token in ("RAM", "CPU", "磁盘", "耗时")),
        "real-data": all(token in text for token in ("500", "56", "780", "30284068.v2")),
        "source-methods": all(token in text for token in ("Prodigal 2.6.3", "KOfam r111", "HMMER 3.3", "10^-10")),
        "inference-contract": all(token in text for token in ("9,999", "2,000", "温泉中位数", "已回收 MAG 池")),
        "ambiguity-audit": all(token in text for token in ("McrABG", "NarGH/NxrAB", "pmoCAB/amoCAB", "DsrAB")),
        "methods-template": "**Element-cycle profiling.**" in text,
        "results-template": "**Results template.**" in text,
        "frozen-input": "data/small/62-element-cycling-frozen" in text,
        "figures": all(f"../figures/{stem}.png" in text for stem in FIGURES),
        "vector-raster-export": all(token in text for token in (".pdf", ".png", ".tiff", 'compression = "lzw"')),
        "no-source-theme": 'source("R/theme_pub.R")' not in text and "source('R/theme_pub.R')" not in text,
        "no-placeholders": re.search(r"__[A-Z0-9_]+__|TODO|TBD|NNN|vX", text) is None,
        "no-planning-comment": "Planned chapter" not in text and "Do not publish" not in text,
        "no-meta-prose": not any(token in text for token in (
            "本篇可独立", "本文可独立", "全系列约定", "接口只学一次",
            "作者代码通常长这样", "（即本文）", "无头服务器",
        )),
    }
    for section in sections:
        checks[f"section-{section}"] = re.search(rf"(?m)^##\s+{re.escape(section)}", text) is not None
    for check, status in checks.items():
        audit.add("Chapter", check, status, check)


def stage_reanalysis(
    root: Path,
    frozen: Path,
    stage: Path,
    python: Path,
    audit: Audit,
) -> None:
    if stage.exists():
        shutil.rmtree(stage)
    cache = stage / "cache"
    analysis = stage / "analysis"
    figures = stage / "figures"
    shutil.copytree(frozen / "inputs", cache)
    environment = os.environ.copy()
    environment.update(
        {
            "MPLCONFIGDIR": str(stage / "matplotlib"),
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    run = subprocess.run(
        [
            str(python),
            str(frozen / "scripts/run_article62_element_cycles.py"),
            "--cache-dir",
            str(cache),
            "--rules",
            str(cache / "62-element-cycle-marker-rules.tsv"),
            "--output-dir",
            str(analysis),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    audit.add("Reanalysis", "analysis-exit", run.returncode == 0, (run.stdout + run.stderr)[-4000:])
    if run.returncode != 0:
        return
    for name in DETERMINISTIC_ANALYSIS:
        observed = sha256(analysis / name) if (analysis / name).is_file() else "MISSING"
        expected = sha256(frozen / name)
        audit.add("Reanalysis", f"hash-{name}", observed == expected, observed)

    plotted = subprocess.run(
        [
            str(python),
            str(frozen / "scripts/plot_article62_element_cycles.py"),
            "--analysis-dir",
            str(analysis),
            "--figure-dir",
            str(figures),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    audit.add("Reanalysis", "plot-exit", plotted.returncode == 0, (plotted.stdout + plotted.stderr)[-4000:])
    if plotted.returncode != 0:
        return
    audit_figures(figures, audit, FIGURES)
    for stem in FIGURES:
        staged = figures / f"{stem}.png"
        published = root / "figures" / f"{stem}.png"
        status = staged.is_file() and published.is_file()
        observed = image_pixel_sha(staged) if staged.is_file() else "MISSING"
        expected = image_pixel_sha(published) if published.is_file() else "MISSING"
        audit.add("Reanalysis", f"pixel-match-{stem}", status and observed == expected, observed)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    frozen = args.frozen_dir.resolve()
    figures = args.figure_dir.resolve()
    audit = Audit()

    audit_checksums(frozen, audit)
    contract = json.loads((frozen / "frozen-contract.json").read_text(encoding="utf-8"))
    audit.add("Contract", "article", contract.get("article") == 62, contract)
    audit.add("Contract", "dataset", contract.get("dataset_doi") == "10.6084/m9.figshare.30284068.v2", contract)
    audit.add("Contract", "dimensions", contract.get("samples") == 500 and contract.get("independent_hot_springs") == 56 and contract.get("mags") == 780, contract)
    audit.add("Contract", "processes", contract.get("processes") == 20, contract)
    audit.add("Contract", "seed", contract.get("seed") == 62001, contract)
    audit.add("Contract", "resampling", contract.get("permutations") == 9999 and contract.get("bootstraps") == 2000, contract)
    audit.add("Contract", "source-included", contract.get("source_tables_included") is True, contract)
    audit.add("Contract", "interpretation", contract.get("interpretation") == "genetic potential; no activity, direction, or rate claim", contract)

    manifest = json.loads((frozen / "inputs/download-manifest.json").read_text(encoding="utf-8"))
    records = {row["local_file"]: row for row in manifest["resources"]}
    audit.add("Input", "manifest-release", manifest.get("dataset_doi") == "10.6084/m9.figshare.30284068.v2" and manifest.get("figshare_version") == 2, manifest)
    for name, (size, digest) in EXPECTED_RESOURCES.items():
        path = frozen / "inputs" / name
        record = records.get(name, {})
        audit.add(
            "Input",
            f"locked-{name}",
            path.is_file()
            and path.stat().st_size == size
            and sha256(path) == digest
            and int(record.get("bytes", size)) == size
            and record.get("sha256") == digest,
            {"bytes": path.stat().st_size if path.is_file() else -1, "sha256": sha256(path) if path.is_file() else "MISSING"},
        )

    rules = read_tsv(frozen / "process-rules.tsv")
    audit.add("Rule", "twenty-processes", len(rules) == 20 and len({row["ProcessID"] for row in rules}) == 20, len(rules))
    audit.add("Rule", "locked-sha", sha256(frozen / "inputs/62-element-cycle-marker-rules.tsv") == "f0e2761ac2220c401be16b6d3d2607ce30fc2eadb8bf34c5c0a1bf9536b28dca", sha256(frozen / "inputs/62-element-cycle-marker-rules.tsv"))
    audit.add("Rule", "caveats-complete", all(row["DirectionCaveat"] and row["StrictCarrierDNF"] for row in rules), "all rows")

    metrics = json.loads((frozen / "analysis-metrics.json").read_text(encoding="utf-8"))
    audit.add("Result", "dimensions", metrics.get("samples") == 500 and metrics.get("hot_springs") == 56 and metrics.get("mags") == 780, metrics)
    audit.add("Result", "source-ko-widths", metrics.get("community_source_kos") == 10011 and metrics.get("mag_source_kos") == 6841 and metrics.get("required_marker_kos") == 68, metrics)
    audit.add("Result", "missing-marker-audit", metrics.get("missing_community_marker_columns") == 9 and metrics.get("missing_mag_marker_columns") == 10, metrics)
    audit.add("Result", "association-counts", metrics.get("significant_temperature_associations_fdr05") == 3 and metrics.get("significant_ph_associations_fdr05") == 7 and metrics.get("significant_community_mag_concordance_fdr05") == 14, metrics)
    audit.add("Result", "community-completeness-gate", metrics.get("community_process_sample_pairs_passing_complete_gate") == 6504, metrics)
    audit.add("Result", "carrier-calls", metrics.get("strict_carrier_calls") == 818 and metrics.get("relaxed_carrier_calls") == 852, metrics)
    audit.add("Result", "chain", metrics.get("complete_denitrification_chain_mags") == 1, metrics)
    audit.add("Result", "recruitment-median", near(metrics.get("median_mag_recruitment_rate"), 0.111735), metrics)

    ledger = read_tsv(frozen / "run-ledger.tsv")
    audit.add("Run", "nine-steps", len(ledger) == 9, len(ledger))
    audit.add("Run", "all-passed", all(row["Status"] == "passed" for row in ledger), ledger)
    associations = read_tsv(frozen / "environment-associations.tsv")
    audit.add("Run", "forty-associations", len(associations) == 40 and {row["Variable"] for row in associations} == {"Temperature", "pH"}, len(associations))
    audit.add("Run", "independent-unit", all(row["IndependentSprings"] == "56" for row in associations), "56 springs")
    audit.add("Run", "resampling-lock", all(row["Permutations"] == "9999" and row["Bootstraps"] == "2000" for row in associations), "9999/2000")
    audit.add("Run", "sample-process-rows", len(read_tsv(frozen / "sample-process-index.tsv.gz")) == 10_000, "500 x 20")
    audit.add("Run", "mag-evidence-rows", len(read_tsv(frozen / "mag-process-evidence.tsv.gz")) == 15_600, "780 x 20")
    audit.add("Run", "spring-process-rows", len(read_tsv(frozen / "spring-process-index.tsv")) == 1_120, "56 x 20")

    audit_figures(figures, audit, FIGURES)
    plot_source = (root / "scripts/plot_article62_element_cycles.py").read_text(encoding="utf-8")
    audit.add("Figure", "english-only-source", re.search(r"[\u4e00-\u9fff]", plot_source) is None, "plot source contains no CJK text")
    audit_chapter_62(args.chapter.resolve(), audit)
    stage_reanalysis(root, frozen, args.stage_dir.resolve(), args.python.resolve(), audit)

    return finish(
        article=62,
        audit=audit,
        output=args.output_dir.resolve(),
        payload={
            "samples": metrics["samples"],
            "hot_springs": metrics["hot_springs"],
            "mags": metrics["mags"],
            "figures": len(FIGURES),
        },
    )


if __name__ == "__main__":
    raise SystemExit(main())
