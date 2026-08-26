#!/usr/bin/env python3
"""Freeze compact, checksum-covered virus-discovery evidence for Article 54."""

from __future__ import annotations

import argparse
from pathlib import Path

from article41_44_utils import freeze_compact_bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/54-virus-discovery-quality"),
    )
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    output = args.output_dir.resolve()
    results = args.results_dir
    if not results.is_absolute():
        results = (root / results).resolve()

    sources = [
        root / "scripts/article41_44_utils.py",
        root / "scripts/article42_44_validation_utils.py",
        root / "scripts/prepare_article54_virus_discovery.py",
        root / "scripts/run_article54_virus_discovery.py",
        root / "scripts/summarize_article54_virus_discovery.py",
        root / "scripts/plot_article54_virus_discovery.R",
        root / "scripts/freeze_article54_virus_discovery.py",
        root / "scripts/validate_article54_virus_discovery.py",
        root / "env/virus-discovery.yml",
        root / "env/virus-discovery-linux-64.lock",
        root / "env/virsorter2.yml",
        root / "env/virsorter2-linux-64.lock",
        root / "env/virsorter2-runtime-linux-64.lock",
        root / "env/condarc-cn-mirror.yml",
        root / "data/small/54-virus-database-manifest.tsv",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing Article 54 source files: " + ", ".join(missing)
        )

    extras = [
        (work / "input/PMC6871006.xml", Path("raw/PMC6871006.xml")),
        (
            work / "input/checkv-test-sequences.fna",
            Path("raw/checkv-test-sequences.fna"),
        ),
        (
            work / "input/checkv-upstream-ground-truth-quality-summary.tsv",
            Path("raw/checkv-upstream-ground-truth-quality-summary.tsv"),
        ),
        (
            results
            / "genomad/checkv-test-sequences_summary/checkv-test-sequences_virus_summary.tsv",
            Path("upstream/genomad-virus-summary.tsv"),
        ),
        (
            results
            / "genomad/checkv-test-sequences_aggregated_classification/checkv-test-sequences_aggregated_classification.tsv",
            Path("upstream/genomad-aggregated-classification.tsv"),
        ),
        (
            results / "virsorter2/final-viral-score.tsv",
            Path("upstream/virsorter2-final-viral-score.tsv"),
        ),
        (
            results / "virsorter2/final-viral-boundary.tsv",
            Path("upstream/virsorter2-final-viral-boundary.tsv"),
        ),
        (
            results / "virsorter2/iter-0/all-fullseq-proba.tsv",
            Path("upstream/virsorter2-all-fullseq-proba.tsv"),
        ),
        (
            results / "checkv/quality_summary.tsv",
            Path("upstream/checkv-quality-summary.tsv"),
        ),
        (
            results / "checkv/complete_genomes.tsv",
            Path("upstream/checkv-complete-genomes.tsv"),
        ),
        (
            results / "checkv/contamination.tsv",
            Path("upstream/checkv-contamination.tsv"),
        ),
        (
            results / "checkv/proviruses.fna",
            Path("upstream/checkv-proviruses.fna"),
        ),
        (
            results / "votu/all-vs-all-blast.tsv",
            Path("upstream/votu-all-vs-all-blast.tsv"),
        ),
        (results / "votu/ani.tsv", Path("upstream/votu-ani.tsv")),
        (
            results / "votu/votu-clusters.tsv",
            Path("upstream/votu-clusters.tsv"),
        ),
    ]
    missing_evidence = [
        str(source)
        for source, _ in extras
        if not source.is_file() or source.stat().st_size == 0
    ]
    if missing_evidence:
        raise FileNotFoundError(
            "Missing Article 54 evidence files: " + ", ".join(missing_evidence)
        )

    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=54,
        slug="virus-discovery-quality",
        source_files=sources,
        extra_files=extras,
    )
    print(f"Frozen Article 54 evidence: {output}")


if __name__ == "__main__":
    main()
