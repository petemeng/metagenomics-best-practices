#!/usr/bin/env python3
"""Freeze compact, checksum-covered taxonomy and abundance evidence for Article 55."""

from __future__ import annotations

import argparse
from pathlib import Path

from article41_44_utils import freeze_compact_bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    results = args.results_dir.resolve()
    output = args.output_dir.resolve()

    sources = [
        root / "scripts/article41_44_utils.py",
        root / "scripts/article42_44_validation_utils.py",
        root / "scripts/prepare_article55_virus_abundance.py",
        root / "scripts/download_article55_illumina_fastq.sh",
        root / "scripts/run_article55_virus_abundance.py",
        root / "scripts/summarize_article55_virus_abundance.py",
        root / "scripts/plot_article55_virus_abundance.R",
        root / "scripts/freeze_article55_virus_abundance.py",
        root / "scripts/validate_article55_virus_abundance.py",
        root / "env/virus-abundance.yml",
        root / "env/virus-abundance-linux-64.lock",
        root / "env/virus-discovery.yml",
        root / "env/virus-discovery-linux-64.lock",
        root / "data/small/55-phage-reference-manifest.tsv",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 55 sources: " + ", ".join(missing))

    extras = [
        (work / "input/PMC10926689.xml", Path("raw/PMC10926689.xml")),
        (work / "input/PMC11521287.xml", Path("raw/PMC11521287.xml")),
        (
            work / "input/mgen-10-1198-s002.xlsx",
            Path("raw/mgen-10-1198-s002.xlsx"),
        ),
        (
            work / "input/ena-prjeb56639-filereport.tsv",
            Path("raw/ena-prjeb56639-filereport.tsv"),
        ),
        (
            work / "input/cook15-phage-reference.fna",
            Path("raw/cook15-phage-reference.fna"),
        ),
        (work / "phage-metadata.tsv", Path("source/phage-metadata.tsv")),
        (
            work / "published-coverage-depth.tsv",
            Path("source/published-coverage-depth.tsv"),
        ),
        (
            work / "illumina-library-metadata.tsv",
            Path("source/illumina-library-metadata.tsv"),
        ),
        (
            work / "ena-illumina-run-manifest.tsv",
            Path("source/ena-illumina-run-manifest.tsv"),
        ),
        (work / "reference-audit.tsv", Path("source/reference-audit.tsv")),
        (
            work / "author-source-assertions.tsv",
            Path("source/author-source-assertions.tsv"),
        ),
        (
            work / "lifecycle-source-assertions.tsv",
            Path("source/lifecycle-source-assertions.tsv"),
        ),
        (
            work / "bbmap-smoke-coverage.tsv",
            Path("source/bbmap-smoke-coverage.tsv"),
        ),
        (
            results
            / "genomad/cook15-phage-reference_summary/cook15-phage-reference_virus_summary.tsv",
            Path("upstream/genomad-virus-summary.tsv"),
        ),
        (
            results
            / "genomad/cook15-phage-reference_annotate/cook15-phage-reference_taxonomy.tsv",
            Path("upstream/genomad-taxonomy.tsv"),
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
            "Missing Article 55 evidence files: " + ", ".join(missing_evidence)
        )

    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=55,
        slug="virus-taxonomy-abundance",
        source_files=sources,
        extra_files=extras,
    )
    print(f"Frozen Article 55 evidence: {output}")


if __name__ == "__main__":
    main()
