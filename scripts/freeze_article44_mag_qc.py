#!/usr/bin/env python3
"""Freeze compact checksum-covered Article 44 evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

from article41_44_utils import freeze_compact_bundle


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root, work, output = args.project_root.resolve(), args.work_dir.resolve(), args.output_dir.resolve()
    sources = [
        root / "scripts/article41_44_utils.py",
        root / "scripts/article42_44_validation_utils.py",
        root / "scripts/prepare_article44_mag_qc.py",
        root / "scripts/run_article44_mag_qc.py",
        root / "scripts/summarize_article44_mag_qc.py",
        root / "scripts/plot_article44_mag_qc.R",
        root / "scripts/freeze_article44_mag_qc.py",
        root / "scripts/validate_article44_mag_qc.py",
        root / "env/assembly.yml", root / "env/assembly-linux-64.lock",
        root / "env/mag-qc.yml", root / "env/mag-qc-linux-64.lock",
        root / "env/gunc.yml", root / "env/gunc-linux-64.lock",
        root / "env/checkm1.yml", root / "env/checkm1-linux-64.lock",
        root / "data/small/44-mag-qc-database-manifest.tsv",
    ]
    maxcss = sorted((work / "qc/gunc").glob("GUNC.*.maxCSS_level.tsv"))
    if len(maxcss) != 1:
        raise RuntimeError(f"Expected one Article 44 GUNC maxCSS table: {maxcss}")
    extra = [
        (work / "qc/checkm2/quality_report.tsv", Path("raw/checkm2-quality-report.tsv")),
        (maxcss[0], Path("raw/gunc-maxcss.tsv")),
        (work / "qc/checkm1-qa.tsv", Path("raw/checkm1-qa.tsv")),
    ]
    for feature in sorted((work / "features").glob("*/*")):
        if feature.is_file() and feature.suffix in {".gff", ".tsv", ".stats"} and feature.stat().st_size < 10_000_000:
            extra.append((feature, Path("raw/features") / feature.parent.name / feature.name))
    freeze_compact_bundle(
        root=root, work=work, output=output, article=44, slug="mag-qc-mimag-graph",
        source_files=sources, extra_files=extra,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
