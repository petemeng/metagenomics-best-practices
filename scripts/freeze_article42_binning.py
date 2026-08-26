#!/usr/bin/env python3
"""Freeze compact checksum-covered Article 42 evidence."""

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
        root / "scripts/prepare_article42_binning.py",
        root / "scripts/convert_article42_kraken_taxonomy.py",
        root / "scripts/run_article42_binning.py",
        root / "scripts/summarize_article42_binning.py",
        root / "scripts/run_article42_candidate_qc.py",
        root / "scripts/plot_article42_binning.R",
        root / "scripts/freeze_article42_binning.py",
        root / "scripts/validate_article42_binning.py",
        root / "env/assembly.yml", root / "env/assembly-linux-64.lock",
        root / "env/binning.yml", root / "env/binning-linux-64.lock",
        root / "env/mag-qc.yml", root / "env/mag-qc-linux-64.lock",
        root / "env/gunc.yml", root / "env/gunc-linux-64.lock",
        root / "data/small/44-mag-qc-database-manifest.tsv",
    ]
    maxcss = sorted((work / "qc/gunc").glob("GUNC.*.maxCSS_level.tsv"))
    if len(maxcss) != 1:
        raise RuntimeError(f"Expected one Article 42 GUNC maxCSS table: {maxcss}")
    freeze_compact_bundle(
        root=root, work=work, output=output, article=42, slug="binning-comparison",
        source_files=sources,
        extra_files=[
            (work / "taxonomy/taxvamb-taxonomy.tsv", Path("raw/taxvamb-taxonomy.tsv")),
            (work / "qc/checkm2/quality_report.tsv", Path("raw/checkm2-quality-report.tsv")),
            (maxcss[0], Path("raw/gunc-maxcss.tsv")),
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
