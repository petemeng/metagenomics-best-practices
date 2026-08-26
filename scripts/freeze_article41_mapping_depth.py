#!/usr/bin/env python3
"""Freeze compact checksum-covered Article 41 evidence."""

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
        root / "scripts/prepare_article41_mapping_depth.py",
        root / "scripts/run_article41_mapping_depth.py",
        root / "scripts/summarize_article41_mapping_depth.py",
        root / "scripts/plot_article41_mapping_depth.R",
        root / "scripts/freeze_article41_mapping_depth.py",
        root / "scripts/validate_article41_mapping_depth.py",
        root / "env/assembly.yml",
        root / "env/assembly-linux-64.lock",
    ]
    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=41,
        slug="read-mapping-depth",
        source_files=sources,
        extra_files=[
            (work / "depth/jgi-depth.tsv", Path("raw/jgi-depth.tsv")),
            (work / "depth/paired-contigs.tsv", Path("raw/paired-contigs.tsv")),
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
