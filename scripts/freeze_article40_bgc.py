#!/usr/bin/env python3
"""Freeze compact checksum-covered Article 40 evidence."""

import argparse
from pathlib import Path

from article37_40_utils import freeze_bundle


parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, required=True)
parser.add_argument("--work-dir", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
args = parser.parse_args()
root, work, output = args.project_root.resolve(), args.work_dir.resolve(), args.output_dir.resolve()
files = [root / "scripts" / name for name in (
    "article37_40_utils.py", "prepare_article40_bgc.py", "run_article40_bgc.py",
    "summarize_article40_bgc.py", "freeze_article40_bgc.py",
    "validate_article37_40.py", "validate_article40_bgc.py",
)] + [
    root / "env/bgc-gecco.yml",
    root / "env/antismash8.yml",
    root / "env/bgc-gecco-linux-64.lock",
    root / "env/antismash8-linux-64.lock",
    root / "data/small/40-bgc-database-manifest.tsv",
]
freeze_bundle(
    root=root,
    work=work,
    output=output,
    article=40,
    slug="bgc-natural-products",
    source_files=files,
    selected_raw=[
        (work / "inputs/salinispora-fragment-map.tsv", Path("raw/salinispora-fragment-map.tsv")),
        (work / "inputs/coassembly-id-map.tsv", Path("raw/coassembly-id-map.tsv")),
    ],
)
