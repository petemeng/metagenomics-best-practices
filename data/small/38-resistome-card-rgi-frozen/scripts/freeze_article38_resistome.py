#!/usr/bin/env python3
"""Freeze compact checksum-covered Article 38 evidence."""

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
    "article37_40_utils.py", "prepare_article38_resistome.py", "run_article38_resistome.py",
    "summarize_article38_resistome.py", "freeze_article38_resistome.py", "validate_article37_40.py", "validate_article38_resistome.py",
)] + [root / "env/resistome.yml", root / "env/resistome-linux-64.lock", root / "env/deeparg-legacy.yml"]
freeze_bundle(root=root, work=work, output=output, article=38, slug="resistome-card-rgi", source_files=files)
