#!/usr/bin/env python3
"""Freeze compact checksum-covered Article 39 evidence."""

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
    "article37_40_utils.py", "build_vfdb_abricate_database.py", "prepare_article39_virulome.py",
    "run_article39_virulome.py", "summarize_article39_virulome.py", "freeze_article39_virulome.py",
    "validate_article37_40.py", "validate_article39_virulome.py",
)] + [root / "env/virulome.yml", root / "env/virulome-linux-64.lock"]
freeze_bundle(root=root, work=work, output=output, article=39, slug="virulome-vfdb-abricate", source_files=files)
