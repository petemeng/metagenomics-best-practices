#!/usr/bin/env python3
"""Freeze compact checksum-covered Article 37 evidence."""

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
    "article37_40_utils.py", "prepare_article37_cazymes.py", "run_article37_cazymes.py",
    "summarize_article37_cazymes.py", "freeze_article37_cazymes.py", "validate_article37_40.py", "validate_article37_cazymes.py",
)] + [root / "env/cazyme.yml", root / "env/cazyme-linux-64.lock", root / "data/small/37-dbcan-database-manifest.tsv"]
freeze_bundle(root=root, work=work, output=output, article=37, slug="cazymes-dbcan", source_files=files,
              selected_raw=[
                  (work / "catalog/overview.tsv", Path("raw/catalog-overview.tsv")),
                  (work / "btheta/cgc_standard_out.tsv", Path("raw/btheta-cgc-standard.tsv")),
                  (work / "btheta/substrate_prediction.tsv", Path("raw/btheta-substrate-prediction.tsv")),
              ])
