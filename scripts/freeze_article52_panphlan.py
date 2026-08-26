#!/usr/bin/env python3
"""Freeze compact, checksum-covered PanPhlAn evidence for Article 52."""

from __future__ import annotations

import argparse
from pathlib import Path

from article41_44_utils import freeze_compact_bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    output = args.output_dir.resolve()
    sources = [
        root / "scripts/article41_44_utils.py",
        root / "scripts/article42_44_validation_utils.py",
        root / "scripts/prepare_article52_panphlan.py",
        root / "scripts/run_article52_panphlan.py",
        root / "scripts/summarize_article52_panphlan.py",
        root / "scripts/plot_article52_panphlan.R",
        root / "scripts/freeze_article52_panphlan.py",
        root / "scripts/validate_article52_panphlan.py",
        root / "env/panphlan.yml",
        root / "env/panphlan-linux-64.lock",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 52 source files: " + ", ".join(missing))
    extras = [
        (work / "pangenome-file-manifest.tsv", Path("pangenome-file-manifest.tsv")),
        (work / "software-provenance.tsv", Path("software-provenance.tsv")),
        (work / "output-shape-audit.tsv", Path("output-shape-audit.tsv")),
        (work / "software/panphlan_profiling.py", Path("pinned-source/panphlan_profiling.py")),
        (work / "software/misc.py", Path("pinned-source/misc.py")),
        (work / "output/primary-index.tsv", Path("raw/primary-index.tsv")),
        (work / "output/sensitive-index.tsv", Path("raw/sensitive-index.tsv")),
    ]
    for source, _ in extras:
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(f"Missing Article 52 evidence: {source}")
    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=52,
        slug="panphlan-pangenome",
        source_files=sources,
        extra_files=extras,
    )
    print(f"Frozen Article 52 evidence: {output}")


if __name__ == "__main__":
    main()
