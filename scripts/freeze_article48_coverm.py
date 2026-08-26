#!/usr/bin/env python3
"""Freeze compact checksum-covered Article 48 CoverM evidence."""

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
    root, work, output = args.project_root.resolve(), args.work_dir.resolve(), args.output_dir.resolve()
    sources = [
        root / "scripts/article41_44_utils.py",
        root / "scripts/article42_44_validation_utils.py",
        root / "scripts/prepare_article48_coverm.py",
        root / "scripts/run_article48_coverm.py",
        root / "scripts/summarize_article48_coverm.py",
        root / "scripts/plot_article48_coverm.R",
        root / "scripts/freeze_article48_coverm.py",
        root / "scripts/validate_article48_coverm.py",
        root / "env/assembly.yml",
        root / "env/assembly-linux-64.lock",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 48 source files: " + ", ".join(missing))
    extra = [
        (work / "genome-ledger.tsv", Path("genome-ledger.tsv")),
        (work / "samples.tsv", Path("samples.tsv")),
        (work / "raw/coverm-identity95.tsv", Path("raw/coverm-identity95.tsv")),
        (work / "raw/coverm-identity97.tsv", Path("raw/coverm-identity97.tsv")),
        (
            root / "data/small/45-drep-dereplication-frozen/file-checksums.sha256",
            Path("upstream/article45-file-checksums.sha256"),
        ),
        (
            root / "data/small/41-read-mapping-depth-frozen/input-audit.tsv",
            Path("upstream/article41-input-audit.tsv"),
        ),
    ]
    freeze_compact_bundle(
        root=root, work=work, output=output, article=48, slug="mag-abundance-coverm",
        source_files=sources, extra_files=extra,
    )
    print(f"Frozen Article 48 evidence: {output}")


if __name__ == "__main__":
    main()
