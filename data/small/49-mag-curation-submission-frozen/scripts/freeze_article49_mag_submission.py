#!/usr/bin/env python3
"""Freeze checksum-covered Article 49 review evidence without fake accessions."""

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
    root, work, output = (
        args.project_root.resolve(), args.work_dir.resolve(), args.output_dir.resolve()
    )
    sources = [
        root / "scripts/article41_44_utils.py",
        root / "scripts/article42_44_validation_utils.py",
        root / "scripts/prepare_article49_mag_submission.py",
        root / "scripts/run_article49_mag_submission.py",
        root / "scripts/plot_article49_mag_submission.R",
        root / "scripts/freeze_article49_mag_submission.py",
        root / "scripts/validate_article49_mag_submission.py",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 49 source files: " + ", ".join(missing))
    extra: list[tuple[Path, Path]] = []
    for path in sorted((work / "review-fasta").glob("*.fsa.gz")):
        extra.append((path, Path("review-fasta") / path.name))
    for number, slug in {
        41: "41-read-mapping-depth-frozen",
        43: "43-bin-refinement-frozen",
        44: "44-mag-qc-mimag-graph-frozen",
        45: "45-drep-dereplication-frozen",
        46: "46-gtdbtk-taxonomy-frozen",
        48: "48-mag-abundance-coverm-frozen",
    }.items():
        extra.append((
            root / f"data/small/{slug}/file-checksums.sha256",
            Path("upstream") / f"article{number}-file-checksums.sha256",
        ))
    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=49,
        slug="mag-curation-submission",
        source_files=sources,
        extra_files=extra,
    )
    print(f"Frozen Article 49 evidence: {output}")


if __name__ == "__main__":
    main()
