#!/usr/bin/env python3
"""Freeze compact, checksum-covered transmission evidence for Article 53."""

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
        root / "scripts/prepare_article53_transmission.py",
        root / "scripts/run_article53_transmission.py",
        root / "scripts/summarize_article53_transmission.py",
        root / "scripts/plot_article53_transmission.R",
        root / "scripts/freeze_article53_transmission.py",
        root / "scripts/validate_article53_transmission.py",
        root / "env/strain-transmission.yml",
        root / "env/strain-transmission-linux-64.lock",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 53 source files: " + ", ".join(missing))
    extras = [
        (work / "source-provenance.tsv", Path("source-provenance.tsv")),
        (work / "input/sys001172080st8.xlsx", Path("raw/sys001172080st8.xlsx")),
        (work / "input/sys001172080st9.xlsx", Path("raw/sys001172080st9.xlsx")),
        (work / "input/PMC5264247.xml", Path("raw/PMC5264247.xml")),
        (
            work / "input/40168_2022_1251_MOESM7_ESM.xlsx",
            Path("raw/40168_2022_1251_MOESM7_ESM.xlsx"),
        ),
        (work / "input/PMC8951724.xml", Path("raw/PMC8951724.xml")),
    ]
    for source, _ in extras:
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(f"Missing Article 53 evidence: {source}")
    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=53,
        slug="strain-transmission",
        source_files=sources,
        extra_files=extras,
    )
    print(f"Frozen Article 53 evidence: {output}")


if __name__ == "__main__":
    main()

