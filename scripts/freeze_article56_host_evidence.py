#!/usr/bin/env python3
"""Freeze compact, checksum-covered host-evidence assets for Article 56."""

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
        root / "scripts/prepare_article56_host_evidence.py",
        root / "scripts/plot_article56_host_evidence.R",
        root / "scripts/freeze_article56_host_evidence.py",
        root / "scripts/validate_article56_host_evidence.py",
        root / "env/virus-host.yml",
        root / "env/virus-host-paper.yml",
        root / "env/virus-host-linux-64.lock",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 56 sources: " + ", ".join(missing))

    extras = [
        (work / "input/PMC10155999.xml", Path("raw/PMC10155999.xml")),
        (work / "input/PMC6871006.xml", Path("raw/PMC6871006.xml")),
    ]
    missing_evidence = [
        str(source)
        for source, _ in extras
        if not source.is_file() or source.stat().st_size == 0
    ]
    if missing_evidence:
        raise FileNotFoundError(
            "Missing Article 56 evidence files: " + ", ".join(missing_evidence)
        )

    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=56,
        slug="virus-host-evidence",
        source_files=sources,
        extra_files=extras,
    )
    print(f"Frozen Article 56 evidence: {output}")


if __name__ == "__main__":
    main()
