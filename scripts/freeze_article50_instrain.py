#!/usr/bin/env python3
"""Freeze compact, checksum-covered Article 50 inStrain evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

from article41_44_utils import freeze_compact_bundle


def one_file(pattern: str, base: Path) -> Path:
    matches = sorted(base.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"Expected one {pattern} under {base}, observed {matches}")
    return matches[0]


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
        root / "scripts/prepare_article50_instrain.py",
        root / "scripts/run_article50_instrain.py",
        root / "scripts/summarize_article50_instrain.py",
        root / "scripts/plot_article50_instrain.R",
        root / "scripts/freeze_article50_instrain.py",
        root / "scripts/validate_article50_instrain.py",
        root / "env/instrain.yml",
        root / "env/instrain-linux-64.lock",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 50 source files: " + ", ".join(missing))

    extra = [
        (work / "genome-ledger.tsv", Path("genome-ledger.tsv")),
        (work / "sample-ledger.tsv", Path("sample-ledger.tsv")),
        (
            one_file("*_genome_info.tsv*", work / "profiles/MOCK1/output"),
            Path("raw/MOCK1-genome_info.tsv.gz"),
        ),
        (
            one_file("*_genome_info.tsv*", work / "profiles/MOCK2/output"),
            Path("raw/MOCK2-genome_info.tsv.gz"),
        ),
        (
            one_file(
                "*_genomeWide_compare.tsv*",
                work / "comparison/MOCK1-vs-MOCK2/output",
            ),
            Path("raw/MOCK1-vs-MOCK2-genomeWide_compare.tsv.gz"),
        ),
        (
            root / "data/small/45-drep-dereplication-frozen/file-checksums.sha256",
            Path("upstream/article45-file-checksums.sha256"),
        ),
        (
            root / "data/small/46-gtdbtk-taxonomy-frozen/file-checksums.sha256",
            Path("upstream/article46-file-checksums.sha256"),
        ),
        (
            root / "data/small/48-mag-abundance-coverm-frozen/file-checksums.sha256",
            Path("upstream/article48-file-checksums.sha256"),
        ),
    ]
    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=50,
        slug="instrain-microdiversity",
        source_files=sources,
        extra_files=extra,
    )
    print(f"Frozen Article 50 evidence: {output}")


if __name__ == "__main__":
    main()
