#!/usr/bin/env python3
"""Freeze compact, checksum-covered Article 51 StrainPhlAn evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

from article41_44_utils import freeze_compact_bundle, read_tsv


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
        root / "scripts/prepare_article51_strainphlan.py",
        root / "scripts/run_article51_strainphlan.py",
        root / "scripts/summarize_article51_strainphlan.py",
        root / "scripts/plot_article51_strainphlan.R",
        root / "scripts/freeze_article51_strainphlan.py",
        root / "scripts/validate_article51_strainphlan.py",
        root / "env/biobakery.yml",
        root / "env/biobakery-linux-64.lock",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 51 source files: " + ", ".join(missing))

    paths = read_tsv(work / "output-paths.tsv")[0]
    extra = [
        (Path(paths["Info"]), Path("raw/strainphlan.info")),
        (Path(paths["Polymorphic"]), Path("raw/strainphlan.polymorphic.tsv")),
        (Path(paths["Alignment"]), Path("raw/strainphlan-concatenated.aln")),
        (Path(paths["Tree"]), Path("raw/strainphlan-raxml.tree")),
        (
            Path(paths["OfficialThresholdInfo"]),
            Path("raw/official-thresholds-strainphlan.info"),
        ),
        (
            Path(paths["OfficialThresholdPolymorphic"]),
            Path("raw/official-thresholds-strainphlan.polymorphic.tsv"),
        ),
        (
            Path(paths["OfficialThresholdAlignment"]),
            Path("raw/official-thresholds-strainphlan-concatenated.aln"),
        ),
        (
            Path(paths["OfficialThresholdTree"]),
            Path("raw/official-thresholds-strainphlan-raxml.tree"),
        ),
        (Path(paths["Configuration"]), Path("raw/phylophlan-seeded.cfg")),
        (
            work / "official_baseline/t__SGB4933_group.info",
            Path("raw/official-baseline.info"),
        ),
        (
            work / "official_baseline/RAxML_bestTree.t__SGB4933_group.StrainPhlAn4.tre",
            Path("raw/official-baseline.tree"),
        ),
    ]
    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=51,
        slug="strainphlan",
        source_files=sources,
        extra_files=extra,
    )
    print(f"Frozen Article 51 evidence: {output}")


if __name__ == "__main__":
    main()
