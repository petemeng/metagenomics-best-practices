#!/usr/bin/env python3
"""Freeze compact, checksum-covered Article 46 GTDB-Tk evidence."""

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
        root / "scripts/prepare_article46_gtdbtk.py",
        root / "scripts/run_article46_gtdbtk.py",
        root / "scripts/summarize_article46_gtdbtk.py",
        root / "scripts/plot_article46_gtdbtk.R",
        root / "scripts/freeze_article46_gtdbtk.py",
        root / "scripts/validate_article46_gtdbtk.py",
        root / "env/gtdbtk.yml",
        root / "env/gtdbtk-linux-64.lock",
        root / "db/download_db.sh",
        root / "data/small/46-gtdbtk-database-manifest.tsv",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 46 source files: " + ", ".join(missing))

    native_root = work / "gtdbtk"
    native_files: list[Path] = []
    for path in sorted(native_root.rglob("*")):
        if not path.is_file() or path.stat().st_size > 25 * 1024 * 1024:
            continue
        name = path.name
        relative = path.relative_to(native_root).as_posix()
        keep = (
            name.endswith((".summary.tsv", ".warnings.tsv", ".failed_genomes.tsv"))
            or name.endswith(".ani_summary.tsv")
            or name.endswith((".markers_summary.tsv", ".translation_table_summary.tsv"))
            or name.endswith(".user_msa.fasta.gz")
            or (name.endswith(".tree") and "classify" in relative)
            or name.endswith(".tree.mapping.tsv")
            or name in {"gtdbtk.log", "gtdbtk.warnings.log", "gtdbtk.json"}
        )
        if keep:
            native_files.append(path)
    if not any(path.name.endswith(".summary.tsv") for path in native_files):
        raise FileNotFoundError("No native GTDB-Tk summary table selected for freezing")

    extra = [
        (work / "genome-ledger.tsv", Path("genome-ledger.tsv")),
        (
            root / "data/small/45-drep-dereplication-frozen/file-checksums.sha256",
            Path("upstream/article45-file-checksums.sha256"),
        ),
    ]
    extra.extend(
        (path, Path("raw/gtdbtk") / path.relative_to(native_root))
        for path in native_files
    )
    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=46,
        slug="gtdbtk-taxonomy",
        source_files=sources,
        extra_files=extra,
    )
    print(f"Frozen Article 46 evidence: {output}")


if __name__ == "__main__":
    main()
