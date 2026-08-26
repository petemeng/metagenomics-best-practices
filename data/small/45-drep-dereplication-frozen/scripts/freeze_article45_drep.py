#!/usr/bin/env python3
"""Freeze compact, checksum-covered Article 45 evidence."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

from article41_44_utils import checksum_tree, freeze_compact_bundle


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
        root / "scripts" / "article41_44_utils.py",
        root / "scripts" / "article42_44_validation_utils.py",
        root / "scripts" / "prepare_article45_drep.py",
        root / "scripts" / "run_article45_drep.py",
        root / "scripts" / "summarize_article45_drep.py",
        root / "scripts" / "plot_article45_drep.R",
        root / "scripts" / "freeze_article45_drep.py",
        root / "scripts" / "validate_article45_drep.py",
        root / "env" / "drep.yml",
        root / "env" / "drep-linux-64.lock",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 45 source files: " + ", ".join(missing))

    extra: list[tuple[Path, Path]] = []
    for branch in ("species95", "nearclone999"):
        tables = work / "drep" / branch / "data_tables"
        for name in ("Cdb.csv", "Wdb.csv", "Ndb.csv", "genomeInformation.csv"):
            extra.append((tables / name, Path("raw") / branch / name))
    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=45,
        slug="drep-dereplication",
        source_files=sources,
        extra_files=extra,
    )
    for stale in (output / "logs").glob("drep-strain99.*"):
        stale.unlink()
    representative_dir = output / "representative-genomes"
    representative_dir.mkdir()
    for source in sorted((work / "drep/species95/dereplicated_genomes").glob("*.fna")):
        target = representative_dir / f"{source.name}.gz"
        with source.open("rb") as input_handle, target.open("wb") as raw_output:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, compresslevel=9, mtime=0) as output_handle:
                shutil.copyfileobj(input_handle, output_handle)
    checksum_tree(output)
    print(f"Frozen Article 45 evidence: {output}")


if __name__ == "__main__":
    main()
