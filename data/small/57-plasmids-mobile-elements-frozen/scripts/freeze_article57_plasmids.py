#!/usr/bin/env python3
"""Create the compact, checksum-covered Article 57 frozen bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    target = args.output_dir.resolve()
    if not (work / ".article57-summary-complete").is_file():
        raise FileNotFoundError("Article 57 summary sentinel is missing")
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)

    for path in sorted((work / "summary").glob("*")):
        if path.is_file():
            copy(path, target / path.name)
    for name in (
        "asset-check-audit.tsv",
        "input-lineage.tsv",
        "reference-replicon-labels.tsv",
        "reference-benchmark-labels.tsv",
        "rgi-primary-coassembly.tsv",
        "rgi-primary-staphylococcus.tsv",
        "run-contract.json",
        "command-log.tsv",
        "tool-versions.tsv",
        "database-audit.tsv",
    ):
        copy(work / name, target / name)
    copy(
        work / "results/coassembly-plasmid-to-reference.paf",
        target / "coassembly-plasmid-to-reference.paf",
    )

    for branch in ("reference-benchmark", "coassembly", "staphylococcus"):
        source_root = work / "results" / branch
        patterns = (
            "*_summary/*_plasmid_summary.tsv",
            "*_summary/*_plasmid_genes.tsv",
            "*_summary/*_plasmid.fna",
            "*_summary/*_summary.json",
        )
        matches: list[Path] = []
        for pattern in patterns:
            current = sorted(source_root.glob(pattern))
            if len(current) != 1:
                raise RuntimeError(f"Expected one {branch}/{pattern}: {current}")
            matches.extend(current)
        for path in matches:
            copy(path, target / "genomad" / branch / path.name)

    for path in sorted((work / "logs").glob("*")):
        if path.is_file() and not path.name.startswith("genomad-references."):
            copy(path, target / "logs" / path.name)

    script_names = (
        "download_article57_plasmid_sources.sh",
        "prepare_article57_plasmids.py",
        "run_article57_plasmids.py",
        "summarize_article57_plasmids.py",
        "plot_article57_plasmids.R",
        "freeze_article57_plasmids.py",
        "validate_article57_plasmids.py",
    )
    for name in script_names:
        copy(root / "scripts" / name, target / "scripts" / name)
    for name in (
        "virus-discovery.yml",
        "virus-discovery-linux-64.lock",
        "resistome.yml",
        "resistome-linux-64.lock",
        "assembly.yml",
        "assembly-linux-64.lock",
    ):
        copy(root / "env" / name, target / "env" / name)
    copy(
        root / "data/small/54-virus-database-manifest.tsv",
        target / "genomad-database-manifest.tsv",
    )

    notice = """Article 57 frozen evidence bundle

The complete PRJEB52977 FASTQ archives, full 292.8-Mb MOCK2 reference FASTA,
84.8-Mb co-assembly, CARD database, and geNomad database are not duplicated in
this bundle. Their public identities, immutable checksums, database releases,
commands, resource logs, and compact derived evidence tables are retained here.

The 43 plasmid labels are literal GenBank header labels in the fixed Meslier
benchmark repository. The 43 negative comparators are deterministically
length-matched sequences without the word "plasmid" in those headers; they are
not asserted to be biologically perfect negatives.
"""
    (target / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    contract = {
        "article": 57,
        "created_from": str(work.relative_to(root)),
        "large_inputs_included": False,
        "random_output_requested": False,
        "seed": 20260757,
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    files = sorted(
        path for path in target.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    lines = [f"{sha256(path)}  {path.relative_to(target).as_posix()}" for path in files]
    (target / "file-checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Frozen Article 57 bundle: {len(files)} payload files in {target}")


if __name__ == "__main__":
    main()
