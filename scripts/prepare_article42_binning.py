#!/usr/bin/env python3
"""Prepare common-coordinate inputs for the Article 42 binner comparison."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from article41_44_utils import dump_json, fasta_records, fasta_summary, sha256, write_tsv


EXPECTED = {
    "MetaBAT2": "2.18",
    "SemiBin2": "2.3.0",
    "Vamb": "5.0.4",
    "Kraken2": "2.17.1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--assembly-env", type=Path, required=True)
    parser.add_argument("--binning-env", type=Path, required=True)
    parser.add_argument("--kraken-db", type=Path, required=True)
    parser.add_argument("--minimum-contig", type=int, default=1500)
    return parser.parse_args()


def command_text(command: list[str], runtime: dict[str, str], allow_nonzero: bool = False) -> str:
    result = subprocess.run(command, capture_output=True, text=True, env=runtime)
    if result.returncode and not allow_nonzero:
        raise RuntimeError(f"Command failed: {command}\n{result.stderr}")
    return result.stdout + result.stderr


def parse_version(tool: str, text: str) -> str:
    patterns = {
        "MetaBAT2": r"version\s+2:(\d+(?:\.\d+)*)",
        "SemiBin2": r"(?m)^(\d+(?:\.\d+)*)$",
        "Vamb": r"Vamb\s+(\d+(?:\.\d+)*)",
        "Kraken2": r"Kraken version\s+(\d+(?:\.\d+)*)",
    }
    match = re.search(patterns[tool], text)
    if not match:
        raise ValueError(f"Cannot parse {tool} version from {text[:500]}")
    return match.group(1)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    assembly_env = args.assembly_env.resolve()
    binning_env = args.binning_env.resolve()
    kraken_db = args.kraken_db.resolve()
    if (work / ".article42-inputs-complete").is_file():
        print(f"Article 42 inputs already prepared: {work}")
        return 0
    if work.exists() and any(path.is_file() for path in work.rglob("*")):
        raise SystemExit(f"Refusing incomplete work directory containing files: {work}")
    for directory in ("inputs", "taxonomy", "bins", "logs", "summary", "tmp"):
        (work / directory).mkdir(parents=True, exist_ok=True)

    runtime = os.environ.copy()
    runtime.update(
        {
            "PATH": f"{binning_env / 'bin'}:{assembly_env / 'bin'}:/usr/bin:/bin",
            "LC_ALL": "C",
            "TZ": "UTC",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "TMPDIR": str(work / "tmp"),
            "OMP_NUM_THREADS": "24",
            "MKL_NUM_THREADS": "24",
        }
    )
    tools = {
        "MetaBAT2": assembly_env / "bin/metabat2",
        "SemiBin2": binning_env / "bin/SemiBin2",
        "Vamb": binning_env / "bin/vamb",
        "Kraken2": binning_env / "bin/kraken2",
    }
    for label, path in tools.items():
        if not path.is_file() or not os.access(path, os.X_OK):
            raise FileNotFoundError(f"Missing {label}: {path}")
    versions = {
        "MetaBAT2": parse_version("MetaBAT2", command_text([str(tools["MetaBAT2"]), "-h"], runtime, True)),
        "SemiBin2": parse_version("SemiBin2", command_text([str(tools["SemiBin2"]), "--version"], runtime)),
        "Vamb": parse_version("Vamb", command_text([str(tools["Vamb"]), "--version"], runtime)),
        "Kraken2": parse_version("Kraken2", command_text([str(tools["Kraken2"]), "--version"], runtime)),
    }
    for tool, expected in EXPECTED.items():
        if versions[tool] != expected:
            raise RuntimeError(f"{tool} version drift: expected {expected}, observed {versions[tool]}")
    write_tsv(
        work / "tool-versions.tsv",
        [
            {
                "Tool": tool,
                "Version": version,
                "Executable": (
                    f"${{ASSEMBLY_ENV_PREFIX}}/bin/metabat2"
                    if tool == "MetaBAT2"
                    else f"${{BINNING_ENV_PREFIX}}/bin/{tools[tool].name}"
                ),
            }
            for tool, version in versions.items()
        ],
    )

    required_db = ["hash.k2d", "opts.k2d", "taxo.k2d", "nodes.dmp", "names.dmp"]
    for name in required_db:
        if not (kraken_db / name).is_file():
            raise FileNotFoundError(kraken_db / name)
    database_rows = list(
        csv.DictReader(
            (root / "data/small/17-kraken-database-confidence-frozen/database-audit.tsv").open(encoding="utf-8"),
            delimiter="\t",
        )
    )
    locked = next(row for row in database_rows if row["DatabaseID"] == "standard8")
    observed_bytes = sum(path.stat().st_size for path in kraken_db.iterdir() if path.is_file())
    if observed_bytes != int(locked["InstalledBytes"]):
        raise RuntimeError(f"Kraken database byte drift: {observed_bytes} != {locked['InstalledBytes']}")
    write_tsv(
        work / "database-audit.tsv",
        [
            {
                "Database": "Kraken2 Standard-8",
                "Release": locked["ReleaseID"],
                "InstalledFiles": locked["InstalledFiles"],
                "InstalledBytes": observed_bytes,
                "ArchiveMD5": locked["ArchiveMD5"],
                "ArchiveSHA256": locked["ArchiveSHA256"],
                "InternalMD5Failures": locked["InternalMD5Failures"],
                "Evidence": "Article 17 checksum-locked database audit",
                "Status": "PASS",
            }
        ],
    )

    source_fasta = root / "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz"
    source_depth = root / "data/small/41-read-mapping-depth-frozen/raw/jgi-depth.tsv"
    if not source_fasta.is_file() or not source_depth.is_file():
        raise FileNotFoundError("Article 30/41 frozen input is missing")
    kept_names: list[str] = []
    output_fasta = work / f"inputs/megahit-coassembly.ge{args.minimum_contig}.fna"
    with output_fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in fasta_records(source_fasta):
            if len(sequence) < args.minimum_contig:
                continue
            kept_names.append(name)
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    kept = set(kept_names)
    assembly_summary, _ = fasta_summary(output_fasta)

    with source_depth.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(source_depth)
        depth_rows = [row for row in reader if row["contigName"] in kept]
        depth_header = reader.fieldnames
    if len(depth_rows) != len(kept):
        raise RuntimeError("Filtered FASTA/depth coordinate mismatch")
    depth_multi = work / "inputs/jgi-depth.ge1500.tsv"
    with depth_multi.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=depth_header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(depth_rows)
    mock1_mean = next(name for name in depth_header if name.startswith("MOCK1") and not name.endswith("-var"))
    mock1_var = next(name for name in depth_header if name.startswith("MOCK1") and name.endswith("-var"))
    write_tsv(
        work / "inputs/jgi-depth.MOCK1-only.tsv",
        [
            {
                "contigName": row["contigName"],
                "contigLen": row["contigLen"],
                "totalAvgDepth": row[mock1_mean],
                "MOCK1": row[mock1_mean],
                "MOCK1-var": row[mock1_var],
            }
            for row in depth_rows
        ],
    )
    mock2_mean = next(name for name in depth_header if name.startswith("MOCK2") and not name.endswith("-var"))
    write_tsv(
        work / "inputs/vamb-abundance.tsv",
        [
            {
                "contigname": row["contigName"],
                "MOCK1": row[mock1_mean],
                "MOCK2": row[mock2_mean],
            }
            for row in depth_rows
        ],
    )

    input_audit = [
        {
            "Role": "Article30-coassembly-ge1000",
            "Path": "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz",
            "Bytes": source_fasta.stat().st_size,
            "SHA256": sha256(source_fasta),
            "Status": "PASS",
        },
        {
            "Role": "Article41-JGI-depth",
            "Path": "data/small/41-read-mapping-depth-frozen/raw/jgi-depth.tsv",
            "Bytes": source_depth.stat().st_size,
            "SHA256": sha256(source_depth),
            "Status": "PASS",
        },
        {
            "Role": "Article42-common-coordinate-fasta",
            "Path": "${ARTICLE42_WORK_DIR}/inputs/megahit-coassembly.ge1500.fna",
            "Bytes": output_fasta.stat().st_size,
            "SHA256": sha256(output_fasta),
            "Status": "PASS",
        },
    ]
    write_tsv(work / "input-audit.tsv", input_audit)
    write_tsv(
        work / "input-lineage.tsv",
        [
            {
                "Output": "Common binner FASTA and depth matrix",
                "ImmediateInput": "Article 30 co-assembly + Article 41 JGI depth",
                "Transformation": f"exact identifier join; contig length >= {args.minimum_contig} bp",
                "TruthUsedDuringBinning": "No",
                "TaxonomyUsedDuringBinning": "Kraken2 Standard-8 only for TaxVAMB",
                "Evidence": "input-audit.tsv; run-contract.json",
            },
            {
                "Output": "Post-hoc mock truth audit",
                "ImmediateInput": "Article 33 MetaQUAST coordinates + 87 known MOCK2 genomes",
                "Transformation": "best non-overlapping aligned query span per genome; never supplied to binners",
                "TruthUsedDuringBinning": "No",
                "TaxonomyUsedDuringBinning": "No",
                "Evidence": "truth-contig-assignment.tsv.gz",
            },
        ],
    )
    dump_json(
        work / "run-contract.json",
        {
            "article": 42,
            "seed": 20260742,
            "minimum_contig_bp": args.minimum_contig,
            "minimum_output_bin_bp": 200000,
            "coordinate_set": assembly_summary,
            "branches": [
                "MetaBAT2-MOCK1-only",
                "MetaBAT2-multisample",
                "SemiBin2-self-supervised",
                "VAMB-taxonomy-free",
                "TaxVAMB-Kraken2",
            ],
            "kraken": {
                "release": locked["ReleaseID"],
                "confidence": 0.05,
                "minimum_hit_groups": 2,
            },
            "truth_blinding": True,
            "prepared_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    (work / ".article42-inputs-complete").write_text("PASS\n", encoding="utf-8")
    print(json.dumps({"work": str(work), "coordinate_set": assembly_summary, "versions": versions}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
