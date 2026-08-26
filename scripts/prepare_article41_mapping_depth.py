#!/usr/bin/env python3
"""Prepare checksum-audited Article 41 mapping/depth inputs."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from article41_44_utils import dump_json, fasta_summary, sha256, write_tsv


EXPECTED = {
    "Bowtie2": "2.5.5",
    "SAMtools": "1.23.1",
    "MetaBAT2 depth": "2.18",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--assembly-env", type=Path, required=True)
    return parser.parse_args()


def run_text(command: list[str], env: dict[str, str]) -> str:
    result = subprocess.run(command, check=True, capture_output=True, text=True, env=env)
    return result.stdout + result.stderr


def parse_version(tool: str, text: str) -> str:
    import re

    patterns = {
        "Bowtie2": r"\bversion\s+(\d+(?:\.\d+){1,3})\b",
        "SAMtools": r"(?m)^samtools\s+(\d+(?:\.\d+){1,3})\b",
        "MetaBAT2 depth": r"jgi_summarize_bam_contig_depths\s+(\d+(?:\.\d+){1,3})\b",
    }
    match = re.search(patterns[tool], text)
    if not match:
        raise ValueError(f"Cannot parse version from: {text[:300]}")
    return match.group(1)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    prefix = args.assembly_env.resolve()
    if (work / ".article41-inputs-complete").is_file():
        print(f"Article 41 inputs already prepared: {work}")
        return 0
    if work.exists():
        unexpected = [path for path in work.rglob("*") if path.is_file()]
        if unexpected:
            raise SystemExit(f"Refusing incomplete work directory containing files: {work}")
    for directory in ("inputs", "index", "bam", "depth", "logs", "summary", "tmp"):
        (work / directory).mkdir(parents=True, exist_ok=True)

    runtime = os.environ.copy()
    runtime.update(
        {
            "PATH": f"{prefix / 'bin'}:/usr/bin:/bin",
            "LC_ALL": "C",
            "TZ": "UTC",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "TMPDIR": str(work / "tmp"),
        }
    )
    executables = {
        "Bowtie2": prefix / "bin/bowtie2",
        "Bowtie2-build": prefix / "bin/bowtie2-build",
        "SAMtools": prefix / "bin/samtools",
        "MetaBAT2 depth": prefix / "bin/jgi_summarize_bam_contig_depths",
    }
    for label, path in executables.items():
        if not path.is_file() or not os.access(path, os.X_OK):
            raise FileNotFoundError(f"Missing {label}: {path}")
    versions = {
        "Bowtie2": parse_version("Bowtie2", run_text([str(executables["Bowtie2"]), "--version"], runtime)),
        "SAMtools": parse_version("SAMtools", run_text([str(executables["SAMtools"]), "--version"], runtime)),
        "MetaBAT2 depth": parse_version("MetaBAT2 depth", run_text([str(executables["MetaBAT2 depth"])], runtime)),
    }
    for tool, expected in EXPECTED.items():
        if versions[tool] != expected:
            raise RuntimeError(f"{tool} version drift: expected {expected}, observed {versions[tool]}")
    write_tsv(
        work / "tool-versions.tsv",
        [
            {
                "Tool": label,
                "Version": versions.get(label, versions.get("Bowtie2", "")),
                "Executable": f"${{ASSEMBLY_ENV_PREFIX}}/bin/{path.name}",
            }
            for label, path in executables.items()
        ],
    )

    source_assembly = root / "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz"
    assembly = work / "inputs/megahit-coassembly.ge1000.fna"
    if not source_assembly.is_file():
        raise FileNotFoundError(source_assembly)
    with gzip.open(source_assembly, "rb") as source, assembly.open("wb") as target:
        shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
    assembly_stats, _ = fasta_summary(assembly)
    if assembly_stats["Contigs"] != 18354 or assembly_stats["MinimumBp"] < 1000:
        raise RuntimeError(f"Unexpected Article 30 assembly contract: {assembly_stats}")

    sample_rows = []
    input_audit = [
        {
            "Role": "coassembly-gzip",
            "Path": "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz",
            "Bytes": source_assembly.stat().st_size,
            "SHA256": sha256(source_assembly),
            "Status": "PASS",
        },
        {
            "Role": "coassembly-fasta",
            "Path": "${ARTICLE41_WORK_DIR}/inputs/megahit-coassembly.ge1000.fna",
            "Bytes": assembly.stat().st_size,
            "SHA256": sha256(assembly),
            "Status": "PASS",
        },
    ]
    for sample, run in (("MOCK1", "ERR9765746"), ("MOCK2", "ERR9765747")):
        fastp = root / f"data/raw/article30/work/fastp/{sample}.json"
        payload = json.loads(fastp.read_text(encoding="utf-8"))
        after = payload["summary"]["after_filtering"]
        pairs = int(after["total_reads"]) // 2
        paths = []
        for mate in ("R1", "R2"):
            path = root / f"data/raw/article30/clean/{run}_clean_{mate}.fastq.gz"
            if not path.is_file():
                raise FileNotFoundError(path)
            with gzip.open(path, "rb") as handle:
                while handle.read(8 * 1024 * 1024):
                    pass
            paths.append(path)
            input_audit.append(
                {
                    "Role": f"{sample}-{mate}",
                    "Path": f"data/raw/article30/clean/{run}_clean_{mate}.fastq.gz",
                    "Bytes": path.stat().st_size,
                    "SHA256": sha256(path),
                    "Status": "PASS",
                }
            )
        sample_rows.append(
            {
                "Sample": sample,
                "RunAccession": run,
                "ReadPairs": pairs,
                "R1": str(paths[0]),
                "R2": str(paths[1]),
                "ReadLengthMean": after["read1_mean_length"],
                "SourceArticle": 30,
            }
        )
    write_tsv(work / "input-audit.tsv", input_audit)
    write_tsv(work / "inputs/samples.tsv", sample_rows)
    write_tsv(
        work / "input-lineage.tsv",
        [
            {
                "Output": "megahit-coassembly.ge1000.fna",
                "ImmediateInput": "Article 30 MEGAHIT co-assembly",
                "Transformation": "deterministic length filter >=1000 bp during Article 30 freezing; gzip decompression only here",
                "SourceRun": "ERR9765746 + ERR9765747",
                "SourceStudy": "PRJEB62467",
                "Evidence": "Article 30 checksum and Article 41 input-audit.tsv",
            },
            {
                "Output": "MOCK1/2 coordinate-sorted BAM",
                "ImmediateInput": "Article 30 fastp-clean paired FASTQ + common co-assembly",
                "Transformation": "Bowtie2 very-sensitive; primary mapped records; coordinate sort",
                "SourceRun": "ERR9765746 / ERR9765747",
                "SourceStudy": "PRJEB62467",
                "Evidence": "Article 41 command-log.tsv and mapping-summary.tsv",
            },
        ],
    )
    dump_json(
        work / "run-contract.json",
        {
            "article": 41,
            "seed": 20260741,
            "coordinate_system": "Article 30 MEGAHIT co-assembly contigs >=1000 bp",
            "alignment": {
                "preset": "--very-sensitive",
                "seed": 20260741,
                "bam_exclude_flag": 3588,
                "minimum_mapq": 0,
            },
            "jgi_depth": {
                "minimum_end_to_end_identity_percent": 97,
                "minimum_mapq": 0,
                "edge_bases_excluded": True,
            },
            "assembly": assembly_stats,
            "prepared_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    (work / ".article41-inputs-complete").write_text("PASS\n", encoding="utf-8")
    print(json.dumps({"work": str(work), "assembly": assembly_stats, "samples": sample_rows}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
