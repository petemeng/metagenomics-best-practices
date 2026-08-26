#!/usr/bin/env python3
"""Freeze compact, checksum-locked Article 33 assembly-QC evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_NAMES = (
    "download_article33_qc_sources.sh",
    "prepare_article33_qc_inputs.py",
    "run_article33_assembly_qc.sh",
    "summarize_article33_assembly_qc.py",
    "freeze_article33_assembly_qc.py",
    "validate_article33_assembly_qc.py",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--benchmark-repo", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_copy(source: Path, target: Path, replacements: dict[str, str]) -> None:
    text = source.read_text(encoding="utf-8", errors="replace")
    for old, new in sorted(replacements.items(), key=lambda item: -len(item[0])):
        text = text.replace(old, new)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    env_prefix = args.env_prefix.resolve()
    repo = args.benchmark_repo.resolve()
    work = args.work_dir.resolve()
    frozen = args.frozen_dir.resolve()
    if frozen.exists() and any(frozen.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty frozen directory: {frozen}")
    frozen.mkdir(parents=True, exist_ok=True)
    for directory in ("env", "logs", "metaquast", "quast", "resources", "scripts", "sources"):
        (frozen / directory).mkdir(parents=True, exist_ok=True)

    required = (
        work / "summary/run-summary.json",
        work / "summary/branch-metrics.tsv",
        work / "summary/per-genome-metaquast.tsv",
        work / "summary/diagnostic-control-effects.tsv",
        work / "quast/.article33-complete",
        work / "metaquast/MOCK1/.article33-complete",
        work / "metaquast/MOCK2/.article33-complete",
        work / "metaquast/MOCK1_MOCK2/.article33-complete",
    )
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    direct = {
        root / "data/small/33-source-manifest.tsv": frozen / "source-manifest.tsv",
        root / "data/small/33-software-releases.tsv": frozen / "software-releases.tsv",
        root / "data/small/33-data-NOTICE.txt": frozen / "data-NOTICE.txt",
        root / "env/hybrid-assembly.yml": frozen / "env/hybrid-assembly.yml",
        root / "env/hybrid-assembly-linux-64.lock": frozen / "env/hybrid-assembly-linux-64.lock",
        repo / "script_r/Supplementary_Table_S1.xlsx": frozen / "sources/Supplementary_Table_S1.xlsx",
        repo / "LICENSE": frozen / "sources/benchmark-repository-LICENSE.txt",
        work / "tool-versions.tsv": frozen / "tool-versions.tsv",
        work / "quast/transposed_report.tsv": frozen / "quast/reference-free-transposed-report.tsv",
    }
    summary_names = (
        "prepare-summary.json",
        "truth-audit.json",
        "input-lineage.tsv",
        "source-bundle-audit.tsv",
        "control-block-audit.tsv",
        "truth-manifest.tsv",
        "branch-metrics.tsv",
        "per-genome-metaquast.tsv",
        "abundance-bin-recovery.tsv",
        "length-threshold-sensitivity.tsv",
        "metric-correlation-audit.tsv",
        "diagnostic-control-effects.tsv",
        "split-scaffold-sensitivity.tsv",
        "resource-usage.tsv",
        "run-summary.json",
    )
    for name in summary_names:
        direct[work / "summary" / name] = frozen / name
    for safe_name in ("MOCK1", "MOCK2", "MOCK1_MOCK2"):
        base = work / "metaquast" / safe_name / "combined_reference"
        direct[base / "transposed_report.tsv"] = frozen / "metaquast" / f"{safe_name}-combined-transposed-report.tsv"
        direct[base / "contigs_reports/transposed_report_misassemblies.tsv"] = frozen / "metaquast" / f"{safe_name}-misassemblies-transposed-report.tsv"

    for source, target in direct.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

    for name in SCRIPT_NAMES:
        source = root / "scripts" / name
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, frozen / "scripts" / name)

    replacements = {
        str(work): "${ARTICLE33_WORK_DIR}",
        str(repo): "${BENCHMARK_REPO}",
        str(root): "${PROJECT_ROOT}",
        str(env_prefix): "${QC_ENV_PREFIX}",
        str(Path.home()): "${HOME}",
    }
    for source in sorted((work / "logs").glob("*.log")):
        normalized_copy(source, frozen / "logs" / source.name, replacements)
    for source in sorted((work / "resources").glob("*.txt")):
        normalized_copy(source, frozen / "resources" / source.name, replacements)

    frozen_contract = {
        "article": 33,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "biological_assemblies": 15,
        "diagnostic_controls": 2,
        "truth_genomes": {"MOCK1": 71, "MOCK2": 87, "MOCK1+MOCK2": 87},
        "per_genome_rows": 1271,
        "large_fasta_policy": "assemblies and per-genome references remain in checksum-locked upstream bundles or Git-ignored source cache",
        "fastq_included": False,
        "assembly_fasta_included": False,
        "routine_validation_requires_network": False,
        "routine_validation_reruns_metaquast": False,
        "upstream_bundle_manifests": [
            "data/small/30-short-read-assembly-frozen/file-checksums.sha256",
            "data/small/31-long-read-assembly-frozen/file-checksums.sha256",
            "data/small/32-hybrid-assembly-polishing-frozen/file-checksums.sha256",
        ],
    }
    (frozen / "frozen-contract.json").write_text(json.dumps(frozen_contract, indent=2) + "\n", encoding="utf-8")

    payloads = sorted(
        path for path in frozen.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    with (frozen / "file-checksums.sha256").open("w", encoding="utf-8", newline="\n") as handle:
        for path in payloads:
            handle.write(f"{sha256(path)}  {path.relative_to(frozen).as_posix()}\n")
    print(json.dumps({"frozen_dir": str(frozen), "payloads": len(payloads)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
