#!/usr/bin/env python3
"""Reconstruct checksum-audited Article 42 bins for DAS Tool and Binette."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from article41_44_utils import dump_json, fasta_records, fasta_summary, read_tsv, sha256, write_tsv


EXPECTED_BRANCHES = (
    "MetaBAT2-MOCK1-only", "MetaBAT2-multisample",
    "SemiBin2-self-supervised", "VAMB-taxonomy-free", "TaxVAMB-Kraken2",
)
CHECKM2_DB_SHA256 = "1b86ef3eac0813c1853f53182c17657045e3763d66f384ec95747261a63ae46f"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--article42-frozen", type=Path, required=True)
    parser.add_argument("--assembly-env", type=Path, required=True)
    parser.add_argument("--mag-qc-env", type=Path, required=True)
    parser.add_argument("--checkm2-db", type=Path, required=True)
    return parser.parse_args()


def output(command: list[str], env: dict[str, str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    if result.returncode:
        raise RuntimeError(f"Command failed: {command}\n{result.stdout}\n{result.stderr}")
    return result.stdout + result.stderr


def main() -> int:
    args = parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    frozen = args.article42_frozen.resolve()
    assembly_env, mag_qc_env = args.assembly_env.resolve(), args.mag_qc_env.resolve()
    checkm2_db = args.checkm2_db.resolve()
    if (work / ".article43-inputs-complete").is_file():
        print(f"Article 43 inputs already prepared: {work}")
        return 0
    if work.exists() and any(path.is_file() for path in work.rglob("*")):
        raise SystemExit(f"Refusing incomplete work directory containing files: {work}")
    required_frozen = (
        "bin-membership.tsv.gz", "bin-quality-truth-audit.tsv", "truth-contig-assignment.tsv.gz",
        "file-checksums.sha256", "run-summary.json",
    )
    for name in required_frozen:
        if not (frozen / name).is_file():
            raise FileNotFoundError(frozen / name)
    if not checkm2_db.is_file() or sha256(checkm2_db) != CHECKM2_DB_SHA256:
        raise RuntimeError("CheckM2 version-3 database is missing or has checksum drift")
    for directory in ("inputs", "tables", "source-bins", "refinement", "qc", "logs", "summary", "tmp"):
        (work / directory).mkdir(parents=True, exist_ok=True)

    runtime = os.environ.copy()
    runtime.update(
        {
            "PATH": f"{mag_qc_env / 'bin'}:{assembly_env / 'bin'}:/usr/bin:/bin",
            "LC_ALL": "C", "TZ": "UTC", "PYTHONHASHSEED": "0", "PYTHONNOUSERSITE": "1",
        }
    )
    dastool = assembly_env / "bin/DAS_Tool"
    prodigal = assembly_env / "bin/prodigal"
    binette = mag_qc_env / "bin/binette"
    checkm2 = mag_qc_env / "bin/checkm2"
    versions = {
        "DAS Tool": re.search(r"DAS Tool\s+([0-9.]+)", output([str(dastool), "--version"], runtime)).group(1),
        "Binette": re.search(r"Binette\s+([0-9.]+)", output([str(binette), "--version"], runtime)).group(1),
        "Prodigal": re.search(r"Prodigal V([0-9.]+)", output([str(prodigal), "-v"], runtime)).group(1),
        "CheckM2": output([str(checkm2), "--version"], runtime).strip(),
    }
    expected = {"DAS Tool": "1.1.7", "Binette": "1.2.1", "Prodigal": "2.6.3", "CheckM2": "1.1.0"}
    if versions != expected:
        raise RuntimeError(f"Article 43 version drift: {versions}")
    write_tsv(
        work / "tool-versions.tsv",
        [
            {"Tool": tool, "Version": version, "Executable": str({"DAS Tool": dastool, "Binette": binette, "Prodigal": prodigal, "CheckM2": checkm2}[tool])}
            for tool, version in versions.items()
        ],
    )

    source = root / "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz"
    if not source.is_file():
        raise FileNotFoundError(source)
    common = work / "inputs/megahit-coassembly.ge1500.fna"
    sequences = {}
    with common.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in fasta_records(source):
            if len(sequence) < 1500:
                continue
            sequences[name] = sequence
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    assembly_summary, _ = fasta_summary(common)
    if assembly_summary != {"Contigs": 10203, "TotalBp": 74932939, "MinimumBp": 1500, "MaximumBp": 1064594, "N50Bp": 23167, "GCPct": 48.517488417209954}:
        raise RuntimeError(f"Article 42 coordinate contract drift: {assembly_summary}")

    membership = read_tsv(frozen / "bin-membership.tsv.gz")
    quality = {row["CandidateID"]: row for row in read_tsv(frozen / "bin-quality-truth-audit.tsv")}
    observed_branches = tuple(dict.fromkeys(row["Branch"] for row in membership))
    if set(observed_branches) != set(EXPECTED_BRANCHES):
        raise RuntimeError(f"Unexpected Article 42 branches: {observed_branches}")
    by_branch_bin: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in membership:
        if row["Contig"] not in sequences:
            raise RuntimeError(f"Article 42 membership has an unknown contig: {row['Contig']}")
        by_branch_bin[row["Branch"]][row["CandidateID"]].append(row["Contig"])
    table_rows = []
    reconstruction_rows = []
    for branch in EXPECTED_BRANCHES:
        seen = set()
        branch_dir = work / "source-bins" / re.sub(r"[^A-Za-z0-9]+", "-", branch).strip("-").lower()
        branch_dir.mkdir()
        rows = []
        for candidate_id, names in sorted(by_branch_bin[branch].items()):
            if seen & set(names):
                raise RuntimeError(f"Within-branch overlapping bins: {branch}")
            seen.update(names)
            target = branch_dir / f"{candidate_id}.fna"
            with target.open("w", encoding="utf-8", newline="\n") as handle:
                for name in sorted(names):
                    handle.write(f">{name}\n")
                    sequence = sequences[name]
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80] + "\n")
                    rows.append({"Contig": name, "Bin": candidate_id})
            expected_sha = quality[candidate_id]["CandidateSHA256"]
            observed_sha = sha256(target)
            if observed_sha != expected_sha:
                raise RuntimeError(f"Candidate reconstruction checksum drift: {candidate_id}")
            reconstruction_rows.append(
                {"Branch": branch, "CandidateID": candidate_id, "Contigs": len(names), "Bytes": target.stat().st_size, "SHA256": observed_sha, "Status": "PASS"}
            )
        table_path = work / "tables" / f"{re.sub(r'[^A-Za-z0-9]+', '-', branch).strip('-').lower()}.contig2bin.tsv"
        # DAS Tool and Binette both expect a headerless two-column mapping.
        # A literal ``Contig\tBin`` header would otherwise be interpreted as
        # one more (non-existent) contig assignment.
        with table_path.open("w", encoding="utf-8", newline="\n") as handle:
            for row in sorted(rows, key=lambda item: item["Contig"]):
                handle.write(f"{row['Contig']}\t{row['Bin']}\n")
        table_rows.append({"Branch": branch, "Table": str(table_path), "Bins": len(by_branch_bin[branch]), "BinnedContigs": len(seen), "SHA256": sha256(table_path)})
    write_tsv(work / "candidate-reconstruction-audit.tsv", reconstruction_rows)
    write_tsv(work / "input-binsets.tsv", table_rows)
    write_tsv(
        work / "input-audit.tsv",
        [
            {"Role": "Article42-frozen-manifest", "Path": str((frozen / "file-checksums.sha256").relative_to(root)), "Bytes": (frozen / "file-checksums.sha256").stat().st_size, "SHA256": sha256(frozen / "file-checksums.sha256"), "Status": "PASS"},
            {"Role": "Article30-coassembly-ge1000", "Path": str(source.relative_to(root)), "Bytes": source.stat().st_size, "SHA256": sha256(source), "Status": "PASS"},
            {"Role": "Article43-common-coordinate-fasta", "Path": "${ARTICLE43_WORK_DIR}/inputs/megahit-coassembly.ge1500.fna", "Bytes": common.stat().st_size, "SHA256": sha256(common), "Status": "PASS"},
            {"Role": "CheckM2-database-v3", "Path": str(checkm2_db), "Bytes": checkm2_db.stat().st_size, "SHA256": sha256(checkm2_db), "Status": "PASS"},
        ],
    )
    write_tsv(
        work / "input-lineage.tsv",
        [
            {
                "Output": "Five reconstructed binner partitions",
                "ImmediateInput": "Article 42 checksum-covered bin membership + Article 30 assembly",
                "Transformation": "exact contig-ID reconstruction at the locked >=1500-bp coordinate set",
                "TruthUsedDuringRefinement": "No",
                "Evidence": "candidate-reconstruction-audit.tsv; input-binsets.tsv",
            },
            {
                "Output": "DAS Tool and Binette refined partitions",
                "ImmediateInput": "same five reconstructed binner partitions",
                "Transformation": "DAS Tool SCG scoring or Binette set operations + CheckM2 scoring",
                "TruthUsedDuringRefinement": "No",
                "Evidence": "command-log.tsv; refinement provenance",
            },
        ],
    )
    dump_json(
        work / "run-contract.json",
        {
            "article": 43,
            "seed": 20260743,
            "coordinate_set": assembly_summary,
            "input_branches": list(EXPECTED_BRANCHES),
            "dastool": {"score_threshold": 0.5, "duplicate_penalty": 0.6, "megabin_penalty": 0.5},
            "binette": {"minimum_completeness": 40, "maximum_contamination": 10, "minimum_length_bp": 200000, "contamination_weight": 2.0},
            "final_method_selection": "more CheckM2>=50%, contamination<10%, GUNC-pass bins; tie by summed completeness-5*contamination among passing bins; final tie favors Binette",
            "truth_blinding": True,
            "prepared_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    (work / ".article43-inputs-complete").write_text("PASS\n", encoding="utf-8")
    print(json.dumps({"work": str(work), "coordinate_set": assembly_summary, "input_bins": len(reconstruction_rows)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
