#!/usr/bin/env python3
"""Run CheckM2 and GUNC on all standardized Article 42 candidate bins."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from article41_44_utils import sha256, write_tsv


CHECKM2_DB_SHA256 = "1b86ef3eac0813c1853f53182c17657045e3763d66f384ec95747261a63ae46f"
GUNC_DB_SHA256 = "2dabe83f2ab7f0b38e78cfdbd8ca33bdc578b330d7501cb42457d331bc8c09d4"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--mag-qc-env", type=Path, required=True)
    parser.add_argument("--gunc-env", type=Path, required=True)
    parser.add_argument("--checkm2-db", type=Path, required=True)
    parser.add_argument("--gunc-db", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=24)
    return parser.parse_args()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_timed(label: str, command: list[str], work: Path, runtime: dict[str, str], rows: list[dict[str, object]]) -> None:
    log_dir = work / "logs"
    stdout_path = log_dir / f"{label}.stdout.log"
    stderr_path = log_dir / f"{label}.stderr.log"
    resource_path = log_dir / f"{label}.time.txt"
    started = datetime.now(timezone.utc)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        result = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(resource_path), *command],
            stdout=stdout,
            stderr=stderr,
            env=runtime,
            text=True,
        )
    ended = datetime.now(timezone.utc)
    rows.append(
        {
            "Label": label,
            "Command": shlex.join(command),
            "ReturnCode": result.returncode,
            "StartedUTC": started.isoformat(),
            "EndedUTC": ended.isoformat(),
            "Stdout": str(stdout_path),
            "Stderr": str(stderr_path),
            "ResourceLog": str(resource_path),
        }
    )
    write_tsv(work / "qc-command-log.tsv", rows)
    if result.returncode:
        tail = stderr_path.read_text(encoding="utf-8", errors="replace")[-8000:]
        raise RuntimeError(f"{label} failed ({result.returncode})\n{tail}")


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    mag_qc = args.mag_qc_env.resolve()
    gunc_env = args.gunc_env.resolve()
    checkm_db = args.checkm2_db.resolve()
    gunc_db = args.gunc_db.resolve()
    if not (work / ".article42-summary-complete").is_file():
        raise FileNotFoundError("Run summarize_article42_binning.py once to standardize candidate bins")
    if (work / ".article42-qc-complete").is_file():
        print(f"Article 42 candidate QC already complete: {work}")
        return 0
    if not checkm_db.is_file() or not gunc_db.is_file():
        raise FileNotFoundError(f"QC database missing: CheckM2={checkm_db}; GUNC={gunc_db}")
    observed_checkm_sha = sha256(checkm_db)
    if observed_checkm_sha != CHECKM2_DB_SHA256:
        raise RuntimeError(f"CheckM2 v3 database SHA-256 drift: {observed_checkm_sha}")
    observed_gunc_sha = sha256(gunc_db)
    if observed_gunc_sha != GUNC_DB_SHA256:
        raise RuntimeError(f"GUNC ProGenomes 2.1 database SHA-256 drift: {observed_gunc_sha}")
    qc = work / "qc"
    flat = qc / "all-candidates"
    logs = work / "logs"
    for directory in (qc, flat, logs, work / "tmp/checkm2", work / "tmp/gunc"):
        directory.mkdir(parents=True, exist_ok=True)
    candidates = sorted((work / "candidates").glob("*/*.fna"))
    if not candidates:
        raise RuntimeError("No standardized Article 42 candidate FASTAs")
    if len({path.name for path in candidates}) != len(candidates):
        raise RuntimeError("Candidate FASTA basenames are not globally unique")
    for source in candidates:
        target = flat / source.name
        if target.exists():
            if sha256(target) != sha256(source):
                raise RuntimeError(f"Existing flat candidate differs: {target}")
            continue
        try:
            os.link(source, target)
        except OSError:
            shutil.copy2(source, target)
    write_tsv(
        work / "qc-database-audit.tsv",
        [
            {
                "Tool": "CheckM2",
                "ToolVersion": "1.1.0",
                "Database": "CheckM2 DIAMOND database version 3",
                "Reference": "Zenodo 14897628",
                "Path": str(checkm_db),
                "Bytes": checkm_db.stat().st_size,
                "MD5": md5(checkm_db),
                "SHA256": observed_checkm_sha,
                "ExpectedSHA256": CHECKM2_DB_SHA256,
                "InstallerIntegrityCheck": "PASS",
            },
            {
                "Tool": "GUNC",
                "ToolVersion": "1.1.0",
                "Database": "ProGenomes 2.1",
                "Reference": "GUNC official download_db endpoint",
                "Path": str(gunc_db),
                "Bytes": gunc_db.stat().st_size,
                "MD5": md5(gunc_db),
                "SHA256": observed_gunc_sha,
                "ExpectedSHA256": GUNC_DB_SHA256,
                "InstallerIntegrityCheck": "PASS (upstream MD5 checked during download)",
            },
        ],
    )
    runtime = os.environ.copy()
    runtime.update(
        {
            "LC_ALL": "C",
            "TZ": "UTC",
            "PYTHONHASHSEED": "0",
            "PYTHONNOUSERSITE": "1",
            "OMP_NUM_THREADS": str(args.threads),
            "OPENBLAS_NUM_THREADS": str(args.threads),
            "MKL_NUM_THREADS": str(args.threads),
        }
    )
    rows: list[dict[str, object]] = []
    checkm_output = qc / "checkm2"
    checkm_command = [
        str(mag_qc / "bin/checkm2"), "predict",
        "--input", str(flat),
        "--output-directory", str(checkm_output),
        "--extension", ".fna",
        "--threads", str(args.threads),
        "--database_path", str(checkm_db),
        "--tmpdir", str(work / "tmp/checkm2"),
        "--remove_intermediates",
    ]
    quality = checkm_output / "quality_report.tsv"
    checkm_done = work / ".article42-checkm2-complete"
    if not checkm_done.is_file():
        for partial in (checkm_output, work / "tmp/checkm2"):
            if partial.exists():
                shutil.rmtree(partial)
            partial.mkdir(parents=True)
        run_timed("checkm2-all-candidates", checkm_command, work, {**runtime, "PATH": f"{mag_qc / 'bin'}:/usr/bin:/bin"}, rows)
        if not quality.is_file():
            raise FileNotFoundError(quality)
        checkm_done.write_text("PASS\n", encoding="utf-8")
    checkm_rows = list(csv.DictReader(quality.open(encoding="utf-8"), delimiter="\t"))
    if len(checkm_rows) != len(candidates):
        raise RuntimeError(f"CheckM2 candidate count mismatch: {len(checkm_rows)} != {len(candidates)}")
    if not rows:
        rows.append(
            {
                "Label": "checkm2-all-candidates",
                "Command": shlex.join(checkm_command),
                "ReturnCode": 0,
                "StartedUTC": "reused",
                "EndedUTC": "reused",
                "Stdout": str(work / "logs/checkm2-all-candidates.stdout.log"),
                "Stderr": str(work / "logs/checkm2-all-candidates.stderr.log"),
                "ResourceLog": str(work / "logs/checkm2-all-candidates.time.txt"),
            }
        )

    gunc_output = qc / "gunc"
    gunc_command = [
        str(gunc_env / "bin/gunc"), "run",
        "--input_dir", str(flat),
        "--file_suffix", ".fna",
        "--db_file", str(gunc_db),
        "--threads", str(args.threads),
        "--out_dir", str(gunc_output),
        "--temp_dir", str(work / "tmp/gunc"),
        "--contig_taxonomy_output",
    ]
    maxcss_candidates = sorted(gunc_output.glob("GUNC.*.maxCSS_level.tsv"))
    gunc_done = work / ".article42-gunc-complete"
    if not gunc_done.is_file():
        for partial in (gunc_output, work / "tmp/gunc"):
            if partial.exists():
                shutil.rmtree(partial)
            partial.mkdir(parents=True)
        run_timed("gunc-all-candidates", gunc_command, work, {**runtime, "PATH": f"{gunc_env / 'bin'}:/usr/bin:/bin"}, rows)
        maxcss_candidates = sorted(gunc_output.glob("GUNC.*.maxCSS_level.tsv"))
        if len(maxcss_candidates) != 1:
            raise RuntimeError(f"Expected one GUNC maxCSS table: {maxcss_candidates}")
        gunc_done.write_text("PASS\n", encoding="utf-8")
    if len(maxcss_candidates) != 1:
        raise RuntimeError(f"Expected one GUNC maxCSS table: {maxcss_candidates}")
    gunc_rows = list(csv.DictReader(maxcss_candidates[0].open(encoding="utf-8"), delimiter="\t"))
    if len(gunc_rows) != len(candidates):
        raise RuntimeError(f"GUNC candidate count mismatch: {len(gunc_rows)} != {len(candidates)}")
    if not any(row["Label"] == "gunc-all-candidates" for row in rows):
        rows.append(
            {
                "Label": "gunc-all-candidates",
                "Command": shlex.join(gunc_command),
                "ReturnCode": 0,
                "StartedUTC": "reused",
                "EndedUTC": "reused",
                "Stdout": str(work / "logs/gunc-all-candidates.stdout.log"),
                "Stderr": str(work / "logs/gunc-all-candidates.stderr.log"),
                "ResourceLog": str(work / "logs/gunc-all-candidates.time.txt"),
            }
        )
    write_tsv(work / "qc-command-log.tsv", rows)
    payload = {
        "article": 42,
        "candidate_bins": len(candidates),
        "checkm2_rows": len(checkm_rows),
        "gunc_rows": len(gunc_rows),
        "checkm2_quality": str(quality),
        "gunc_maxcss": str(maxcss_candidates[0]),
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    (work / ".article42-qc-complete").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
