#!/usr/bin/env python3
"""Run truth-blinded DAS Tool/Binette refinement and independent QC for Article 43."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from article41_44_utils import fasta_records, sha256, write_tsv


CHECKM2_DB_SHA256 = "1b86ef3eac0813c1853f53182c17657045e3763d66f384ec95747261a63ae46f"
GUNC_DB_SHA256 = "2dabe83f2ab7f0b38e78cfdbd8ca33bdc578b330d7501cb42457d331bc8c09d4"
FASTA_SUFFIXES = {".fa", ".fna", ".fasta"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--assembly-env", type=Path, required=True)
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


def load_log(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    return [
        row
        for row in csv.DictReader(path.open(encoding="utf-8"), delimiter="\t")
        if row.get("ReturnCode") == "0"
    ]


def save_log(path: Path, rows: list[dict[str, object]]) -> None:
    if rows:
        write_tsv(path, rows)


def run_timed(
    label: str,
    command: list[str],
    work: Path,
    runtime: dict[str, str],
    rows: list[dict[str, object]],
) -> None:
    stdout_path = work / "logs" / f"{label}.stdout.log"
    stderr_path = work / "logs" / f"{label}.stderr.log"
    resource_path = work / "logs" / f"{label}.time.txt"
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
    save_log(work / "command-log.tsv", rows)
    if result.returncode:
        tail = stderr_path.read_text(encoding="utf-8", errors="replace")[-8000:]
        raise RuntimeError(f"{label} failed ({result.returncode})\n{tail}")


def fasta_paths(directory: Path) -> list[Path]:
    return sorted(
        path for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES and path.stat().st_size > 0
    )


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    assembly_env = args.assembly_env.resolve()
    mag_qc_env = args.mag_qc_env.resolve()
    gunc_env = args.gunc_env.resolve()
    checkm2_db = args.checkm2_db.resolve()
    gunc_db = args.gunc_db.resolve()
    if not (work / ".article43-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article43_refinement.py first")
    if (work / ".article43-run-complete").is_file():
        print(f"Article 43 refinement already complete: {work}")
        return 0
    observed_checkm2_sha = sha256(checkm2_db) if checkm2_db.is_file() else ""
    observed_gunc_sha = sha256(gunc_db) if gunc_db.is_file() else ""
    if observed_checkm2_sha != CHECKM2_DB_SHA256:
        raise RuntimeError("CheckM2 version-3 database is missing or has checksum drift")
    if observed_gunc_sha != GUNC_DB_SHA256:
        raise RuntimeError("GUNC ProGenomes 2.1 database is missing or has checksum drift")
    for directory in (work / "logs", work / "refinement", work / "qc", work / "tmp"):
        directory.mkdir(parents=True, exist_ok=True)

    common = work / "inputs/megahit-coassembly.ge1500.fna"
    proteins = work / "inputs/megahit-coassembly.ge1500.proteins.faa"
    tables = sorted((work / "tables").glob("*.contig2bin.tsv"))
    if len(tables) != 5:
        raise RuntimeError(f"Expected five input bin tables, observed {len(tables)}")
    command_rows: list[dict[str, object]] = load_log(work / "command-log.tsv")
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

    prodigal_done = work / ".article43-prodigal-complete"
    if not prodigal_done.is_file():
        if proteins.exists():
            proteins.unlink()
        prodigal_command = [
            str(assembly_env / "bin/prodigal"), "-i", str(common), "-a", str(proteins),
            "-p", "meta", "-f", "gff", "-o", str(work / "inputs/megahit-coassembly.ge1500.genes.gff"),
        ]
        run_timed(
            "prodigal-shared-proteins",
            prodigal_command,
            work,
            {**runtime, "PATH": f"{assembly_env / 'bin'}:/usr/bin:/bin"},
            command_rows,
        )
        if not proteins.is_file() or proteins.stat().st_size == 0:
            raise RuntimeError("Prodigal did not create shared proteins")
        prodigal_done.write_text("PASS\n", encoding="utf-8")

    dastool_base = work / "refinement/dastool/article43"
    dastool_bins = work / "refinement/dastool/article43_DASTool_bins"
    dastool_done = work / ".article43-dastool-complete"
    if not dastool_done.is_file():
        if dastool_bins.exists():
            shutil.rmtree(dastool_bins)
        for old in dastool_base.parent.glob(f"{dastool_base.name}_*"):
            if old.is_file():
                old.unlink()
        labels = [path.name.removesuffix(".contig2bin.tsv") for path in tables]
        dastool_command = [
            str(assembly_env / "bin/DAS_Tool"),
            "-i", ",".join(map(str, tables)),
            "-l", ",".join(labels),
            "-c", str(common),
            "-p", str(proteins),
            "-o", str(dastool_base),
            "-t", str(args.threads),
            "--score_threshold", "0.5",
            "--duplicate_penalty", "0.6",
            "--megabin_penalty", "0.5",
            "--write_bin_evals",
            "--write_bins",
        ]
        run_timed(
            "dastool-refinement",
            dastool_command,
            work,
            {**runtime, "PATH": f"{assembly_env / 'bin'}:/usr/bin:/bin"},
            command_rows,
        )
        if not dastool_bins.is_dir() or not fasta_paths(dastool_bins):
            raise RuntimeError("DAS Tool produced no refined FASTA bins")
        dastool_done.write_text("PASS\n", encoding="utf-8")

    binette_dir = work / "refinement/binette"
    binette_bins = binette_dir / "final_bins"
    binette_done = work / ".article43-binette-complete"
    if not binette_done.is_file():
        if binette_dir.exists():
            shutil.rmtree(binette_dir)
        binette_command = [str(mag_qc_env / "bin/binette"), "-b"]
        binette_command.extend(map(str, tables))
        binette_command.extend(
            [
                "-c", str(common), "-p", str(proteins), "-o", str(binette_dir),
                "--prefix", "article43", "--threads", str(args.threads),
                "--min_completeness", "40", "--max_contamination", "10",
                "--min_length", "200000", "--contamination_weight", "2.0",
                "--checkm2_db", str(checkm2_db), "--no-progress",
            ]
        )
        run_timed(
            "binette-refinement",
            binette_command,
            work,
            {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"},
            command_rows,
        )
        if not binette_bins.is_dir() or not fasta_paths(binette_bins):
            raise RuntimeError("Binette produced no refined FASTA bins")
        binette_done.write_text("PASS\n", encoding="utf-8")

    flat = work / "qc/all-refined"
    flat_done = work / ".article43-flat-complete"
    if not flat_done.is_file():
        if flat.exists():
            shutil.rmtree(flat)
        flat.mkdir(parents=True)
        common_sequences = dict(fasta_records(common))
        source_rows = []
        for method, directory in (("DAS Tool", dastool_bins), ("Binette", binette_bins)):
            seen: set[str] = set()
            slug = "dastool" if method == "DAS Tool" else "binette"
            for index, source in enumerate(fasta_paths(directory), start=1):
                records = list(fasta_records(source))
                names = [name for name, _ in records]
                unknown = set(names) - set(common_sequences)
                overlap = seen & set(names)
                if unknown or overlap:
                    raise RuntimeError(f"Invalid {method} partition: unknown={len(unknown)}, overlap={len(overlap)}")
                seen.update(names)
                refined_id = f"{slug}__{index:03d}"
                target = flat / f"{refined_id}.fna"
                with target.open("w", encoding="utf-8", newline="\n") as handle:
                    for name in sorted(names):
                        sequence = common_sequences[name]
                        handle.write(f">{name}\n")
                        for start in range(0, len(sequence), 80):
                            handle.write(sequence[start : start + 80] + "\n")
                source_rows.append(
                    {
                        "Method": method,
                        "RefinedID": refined_id,
                        "SourceBin": source.name,
                        "Contigs": len(names),
                        "FASTA": f"${{ARTICLE43_WORK_DIR}}/qc/all-refined/{target.name}",
                        "SHA256": sha256(target),
                    }
                )
        write_tsv(work / "refined-bin-source.tsv", source_rows)
        flat_done.write_text("PASS\n", encoding="utf-8")

    refined = fasta_paths(flat)
    qc_done = work / ".article43-qc-complete"
    if not qc_done.is_file():
        for partial in (work / "qc/checkm2", work / "qc/gunc", work / "tmp/checkm2", work / "tmp/gunc"):
            if partial.exists():
                shutil.rmtree(partial)
            partial.mkdir(parents=True)
        checkm2_command = [
            str(mag_qc_env / "bin/checkm2"), "predict",
            "--input", str(flat), "--output-directory", str(work / "qc/checkm2"),
            "--extension", ".fna", "--threads", str(args.threads),
            "--database_path", str(checkm2_db), "--tmpdir", str(work / "tmp/checkm2"),
            "--remove_intermediates",
        ]
        run_timed(
            "checkm2-refined-bins",
            checkm2_command,
            work,
            {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"},
            command_rows,
        )
        quality = work / "qc/checkm2/quality_report.tsv"
        if not quality.is_file() or len(list(csv.DictReader(quality.open(encoding="utf-8"), delimiter="\t"))) != len(refined):
            raise RuntimeError("CheckM2 refined-bin row count mismatch")
        gunc_command = [
            str(gunc_env / "bin/gunc"), "run",
            "--input_dir", str(flat), "--file_suffix", ".fna",
            "--db_file", str(gunc_db), "--threads", str(args.threads),
            "--out_dir", str(work / "qc/gunc"), "--temp_dir", str(work / "tmp/gunc"),
            "--contig_taxonomy_output",
        ]
        run_timed(
            "gunc-refined-bins",
            gunc_command,
            work,
            {**runtime, "PATH": f"{gunc_env / 'bin'}:/usr/bin:/bin"},
            command_rows,
        )
        maxcss = sorted((work / "qc/gunc").glob("GUNC.*.maxCSS_level.tsv"))
        if len(maxcss) != 1 or len(list(csv.DictReader(maxcss[0].open(encoding="utf-8"), delimiter="\t"))) != len(refined):
            raise RuntimeError("GUNC refined-bin row count mismatch")
        qc_done.write_text("PASS\n", encoding="utf-8")

    maxcss = sorted((work / "qc/gunc").glob("GUNC.*.maxCSS_level.tsv"))
    write_tsv(
        work / "qc-database-audit.tsv",
        [
            {
                "Tool": "CheckM2", "ToolVersion": "1.1.0", "Database": "DIAMOND database version 3",
                "Reference": "Zenodo 14897628", "Path": str(checkm2_db), "Bytes": checkm2_db.stat().st_size,
                "MD5": md5(checkm2_db), "SHA256": observed_checkm2_sha, "ExpectedSHA256": CHECKM2_DB_SHA256,
                "InstallerIntegrityCheck": "PASS",
            },
            {
                "Tool": "GUNC", "ToolVersion": "1.1.0", "Database": "ProGenomes 2.1",
                "Reference": "GUNC official download_db endpoint", "Path": str(gunc_db), "Bytes": gunc_db.stat().st_size,
                "MD5": md5(gunc_db), "SHA256": observed_gunc_sha, "ExpectedSHA256": GUNC_DB_SHA256,
                "InstallerIntegrityCheck": "PASS (upstream MD5 checked during download)",
            },
        ],
    )
    payload = {
        "article": 43,
        "input_binsets": len(tables),
        "dastool_bins": len(fasta_paths(dastool_bins)),
        "binette_bins": len(fasta_paths(binette_bins)),
        "refined_bins_qc": len(refined),
        "checkm2_quality": str(work / "qc/checkm2/quality_report.tsv"),
        "gunc_maxcss": str(maxcss[0]),
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    (work / ".article43-run-complete").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
