#!/usr/bin/env python3
"""Build the local dbCAN and MEROPS indexes required by METABOLIC-G v4.0."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shlex
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--metabolic-dir", type=Path, required=True)
    parser.add_argument("--metabolic-env", type=Path, required=True)
    return parser.parse_args()


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def timed(
    label: str,
    command: list[str],
    cwd: Path,
    logs: Path,
    environment: dict[str, str],
) -> dict[str, object]:
    stdout = logs / f"{label}.stdout.log"
    stderr = logs / f"{label}.stderr.log"
    timing = logs / f"{label}.time.txt"
    wrapped = ["/usr/bin/time", "-v", "-o", str(timing), *command]
    with stdout.open("w", encoding="utf-8") as out_handle, stderr.open(
        "w", encoding="utf-8"
    ) as err_handle:
        completed = subprocess.run(
            wrapped,
            cwd=cwd,
            env=environment,
            stdout=out_handle,
            stderr=err_handle,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed; inspect {stderr}")
    return {
        "Label": label,
        "Command": shlex.join(command),
        "ExitStatus": completed.returncode,
        "Stdout": str(stdout),
        "Stderr": str(stderr),
        "Timing": str(timing),
    }


def normalize_merops(source: Path, destination: Path) -> dict[str, object]:
    """Mirror METABOLIC's whitespace cleanup without mutating the source asset."""
    headers = 0
    sequences = 0
    removed_whitespace = 0
    with source.open("rb") as input_handle, destination.open("wb") as output_handle:
        for raw in input_handle:
            line = raw.rstrip(b"\r\n")
            if line.startswith(b">"):
                headers += 1
                output_handle.write(line + b"\n")
            else:
                cleaned = b"".join(line.split())
                removed_whitespace += len(line) - len(cleaned)
                if cleaned:
                    sequences += 1
                    output_handle.write(cleaned + b"\n")
    if headers == 0 or sequences == 0 or removed_whitespace == 0:
        raise RuntimeError("Unexpected MEROPS normalization result")
    return {
        "Headers": headers,
        "SequenceLines": sequences,
        "WhitespaceBytesRemoved": removed_whitespace,
        "NormalizedBytes": destination.stat().st_size,
        "NormalizedSHA256": digest(destination),
    }


def main() -> None:
    args = parse_args()
    work = args.work_dir.resolve()
    metabolic = args.metabolic_dir.resolve()
    environment = args.metabolic_env.resolve()
    if not (work / ".article59-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article59_metabolism.py first")
    marker = work / ".article59-metabolic-database-complete"
    dbcan = metabolic / "dbCAN2/dbCAN-fam-HMMs.txt"
    dbcan_indexes = [Path(f"{dbcan}.{suffix}") for suffix in ("h3f", "h3i", "h3m", "h3p")]
    merops = metabolic / "MEROPS/pepunit.db.dmnd"
    expected = [*dbcan_indexes, merops]
    if marker.is_file():
        if not all(path.is_file() and path.stat().st_size > 0 for path in expected):
            raise RuntimeError("METABOLIC database marker exists but indexes are missing")
        print("Reusing completed METABOLIC indexes")
        return
    present_dbcan = [path for path in dbcan_indexes if path.exists()]
    if present_dbcan and len(present_dbcan) != len(dbcan_indexes):
        raise RuntimeError(f"Partial dbCAN indexes require inspection: {present_dbcan}")
    if merops.exists():
        raise RuntimeError(f"MEROPS index exists without completion marker: {merops}")

    profile_count = sum(1 for _ in (metabolic / "kofam_database/profiles").glob("K*.hmm"))
    if profile_count != 2643:
        raise RuntimeError(f"Expected 2,643 compatible KOfam profiles, observed {profile_count}")
    env = os.environ.copy()
    env["PATH"] = str(environment / "bin") + os.pathsep + env.get("PATH", "")
    logs = work / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    commands: list[dict[str, object]] = []
    if present_dbcan:
        commands.append(
            {
                "Label": "metabolic-dbcan-hmmpress",
                "Command": "reused complete indexes from prior audited setup attempt",
                "ExitStatus": 0,
                "Stdout": str(logs / "metabolic-dbcan-hmmpress.stdout.log"),
                "Stderr": str(logs / "metabolic-dbcan-hmmpress.stderr.log"),
                "Timing": str(logs / "metabolic-dbcan-hmmpress.time.txt"),
            }
        )
    else:
        commands.append(
            timed(
                "metabolic-dbcan-hmmpress",
                [str(environment / "bin/hmmpress"), "-f", str(dbcan)],
                metabolic,
                logs,
                env,
            )
        )
    normalized_merops = metabolic / "MEROPS/pepunit.normalized.lib"
    normalization = normalize_merops(
        metabolic / "MEROPS/pepunit.lib", normalized_merops
    )
    commands.append(
        timed(
            "metabolic-merops-diamond",
            [
                str(environment / "bin/diamond"),
                "makedb",
                "--in",
                str(normalized_merops),
                "--db",
                str(metabolic / "MEROPS/pepunit.db"),
            ],
            metabolic,
            logs,
            env,
        )
    )
    if not all(path.is_file() and path.stat().st_size > 0 for path in expected):
        raise RuntimeError("METABOLIC index construction did not produce all required files")
    write_tsv(work / "metabolic-database-build.tsv", commands)
    write_tsv(work / "merops-normalization-audit.tsv", [normalization])
    audit = [
        {"File": str(path.relative_to(metabolic)), "Bytes": path.stat().st_size, "SHA256": digest(path)}
        for path in expected
    ]
    write_tsv(work / "metabolic-database-file-audit.tsv", audit)
    marker.write_text("METABOLIC dbCAN and MEROPS indexes completed\n", encoding="utf-8")
    print("Built and audited METABOLIC dbCAN and MEROPS indexes")


if __name__ == "__main__":
    main()
