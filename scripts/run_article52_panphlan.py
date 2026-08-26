#!/usr/bin/env python3
"""Run checksum-gated PanPhlAn profiling and a deterministic replay for Article 52."""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path

from article41_44_utils import parse_time, read_tsv, sha256, write_tsv


SOURCE_SHA256 = {
    "panphlan_profiling.py": "901b4dc3710a26145dca064216d62769772b743db55683a956b66253565e7480",
    "misc.py": "d8f7283e847c506178205d4e05da878c5d6a06e1eb767118ca509be4a4af2c15",
}


def run_timed(label: str, command: list[str], work: Path, env: dict[str, str]) -> dict[str, object]:
    stdout_path = work / "logs" / f"{label}.stdout.log"
    stderr_path = work / "logs" / f"{label}.stderr.log"
    time_path = work / "logs" / f"{label}.time.txt"
    timed = ["/usr/bin/time", "-v", "-o", str(time_path), *command]
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        completed = subprocess.run(timed, stdout=stdout, stderr=stderr, env=env, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed; inspect {stderr_path}")
    return {
        "Label": label,
        "ExitStatus": completed.returncode,
        "Command": shlex.join(command),
        "Stdout": str(stdout_path),
        "Stderr": str(stderr_path),
        "TimeLog": str(time_path),
    }


def version(command: list[str], env: dict[str, str]) -> str:
    completed = subprocess.run(command, text=True, capture_output=True, env=env, check=False)
    lines = [
        line.strip()
        for line in (completed.stdout + "\n" + completed.stderr).splitlines()
        if line.strip() and not line.startswith("[WARNING]")
    ]
    return " | ".join(lines[:3]) if lines else "VERSION_NOT_REPORTED"


def matrix_shape(path: Path) -> tuple[int, int, list[str]]:
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = sum(1 for line in handle if line.strip())
    return rows, len(header) - 1, header[1:]


def verify_assets(work: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    archive_locations = {
        "Eubacterium rectale pangenome archive": work / "downloads/Eubacterium_rectale.tar.bz2",
        "25 official tutorial mapping profiles": work / "downloads/panphlan_tutorial_map_results.tar.bz2",
    }
    for item in read_tsv(work / "asset-manifest.tsv"):
        path = archive_locations[item["Asset"]]
        observed_bytes = path.stat().st_size if path.is_file() else -1
        observed_sha = sha256(path) if path.is_file() else "MISSING"
        rows.append({
            "AssetType": "official archive",
            "Asset": item["Asset"],
            "Path": str(path),
            "ExpectedBytes": int(item["Bytes"]),
            "ObservedBytes": observed_bytes,
            "ExpectedSHA256": item["SHA256"],
            "ObservedSHA256": observed_sha,
            "ChecksumPass": observed_bytes == int(item["Bytes"]) and observed_sha == item["SHA256"],
        })

    metadata = read_tsv(work / "sample-metadata.tsv")
    if len(metadata) != 25 or len({row["Sample"] for row in metadata}) != 25:
        raise ValueError("Expected 25 unique official tutorial samples")
    for item in metadata:
        path = work / "decoded-maps" / item["Sample"]
        observed_bytes = path.stat().st_size if path.is_file() else -1
        observed_sha = sha256(path) if path.is_file() else "MISSING"
        rows.append({
            "AssetType": "decoded mapping profile",
            "Asset": item["Sample"],
            "Path": str(path),
            "ExpectedBytes": int(item["DecodedBytes"]),
            "ObservedBytes": observed_bytes,
            "ExpectedSHA256": item["DecodedSHA256"],
            "ObservedSHA256": observed_sha,
            "ChecksumPass": observed_bytes == int(item["DecodedBytes"]) and observed_sha == item["DecodedSHA256"],
        })

    for item in read_tsv(work / "pangenome-file-manifest.tsv"):
        path = work / "pangenome/Eubacterium_rectale" / item["File"]
        observed_bytes = path.stat().st_size if path.is_file() else -1
        observed_sha = sha256(path) if path.is_file() else "MISSING"
        rows.append({
            "AssetType": "pangenome member",
            "Asset": item["File"],
            "Path": str(path),
            "ExpectedBytes": int(item["Bytes"]),
            "ObservedBytes": observed_bytes,
            "ExpectedSHA256": item["SHA256"],
            "ObservedSHA256": observed_sha,
            "ChecksumPass": observed_bytes == int(item["Bytes"]) and observed_sha == item["SHA256"],
        })

    for name, expected in SOURCE_SHA256.items():
        path = work / "software" / name
        observed_sha = sha256(path) if path.is_file() else "MISSING"
        rows.append({
            "AssetType": "pinned source",
            "Asset": name,
            "Path": str(path),
            "ExpectedBytes": path.stat().st_size if path.is_file() else -1,
            "ObservedBytes": path.stat().st_size if path.is_file() else -1,
            "ExpectedSHA256": expected,
            "ObservedSHA256": observed_sha,
            "ChecksumPass": observed_sha == expected,
        })
    write_tsv(work / "asset-check-audit.tsv", rows)
    failed = [row["Asset"] for row in rows if not row["ChecksumPass"]]
    if failed:
        raise ValueError(f"Refusing to execute checksum-failed assets: {failed}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--environment", default="metagenome-panphlan-2026.07")
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article52-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article52_panphlan.py first")
    contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))
    if contract["seed"] != 20260752 or contract["random_output_requested"]:
        raise ValueError("Unexpected deterministic run contract")
    verify_assets(work)

    output = work / "output"
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)
    (work / "logs").mkdir(exist_ok=True)
    for path in (work / "logs").glob("panphlan-*"):
        path.unlink()

    env = os.environ.copy()
    env.update({
        "PYTHONHASHSEED": "0",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    })
    profiler = work / "software/panphlan_profiling.py"
    pangenome = work / "pangenome/Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv"
    prefix = ["conda", "run", "-n", args.environment, "python", str(profiler)]
    primary = [
        *prefix,
        "-i", str(work / "decoded-maps"),
        "-p", str(pangenome),
        "--o_matrix", str(output / "primary-matrix.tsv"),
        "--o_covmat", str(output / "coverage-matrix.tsv"),
        "--o_idx", str(output / "primary-index.tsv"),
        "--min_coverage", "2", "--left_max", "1.25", "--right_min", "0.75",
        "--th_non_present", "0.25", "--th_present", "0.5", "--th_multicopy", "1.5",
        "--add_ref", "-v",
    ]
    sensitivity = [
        *prefix,
        "--i_covmat", str(output / "coverage-matrix.tsv"),
        "-p", str(pangenome),
        "--o_matrix", str(output / "sensitive-matrix.tsv"),
        "--o_idx", str(output / "sensitive-index.tsv"),
        "--min_coverage", "1", "--left_max", "1.70", "--right_min", "0.30",
        "--th_non_present", "0.25", "--th_present", "0.5", "--th_multicopy", "1.5",
        "--add_ref", "-v",
    ]
    replay = list(primary)
    replacements = {
        str(output / "primary-matrix.tsv"): str(output / "replay-matrix.tsv"),
        str(output / "coverage-matrix.tsv"): str(output / "replay-coverage-matrix.tsv"),
        str(output / "primary-index.tsv"): str(output / "replay-index.tsv"),
    }
    replay = [replacements.get(token, token) for token in replay]

    commands = [
        run_timed("panphlan-primary", primary, work, env),
        run_timed("panphlan-sensitivity", sensitivity, work, env),
        run_timed("panphlan-primary-replay", replay, work, env),
    ]
    if any("--o_covplot_normed" in row["Command"] for row in commands):
        raise ValueError("Random native coverage plotting must remain disabled")

    expected_shapes = {
        "primary-matrix.tsv": (11069, 37),
        "primary-index.tsv": (11069, 22),
        "coverage-matrix.tsv": (10703, 25),
        "sensitive-matrix.tsv": (11069, 40),
        "sensitive-index.tsv": (11069, 25),
    }
    shape_rows = []
    for name, expected in expected_shapes.items():
        path = output / name
        rows, columns, names = matrix_shape(path)
        if (rows, columns) != expected:
            raise ValueError(f"Unexpected {name} shape: {(rows, columns)} != {expected}")
        shape_rows.append({
            "File": name, "Rows": rows, "DataColumns": columns,
            "MetagenomeSamples": sum(not value.startswith("REF_") for value in names),
            "ReferenceGenomes": sum(value.startswith("REF_") for value in names),
        })
    write_tsv(work / "output-shape-audit.tsv", shape_rows)

    replay_pairs = (
        ("presence-absence", output / "primary-matrix.tsv", output / "replay-matrix.tsv"),
        ("plateau-index", output / "primary-index.tsv", output / "replay-index.tsv"),
        ("coverage-matrix", output / "coverage-matrix.tsv", output / "replay-coverage-matrix.tsv"),
    )
    deterministic_rows = []
    for label, original, duplicate in replay_pairs:
        original_sha, duplicate_sha = sha256(original), sha256(duplicate)
        deterministic_rows.append({
            "Artifact": label,
            "PrimarySHA256": original_sha,
            "ReplaySHA256": duplicate_sha,
            "ByteIdentical": original_sha == duplicate_sha,
            "Seed": contract["seed"],
            "RandomCoveragePlotRequested": False,
        })
    if not all(row["ByteIdentical"] for row in deterministic_rows):
        raise ValueError("PanPhlAn deterministic replay mismatch")
    write_tsv(work / "determinism-audit.tsv", deterministic_rows)
    write_tsv(work / "command-log.tsv", commands)
    write_tsv(work / "resource-summary.tsv", [parse_time(Path(row["TimeLog"])) for row in commands])

    conda = ["conda", "run", "-n", args.environment]
    source_text = profiler.read_text(encoding="utf-8")
    source_version = re.search(r"__version__\s*=\s*['\"]([^'\"]+)", source_text)
    versions = [
        {"Software": "PanPhlAn profiling source", "Version": source_version.group(1) if source_version else "UNKNOWN"},
        {"Software": "Python", "Version": version([*conda, "python", "-c", "import sys; print(sys.version.split()[0])"], env)},
        {"Software": "NumPy", "Version": version([*conda, "python", "-c", "import numpy; print(numpy.__version__)"], env)},
        {"Software": "Pandas", "Version": version([*conda, "python", "-c", "import pandas; print(pandas.__version__)"], env)},
        {"Software": "Bowtie2", "Version": version([*conda, "bowtie2", "--version"], env)},
        {"Software": "SAMtools", "Version": version([*conda, "samtools", "--version"], env)},
    ]
    write_tsv(work / "tool-versions.tsv", versions)
    write_tsv(work / "output-paths.tsv", [{
        "PrimaryMatrix": str(output / "primary-matrix.tsv"),
        "PrimaryIndex": str(output / "primary-index.tsv"),
        "CoverageMatrix": str(output / "coverage-matrix.tsv"),
        "SensitiveMatrix": str(output / "sensitive-matrix.tsv"),
        "SensitiveIndex": str(output / "sensitive-index.tsv"),
        "Annotation": str(work / "pangenome/Eubacterium_rectale/panphlan_Eubacterium_rectale_annot.tsv"),
    }])
    (work / ".article52-run-complete").write_text("complete\n", encoding="utf-8")
    print("PanPhlAn Article 52 run complete: 22 primary and 25 sensitivity samples")


if __name__ == "__main__":
    main()
