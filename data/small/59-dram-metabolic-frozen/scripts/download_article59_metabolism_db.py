#!/usr/bin/env python3
"""Download and verify the checksum-locked Article 59 metabolism resources."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
from pathlib import Path


DRAM_COMMIT = "fe61d759303f30db058d5d505c448b28e41b03f1"
METABOLIC_COMMIT = "97236332519180f1d76a242dedb0aaa8191fdbb3"
REPOSITORIES = {
    "DRAM": ("https://github.com/WrightonLabCSU/DRAM.git", DRAM_COMMIT),
    "METABOLIC": ("https://github.com/AnantharamanLab/METABOLIC.git", METABOLIC_COMMIT),
}
EXTERNAL_TARGETS = {
    "profiles.tar.gz": "kofam-2026-06-01/profiles.tar.gz",
    "ko_list.gz": "kofam-2026-06-01/ko_list.gz",
    "dbCAN.hmm": "dbcan-5.2.9/dbCAN.hmm",
    "pepunit.lib": "merops-2023-02-22/pepunit.lib",
}
SOURCE_TARGETS = {
    "METABOLIC_hmm_db.tgz": "sources/METABOLIC/METABOLIC_hmm_db.tgz",
    "METABOLIC_template_and_database.tgz": "sources/METABOLIC/METABOLIC_template_and_database.tgz",
    "Accessory_scripts.tgz": "sources/METABOLIC/Accessory_scripts.tgz",
    "genome_summary_form.tsv": "sources/DRAM/data/genome_summary_form.tsv",
    "module_step_form.tsv": "sources/DRAM/data/module_step_form.tsv",
    "etc_module_database.tsv": "sources/DRAM/data/etc_module_database.tsv",
    "function_heatmap_form.tsv": "sources/DRAM/data/function_heatmap_form.tsv",
    "amg_database.tsv": "sources/DRAM/data/amg_database.tsv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("mode", choices=("download", "verify"))
    return parser.parse_args()


def run(command: list[str], cwd: Path | None = None) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(command)}\n{completed.stdout}")
    return completed.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 12 or len({row["Asset"] for row in rows}) != 12:
        raise RuntimeError("Article 59 database manifest must contain 12 unique assets")
    return rows


def clone_or_verify(cache: Path, name: str, url: str, commit: str, download: bool) -> None:
    target = cache / "sources" / name
    if not target.exists():
        if not download:
            raise FileNotFoundError(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        run(["git", "clone", "--no-checkout", url, str(target)])
    if not (target / ".git").is_dir():
        raise RuntimeError(f"Source directory is not a Git checkout: {target}")
    if download:
        run(["git", "fetch", "--depth", "1", "origin", commit], cwd=target)
        run(["git", "checkout", "--detach", commit], cwd=target)
    observed = run(["git", "rev-parse", "HEAD"], cwd=target)
    if observed != commit:
        raise RuntimeError(f"{name} commit mismatch: expected {commit}, observed {observed}")


def download_external(cache: Path, row: dict[str, str]) -> None:
    target = cache / EXTERNAL_TARGETS[row["Asset"]]
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_file() and target.stat().st_size == int(row["ExpectedBytes"]) and sha256(target) == row["SHA256"]:
        return
    part = target.with_name(target.name + ".part")
    run(
        [
            "curl", "--fail", "--location", "--retry", "50", "--retry-all-errors",
            "--retry-delay", "5", "--connect-timeout", "30", "--speed-time", "180",
            "--speed-limit", "1024", "--continue-at", "-", "--output", str(part), row["Source"],
        ]
    )
    if part.stat().st_size != int(row["ExpectedBytes"]) or sha256(part) != row["SHA256"]:
        raise RuntimeError(f"Downloaded asset failed identity gate: {row['Asset']}")
    part.replace(target)


def extract_resources(cache: Path) -> None:
    kofam = cache / "kofam-2026-06-01"
    profiles = kofam / "profiles"
    if not profiles.is_dir() or len(list(profiles.glob("K*.hmm"))) != 27_754:
        if profiles.exists():
            shutil.rmtree(profiles)
        run(["tar", "-xzf", str(kofam / "profiles.tar.gz"), "-C", str(kofam)])
    metabolic = cache / "sources/METABOLIC"
    for archive in (
        "METABOLIC_hmm_db.tgz",
        "METABOLIC_template_and_database.tgz",
        "Accessory_scripts.tgz",
        "Motif.tgz",
    ):
        run(["tar", "-xzf", str(metabolic / archive), "-C", str(metabolic)])
    for directory in ("kofam_database", "dbCAN2", "MEROPS"):
        (metabolic / directory).mkdir(exist_ok=True)


def asset_path(cache: Path, asset: str) -> Path:
    if asset in EXTERNAL_TARGETS:
        return cache / EXTERNAL_TARGETS[asset]
    return cache / SOURCE_TARGETS[asset]


def verify(cache: Path, rows: list[dict[str, str]]) -> None:
    for name, (_, commit) in REPOSITORIES.items():
        target = cache / "sources" / name
        observed = run(["git", "rev-parse", "HEAD"], cwd=target)
        if observed != commit:
            raise RuntimeError(f"{name} commit mismatch: expected {commit}, observed {observed}")
    for row in rows:
        path = asset_path(cache, row["Asset"])
        if not path.is_file():
            raise FileNotFoundError(path)
        observed_bytes = path.stat().st_size
        observed_sha = sha256(path)
        if observed_bytes != int(row["ExpectedBytes"]) or observed_sha != row["SHA256"]:
            raise RuntimeError(
                f"Identity gate failed for {row['Asset']}: bytes={observed_bytes}, sha256={observed_sha}"
            )
        print(f"PASS\t{row['Tool']}\t{row['Release']}\t{row['Asset']}")
    profiles = cache / "kofam-2026-06-01/profiles"
    count = len(list(profiles.glob("K*.hmm")))
    if count != 27_754:
        raise RuntimeError(f"Expected 27,754 extracted KOfam profiles, observed {count}")
    required_metabolic = (
        "METABOLIC_hmm_db", "METABOLIC_template_and_database", "Accessory_scripts",
        "Motif", "kofam_database", "dbCAN2", "MEROPS",
    )
    missing = [name for name in required_metabolic if not (cache / "sources/METABOLIC" / name).is_dir()]
    if missing:
        raise RuntimeError(f"METABOLIC setup directories are missing: {missing}")
    print(f"PASS\tKOfam extracted profiles\t{count}")
    print("PASS\tArticle 59 metabolism database/source audit\t12 assets")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    cache = args.cache_dir.resolve()
    rows = read_manifest(root / "data/small/59-metabolism-database-manifest.tsv")
    if args.mode == "download":
        cache.mkdir(parents=True, exist_ok=True)
        for name, (url, commit) in REPOSITORIES.items():
            clone_or_verify(cache, name, url, commit, True)
        for row in rows:
            if row["Asset"] in EXTERNAL_TARGETS:
                download_external(cache, row)
        extract_resources(cache)
    verify(cache, rows)


if __name__ == "__main__":
    main()
