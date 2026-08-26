#!/usr/bin/env python3
"""Build the checksum-gated, KOfam-only DRAM 1.5.0 database for Article 59."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--kofam-dir", type=Path, required=True)
    parser.add_argument("--dram-source", type=Path, required=True)
    parser.add_argument("--dram-env", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=32)
    return parser.parse_args()


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty table")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    work = args.work_dir.resolve()
    kofam = args.kofam_dir.resolve()
    source = args.dram_source.resolve()
    environment = args.dram_env.resolve()
    if not (work / ".article59-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article59_metabolism.py first")

    database = work / "database/dram-kofam-2026-06-01"
    config = work / "dram-config.json"
    logs = work / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    marker = work / ".article59-dram-database-complete"
    expected = [
        database / "kofam_profiles.hmm",
        database / "kofam_profiles.hmm.h3f",
        database / "kofam_profiles.hmm.h3i",
        database / "kofam_profiles.hmm.h3m",
        database / "kofam_profiles.hmm.h3p",
        database / "kofam_ko_list.tsv",
        database / "description_db.sqlite",
        database / "genome_summary_form.tsv",
        database / "module_step_form.tsv",
        database / "etc_module_database.tsv",
        database / "function_heatmap_form.tsv",
    ]
    if marker.is_file():
        if not all(path.is_file() and path.stat().st_size > 0 for path in expected):
            raise RuntimeError("DRAM database marker exists but required files are missing")
        print(f"Reusing completed DRAM database: {database}")
        return
    if database.exists() or config.exists():
        raise RuntimeError(
            "Partial DRAM database state exists without its completion marker; "
            "inspect it before choosing a new work directory"
        )

    installed_config = (
        environment
        / "lib/python3.10/site-packages/mag_annotator/CONFIG"
    )
    if not installed_config.is_file():
        raise FileNotFoundError(installed_config)
    shutil.copy2(installed_config, config)

    forms = source / "data"
    command = [
        str(environment / "bin/DRAM-setup.py"),
        "prepare_databases",
        "--output_dir",
        str(database),
        "--kofam_hmm_loc",
        str(kofam / "profiles.tar.gz"),
        "--kofam_ko_list_loc",
        str(kofam / "ko_list.gz"),
        "--genome_summary_form_loc",
        str(forms / "genome_summary_form.tsv"),
        "--module_step_form_loc",
        str(forms / "module_step_form.tsv"),
        "--etc_module_database_loc",
        str(forms / "etc_module_database.tsv"),
        "--function_heatmap_form_loc",
        str(forms / "function_heatmap_form.tsv"),
        "--amg_database_loc",
        str(forms / "amg_database.tsv"),
        "--threads",
        str(args.threads),
        "--clear_config",
        "--verbose",
    ]
    for database_name in (
        "kofam_hmm",
        "kofam_ko_list",
        "genome_summary_form",
        "module_step_form",
        "etc_module_database",
        "function_heatmap_form",
        "amg_database",
    ):
        command.extend(["--select_db", database_name])

    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"
    env["DRAM_CONFIG_LOCATION"] = str(config)
    stdout = logs / "dram-database.stdout.log"
    stderr = logs / "dram-database.stderr.log"
    timing = logs / "dram-database.time.txt"
    wrapped = ["/usr/bin/time", "-v", "-o", str(timing), *command]
    with stdout.open("w", encoding="utf-8") as out_handle, stderr.open(
        "w", encoding="utf-8"
    ) as err_handle:
        completed = subprocess.run(
            wrapped,
            cwd=work,
            env=env,
            stdout=out_handle,
            stderr=err_handle,
            check=False,
        )
    write_tsv(
        work / "dram-database-build.tsv",
        [
            {
                "Tool": "DRAM-setup.py",
                "Version": "1.5.0",
                "Database": "KOfam 2026-06-01",
                "Command": shlex.join(command),
                "ExitStatus": completed.returncode,
                "Stdout": str(stdout),
                "Stderr": str(stderr),
                "Timing": str(timing),
                "PythonUserSiteDisabled": "true",
            }
        ],
    )
    if completed.returncode != 0:
        raise RuntimeError(f"DRAM database build failed; inspect {stderr}")
    if not all(path.is_file() and path.stat().st_size > 0 for path in expected):
        missing = [str(path) for path in expected if not path.is_file()]
        raise RuntimeError(f"DRAM database build is incomplete: {missing}")

    audit = [
        {
            "File": path.name,
            "Bytes": path.stat().st_size,
            "SHA256": digest(path),
        }
        for path in expected
    ]
    write_tsv(work / "dram-database-file-audit.tsv", audit)
    config_data = json.loads(config.read_text(encoding="utf-8"))
    configured = config_data.get("search_databases", {})
    if not configured.get("kofam_hmm") or not configured.get("kofam_ko_list"):
        raise RuntimeError("DRAM configuration lacks KOfam paths")
    marker.write_text(
        "checksum-gated KOfam-only DRAM database completed\n", encoding="utf-8"
    )
    print(f"Built and audited DRAM database: {database}")


if __name__ == "__main__":
    main()
