#!/usr/bin/env python3
"""Run checksum-gated gapseq and CarveMe reconstructions for Article 60."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--gapseq-env", type=Path, required=True)
    parser.add_argument("--carveme-env", type=Path, required=True)
    parser.add_argument("--seqdb", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--threads-per-job", type=int, default=4)
    parser.add_argument("--only", nargs="*", default=None)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def deterministic_environment(prefix: Path, home: Path) -> dict[str, str]:
    home.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{prefix / 'bin'}:/usr/bin:/bin",
            "HOME": str(home),
            "TMPDIR": "/tmp",
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    return env


def run_step(
    *,
    genome: str,
    tool: str,
    step: str,
    command: list[str],
    required: list[Path],
    log_dir: Path,
    env: dict[str, str],
    lock: threading.Lock,
    ledger: list[dict[str, object]],
) -> None:
    if required and all(path.is_file() and path.stat().st_size > 0 for path in required):
        with lock:
            ledger.append(
                {
                    "Genome": genome,
                    "Tool": tool,
                    "Step": step,
                    "Status": "skipped_nonempty_existing",
                    "ElapsedSeconds": "0.000",
                    "ReturnCode": 0,
                    "Command": shlex.join(command),
                }
            )
        return
    for path in required:
        path.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    prefix = log_dir / f"{tool}-{step}-{genome}"
    stdout = prefix.with_suffix(".stdout.log")
    stderr = prefix.with_suffix(".stderr.log")
    usage = prefix.with_suffix(".time.txt")
    wrapped = ["/usr/bin/time", "-v", "-o", str(usage), *command]
    started = time.perf_counter()
    with stdout.open("w", encoding="utf-8") as out, stderr.open(
        "w", encoding="utf-8"
    ) as err:
        completed = subprocess.run(wrapped, env=env, stdout=out, stderr=err, check=False)
    elapsed = time.perf_counter() - started
    with lock:
        ledger.append(
            {
                "Genome": genome,
                "Tool": tool,
                "Step": step,
                "Status": "completed" if completed.returncode == 0 else "failed",
                "ElapsedSeconds": f"{elapsed:.3f}",
                "ReturnCode": completed.returncode,
                "Command": shlex.join(command),
            }
        )
    if completed.returncode != 0:
        tail = "\n".join(stderr.read_text(encoding="utf-8").splitlines()[-30:])
        raise RuntimeError(
            f"{tool} {step} failed for {genome} ({completed.returncode})\n{tail}"
        )
    missing = [path for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"{tool} {step} produced missing/empty outputs for {genome}: {missing}")


def run_gapseq(
    args: argparse.Namespace,
    row: dict[str, str],
    env: dict[str, str],
    lock: threading.Lock,
    ledger: list[dict[str, object]],
) -> None:
    genome = row["Genome"]
    gapseq = args.gapseq_env / "bin/gapseq"
    protein = args.work_dir / "inputs/proteins" / f"{genome}.faa"
    output = args.work_dir / "gapseq" / genome
    output.mkdir(parents=True, exist_ok=True)
    reactions = output / f"{genome}-all-Reactions.tbl"
    pathways = output / f"{genome}-all-Pathways.tbl"
    transport = output / f"{genome}-Transporter.tbl"
    draft_rds = output / f"{genome}-draft.RDS"
    draft_xml = output / f"{genome}-draft.xml"
    taxonomy = row["Domain"]
    biomass = "Archaea" if taxonomy == "Archaea" else "Bacteria"
    common = {
        "genome": genome,
        "log_dir": args.work_dir / "logs/gapseq",
        "env": env,
        "lock": lock,
        "ledger": ledger,
    }
    find_command = [
        str(gapseq), "find", "-p", "all", "-t", taxonomy,
    ]
    if taxonomy == "Archaea":
        # The official archaeal tutorial fixes both the sequence domain (-t)
        # and the pathway taxonomic range (-m). Leaving -m at auto can omit
        # archaeal pathways even when the reference sequence search is correct.
        find_command.extend(["-m", "Archaea"])
    find_command.extend(
        [
            "-b", "200", "-i", "0", "-c", "75", "-K",
            str(args.threads_per_job), "-A", "diamond", "-O", "-D",
            str(args.seqdb), "-v", "1", "-f", str(output), str(protein),
        ]
    )
    run_step(
        **common,
        tool="gapseq",
        step="find",
        command=find_command,
        required=[reactions, pathways],
    )
    run_step(
        **common,
        tool="gapseq",
        step="transport",
        command=[
            str(gapseq), "find-transport", "-b", "50", "-i", "0", "-c", "75",
            "-K", str(args.threads_per_job), "-A", "diamond", "-v", "1",
            "-f", str(output), str(protein),
        ],
        required=[transport],
    )
    run_step(
        **common,
        tool="gapseq",
        step="draft",
        command=[
            str(gapseq), "draft", "-r", str(reactions), "-t", str(transport),
            "-p", str(pathways), "-b", biomass, "-n", genome, "-u", "200",
            "-l", "100", "-f", str(output),
        ],
        required=[draft_rds, draft_xml],
    )
    rich_dir = output / "filled-permissive"
    rich_rds = rich_dir / f"{genome}.RDS"
    rich_xml = rich_dir / f"{genome}.xml"
    profile_by_genome = {
        # The three archaeal stress cases require ecological media from the
        # official gapseq distribution; a bacterial rich proxy fails even
        # when the entire archaeal reaction universe is offered.
        "SGB_008": ("autotrophic.csv", "highH2"),
        "SGB_010": ("MM_anaerobic_CO2_H2.csv", "highH2"),
        "SGB_018": ("meerwasser.csv", None),
    }
    medium_name, environmental_condition = profile_by_genome.get(
        genome, ("ALLmed.csv", None)
    )
    medium = args.gapseq_env / "share/gapseq/dat/media" / medium_name
    fill_command = [
        str(gapseq), "fill", "-m", str(draft_rds), "-n", str(medium),
        "-b", "100", "-k", "0.01", "-z", "glpk", "-v",
    ]
    if environmental_condition:
        fill_command.extend(["-e", environmental_condition])
    fill_command.extend(["-f", str(rich_dir)])
    run_step(
        **common,
        tool="gapseq",
        step="gapfill-construction-medium",
        command=fill_command,
        required=[rich_rds, rich_xml],
    )


def run_carveme(
    args: argparse.Namespace,
    row: dict[str, str],
    env: dict[str, str],
    package_root: Path,
    lock: threading.Lock,
    ledger: list[dict[str, object]],
) -> None:
    genome = row["Genome"]
    carve = args.carveme_env / "bin/carve"
    gapfill = args.carveme_env / "bin/gapfill"
    diamond = args.carveme_env / "bin/diamond"
    protein = args.work_dir / "inputs/proteins" / f"{genome}.faa"
    output = args.work_dir / "carveme" / genome
    output.mkdir(parents=True, exist_ok=True)
    annotations = output / f"{genome}-diamond.tsv"
    universe = "archaea" if row["Domain"] == "Archaea" else "bacteria"
    database = package_root / "data/generated/bigg_proteins.dmnd"
    draft = output / f"{genome}-draft.xml"
    rich = output / f"{genome}-filled-LB.xml"
    common = {
        "genome": genome,
        "log_dir": args.work_dir / "logs/carveme",
        "env": env,
        "lock": lock,
        "ledger": ledger,
    }
    run_step(
        **common,
        tool="CarveMe",
        step="diamond",
        command=[
            str(diamond), "blastp", "-d", str(database), "-q", str(protein),
            "-o", str(annotations), "--more-sensitive", "--top", "10", "--quiet",
            "--threads", str(args.threads_per_job),
        ],
        required=[annotations],
    )
    run_step(
        **common,
        tool="CarveMe",
        step="draft",
        command=[
            str(carve), "--diamond", "-u", universe, "--fbc2", "--solver", "scip",
            "-v", "-o", str(draft), str(annotations),
        ],
        required=[draft],
    )
    run_step(
        **common,
        tool="CarveMe",
        step="gapfill-rich",
        command=[
            str(gapfill), "-m", "LB", "-u", universe, "--fbc2", "-v",
            "-o", str(rich), str(draft),
        ],
        required=[rich],
    )


def package_root(environment: Path) -> Path:
    paths = sorted(
        {path.resolve() for path in (environment / "lib").glob("python*/site-packages/carveme")}
    )
    if len(paths) != 1:
        raise RuntimeError(f"Expected one CarveMe package path, observed {paths}")
    return paths[0]


def main() -> None:
    args = parse_args()
    for attribute in ("project_root", "work_dir", "gapseq_env", "carveme_env", "seqdb"):
        setattr(args, attribute, getattr(args, attribute).resolve())
    if args.jobs < 1 or args.threads_per_job < 1:
        raise ValueError("--jobs and --threads-per-job must both be positive")
    marker = args.work_dir / ".article60-inputs-complete"
    if not marker.is_file():
        raise FileNotFoundError("Run prepare_article60_gem.py first")
    for domain in ("Bacteria", "Archaea"):
        version = args.seqdb / domain / "version_seqDB.json"
        metadata = json.loads(version.read_text(encoding="utf-8"))
        if metadata.get("zenodoID") != 20446806 or metadata.get("version") != "1.5":
            raise RuntimeError(f"Wrong gapseq sequence database identity: {version}")

    rows = read_tsv(args.work_dir / "input-mag-ledger.tsv")
    if args.only:
        requested = set(args.only)
        rows = [row for row in rows if row["Genome"] in requested]
        missing = requested - {row["Genome"] for row in rows}
        if missing:
            raise RuntimeError(f"Unknown --only genome(s): {sorted(missing)}")
    if not rows:
        raise RuntimeError("No genomes selected")
    gapseq_env = deterministic_environment(
        args.gapseq_env, args.work_dir / "home/gapseq"
    )
    carveme_env = deterministic_environment(
        args.carveme_env, args.work_dir / "home/carveme"
    )
    carveme_root = package_root(args.carveme_env)
    lock = threading.Lock()
    ledger: list[dict[str, object]] = []

    futures = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        for row in rows:
            future = executor.submit(
                run_gapseq, args, row, gapseq_env, lock, ledger
            )
            futures[future] = (row["Genome"], "gapseq")
            future = executor.submit(
                run_carveme, args, row, carveme_env, carveme_root, lock, ledger
            )
            futures[future] = (row["Genome"], "CarveMe")
        failures: list[str] = []
        for future in as_completed(futures):
            genome, tool = futures[future]
            try:
                future.result()
                print(f"PASS\t{tool}\t{genome}", flush=True)
            except Exception as error:  # retain every completed command in the ledger
                failures.append(f"{tool} {genome}: {error}")
                print(f"FAIL\t{tool}\t{genome}\t{error}", flush=True)

    ordered = sorted(ledger, key=lambda row: (str(row["Genome"]), str(row["Tool"]), str(row["Step"])))
    write_tsv(args.work_dir / "command-log.tsv", ordered)
    gapseq_profiles = {
        "SGB_008": ("autotrophic.csv", "highH2", "autotrophic hyperthermophile proxy"),
        "SGB_010": ("MM_anaerobic_CO2_H2.csv", "highH2", "official methanogen tutorial medium"),
        "SGB_018": ("meerwasser.csv", "none", "marine host-rich proxy for an obligate symbiont"),
    }
    plan_rows = []
    for row in rows:
        medium, condition, rationale = gapseq_profiles.get(
            row["Genome"],
            ("ALLmed.csv", "none", "common permissive bacterial reconstruction proxy"),
        )
        plan_rows.append({
            "Genome": row["Genome"],
            "Species": row["Species"],
            "Domain": row["Domain"],
            "gapseqTaxonomy": row["Domain"],
            "gapseqPathwayTaxonomicRange": "Archaea" if row["Domain"] == "Archaea" else "auto",
            "gapseqGapfillMedium": medium,
            "gapseqEnvironmentalCondition": condition,
            "gapseqMediumRationale": rationale,
            "gapseqAuditMedium": "ALLmed.csv and MM_glu.csv",
            "CarveMeUniverse": "archaea" if row["Domain"] == "Archaea" else "bacteria",
            "CarveMeGapfillMedium": "LB",
            "CarveMeAuditMedium": "M9",
            "MediumPairingRule": "construction versus common-rich/minimal sensitivity; not chemically identical across tools or ecologies",
        })
    write_tsv(args.work_dir / "model-plan.tsv", plan_rows)
    contract = {
        "article": 60,
        "genomes": len(rows),
        "jobs": args.jobs,
        "threads_per_job": args.threads_per_job,
        "gapseq_solver": "GLPK 5.0",
        "carveme_solver": "SCIP 10.0.3",
        "carveme_scip_time_limit_seconds": 600,
        "carveme_scip_relative_gap": 0.001,
        "carveme_ensemble_size": 1,
        "carveme_gapfill_mode": "independent gapfill CLI applied to the fixed draft model",
        "python_hash_seed": 0,
        "gapseq_sequence_db": "v1.5 / Zenodo 20446806",
        "gapseq_bacterial_gapfill_medium": "ALLmed.csv",
        "gapseq_archaeal_gapfill_profiles": {
            "SGB_008": "autotrophic.csv + highH2",
            "SGB_010": "MM_anaerobic_CO2_H2.csv + highH2",
            "SGB_018": "meerwasser.csv",
        },
        "failures": failures,
    }
    (args.work_dir / "run-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if failures:
        raise RuntimeError("Article 60 model run failed:\n" + "\n".join(failures))
    (args.work_dir / ".article60-models-complete").write_text(
        "Article 60 gapseq and CarveMe models complete\n", encoding="utf-8"
    )
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
