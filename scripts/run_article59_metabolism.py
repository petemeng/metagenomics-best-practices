#!/usr/bin/env python3
"""Run DRAM 1.5.0 and METABOLIC-G v4.0 on the Article 59 MAG panel."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path


EXPECTED_GENOMES = 28


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--metabolic-dir", type=Path, required=True)
    parser.add_argument("--dram-env", type=Path, required=True)
    parser.add_argument("--metabolic-env", type=Path, required=True)
    parser.add_argument("--dram-jobs", type=int, default=24)
    parser.add_argument("--dram-threads-per-job", type=int, default=3)
    parser.add_argument("--metabolic-threads", type=int, default=32)
    return parser.parse_args()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty command table")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def timed_run(
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
    row = {
        "Label": label,
        "Command": shlex.join(command),
        "ExitStatus": completed.returncode,
        "Stdout": str(stdout),
        "Stderr": str(stderr),
        "Timing": str(timing),
    }
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed; inspect {stderr}")
    return row


def existing_run_row(
    label: str, command: list[str], logs: Path
) -> dict[str, object]:
    """Reconstruct a command row after a completed branch outlives a sibling failure."""
    stdout = logs / f"{label}.stdout.log"
    stderr = logs / f"{label}.stderr.log"
    timing = logs / f"{label}.time.txt"
    required = (stdout, stderr, timing)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Cannot resume {label}; run evidence is missing: {missing}")
    return {
        "Label": label,
        "Command": shlex.join(command),
        "ExitStatus": 0,
        "Stdout": str(stdout),
        "Stderr": str(stderr),
        "Timing": str(timing),
    }


def dram_branch_is_complete(work: Path) -> bool:
    required = (
        work / "dram-annotation/annotations.tsv",
        work / "dram-annotation/genes.faa",
        work / "dram-annotation/genes.fna",
        work / "dram-annotation/genes.gff",
        work / "dram-distill/genome_stats.tsv",
        work / "dram-distill/metabolism_summary.xlsx",
        work / "dram-distill/product.tsv",
        work / "dram-distill/product.html",
    )
    if any(not path.is_file() or path.stat().st_size == 0 for path in required):
        return False
    return (
        row_count(work / "dram-annotation/annotations.tsv") >= 30_000
        and row_count(work / "dram-distill/genome_stats.tsv") == EXPECTED_GENOMES
    )


def run_dram(args: argparse.Namespace, work: Path, logs: Path) -> list[dict[str, object]]:
    output = work / "dram-annotation"
    distill = work / "dram-distill"
    shards = work / "dram-shards"
    if output.exists() or distill.exists():
        if dram_branch_is_complete(work):
            return [
                existing_run_row("dram-annotate", [], logs),
                existing_run_row("dram-merge", [], logs),
                existing_run_row("dram-distill", [], logs),
            ]
        raise RuntimeError("Partial merged DRAM output exists without the Article 59 run marker")
    genome_files = sorted((work / "inputs/dram-mags").glob("*.fna"))
    if len(genome_files) != EXPECTED_GENOMES:
        raise RuntimeError(f"Expected {EXPECTED_GENOMES} DRAM inputs, observed {len(genome_files)}")
    if args.dram_jobs < 1 or args.dram_threads_per_job < 1:
        raise ValueError("DRAM jobs and threads per job must be positive")
    shards.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["PATH"] = str(args.dram_env / "bin") + os.pathsep + environment.get("PATH", "")
    environment["PYTHONNOUSERSITE"] = "1"
    environment["DRAM_CONFIG_LOCATION"] = str(work / "dram-config.json")
    shard_commands: list[str] = []
    reused_shards: list[str] = []
    for genome in genome_files:
        shard = shards / genome.stem
        annotation_table = shard / "annotations.tsv"
        required_shard_files = (
            annotation_table,
            shard / "genes.faa",
            shard / "genes.fna",
            shard / "genes.gff",
            shard / "scaffolds.fna",
        )
        if all(path.is_file() and path.stat().st_size > 0 for path in required_shard_files):
            reused_shards.append(genome.stem)
            continue
        if shard.exists():
            shutil.rmtree(shard)
        annotate = [
            str(args.dram_env / "bin/DRAM.py"),
            "annotate",
            "-i",
            str(genome),
            "-o",
            str(shards / genome.stem),
            "--min_contig_size",
            "1500",
            "--prodigal_mode",
            "meta",
            "--low_mem_mode",
            "--skip_trnascan",
            "--threads",
            str(args.dram_threads_per_job),
            "--config_loc",
            str(work / "dram-config.json"),
        ]
        shard_commands.append(shlex.join(annotate))
    command_file = work / "dram-annotate-commands.txt"
    command_file.write_text("\n".join(shard_commands) + "\n", encoding="utf-8")
    parallel_annotate = [
        str(args.metabolic_env / "bin/parallel"),
        "--jobs",
        str(args.dram_jobs),
        "--halt",
        "now,fail=1",
        "--joblog",
        str(logs / "dram-annotate.joblog.tsv"),
        "::::",
        str(command_file),
    ]
    if shard_commands:
        previous_joblog = logs / "dram-annotate.joblog.tsv"
        if previous_joblog.is_file():
            archived_joblog = logs / "dram-annotate.initial-interrupted.joblog.tsv"
            if not archived_joblog.exists():
                previous_joblog.replace(archived_joblog)
        annotate_row = timed_run("dram-annotate", parallel_annotate, work, logs, environment)
    else:
        annotate_row = {
            "Label": "dram-annotate",
            "Command": "all DRAM shards reused after output validation",
            "ExitStatus": 0,
            "Stdout": "not rerun",
            "Stderr": "not rerun",
            "Timing": "not rerun",
        }
    annotate_row["ReusedShards"] = ",".join(reused_shards)
    shard_tables = sorted(shards.glob("*/annotations.tsv"))
    if len(shard_tables) != EXPECTED_GENOMES or any(path.stat().st_size == 0 for path in shard_tables):
        raise RuntimeError("Parallel DRAM shard outputs are incomplete")
    merge_command = [
        str(args.dram_env / "bin/DRAM.py"),
        "merge_annotations",
        "-i",
        str(shards / "*"),
        "-o",
        str(output),
    ]
    merge_row = timed_run("dram-merge", merge_command, work, logs, environment)
    annotation_table = output / "annotations.tsv"
    if not annotation_table.is_file() or annotation_table.stat().st_size == 0:
        raise RuntimeError("DRAM annotation table is missing")
    distill_command = [
        str(args.dram_env / "bin/DRAM.py"),
        "distill",
        "-i",
        str(annotation_table),
        "-o",
        str(distill),
        "--config_loc",
        str(work / "dram-config.json"),
    ]
    distill_row = timed_run(
        "dram-distill", distill_command, work, logs, environment
    )
    return [annotate_row, merge_row, distill_row]


def metabolic_branch_is_complete(work: Path) -> bool:
    output = work / "metabolic-output"
    proteins = work / "metabolic-proteins"
    if not output.is_dir() or not proteins.is_dir():
        return False
    sheets = [
        output / "METABOLIC_result_each_spreadsheet" / f"METABOLIC_result_worksheet{index}.tsv"
        for index in range(1, 7)
    ]
    required = sheets + [output / "METABOLIC_result.xlsx", output / "METABOLIC_log.log"]
    if any(not path.is_file() or path.stat().st_size == 0 for path in required):
        return False
    if len(list(proteins.glob("*.faa"))) != EXPECTED_GENOMES:
        return False
    for protein in proteins.glob("*.faa"):
        prefix = f">{protein.stem}__"
        headers = [
            line for line in protein.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        ]
        if not headers or any(not header.startswith(prefix) for header in headers):
            return False
    result_dir = output / "KEGG_identifier_result"
    result_files = sorted(result_dir.glob("*.result.txt"))
    if len(result_files) != EXPECTED_GENOMES:
        return False
    if len(list(result_dir.glob("*.hits.txt"))) != EXPECTED_GENOMES:
        return False
    for result in result_files:
        rows = result.read_text(encoding="utf-8").splitlines()
        if len(rows) != 2_678:
            return False
        if not any(len(line.split("\t")) >= 2 and line.split("\t")[1].strip() for line in rows):
            return False
    log_text = (output / "METABOLIC_log.log").read_text(encoding="utf-8", errors="replace")
    return all(
        message in log_text
        for message in (
            "The hmmsearch is finished",
            "dbCAN2 searching is done",
            "MEROPS peptidase searching is done",
            "METABOLIC-G was done",
        )
    )


def enforce_unique_protein_ids(protein_files: list[Path]) -> dict[str, object]:
    """Prefix Prodigal IDs with the genome name and reject cross-genome collisions."""
    observed: set[str] = set()
    records = 0
    already_prefixed = 0
    for protein in protein_files:
        prefix = f"{protein.stem}__"
        temporary = protein.with_name(protein.name + ".prefixed.tmp")
        headers = 0
        with protein.open(encoding="utf-8") as source, temporary.open(
            "w", encoding="utf-8", newline="\n"
        ) as target:
            for raw in source:
                if raw.startswith(">"):
                    headers += 1
                    records += 1
                    identifier, separator, description = raw[1:].rstrip("\n").partition(" ")
                    if identifier.startswith(prefix):
                        already_prefixed += 1
                    else:
                        identifier = prefix + identifier
                    if identifier in observed:
                        raise RuntimeError(f"Duplicate protein identifier after prefixing: {identifier}")
                    observed.add(identifier)
                    target.write(">" + identifier)
                    if separator:
                        target.write(" " + description)
                    target.write("\n")
                else:
                    target.write(raw)
        if headers == 0:
            temporary.unlink(missing_ok=True)
            raise RuntimeError(f"No protein records in {protein}")
        temporary.replace(protein)
    return {
        "ProteinFiles": len(protein_files),
        "ProteinRecords": records,
        "InputIDsAlreadyGenomePrefixed": already_prefixed,
        "OutputIDsGenomePrefixed": records,
        "UniqueOutputProteinIDs": len(observed),
        "DuplicateOutputProteinIDs": records - len(observed),
        "PrefixPattern": "genome_id__prodigal_id",
    }


def run_metabolic(
    args: argparse.Namespace, work: Path, logs: Path
) -> list[dict[str, object]]:
    output = work / "metabolic-output"
    proteins = work / "metabolic-proteins"
    environment = os.environ.copy()
    environment["PATH"] = str(args.metabolic_env / "bin") + os.pathsep + environment.get("PATH", "")
    genome_files = sorted((work / "inputs/metabolic-mags").glob("*.fasta"))
    if len(genome_files) != EXPECTED_GENOMES:
        raise RuntimeError(
            f"Expected {EXPECTED_GENOMES} METABOLIC genome inputs, observed {len(genome_files)}"
        )
    command_file = work / "metabolic-prodigal-commands.txt"
    prodigal = args.metabolic_env / "bin/prodigal"
    prodigal_commands = []
    for genome in genome_files:
        stem = genome.stem
        command = [
            str(prodigal),
            "-i",
            str(genome),
            "-a",
            str(proteins / f"{stem}.faa"),
            "-o",
            str(proteins / f"{stem}.gff"),
            "-f",
            "gff",
            "-p",
            "meta",
            "-q",
        ]
        prodigal_commands.append(shlex.join(command))
    if not command_file.is_file():
        command_file.write_text("\n".join(prodigal_commands) + "\n", encoding="utf-8")
    parallel_command = [
        str(args.metabolic_env / "bin/parallel"),
        "--jobs",
        str(args.metabolic_threads),
        "--halt",
        "now,fail=1",
        "--joblog",
        str(logs / "metabolic-prodigal.joblog.tsv"),
        "::::",
        str(command_file),
    ]
    command = [
        str(args.metabolic_env / "bin/perl"),
        str(args.metabolic_dir / "METABOLIC-G.pl"),
        "-t",
        str(args.metabolic_threads),
        "-m-cutoff",
        "0.75",
        "-in",
        str(proteins),
        "-kofam-db",
        "small",
        "-o",
        str(output),
    ]
    if output.exists() or proteins.exists():
        if metabolic_branch_is_complete(work):
            return [
                existing_run_row("metabolic-prodigal", parallel_command, logs),
                existing_run_row("metabolic-g", command, logs),
            ]
        raise RuntimeError("Partial METABOLIC output exists without a complete branch")
    proteins.mkdir(parents=True)
    prodigal_row = timed_run("metabolic-prodigal", parallel_command, work, logs, environment)
    protein_files = sorted(proteins.glob("*.faa"))
    if len(protein_files) != EXPECTED_GENOMES or any(
        path.stat().st_size == 0 for path in protein_files
    ):
        raise RuntimeError("Standalone Prodigal protein outputs are incomplete")
    protein_id_audit = enforce_unique_protein_ids(protein_files)
    write_tsv(work / "metabolic-protein-id-audit.tsv", [protein_id_audit])
    metabolic_row = timed_run("metabolic-g", command, work, logs, environment)
    return [prodigal_row, metabolic_row]


def row_count(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def verify_outputs(work: Path) -> dict[str, object]:
    dram_required = [
        work / "dram-annotation/annotations.tsv",
        work / "dram-annotation/genes.faa",
        work / "dram-annotation/genes.fna",
        work / "dram-annotation/genes.gff",
        work / "dram-distill/genome_stats.tsv",
        work / "dram-distill/metabolism_summary.xlsx",
        work / "dram-distill/product.tsv",
        work / "dram-distill/product.html",
    ]
    metabolic_sheet_dir = work / "metabolic-output/METABOLIC_result_each_spreadsheet"
    metabolic_required = [
        metabolic_sheet_dir / f"METABOLIC_result_worksheet{index}.tsv"
        for index in range(1, 7)
    ] + [
        work / "metabolic-output/METABOLIC_result.xlsx",
        work / "metabolic-output/METABOLIC_run.log",
        work / "metabolic-output/METABOLIC_log.log",
    ]
    required = dram_required + metabolic_required
    incomplete = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if incomplete:
        raise RuntimeError(f"Required outputs are incomplete: {incomplete}")

    dram_annotations = row_count(work / "dram-annotation/annotations.tsv")
    dram_genomes = row_count(work / "dram-distill/genome_stats.tsv")
    result_dir = work / "metabolic-output/KEGG_identifier_result"
    metabolic_result_files = sorted(result_dir.glob("*.result.txt"))
    metabolic_hit_files = sorted(result_dir.glob("*.hits.txt"))
    protein_files = sorted((work / "metabolic-proteins").glob("*.faa"))
    if dram_annotations < 30_000:
        raise RuntimeError(f"Implausibly few DRAM annotations: {dram_annotations}")
    if dram_genomes != EXPECTED_GENOMES:
        raise RuntimeError(f"Expected {EXPECTED_GENOMES} DRAM genomes, observed {dram_genomes}")
    if len(metabolic_result_files) != EXPECTED_GENOMES or len(metabolic_hit_files) != EXPECTED_GENOMES:
        raise RuntimeError("METABOLIC did not produce one KO result/hit pair per genome")
    if len(protein_files) != EXPECTED_GENOMES or any(path.stat().st_size == 0 for path in protein_files):
        raise RuntimeError("METABOLIC Prodigal outputs are incomplete")
    protein_ids: set[str] = set()
    for protein in protein_files:
        prefix = f">{protein.stem}__"
        for line in protein.read_text(encoding="utf-8").splitlines():
            if not line.startswith(">"):
                continue
            if not line.startswith(prefix):
                raise RuntimeError(f"METABOLIC protein ID lacks a genome prefix: {protein.name}")
            identifier = line[1:].split(maxsplit=1)[0]
            if identifier in protein_ids:
                raise RuntimeError(f"Duplicate METABOLIC protein identifier: {identifier}")
            protein_ids.add(identifier)
    protein_id_audit = work / "metabolic-protein-id-audit.tsv"
    if not protein_id_audit.is_file() or protein_id_audit.stat().st_size == 0:
        raise RuntimeError("METABOLIC global protein-ID audit is missing")
    # METABOLIC ``*.result.txt`` tables are headerless: every physical line is
    # one KO-grid row.  Do not use row_count(), which intentionally subtracts
    # a header for DRAM/worksheet tables.
    minimum_ko_rows = min(
        sum(1 for _ in path.open(encoding="utf-8", errors="replace"))
        for path in metabolic_result_files
    )
    if minimum_ko_rows < 2_500:
        raise RuntimeError(f"METABOLIC KO tables are unexpectedly short: {minimum_ko_rows}")
    empty_hit_tables = []
    for path in metabolic_result_files:
        rows = path.read_text(encoding="utf-8").splitlines()
        if len(rows) != 2_678 or not any(
            len(line.split("\t")) >= 2 and line.split("\t")[1].strip()
            for line in rows
        ):
            empty_hit_tables.append(path.name)
    if empty_hit_tables:
        raise RuntimeError(
            "METABOLIC KO tables lack genome-specific hits; check globally unique protein IDs: "
            f"{empty_hit_tables}"
        )

    log_text = (work / "metabolic-output/METABOLIC_log.log").read_text(
        encoding="utf-8", errors="replace"
    )
    required_messages = (
        "The hmmsearch is finished",
        "dbCAN2 searching is done",
        "MEROPS peptidase searching is done",
        "METABOLIC-G was done",
    )
    missing_messages = [message for message in required_messages if message not in log_text]
    if missing_messages:
        raise RuntimeError(f"METABOLIC completion messages are missing: {missing_messages}")
    fatal_tokens = ("No such file or directory", "command not found", "Segmentation fault")
    observed_fatal = [token for token in fatal_tokens if token in log_text]
    if observed_fatal:
        raise RuntimeError(f"METABOLIC log contains fatal tokens: {observed_fatal}")

    return {
        "expected_genomes": EXPECTED_GENOMES,
        "dram_annotation_rows": dram_annotations,
        "dram_genome_rows": dram_genomes,
        "metabolic_result_files": len(metabolic_result_files),
        "metabolic_minimum_ko_rows": minimum_ko_rows,
        "metabolic_worksheet1_rows": row_count(metabolic_required[0]),
        "metabolic_worksheet2_rows": row_count(metabolic_required[1]),
        "metabolic_worksheet3_rows": row_count(metabolic_required[2]),
        "metabolic_worksheet4_rows": row_count(metabolic_required[3]),
        "metabolic_completion_messages": list(required_messages),
    }


def main() -> None:
    args = parse_args()
    args.work_dir = args.work_dir.resolve()
    args.metabolic_dir = args.metabolic_dir.resolve()
    args.dram_env = args.dram_env.resolve()
    args.metabolic_env = args.metabolic_env.resolve()
    work = args.work_dir
    required_markers = (
        ".article59-inputs-complete",
        ".article59-dram-database-complete",
        ".article59-metabolic-database-complete",
    )
    missing_markers = [name for name in required_markers if not (work / name).is_file()]
    if missing_markers:
        raise FileNotFoundError(f"Run setup steps first: {missing_markers}")
    marker = work / ".article59-runs-complete"
    if marker.is_file():
        contract = verify_outputs(work)
        contract_path = work / "run-contract.json"
        if contract_path.is_file():
            recorded = json.loads(contract_path.read_text(encoding="utf-8"))
            recorded.update(contract)
            contract = recorded
            contract_path.write_text(
                json.dumps(contract, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        print(json.dumps(contract, indent=2, sort_keys=True))
        return

    logs = work / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures = {
            pool.submit(run_dram, args, work, logs): "DRAM",
            pool.submit(run_metabolic, args, work, logs): "METABOLIC",
        }
        command_rows: list[dict[str, object]] = []
        for future in concurrent.futures.as_completed(futures):
            label = futures[future]
            try:
                command_rows.extend(future.result())
            except Exception as error:
                raise RuntimeError(f"{label} branch failed") from error

    command_rows.sort(key=lambda row: str(row["Label"]))
    write_tsv(work / "command-log.tsv", command_rows)
    contract = verify_outputs(work)
    contract.update(
        {
            "article": 59,
            "dram_jobs": args.dram_jobs,
            "dram_threads_per_job": args.dram_threads_per_job,
            "dram_total_thread_ceiling": args.dram_jobs * args.dram_threads_per_job,
            "metabolic_threads": args.metabolic_threads,
            "parallel_tool_branches": True,
            "module_cutoff": 0.75,
            "minimum_contig_bp": 1500,
            "prodigal_mode": "meta",
            "metabolic_input_mode": "supported-protein-input-after-pinned-prodigal",
            "metabolic_protein_ids_globally_prefixed": True,
            "kofam_mode_metabolic": "small-2643-compatible-profiles",
            "python_user_site_disabled_for_dram": True,
            "random_processes": "none after deterministic truncation seed 59002",
        }
    )
    (work / "run-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    marker.write_text("DRAM and METABOLIC Article 59 runs completed\n", encoding="utf-8")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
