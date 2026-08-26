#!/usr/bin/env python3
"""Prepare checksum-gated MAGs and databases for Article 59."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import random
import subprocess
from pathlib import Path
from typing import Iterable, Iterator


PRIMARY_MAG_COUNT = 24
TRUNCATION_SEED = 59002
TRUNCATION_FRACTIONS = (0.50, 0.70, 0.90, 1.00)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--kofam-dir", type=Path, required=True)
    parser.add_argument("--metabolic-dir", type=Path, required=True)
    parser.add_argument("--merops-file", type=Path, required=True)
    parser.add_argument("--dbcan-hmm", type=Path, required=True)
    parser.add_argument("--dram-source", type=Path, required=True)
    parser.add_argument("--dram-env", type=Path, required=True)
    parser.add_argument("--metabolic-env", type=Path, required=True)
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


def digest(path: Path, algorithm: str = "sha256") -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def digest_gzip_payload(path: Path) -> str:
    hasher = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def digest_record_set(records: Iterable[tuple[str, str]]) -> str:
    """Hash a FASTA record set independently of record order and line wrapping."""
    hasher = hashlib.sha256()
    for header, sequence in sorted(records, key=lambda record: record[0]):
        hasher.update(header.encode("utf-8"))
        hasher.update(b"\0")
        hasher.update(sequence.encode("ascii"))
        hasher.update(b"\0")
    return hasher.hexdigest()


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    parts: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> tuple[int, int]:
    count = 0
    bases = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            count += 1
            bases += len(sequence)
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    if count == 0 or bases == 0:
        raise ValueError(f"Empty FASTA output: {path}")
    return count, bases


def link_exact(source: Path, destination: Path) -> None:
    source = source.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        if destination.resolve() != source:
            raise RuntimeError(f"Wrong symlink target: {destination}")
        return
    if destination.exists():
        raise FileExistsError(f"Expected a symlink, found another object: {destination}")
    destination.symlink_to(source)


def capture(command: list[str], environment: dict[str, str] | None = None) -> str:
    completed = subprocess.run(
        command,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if completed.returncode != 0 or not lines:
        raise RuntimeError(f"Version command failed: {' '.join(command)}\n{completed.stdout}")
    return lines[0]


def asset_paths(args: argparse.Namespace) -> dict[str, Path]:
    forms = args.dram_source / "data"
    return {
        "profiles.tar.gz": args.kofam_dir / "profiles.tar.gz",
        "ko_list.gz": args.kofam_dir / "ko_list.gz",
        "METABOLIC_hmm_db.tgz": args.metabolic_dir / "METABOLIC_hmm_db.tgz",
        "METABOLIC_template_and_database.tgz": args.metabolic_dir
        / "METABOLIC_template_and_database.tgz",
        "Accessory_scripts.tgz": args.metabolic_dir / "Accessory_scripts.tgz",
        "dbCAN.hmm": args.dbcan_hmm,
        "pepunit.lib": args.merops_file,
        "genome_summary_form.tsv": forms / "genome_summary_form.tsv",
        "module_step_form.tsv": forms / "module_step_form.tsv",
        "etc_module_database.tsv": forms / "etc_module_database.tsv",
        "function_heatmap_form.tsv": forms / "function_heatmap_form.tsv",
        "amg_database.tsv": forms / "amg_database.tsv",
    }


def audit_databases(args: argparse.Namespace, work: Path) -> list[dict[str, object]]:
    manifest = read_tsv(
        args.project_root / "data/small/59-metabolism-database-manifest.tsv"
    )
    paths = asset_paths(args)
    rows: list[dict[str, object]] = []
    for row in manifest:
        path = paths[row["Asset"]]
        if not path.is_file():
            raise FileNotFoundError(path)
        observed_bytes = path.stat().st_size
        observed_sha = digest(path)
        passed = (
            observed_bytes == int(row["ExpectedBytes"])
            and observed_sha == row["SHA256"]
        )
        rows.append(
            {
                "Asset": row["Asset"],
                "Tool": row["Tool"],
                "Release": row["Release"],
                "ExpectedBytes": row["ExpectedBytes"],
                "ObservedBytes": observed_bytes,
                "ExpectedSHA256": row["SHA256"],
                "ObservedSHA256": observed_sha,
                "ChecksumPass": str(passed).lower(),
                "Redistribution": row["Redistribution"],
            }
        )
        if not passed:
            raise RuntimeError(f"Database identity gate failed for {path}")
    write_tsv(work / "database-audit.tsv", rows)
    return rows


def prepare_mags(root: Path, work: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    taxonomy_path = (
        root / "data/small/46-gtdbtk-taxonomy-frozen/taxonomy-summary.tsv"
    )
    genome_dir = (
        root / "data/small/45-drep-dereplication-frozen/representative-genomes"
    )
    taxonomy = read_tsv(taxonomy_path)
    if len(taxonomy) != PRIMARY_MAG_COUNT:
        raise RuntimeError(f"Expected {PRIMARY_MAG_COUNT} taxonomy rows")

    dram_input = work / "inputs/dram-mags"
    metabolic_input = work / "inputs/metabolic-mags"
    primary_rows: list[dict[str, object]] = []
    source_records: dict[str, list[tuple[str, str]]] = {}
    for row in taxonomy:
        sgb = row["SGB"]
        representative = row["Representative"]
        source = genome_dir / f"{representative}.gz"
        if not source.is_file():
            raise FileNotFoundError(source)
        payload_sha = digest_gzip_payload(source)
        if payload_sha != row["RepresentativeSHA256"]:
            raise RuntimeError(f"Representative SHA mismatch for {sgb}")
        records = list(fasta_records(source))
        if not records:
            raise RuntimeError(f"No sequences in {source}")
        source_records[sgb] = records
        dram_path = dram_input / f"{sgb}.fna"
        metabolic_path = metabolic_input / f"{sgb}.fasta"
        contigs, bases = write_fasta(dram_path, records)
        write_fasta(
            metabolic_path,
            ((f"{sgb}__{header}", sequence) for header, sequence in records),
        )
        primary_rows.append(
            {
                "Genome": sgb,
                "AnalysisSet": "Primary real MAG",
                "Representative": representative,
                "RepresentativeSHA256": payload_sha,
                "SequenceSetSHA256": digest_record_set(records),
                "DRAMInputFileSHA256": digest(dram_path),
                "GTDBRelease": row["GTDBRelease"],
                "Species": row["Species"].removeprefix("s__") or "Unresolved",
                "Domain": row["Domain"].removeprefix("d__"),
                "Phylum": row["Phylum"].removeprefix("p__"),
                "CompletenessPct": row["Completeness"],
                "ContaminationPct": row["Contamination"],
                "Contigs": contigs,
                "GenomeBp": bases,
                "RetentionTargetPct": "100",
                "RetentionObservedPct": "100",
                "ParentGenome": sgb,
                "DeterministicSeed": "NA",
            }
        )

    parent = "SGB_002"
    records = list(source_records[parent])
    total_bp = sum(len(sequence) for _, sequence in records)
    order = list(range(len(records)))
    random.Random(TRUNCATION_SEED).shuffle(order)
    truncation_rows: list[dict[str, object]] = []
    parent_meta = next(row for row in primary_rows if row["Genome"] == parent)
    for fraction in TRUNCATION_FRACTIONS:
        target_bp = total_bp * fraction
        retained: list[tuple[str, str]] = []
        retained_bp = 0
        for index in order:
            retained.append(records[index])
            retained_bp += len(records[index][1])
            if retained_bp >= target_bp:
                break
        label = f"TRUNC_{int(fraction * 100):03d}"
        dram_path = dram_input / f"{label}.fna"
        metabolic_path = metabolic_input / f"{label}.fasta"
        contigs, bases = write_fasta(dram_path, retained)
        write_fasta(
            metabolic_path,
            ((f"{label}__{header}", sequence) for header, sequence in retained),
        )
        truncation_rows.append(
            {
                "Genome": label,
                "AnalysisSet": "Deterministic truncation sensitivity",
                "Representative": parent_meta["Representative"],
                "RepresentativeSHA256": parent_meta["RepresentativeSHA256"],
                "SequenceSetSHA256": digest_record_set(retained),
                "DRAMInputFileSHA256": digest(dram_path),
                "GTDBRelease": parent_meta["GTDBRelease"],
                "Species": parent_meta["Species"],
                "Domain": parent_meta["Domain"],
                "Phylum": parent_meta["Phylum"],
                "CompletenessPct": "NA",
                "ContaminationPct": "NA",
                "Contigs": contigs,
                "GenomeBp": bases,
                "RetentionTargetPct": f"{100 * fraction:.0f}",
                "RetentionObservedPct": f"{100 * bases / total_bp:.6f}",
                "ParentGenome": parent,
                "DeterministicSeed": TRUNCATION_SEED,
            }
        )

    all_rows = primary_rows + truncation_rows
    write_tsv(work / "input-mag-ledger.tsv", all_rows)
    write_tsv(work / "truncation-ledger.tsv", truncation_rows)
    return primary_rows, truncation_rows


def prepare_metabolic_links(
    args: argparse.Namespace, work: Path
) -> list[dict[str, object]]:
    requested_path = args.metabolic_dir / "All_Module_KO_ids.txt"
    requested = [
        line.strip()
        for line in requested_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(requested) != len(set(requested)):
        raise RuntimeError("METABOLIC All_Module_KO_ids contains duplicates")
    source_profiles = args.kofam_dir / "profiles"
    target_profiles = args.metabolic_dir / "kofam_database/profiles"
    target_profiles.mkdir(parents=True, exist_ok=True)
    compatibility: list[dict[str, object]] = []
    retained: list[str] = []
    for profile in requested:
        source = source_profiles / profile
        present = source.is_file()
        status = "included" if present else "excluded_missing_in_locked_snapshot"
        compatibility.append(
            {
                "Profile": profile,
                "RequestedByMETABOLICv4": "true",
                "PresentInKOfam2026_06_01": str(present).lower(),
                "CompatibilityAction": status,
            }
        )
        if present:
            retained.append(profile)
            link_exact(source, target_profiles / profile)
    if not retained or len(requested) - len(retained) > 5:
        raise RuntimeError("Unexpectedly large KOfam/METABOLIC profile mismatch")
    (target_profiles / "All_Module_KO_ids.txt").write_text(
        "\n".join(retained) + "\n", encoding="utf-8"
    )
    with gzip.open(args.kofam_dir / "ko_list.gz", "rt", encoding="utf-8") as source:
        ko_text = source.read()
    (args.metabolic_dir / "kofam_database/ko_list").write_text(
        ko_text, encoding="utf-8"
    )
    link_exact(
        args.dbcan_hmm,
        args.metabolic_dir / "dbCAN2/dbCAN-fam-HMMs.txt",
    )
    link_exact(args.merops_file, args.metabolic_dir / "MEROPS/pepunit.lib")
    write_tsv(work / "kofam-compatibility-audit.tsv", compatibility)
    return compatibility


def smoke_versions(args: argparse.Namespace, work: Path) -> list[dict[str, object]]:
    dram_env = os.environ.copy()
    dram_env["PYTHONNOUSERSITE"] = "1"
    metabolic_env = os.environ.copy()
    metabolic_env["PATH"] = str(args.metabolic_env / "bin") + os.pathsep + metabolic_env.get("PATH", "")
    rows = [
        {
            "Tool": "DRAM",
            "VersionEvidence": capture(
                [str(args.dram_env / "bin/DRAM-setup.py"), "version"], dram_env
            ),
        },
        {
            "Tool": "Python isolation",
            "VersionEvidence": capture(
                [
                    str(args.dram_env / "bin/python"),
                    "-c",
                    "import numpy,skbio; print(f'numpy={numpy.__version__};skbio={skbio.__version__};user_site_disabled')",
                ],
                dram_env,
            ),
        },
        {
            "Tool": "METABOLIC-G",
            "VersionEvidence": capture(
                [
                    str(args.metabolic_env / "bin/perl"),
                    str(args.metabolic_dir / "METABOLIC-G.pl"),
                    "--version",
                ],
                metabolic_env,
            ),
        },
        {
            "Tool": "HMMER",
            "VersionEvidence": capture(
                [str(args.metabolic_env / "bin/hmmsearch"), "-h"], metabolic_env
            ),
        },
        {
            "Tool": "DIAMOND",
            "VersionEvidence": capture(
                [str(args.metabolic_env / "bin/diamond"), "version"], metabolic_env
            ),
        },
        {
            "Tool": "Prodigal",
            "VersionEvidence": capture(
                [str(args.metabolic_env / "bin/prodigal"), "-v"], metabolic_env
            ),
        },
    ]
    write_tsv(work / "tool-smoke.tsv", rows)
    return rows


def main() -> None:
    args = parse_args()
    args.project_root = args.project_root.resolve()
    args.work_dir = args.work_dir.resolve()
    args.kofam_dir = args.kofam_dir.resolve()
    args.metabolic_dir = args.metabolic_dir.resolve()
    args.merops_file = args.merops_file.resolve()
    args.dbcan_hmm = args.dbcan_hmm.resolve()
    args.dram_source = args.dram_source.resolve()
    args.dram_env = args.dram_env.resolve()
    args.metabolic_env = args.metabolic_env.resolve()
    work = args.work_dir
    work.mkdir(parents=True, exist_ok=True)

    database_rows = audit_databases(args, work)
    primary, truncations = prepare_mags(args.project_root, work)
    compatibility = prepare_metabolic_links(args, work)
    versions = smoke_versions(args, work)

    contract = {
        "article": 59,
        "primary_real_mags": len(primary),
        "sensitivity_genomes": len(truncations),
        "truncation_parent": "SGB_002",
        "truncation_seed": TRUNCATION_SEED,
        "kofam_release": "2026-06-01",
        "metabolic_commit": "97236332519180f1d76a242dedb0aaa8191fdbb3",
        "dram_release": "1.5.0",
        "dram_source_commit": "fe61d759303f30db058d5d505c448b28e41b03f1",
        "database_assets_passed": len(database_rows),
        "metabolic_requested_kofam_profiles": len(compatibility),
        "metabolic_available_kofam_profiles": sum(
            row["PresentInKOfam2026_06_01"] == "true" for row in compatibility
        ),
        "tool_smoke_rows": len(versions),
        "python_user_site_disabled_for_dram": True,
    }
    (work / "preparation-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (work / ".article59-inputs-complete").write_text(
        "checksum-gated Article 59 inputs prepared\n", encoding="utf-8"
    )
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
