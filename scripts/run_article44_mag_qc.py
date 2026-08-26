#!/usr/bin/env python3
"""Run the complete Article 44 MAG quality, marker, and assembly-graph audit."""

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
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--assembly-env", type=Path, required=True)
    parser.add_argument("--mag-qc-env", type=Path, required=True)
    parser.add_argument("--gunc-env", type=Path, required=True)
    parser.add_argument("--checkm1-env", type=Path, required=True)
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


def save_log(work: Path, rows: list[dict[str, object]]) -> None:
    if rows:
        write_tsv(work / "command-log.tsv", rows)


def run_timed(
    label: str,
    command: list[str],
    work: Path,
    runtime: dict[str, str],
    rows: list[dict[str, object]],
    stdout_target: Path | None = None,
) -> None:
    stdout_path = stdout_target or (work / "logs" / f"{label}.stdout.log")
    stderr_path = work / "logs" / f"{label}.stderr.log"
    resource_path = work / "logs" / f"{label}.time.txt"
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        result = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(resource_path), *command],
            stdout=stdout, stderr=stderr, env=runtime, text=True,
        )
    ended = datetime.now(timezone.utc)
    rows.append({
        "Label": label, "Command": shlex.join(command), "ReturnCode": result.returncode,
        "StartedUTC": started.isoformat(), "EndedUTC": ended.isoformat(),
        "Stdout": str(stdout_path), "Stderr": str(stderr_path), "ResourceLog": str(resource_path),
    })
    save_log(work, rows)
    if result.returncode:
        tail = stderr_path.read_text(encoding="utf-8", errors="replace")[-8000:]
        raise RuntimeError(f"{label} failed ({result.returncode})\n{tail}")


RANK_CODES = {
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


def parse_checkm_domains(path: Path, taxon_marker_sets: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    taxon_domains: dict[tuple[str, str], str] = {}
    with taxon_marker_sets.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 3:
                continue
            rank, taxon, lineage = row[:3]
            domain = lineage.split(";", 1)[0]
            if domain not in {"Archaea", "Bacteria"}:
                continue
            key = (rank, taxon)
            previous = taxon_domains.get(key)
            if previous is not None and previous != domain:
                raise RuntimeError(f"Ambiguous CheckM taxon domain for {key}: {previous}, {domain}")
            taxon_domains[key] = domain
    rows = list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))
    result: dict[str, str] = {}
    audit: list[dict[str, str]] = []
    for row in rows:
        name = row.get("Bin Id") or row.get("Bin ID") or row.get("Name")
        lineage = row.get("Marker lineage") or row.get("Marker Lineage") or ""
        if not name:
            raise RuntimeError(f"Cannot identify CheckM1 bin column: {list(row)}")
        marker = lineage.split(" (UID", 1)[0]
        if marker == "k__Archaea":
            rank, taxon, domain = "domain", "Archaea", "Archaea"
        elif marker == "k__Bacteria":
            rank, taxon, domain = "domain", "Bacteria", "Bacteria"
        elif len(marker) > 3 and marker[1:3] == "__" and marker[0] in RANK_CODES:
            rank, taxon = RANK_CODES[marker[0]], marker[3:]
            domain = taxon_domains.get((rank, taxon), "")
            if not domain:
                raise RuntimeError(
                    f"Cannot infer domain for CheckM marker lineage {lineage!r}; "
                    f"missing ({rank}, {taxon}) in {taxon_marker_sets}"
                )
        else:
            raise RuntimeError(f"Unsupported CheckM marker lineage for domain selection: {lineage!r}")
        mag = name.removesuffix(".fna")
        result[mag] = domain
        audit.append({
            "MAG": mag,
            "MarkerLineage": lineage,
            "MarkerRank": rank,
            "MarkerTaxon": taxon,
            "Domain": domain,
            "DomainSource": "CheckM1 marker lineage plus taxon_marker_sets.tsv",
        })
    return result, audit


def main() -> int:
    args = parse_args()
    root, work = args.project_root.resolve(), args.work_dir.resolve()
    assembly_env, mag_qc_env = args.assembly_env.resolve(), args.mag_qc_env.resolve()
    gunc_env, checkm1_env = args.gunc_env.resolve(), args.checkm1_env.resolve()
    checkm2_db, gunc_db = args.checkm2_db.resolve(), args.gunc_db.resolve()
    if not (work / ".article44-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article44_mag_qc.py first")
    if (work / ".article44-run-complete").is_file():
        print(f"Article 44 run already complete: {work}")
        return 0
    observed_checkm2_sha = sha256(checkm2_db) if checkm2_db.is_file() else ""
    observed_gunc_sha = sha256(gunc_db) if gunc_db.is_file() else ""
    if observed_checkm2_sha != CHECKM2_DB_SHA256:
        raise RuntimeError("CheckM2 version-3 database is missing or has checksum drift")
    if observed_gunc_sha != GUNC_DB_SHA256:
        raise RuntimeError("GUNC ProGenomes 2.1 database is missing or has checksum drift")
    bins = work / "bins"
    bin_fastas = sorted(bins.glob("*.fna"))
    if not bin_fastas:
        raise RuntimeError("No selected Article 44 MAG candidates")
    for directory in (work / "logs", work / "qc", work / "features", work / "graph", work / "tmp"):
        directory.mkdir(parents=True, exist_ok=True)

    runtime = os.environ.copy()
    for key in ("PERL5LIB", "PERL_LOCAL_LIB_ROOT", "PERL_MB_OPT", "PERL_MM_OPT"):
        runtime.pop(key, None)
    runtime.update({
        "LC_ALL": "C", "TZ": "UTC", "PYTHONHASHSEED": "0", "PYTHONNOUSERSITE": "1",
        "OMP_NUM_THREADS": str(args.threads), "OPENBLAS_NUM_THREADS": str(args.threads),
        "MKL_NUM_THREADS": str(args.threads), "MPLCONFIGDIR": str(work / "tmp/matplotlib"),
    })
    (work / "tmp/matplotlib").mkdir(exist_ok=True)
    rows: list[dict[str, object]] = load_log(work / "command-log.tsv")

    primary_done = work / ".article44-primary-qc-complete"
    if not primary_done.is_file():
        for partial in (work / "qc/checkm2", work / "qc/gunc", work / "tmp/checkm2", work / "tmp/gunc"):
            if partial.exists():
                shutil.rmtree(partial)
            partial.mkdir(parents=True)
        run_timed(
            "checkm2-selected-mags",
            [
                str(mag_qc_env / "bin/checkm2"), "predict", "--input", str(bins),
                "--output-directory", str(work / "qc/checkm2"), "--extension", ".fna",
                "--threads", str(args.threads), "--database_path", str(checkm2_db),
                "--tmpdir", str(work / "tmp/checkm2"), "--remove_intermediates",
            ],
            work, {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"}, rows,
        )
        quality_rows = list(csv.DictReader((work / "qc/checkm2/quality_report.tsv").open(encoding="utf-8"), delimiter="\t"))
        if len(quality_rows) != len(bin_fastas):
            raise RuntimeError("Article 44 CheckM2 row count mismatch")
        run_timed(
            "gunc-selected-mags",
            [
                str(gunc_env / "bin/gunc"), "run", "--input_dir", str(bins),
                "--file_suffix", ".fna", "--db_file", str(gunc_db), "--threads", str(args.threads),
                "--out_dir", str(work / "qc/gunc"), "--temp_dir", str(work / "tmp/gunc"),
                "--contig_taxonomy_output",
            ],
            work, {**runtime, "PATH": f"{gunc_env / 'bin'}:/usr/bin:/bin"}, rows,
        )
        maxcss = sorted((work / "qc/gunc").glob("GUNC.*.maxCSS_level.tsv"))
        if len(maxcss) != 1 or len(list(csv.DictReader(maxcss[0].open(encoding="utf-8"), delimiter="\t"))) != len(bin_fastas):
            raise RuntimeError("Article 44 GUNC row count mismatch")
        primary_done.write_text("PASS\n", encoding="utf-8")

    checkm1_done = work / ".article44-checkm1-complete"
    checkm1_table = work / "qc/checkm1-qa.tsv"
    if not checkm1_done.is_file():
        for partial in (work / "qc/checkm1", work / "tmp/checkm1"):
            if partial.exists():
                shutil.rmtree(partial)
            partial.mkdir(parents=True)
        run_timed(
            "checkm1-lineage-selected-mags",
            [
                str(checkm1_env / "bin/checkm"), "lineage_wf", "-x", "fna", "-t", str(args.threads),
                "--pplacer_threads", "1", "--aai_strain", "0.9", "--tab_table", "-f", str(checkm1_table),
                "--tmpdir", str(work / "tmp/checkm1"), str(bins), str(work / "qc/checkm1"),
            ],
            work, {**runtime, "PATH": f"{checkm1_env / 'bin'}:/usr/bin:/bin"}, rows,
        )
        domains, _ = parse_checkm_domains(checkm1_table, checkm1_env / "checkm_data/taxon_marker_sets.tsv")
        if len(domains) != len(bin_fastas):
            raise RuntimeError("Article 44 CheckM1 row count mismatch")
        checkm1_done.write_text("PASS\n", encoding="utf-8")
    domains, domain_audit = parse_checkm_domains(
        checkm1_table, checkm1_env / "checkm_data/taxon_marker_sets.tsv"
    )
    if set(domains) != {path.stem for path in bin_fastas}:
        raise RuntimeError("Article 44 CheckM1 domain coordinate mismatch")
    write_tsv(work / "qc-domain-audit.tsv", domain_audit)

    feature_done = work / ".article44-feature-audit-complete"
    if not feature_done.is_file():
        feature_root = work / "features"
        if feature_root.exists():
            shutil.rmtree(feature_root)
        feature_root.mkdir(parents=True)
        feature_prefixes = ("barrnap-", "trnascan-", "prodigal-")
        rows = [row for row in rows if not str(row.get("Label", "")).startswith(feature_prefixes)]
        for path in (work / "logs").glob("*"):
            if path.is_file() and path.name.startswith(feature_prefixes):
                path.unlink()
        save_log(work, rows)
        for fasta in bin_fastas:
            mag = fasta.stem
            mag_dir = feature_root / mag
            mag_dir.mkdir(parents=True, exist_ok=True)
            domain = domains[mag]
            kingdom = "arc" if domain == "Archaea" else "bac"
            trna_mode = "-A" if domain == "Archaea" else "-B"
            barrnap_gff = mag_dir / "barrnap.gff"
            run_timed(
                f"barrnap-{mag}",
                [str(mag_qc_env / "bin/barrnap"), "--kingdom", kingdom, "--threads", str(min(args.threads, 8)), "--addids", str(fasta)],
                work, {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"}, rows,
                stdout_target=barrnap_gff,
            )
            trna_table = mag_dir / "trnascan.tsv"
            run_timed(
                f"trnascan-{mag}",
                [
                    str(mag_qc_env / "bin/tRNAscan-SE"), trna_mode, "--thread", str(min(args.threads, 8)),
                    "--output", str(trna_table), "--stats", str(mag_dir / "trnascan.stats"),
                    "--gff", str(mag_dir / "trnascan.gff"), "--log", str(mag_dir / "trnascan.log"),
                    "--forceow", "--quiet", str(fasta),
                ],
                work, {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"}, rows,
            )
            run_timed(
                f"prodigal-{mag}",
                [
                    str(mag_qc_env / "bin/prodigal"), "-i", str(fasta), "-p", "meta", "-f", "gff",
                    "-o", str(mag_dir / "prodigal.gff"), "-a", str(mag_dir / "proteins.faa"),
                ],
                work, {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"}, rows,
            )
            if not barrnap_gff.is_file() or not trna_table.is_file() or not (mag_dir / "prodigal.gff").is_file():
                raise RuntimeError(f"Incomplete marker feature output: {mag}")
        feature_done.write_text("PASS\n", encoding="utf-8")

    graph_done = work / ".article44-graph-complete"
    graph_fastg = work / "graph/megahit-k141.fastg"
    if not graph_done.is_file():
        k141 = root / "data/raw/article30/work/assemblies/megahit-coassembly/intermediate_contigs/k141.contigs.fa"
        run_timed(
            "megahit-k141-fastg",
            [str(assembly_env / "bin/megahit_toolkit"), "contig2fastg", "141", str(k141)],
            work, {**runtime, "PATH": f"{assembly_env / 'bin'}:/usr/bin:/bin"}, rows,
            stdout_target=graph_fastg,
        )
        if not graph_fastg.is_file() or graph_fastg.stat().st_size == 0:
            raise RuntimeError("MEGAHIT graph conversion produced no FASTG")
        graph_done.write_text("PASS\n", encoding="utf-8")

    maxcss = sorted((work / "qc/gunc").glob("GUNC.*.maxCSS_level.tsv"))
    checkm_manifest = checkm1_env / "checkm_data/.dmanifest"
    write_tsv(work / "qc-database-audit.tsv", [
        {"Tool": "CheckM2", "ToolVersion": "1.1.0", "Database": "DIAMOND database version 3", "Reference": "Zenodo 14897628", "Path": str(checkm2_db), "Bytes": checkm2_db.stat().st_size, "MD5": md5(checkm2_db), "SHA256": observed_checkm2_sha, "Integrity": "PASS"},
        {"Tool": "GUNC", "ToolVersion": "1.1.0", "Database": "ProGenomes 2.1", "Reference": "GUNC official download_db endpoint", "Path": str(gunc_db), "Bytes": gunc_db.stat().st_size, "MD5": md5(gunc_db), "SHA256": observed_gunc_sha, "Integrity": "PASS (upstream MD5 checked during download)"},
        {"Tool": "CheckM1", "ToolVersion": "1.2.5", "Database": "CheckM reference data 2015-01-16", "Reference": "CheckM post-link installer", "Path": str(checkm_manifest), "Bytes": checkm_manifest.stat().st_size, "MD5": md5(checkm_manifest), "SHA256": sha256(checkm_manifest), "Integrity": "PASS (.dmanifest retained)"},
    ])
    payload = {
        "article": 44, "selected_bins": len(bin_fastas),
        "checkm2_quality": str(work / "qc/checkm2/quality_report.tsv"), "gunc_maxcss": str(maxcss[0]),
        "checkm1_quality": str(checkm1_table), "feature_directories": len(list((work / "features").iterdir())),
        "graph_fastg": str(graph_fastg), "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    (work / ".article44-run-complete").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
