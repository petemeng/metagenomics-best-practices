#!/usr/bin/env python3
"""Create a compact, checksum-covered Article 60 GEM evidence bundle."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
from pathlib import Path


SUMMARY_FILES = (
    "model-structure-summary.tsv",
    "medium-feasibility.tsv",
    "gapfill-burden.tsv",
    "determinism-control.tsv",
    "truncation-sensitivity.tsv",
    "resource-usage.tsv",
    "evidence-ladder.tsv",
    "summary.json",
)

LEDGER_FILES = (
    "input-mag-ledger.tsv",
    "truncation-ledger.tsv",
    "protein-id-audit.tsv",
    "prodigal-command-log.tsv",
    "tool-smoke.tsv",
    "preparation-contract.json",
    "run-contract.json",
    "command-log.tsv",
    "model-plan.tsv",
)

SCRIPT_FILES = (
    "download_article60_gem_db.py",
    "prepare_article60_gem.py",
    "run_article60_gem.py",
    "summarize_article60_gem.py",
    "plot_article60_gem.py",
    "freeze_article60_gem.py",
    "validate_article60_gem.py",
    "article42_44_validation_utils.py",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gapseq-env", type=Path, required=True)
    parser.add_argument("--carveme-env", type=Path, required=True)
    parser.add_argument("--database-cache", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def gzip_copy(source: Path, target: Path, *, allow_empty: bool = False) -> None:
    if not source.is_file() or (source.stat().st_size == 0 and not allow_empty):
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as source_handle, target.open("wb") as raw_target:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_target, mtime=0
        ) as gzip_handle:
            shutil.copyfileobj(source_handle, gzip_handle, length=8 * 1024 * 1024)


def carveme_root(environment: Path) -> Path:
    candidates = sorted(
        {
            path.resolve()
            for path in (environment / "lib").glob("python*/site-packages/carveme")
        }
    )
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one CarveMe package root, observed {candidates}")
    return candidates[0]


def model_source(work: Path, genome: str, tool: str, stage: str) -> Path:
    if tool == "gapseq":
        base = work / "gapseq" / genome
        if stage == "draft":
            return base / f"{genome}-draft.xml"
        return base / "filled-permissive" / f"{genome}.xml"
    base = work / "carveme" / genome
    if stage == "draft":
        return base / f"{genome}-draft.xml"
    return base / f"{genome}-filled-LB.xml"


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    target = args.output_dir.resolve()
    gapseq_env = args.gapseq_env.resolve()
    carveme_env = args.carveme_env.resolve()
    database_cache = args.database_cache.resolve()
    if not (work / ".article60-models-complete").is_file():
        raise FileNotFoundError("Article 60 model sentinel is missing")
    for name in SUMMARY_FILES:
        if not (work / name).is_file():
            raise FileNotFoundError(f"Run the Article 60 summarizer first: {name}")
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)

    for name in (*SUMMARY_FILES, *LEDGER_FILES):
        copy(work / name, target / name)

    ledger = read_tsv(work / "input-mag-ledger.tsv")
    genomes = [row["Genome"] for row in ledger]
    if len(genomes) != 12 or len(set(genomes)) != 12:
        raise RuntimeError(f"Expected 12 unique Article 60 inputs, observed {genomes}")

    model_count = 0
    for genome in genomes:
        for tool in ("gapseq", "CarveMe"):
            for stage in ("draft", "gapfilled"):
                source = model_source(work, genome, tool, stage)
                destination = target / "models" / tool / genome / f"{stage}.xml.gz"
                gzip_copy(source, destination)
                model_count += 1
        gzip_copy(
            work / "inputs/proteins" / f"{genome}.faa",
            target / "inputs/proteins" / f"{genome}.faa.gz",
        )

    for path in sorted((work / "logs").rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(work / "logs")
        if path.suffix == ".log":
            # Empty stderr/stdout logs are valid evidence for successful tools;
            # retain them as deterministic empty gzip streams rather than
            # treating them as missing model payloads.
            gzip_copy(
                path,
                target / "logs" / relative.with_suffix(path.suffix + ".gz"),
                allow_empty=True,
            )
        else:
            copy(path, target / "logs" / relative)

    for name in SCRIPT_FILES:
        copy(root / "scripts" / name, target / "scripts" / name)
    for name in (
        "gapseq.yml",
        "gapseq-linux-64.lock",
        "carveme.yml",
        "carveme-linux-64.lock",
        "drep.yml",
        "drep-linux-64.lock",
    ):
        copy(root / "env" / name, target / "env" / name)
    copy(
        root / "data/small/60-gem-database-manifest.tsv",
        target / "database-manifest.tsv",
    )
    copy(
        database_cache / "gapseq-seqdb-v1.5/database-audit.json",
        target / "database-audit.json",
    )
    copy(
        database_cache / "gapseq-seqdb-v1.5/archives/md5sums.txt",
        target / "database/gapseq-seqdb-v1.5-md5sums.txt",
    )
    for name in (
        "ALLmed.csv",
        "LBmed.csv",
        "MM_glu.csv",
        "autotrophic.csv",
        "MM_anaerobic_CO2_H2.csv",
        "meerwasser.csv",
    ):
        copy(
            gapseq_env / "share/gapseq/dat/media" / name,
            target / "database/gapseq-media" / name,
        )
    copy(
        carveme_root(carveme_env) / "data/input/media_db.tsv",
        target / "database/carveme-media-db.tsv",
    )

    notice = """Article 60 frozen evidence bundle

Eight primary inputs are real dRep representatives reconstructed from the
PRJEB52977 uneven mock study and classified with GTDB-Tk R232. Four additional
inputs are deterministic SGB_002 contig-retention sensitivities generated with
seed 59002. Their shared Prodigal 2.6.3 protein FASTAs are included in
compressed form, while the source MAG FASTAs remain in the checksum-covered
Article 45 representative-genome bundle.

The bundle contains all 48 SBML models: fixed drafts and gap-filled models from
gapseq 2.1.0 and CarveMe 1.6.6 for 12 inputs. Bacterial gapseq models use
ALLmed; three archaeal stress cases use versioned ecological profiles from the
official distribution (autotrophic plus highH2, MM_anaerobic_CO2_H2 plus
highH2, or meerwasser). ALLmed and MM_glu are common post-hoc sensitivities.
CarveMe applies its independent gapfill command to a fixed draft under the
built-in LB definition and audits M9 post hoc. The medium files are included
exactly as used.

The 1.5-GB extracted gapseq sequence database and CarveMe universes/DIAMOND
database are not duplicated. Their immutable source identities, byte counts,
checksums, extracted-file counts and environment locks are retained. Gap-filled
reactions without GPR support are optimization hypotheses, and a positive FBA
objective is constraint feasibility rather than measured growth or flux.
"""
    (target / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    contract = {
        "article": 60,
        "created_from": str(work.relative_to(root)),
        "primary_real_mags": 8,
        "deterministic_sensitivity_genomes": 4,
        "shared_protein_fastas_included": len(genomes),
        "sbml_models_included": model_count,
        "large_databases_included": False,
        "source_mag_fastas_duplicated": False,
        "gapseq_bacterial_gapfill_medium": "ALLmed.csv",
        "gapseq_archaeal_gapfill_profiles": {
            "SGB_008": "autotrophic.csv + highH2",
            "SGB_010": "MM_anaerobic_CO2_H2.csv + highH2",
            "SGB_018": "meerwasser.csv",
        },
        "carveme_gapfill_medium": "LB",
        "truncation_seed": 59002,
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    files = sorted(
        path
        for path in target.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    lines = [f"{sha256(path)}  {path.relative_to(target).as_posix()}" for path in files]
    (target / "file-checksums.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(
        f"Frozen Article 60 bundle: {len(files)} payload files, "
        f"{model_count} compressed SBML models in {target}"
    )


if __name__ == "__main__":
    main()
