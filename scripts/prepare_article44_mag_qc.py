#!/usr/bin/env python3
"""Reconstruct the Article 43 selected MAG candidates for the Article 44 QC audit."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from article41_44_utils import dump_json, fasta_records, fasta_summary, read_tsv, sha256, write_tsv


CHECKM2_DB_SHA256 = "1b86ef3eac0813c1853f53182c17657045e3763d66f384ec95747261a63ae46f"
GUNC_DB_SHA256 = "2dabe83f2ab7f0b38e78cfdbd8ca33bdc578b330d7501cb42457d331bc8c09d4"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--article43-frozen", type=Path, required=True)
    parser.add_argument("--assembly-env", type=Path, required=True)
    parser.add_argument("--mag-qc-env", type=Path, required=True)
    parser.add_argument("--gunc-env", type=Path, required=True)
    parser.add_argument("--checkm1-env", type=Path, required=True)
    parser.add_argument("--checkm2-db", type=Path, required=True)
    parser.add_argument("--gunc-db", type=Path, required=True)
    return parser.parse_args()


def output(command: list[str], env: dict[str, str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    if result.returncode:
        raise RuntimeError(f"Command failed: {command}\n{result.stdout}\n{result.stderr}")
    return result.stdout + result.stderr


def conda_package_version(prefix: Path, package: str) -> str:
    matches = sorted((prefix / "conda-meta").glob(f"{package}-*.json"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {package} conda record in {prefix}: {matches}")
    return str(json.loads(matches[0].read_text(encoding="utf-8"))["version"])


def main() -> int:
    args = parse_args()
    root, work, frozen43 = args.project_root.resolve(), args.work_dir.resolve(), args.article43_frozen.resolve()
    assembly_env, mag_qc_env = args.assembly_env.resolve(), args.mag_qc_env.resolve()
    gunc_env, checkm1_env = args.gunc_env.resolve(), args.checkm1_env.resolve()
    checkm2_db, gunc_db = args.checkm2_db.resolve(), args.gunc_db.resolve()
    if (work / ".article44-inputs-complete").is_file():
        print(f"Article 44 inputs already prepared: {work}")
        return 0
    if work.exists() and any(path.is_file() for path in work.rglob("*")):
        raise SystemExit(f"Refusing incomplete work directory containing files: {work}")
    required = (
        "selected-mag-candidates.tsv", "selected-refinement-membership.tsv.gz",
        "final-method-selection.tsv", "file-checksums.sha256", "run-summary.json",
    )
    for name in required:
        if not (frozen43 / name).is_file():
            raise FileNotFoundError(frozen43 / name)
    observed_checkm2_sha = sha256(checkm2_db) if checkm2_db.is_file() else ""
    observed_gunc_sha = sha256(gunc_db) if gunc_db.is_file() else ""
    if observed_checkm2_sha != CHECKM2_DB_SHA256:
        raise RuntimeError("CheckM2 version-3 database is missing or has checksum drift")
    if observed_gunc_sha != GUNC_DB_SHA256:
        raise RuntimeError("GUNC ProGenomes 2.1 database is missing or has checksum drift")
    for directory in ("inputs", "bins", "qc", "features", "graph", "logs", "summary", "tmp"):
        (work / directory).mkdir(parents=True, exist_ok=True)

    runtime = os.environ.copy()
    for key in ("PERL5LIB", "PERL_LOCAL_LIB_ROOT", "PERL_MB_OPT", "PERL_MM_OPT"):
        runtime.pop(key, None)
    runtime.update({"LC_ALL": "C", "TZ": "UTC", "PYTHONHASHSEED": "0", "PYTHONNOUSERSITE": "1"})
    tools = {
        "CheckM2": mag_qc_env / "bin/checkm2",
        "GUNC": gunc_env / "bin/gunc",
        "CheckM1": checkm1_env / "bin/checkm",
        "barrnap": mag_qc_env / "bin/barrnap",
        "tRNAscan-SE": mag_qc_env / "bin/tRNAscan-SE",
        "Prodigal": mag_qc_env / "bin/prodigal",
        "MEGAHIT toolkit": assembly_env / "bin/megahit_toolkit",
    }
    versions = {
        "CheckM2": output([str(tools["CheckM2"]), "--version"], {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"}).strip(),
        "GUNC": output([str(tools["GUNC"]), "--version"], {**runtime, "PATH": f"{gunc_env / 'bin'}:/usr/bin:/bin"}).strip(),
        "CheckM1": conda_package_version(checkm1_env, "checkm-genome"),
        "barrnap": re.search(r"barrnap\s+([0-9.]+)", output([str(tools["barrnap"]), "--version"], {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"})).group(1),
        "tRNAscan-SE": re.search(r"tRNAscan-SE\s+([0-9.]+)", output([str(tools["tRNAscan-SE"]), "--help"], {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"})).group(1),
        "Prodigal": re.search(r"Prodigal V([0-9.]+)", output([str(tools["Prodigal"]), "-v"], {**runtime, "PATH": f"{mag_qc_env / 'bin'}:/usr/bin:/bin"})).group(1),
        "MEGAHIT toolkit": output([str(tools["MEGAHIT toolkit"]), "dumpversion"], {**runtime, "PATH": f"{assembly_env / 'bin'}:/usr/bin:/bin"}).strip().removeprefix("v"),
    }
    expected = {
        "CheckM2": "1.1.0", "GUNC": "1.1.0", "CheckM1": "1.2.5", "barrnap": "1.10.5",
        "tRNAscan-SE": "2.0.13", "Prodigal": "2.6.3", "MEGAHIT toolkit": "1.2.9",
    }
    if versions != expected:
        raise RuntimeError(f"Article 44 version drift: {versions}")
    write_tsv(work / "tool-versions.tsv", [
        {"Tool": name, "Version": versions[name], "Executable": str(tools[name])} for name in tools
    ])

    source = root / "data/small/30-short-read-assembly-frozen/contigs/megahit-coassembly.ge1000.fna.gz"
    full_sequences = dict(fasta_records(source))
    common_sequences = {name: sequence for name, sequence in full_sequences.items() if len(sequence) >= 1500}
    if len(common_sequences) != 10203:
        raise RuntimeError(f"Article 44 common-coordinate drift: {len(common_sequences)} contigs")
    common = work / "inputs/megahit-coassembly.ge1500.fna"
    with common.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in common_sequences.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    coordinate_summary, _ = fasta_summary(common)

    selected = read_tsv(frozen43 / "selected-mag-candidates.tsv")
    membership = read_tsv(frozen43 / "selected-refinement-membership.tsv.gz")
    selected_ids = {row["RefinedID"] for row in selected}
    if {row["RefinedID"] for row in membership} != selected_ids:
        raise RuntimeError("Article 43 selected-candidate membership mismatch")
    by_bin: dict[str, list[str]] = defaultdict(list)
    for row in membership:
        if row["Contig"] not in common_sequences:
            raise RuntimeError(f"Unknown selected contig: {row['Contig']}")
        by_bin[row["RefinedID"]].append(row["Contig"])
    expected_sha = {row["RefinedID"]: row["RefinedSHA256"] for row in selected}
    reconstruction = []
    seen = set()
    for refined_id, names in sorted(by_bin.items()):
        if seen & set(names):
            raise RuntimeError("Selected Article 43 bins overlap")
        seen.update(names)
        target = work / "bins" / f"{refined_id}.fna"
        with target.open("w", encoding="utf-8", newline="\n") as handle:
            for name in sorted(names):
                sequence = common_sequences[name]
                handle.write(f">{name}\n")
                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start : start + 80] + "\n")
        observed = sha256(target)
        if observed != expected_sha[refined_id]:
            raise RuntimeError(f"Selected MAG reconstruction checksum drift: {refined_id}")
        reconstruction.append({"MAG": refined_id, "Contigs": len(names), "BinBp": sum(len(common_sequences[name]) for name in names), "Bytes": target.stat().st_size, "SHA256": observed, "Status": "PASS"})
    write_tsv(work / "selected-mag-reconstruction-audit.tsv", reconstruction)

    k141 = root / "data/raw/article30/work/assemblies/megahit-coassembly/intermediate_contigs/k141.contigs.fa"
    paired = root / "data/small/41-read-mapping-depth-frozen/raw/paired-contigs.tsv"
    checkm_manifest = checkm1_env / "checkm_data/.dmanifest"
    for path in (k141, paired, checkm_manifest):
        if not path.is_file():
            raise FileNotFoundError(path)
    write_tsv(work / "input-audit.tsv", [
        {"Role": "Article43-frozen-manifest", "Path": str((frozen43 / "file-checksums.sha256").relative_to(root)), "Bytes": (frozen43 / "file-checksums.sha256").stat().st_size, "SHA256": sha256(frozen43 / "file-checksums.sha256"), "Status": "PASS"},
        {"Role": "Article30-coassembly-ge1000", "Path": str(source.relative_to(root)), "Bytes": source.stat().st_size, "SHA256": sha256(source), "Status": "PASS"},
        {"Role": "MEGAHIT-k141-graph-source", "Path": str(k141.relative_to(root)), "Bytes": k141.stat().st_size, "SHA256": sha256(k141), "Status": "PASS"},
        {"Role": "Article41-read-pair-links", "Path": str(paired.relative_to(root)), "Bytes": paired.stat().st_size, "SHA256": sha256(paired), "Status": "PASS"},
        {"Role": "CheckM1-2015-dmanifest", "Path": str(checkm_manifest), "Bytes": checkm_manifest.stat().st_size, "SHA256": sha256(checkm_manifest), "Status": "PASS"},
        {"Role": "CheckM2-database-v3", "Path": str(checkm2_db), "Bytes": checkm2_db.stat().st_size, "SHA256": observed_checkm2_sha, "Status": "PASS"},
        {"Role": "GUNC-ProGenomes2.1", "Path": str(gunc_db), "Bytes": gunc_db.stat().st_size, "SHA256": observed_gunc_sha, "Status": "PASS"},
    ])
    write_tsv(work / "input-lineage.tsv", [
        {"Output": "Selected MAG candidate FASTAs", "ImmediateInput": "Article 43 reference-free selected membership + Article 30 common assembly", "Transformation": "checksum-exact reconstruction", "TruthUsed": "No", "Evidence": "selected-mag-reconstruction-audit.tsv"},
        {"Output": "MIMAG quality tier", "ImmediateInput": "CheckM2, GUNC, barrnap, tRNAscan-SE outputs", "Transformation": ">90/<5 + GUNC + complete 5S/16S/23S + >=18 tRNA isotypes for HQ; >=50/<10 + GUNC for MQ", "TruthUsed": "No", "Evidence": "summary/mag-quality-summary.tsv"},
        {"Output": "Assembly graph audit", "ImmediateInput": "MEGAHIT k141 graph + Article 41 paired-contig links", "Transformation": "within-bin components and boundary-link summaries", "TruthUsed": "No", "Evidence": "summary/assembly-graph-audit.tsv"},
    ])
    dump_json(work / "run-contract.json", {
        "article": 44, "seed": 20260744, "coordinate_set": coordinate_summary,
        "selected_bins": len(selected_ids),
        "mimag_high": {"completeness": ">90", "contamination": "<5", "gunc_pass": True, "complete_rrna": ["5S", "16S", "23S"], "minimum_trna_isotypes": 18},
        "mimag_medium": {"completeness": ">=50", "contamination": "<10", "gunc_pass": True},
        "checkm1_role": "supplementary lineage and strain-heterogeneity audit only",
        "graph_role": "triage evidence; not part of the MIMAG tier and not proof of circularity",
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
    })
    (work / ".article44-inputs-complete").write_text("PASS\n", encoding="utf-8")
    print(json.dumps({"work": str(work), "selected_bins": len(selected_ids), "coordinate_set": coordinate_summary}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
