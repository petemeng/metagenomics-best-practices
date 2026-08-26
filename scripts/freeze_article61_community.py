#!/usr/bin/env python3
"""Create a compact, checksum-covered Article 61 evidence bundle."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
from pathlib import Path


RAW_FILES = (
    "resource-audit.tsv",
    "sample-selection-audit.tsv",
    "selected-samples.tsv",
    "micom-taxonomy.tsv",
    "micom-taxonomy-equal.tsv",
    "taxon-match-audit.tsv",
    "model-coverage.tsv",
    "database-manifest.tsv",
    "medium.tsv",
    "smetana-media.tsv",
    "smetana-medium-compatibility.tsv",
    "smetana-subcommunity.tsv",
    "smetana-model-audit.tsv",
    "smetana-compatibility-audit.tsv",
    "smetana-compatibility-summary.json",
    "smetana-run.log",
    "preparation-contract.json",
    "software-versions.json",
    "micom-observed-manifest.tsv",
    "micom-equal-manifest.tsv",
    "micom-tradeoff-growth.tsv",
    "micom-primary-growth.tsv",
    "micom-primary-exchanges.tsv",
    "micom-exchange-annotations.tsv",
    "micom-medium-sensitivity.tsv",
    "micom-equal-abundance-growth.tsv",
    "smetana-global.tsv",
    "smetana-detailed.tsv",
    "run-ledger.tsv",
)

SUMMARY_FILES = (
    "sample-selection-audit.tsv",
    "selected-samples.tsv",
    "model-coverage.tsv",
    "taxon-match-audit.tsv",
    "tradeoff-summary.tsv",
    "primary-growth.tsv",
    "primary-growth-summary.tsv",
    "medium-sensitivity-summary.tsv",
    "abundance-sensitivity-summary.tsv",
    "primary-exchanges.tsv",
    "net-community-flux.tsv",
    "micom-crossfeeding-potential.tsv",
    "focal-micom-fluxes.tsv",
    "focal-micom-potential-edges.tsv",
    "smetana-global.tsv",
    "smetana-detailed.tsv",
    "smetana-compatibility-audit.tsv",
    "smetana-component-summary.tsv",
    "smetana-pair-summary.tsv",
    "cross-method-concordance.tsv",
    "run-ledger.tsv",
    "analysis-metrics.json",
)

SCRIPT_FILES = (
    "download_article61_community_db.py",
    "prepare_article61_community.py",
    "run_article61_community.py",
    "rerun_article61_smetana.py",
    "audit_article61_smetana_compatibility.py",
    "summarize_article61_community.py",
    "plot_article61_community.py",
    "freeze_article61_community.py",
    "validate_article61_community.py",
    "article42_44_validation_utils.py",
)

ENV_FILES = (
    "community-metabolism.yml",
    "community-metabolism-linux-64.lock",
    "community-metabolism-pip.lock",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--summary-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


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


def gzip_copy(source: Path, target: Path) -> None:
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as source_handle, target.open("wb") as raw_target:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_target, mtime=0) as target_handle:
            shutil.copyfileobj(source_handle, target_handle, length=8 * 1024 * 1024)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    summary = args.summary_dir.resolve()
    target = args.output_dir.resolve()

    for name in RAW_FILES:
        if not (work / name).is_file():
            raise FileNotFoundError(f"Article 61 raw output missing: {name}")
    for name in SUMMARY_FILES:
        if not (summary / name).is_file():
            raise FileNotFoundError(f"Run the Article 61 summarizer first: {name}")

    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)

    for name in RAW_FILES:
        copy(work / name, target / "raw" / name)
    for name in SUMMARY_FILES:
        copy(summary / name, target / name)

    model_audit = (work / "smetana-model-audit.tsv").read_text(encoding="utf-8").splitlines()
    if len(model_audit) != 7:
        raise RuntimeError("Expected six SMETANA model rows")
    model_count = 0
    for source in sorted((work / "smetana/models").glob("*.xml")):
        gzip_copy(source, target / "models/smetana" / f"{source.name}.gz")
        model_count += 1
    if model_count != 6:
        raise RuntimeError(f"Expected six SMETANA SBML files, observed {model_count}")

    for name in SCRIPT_FILES:
        copy(root / "scripts" / name, target / "scripts" / name)
    for name in ENV_FILES:
        copy(root / "env" / name, target / "env" / name)
    copy(
        root / "data/small/61-community-database-manifest.tsv",
        target / "database-resource-manifest.tsv",
    )
    copy(
        root / "data/small/22-diversity-inputs/human-sample-metadata.tsv",
        target / "inputs/human-sample-metadata.tsv",
    )
    copy(
        root / "data/small/22-diversity-inputs/species-relative-abundance.tsv.gz",
        target / "inputs/species-relative-abundance.tsv.gz",
    )
    copy(
        work / "database/agora201/.identity.json",
        target / "database-extraction-identity.json",
    )

    notice = """Article 61 frozen evidence bundle

The biological inputs are the checksum-locked AsnicarF_2017 species profile and
metadata distributed with curatedMetagenomicData 3.12.0. Six independent adult
subjects pass the prespecified read-depth and visit-selection gates. Abundance
coverage is always reported against the original whole-profile denominator.

The large 264-MB AGORA 2.01 QIIME artifact, its approximately 2-GB extraction,
and MICOM community pickle files are not duplicated. Their exact release,
checksums, manifest, model counts, taxonomy and extraction UUID are retained.
The six uniquely named SBML files used by SMETANA are compressed in this bundle.

MICOM pFBA fluxes are constraint-dependent metabolic hypotheses. SMETANA 1.2.1
returned 7,853 detailed component rows, but global MIP/MRO was not estimable:
the legacy non-interacting merge had zero complete-medium growth for all six
exported AGORA2 SBML members, and the MRO member solve failed. Zero composite
SMETANA scores are therefore not evidence that biological exchange is absent.
None of these outputs are measured in vivo fluxes, metabolite concentrations
or proof of mutualism.
"""
    (target / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    metrics = json.loads((summary / "analysis-metrics.json").read_text(encoding="utf-8"))
    contract = {
        "article": 61,
        "created_from": str(work.relative_to(root)),
        "selected_samples": metrics["selected_samples"],
        "independent_subjects": metrics["independent_subjects"],
        "smetana_subcommunity_size": metrics["smetana_subcommunity_size"],
        "smetana_sbml_models_included": model_count,
        "smetana_global_estimable": metrics["smetana_global_estimable"],
        "smetana_cross_method_comparison_estimable": metrics[
            "smetana_cross_method_comparison_estimable"
        ],
        "agora_qza_included": False,
        "micom_pickles_included": False,
        "source_profile_included": True,
        "seed": 61001,
        "primary_tradeoff": 0.5,
        "relative_abundance_cutoff": 0.001,
        "minimum_reads": 1_000_000,
        "checksum_manifest": "file-checksums.sha256",
    }
    (target / "frozen-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    payloads = sorted(
        path for path in target.rglob("*")
        if path.is_file() and path.name != "file-checksums.sha256"
    )
    lines = [f"{sha256(path)}  {path.relative_to(target).as_posix()}" for path in payloads]
    (target / "file-checksums.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(
        f"Frozen Article 61 bundle: {len(payloads)} payload files, "
        f"{model_count} compressed SMETANA models in {target}"
    )


if __name__ == "__main__":
    main()
