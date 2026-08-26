#!/usr/bin/env python3
"""Prepare checksum-gated open data for the Article 53 transmission audit."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import zipfile
from pathlib import Path

from article41_44_utils import dump_json, sha256, write_tsv


MOTHER_SUPPLEMENT_URL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/"
    "PMC5264247/supplementaryFiles"
)
MOTHER_XML_URL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5264247/fullTextXML"
)
FMT_WORKBOOK_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1186%2Fs40168-022-01251-w/MediaObjects/"
    "40168_2022_1251_MOESM7_ESM.xlsx"
)
FMT_XML_URL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8951724/fullTextXML"
)

EXPECTED = {
    "sys001172080st8.xlsx": (
        10_992,
        "8aa0eb1e386eb78b983fac927a1ba1ab2099956570688704d3f9a8a6ca7299a3",
    ),
    "sys001172080st9.xlsx": (
        118_131,
        "12b5dc05b3237b65512c3221812dffef9b749287f46af32ca9d983835c0e9144",
    ),
    "PMC5264247.xml": (
        201_829,
        "1891a7c22bcb3616e77dd94757f1ef6be6b2a04bd7271a0deea1c4a87a8d17bf",
    ),
    "40168_2022_1251_MOESM7_ESM.xlsx": (
        168_598,
        "0f51d05906a2a574970070cab4302204c0f04c2a040bce5d923cd921d91ce757",
    ),
    "PMC8951724.xml": (
        141_847,
        "382a06a553b33b4e6edbcdba900f476986794d699eced2d00ad96caa62219411",
    ),
}


def download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "curl",
            "-fL",
            "--retry",
            "8",
            "--retry-all-errors",
            "--connect-timeout",
            "30",
            "--max-time",
            "900",
            "--output",
            str(target),
            url,
        ],
        check=True,
    )


def stage(source: Path | None, url: str, target: Path) -> None:
    if source is None:
        download(url, target)
    else:
        source = source.resolve()
        if not source.is_file():
            raise FileNotFoundError(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def gate(path: Path) -> None:
    expected_bytes, expected_hash = EXPECTED[path.name]
    if path.stat().st_size != expected_bytes or sha256(path) != expected_hash:
        raise ValueError(
            f"Pinned source mismatch for {path.name}: "
            f"observed {path.stat().st_size} bytes / {sha256(path)}"
        )


def extract_tables(archive: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    required = ("sys001172080st8.xlsx", "sys001172080st9.xlsx")
    with zipfile.ZipFile(archive) as handle:
        members = {Path(name).name: name for name in handle.namelist()}
        for name in required:
            if name not in members:
                raise FileNotFoundError(f"Missing {name} in Europe PMC supplement bundle")
            target = output / name
            with handle.open(members[name]) as source, target.open("wb") as sink:
                shutil.copyfileobj(source, sink)
            gate(target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--mother-supplement", type=Path)
    parser.add_argument("--mother-xml", type=Path)
    parser.add_argument("--fmt-workbook", type=Path)
    parser.add_argument("--fmt-xml", type=Path)
    args = parser.parse_args()

    work = args.work_dir.resolve()
    if work.exists():
        shutil.rmtree(work)
    for folder in ("downloads", "input", "logs", "summary"):
        (work / folder).mkdir(parents=True)

    mother_archive = work / "downloads/PMC5264247_SupplementaryFiles.zip"
    stage(args.mother_supplement, MOTHER_SUPPLEMENT_URL, mother_archive)
    extract_tables(mother_archive, work / "input")

    mother_xml = work / "input/PMC5264247.xml"
    fmt_workbook = work / "input/40168_2022_1251_MOESM7_ESM.xlsx"
    fmt_xml = work / "input/PMC8951724.xml"
    stage(args.mother_xml, MOTHER_XML_URL, mother_xml)
    stage(args.fmt_workbook, FMT_WORKBOOK_URL, fmt_workbook)
    stage(args.fmt_xml, FMT_XML_URL, fmt_xml)
    for path in (mother_xml, fmt_workbook, fmt_xml):
        gate(path)

    assets = [
        {
            "Asset": "Mother-infant sample design (Table S1)",
            "Study": "Asnicar et al. 2017",
            "Identifier": "PMC5264247 / sys001172080st8.xlsx",
            "DOI": "10.1128/mSystems.00164-16",
            "URL": MOTHER_SUPPLEMENT_URL,
            "SourceFile": "input/sys001172080st8.xlsx",
            "Bytes": EXPECTED["sys001172080st8.xlsx"][0],
            "SHA256": EXPECTED["sys001172080st8.xlsx"][1],
            "License": "CC BY 4.0",
            "GateScope": "extracted member; dynamic ZIP envelope is not the identity",
        },
        {
            "Asset": "Mother-infant MetaPhlAn2 profile (Table S2)",
            "Study": "Asnicar et al. 2017",
            "Identifier": "PMC5264247 / sys001172080st9.xlsx",
            "DOI": "10.1128/mSystems.00164-16",
            "URL": MOTHER_SUPPLEMENT_URL,
            "SourceFile": "input/sys001172080st9.xlsx",
            "Bytes": EXPECTED["sys001172080st9.xlsx"][0],
            "SHA256": EXPECTED["sys001172080st9.xlsx"][1],
            "License": "CC BY 4.0",
            "GateScope": "extracted member; dynamic ZIP envelope is not the identity",
        },
        {
            "Asset": "Mother-infant full-text methods and reported strain evidence",
            "Study": "Asnicar et al. 2017",
            "Identifier": "PMC5264247 JATS XML",
            "DOI": "10.1128/mSystems.00164-16",
            "URL": MOTHER_XML_URL,
            "SourceFile": "input/PMC5264247.xml",
            "Bytes": EXPECTED["PMC5264247.xml"][0],
            "SHA256": EXPECTED["PMC5264247.xml"][1],
            "License": "CC BY 4.0",
            "GateScope": "complete file",
        },
        {
            "Asset": "SameStr FMT supplementary workbook (Tables S1-S11)",
            "Study": "Podlesny et al. 2022",
            "Identifier": "PMC8951724 / MOESM7",
            "DOI": "10.1186/s40168-022-01251-w",
            "URL": FMT_WORKBOOK_URL,
            "SourceFile": "input/40168_2022_1251_MOESM7_ESM.xlsx",
            "Bytes": EXPECTED["40168_2022_1251_MOESM7_ESM.xlsx"][0],
            "SHA256": EXPECTED["40168_2022_1251_MOESM7_ESM.xlsx"][1],
            "License": "CC BY 4.0",
            "GateScope": "complete file",
        },
        {
            "Asset": "SameStr full-text methods",
            "Study": "Podlesny et al. 2022",
            "Identifier": "PMC8951724 JATS XML",
            "DOI": "10.1186/s40168-022-01251-w",
            "URL": FMT_XML_URL,
            "SourceFile": "input/PMC8951724.xml",
            "Bytes": EXPECTED["PMC8951724.xml"][0],
            "SHA256": EXPECTED["PMC8951724.xml"][1],
            "License": "CC BY 4.0",
            "GateScope": "complete file",
        },
    ]
    write_tsv(work / "asset-manifest.tsv", assets)
    write_tsv(
        work / "input-lineage.tsv",
        [
            {
                "Output": "mother-infant species-sharing negative-control table",
                "ImmediateInput": "Asnicar Table S2 MetaPhlAn2 percentage profile",
                "Transformation": "named species only; >0.1% presence; time-matched mother comparisons",
                "Evidence": "asset-manifest.tsv and pairwise-sharing output",
            },
            {
                "Output": "mother-infant strain-divergence evidence",
                "ImmediateInput": "Asnicar article JATS XML and cited Figures 2/S3-S5",
                "Transformation": "assert exact reported values; no recomputation from absent marker alignments",
                "Evidence": "published-strain-evidence.tsv",
            },
            {
                "Output": "FMT sharing, relatedness and source summaries",
                "ImmediateInput": "Podlesny Tables S6-S8",
                "Transformation": "lossless field normalization and descriptive aggregation",
                "Evidence": "classifier-performance, casewise-sharing and source-event tables",
            },
        ],
    )
    dump_json(
        work / "run-contract.json",
        {
            "article": 53,
            "seed": 20260753,
            "mother_infant": {
                "species_presence_threshold_percent": 0.1,
                "comparison": "infant versus mothers sampled at the same nominal time point",
                "named_species_only": True,
                "formal_pairwise_test": False,
            },
            "same_str": {
                "mvs_threshold": 0.999,
                "minimum_overlap_bp": 5000,
                "minimum_allele_fraction": 0.10,
                "marker_database": "db_v20 / mpa_v20_m200",
            },
            "random_output_requested": False,
        },
    )
    write_tsv(
        work / "source-provenance.tsv",
        [
            {
                "Study": "Asnicar et al. 2017",
                "PMCID": "PMC5264247",
                "BioProject": "PRJNA339914",
                "DataRole": "mother-infant species profiles and reported strain evidence",
            },
            {
                "Study": "Podlesny et al. 2022",
                "PMCID": "PMC8951724",
                "BioProject": "PRJEB39023 plus cited public cohorts",
                "DataRole": "SameStr validation and FMT engraftment outputs",
            },
        ],
    )
    (work / ".article53-inputs-complete").write_text("complete\n", encoding="utf-8")
    print("Prepared five checksum-gated Article 53 source objects")


if __name__ == "__main__":
    main()

