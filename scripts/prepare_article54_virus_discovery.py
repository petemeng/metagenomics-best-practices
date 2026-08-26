#!/usr/bin/env python3
"""Download and checksum-gate the compact public inputs for Article 54."""

from __future__ import annotations

import argparse
import shutil
import time
import urllib.request
from pathlib import Path

from article41_44_utils import dump_json, sha256, write_tsv


CHECKV_COMMIT = "6a118f20e895105ce0e4f10257955494c60f1293"
ASSETS = (
    {
        "AssetID": "miuvig-jats",
        "File": "PMC6871006.xml",
        "URL": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6871006/fullTextXML",
        "Bytes": 183_966,
        "SHA256": "289c221ece18f17a0a8930a354cfc3e01b712751abfdd52ff7b054c6f31c36fa",
        "DOI": "10.1038/nbt.4306",
        "License": "CC BY 4.0",
        "Role": "MIUViG definitions and vOTU threshold",
    },
    {
        "AssetID": "checkv-regression-input",
        "File": "checkv-test-sequences.fna",
        "URL": (
            "https://bitbucket.org/berkeleylab/checkv/raw/"
            f"{CHECKV_COMMIT}/test/test_sequences.fna"
        ),
        "Bytes": 910_376,
        "SHA256": "4347d1f5f52d4e2b6986845afb8e401fd146bb5a9628768e0ba125ec380144f5",
        "DOI": "10.1038/s41587-020-00774-7",
        "License": "BSD-3-Clause",
        "Role": "46-sequence official CheckV regression input",
    },
    {
        "AssetID": "checkv-regression-truth",
        "File": "checkv-upstream-ground-truth-quality-summary.tsv",
        "URL": (
            "https://bitbucket.org/berkeleylab/checkv/raw/"
            f"{CHECKV_COMMIT}/test/ground_truth/quality_summary.tsv"
        ),
        "Bytes": 5_518,
        "SHA256": "f0d692ab02446bca007722b0f97a928b4cc893c94bcd311e151196bb9afa5a76",
        "DOI": "10.1038/s41587-020-00774-7",
        "License": "BSD-3-Clause",
        "Role": "official expected CheckV quality summary",
    },
)


def retrieve(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    for attempt in range(1, 6):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "metagenomics-best-practices/article54"}
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                with partial.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
            partial.replace(target)
            return
        except Exception:
            if attempt == 5:
                raise
            time.sleep(2**attempt)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    raw = root / "data/raw/article54"
    inputs = work / "input"
    inputs.mkdir(parents=True, exist_ok=True)

    audit_rows = []
    manifest_rows = []
    for asset in ASSETS:
        source = raw / str(asset["File"])
        if not source.is_file():
            retrieve(str(asset["URL"]), source)
        observed_bytes = source.stat().st_size
        observed_hash = sha256(source)
        passed = (
            observed_bytes == int(asset["Bytes"])
            and observed_hash == str(asset["SHA256"])
        )
        audit_rows.append(
            {
                "AssetID": asset["AssetID"],
                "File": asset["File"],
                "ExpectedBytes": asset["Bytes"],
                "ObservedBytes": observed_bytes,
                "ExpectedSHA256": asset["SHA256"],
                "ObservedSHA256": observed_hash,
                "ChecksumPass": passed,
            }
        )
        if not passed:
            raise RuntimeError(f"Checksum gate failed for {source}")
        shutil.copy2(source, inputs / str(asset["File"]))
        manifest_rows.append(
            {
                **asset,
                "SourceIdentity": (
                    f"CheckV commit {CHECKV_COMMIT}"
                    if str(asset["AssetID"]).startswith("checkv-")
                    else "Europe PMC PMCID PMC6871006"
                ),
            }
        )

    write_tsv(work / "asset-manifest.tsv", manifest_rows)
    write_tsv(work / "asset-check-audit.tsv", audit_rows)
    write_tsv(
        work / "input-lineage.tsv",
        [
            {
                "Input": "checkv-test-sequences.fna",
                "Records": 46,
                "Origin": "CheckV upstream regression fixture",
                "GeneratedBy": f"CheckV commit {CHECKV_COMMIT}",
                "NextStep": "geNomad, VirSorter2, CheckV, all-vs-all BLASTN",
            }
        ],
    )
    dump_json(
        work / "run-contract.json",
        {
            "article": 54,
            "seed": 20260754,
            "input_records": 46,
            "checkv_commit": CHECKV_COMMIT,
            "votu_min_ani_pct": 95,
            "votu_min_alignment_fraction_shorter_pct": 85,
            "virsorter2_default_score": 0.5,
            "virsorter2_high_confidence_rule": (
                "max_score>=0.9 OR (max_score>=0.7 AND hallmark>=1)"
            ),
            "random_output_requested": False,
        },
    )
    print(f"Article 54 inputs verified: {len(ASSETS)} assets in {inputs}")


if __name__ == "__main__":
    main()
