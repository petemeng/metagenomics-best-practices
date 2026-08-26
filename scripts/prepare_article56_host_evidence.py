#!/usr/bin/env python3
"""Build source-locked evidence tables for Article 56.

The numerical results in this chapter are extracted from the primary iPHoP
and MIUViG full-text XML files.  No missing supplementary spreadsheet value is
reconstructed or simulated.
"""

from __future__ import annotations

import argparse
import platform
import shutil
import sys
import time
import urllib.request
from pathlib import Path

from lxml import etree

from article41_44_utils import dump_json, sha256, write_tsv


SEED = 20260756
ASSETS = (
    {
        "AssetID": "iphop-primary-xml",
        "File": "PMC10155999.xml",
        "Source": "data/raw/article56/PMC10155999.xml",
        "URL": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10155999/fullTextXML",
        "Bytes": 260_284,
        "SHA256": "c2194569972cc1572a0d709e58916c69fd7f1603c7f4645098194d641319e378",
        "DOI": "10.1371/journal.pbio.3002083",
        "License": "CC BY 4.0",
        "Role": "iPHoP benchmark, calibration, methods, and limitations",
    },
    {
        "AssetID": "miuvig-primary-xml",
        "File": "PMC6871006.xml",
        "Source": "data/raw/article54/PMC6871006.xml",
        "URL": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6871006/fullTextXML",
        "Bytes": 183_966,
        "SHA256": "289c221ece18f17a0a8930a354cfc3e01b712751abfdd52ff7b054c6f31c36fa",
        "DOI": "10.1038/nbt.4306",
        "License": "CC BY 4.0",
        "Role": "evidence classes and claim limitations for uncultivated viruses",
    },
)


def retrieve(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    for attempt in range(1, 6):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "metagenomics-best-practices/article56"},
            )
            with urllib.request.urlopen(request, timeout=180) as response:
                with partial.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
            partial.replace(target)
            return
        except Exception:
            partial.unlink(missing_ok=True)
            if attempt == 5:
                raise
            time.sleep(2**attempt)


def normalized_text(element: etree._Element) -> str:
    return " ".join("".join(element.itertext()).split())


def paragraphs(tree: etree._ElementTree, *, main_body_only: bool) -> list[str]:
    if main_body_only:
        nodes = tree.xpath(
            "/*[local-name()='article']/*[local-name()='body']"
            "//*[local-name()='p']"
        )
    else:
        nodes = tree.xpath("//*[local-name()='p']")
    return [normalized_text(node) for node in nodes]


def locate(paragraph_text: list[str], needle: str) -> tuple[int, str]:
    hits = [(index, text) for index, text in enumerate(paragraph_text, 1) if needle in text]
    if len(hits) > 1 and len({text for _, text in hits}) == 1:
        # Europe PMC may repeat an identical supplementary caption in the
        # manuscript and review-history subarticle.  Treat byte-identical
        # paragraph text as one assertion, retaining the first coordinate.
        return hits[0]
    if len(hits) != 1:
        raise RuntimeError(
            f"Expected one paragraph containing {needle!r}; observed {len(hits)}"
        )
    return hits[0]


def mirror_tsv(work: Path, name: str, rows: list[dict[str, object]]) -> None:
    write_tsv(work / name, rows)
    write_tsv(work / "summary" / name, rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    work = args.work_dir.resolve()
    inputs = work / "input"
    summary_dir = work / "summary"
    logs = work / "logs"
    for directory in (inputs, summary_dir, logs):
        directory.mkdir(parents=True, exist_ok=True)

    asset_manifest: list[dict[str, object]] = []
    asset_audit: list[dict[str, object]] = []
    input_paths: dict[str, Path] = {}
    for asset in ASSETS:
        source = root / str(asset["Source"])
        if not source.is_file():
            retrieve(str(asset["URL"]), source)
        observed_bytes = source.stat().st_size
        observed_hash = sha256(source)
        passed = (
            observed_bytes == int(asset["Bytes"])
            and observed_hash == str(asset["SHA256"])
        )
        asset_audit.append(
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
        target = inputs / str(asset["File"])
        shutil.copy2(source, target)
        input_paths[str(asset["AssetID"])] = target
        asset_manifest.append({key: value for key, value in asset.items() if key != "Source"})
    mirror_tsv(work, "asset-check-audit.tsv", asset_audit)
    mirror_tsv(work, "asset-manifest.tsv", asset_manifest)

    iphop_tree = etree.parse(str(input_paths["iphop-primary-xml"]))
    miuvig_tree = etree.parse(str(input_paths["miuvig-primary-xml"]))
    iphop_main = paragraphs(iphop_tree, main_body_only=True)
    iphop_all = paragraphs(iphop_tree, main_body_only=False)
    miuvig_main = paragraphs(miuvig_tree, main_body_only=True)

    assertion_specs = (
        (
            "IPHOP-TEST-SCOPE",
            "iPHoP",
            iphop_main,
            "This test dataset contained 1,870 genomes, spanning 170 host genera",
            "Benchmark scope",
            "1870 viruses; 170 host genera",
        ),
        (
            "IPHOP-ALIGNMENT-PRECISION",
            "iPHoP",
            iphop_main,
            "alignment-based tools, both phage-based and host-based, were able to reach high (>80%) PPV",
            "Method performance",
            "Score-filtered alignment evidence exceeded 80% PPV",
        ),
        (
            "IPHOP-RAFAH",
            "iPHoP",
            iphop_main,
            "RaFAH [25] in particular able to maintain a very low FDR (<5%)",
            "Method performance",
            "RaFAH FDR below 5% in the benchmark, with reference bias",
        ),
        (
            "IPHOP-SCORE",
            "iPHoP",
            iphop_main,
            "By default, only predictions with a confidence score ≥90",
            "Calibration",
            "Default score >=90; estimated FDR <10%",
        ),
        (
            "IPHOP-IMGVR",
            "iPHoP",
            iphop_main,
            "216,015 high-quality genomes from the IMG/VR v3 database",
            "External application",
            "216015 high-quality IMG/VR genomes",
        ),
        (
            "IPHOP-NOVEL-HOST",
            "iPHoP",
            iphop_main,
            "restricted to “high” score (i.e., score ≥90) and interpreted at the family rank",
            "Claim ceiling",
            "When the host genus is absent, use >=90 and interpret at family rank",
        ),
        (
            "IPHOP-EUK-CONTROL",
            "iPHoP",
            iphop_main,
            "when applied to 8,128 eukaryotic virus genomes from RefSeq",
            "Negative control",
            "8128 eukaryotic viruses; 1018 erroneous prokaryotic-host calls",
        ),
        (
            "IPHOP-RUNTIME",
            "iPHoP",
            iphop_main,
            "approximately 12 minutes for a test set of 5 complete phage genomes",
            "Resource observation",
            "12 minutes; 5 genomes; 6 CPUs; Sept_2021_pub",
        ),
        (
            "IPHOP-LIFESTYLE-STRATA",
            "iPHoP",
            iphop_all,
            "predicted as temperate, either via BacPhlip or based on the genome annotation (n = 949)",
            "Benchmark composition",
            "949 temperate; 663 virulent; 258 not in either displayed stratum",
        ),
        (
            "MIUVIG-HOST-EVIDENCE",
            "MIUViG",
            miuvig_main,
            "These sequence similarities can range from short exact matches",
            "Evidence hierarchy",
            "CRISPR/long matches are more reliable than composition and abundance profiles",
        ),
    )
    source_assertions: list[dict[str, object]] = []
    for assertion_id, source, corpus, needle, evidence_type, interpretation in assertion_specs:
        index, text = locate(corpus, needle)
        source_assertions.append(
            {
                "AssertionID": assertion_id,
                "Source": source,
                "ParagraphIndex": index,
                "EvidenceType": evidence_type,
                "ValidationNeedle": needle,
                "ParagraphSHA256": __import__("hashlib").sha256(text.encode()).hexdigest(),
                "Interpretation": interpretation,
                "SourceCheckPass": True,
            }
        )
    mirror_tsv(work, "source-assertions.tsv", source_assertions)

    benchmark_scope = [
        {
            "Dataset": "Held-out GenBank test set",
            "Category": "All test viruses",
            "Count": 1870,
            "Unit": "viral genomes",
            "ValueStatus": "reported",
            "SourceAssertion": "IPHOP-TEST-SCOPE",
        },
        {
            "Dataset": "Held-out GenBank test set",
            "Category": "Host genera",
            "Count": 170,
            "Unit": "host genera",
            "ValueStatus": "reported",
            "SourceAssertion": "IPHOP-TEST-SCOPE",
        },
        {
            "Dataset": "Held-out GenBank test set",
            "Category": "Temperate stratum",
            "Count": 949,
            "Unit": "viral genomes",
            "ValueStatus": "reported",
            "SourceAssertion": "IPHOP-LIFESTYLE-STRATA",
        },
        {
            "Dataset": "Held-out GenBank test set",
            "Category": "Virulent stratum",
            "Count": 663,
            "Unit": "viral genomes",
            "ValueStatus": "reported",
            "SourceAssertion": "IPHOP-LIFESTYLE-STRATA",
        },
        {
            "Dataset": "Held-out GenBank test set",
            "Category": "Not in either displayed lifestyle stratum",
            "Count": 258,
            "Unit": "viral genomes",
            "ValueStatus": "derived: 1870-949-663",
            "SourceAssertion": "IPHOP-LIFESTYLE-STRATA",
        },
        {
            "Dataset": "IMG/VR v3 application",
            "Category": "High-quality prokaryotic-virus genomes",
            "Count": 216015,
            "Unit": "viral genomes",
            "ValueStatus": "reported",
            "SourceAssertion": "IPHOP-IMGVR",
        },
        {
            "Dataset": "RefSeq r214 negative control",
            "Category": "Eukaryotic-virus genomes",
            "Count": 8128,
            "Unit": "viral genomes",
            "ValueStatus": "reported",
            "SourceAssertion": "IPHOP-EUK-CONTROL",
        },
    ]
    mirror_tsv(work, "benchmark-scope.tsv", benchmark_scope)

    evidence_hierarchy = [
        {
            "Rank": 1,
            "Evidence": "Integrated prophage + host flanks",
            "Signal": "Virus sequence bounded by cellular host sequence",
            "Directness": 5.0,
            "iPHoPComponent": "BLAST may recover this signal",
            "ClaimCeiling": "Host lineage of the linked cellular sequence",
            "MainFalsePositive": "Viral contig merely co-binned; chimeric assembly; contaminated MAG",
            "RequiredAudit": "Inspect both junctions, host markers, coverage, and assembly graph",
        },
        {
            "Rank": 2,
            "Evidence": "CRISPR spacer match",
            "Signal": "Host immune-memory spacer matches viral sequence",
            "Directness": 4.2,
            "iPHoPComponent": "CRISPR",
            "ClaimCeiling": "Host genus when spacer ownership and mismatch rule are credible",
            "MainFalsePositive": "Loose mismatch rule; short/low-complexity spacer; array misassignment",
            "RequiredAudit": "Report spacer length, mismatches, complexity, and host-contig identity",
        },
        {
            "Rank": 3,
            "Evidence": "Long-read / Hi-C linkage",
            "Signal": "Single-molecule span or proximity-ligation contacts",
            "Directness": 3.6,
            "iPHoPComponent": "Not integrated in iPHoP v1.0",
            "ClaimCeiling": "Physical linkage to a host genome under a calibrated contact model",
            "MainFalsePositive": "Index hopping, crosslink/contact bias, abundance-driven contacts",
            "RequiredAudit": "Use negative controls, replicate support, contact normalization, and MAPQ gates",
        },
        {
            "Rank": 4,
            "Evidence": "Reference-phage protein content",
            "Signal": "Related phages with known hosts share protein content",
            "Directness": 3.0,
            "iPHoPComponent": "RaFAH",
            "ClaimCeiling": "Predicted host genus; strongest for reference-like viruses",
            "MainFalsePositive": "Reference-host annotation error and novelty/reference bias",
            "RequiredAudit": "Report AAI to closest viral reference and main prediction method",
        },
        {
            "Rank": 5,
            "Evidence": "k-mer / composition",
            "Signal": "Long-term nucleotide-usage adaptation to a host lineage",
            "Directness": 2.0,
            "iPHoPComponent": "WIsH, VHM, PHP, combined-hosts RF",
            "ClaimCeiling": "Predicted host genus or family after score calibration",
            "MainFalsePositive": "Shared ecology/composition; short contigs; eukaryotic-virus leakage",
            "RequiredAudit": "Apply length/QC gates, host-domain screen, calibrated score, and negatives",
        },
        {
            "Rank": 6,
            "Evidence": "Co-abundance",
            "Signal": "Virus and candidate host covary across samples or time",
            "Directness": 1.0,
            "iPHoPComponent": "Not integrated in iPHoP v1.0",
            "ClaimCeiling": "Candidate ecological association only",
            "MainFalsePositive": "Shared environmental response, compositionality, repeated measures",
            "RequiredAudit": "Use enough independent samples, covariates, FDR, and lag/sensitivity analyses",
        },
    ]
    mirror_tsv(work, "evidence-hierarchy.tsv", evidence_hierarchy)

    confidence_contract = [
        {
            "MinimumScore": 75,
            "NominalMaximumFDRPct": 25,
            "Use": "Exploratory candidate list",
            "Report": "Hypothesis only; require orthogonal evidence",
            "Caveat": "Lowest selectable score; calibration is dataset dependent",
        },
        {
            "MinimumScore": 90,
            "NominalMaximumFDRPct": 10,
            "Use": "Default high-confidence screen",
            "Report": "Predicted host genus; inspect main method",
            "Caveat": "Approximate FDR, not a per-link experimental probability",
        },
        {
            "MinimumScore": 95,
            "NominalMaximumFDRPct": 5,
            "Use": "Stringent sensitivity analysis",
            "Report": "Very-high-confidence prediction with lower recall",
            "Caveat": "Still database- and benchmark-dependent",
        },
    ]
    mirror_tsv(work, "confidence-contract.tsv", confidence_contract)

    false_rate = 100 * 1018 / 8128
    negative_control = [
        {
            "Metric": "Eukaryotic viruses tested",
            "Count": 8128,
            "Percent": 100.0,
            "Denominator": "all eukaryotic-virus controls",
            "ValueStatus": "reported",
        },
        {
            "Metric": "Erroneous prokaryotic-host predictions",
            "Count": 1018,
            "Percent": round(false_rate, 6),
            "Denominator": "all eukaryotic-virus controls",
            "ValueStatus": "count reported; percent derived",
        },
        {
            "Metric": "Errors originating from k-mer comparison",
            "Count": "",
            "Percent": 85.0,
            "Denominator": "erroneous predictions",
            "ValueStatus": "reported percentage",
        },
        {
            "Metric": "Errors with iPHoP score below 90",
            "Count": "",
            "Percent": 90.0,
            "Denominator": "erroneous predictions",
            "ValueStatus": "reported percentage",
        },
        {
            "Metric": "Riboviria among errors",
            "Count": 640,
            "Percent": round(100 * 640 / 1018, 6),
            "Denominator": "erroneous predictions",
            "ValueStatus": "count reported; percent derived",
        },
        {
            "Metric": "Monodnaviria among errors",
            "Count": 155,
            "Percent": round(100 * 155 / 1018, 6),
            "Denominator": "erroneous predictions",
            "ValueStatus": "count reported; percent derived",
        },
    ]
    mirror_tsv(work, "negative-control.tsv", negative_control)

    integration = [
        {
            "Component": "BLAST",
            "EvidenceFamily": "Virus-host sequence alignment",
            "PublishedGate": "identity >=80%; hit length >=500 nt; e-value <=1e-3",
            "IncludedInComposite": True,
            "Interpretation": "Can capture integrated prophage or other shared sequence",
        },
        {
            "Component": "CRISPR",
            "EvidenceFamily": "Spacer-virus alignment",
            "PublishedGate": "spacer >=25 nt; <8 mismatches; complexity <0.6",
            "IncludedInComposite": True,
            "Interpretation": "Fig. 1 strict benchmark used <=2 mismatches",
        },
        {
            "Component": "Combined-hosts-RF",
            "EvidenceFamily": "Ten calibrated host-based classifiers",
            "PublishedGate": "empirical classifier score",
            "IncludedInComposite": True,
            "Interpretation": "Includes alignment-free WIsH, VHM, and PHP signals",
        },
        {
            "Component": "RaFAH",
            "EvidenceFamily": "Reference-phage protein content",
            "PublishedGate": "empirical PPV from RaFAH score",
            "IncludedInComposite": True,
            "Interpretation": "High precision but biased toward reference-like viruses",
        },
        {
            "Component": "Long-read / Hi-C",
            "EvidenceFamily": "Physical linkage",
            "PublishedGate": "not applicable",
            "IncludedInComposite": False,
            "Interpretation": "Must be joined as an orthogonal evidence channel",
        },
        {
            "Component": "Co-abundance",
            "EvidenceFamily": "Cross-sample ecological association",
            "PublishedGate": "not applicable",
            "IncludedInComposite": False,
            "Interpretation": "Must be joined separately and labelled predictive",
        },
    ]
    mirror_tsv(work, "iphop-component-ledger.tsv", integration)

    statuses = {
        "Allowed": 3,
        "Conditional": 2,
        "Avoid": 1,
    }
    claim_matrix = {
        "Integrated prophage + host flanks": [
            "Allowed", "Allowed", "Allowed", "Conditional", "Avoid", "Avoid"
        ],
        "CRISPR spacer match": [
            "Allowed", "Allowed", "Conditional", "Avoid", "Avoid", "Avoid"
        ],
        "Long-read / Hi-C linkage": [
            "Allowed", "Conditional", "Conditional", "Conditional", "Conditional", "Avoid"
        ],
        "iPHoP score >=90": [
            "Conditional", "Conditional", "Conditional", "Avoid", "Avoid", "Avoid"
        ],
        "k-mer / composition alone": [
            "Conditional", "Conditional", "Conditional", "Avoid", "Avoid", "Avoid"
        ],
        "Co-abundance alone": [
            "Avoid", "Conditional", "Conditional", "Avoid", "Avoid", "Avoid"
        ],
    }
    claims = [
        "Host domain",
        "Host family",
        "Host genus",
        "Host species/strain",
        "Active infection now",
        "Causal ecological effect",
    ]
    claim_ceiling = []
    for evidence, labels in claim_matrix.items():
        for claim, label in zip(claims, labels):
            claim_ceiling.append(
                {
                    "Evidence": evidence,
                    "Claim": claim,
                    "Status": label,
                    "StatusCode": statuses[label],
                }
            )
    mirror_tsv(work, "claim-ceiling.tsv", claim_ceiling)

    resource_contract = [
        {
            "Mode": "Frozen evidence and figure regeneration",
            "CPU": 4,
            "RAMGiB": 4,
            "DiskGiB": 0.1,
            "Elapsed": "about 1 minute",
            "MeasurementStatus": "workflow planning bound",
        },
        {
            "Mode": "Published iPHoP benchmark example",
            "CPU": 6,
            "RAMGiB": "not reported",
            "DiskGiB": "not reported",
            "Elapsed": "about 12 minutes for 5 complete phage genomes",
            "MeasurementStatus": "reported by Roux et al.; Sept_2021_pub",
        },
        {
            "Mode": "Cohort-scale current run",
            "CPU": "16-32",
            "RAMGiB": "32-64",
            "DiskGiB": ">=400 including database and work files",
            "Elapsed": "hours to days, input dependent",
            "MeasurementStatus": "conservative planning recommendation, not a benchmark",
        },
    ]
    mirror_tsv(work, "resource-contract.tsv", resource_contract)

    tool_versions = [
        {
            "Tool": "iPHoP",
            "Version": "1.0",
            "Role": "published IMG/VR application",
            "Evidence": "primary article Methods",
        },
        {
            "Tool": "iPHoP",
            "Version": "1.4.2",
            "Role": "current reproducibility environment",
            "Evidence": "Bioconda recipe checked 2026-07-30",
        },
        {
            "Tool": "iPHoP database",
            "Version": "Sept_2021_pub / iPHoP_db_Sept21",
            "Role": "published benchmark database",
            "Evidence": "GTDB r202; IMG 2021-07-07; GEM catalog",
        },
        {
            "Tool": "lxml",
            "Version": etree.LXML_VERSION,
            "Role": "source XML parsing",
            "Evidence": "local runtime",
        },
        {
            "Tool": "Python",
            "Version": platform.python_version(),
            "Role": "deterministic table construction",
            "Evidence": "local runtime",
        },
    ]
    mirror_tsv(work, "tool-versions.tsv", tool_versions)

    determinism = [
        {
            "Step": "Source extraction",
            "Random": False,
            "Seed": SEED,
            "Control": "exact XML SHA-256 plus unique paragraph needles",
        },
        {
            "Step": "Derived percentages",
            "Random": False,
            "Seed": SEED,
            "Control": "fixed arithmetic and six-decimal rounding",
        },
        {
            "Step": "Figure layout",
            "Random": False,
            "Seed": SEED,
            "Control": "fixed factor order, palette, dimensions, and export DPI",
        },
    ]
    mirror_tsv(work, "determinism-audit.tsv", determinism)

    run_contract = {
        "article": 56,
        "seed": SEED,
        "primary_minimum_iphop_score": 90,
        "sensitivity_minimum_scores": [75, 95],
        "primary_claim_rank": "genus only when represented; family under novel-host audit",
        "host_domain_prescreen_required": True,
        "same_mag_bin_is_not_prophage_proof": True,
        "supplementary_spreadsheets_used": False,
        "random_output_requested": False,
    }
    dump_json(work / "run-contract.json", run_contract)
    dump_json(summary_dir / "run-contract.json", run_contract)

    summary = {
        "article": 56,
        "seed": SEED,
        "test_viral_genomes": 1870,
        "test_host_genera": 170,
        "temperate_stratum": 949,
        "virulent_stratum": 663,
        "neither_displayed_lifestyle_stratum": 258,
        "imgvr_hq_genomes": 216015,
        "eukaryotic_negative_controls": 8128,
        "eukaryotic_false_host_calls": 1018,
        "eukaryotic_false_host_call_pct": round(false_rate, 6),
        "false_calls_from_kmer_pct": 85.0,
        "false_calls_below_score90_pct": 90.0,
        "riboviria_false_calls": 640,
        "monodnaviria_false_calls": 155,
        "default_minimum_score": 90,
        "default_estimated_max_fdr_pct": 10,
        "published_runtime_minutes": 12,
        "published_runtime_genomes": 5,
        "published_runtime_cpus": 6,
        "published_iphop_version": "1.0",
        "published_database": "iPHoP_db_Sept21",
        "current_locked_iphop_version": "1.4.2",
        "evidence_tiers": 6,
        "source_assertions": len(source_assertions),
        "supplementary_spreadsheets_used": False,
        "random_output_requested": False,
    }
    dump_json(work / "summary.json", summary)
    dump_json(summary_dir / "summary.json", summary)

    command_log = [
        {
            "Step": "prepare",
            "Command": "python3 scripts/prepare_article56_host_evidence.py --project-root ${PROJECT_ROOT} --work-dir ${ARTICLE56_WORK_DIR}",
            "ExitStatus": 0,
            "Output": "checksum-gated source assertions and deterministic evidence tables",
        }
    ]
    write_tsv(work / "command-log.tsv", command_log)
    (logs / "prepare.log").write_text(
        "Article 56 source extraction completed.\n"
        f"iPHoP paragraphs (main/all): {len(iphop_main)}/{len(iphop_all)}\n"
        f"MIUViG main paragraphs: {len(miuvig_main)}\n"
        f"Assertions validated: {len(source_assertions)}\n"
        "Supplementary spreadsheets used: false\n",
        encoding="utf-8",
    )
    (work / ".article56-summary-complete").write_text(
        "Article 56 summary complete\n", encoding="utf-8"
    )
    print(f"Prepared Article 56 evidence: {work}")


if __name__ == "__main__":
    main()
