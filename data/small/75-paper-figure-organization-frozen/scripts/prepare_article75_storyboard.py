#!/usr/bin/env python3
"""Build Article 75's paper-figure storyboard from public, checksum-locked evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from PIL import Image


ARTICLE = 75
ANALYSIS_SEED = 75_001
PLOT_SEED = 20_260_775
SNAPSHOT_DATE = "2026-08-23"
PMCID = "PMC7984229"
DOI = "10.1038/s41591-019-0406-6"
ANCHOR_MEMBER = "emss-81948-f001.jpg"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n", float_format="%.10g")


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def normalized_text(node: ET.Element) -> str:
    return " ".join("".join(node.itertext()).split())


def paper_sources(evidence: Path, output: Path) -> pd.DataFrame:
    manifest = read_json(evidence / "download-manifest.json")
    if manifest.get("article") != ARTICLE or manifest["paper"]["doi"] != DOI:
        raise ValueError("Unexpected Article 75 download manifest")
    xml_path = evidence / "PMC7984229-fulltext.xml"
    zip_path = evidence / "PMC7984229-supplementary.zip"
    xml_record = manifest["full_text_xml"]
    if xml_path.stat().st_size != xml_record["bytes"] or sha256(xml_path) != xml_record["sha256"]:
        raise ValueError("Full-text XML failed checksum verification")

    root = ET.parse(xml_path).getroot()
    if root.findtext('.//article-id[@pub-id-type="doi"]') != DOI:
        raise ValueError("Full-text XML DOI mismatch")
    xlink = "{http://www.w3.org/1999/xlink}href"
    captions: dict[int, tuple[str, str]] = {}
    with zipfile.ZipFile(zip_path) as archive:
        for figure in root.findall(".//fig"):
            label = figure.find("label")
            caption = figure.find("caption")
            graphic = figure.find("graphic")
            if label is None or caption is None or graphic is None:
                continue
            label_text = normalized_text(label)
            if label_text not in {f"Figure {number}" for number in range(1, 6)}:
                continue
            number = int(label_text.split()[1])
            member = graphic.attrib.get(xlink, "")
            expected_member = f"emss-81948-f{number:03d}.jpg"
            if member != expected_member:
                raise ValueError(f"Unexpected Figure {number} member: {member}")
            payload = archive.read(member)
            record = manifest["selected_members"][member]
            if len(payload) != record["bytes"] or sha256_bytes(payload) != record["sha256"]:
                raise ValueError(f"Figure member checksum mismatch: {member}")
            captions[number] = (normalized_text(caption), member)
        anchor = archive.read(ANCHOR_MEMBER)

    if sorted(captions) != [1, 2, 3, 4, 5]:
        raise ValueError(f"Expected five main figures, observed {sorted(captions)}")
    anchor_path = output / "wirbel-figure1-original.jpg"
    anchor_path.write_bytes(anchor)
    with Image.open(anchor_path) as image:
        width, height = image.size
    if (width, height) != (800, 941):
        raise ValueError(f"Unexpected anchor dimensions: {(width, height)}")

    roles = {
        1: (
            4,
            "Cross-study discovery",
            "574 independent observations",
            "Blocked meta-analysis plus explicit study-effect audit",
            "A core set of species is reproducibly associated with CRC despite study heterogeneity.",
        ),
        2: (
            4,
            "Phenotypic structure",
            "285 independent CRC samples",
            "Study-blocked subgroup tests",
            "Co-occurring marker clusters expose heterogeneity among CRC cases.",
        ),
        3: (
            5,
            "Transportability",
            "Complete studies as train/test domains",
            "Study-to-study transfer, LOSO, and other-disease specificity",
            "Taxonomic and functional models transfer better when trained across studies.",
        ),
        4: (
            6,
            "Function and orthogonal assay",
            "574 profiles; 47-sample qPCR subset",
            "Functional meta-analysis plus DNA/RNA qPCR",
            "Functional associations are narrowed to a bai-operon signal with an orthogonal assay.",
        ),
        5: (
            3,
            "Independent validation",
            "193 independent observations in three populations",
            "Untouched external cohorts",
            "Taxonomic, functional, and selected virulence signals are tested outside discovery.",
        ),
    }
    rows = []
    for number in range(1, 6):
        caption, member = captions[number]
        panels, role, unit, validation, takeaway = roles[number]
        rows.append(
            {
                "Figure": number,
                "SourceLabel": f"Figure {number}",
                "SourceImage": member,
                "Panels": panels,
                "NarrativeRole": role,
                "StatisticalUnit": unit,
                "ValidationRole": validation,
                "TakeHomeParaphrase": takeaway,
                "CaptionSHA256": hashlib.sha256(caption.encode("utf-8")).hexdigest(),
                "ImageBytes": manifest["selected_members"][member]["bytes"],
                "ImageSHA256": manifest["selected_members"][member]["sha256"],
            }
        )
    ledger = pd.DataFrame(rows)
    write_tsv(ledger, output / "wirbel-main-figure-ledger.tsv")

    source_manifest = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "paper": manifest["paper"],
        "pmcid": PMCID,
        "xml": {
            "url": xml_record["url"],
            "bytes": xml_record["bytes"],
            "sha256": xml_record["sha256"],
        },
        "figure_endpoint": manifest["supplementary_endpoint"],
        "selected_figure_members": manifest["selected_members"],
        "anchor": {
            "member": ANCHOR_MEMBER,
            "bytes": len(anchor),
            "sha256": sha256_bytes(anchor),
            "width": width,
            "height": height,
        },
        "rights_boundary": (
            "The original anchor remains under the publisher/rightsholder terms and "
            "is reproduced only for attributed scholarly commentary; it is excluded "
            "from the repository CC BY/MIT grant."
        ),
    }
    (output / "paper-source-manifest.json").write_text(
        json.dumps(source_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return ledger


def source_artifacts(root: Path, output: Path) -> dict[str, Path]:
    selected = {
        "A28_cohorts": root / "data/small/28-cross-cohort/cohort-summary.tsv",
        "A28_lodo": root / "results/28-cross-cohort/lodo-performance.tsv",
        "A28_meta": root / "results/28-cross-cohort/performance-meta-analysis.tsv",
        "A28_bootstrap": root / "results/28-cross-cohort/hierarchical-bootstrap-summary.tsv",
        "A44_mimag": root / "data/small/44-mag-qc-mimag-graph-frozen/mimag-tier-counts.tsv",
        "A44_summary": root / "data/small/44-mag-qc-mimag-graph-frozen/run-summary.json",
        "A50_summary": root / "data/small/50-instrain-microdiversity-frozen/run-summary.json",
        "A51_summary": root / "data/small/51-strainphlan-frozen/run-summary.json",
        "A54_summary": root / "data/small/54-virus-discovery-quality-frozen/summary.json",
        "A55_summary": root / "data/small/55-virus-taxonomy-abundance-frozen/summary.json",
        "A56_summary": root / "data/small/56-virus-host-evidence-frozen/summary.json",
        "A60_summary": root / "data/small/60-genome-scale-models-frozen/summary.json",
        "A61_metrics": root / "data/small/61-community-metabolism-frozen/analysis-metrics.json",
        "A72_metrics": root / "data/small/72-causal-evidence-frozen/analysis-metrics.json",
        "A74_databases": root / "data/small/74-nfcore-mag-workflows-frozen/database-lock.tsv",
    }
    missing = [str(path) for path in selected.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing source artifacts: {missing}")
    rows = []
    for key, path in selected.items():
        rows.append(
            {
                "ArtifactID": key,
                "SourceArticle": int(key[1:3]),
                "RelativePath": path.relative_to(root).as_posix(),
                "Bytes": path.stat().st_size,
                "SHA256": sha256(path),
                "Role": {
                    "A28_cohorts": "Cross-cohort design",
                    "A28_lodo": "External validation",
                    "A28_meta": "Performance heterogeneity",
                    "A28_bootstrap": "External uncertainty",
                    "A44_mimag": "MAG quality tiers",
                    "A44_summary": "MAG acceptance audit",
                    "A50_summary": "Within-species comparison",
                    "A51_summary": "Strain phylogeny",
                    "A54_summary": "Virus discovery and quality",
                    "A55_summary": "vOTU abundance and taxonomy",
                    "A56_summary": "Virus-host claim ceiling",
                    "A60_summary": "Genome-scale model QC",
                    "A61_metrics": "Community metabolism sensitivity",
                    "A72_metrics": "Causal evidence ceiling",
                    "A74_databases": "Database version ledger",
                }[key],
            }
        )
    write_tsv(pd.DataFrame(rows), output / "source-artifact-manifest.tsv")
    return selected


def evidence_metrics(paths: dict[str, Path], output: Path) -> pd.DataFrame:
    cohorts = pd.read_csv(paths["A28_cohorts"], sep="\t")
    lodo = pd.read_csv(paths["A28_lodo"], sep="\t")
    meta = pd.read_csv(paths["A28_meta"], sep="\t").iloc[0]
    boot = pd.read_csv(paths["A28_bootstrap"], sep="\t").iloc[0]
    mag = read_json(paths["A44_summary"])
    instrain = read_json(paths["A50_summary"])
    strainphlan = read_json(paths["A51_summary"])
    virus = read_json(paths["A54_summary"])
    votu = read_json(paths["A55_summary"])
    host = read_json(paths["A56_summary"])
    models = read_json(paths["A60_summary"])
    community = read_json(paths["A61_metrics"])
    causal = read_json(paths["A72_metrics"])
    rows = [
        ("E01", "Cross-cohort", "Cohorts", len(cohorts), "cohorts", "A28_cohorts", "Eight complete domains; not random folds."),
        ("E02", "Cross-cohort", "Independent subjects", int(cohorts["Samples"].sum()), "subjects", "A28_cohorts", "CRC/control profiles after cohort-specific exclusions."),
        ("E03", "External validation", "Macro AUROC", float(boot["Estimate"]), "AUROC", "A28_bootstrap", f"95% CI {boot['Lower95']:.3f}-{boot['Upper95']:.3f}; no clinical utility claim."),
        ("E04", "External validation", "Performance heterogeneity", float(meta["I2Percent"]), "I2 percent", "A28_meta", "LODO heterogeneity is a result, not noise to hide."),
        ("E05", "MAG", "Accepted MAGs", int(mag["selected_bins"]), "MAGs", "A44_summary", f"{mag['mimag_counts']['High quality']} HQ and {mag['mimag_counts']['Medium quality']} MQ; no circular genomes claimed."),
        ("E06", "Strain", "Compared genomes", int(instrain["genomes_compared"]), "genomes", "A50_summary", f"{instrain['same_strain_calls']} same-strain calls in a two-sample benchmark."),
        ("E07", "Strain", "Phylogeny alignment", int(strainphlan["alignment_sites"]), "sites", "A51_summary", f"{strainphlan['retained_samples']} samples and {strainphlan['selected_markers']} markers; geography withheld from inference."),
        ("E08", "Virus", "CheckV input", int(virus["input_sequences"]), "viral sequences", "A54_summary", f"Quality: {virus['checkv_quality_counts']['High-quality']} HQ, {virus['checkv_quality_counts']['Medium-quality']} MQ, {virus['checkv_quality_counts']['Low-quality']} LQ, {virus['checkv_quality_counts']['Not-determined']} ND."),
        ("E09", "Virus", "Mock-community vOTUs", int(votu["votu_clusters"]), "vOTUs", "A55_summary", f"Derived from {votu['mock_phages']} mock phages; not a clinical cohort."),
        ("E10", "Virus-host", "Negative-control genomes", int(host["eukaryotic_negative_controls"]), "genomes", "A56_summary", f"{host['eukaryotic_false_host_calls']} false calls demonstrate the claim ceiling."),
        ("E11", "Metabolic models", "Reconstructed models", int(models["models"]), "models", "A60_summary", f"{models['high_gapfill_burden_flags']} high gap-fill burden flags across {models['input_genomes']} input genomes."),
        ("E12", "Community metabolism", "Modeled subjects", int(community["independent_subjects"]), "subjects", "A61_metrics", f"{community['primary_positive_taxa']} of {community['primary_total_taxa']} modeled taxa had positive primary growth."),
        ("E13", "Causal evidence", "CRC randomized targeted interventions", int(causal["crc_human_randomized_targeted_interventions"]), "studies", "A72_metrics", "Zero: association and models do not become human causality."),
    ]
    frame = pd.DataFrame(
        rows,
        columns=["EvidenceID", "Domain", "Metric", "Value", "Unit", "ArtifactID", "EvidenceBoundary"],
    )
    if int(lodo["Samples"].sum()) != 771 or len(lodo) != 8:
        raise ValueError("Article 28 LODO evidence changed unexpectedly")
    write_tsv(frame, output / "series-evidence-metrics.tsv")
    return frame


def storyboard_tables(output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    storyboard = pd.DataFrame(
        [
            (1, "Design and global map", "Who and what was measured?", "subject / sample / cohort", "Flow, missingness, sequencing depth, global effect", 4, "Cohort and negative-control balance", "Stop if exclusions or batch structure are unresolved"),
            (2, "Reproducible taxa and functions", "Which associations repeat across domains?", "feature within independent subject", "Effect size + 95% CI + FDR + prevalence", 5, "Direction and magnitude in held-out cohorts", "Do not promote a single-cohort hit"),
            (3, "MAG, strain, and virus discovery", "What genome-resolved entities were recovered?", "MAG / strain comparison / vOTU", "Quality gate, breadth, ANI, tree, host-evidence tier", 6, "Threshold and database sensitivity", "Split the figure if one branch lacks validation"),
            (4, "Mechanism or metabolic model", "What process is compatible with the observations?", "genome model / metabolite / experimental unit", "Flux, gap-fill burden, medium sensitivity, orthogonal assay", 5, "Perturbation, rescue, or independent assay", "Use compatibility language without intervention evidence"),
            (5, "External or experimental validation", "Does the central claim survive new data?", "untouched cohort or experimental unit", "Effect, calibration, AUROC, assay endpoint", 5, "Locked pipeline and prespecified endpoint", "Keep this dataset out of all selection"),
        ],
        columns=["Figure", "ShortTitle", "Question", "PrimaryUnit", "PrimaryStatistic", "Panels", "ValidationGate", "StopRule"],
    )
    write_tsv(storyboard, output / "main-figure-storyboard.tsv")

    panel_rows = [
        (1, "A", "Context", "Enrollment and exclusions", "subject", "Count and reason", "Sample ledger frozen"),
        (1, "B", "Context", "Cohort and covariate map", "cohort", "Balance and missingness", "Metadata dictionary frozen"),
        (1, "C", "Diagnostic", "Read, host, and negative-control QC", "sample", "Depth and contamination", "QC threshold prespecified"),
        (1, "D", "Primary", "Global community structure", "subject", "Effect size with cohort adjustment", "No joint train-test correction"),
        (2, "A", "Primary", "Species effects", "species", "Effect + CI + FDR", "Cohort-aware model"),
        (2, "B", "Primary", "Functional effects", "pathway / gene family", "Effect + CI + FDR", "Independent denominator"),
        (2, "C", "Validation", "Cross-cohort concordance", "feature by cohort", "Direction and heterogeneity", "Complete cohort holdout"),
        (2, "D", "Secondary", "Taxon-function linkage", "feature pair", "Association with multiplicity control", "No causal arrow"),
        (2, "E", "Sensitivity", "Prevalence and transform sensitivity", "feature", "Stable sign and rank", "Prespecified grid"),
        (3, "A", "Diagnostic", "Assembly and MAG acceptance", "MAG", "Completeness, contamination, GUNC", "MIMAG table complete"),
        (3, "B", "Primary", "Genome phylogeny and novelty", "MAG / species", "Tree and ANI", "Reference release locked"),
        (3, "C", "Primary", "Within-species structure", "sample pair", "Breadth, popANI, p-distance", "Coverage gate passed"),
        (3, "D", "Primary", "Virus discovery and vOTUs", "viral contig / vOTU", "CheckV tier, ANI, alignment fraction", "MIUViG fields complete"),
        (3, "E", "Validation", "Virus-host linkage", "vOTU-host pair", "Evidence tier and score", "Negative control reported"),
        (3, "F", "Sensitivity", "Genome-resolved threshold audit", "entity", "Yield versus confidence", "Alternative thresholds shown"),
        (4, "A", "Context", "Mechanistic hypothesis", "claim", "Directed evidence chain", "No arrow without evidence"),
        (4, "B", "Primary", "Metabolic prediction", "model / metabolite", "Flux or growth", "Medium defined"),
        (4, "C", "Sensitivity", "Gap-fill and medium dependence", "model", "Range across settings", "No single-medium conclusion"),
        (4, "D", "Validation", "Orthogonal assay", "experimental unit", "Agreement and uncertainty", "Blinded or held-out where possible"),
        (4, "E", "Validation", "Perturbation or rescue", "experimental unit", "Endpoint effect", "Comparator and time fixed"),
        (5, "A", "Validation", "Untouched cohort effect", "cohort / subject", "Effect + CI", "No refitting"),
        (5, "B", "Validation", "External prediction", "cohort / subject", "AUROC, AUPRC, calibration", "Threshold locked"),
        (5, "C", "Sensitivity", "Domain-shift audit", "cohort", "Performance heterogeneity", "Failure cohort retained"),
        (5, "D", "Validation", "Independent experiment", "experimental unit", "Prespecified endpoint", "Biological replicates"),
        (5, "E", "Secondary", "Claim boundary", "claim", "Highest supported evidence rung", "Forbidden leap stated"),
    ]
    panels = pd.DataFrame(
        panel_rows,
        columns=["Figure", "Panel", "PanelRole", "Question", "AnalysisUnit", "Statistic", "Gate"],
    )
    panels["ResultID"] = [f"F{row.Figure}{row.Panel}" for row in panels.itertuples(index=False)]
    panels["ClaimID"] = [f"C{number}" for number in panels["Figure"]]
    write_tsv(panels, output / "panel-register.tsv")
    return storyboard, panels


def supplement_table(output: Path) -> pd.DataFrame:
    rows = [
        ("S01", "Study flow and exclusion reasons", 2, 2, "Main flow; full sample ledger and missingness in supplement"),
        ("S02", "Read QC, host removal, controls", 1, 2, "Main diagnostic summary; all per-sample metrics retained"),
        ("S03", "All taxa and functional effects", 1, 2, "Selected effects in main; complete effect table in supplement"),
        ("S04", "Compositional and prevalence sensitivity", 1, 2, "One stability panel in main; full grid in supplement"),
        ("S05", "Cross-cohort transfer and calibration", 2, 2, "External performance is a main result; all folds and curves are archived"),
        ("S06", "Assembly and binning diagnostics", 1, 2, "Yield summary in main; contig/bin diagnostics in supplement"),
        ("S07", "MAG acceptance table", 1, 2, "Representative genomes in main; every MAG and reason in supplement"),
        ("S08", "Strain breadth, ANI, and tree robustness", 1, 2, "Key lineage in main; thresholds and pairwise matrices in supplement"),
        ("S09", "Virus detection, CheckV, and vOTU audit", 1, 2, "Validated discoveries in main; complete MIUViG ledger in supplement"),
        ("S10", "Virus-host evidence and negatives", 1, 2, "High-confidence links in main; all evidence tiers and negatives retained"),
        ("S11", "Metabolic models and gap-fill burden", 1, 2, "Mechanistic result in main; every model and reaction audit in supplement"),
        ("S12", "Medium, abundance, and solver sensitivity", 1, 2, "One sensitivity panel in main; full parameter surface retained"),
        ("S13", "Software, database, container, checksum", 0, 2, "Version ledger belongs in Methods and supplement, not a decorative panel"),
        ("S14", "Negative, null, and failed analyses", 0, 2, "Failures stay visible in supplement and result ledger"),
    ]
    frame = pd.DataFrame(rows, columns=["SectionID", "EvidenceBlock", "MainSpace", "SupplementDetail", "PlacementRule"])
    write_tsv(frame, output / "main-supplement-map.tsv")
    return frame


def claim_tables(output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    layers = [
        "Study design and QC",
        "Adjusted association",
        "Independent replication",
        "Genome-resolved evidence",
        "Model or orthogonal assay",
        "Perturbation or intervention",
    ]
    claims = {
        "C1 Cohort-level pattern": [2, 2, 1, -1, -1, -1],
        "C2 Taxon/function biomarker": [2, 2, 2, 1, 1, 0],
        "C3 MAG/strain/virus discovery": [2, 1, 2, 2, 1, 0],
        "C4 Mechanistic compatibility": [2, 1, 2, 2, 2, 1],
        "C5 Clinical or causal effect": [2, 1, 2, 1, 1, 2],
    }
    labels = {-1: "Not applicable", 0: "Insufficient alone", 1: "Supportive", 2: "Required"}
    rows = []
    for claim, values in claims.items():
        for layer, value in zip(layers, values, strict=True):
            rows.append((claim.split()[0], claim.split(" ", 1)[1], layer, value, labels[value]))
    matrix = pd.DataFrame(rows, columns=["ClaimID", "Claim", "EvidenceLayer", "RequirementCode", "Requirement"])
    write_tsv(matrix, output / "claim-evidence-matrix.tsv")

    ladder = pd.DataFrame(
        [
            (1, "Description", "Observed in this dataset", "causes, protects, or treats"),
            (2, "Adjusted association", "Associated after prespecified adjustment", "independent biomarker or mechanism"),
            (3, "Independent replication", "Replicated across untouched domains", "general population utility"),
            (4, "Genome-resolved compatibility", "Compatible with a strain, host, or pathway hypothesis", "active interaction or flux in vivo"),
            (5, "Perturbation and rescue", "Contributes in the tested model system", "randomized human efficacy"),
            (6, "Randomized human intervention", "Changes the prespecified outcome in the trial population", "universal mechanism beyond the tested intervention"),
        ],
        columns=["Rung", "Evidence", "AllowedLanguage", "ForbiddenLeap"],
    )
    write_tsv(ladder, output / "evidence-language-ladder.tsv")
    return matrix, ladder


def sensitivity_table(output: Path) -> pd.DataFrame:
    axes = {
        "Cohort and batch specification": [2, 2, 1, 1, 2],
        "Prevalence and abundance filter": [1, 2, 1, 1, 1],
        "Compositional transform": [1, 2, 0, 1, 2],
        "Outer-domain split": [1, 2, 1, 1, 2],
        "MAG quality threshold": [0, 0, 2, 2, 1],
        "Strain breadth and ANI": [0, 0, 2, 1, 1],
        "Virus, vOTU, and host threshold": [0, 0, 2, 1, 1],
        "Medium, abundance, and gap filling": [0, 0, 1, 2, 1],
        "Database and software release": [1, 1, 2, 2, 2],
    }
    label = {0: "Not primary", 1: "Conditional", 2: "Required"}
    rows = []
    for axis, values in axes.items():
        for figure, value in enumerate(values, start=1):
            rows.append((axis, figure, value, label[value]))
    frame = pd.DataFrame(rows, columns=["SensitivityAxis", "Figure", "RequirementCode", "Requirement"])
    write_tsv(frame, output / "sensitivity-matrix.tsv")
    return frame


def audit_and_style_tables(panels: pd.DataFrame, output: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    attack = pd.DataFrame(
        [
            ("Outcome leakage", 5, 5, "Outer-domain split and frozen preprocessing", "Main Fig 5 + supplement"),
            ("Cohort confounding", 5, 5, "Cohort-adjusted effects and study-effect panel", "Main Figs 1-2"),
            ("Missing external validation", 5, 5, "Untouched cohort or experiment", "Main Fig 5"),
            ("Compositional artifact", 4, 4, "Alternative denominator and transform grid", "Main Fig 2 + supplement"),
            ("Incomplete MAG QC", 4, 5, "CheckM2, GUNC, MIMAG acceptance table", "Main Fig 3 + supplement"),
            ("Strain threshold drift", 3, 4, "Breadth, ANI, marker, and topology sensitivity", "Main Fig 3 + supplement"),
            ("Virus false positive", 4, 5, "Two callers, CheckV, vOTU audit", "Main Fig 3 + supplement"),
            ("Unsupported host link", 4, 5, "Evidence tiers and negative controls", "Main Fig 3 + supplement"),
            ("Gap-fill dependence", 4, 4, "Medium and gap-fill burden surface", "Main Fig 4 + supplement"),
            ("Version drift", 4, 4, "Database, software, container, and checksum ledger", "Methods + supplement"),
        ],
        columns=["Attack", "Likelihood", "Impact", "Defense", "Location"],
    )
    write_tsv(attack, output / "reviewer-attack-map.tsv")

    style = pd.DataFrame(
        [
            ("Control", "Color", "#0072B2", "circle", "solid", "Reference group only"),
            ("Case", "Color", "#D55E00", "circle", "solid", "Case group only"),
            ("Discovery cohort", "Fill", "#009E73", "circle", "solid", "Training/discovery role"),
            ("External cohort", "Outline", "#CC79A7", "diamond", "solid", "Untouched validation role"),
            ("Accepted entity", "Status", "#2A9D8F", "circle", "solid", "Passed prespecified gate"),
            ("Excluded entity", "Status", "#8D99AE", "cross", "solid", "Excluded with reason"),
            ("95% confidence interval", "Line", "#264653", "none", "solid", "Uncertainty only"),
            ("Prespecified threshold", "Line", "#E9C46A", "none", "dashed", "Decision threshold only"),
        ],
        columns=["Semantic", "Channel", "Hex", "Marker", "LineType", "DoNotReuseFor"],
    )
    write_tsv(style, output / "figure-style-contract.tsv")

    trace = panels[["ResultID", "ClaimID", "Figure", "Panel", "AnalysisUnit", "Statistic", "Gate"]].copy()
    trace["SourceArtifact"] = trace["Figure"].map({1: "A28_cohorts", 2: "A28_meta", 3: "A44/A50/A51/A54/A55/A56", 4: "A60/A61/A72", 5: "A28_lodo/A28_bootstrap"})
    trace["CodeTarget"] = trace["ResultID"].map(lambda value: f"analysis/{value.lower()}.R")
    trace["OutputTarget"] = trace["ResultID"].map(lambda value: f"results/{value.lower()}.tsv")
    trace["UnitRecorded"] = True
    trace["StatisticRecorded"] = True
    trace["ValidationRecorded"] = True
    trace["ChecksumRequired"] = True
    write_tsv(trace, output / "result-traceability-ledger.tsv")
    return attack, style, trace


def version_ledger(paths: dict[str, Path], output: Path) -> pd.DataFrame:
    databases = pd.read_csv(paths["A74_databases"], sep="\t")
    rows = [
        {
            "Layer": "Read-based profile",
            "Resource": "curatedMetagenomicData",
            "Release": "3.12.0; resource snapshot 2021-03-31",
            "Artifact": "MetaPhlAn 3 species profiles; CHOCOPhlAn 201901",
            "Checksum": "Per-cohort SHA256 in Article 28 resource manifest",
            "Use": "Cross-cohort tutorial example",
        }
    ]
    for row in databases.itertuples(index=False):
        rows.append(
            {
                "Layer": "Genome-resolved database",
                "Resource": row.Database,
                "Release": row.Release,
                "Artifact": row.Artifact,
                "Checksum": f"{row.ChecksumType}:{row.Checksum}",
                "Use": "Version-ledger example; not a shared biological cohort",
            }
        )
    frame = pd.DataFrame(rows)
    write_tsv(frame, output / "version-ledger-example.tsv")
    return frame


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    evidence = args.evidence_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    paper = paper_sources(evidence, output)
    artifacts = source_artifacts(root, output)
    metrics = evidence_metrics(artifacts, output)
    storyboard, panels = storyboard_tables(output)
    supplement = supplement_table(output)
    claims, ladder = claim_tables(output)
    sensitivity = sensitivity_table(output)
    attacks, style, trace = audit_and_style_tables(panels, output)
    versions = version_ledger(artifacts, output)

    notice = """Article 75 data notice

The anchor image is Wirbel et al. 2019, Nature Medicine, Figure 1, obtained from
Europe PMC for attributed scholarly commentary. Copyright and reuse rights remain
with the original rightsholder; this image is excluded from this repository's
CC BY 4.0 / MIT license grant. The full XML and transport ZIP are download-time
evidence and are not redistributed in the frozen bundle. Only hashes, paraphrased
figure roles, and the attributed anchor are retained.

The numerical examples from Articles 28, 44, 50, 51, 54, 55, 56, 60, 61, 72,
and 74 come from different public datasets and validation fixtures. They illustrate
what evidence belongs in a figure or supplement. They must not be combined into a
single biological claim or represented as one cohort.
"""
    (output / "data-NOTICE.txt").write_text(notice, encoding="utf-8")

    methods = {
        "article": ARTICLE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "primary_unit_rule": "Every panel declares its independent analysis unit.",
        "main_figure_rule": "One claim, one primary result, one validation gate per main figure.",
        "supplement_rule": "Complete QC, entity ledgers, versions, negatives, and sensitivity grids remain visible.",
        "validation_rule": "External cohorts or experiments remain untouched by selection and threshold tuning.",
        "language_rule": "Claim wording cannot exceed the highest supported evidence rung.",
        "source_boundary": "Series metrics are heterogeneous tutorial examples, not a combined biological study.",
    }
    (output / "methods-contract.json").write_text(
        json.dumps(methods, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    analysis = {
        "article": ARTICLE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "snapshot_date": SNAPSHOT_DATE,
        "anchor_paper_figures": len(paper),
        "anchor_paper_panels": int(paper["Panels"].sum()),
        "storyboard_figures": len(storyboard),
        "storyboard_panels": len(panels),
        "supplement_blocks": len(supplement),
        "claim_matrix_cells": len(claims),
        "evidence_rungs": len(ladder),
        "sensitivity_cells": len(sensitivity),
        "reviewer_attacks": len(attacks),
        "traceability_results": len(trace),
        "series_evidence_metrics": len(metrics),
        "version_records": len(versions),
        "python": platform.python_version(),
        "pandas": pd.__version__,
    }
    (output / "analysis-metrics.json").write_text(
        json.dumps(analysis, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"prepared\t{output}\t{len(panels)} panel contracts\t{len(metrics)} evidence metrics")


if __name__ == "__main__":
    main()
