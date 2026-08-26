#!/usr/bin/env python3
"""Build the claim-specific causal-evidence ledger for Article 72."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import platform
import re
import shutil
from pathlib import Path

import pandas as pd


ARTICLE = 72
ANALYSIS_SEED = 72_001
PLOT_SEED = 20_260_772


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_text(value: object) -> str:
    text = html.unescape(str(value))
    text = re.sub(r"<[^>]+>", "", text)
    return re.sub(r"\s+", " ", text).strip()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
        float_format="%.15g",
    )


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((cache / "download-manifest.json").read_text())
    if manifest["article"] != ARTICLE or manifest["publication_count"] != 11:
        raise ValueError("Unexpected Article 72 download manifest")
    metadata_path = cache / "publication-metadata.json"
    anchor_path = cache / "buffie-figure4-original.jpg"
    metadata_record = manifest["resources"]["publication-metadata.json"]
    anchor_record = manifest["anchor"]
    for path, record in ((metadata_path, metadata_record), (anchor_path, anchor_record)):
        if path.stat().st_size != record["bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"Checksum mismatch for {path.name}")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata_index = {record["citation_key"]: record for record in metadata}
    if len(metadata_index) != 11:
        raise ValueError("Expected 11 unique publication records")
    publication_rows = []
    for record in metadata:
        authors = record.get("authors", [])
        publication_rows.append(
            {
                "CitationKey": record["citation_key"],
                "DOI": record["doi"],
                "Year": record["year"],
                "Title": clean_text(record["title"]),
                "Journal": clean_text(record["journal"]),
                "Volume": record["volume"],
                "Issue": record["issue"],
                "Pages": record["pages"],
                "FirstAuthor": clean_text(authors[0]["family"]) if authors else "",
                "DOIURL": record["source_url"],
            }
        )
    publications = pd.DataFrame(publication_rows).sort_values(
        ["Year", "CitationKey"], kind="stable"
    )
    write_tsv(publications, output / "publication-metadata.tsv")

    rungs = pd.DataFrame(
        [
            (1, "Replicated association", "Does the feature travel across cohorts?", "Independent human cohorts with harmonized measurement", "Population relevance and reproducibility", "Direction, intervention effect, or mechanism"),
            (2, "Temporality and reversibility", "Does the change precede the outcome and reverse as predicted?", "Prospective sampling, time-varying exposure, or controlled perturbation", "Ordering and within-unit change", "Exchangeability when exposure is not randomized"),
            (3, "Human causal bridge", "Do mediation, MR, negative controls, or natural experiments survive their assumptions?", "Design-specific causal estimand with sensitivity analyses", "Triangulation under a different bias structure", "A universal causal verdict or a treatment effect"),
            (4, "Randomized human intervention", "Does changing the intervention change a patient outcome?", "Random allocation, prespecified endpoint, intention-to-treat", "Effect of the tested intervention bundle", "Which strain, gene, or metabolite produced the effect"),
            (5, "Isolation or defined product", "Can the candidate strain/function be physically specified?", "Pure culture, defined consortium, or composition-controlled product", "Entity and dose specificity", "Efficacy in the target population"),
            (6, "Host transfer or perturbation", "Can the phenotype be transferred or modified in a susceptible host?", "FMT, gnotobiotic host, colonization, depletion, or rescue", "Experimental sufficiency in the model", "Human transportability or natural exposure equivalence"),
            (7, "Molecular perturbation and rescue", "Does a gene/metabolite perturbation alter the phenotype and can it be rescued?", "Knockout, inhibitor, complementation, metabolite add-back, or receptor block", "Mechanistic specificity within the model", "Population effect size or clinical benefit"),
        ],
        columns=[
            "RungOrder",
            "Rung",
            "Question",
            "MinimumDesign",
            "WhatItAdds",
            "DoesNotProve",
        ],
    )
    write_tsv(rungs, output / "rung-definitions.tsv")

    claims = pd.DataFrame(
        [
            {
                "ClaimID": "C1",
                "Case": "Recurrent CDI: clinical restoration",
                "ExactClaim": "After symptom resolution with standard-of-care antibiotics, a tested donor-derived microbiota intervention can reduce short-term recurrent CDI in trial-eligible adults.",
                "TargetPopulation": "Adults with recurrent CDI after response to standard-of-care antibiotics",
                "InterventionOrExposure": "Protocol-specific donor FMT or purified donor-derived Firmicutes spores",
                "Comparator": "Protocol-specific antibiotic or placebo control",
                "Outcome": "Sustained clinical response / recurrent CDI",
                "TimeHorizon": "8 to 10 weeks",
                "AllowedConclusion": "The tested intervention bundle changes recurrence risk in its trial population.",
                "ForbiddenLeap": "FMT efficacy proves that one named strain or bile acid is the active human mechanism.",
            },
            {
                "ClaimID": "C2",
                "Case": "C. scindens: bile-acid resistance",
                "ExactClaim": "Cultured bile-acid 7alpha-dehydroxylating C. scindens can increase colonization resistance to C. difficile in specified antibiotic-perturbed mouse and ex-vivo systems.",
                "TargetPopulation": "Specified antibiotic-perturbed experimental systems",
                "InterventionOrExposure": "C. scindens administration / bile-acid-dependent activity",
                "Comparator": "Vehicle, susceptible community, or bile-acid sequestration control",
                "Outcome": "C. difficile burden, mortality, growth inhibition, and bile-acid profile",
                "TimeHorizon": "Acute experimental challenge",
                "AllowedConclusion": "C. scindens is sufficient to alter resistance in those model systems through bile-acid-dependent activity.",
                "ForbiddenLeap": "C. scindens alone is the clinically proven active ingredient of FMT in recurrent CDI.",
            },
            {
                "ClaimID": "C3",
                "Case": "F. nucleatum: colorectal cancer phenotypes",
                "ExactClaim": "Tumor-associated F. nucleatum can promote colorectal-cancer phenotypes through strain- and context-dependent mechanisms in specified cell, mouse, and xenograft systems.",
                "TargetPopulation": "Human CRC cohorts plus specified experimental models",
                "InterventionOrExposure": "F. nucleatum, FadA signaling, or antibiotic depletion in model systems",
                "Comparator": "Non-tumor tissue, vehicle/control bacteria, peptide block, or antibiotic control",
                "Outcome": "Tumor burden, immune signaling, growth, persistence, or xenograft growth",
                "TimeHorizon": "Study-specific",
                "AllowedConclusion": "Human association and multiple experimental systems are mechanistically compatible with a tumor-promoting role.",
                "ForbiddenLeap": "F. nucleatum has been shown by a randomized human intervention to cause or treat CRC.",
            },
        ]
    )
    write_tsv(claims, output / "claim-contracts.tsv")

    evidence_rows = [
        {
            "EvidenceID": "E01", "ClaimID": "C1", "CitationKey": "vannood2013fmt", "StudyDesign": "Open-label randomized controlled trial", "PopulationModel": "Adults with recurrent CDI", "NLabel": "42 evaluable: 16 FMT, 13 vancomycin, 13 vancomycin plus lavage", "ShotgunStrainResolved": 0, "IndependentReplication": 0, "EvidenceRole": "Direct", "MainObservation": "Overall cure without relapse was 15/16 after donor FMT versus 4/13 and 3/13 in the two control groups.", "Boundary": "Tests a protocol-level donor-feces intervention; it does not isolate the active microbial component.",
            "Replicated association": "Supporting", "Temporality and reversibility": "Direct", "Human causal bridge": "Not addressed", "Randomized human intervention": "Direct", "Isolation or defined product": "Not addressed", "Host transfer or perturbation": "Direct", "Molecular perturbation and rescue": "Not addressed",
        },
        {
            "EvidenceID": "E02", "ClaimID": "C1", "CitationKey": "feuerstadt2022ser109", "StudyDesign": "Phase 3 double-blind randomized placebo-controlled trial", "PopulationModel": "Adults with recurrent CDI after standard-of-care antibiotics", "NLabel": "182: 89 SER-109, 93 placebo", "ShotgunStrainResolved": 0, "IndependentReplication": 1, "EvidenceRole": "Direct", "MainObservation": "Week-8 recurrence was 11/89 with SER-109 and 37/93 with placebo (RR 0.32, 95% CI 0.18-0.58).", "Boundary": "Tests a purified donor-derived Firmicutes spore product, not C. scindens as a single active ingredient.",
            "Replicated association": "Not addressed", "Temporality and reversibility": "Direct", "Human causal bridge": "Not addressed", "Randomized human intervention": "Direct", "Isolation or defined product": "Direct", "Host transfer or perturbation": "Direct", "Molecular perturbation and rescue": "Supporting",
        },
        {
            "EvidenceID": "E03", "ClaimID": "C1", "CitationKey": "smillie2018engraftment", "StudyDesign": "Longitudinal shotgun strain-tracking clinical experiment", "PopulationModel": "19 recurrent-CDI FMT recipients and four donors", "NLabel": "19 recipients; 4 donors", "ShotgunStrainResolved": 1, "IndependentReplication": 1, "EvidenceRole": "Supporting", "MainObservation": "Donor and recipient strain signatures were tracked over time; engraftment depended on abundance, phylogeny, and recipient context.", "Boundary": "Engraftment is a post-treatment event and is not, by itself, a randomized mediator of cure.",
            "Replicated association": "Supporting", "Temporality and reversibility": "Direct", "Human causal bridge": "Not addressed", "Randomized human intervention": "Not addressed", "Isolation or defined product": "Not addressed", "Host transfer or perturbation": "Direct", "Molecular perturbation and rescue": "Not addressed",
        },
        {
            "EvidenceID": "E04", "ClaimID": "C1", "CitationKey": "li2016fmtstrains", "StudyDesign": "Longitudinal shotgun SNV tracking after FMT", "PopulationModel": "Metabolic-syndrome FMT recipients", "NLabel": "Study-specific longitudinal recipients", "ShotgunStrainResolved": 1, "IndependentReplication": 1, "EvidenceRole": "Supporting", "MainObservation": "Donor and recipient strains coexisted for at least three months with recipient-specific transfer patterns.", "Boundary": "Shows strain persistence in another indication; it does not estimate recurrent-CDI treatment efficacy.",
            "Replicated association": "Supporting", "Temporality and reversibility": "Direct", "Human causal bridge": "Not addressed", "Randomized human intervention": "Not addressed", "Isolation or defined product": "Not addressed", "Host transfer or perturbation": "Direct", "Molecular perturbation and rescue": "Not addressed",
        },
        {
            "EvidenceID": "E05", "ClaimID": "C2", "CitationKey": "buffie2015cscindens", "StudyDesign": "Cross-species association, isolate transfer, metabolomics, and ex-vivo perturbation", "PopulationModel": "Antibiotic-perturbed mice, hospitalized patients, and ex-vivo intestinal content", "NLabel": "Multiple linked mouse, human, and ex-vivo experiments", "ShotgunStrainResolved": 0, "IndependentReplication": 0, "EvidenceRole": "Direct", "MainObservation": "C. scindens administration increased resistance, tracked with secondary bile acids, and bile sequestration altered ex-vivo inhibition.", "Boundary": "Demonstrates model-system sufficiency; it does not test C. scindens monotherapy in a randomized human trial.",
            "Replicated association": "Direct", "Temporality and reversibility": "Direct", "Human causal bridge": "Not addressed", "Randomized human intervention": "Not addressed", "Isolation or defined product": "Direct", "Host transfer or perturbation": "Direct", "Molecular perturbation and rescue": "Direct",
        },
        {
            "EvidenceID": "E06", "ClaimID": "C2", "CitationKey": "theriot2014metabolome", "StudyDesign": "Controlled antibiotic perturbation with microbiome-metabolome profiling and pathogen challenge", "PopulationModel": "Mouse C. difficile susceptibility model", "NLabel": "Study-specific mouse experiments", "ShotgunStrainResolved": 0, "IndependentReplication": 1, "EvidenceRole": "Supporting", "MainObservation": "Antibiotic-specific microbiome and bile-acid shifts preceded altered susceptibility to C. difficile challenge.", "Boundary": "Supports a bile-acid ecological mechanism but does not isolate C. scindens as the sole causal component.",
            "Replicated association": "Supporting", "Temporality and reversibility": "Direct", "Human causal bridge": "Not addressed", "Randomized human intervention": "Not addressed", "Isolation or defined product": "Not addressed", "Host transfer or perturbation": "Supporting", "Molecular perturbation and rescue": "Supporting",
        },
        {
            "EvidenceID": "E07", "ClaimID": "C3", "CitationKey": "wirbel2019crc", "StudyDesign": "Cross-cohort fecal-shotgun meta-analysis and external prediction", "PopulationModel": "768 fecal metagenomes from eight CRC studies", "NLabel": "768 profiles; 8 cohorts", "ShotgunStrainResolved": 0, "IndependentReplication": 1, "EvidenceRole": "Supporting", "MainObservation": "CRC-associated microbial signatures, including Fusobacterium-related signals, transferred across geographically distinct cohorts.", "Boundary": "Replicated case-control association remains vulnerable to disease effects, screening, treatment, diet, and stool-versus-tissue differences.",
            "Replicated association": "Direct", "Temporality and reversibility": "Not addressed", "Human causal bridge": "Not addressed", "Randomized human intervention": "Not addressed", "Isolation or defined product": "Not addressed", "Host transfer or perturbation": "Not addressed", "Molecular perturbation and rescue": "Not addressed",
        },
        {
            "EvidenceID": "E08", "ClaimID": "C3", "CitationKey": "kostic2013fusobacterium", "StudyDesign": "Human association plus cultured-bacterium administration in a tumor-prone mouse model", "PopulationModel": "Human adenoma/CRC samples and ApcMin/+ mice", "NLabel": "Human cohorts plus controlled mouse experiments", "ShotgunStrainResolved": 0, "IndependentReplication": 1, "EvidenceRole": "Direct", "MainObservation": "F. nucleatum was enriched in human lesions and accelerated tumorigenesis with a pro-inflammatory myeloid signature in ApcMin/+ mice.", "Boundary": "Mouse sufficiency does not estimate the effect of removing F. nucleatum from people at risk of CRC.",
            "Replicated association": "Supporting", "Temporality and reversibility": "Supporting", "Human causal bridge": "Not addressed", "Randomized human intervention": "Not addressed", "Isolation or defined product": "Direct", "Host transfer or perturbation": "Direct", "Molecular perturbation and rescue": "Supporting",
        },
        {
            "EvidenceID": "E09", "ClaimID": "C3", "CitationKey": "rubinstein2013fada", "StudyDesign": "Cultured-bacterium, adhesin, receptor-binding, and peptide-block experiments", "PopulationModel": "CRC cells and human colon tissues", "NLabel": "Study-specific cell and tissue experiments", "ShotgunStrainResolved": 0, "IndependentReplication": 1, "EvidenceRole": "Direct", "MainObservation": "FadA bound E-cadherin, activated beta-catenin signaling, and an E-cadherin-derived peptide blocked induced growth and signaling.", "Boundary": "Mechanistic specificity in cells/tissues does not establish population attributable risk or clinical benefit from targeting FadA.",
            "Replicated association": "Supporting", "Temporality and reversibility": "Not addressed", "Human causal bridge": "Not addressed", "Randomized human intervention": "Not addressed", "Isolation or defined product": "Direct", "Host transfer or perturbation": "Not addressed", "Molecular perturbation and rescue": "Direct",
        },
        {
            "EvidenceID": "E10", "ClaimID": "C3", "CitationKey": "bullman2017fusobacterium", "StudyDesign": "Paired primary-metastasis analysis, culture, xenograft transfer, and antibiotic perturbation", "PopulationModel": "Human CRC tissues and patient-derived xenografts", "NLabel": "Paired human tumors plus xenograft experiments", "ShotgunStrainResolved": 1, "IndependentReplication": 1, "EvidenceRole": "Direct", "MainObservation": "Concordant Fusobacterium strains persisted in primary and metastatic tumors; metronidazole reduced Fusobacterium load and xenograft growth.", "Boundary": "A broad anaerobe-active antibiotic in xenografts is not a targeted randomized treatment in patients with CRC.",
            "Replicated association": "Direct", "Temporality and reversibility": "Supporting", "Human causal bridge": "Not addressed", "Randomized human intervention": "Not addressed", "Isolation or defined product": "Direct", "Host transfer or perturbation": "Direct", "Molecular perturbation and rescue": "Supporting",
        },
    ]
    evidence = pd.DataFrame(evidence_rows)
    evidence = evidence.merge(
        publications[["CitationKey", "DOI", "Year", "Title", "Journal"]],
        on="CitationKey",
        how="left",
        validate="many_to_one",
    )
    if evidence[["DOI", "Year", "Title", "Journal"]].isna().any().any():
        raise ValueError("Evidence ledger contains an unknown citation key")
    ordered_columns = [
        "EvidenceID", "ClaimID", "CitationKey", "DOI", "Year", "Title", "Journal",
        "StudyDesign", "PopulationModel", "NLabel", "ShotgunStrainResolved",
        "IndependentReplication", "EvidenceRole", "MainObservation", "Boundary",
    ] + rungs["Rung"].tolist()
    evidence = evidence[ordered_columns]
    write_tsv(evidence, output / "evidence-ledger.tsv")

    coverage = evidence.melt(
        id_vars=["EvidenceID", "ClaimID", "CitationKey", "Year", "EvidenceRole"],
        value_vars=rungs["Rung"].tolist(),
        var_name="Rung",
        value_name="Coverage",
    ).merge(rungs[["RungOrder", "Rung"]], on="Rung", how="left", validate="many_to_one")
    coverage["CoverageScore"] = coverage["Coverage"].map(
        {"Not addressed": 0, "Supporting": 1, "Direct": 2}
    )
    if coverage["CoverageScore"].isna().any():
        raise ValueError("Unknown evidence-coverage code")
    coverage = coverage.sort_values(["ClaimID", "EvidenceID", "RungOrder"], kind="stable")
    write_tsv(coverage, output / "study-domain-coverage.tsv")

    outcomes = pd.DataFrame(
        [
            ("van Nood 2013", "Donor FMT", "Cure without relapse", 15, 16, 10, "Open-label randomized trial", "Overall response after up to two infusions"),
            ("van Nood 2013", "Vancomycin", "Cure without relapse", 4, 13, 10, "Open-label randomized trial", "Standard vancomycin arm"),
            ("van Nood 2013", "Vancomycin + lavage", "Cure without relapse", 3, 13, 10, "Open-label randomized trial", "Standard vancomycin plus bowel-lavage arm"),
            ("Feuerstadt 2022", "SER-109", "No recurrent CDI", 78, 89, 8, "Double-blind randomized trial", "Complement of 11/89 recurrent CDI"),
            ("Feuerstadt 2022", "Placebo", "No recurrent CDI", 56, 93, 8, "Double-blind randomized trial", "Complement of 37/93 recurrent CDI"),
        ],
        columns=[
            "Study", "Arm", "FavorableEndpoint", "FavorableEvents", "Total",
            "TimeWeeks", "Design", "EndpointNote",
        ],
    )
    outcomes["FavorableRate"] = outcomes["FavorableEvents"] / outcomes["Total"]
    write_tsv(outcomes, output / "human-intervention-outcomes.tsv")

    downgrades = pd.DataFrame(
        [
            ("FMT works because C. scindens restores secondary bile acids in patients.", "Donor FMT and a purified spore product reduce recurrent-CDI risk in their trial populations; C. scindens has bile-acid-dependent effects in specified experimental systems.", "Intervention bundle and single-strain mechanism come from different estimands."),
            ("Donor strains engrafted, therefore those strains caused the cure.", "Shotgun strain tracking establishes donor-origin persistence after FMT; a mediator effect on cure requires an additional identified design.", "Post-treatment association can reflect treatment delivery, host selection, or outcome-related survival."),
            ("F. nucleatum causes human colorectal cancer.", "F. nucleatum is reproducibly associated with human CRC and can promote CRC phenotypes in specified mouse, cell, and xenograft systems.", "No randomized targeted human intervention establishes prevention or treatment benefit."),
            ("A significant mediation or MR result proves a microbial treatment will work.", "Mediation and MR provide design-specific causal bridges under explicit assumptions; their estimands need not equal a short-term microbial intervention effect.", "Unmeasured mediator-outcome confounding and horizontal pleiotropy remain design-specific threats."),
            ("Mouse FMT reproduced the phenotype, so the same mechanism operates in people.", "Host-transfer experiments establish sufficiency in that model and intervention context; human transport requires matched exposure, outcome, ecology, and dose.", "Species, diet, housing, immune state, and donor-recipient compatibility can modify the effect."),
        ],
        columns=["Overclaim", "SupportedClaim", "WhyDowngraded"],
    )
    write_tsv(downgrades, output / "claim-downgrade-examples.tsv")

    packet = pd.DataFrame(
        [
            (1, "Claim contract", "Population, intervention/exposure, comparator, outcome, horizon, microbial resolution", "Before literature synthesis"),
            (2, "Human lineage", "Cohorts, eligibility, sampling time, profiling/database release, attrition", "Association and longitudinal evidence"),
            (3, "Bias ledger", "Confounding, reverse causation, selection, measurement error, multiple testing", "Every observational design"),
            (4, "Perturbation ledger", "Randomization, allocation, dose, adherence, contamination, intention-to-treat", "Human and model interventions"),
            (5, "Entity ledger", "Strain/isolate/consortium, gene, metabolite, receptor, viability and dose", "Culture and mechanism"),
            (6, "Transport ledger", "Host species, background community, diet, antibiotic state, tissue, outcome scale", "FMT/gnotobiotic/model-to-human"),
            (7, "Contradiction ledger", "Null results, sign changes, failed replication, boundary conditions", "Before wording the conclusion"),
            (8, "Reproducibility bundle", "DOIs, accessions, tables, code, versions, checksums, seeds", "Submission and review"),
        ],
        columns=["Order", "Packet", "RequiredFields", "UsedFor"],
    )
    write_tsv(packet, output / "evidence-packet.tsv")

    shutil.copy2(anchor_path, output / "buffie-figure4-original.jpg")
    source_manifest = dict(manifest)
    source_manifest.update(
        {
            "article": ARTICLE,
            "publication_rows": len(publications),
            "evidence_rows": len(evidence),
            "claim_count": len(claims),
            "rung_count": len(rungs),
            "intervention_arm_rows": len(outcomes),
            "framework_citation": "10.1016/j.mib.2017.10.001",
        }
    )
    (output / "source-manifest.json").write_text(
        json.dumps(source_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    methods_contract = {
        "article": ARTICLE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "unit": "one primary publication by one prespecified claim",
        "synthesis": "claim-specific evidence map without an additive causal score",
        "coverage_codes": ["Not addressed", "Supporting", "Direct"],
        "claims": claims["ClaimID"].tolist(),
        "rungs": rungs["Rung"].tolist(),
        "interpretation_limit": (
            "Coverage across evidence domains is not a quality-weighted meta-analysis; "
            "each study supports only the claim, intervention, host, and outcome it tested."
        ),
    }
    (output / "methods-contract.json").write_text(
        json.dumps(methods_contract, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    metrics = {
        "article": ARTICLE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "publication_count": len(publications),
        "evidence_rows": len(evidence),
        "claim_count": len(claims),
        "rung_count": len(rungs),
        "direct_cells": int((coverage["Coverage"] == "Direct").sum()),
        "supporting_cells": int((coverage["Coverage"] == "Supporting").sum()),
        "not_addressed_cells": int((coverage["Coverage"] == "Not addressed").sum()),
        "shotgun_strain_resolved_studies": int(evidence["ShotgunStrainResolved"].sum()),
        "human_randomized_studies": 2,
        "crc_human_randomized_targeted_interventions": 0,
        "python": platform.python_version(),
        "pandas": pd.__version__,
    }
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
