#!/usr/bin/env python3
"""Build Article 77's local release packet and submission-readiness audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from PIL import Image


ARTICLE = 77
ANALYSIS_SEED = 77_001
PLOT_SEED = 20_260_777
SNAPSHOT_DATE = "2026-08-23"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n", float_format="%.10g")


def git_value(root: Path, *arguments: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), *arguments],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout.strip()


def policy_sources() -> pd.DataFrame:
    rows = [
        ("LIFECYCLE", "GigaScience", "Metagenomic data object model", "https://doi.org/10.1093/gigascience/gix047", "2017 article", "No"),
        ("NCBI-SRA-START", "NCBI", "Raw-read scope and human metagenomes", "https://www.ncbi.nlm.nih.gov/sra/docs/submit/", "Live documentation", "Yes"),
        ("NCBI-SRA-PORTAL", "NCBI", "Portal, accepted objects and submission IDs", "https://www.ncbi.nlm.nih.gov/sra/docs/submitportal/", "Live documentation", "Yes"),
        ("NCBI-SRA-META", "NCBI", "Study-Sample-Experiment-Run model", "https://www.ncbi.nlm.nih.gov/sra/docs/submitmeta/", "Live documentation", "Yes"),
        ("NCBI-SRA-FORMAT", "NCBI", "Read formats and quality scores", "https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/", "Live documentation", "Yes"),
        ("NCBI-BIO", "NCBI", "BioProject and BioSample prerequisites", "https://www.ncbi.nlm.nih.gov/sra/docs/submitbio/", "Live documentation", "Yes"),
        ("NCBI-MIMAG", "NCBI", "MIMAG BioSample package", "https://www.ncbi.nlm.nih.gov/biosample/docs/packages/MIMAG.5.0/", "Package 6.0 displayed", "Yes"),
        ("NCBI-MIUVIG", "NCBI", "MIUVIG BioSample package", "https://www.ncbi.nlm.nih.gov/biosample/docs/packages/MIUVIG.5.0/", "Package 6.0 displayed", "Yes"),
        ("ENA-READS", "ENA", "Raw-read Webin-CLI submission", "https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html", "Live documentation", "Yes"),
        ("ENA-PRIMARY", "ENA", "Primary metagenome assembly", "https://ena-docs.readthedocs.io/en/latest/submit/assembly/metagenome/primary.html", "Live documentation", "Yes"),
        ("ENA-MAG", "ENA", "MAG submission", "https://ena-docs.readthedocs.io/en/latest/submit/assembly/metagenome/mag.html", "Live documentation", "Yes"),
        ("ZENODO-RECORD", "Zenodo", "Records, DOI and immutability", "https://help.zenodo.org/docs/deposit/about-records/", "Live documentation", "Yes"),
        ("ZENODO-VERSION", "Zenodo", "Version-specific records and DOIs", "https://help.zenodo.org/docs/deposit/manage-versions/", "Live documentation", "Yes"),
        ("GITHUB-CITATION", "GitHub", "CITATION.cff", "https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files", "Live documentation", "Yes"),
        ("GITHUB-ZENODO", "GitHub", "Release archival with Zenodo", "https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content", "Live documentation", "Yes"),
        ("GITHUB-IMMUTABLE", "GitHub", "Immutable releases", "https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases", "Live documentation", "Yes"),
        ("DOCKER-DIGEST", "Docker", "Image digest semantics", "https://docs.docker.com/dhi/explore/security-concepts/digests/", "Live documentation", "Yes"),
        ("NATURE-DAS", "Springer Nature", "Data Availability statement", "https://support.nature.com/en/support/solutions/articles/6000237611-write-a-data-availability-statement-for-a-paper", "Live documentation", "Yes"),
    ]
    frame = pd.DataFrame(
        rows,
        columns=["SourceID", "Provider", "Scope", "URL", "VersionBoundary", "RecheckAtSubmission"],
    )
    frame.insert(4, "RetrievedDate", SNAPSHOT_DATE)
    return frame


def policy_assertions() -> pd.DataFrame:
    rows = [
        ("SRA-01", "Raw reads", "SRA accepts raw or minimally processed sequence reads, not assembled contigs.", "NCBI-SRA-START; NCBI-SRA-PORTAL"),
        ("SRA-02", "Quality", "Raw read archives require per-base quality information; FASTA alone is insufficient.", "NCBI-SRA-FORMAT"),
        ("SRA-03", "Object model", "Study/BioProject, Sample/BioSample, Experiment and Run are distinct linked objects.", "NCBI-SRA-META"),
        ("SRA-04", "Paired data", "Paired files from one library belong to one Run; different samples are not grouped into one Run.", "NCBI-SRA-META"),
        ("SRA-05", "Identifier", "A submission tracking ID such as SUB is not the accession cited in a paper.", "NCBI-SRA-PORTAL"),
        ("SRA-06", "Governance", "Human metagenomes require consent, access and host-sequence review before release.", "NCBI-SRA-START"),
        ("ENA-01", "Read submission", "A read submission requires pre-registered study and sample records plus experiment/library metadata.", "ENA-READS"),
        ("ENA-02", "Assembly linkage", "A primary metagenome assembly should be linked to its raw runs.", "ENA-PRIMARY"),
        ("ENA-03", "MAG linkage", "A MAG uses a derived sample record and should retain links to lower-level reads and assemblies.", "ENA-MAG"),
        ("ENA-04", "Validation", "Webin-CLI validation and submission are separate operations.", "ENA-PRIMARY; ENA-MAG"),
        ("ENA-05", "Identifier", "Internal submission identifiers are not the final public sequence accessions.", "ENA-PRIMARY; ENA-MAG"),
        ("DOI-01", "Repository", "A Git branch or mutable tag is not a preservation DOI.", "GITHUB-ZENODO; ZENODO-RECORD"),
        ("DOI-02", "Versioning", "A new Zenodo version has new files and a version-specific persistent identifier linked to prior versions.", "ZENODO-VERSION"),
        ("DOI-03", "Citation", "CITATION.cff exposes preferred citation metadata but does not itself preserve a release.", "GITHUB-CITATION; GITHUB-ZENODO"),
        ("OCI-01", "Container", "An image tag is mutable while a digest identifies content.", "DOCKER-DIGEST"),
        ("OCI-02", "Container", "A multi-platform index digest and a platform manifest digest identify different objects.", "DOCKER-DIGEST"),
        ("DB-01", "Database", "A reproducible database record needs release, artifact, checksum and source, not a database name alone.", "LIFECYCLE"),
        ("DAS-01", "Availability", "A Data Availability statement covers both generated and reused data and gives persistent identifiers or restrictions.", "NATURE-DAS"),
        ("DAS-02", "Availability", "Code availability is recorded separately from scientific data availability.", "NATURE-DAS; GITHUB-ZENODO"),
        ("ETHICS-01", "Access", "FAIR does not require sensitive human files to be openly downloadable; metadata and access procedures can be public.", "NCBI-SRA-START; NATURE-DAS"),
        ("OWN-01", "Ownership", "Existing public tutorial data must be cited and reused, not resubmitted as a new investigator-owned study.", "NCBI-SRA-PORTAL"),
        ("LINK-01", "Lineage", "Derived tables, catalogs, MAGs and UViGs must point back to samples/runs and analysis versions.", "LIFECYCLE; ENA-PRIMARY; ENA-MAG"),
    ]
    return pd.DataFrame(rows, columns=["AssertionID", "Domain", "Assertion", "SourceIDs"])


def routing_matrix() -> pd.DataFrame:
    rows = [
        ("Raw reads", "SRA / ENA", "Run", "Study + Sample + Experiment", "FASTQ/BAM/CRAM with quality", "Existing third-party accessions; cite, do not resubmit"),
        ("Study description", "BioProject / ENA Study", "Study", "Project scope and contributors", "Metadata record", "Use existing PRJEB52977 for the tutorial inputs"),
        ("Source material", "BioSample / ENA Sample", "Sample", "Collection and source metadata", "One biological material record", "Use existing SAMEA records for the tutorial inputs"),
        ("Library and instrument", "SRA / ENA", "Experiment", "Sample + library + platform", "Metadata record", "Resolve existing experiment accessions before manuscript freeze"),
        ("Primary metagenome assembly", "GenBank / ENA", "Assembly analysis", "Raw runs + environmental sample", "FASTA + manifest", "Local tutorial assembly is not submitted"),
        ("MAG", "GenBank / ENA", "Genome assembly + derived BioSample", "Primary assembly + source runs", "One qualified nonredundant MAG", "Blocked by ownership, metadata and manual review"),
        ("UViG / vOTU", "INSDC analysis archive + durable data repository", "Virus sequence + MIUViG record", "Source runs + assembly + detection", "FASTA + per-UViG metadata", "Fixture is incomplete and not a submitter-owned study"),
        ("Gene catalog", "Zenodo / institutional repository", "Versioned dataset", "Source runs + assembly + clustering", "Representative FNA/FAA + membership + annotation", "Sequence payload is not staged in this article"),
        ("Functional annotation", "Zenodo / institutional repository", "Versioned dataset", "Gene catalog + database manifest", "Long table + dictionary", "Local only; needs a release DOI"),
        ("Taxonomic/functional profiles", "Zenodo / institutional repository", "Versioned dataset", "Samples + tool/database versions", "Machine-readable tables", "Local only; needs a release DOI"),
        ("Figure source tables", "Zenodo / institutional repository", "Versioned dataset", "Claim + script + upstream artifact", "Tidy source tables + README", "Local checksum bundle, not published"),
        ("Workflow source", "GitHub + Zenodo", "Tagged release", "Commit + tests + license", "Repository archive", "No public release tag or DOI recorded"),
        ("Container image", "OCI registry", "Image manifest", "Build recipe + source commit", "Index/platform digest", "No digest recorded; tag alone would be insufficient"),
        ("Database manifest", "Data repository + code release", "Version ledger", "Every reference database", "Release + artifact + checksum + source", "Three database locks are complete locally"),
        ("Restricted human data", "Controlled-access repository", "Study/Sample/Run metadata + access record", "Consent and governance", "Public metadata; controlled files", "Not applicable to this public non-human mock example"),
    ]
    frame = pd.DataFrame(
        rows,
        columns=["AssetClass", "Destination", "ArchiveObject", "RequiredLinks", "DepositUnit", "WorkedExampleState"],
    )
    frame["RouteOrder"] = range(1, len(frame) + 1)
    return frame


def object_relationships() -> pd.DataFrame:
    rows = [
        ("STUDY-01", "Study", "PRJEB52977", "None", "EXISTING_THIRD_PARTY", "Public source project; cite and reuse"),
        ("SAMPLE-01", "Sample", "SAMEA14435832", "PRJEB52977", "EXISTING_THIRD_PARTY", "Source of ERR9765746"),
        ("SAMPLE-02", "Sample", "SAMEA14435833", "PRJEB52977", "EXISTING_THIRD_PARTY", "Source of ERR9765747"),
        ("EXPERIMENT-01", "Experiment", "EXISTING_NOT_CAPTURED", "SAMEA14435832", "RESOLVE_BEFORE_RELEASE", "Resolve from ENA; do not infer an ERX identifier"),
        ("EXPERIMENT-02", "Experiment", "EXISTING_NOT_CAPTURED", "SAMEA14435833", "RESOLVE_BEFORE_RELEASE", "Resolve from ENA; do not infer an ERX identifier"),
        ("RUN-01", "Run", "ERR9765746", "SAMEA14435832", "EXISTING_THIRD_PARTY", "Paired FASTQ, ENA MD5 recorded"),
        ("RUN-02", "Run", "ERR9765747", "SAMEA14435833", "EXISTING_THIRD_PARTY", "Paired FASTQ, ENA MD5 recorded"),
        ("ANALYSIS-ASM", "Primary assembly", "LOCAL_ARTICLE30_ASSEMBLY", "ERR9765746;ERR9765747", "LOCAL_ONLY", "Not an INSDC accession"),
        ("ANALYSIS-CATALOG", "Gene catalog", "LOCAL_ARTICLE34_CATALOG", "ERR9765746;ERR9765747", "LOCAL_ONLY", "93,782 representative genes; files require deposition"),
        ("ANALYSIS-MAG", "MAG collection", "PENDING_NOT_INVENTED", "LOCAL_ARTICLE30_ASSEMBLY", "BLOCKED", "No accession: public tutorial mock is not submitter-owned"),
        ("ANALYSIS-UVIG", "UViG collection", "PENDING_NOT_INVENTED", "INDEPENDENT_CHECKV_FIXTURE", "BLOCKED", "Independent fixture; incomplete MIUViG metadata"),
        ("RELEASE-DATA", "Derived-data release", "PENDING_NOT_INVENTED", "LOCAL_ARTICLE34_CATALOG", "BLOCKED", "No Zenodo DOI reserved or published"),
        ("RELEASE-CODE", "Code release", "PENDING_NOT_INVENTED", "LOCAL_GIT_COMMIT", "BLOCKED", "No public tag or archive DOI recorded"),
        ("RELEASE-IMAGE", "Container image", "PENDING_NOT_INVENTED", "LOCAL_DOCKERFILE", "BLOCKED", "No OCI digest recorded"),
    ]
    return pd.DataFrame(rows, columns=["ObjectID", "ObjectType", "Identifier", "Parent", "IdentifierState", "EvidenceOrAction"])


def identifier_registry(git_commit: str) -> pd.DataFrame:
    rows = [
        ("BioProject / Study", "PRJEB52977", "EXISTING_THIRD_PARTY", "ENA", "Cite source; never claim ownership"),
        ("BioSample / Sample 1", "SAMEA14435832", "EXISTING_THIRD_PARTY", "ENA", "Links to ERR9765746"),
        ("BioSample / Sample 2", "SAMEA14435833", "EXISTING_THIRD_PARTY", "ENA", "Links to ERR9765747"),
        ("Run 1", "ERR9765746", "EXISTING_THIRD_PARTY", "ENA", "Paired FASTQ source"),
        ("Run 2", "ERR9765747", "EXISTING_THIRD_PARTY", "ENA", "Paired FASTQ source"),
        ("Experiment accessions", "EXISTING_NOT_CAPTURED", "RESOLVE_BEFORE_RELEASE", "ENA", "Resolve; do not derive names from run strings"),
        ("Primary assembly accession", "PENDING_NOT_INVENTED", "BLOCKED", "INSDC", "Submitter ownership and metadata required"),
        ("MAG accessions", "PENDING_NOT_INVENTED", "BLOCKED", "INSDC", "One approved MAG submission unit at a time"),
        ("UViG accessions", "PENDING_NOT_INVENTED", "BLOCKED", "INSDC / repository", "MIUViG and ownership gates unresolved"),
        ("Derived-data DOI", "PENDING_NOT_INVENTED", "BLOCKED", "Zenodo", "Reserve or publish only after payload approval"),
        ("Code release tag", "PENDING_NOT_INVENTED", "BLOCKED", "GitHub", "Clean, test and tag a public release"),
        ("Code version DOI", "PENDING_NOT_INVENTED", "BLOCKED", "Zenodo", "Archive the exact GitHub release"),
        ("Local Git commit", git_commit, "LOCAL_ONLY", "Git", "Provenance only; not a public release identifier"),
        ("OCI index digest", "PENDING_NOT_INVENTED", "BLOCKED", "OCI registry", "Push and inspect the immutable index digest"),
        ("OCI platform digest", "PENDING_NOT_INVENTED", "BLOCKED", "OCI registry", "Record the tested platform manifest digest"),
    ]
    frame = pd.DataFrame(rows, columns=["IdentifierType", "Value", "State", "Authority", "UseOrNextAction"])
    frame["StateCode"] = frame["State"].map(
        {"BLOCKED": 0, "RESOLVE_BEFORE_RELEASE": 1, "LOCAL_ONLY": 2, "EXISTING_THIRD_PARTY": 3}
    )
    return frame


def release_readiness() -> pd.DataFrame:
    rows = [
        ("Reused raw reads", "Existing third-party", "Cite PRJEB52977, ERR9765746 and ERR9765747; no new submission", False),
        ("Study and sample metadata", "Existing third-party", "Cite PRJEB52977 and two SAMEA accessions; resolve experiment accessions", False),
        ("Primary assembly", "Local only", "Package manifest and FASTA only after ownership and release scope are approved", False),
        ("MAG collection", "Blocked", "Manual sign-off, ownership, taxonomy names, BioSamples and source metadata are unresolved", False),
        ("UViG collection", "Blocked", "Independent fixture is not a study deposit and MIUViG has two missing fields", False),
        ("Gene catalog", "Local only", "Summary records 93,782 genes, but representative files are not staged here", False),
        ("Annotations and profiles", "Local only", "Deposit long tables, dictionaries and database links as a versioned dataset", False),
        ("Figure source tables", "Local package ready", "Checksum-covered tables can be bundled after manuscript scope freezes", False),
        ("Workflow code", "Blocked", "Working tree is not a clean public release and no tag is recorded", False),
        ("Code preservation DOI", "Blocked", "No GitHub release has been archived by Zenodo", False),
        ("Container image", "Missing", "No registry, index digest or tested platform digest is recorded", False),
        ("Database manifest", "Local package ready", "GTDB R232, CheckM2 v3 and GUNC ProGenomes 2.1 have checksums", False),
        ("Data Availability", "Draft ready", "Current statement names reused data and honestly marks generated deposits pending", False),
        ("Code Availability", "Draft ready", "Current statement separates local commit evidence from a future public DOI", False),
        ("External release", "Blocked", "No external-ready claim is made by this tutorial packet", False),
    ]
    frame = pd.DataFrame(rows, columns=["ReleaseObject", "Status", "EvidenceOrNextAction", "ExternalReady"])
    frame["StatusCode"] = frame["Status"].map(
        {"Missing": 0, "Blocked": 0, "Local only": 1, "Draft ready": 1, "Local package ready": 2, "Existing third-party": 3}
    )
    return frame


def release_gates(code_clean: bool) -> pd.DataFrame:
    rows = [
        (1, "Release scope frozen", "Review", "Decide which assemblies, catalogs, profiles and figure tables are release payloads", "Corresponding author"),
        (2, "Ownership and non-duplicate", "Blocked", "Tutorial source data and mock MAGs are third-party public examples", "Principal investigator"),
        (3, "Human-data governance", "Not applicable", "Worked example uses public non-human mock-community reads", "Data steward"),
        (4, "Stable sample/run crosswalk", "Pass", "PRJEB52977, SAMEA14435832/3 and ERR9765746/7 are recorded", "Data steward"),
        (5, "Experiment accessions resolved", "Blocked", "Existing ERX identifiers were not copied and are not inferred", "Data steward"),
        (6, "MAG manual review", "Blocked", "Article 49 investigator sign-off remains pending", "Genome analyst"),
        (7, "MAG source metadata", "Blocked", "BioProject, per-MAG BioSamples and mandatory package fields are absent", "Genome analyst"),
        (8, "UViG MIUViG record", "Blocked", "Assembly software and predicted genome type are missing", "Virus analyst"),
        (9, "Payload checksums", "Pass", "This local packet and source artifacts have SHA-256 coverage", "Data steward"),
        (10, "Repository format validation", "Blocked", "No Webin-CLI or NCBI submission validation report exists", "Data steward"),
        (11, "Accession resolvability", "Blocked", "No new sequence accession or derived-data DOI exists", "Data steward"),
        (12, "Clean tested code tree", "Pass" if code_clean else "Blocked", "Local Git state audited; a public release must be clean and tested", "Pipeline lead"),
        (13, "Tagged code release", "Blocked", "No public release tag is recorded", "Pipeline lead"),
        (14, "Code archive DOI", "Blocked", "No version DOI is recorded", "Pipeline lead"),
        (15, "Container digest", "Blocked", "No OCI index/platform digest is recorded", "Pipeline lead"),
        (16, "Database manifest", "Pass", "Release, artifact, checksum and source URL are recorded for three databases", "Pipeline lead"),
        (17, "License and citation files", "Review", "LICENSE exists; CITATION.cff is absent", "Corresponding author"),
        (18, "Availability statements", "Review", "Draft clauses exist but pending identifiers cannot be published as resolved", "Corresponding author"),
    ]
    frame = pd.DataFrame(rows, columns=["GateOrder", "Gate", "Status", "EvidenceOrNextAction", "Owner"])
    frame["StatusCode"] = frame["Status"].map({"Blocked": 0, "Review": 1, "Not applicable": 2, "Pass": 3})
    return frame


def container_ledger(root: Path) -> pd.DataFrame:
    rows = [
        ("Build recipe", "Dockerfile", "LOCAL_ONLY", "File exists" if (root / "Dockerfile").is_file() else "Missing"),
        ("Image tag", "PENDING_NOT_INVENTED", "BLOCKED", "A tag may move and is not sufficient for citation"),
        ("OCI registry", "PENDING_NOT_INVENTED", "BLOCKED", "No public registry location recorded"),
        ("Multi-platform index digest", "PENDING_NOT_INVENTED", "BLOCKED", "Record sha256 after push"),
        ("Tested platform", "linux/amd64", "PLANNED", "Declare the platform used for manuscript results"),
        ("Platform manifest digest", "PENDING_NOT_INVENTED", "BLOCKED", "Record sha256 for the tested platform"),
        ("Source commit link", "LOCAL_GIT_COMMIT", "LOCAL_ONLY", "Replace with public tagged commit URL at release"),
        ("Runtime verification", "PENDING_NOT_INVENTED", "BLOCKED", "Run smoke test from digest-pinned image"),
    ]
    frame = pd.DataFrame(rows, columns=["ContainerField", "Value", "State", "EvidenceOrAction"])
    frame["StateCode"] = frame["State"].map({"BLOCKED": 0, "PLANNED": 1, "LOCAL_ONLY": 2})
    return frame


def availability_components() -> tuple[pd.DataFrame, pd.DataFrame]:
    data_rows = [
        ("Reused study", "Complete", "PRJEB52977", "State that these are third-party public data"),
        ("Reused samples", "Complete", "SAMEA14435832; SAMEA14435833", "Map samples to runs"),
        ("Reused runs", "Complete", "ERR9765746; ERR9765747", "Give archive and accession"),
        ("Primary assembly", "Pending", "PENDING_NOT_INVENTED", "Deposit only with ownership and source links"),
        ("MAGs", "Blocked", "PENDING_NOT_INVENTED", "Do not cite local MAG names as accessions"),
        ("UViGs", "Blocked", "PENDING_NOT_INVENTED", "Complete MIUViG and ownership review first"),
        ("Gene catalog", "Pending", "PENDING_NOT_INVENTED", "Deposit representative FNA/FAA, membership and annotation"),
        ("Figure source data", "Local package ready", "LOCAL_CHECKSUM_BUNDLE", "Publish a DOI-bearing release before final statement"),
        ("Restrictions", "Complete", "No controlled human data in worked example", "For human studies, state repository and access procedure"),
    ]
    code_rows = [
        ("Repository URL", "Blocked", "PENDING_NOT_INVENTED", "Give a public repository URL"),
        ("Release tag", "Blocked", "PENDING_NOT_INVENTED", "Use the exact manuscript release"),
        ("Commit", "Local only", "LOCAL_GIT_COMMIT", "Retain for provenance; link publicly at release"),
        ("Version DOI", "Blocked", "PENDING_NOT_INVENTED", "Archive the tagged release"),
        ("License", "Complete", "MIT code / CC BY 4.0 content", "Retain repository license notices"),
        ("CITATION.cff", "Missing", "PENDING_NOT_INVENTED", "Add preferred citation metadata before release"),
        ("Container digest", "Blocked", "PENDING_NOT_INVENTED", "Give registry plus index and tested-platform digest"),
        ("Database manifest", "Complete", "database-manifest-public.tsv", "Release with code and results"),
    ]
    return (
        pd.DataFrame(data_rows, columns=["Component", "Status", "IdentifierOrEvidence", "RequiredWording"]),
        pd.DataFrame(code_rows, columns=["Component", "Status", "IdentifierOrEvidence", "RequiredWording"]),
    )


def package_index() -> pd.DataFrame:
    rows = [
        ("README", "README.md", "Explain scope, versions, licenses and reuse"),
        ("Data dictionary", "metadata/data-dictionary.tsv", "Define every field, unit and missing-value code"),
        ("Sample crosswalk", "metadata/sample-run-crosswalk.tsv", "Link internal IDs to BioSample and Run"),
        ("Raw-read manifest", "manifests/raw-reads.tsv", "URL, bytes and archive checksum per mate"),
        ("Assembly manifest", "manifests/assemblies.tsv", "Source runs, assembler, parameters and checksum"),
        ("MAG manifest", "manifests/mags.tsv", "MIMAG, taxonomy, source assembly and accession"),
        ("UViG manifest", "manifests/uvigs.tsv", "MIUViG, vOTU membership and accession"),
        ("Gene catalog manifest", "manifests/gene-catalog.tsv", "FNA/FAA/membership/annotation hashes"),
        ("Database manifest", "manifests/databases.tsv", "Release, artifact, checksum and source"),
        ("Workflow parameters", "workflow/params.json", "Non-default parameters and seed"),
        ("Software lock", "workflow/software-versions.tsv", "Workflow, tool and plugin identities"),
        ("Container ledger", "workflow/container-ledger.tsv", "Registry and immutable digests"),
        ("Figure sources", "results/figure-source/", "One tidy table per final figure/panel"),
        ("Analysis outputs", "results/tables/", "Machine-readable derived tables"),
        ("Checksums", "SHA256SUMS", "Cover every released payload file"),
        ("CITATION", "CITATION.cff", "Preferred citation and version"),
        ("License", "LICENSES/", "Content, code and third-party boundaries"),
        ("Availability statements", "manuscript/availability-statements.md", "Exact identifiers and access conditions"),
    ]
    frame = pd.DataFrame(rows, columns=["PackageObject", "SuggestedPath", "Purpose"])
    frame["Order"] = range(1, len(frame) + 1)
    return frame


def artifact_sources(root: Path) -> list[tuple[str, str, str]]:
    return [
        ("A30-RAW", "data/small/30-short-read-assembly-frozen/source-manifest.tsv", "Public ENA raw-read manifest"),
        ("A34-CATALOG", "data/small/34-nonredundant-gene-catalog-frozen/catalog-summary.tsv", "Gene catalog summary and logical sequence hashes"),
        ("A34-RUN", "data/small/34-nonredundant-gene-catalog-frozen/run-summary.json", "Gene catalog run contract"),
        ("A49-GATES", "data/small/49-mag-curation-submission-frozen/submission-readiness-checklist.tsv", "MAG submission gates"),
        ("A49-MIMAG", "data/small/49-mag-curation-submission-frozen/mimag-quality-supplement.tsv", "Per-MAG quality metadata"),
        ("A54-LINEAGE", "data/small/54-virus-discovery-quality-frozen/input-lineage.tsv", "Virus fixture lineage"),
        ("A54-MIUVIG", "data/small/54-virus-discovery-quality-frozen/miuvig-source-assertions.tsv", "MIUViG source assertions"),
        ("A74-RELEASE", "data/small/74-nfcore-mag-workflows-frozen/release-lock.tsv", "Workflow release locks"),
        ("A74-DATABASE", "data/small/74-nfcore-mag-workflows-frozen/database-lock.tsv", "Database locks; local paths excluded downstream"),
        ("A74-RUNTIME", "data/small/74-nfcore-mag-workflows-frozen/runtime-environment.tsv", "Runtime evidence boundary"),
        ("A75-TRACE", "data/small/75-paper-figure-organization-frozen/result-traceability-ledger.tsv", "Claim-to-result traceability"),
        ("A75-VERSION", "data/small/75-paper-figure-organization-frozen/version-ledger-example.tsv", "Public version-ledger example"),
        ("A76-READY", "data/small/76-reporting-standards-frozen/submission-readiness.tsv", "Reporting-standard readiness"),
        ("RECIPE", "Dockerfile", "Local container build recipe"),
        ("LICENSE", "LICENSE", "Repository licensing boundary"),
    ]


def source_artifact_manifest(root: Path) -> pd.DataFrame:
    rows = []
    for artifact_id, relative, role in artifact_sources(root):
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append((artifact_id, relative, path.stat().st_size, sha256(path), role))
    return pd.DataFrame(rows, columns=["ArtifactID", "RelativePath", "Bytes", "SHA256", "Role"])


def release_artifact_manifest(root: Path, raw: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    primary = catalog.loc[catalog["CatalogID"].eq("megahit-mix-primary")].iloc[0]
    rows = []
    for row in raw.itertuples(index=False):
        rows.append(
            (
                f"{row.RunAccession}-{row.Mate}",
                "Raw read",
                "ENA remote",
                row.RunAccession,
                "MD5",
                row.ENAReportedMD5,
                int(row.ENABytes),
                "EXISTING_THIRD_PARTY",
            )
        )
    rows.extend(
        [
            ("megahit-mix-primary.fna", "Gene catalog nucleotide representatives", "Not staged", "PENDING_NOT_INVENTED", "SHA256", primary["RepresentativeFNA_SHA256"], "NOT_STAGED", "LOCAL_LOGICAL_RECORD"),
            ("megahit-mix-primary.faa", "Gene catalog protein representatives", "Not staged", "PENDING_NOT_INVENTED", "SHA256", primary["RepresentativeFAA_SHA256"], "NOT_STAGED", "LOCAL_LOGICAL_RECORD"),
            ("MAG-collection", "MAG sequences", "Not staged", "PENDING_NOT_INVENTED", "None", "PENDING_NOT_INVENTED", "NOT_STAGED", "BLOCKED"),
            ("UViG-collection", "UViG sequences", "Not staged", "PENDING_NOT_INVENTED", "None", "PENDING_NOT_INVENTED", "NOT_STAGED", "BLOCKED"),
            ("Figure-source-tables", "Derived tables", "Article 77 frozen packet", "LOCAL_ONLY", "SHA256", "See file-checksums.sha256", "LOCAL", "LOCAL_PACKAGE_READY"),
            ("Workflow-release", "Code archive", "Not published", "PENDING_NOT_INVENTED", "None", "PENDING_NOT_INVENTED", "NOT_PUBLISHED", "BLOCKED"),
            ("Container-image", "OCI image", "Not published", "PENDING_NOT_INVENTED", "SHA256", "PENDING_NOT_INVENTED", "NOT_PUBLISHED", "BLOCKED"),
        ]
    )
    return pd.DataFrame(
        rows,
        columns=["Artifact", "AssetClass", "LocationClass", "PersistentIdentifier", "ChecksumType", "Checksum", "Bytes", "State"],
    )


def availability_text(git_commit: str) -> str:
    return f"""# Availability statement drafts

## Worked-example Data Availability

The raw sequencing data reused in this tutorial are available from the European Nucleotide Archive under Study PRJEB52977, Samples SAMEA14435832 and SAMEA14435833, and Runs ERR9765746 and ERR9765747. The local assembly, MAG, UViG, gene-catalog and figure-source artifacts in this worked example have not been deposited as new research outputs and have no new accession or DOI. They must not be described as publicly released datasets.

## Production Data Availability template

Raw sequence reads are available from [SRA_OR_ENA] under BioProject/Study [PROJECT_ACCESSION], BioSamples [BIOSAMPLE_ACCESSIONS] and Runs [RUN_ACCESSIONS]. The primary metagenome assemblies and MAG/UViG records are available under [ASSEMBLY_AND_GENOME_ACCESSIONS]. Version [DATA_RELEASE] of the gene catalog, annotations, abundance tables and figure-source data is archived at [DATA_DOI]. Controlled files, if any, are available through [CONTROLLED_REPOSITORY] under [ACCESS_PROCEDURE]; public metadata remain available at [METADATA_ACCESSION].

## Worked-example Code Availability

The local analysis state was audited at Git commit `{git_commit}` but is not a clean, tagged public release. No public repository URL, archived version DOI or OCI image digest is claimed. The database ledger records GTDB R232, CheckM2 database v3 and GUNC ProGenomes 2.1 with source checksums.

## Production Code Availability template

Analysis code is available at [PUBLIC_REPOSITORY_URL], release [RELEASE_TAG], commit [COMMIT_SHA], and archived as version DOI [CODE_VERSION_DOI]. The tested container is [REGISTRY/IMAGE] pinned by multi-platform index digest [OCI_INDEX_DIGEST] and tested-platform manifest digest [OCI_PLATFORM_DIGEST] for [PLATFORM]. Workflow parameters, software locks and database manifests are included in the archived release.
"""


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    evidence = args.evidence_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    download_manifest = json.loads((evidence / "download-manifest.json").read_text(encoding="utf-8"))
    xml = evidence / "PMC5737865-fulltext.xml"
    anchor = evidence / "gix047fig2.jpg"
    if not xml.is_file() or not anchor.is_file():
        raise FileNotFoundError("Article 77 primary evidence is incomplete")
    shutil.copy2(anchor, output / "tenhoopen-figure2-original.jpg")

    git_commit = git_value(root, "rev-parse", "HEAD")
    project_status = git_value(root, "status", "--porcelain", "--", ".")
    code_clean = not bool(project_status)

    raw = pd.read_csv(root / "data/small/30-short-read-assembly-frozen/source-manifest.tsv", sep="\t")
    catalog = pd.read_csv(root / "data/small/34-nonredundant-gene-catalog-frozen/catalog-summary.tsv", sep="\t")
    mag_gates = pd.read_csv(root / "data/small/49-mag-curation-submission-frozen/submission-readiness-checklist.tsv", sep="\t")
    databases = pd.read_csv(root / "data/small/74-nfcore-mag-workflows-frozen/database-lock.tsv", sep="\t")
    database_public = databases.drop(columns=["LocalPath"])
    database_public.insert(2, "SnapshotDate", SNAPSHOT_DATE)

    policies = policy_sources()
    assertions = policy_assertions()
    routes = routing_matrix()
    objects = object_relationships()
    identifiers = identifier_registry(git_commit)
    readiness = release_readiness()
    gates = release_gates(code_clean)
    containers = container_ledger(root)
    data_components, code_components = availability_components()
    packet = package_index()
    source_artifacts = source_artifact_manifest(root)
    release_artifacts = release_artifact_manifest(root, raw, catalog)

    tables = {
        "policy-source-registry.tsv": policies,
        "policy-assertions.tsv": assertions,
        "repository-routing-matrix.tsv": routes,
        "object-relationship.tsv": objects,
        "identifier-registry.tsv": identifiers,
        "release-readiness.tsv": readiness,
        "release-gate-ledger.tsv": gates,
        "container-ledger.tsv": containers,
        "database-manifest-public.tsv": database_public,
        "data-availability-components.tsv": data_components,
        "code-availability-components.tsv": code_components,
        "package-index.tsv": packet,
        "source-artifact-manifest.tsv": source_artifacts,
        "release-artifact-manifest.tsv": release_artifacts,
    }
    for filename, frame in tables.items():
        write_tsv(frame, output / filename)

    (output / "availability-statements.md").write_text(
        availability_text(git_commit) + "\n", encoding="utf-8"
    )
    source_manifest = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "records": [
            {
                "File": xml.name,
                "Bytes": xml.stat().st_size,
                "SHA256": sha256(xml),
                "Identity": "ten Hoopen et al. 2017 full-text XML",
            },
            {
                "File": anchor.name,
                "Bytes": anchor.stat().st_size,
                "SHA256": sha256(anchor),
                "Identity": "ten Hoopen et al. 2017 Figure 2 selected from Europe PMC archive",
            },
        ],
        "anchor": {
            "file": "tenhoopen-figure2-original.jpg",
            "source_member": "gix047fig2.jpg",
            "license": "CC BY 4.0",
            "attribution": "ten Hoopen et al., GigaScience 2017, DOI 10.1093/gigascience/gix047",
            "sha256": sha256(anchor),
            "dimensions": list(Image.open(anchor).size),
        },
        "policy_snapshot": {
            "date": SNAPSHOT_DATE,
            "mutable_pages": len(policies) - 1,
            "rule": "Recheck live repository and journal documentation immediately before submission.",
        },
        "download_manifest": download_manifest,
    }
    (output / "source-manifest.json").write_text(
        json.dumps(source_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    primary = catalog.loc[catalog["CatalogID"].eq("megahit-mix-primary")].iloc[0]
    metrics = {
        "article": ARTICLE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "snapshot_date": SNAPSHOT_DATE,
        "policy_sources": len(policies),
        "policy_assertions": len(assertions),
        "routing_assets": len(routes),
        "object_records": len(objects),
        "identifier_records": len(identifiers),
        "existing_third_party_identifiers": int(identifiers["State"].eq("EXISTING_THIRD_PARTY").sum()),
        "pending_or_blocked_identifiers": int(identifiers["State"].isin(["BLOCKED", "RESOLVE_BEFORE_RELEASE"]).sum()),
        "external_ready_records": int(readiness["ExternalReady"].sum()),
        "release_gates": len(gates),
        "blocked_release_gates": int(gates["Status"].eq("Blocked").sum()),
        "passed_release_gates": int(gates["Status"].eq("Pass").sum()),
        "raw_runs": int(raw["RunAccession"].nunique()),
        "raw_fastq_files": len(raw),
        "raw_bytes": int(raw["ENABytes"].sum()),
        "primary_catalog_genes": int(primary["CatalogGenes"]),
        "article49_blocked_gates": int(mag_gates["Status"].eq("BLOCKED").sum()),
        "database_records": len(database_public),
        "source_artifacts": len(source_artifacts),
        "release_artifacts": len(release_artifacts),
        "code_worktree_clean": code_clean,
        "new_accessions_or_dois_claimed": 0,
    }
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    methods = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "seeds": {"analysis": ANALYSIS_SEED, "plot": PLOT_SEED},
        "object_model": "Study/BioProject -> Sample/BioSample -> Experiment -> Run -> Analysis/derived entities",
        "repository_policy": "SRA/ENA for primary reads and eligible sequence assemblies; durable DOI repository for processed tables and catalogs; GitHub release plus Zenodo DOI for code; OCI registry plus digest for images.",
        "database_policy": "Release, artifact, checksum and source URL; local filesystem paths excluded from public tables.",
        "identifier_policy": "No accession, DOI, tag or digest is invented. PENDING_NOT_INVENTED is a blocking token, never a citeable identifier.",
        "ownership_policy": "Third-party public tutorial inputs are cited and reused, not resubmitted as investigator-owned outputs.",
        "availability_policy": "Data and Code Availability statements are separate and must cover reused and generated assets.",
        "external_action": "None; this script creates a local readiness packet only.",
        "python": platform.python_version(),
        "pandas": pd.__version__,
    }
    (output / "methods-contract.json").write_text(
        json.dumps(methods, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    notice = """Article 77 data notice

This directory is a local, checksum-covered release-readiness example. It is not an SRA/ENA, Zenodo, GitHub or OCI publication and confers no new accession, DOI, release tag or image digest.

Raw inputs PRJEB52977, SAMEA14435832, SAMEA14435833, ERR9765746 and ERR9765747 are third-party public ENA records and must be cited to their original study. They must not be resubmitted or presented as investigator-owned data.

The anchor image is Figure 2 from ten Hoopen et al. (2017), DOI 10.1093/gigascience/gix047, reused under CC BY 4.0 with attribution. Other referenced repository documentation remains subject to its providers' terms. Local tutorial prose is CC BY 4.0 and code is MIT under the repository LICENSE.
"""
    (output / "data-NOTICE.txt").write_text(notice, encoding="utf-8")
    print(f"prepared\t{output}\t{len(tables)} tables")


if __name__ == "__main__":
    main()
