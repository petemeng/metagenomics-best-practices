#!/usr/bin/env python3
"""Prepare the frozen public-resource and environmental-MAG tables for Article 73."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


ARTICLE = 73
ANALYSIS_SEED = 73_001
PLOT_SEED = 20_260_773
SNAPSHOT_DATE = "2026-08-23"


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


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
        float_format="%.10g",
    )


def load_pages(cache: Path, prefix: str, pages: tuple[int, ...]) -> list[dict]:
    items: list[dict] = []
    counts: set[int] = set()
    for page in pages:
        payload = json.loads((cache / f"{prefix}-page{page}.json").read_text())
        counts.add(int(payload["count"]))
        items.extend(payload["items"])
    if len(counts) != 1:
        raise ValueError(f"Pagination count drift for {prefix}: {counts}")
    return items


def nonempty(value: object) -> bool:
    if value is None or pd.isna(value):
        return False
    return bool(str(value).strip())


def infer_nucleic_acid(protocol: str, title: str) -> str:
    combined = f"{protocol} {title}".upper()
    if "NUC-DNA" in combined:
        return "DNA"
    if "NUC-RNA" in combined:
        return "RNA"
    return "Unresolved"


def infer_size_fraction(protocol: str) -> str:
    match = re.search(r"_W([^_]+)$", protocol)
    return match.group(1) if match else "Unresolved"


def main() -> None:
    args = parse_args()
    cache = args.cache_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_dir = output / "source"
    source_dir.mkdir(exist_ok=True)

    source_manifest = json.loads((cache / "download-manifest.json").read_text())
    if source_manifest.get("article") != ARTICLE:
        raise ValueError("Unexpected Article 73 source manifest")
    if source_manifest.get("snapshot_date") != SNAPSHOT_DATE:
        raise ValueError("Unexpected public-resource snapshot date")
    for name, record in source_manifest["resources"].items():
        path = cache / name
        if path.stat().st_size != record["bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"Source checksum mismatch: {name}")
        shutil.copy2(path, source_dir / name)
    shutil.copy2(cache / "download-manifest.json", output / "source-manifest.json")
    shutil.copy2(cache / "gem-figure1-original.png", output / "gem-figure1-original.png")

    study = json.loads((cache / "mgnify-study.json").read_text())
    samples_raw = load_pages(cache, "mgnify-samples", (1, 2))
    analyses_raw = load_pages(cache, "mgnify-analyses", (1, 2, 3))
    if len(samples_raw) != 136 or len(analyses_raw) != 249:
        raise ValueError("Unexpected MGnify relationship counts")

    sample_rows = []
    for sample in samples_raw:
        metadata = sample.get("metadata", {})
        protocol = str(metadata.get("protocol_label", ""))
        title = str(sample.get("sample_title", ""))
        depth_annotation = metadata.get("depth", "")
        sample_rows.append(
            {
                "SampleAccession": sample["accession"],
                "ENAAccessions": ";".join(sample.get("ena_accessions", [])),
                "SampleTitle": title,
                "Latitude": pd.to_numeric(metadata.get("lat", ""), errors="coerce"),
                "Longitude": pd.to_numeric(metadata.get("lon", ""), errors="coerce"),
                "DepthAnnotation": depth_annotation,
                "DepthM": pd.to_numeric(depth_annotation, errors="coerce"),
                "SamplingSite": metadata.get("sampling_site", ""),
                "ProtocolLabel": protocol,
                "NucleicAcid": infer_nucleic_acid(protocol, title),
                "SizeFractionMicrometre": infer_size_fraction(protocol),
                "CollectionDateStart": metadata.get("collection_date_start", ""),
                "CollectionDateEnd": metadata.get("collection_date_end", ""),
                "EnvironmentFeature": metadata.get("environment_feature", ""),
                "EnvironmentMaterial": metadata.get("environment_material", ""),
                "BroadContext": metadata.get("broad_scale_environmental_context", ""),
                "LocalContext": metadata.get("local_environmental_context", ""),
                "Temperature": metadata.get("temperature", ""),
                "Salinity": metadata.get("salinity", ""),
                "Status": metadata.get("status", ""),
                "MetadataStudyAccessions": metadata.get("study_accession", ""),
            }
        )
    samples = pd.DataFrame(sample_rows).sort_values("SampleAccession", kind="stable")
    if not samples["SampleAccession"].is_unique:
        raise ValueError("Duplicate sample accessions")
    samples["ValidCoordinates"] = (
        samples["Latitude"].between(-90, 90)
        & samples["Longitude"].between(-180, 180)
    )
    write_tsv(samples, output / "tara-samples.tsv")

    completeness_fields = [
        ("Coordinates", ["Latitude", "Longitude"], "Geographic coverage and spatial strata"),
        ("Numeric depth", ["DepthM"], "Scalar depth for modeling and comparability"),
        ("Sampling site", ["SamplingSite"], "Station-level blocking and repeated visits"),
        ("Protocol label", ["ProtocolLabel"], "DNA/RNA and size-fraction eligibility"),
        ("Collection date", ["CollectionDateStart"], "Seasonality and temporal ordering"),
        ("Environment feature", ["EnvironmentFeature"], "ENVO-informed habitat layer"),
        ("Environment material", ["EnvironmentMaterial"], "What physical material was sampled"),
        ("Temperature", ["Temperature"], "Environmental covariate"),
        ("Salinity", ["Salinity"], "Environmental covariate"),
    ]
    completeness_rows = []
    for field, columns, use in completeness_fields:
        if field == "Coordinates":
            available = samples["ValidCoordinates"]
        else:
            available = samples[columns].apply(
                lambda row: all(nonempty(value) for value in row), axis=1
            )
        count = int(available.sum())
        completeness_rows.append(
            {
                "Field": field,
                "Available": count,
                "Total": len(samples),
                "CompletenessPct": 100 * count / len(samples),
                "AnalysisUse": use,
            }
        )
    metadata_completeness = pd.DataFrame(completeness_rows)
    write_tsv(metadata_completeness, output / "tara-metadata-completeness.tsv")

    analysis_rows = []
    for analysis in analyses_raw:
        run = analysis.get("run") or {}
        sample = analysis.get("sample") or {}
        title = str(sample.get("sample_title", ""))
        analysis_rows.append(
            {
                "AnalysisAccession": analysis["accession"],
                "PipelineVersion": analysis.get("pipeline_version", ""),
                "ExperimentType": analysis.get("experiment_type", ""),
                "StudyAccession": analysis.get("study_accession", ""),
                "RunAccession": run.get("accession", ""),
                "SampleAccession": sample.get("accession", ""),
                "SampleTitle": title,
                "NucleicAcid": infer_nucleic_acid("", title),
                "AssemblyAccession": (analysis.get("assembly") or {}).get("accession", ""),
            }
        )
    analyses = pd.DataFrame(analysis_rows).sort_values(
        ["PipelineVersion", "AnalysisAccession"], kind="stable"
    )
    if not analyses["AnalysisAccession"].is_unique or not analyses["RunAccession"].is_unique:
        raise ValueError("Expected one analysis per run in the frozen response")
    write_tsv(analyses, output / "tara-analyses.tsv")

    analysis_lineage = (
        analyses.groupby("PipelineVersion", sort=True)
        .agg(
            Analyses=("AnalysisAccession", "size"),
            Runs=("RunAccession", "nunique"),
            Samples=("SampleAccession", "nunique"),
            DNAAnalyses=("NucleicAcid", lambda values: int((values == "DNA").sum())),
            RNAAnalyses=("NucleicAcid", lambda values: int((values == "RNA").sum())),
        )
        .reset_index()
    )
    write_tsv(analysis_lineage, output / "tara-analysis-lineage.tsv")

    downloads = pd.DataFrame(
        [
            {
                "Alias": item["alias"],
                "DownloadType": item["download_type"],
                "DownloadGroup": item["download_group"],
                "ShortDescription": item["short_description"],
                "URL": item["url"],
            }
            for item in study["downloads"]
        ]
    )
    write_tsv(downloads, output / "tara-study-downloads.tsv")

    catalogues_payload = json.loads((cache / "mgnify-catalogues.json").read_text())
    catalogue_rows = []
    for item in catalogues_payload["items"]:
        biome = item.get("biome") or {}
        lineage = str(biome.get("lineage", ""))
        if lineage.startswith("root:Environmental"):
            scope = "Environmental"
        elif ":Plants:" in lineage:
            scope = "Plant-associated"
        else:
            scope = "Host-associated"
        catalogue_rows.append(
            {
                "CatalogueID": item["catalogue_id"],
                "Name": item["name"],
                "Version": item["version"],
                "Status": item["status"],
                "CatalogueType": item["catalogue_type"],
                "Scope": scope,
                "Biome": item.get("catalogue_biome_label", ""),
                "BiomeLineage": lineage,
                "RepresentativeClusters": item.get("genome_count", 0),
                "InputGenomes": item.get("unclustered_genome_count", 0),
                "PipelineVersion": item.get("pipeline_version_tag", ""),
                "UpdatedAt": item.get("updated_at", ""),
                "FTPURL": item.get("ftp_url", ""),
            }
        )
    catalogues = pd.DataFrame(catalogue_rows).sort_values(
        "RepresentativeClusters", ascending=False, kind="stable"
    )
    write_tsv(catalogues, output / "mgnify-catalogues.tsv")

    gem = pd.read_csv(cache / "gem-genome-metadata.tsv", sep="\t")
    if len(gem) != 52_515 or gem["genome_id"].duplicated().any():
        raise ValueError("Unexpected GEM genome table")
    identity_error = (
        gem["quality_score"]
        - (gem["completeness"] - 5 * gem["contamination"])
    ).abs()
    if identity_error.max() > 1e-9:
        raise ValueError("GEM quality-score identity failed")
    if not (
        gem["completeness"].ge(50).all()
        and gem["contamination"].le(5).all()
        and gem["quality_score"].ge(50).all()
    ):
        raise ValueError("GEM minimum quality contract failed")

    gem_biome = (
        gem.groupby("ecosystem_category", sort=False)
        .agg(
            MAGs=("genome_id", "size"),
            Metagenomes=("metagenome_id", "nunique"),
            SpeciesOTUs=("otu_id", "nunique"),
            HighQuality=("mimag_quality", lambda values: int((values == "HQ").sum())),
            MeanCompleteness=("completeness", "mean"),
            MeanContamination=("contamination", "mean"),
        )
        .reset_index()
        .rename(columns={"ecosystem_category": "EcosystemCategory"})
        .sort_values("MAGs", ascending=False, kind="stable")
    )
    gem_biome["HighQualityPct"] = 100 * gem_biome["HighQuality"] / gem_biome["MAGs"]
    write_tsv(gem_biome, output / "gem-biome-summary.tsv")

    rng = np.random.default_rng(ANALYSIS_SEED)
    quality_parts = []
    for category in gem_biome.head(8)["EcosystemCategory"]:
        subset = gem.loc[gem["ecosystem_category"].eq(category)]
        choose = min(len(subset), 700)
        quality_parts.append(subset.iloc[rng.choice(len(subset), choose, replace=False)])
    gem_quality = pd.concat(quality_parts, ignore_index=True)[
        [
            "genome_id",
            "completeness",
            "contamination",
            "quality_score",
            "mimag_quality",
            "ecosystem_category",
        ]
    ].rename(
        columns={
            "genome_id": "GenomeID",
            "completeness": "Completeness",
            "contamination": "Contamination",
            "quality_score": "QualityScore",
            "mimag_quality": "MIMAGQuality",
            "ecosystem_category": "EcosystemCategory",
        }
    )
    write_tsv(gem_quality, output / "gem-quality-visualization-sample.tsv")

    valid_geo = gem.loc[
        gem["longitude"].between(-180, 180)
        & gem["latitude"].between(-90, 90)
    ].copy()
    map_n = min(12_000, len(valid_geo))
    gem_map = valid_geo.iloc[rng.choice(len(valid_geo), map_n, replace=False)][
        ["genome_id", "ecosystem_category", "longitude", "latitude", "mimag_quality"]
    ].rename(
        columns={
            "genome_id": "GenomeID",
            "ecosystem_category": "EcosystemCategory",
            "longitude": "Longitude",
            "latitude": "Latitude",
            "mimag_quality": "MIMAGQuality",
        }
    )
    write_tsv(gem_map, output / "gem-map-visualization-sample.tsv")
    write_tsv(gem.head(12), output / "gem-metadata-preview.tsv")

    gtdb_history = pd.DataFrame(
        [
            {
                "Release": "R10-RS226",
                "ReleaseDate": "2025-04-16",
                "Genomes": 732_475,
                "SpeciesClusters": 143_614,
                "VersionFile": "v226",
            },
            {
                "Release": "R11-RS232",
                "ReleaseDate": "2026-04-15",
                "Genomes": 901_341,
                "SpeciesClusters": 199_923,
                "VersionFile": "v232",
            },
        ]
    )
    write_tsv(gtdb_history, output / "gtdb-release-history.tsv")

    md5_index = {}
    for line in (cache / "gtdb-r232-md5.txt").read_text().splitlines():
        if not line.strip():
            continue
        digest, relative = line.split(maxsplit=1)
        md5_index[relative.removeprefix("./")] = digest
    base = "https://data.gtdb.ecogenomic.org/releases/release232/232.0/"
    selected = [
        ("Archaea metadata", "ar53_metadata_r232.tsv.gz", 7_015_869),
        ("Archaea taxonomy", "ar53_taxonomy_r232.tsv.gz", 311_583),
        ("Bacteria metadata", "bac120_metadata_r232.tsv.gz", 289_043_801),
        ("Bacteria taxonomy", "bac120_taxonomy_r232.tsv.gz", 9_886_191),
        (
            "GTDB-Tk full reference package",
            "auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz",
            60_806_405_195,
        ),
    ]
    gtdb_files = pd.DataFrame(
        [
            {
                "Purpose": purpose,
                "RelativePath": relative,
                "Bytes": size,
                "MD5": md5_index[relative],
                "URL": base + relative,
            }
            for purpose, relative, size in selected
        ]
    )
    write_tsv(gtdb_files, output / "gtdb-selected-files.tsv")

    resource_registry = pd.DataFrame(
        [
            (
                "Archive identity",
                "ENA",
                "PRJEB1787 / ERP001736",
                "Study, sample, experiment, run, raw reads",
                "Primary sequence archive and accession lineage",
                "An ENA project is not a harmonized analysis cohort",
                "https://www.ebi.ac.uk/ena/browser/view/PRJEB1787",
            ),
            (
                "Standardized analysis",
                "MGnify API v2",
                "MGYS00000410",
                "Study, sample, run, MGYA analysis, download",
                "Searchable metadata and versioned pipeline output",
                "API version, pipeline version and database release are different fields",
                "https://www.ebi.ac.uk/metagenomics/api/v2/studies/MGYS00000410",
            ),
            (
                "Cross-biome MAG snapshot",
                "GEM",
                "DOI 10.1038/s41587-020-0718-6",
                "52,515 MAGs from a published fixed snapshot",
                "Global comparative genomics and historical benchmark",
                "The paper-era taxonomy is not current GTDB and the catalogue is not exhaustive",
                "https://portal.nersc.gov/GEM/",
            ),
            (
                "Biome-specific genome catalogue",
                "MGnify Genomes Marine v2.0",
                "marine-v2-0",
                "50,866 input genomes; 13,223 representative species clusters",
                "Biome-matched novelty search, genomes and protein catalogues",
                "Representative-cluster counts across catalogues are not globally dereplicated",
                "https://www.ebi.ac.uk/metagenomics/api/v2/genomes/catalogues/marine-v2-0",
            ),
            (
                "Taxonomy backbone",
                "GTDB",
                "R11-RS232",
                "901,341 genomes; 199,923 species clusters",
                "Release-consistent genome taxonomy and GTDB-Tk reference",
                "GTDB is not an abundance table, ecological catalogue or functional annotation",
                "https://data.gtdb.ecogenomic.org/releases/release232/232.0/",
            ),
        ],
        columns=[
            "Layer",
            "Resource",
            "StableIdentity",
            "Unit",
            "Use",
            "Boundary",
            "URL",
        ],
    )
    write_tsv(resource_registry, output / "resource-registry.tsv")

    accession_crosswalk = pd.DataFrame(
        [
            ("MGnify study", "MGYS00000410", 1, "Discovery and grouped downloads"),
            ("ENA project IDs", "PRJEB1787; ERP001736", 2, "Archive identity"),
            ("Related BioSamples", "SAMEA / ERS", samples["SampleAccession"].nunique(), "Specimen metadata"),
            ("Sequencing runs", "ERR", analyses["RunAccession"].nunique(), "Technical sequencing unit"),
            ("MGnify analyses", "MGYA", analyses["AnalysisAccession"].nunique(), "Pipeline result unit"),
            ("DNA-labelled samples", "protocol NUC-DNA", int(samples["NucleicAcid"].eq("DNA").sum()), "Shotgun-DNA eligibility"),
            ("RNA-labelled samples", "protocol NUC-RNA", int(samples["NucleicAcid"].eq("RNA").sum()), "Exclude or analyze separately"),
            ("Study summary downloads", "ERP001736 aliases", len(downloads), "Legacy V2/V3 summary outputs"),
        ],
        columns=["Object", "IdentifierPattern", "Count", "Role"],
    )
    write_tsv(accession_crosswalk, output / "accession-crosswalk.tsv")

    decisions = pd.DataFrame(
        [
            ("Find environmental studies", "MGnify v2 study/sample JSON", "Kilobytes to megabytes", "Freeze accessions and snapshot date"),
            ("Reanalyse raw reads", "ENA run table + selected FASTQ", "Study-dependent", "Keep BioSample and run crosswalk"),
            ("Reuse standardized output", "Exact MGYA download file", "Megabytes", "Record pipeline and reference-database versions"),
            ("Screen MAG novelty in one biome", "MGnify catalogue representatives", "Gigabytes", "Do not sum species across catalogues"),
            ("Summarize GEM geography/quality", "GEM genome_metadata.tsv", "12.35 MB", "No FASTA is needed for metadata questions"),
            ("Retaxonomize your MAGs", "GTDB-Tk R232 reference package", "60.81 GB compressed", "Verify official MD5 and match GTDB-Tk"),
            ("Join existing GTDB labels", "R232 taxonomy TSVs", "0.31 + 9.89 MB compressed", "Join by genome accession and retain release"),
        ],
        columns=["Question", "MinimumDownload", "Scale", "AuditRequirement"],
    )
    write_tsv(decisions, output / "download-decisions.tsv")

    metrics = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "mgnify_api": "v2",
        "mgnify_study": "MGYS00000410",
        "mgnify_samples": len(samples),
        "mgnify_dna_samples": int(samples["NucleicAcid"].eq("DNA").sum()),
        "mgnify_rna_samples": int(samples["NucleicAcid"].eq("RNA").sum()),
        "mgnify_runs": int(analyses["RunAccession"].nunique()),
        "mgnify_analyses": len(analyses),
        "mgnify_v2_analyses": int(analyses["PipelineVersion"].eq("V2").sum()),
        "mgnify_v3_analyses": int(analyses["PipelineVersion"].eq("V3").sum()),
        "mgnify_catalogues": len(catalogues),
        "marine_input_genomes": 50_866,
        "marine_representative_clusters": 13_223,
        "gem_mags": len(gem),
        "gem_source_metagenomes_with_mags": int(gem["metagenome_id"].nunique()),
        "gem_species_otus": int(gem["otu_id"].nunique()),
        "gem_high_quality": int(gem["mimag_quality"].eq("HQ").sum()),
        "gem_medium_quality": int(gem["mimag_quality"].eq("MQ").sum()),
        "gem_georeferenced": len(valid_geo),
        "gem_mean_completeness": float(gem["completeness"].mean()),
        "gem_mean_contamination": float(gem["contamination"].mean()),
        "gtdb_release": "R11-RS232",
        "gtdb_genomes": 901_341,
        "gtdb_species_clusters": 199_923,
        "gtdbtk_reference_bytes": 60_806_405_195,
    }
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    contract = {
        "article": ARTICLE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "snapshot_date": SNAPSHOT_DATE,
        "unit_contract": {
            "archive": "study -> sample -> run",
            "analysis": "one MGYA per frozen run response",
            "gem": "one published MAG row",
            "mgnify_catalogue": "one catalogue-specific representative cluster",
            "gtdb": "one release-specific species cluster",
        },
        "eligibility_rule": "Use protocol NUC-DNA for the shotgun-DNA subset; do not infer eligibility from the study title alone.",
        "catalogue_limit": "Catalogue representative counts are not abundance and are not globally dereplicated across resources.",
        "taxonomy_limit": "Taxonomic names are release-bound; retain original and retaxonomized labels side by side.",
        "visualization_sampling": {
            "gem_quality": "up to 700 MAGs from each of the eight largest ecosystem categories",
            "gem_map": "12,000 checksum-frozen georeferenced MAG rows",
        },
    }
    (output / "methods-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
