#!/usr/bin/env python3
"""Prepare the checksum-gated real cohort and AGORA inputs for Article 61."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from cobra.io import load_json_model, write_sbml_model
from micom.qiime_formats import load_qiime_medium, metadata


SEED = 61001
MIN_READS = 1_000_000
RELATIVE_ABUNDANCE_CUTOFF = 0.001  # proportion; 0.1%
ARTICLE22_SHA256 = {
    "human-sample-metadata.tsv": "2bac387f2e9dbd6786a79de98e2f2c032b4d249360daacb7b65f2ff53ef32f2a",
    "species-relative-abundance.tsv.gz": "312e55c506279ef5e5caa0e20f9bb5f14be06bf5b22feabf682e9314f29951dc",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def digest(path: Path, algorithm: str = "sha256") -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def read_resource_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 3 or len({row["Asset"] for row in rows}) != 3:
        raise RuntimeError("Article 61 resource manifest must contain three unique assets")
    return {row["Asset"]: row for row in rows}


def verify_resource(path: Path, row: dict[str, str]) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    observed = {
        "Asset": row["Asset"],
        "LocalFile": path.name,
        "ObservedBytes": path.stat().st_size,
        "ObservedMD5": digest(path, "md5"),
        "ObservedSHA256": digest(path),
    }
    if observed["ObservedBytes"] != int(row["ExpectedBytes"]):
        raise RuntimeError(f"Byte-count gate failed for {path}")
    if observed["ObservedSHA256"] != row["SHA256"]:
        raise RuntimeError(f"SHA-256 gate failed for {path}")
    if row["PublisherMD5"] != "NA" and observed["ObservedMD5"] != row["PublisherMD5"]:
        raise RuntimeError(f"Publisher MD5 gate failed for {path}")
    return observed


def verify_article22(source: Path) -> list[dict[str, object]]:
    audit: list[dict[str, object]] = []
    for name, expected in ARTICLE22_SHA256.items():
        path = source / name
        observed = digest(path)
        if observed != expected:
            raise RuntimeError(f"Article 22 checksum gate failed for {path}")
        audit.append(
            {
                "Asset": f"Article22:{name}",
                "LocalFile": str(path),
                "ObservedBytes": path.stat().st_size,
                "ObservedMD5": digest(path, "md5"),
                "ObservedSHA256": observed,
            }
        )
    return audit


def qza_data_member(qza: Path, basename: str) -> str:
    with zipfile.ZipFile(qza) as archive:
        matches = [name for name in archive.namelist() if name.endswith(f"/data/{basename}")]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {basename} in {qza}, observed {matches}")
    return matches[0]


def extract_database(qza: Path, destination: Path) -> tuple[Path, pd.DataFrame]:
    meta = metadata(str(qza))
    if meta.get("type") != "MetabolicModels[JSON]":
        raise RuntimeError(f"Wrong QIIME artifact type: {meta}")
    uuid = meta["uuid"]
    database = destination / uuid / "data"
    stamp = destination / ".identity.json"
    identity = {"uuid": uuid, "sha256": digest(qza), "bytes": qza.stat().st_size}
    if stamp.is_file():
        if json.loads(stamp.read_text(encoding="utf-8")) != identity:
            raise RuntimeError("Existing extracted database has a different archive identity")
    else:
        if destination.exists() and any(destination.iterdir()):
            raise RuntimeError(f"Refusing to overlay a non-empty extraction directory: {destination}")
        destination.mkdir(parents=True, exist_ok=True)
        root = destination.resolve()
        with zipfile.ZipFile(qza) as archive:
            for member in archive.infolist():
                target = (destination / member.filename).resolve()
                if root != target and root not in target.parents:
                    raise RuntimeError(f"Unsafe QIIME archive member: {member.filename}")
            archive.extractall(destination)
        stamp.write_text(json.dumps(identity, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_path = database / "manifest.csv"
    manifest = pd.read_csv(manifest_path)
    required = {"id", "species", "file", "summary_rank"}
    if not required.issubset(manifest.columns):
        raise RuntimeError(f"AGORA manifest lacks columns: {required - set(manifest.columns)}")
    if set(manifest.summary_rank) != {"species"} or manifest.species.duplicated().any():
        raise RuntimeError("AGORA artifact is not a unique species-level database")
    return database, manifest


def slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    if not slug:
        raise RuntimeError(f"Could not create model identifier from {value!r}")
    return slug


def sample_selection(metadata_frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    adult = metadata_frame.loc[metadata_frame.age_category.eq("adult")].copy()
    adult["SourceOrder"] = np.arange(len(adult))
    adult["PassReadDepth"] = adult.number_reads.ge(MIN_READS)
    eligible = adult.loc[adult.PassReadDepth].sort_values(["subject_id", "SourceOrder"])
    selected_ids = set(eligible.drop_duplicates("subject_id", keep="first").sample_id)
    adult["Selected"] = adult.sample_id.isin(selected_ids)

    first_eligible_order = eligible.groupby("subject_id").SourceOrder.min().to_dict()

    def reason(row: pd.Series) -> str:
        if row.Selected:
            return "Selected: earliest visit passing read-depth gate"
        if not row.PassReadDepth:
            return "Excluded: fewer than 1,000,000 reads"
        if row.SourceOrder > first_eligible_order.get(row.subject_id, -1):
            return "Excluded: later eligible visit from selected subject"
        return "Excluded: not selected"

    adult["Decision"] = adult.apply(reason, axis=1)
    selected = adult.loc[adult.Selected].sort_values("SourceOrder").copy()
    if selected.subject_id.duplicated().any() or len(selected) != 6:
        raise RuntimeError(f"Expected six independent adult subjects, observed {len(selected)}")
    if not selected.number_reads.ge(MIN_READS).all():
        raise RuntimeError("A selected sample failed the read-depth gate")
    return adult, selected


def prepare_taxonomy(
    profiles: pd.DataFrame,
    selected: pd.DataFrame,
    database_manifest: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    model_species = set(database_manifest.species.astype(str))
    records: list[dict[str, object]] = []
    taxonomy: list[dict[str, object]] = []
    coverage: list[dict[str, object]] = []
    for sample in selected.itertuples(index=False):
        column = sample.sample_id
        frame = profiles.loc[profiles[column].gt(0), ["Feature", "Species", column]].copy()
        frame.rename(columns={column: "Percent"}, inplace=True)
        total = float(frame.Percent.sum())
        if total <= 0 or total > 100.1:
            raise RuntimeError(f"Invalid species-resolved abundance for {column}: {total}")
        # Keep the original whole-profile denominator. MICOM normalizes the
        # modeled members internally, while this audit retains unresolved
        # higher-rank abundance as an explicit coverage loss.
        frame["Abundance"] = frame.Percent / 100.0
        frame["AboveCutoff"] = frame.Abundance.ge(RELATIVE_ABUNDANCE_CUTOFF)
        frame["InDatabase"] = frame.Species.isin(model_species)
        frame["IncludedInModel"] = frame.AboveCutoff & frame.InDatabase
        for row in frame.itertuples(index=False):
            records.append(
                {
                    "SampleID": column,
                    "SubjectID": sample.subject_id,
                    "Feature": row.Feature,
                    "Species": row.Species,
                    "Percent": row.Percent,
                    "Abundance": row.Abundance,
                    "AboveCutoff": row.AboveCutoff,
                    "InDatabase": row.InDatabase,
                    "IncludedInModel": row.IncludedInModel,
                }
            )
            taxonomy.append(
                {
                    "sample_id": column,
                    "id": slugify(row.Species),
                    "species": row.Species,
                    "abundance": row.Abundance,
                }
            )
        coverage.append(
            {
                "SampleID": column,
                "SubjectID": sample.subject_id,
                "Reads": int(sample.number_reads),
                "ObservedTaxa": len(frame),
                "TaxaAboveCutoff": int(frame.AboveCutoff.sum()),
                "MatchedTaxaAboveCutoff": int(frame.IncludedInModel.sum()),
                "SpeciesResolvedAbundance": total / 100.0,
                "AbundanceRetained": float(frame.loc[frame.AboveCutoff, "Abundance"].sum()),
                "ModeledAbundance": float(frame.loc[frame.IncludedInModel, "Abundance"].sum()),
                "AllMatchedAbundance": float(frame.loc[frame.InDatabase, "Abundance"].sum()),
            }
        )
    taxonomy_frame = pd.DataFrame(taxonomy)
    audit_frame = pd.DataFrame(records)
    coverage_frame = pd.DataFrame(coverage)
    if coverage_frame.ModeledAbundance.min() < 0.5:
        raise RuntimeError("A selected sample has less than 50% modeled abundance")

    equal_rows: list[dict[str, object]] = []
    for sample_id, group in audit_frame.loc[audit_frame.IncludedInModel].groupby("SampleID", sort=False):
        n_taxa = len(group)
        for species in group.Species:
            equal_rows.append(
                {
                    "sample_id": sample_id,
                    "id": slugify(species),
                    "species": species,
                    "abundance": 1.0 / n_taxa,
                }
            )
    return taxonomy_frame, pd.DataFrame(equal_rows), audit_frame, coverage_frame


def prepare_smetana_models(
    work: Path,
    database: Path,
    database_manifest: pd.DataFrame,
    audit: pd.DataFrame,
    focal_sample: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    focal = audit.loc[
        audit.SampleID.eq(focal_sample) & audit.IncludedInModel
    ].sort_values(["Abundance", "Species"], ascending=[False, True]).head(6)
    if len(focal) != 6:
        raise RuntimeError("The focal sample has fewer than six matched abundant taxa")
    merged = focal.merge(
        database_manifest[["id", "species", "file"]],
        left_on="Species",
        right_on="species",
        validate="one_to_one",
    )
    model_dir = work / "smetana/models"
    model_dir.mkdir(parents=True, exist_ok=True)
    model_audit: list[dict[str, object]] = []
    subcommunity: list[dict[str, object]] = []
    for rank, row in enumerate(merged.itertuples(index=False), start=1):
        slug = slugify(row.Species)
        json_path = database / row.file
        sbml_path = model_dir / f"{slug}.xml"
        model = load_json_model(str(json_path))
        model.id = slug
        # AGORA JSON uses metabolite IDs such as ``ac[e]`` and exchange IDs
        # such as ``EX_ac(e)``. COBRA writes the brackets as SBML escape
        # sequences, whereas SMETANA 1.2.1 constructs pooled exchanges from
        # the conventional ``ac_e`` form. Normalize only terminal compartment
        # tokens before export so the declared medium actually reaches the
        # community models.
        metabolite_ids = {
            metabolite.id: re.sub(r"\[([A-Za-z0-9]+)\]$", r"_\1", metabolite.id)
            for metabolite in model.metabolites
        }
        reaction_ids = {
            reaction.id: re.sub(r"\(([A-Za-z0-9]+)\)$", r"_\1", reaction.id)
            for reaction in model.reactions
        }
        if len(set(metabolite_ids.values())) != len(metabolite_ids):
            raise RuntimeError(f"SMETANA metabolite normalization collides for {row.Species}")
        if len(set(reaction_ids.values())) != len(reaction_ids):
            raise RuntimeError(f"SMETANA reaction normalization collides for {row.Species}")
        for metabolite in model.metabolites:
            metabolite.id = metabolite_ids[metabolite.id]
        for reaction in model.reactions:
            reaction.id = reaction_ids[reaction.id]
        model.repair()
        external_metabolites = [
            metabolite.id for metabolite in model.metabolites
            if metabolite.compartment == "e"
        ]
        if not external_metabolites or not all(
            metabolite.endswith("_e") for metabolite in external_metabolites
        ):
            raise RuntimeError(f"SMETANA extracellular ID normalization failed for {row.Species}")
        write_sbml_model(model, str(sbml_path))
        model_audit.append(
            {
                "Species": row.Species,
                "ModelID": slug,
                "SourceJSON": row.file,
                "SBML": sbml_path.name,
                "Reactions": len(model.reactions),
                "Metabolites": len(model.metabolites),
                "Genes": len(model.genes),
                "Exchanges": len(model.exchanges),
                "SMETANAIDConvention": "metabolite_compartment",
                "SBMLBytes": sbml_path.stat().st_size,
                "SBMLSHA256": digest(sbml_path),
            }
        )
        subcommunity.append(
            {
                "SampleID": focal_sample,
                "Rank": rank,
                "Species": row.Species,
                "ModelID": slug,
                "ObservedAbundance": row.Abundance,
                "SBML": str(sbml_path),
            }
        )
    if len({row["ModelID"] for row in subcommunity}) != 6:
        raise RuntimeError("SMETANA model identifiers are not unique")
    return pd.DataFrame(subcommunity), pd.DataFrame(model_audit)


def main() -> None:
    args = parse_args()
    random.seed(SEED)
    np.random.seed(SEED)
    root = args.project_root.resolve()
    cache = args.cache_dir.resolve()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)

    resources = read_resource_manifest(root / "data/small/61-community-database-manifest.tsv")
    resource_paths = {
        "AGORA201RefSeq216Species": cache / "agora201_refseq216_species_1.qza",
        "AGORA2RefSeqSpeciesManifest": cache / "agora2_refseq_species.tsv",
        "WesternDietGutAGORA": cache / "western_diet_gut_agora.qza",
    }
    resource_audit = [verify_resource(resource_paths[key], resources[key]) for key in resource_paths]
    article22 = root / "data/small/22-diversity-inputs"
    resource_audit.extend(verify_article22(article22))

    database, database_manifest = extract_database(
        resource_paths["AGORA201RefSeq216Species"], work / "database/agora201"
    )
    external_manifest = pd.read_csv(
        resource_paths["AGORA2RefSeqSpeciesManifest"], sep="\t", low_memory=False
    )
    external_species = set(external_manifest.species.dropna().astype(str))
    internal_species = set(database_manifest.species.astype(str))
    if not internal_species.issubset(external_species):
        raise RuntimeError("QIIME database species are not a subset of the locked source manifest")

    sample_metadata = pd.read_csv(article22 / "human-sample-metadata.tsv", sep="\t")
    profiles = pd.read_csv(article22 / "species-relative-abundance.tsv.gz", sep="\t")
    selection_audit, selected = sample_selection(sample_metadata)
    taxonomy, equal_taxonomy, taxon_audit, coverage = prepare_taxonomy(
        profiles, selected, database_manifest
    )
    focal_sample = selected.sort_values("SourceOrder").sample_id.iloc[0]

    medium = load_qiime_medium(str(resource_paths["WesternDietGutAGORA"])).copy()
    if medium.reaction.duplicated().any() or (medium.flux < 0).any():
        raise RuntimeError("Growth medium contains duplicated reactions or negative bounds")
    medium["PositiveFlux"] = medium.flux.gt(0)
    medium["Compound"] = medium.metabolite.str.replace(r"_m$", "", regex=True)
    medium["SMETANACompound"] = medium.Compound.str.replace("__", "_", regex=False)
    medium["SMETANAMembership"] = medium.PositiveFlux & medium.Compound.ne("o2")
    smetana_medium = (
        medium.loc[medium.SMETANAMembership, ["SMETANACompound"]]
        .rename(columns={"SMETANACompound": "compound"})
        .assign(medium="western_diet_gut_anoxic_membership")[["medium", "compound"]]
        .drop_duplicates()
        .sort_values("compound")
    )
    if smetana_medium.compound.eq("o2").any() or not medium.Compound.eq("o2").any():
        raise RuntimeError("The declared MICOM-to-SMETANA oxygen adaptation failed")

    subcommunity, model_audit = prepare_smetana_models(
        work, database, database_manifest, taxon_audit, focal_sample
    )

    selected_out = selected[
        [
            "study_name", "sample_id", "subject_id", "age_category", "gender",
            "NCBI_accession", "PMID", "number_reads", "SourceOrder",
        ]
    ].rename(
        columns={
            "sample_id": "SampleID", "subject_id": "SubjectID",
            "number_reads": "Reads",
        }
    )
    write_tsv(work / "resource-audit.tsv", pd.DataFrame(resource_audit))
    write_tsv(work / "sample-selection-audit.tsv", selection_audit)
    write_tsv(work / "selected-samples.tsv", selected_out)
    write_tsv(work / "micom-taxonomy.tsv", taxonomy)
    write_tsv(work / "micom-taxonomy-equal.tsv", equal_taxonomy)
    write_tsv(work / "taxon-match-audit.tsv", taxon_audit)
    write_tsv(work / "model-coverage.tsv", coverage)
    write_tsv(work / "database-manifest.tsv", database_manifest)
    write_tsv(work / "medium.tsv", medium.reset_index(drop=True))
    write_tsv(work / "smetana-media.tsv", smetana_medium)
    write_tsv(work / "smetana-subcommunity.tsv", subcommunity)
    write_tsv(work / "smetana-model-audit.tsv", model_audit)
    (work / "database-path.txt").write_text(str(database) + "\n", encoding="utf-8")
    contract = {
        "article": 61,
        "seed": SEED,
        "minimum_reads": MIN_READS,
        "relative_abundance_cutoff": RELATIVE_ABUNDANCE_CUTOFF,
        "selected_samples": selected_out.SampleID.tolist(),
        "selected_subjects": selected_out.SubjectID.tolist(),
        "focal_smetana_sample": focal_sample,
        "smetana_subcommunity_size": len(subcommunity),
        "agora_species_models": len(database_manifest),
        "agora_source_species": len(external_species),
        "minimum_modeled_abundance": float(coverage.ModeledAbundance.min()),
        "medium_reactions": len(medium),
        "positive_medium_reactions": int(medium.PositiveFlux.sum()),
        "micom_oxygen_flux": float(medium.loc[medium.Compound.eq("o2"), "flux"].iloc[0]),
        "smetana_medium_components": len(smetana_medium),
        "smetana_medium_adaptation": "same positive-flux compound set except oxygen; membership only",
    }
    (work / "preparation-contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
