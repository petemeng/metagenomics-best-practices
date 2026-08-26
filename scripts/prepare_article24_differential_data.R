#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
Sys.setenv(TZ = "UTC")

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
    }
    out[[substring(args[[i]], 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("species-rda", "pathway-rda", "output-dir")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required.", call. = FALSE)
}

species_path <- normalizePath(args[["species-rda"]], mustWork = TRUE)
pathway_path <- normalizePath(args[["pathway-rda"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

md5_file <- function(path) {
  digest::digest(file = path, algo = "md5", serialize = FALSE)
}

load_one <- function(path) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (length(loaded) != 1L) {
    stop("Expected exactly one object in ", path, call. = FALSE)
  }
  list(name = loaded[[1L]], value = get(loaded[[1L]], envir = env, inherits = FALSE))
}

write_wide_gz <- function(x, path, annotations) {
  stopifnot(nrow(x) == nrow(annotations), ncol(x) > 0L)
  con <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(con), add = TRUE)
  out <- data.frame(annotations, x, check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.table(
    out,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  invisible(path)
}

expected_sources <- data.frame(
  resource = c("relative_abundance", "pathway_abundance"),
  object = c(
    "2021-03-31.ZellerG_2014.relative_abundance",
    "2021-03-31.ZellerG_2014.pathway_abundance"
  ),
  bytes = c(147044, 4503559),
  md5 = c(
    "a03686e5adc99d03ca1f026591bfcc84",
    "43e0f9599d8ec62aec6ffa0350351e7b"
  ),
  sha256 = c(
    "d8e0f3fd00b2339b1aa929197ca0869c43990ff885a04fc675e70d4aff5604b2",
    "bf333a90b36cf875960d584500be6ce1dabb0240981f797493fd98ff1dc60bf6"
  ),
  url = c(
    paste0(
      "https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/",
      "curatedMetagenomicData/2021-03-31/ZellerG_2014/",
      "2021-03-31.ZellerG_2014.relative_abundance.rda"
    ),
    paste0(
      "https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/",
      "curatedMetagenomicData/2021-03-31/ZellerG_2014/",
      "2021-03-31.ZellerG_2014.pathway_abundance.rda"
    )
  ),
  stringsAsFactors = FALSE
)
source_paths <- c(species_path, pathway_path)
observed_sources <- expected_sources
observed_sources$observed_bytes <- as.numeric(file.info(source_paths)$size)
observed_sources$observed_md5 <- vapply(source_paths, md5_file, character(1L))
observed_sources$observed_sha256 <- vapply(source_paths, sha256_file, character(1L))
stopifnot(
  identical(observed_sources$observed_bytes, as.numeric(expected_sources$bytes)),
  identical(observed_sources$observed_md5, expected_sources$md5),
  identical(observed_sources$observed_sha256, expected_sources$sha256)
)

species_object <- load_one(species_path)
pathway_object <- load_one(pathway_path)
stopifnot(
  identical(species_object$name, expected_sources$object[[1L]]),
  identical(pathway_object$name, expected_sources$object[[2L]])
)
species_all <- species_object$value
pathway_all <- pathway_object$value
stopifnot(
  is.matrix(species_all),
  is.matrix(pathway_all),
  identical(dim(species_all), c(1019L, 156L)),
  identical(dim(pathway_all), c(22620L, 156L)),
  identical(colnames(species_all), colnames(pathway_all)),
  !anyNA(species_all),
  !anyNA(pathway_all),
  min(species_all) >= 0,
  min(pathway_all) >= 0,
  max(abs(colSums(species_all) - 700)) < 1e-3
)

species_keep <- grepl("|s__", rownames(species_all), fixed = TRUE)
species_percent <- species_all[species_keep, , drop = FALSE]
species_label <- sub("^.*\\|s__", "", rownames(species_percent))
species_label <- gsub("_", " ", species_label, fixed = TRUE)
stopifnot(
  nrow(species_percent) == 661L,
  max(abs(colSums(species_percent) - 100)) < 1e-3,
  !anyDuplicated(rownames(species_percent)),
  !anyDuplicated(species_label)
)

pathway_keep <- !grepl("|", rownames(pathway_all), fixed = TRUE) &
  !rownames(pathway_all) %in% c("UNMAPPED", "UNGROUPED", "UNINTEGRATED")
pathway_ordinary <- pathway_all[pathway_keep, , drop = FALSE]
pathway_denominator <- colSums(pathway_ordinary)
pathway_fraction <- sweep(pathway_ordinary, 2L, pathway_denominator, "/")
stopifnot(
  nrow(pathway_fraction) == 493L,
  all(pathway_denominator > 0),
  max(abs(colSums(pathway_fraction) - 1)) < 1e-10,
  !anyNA(pathway_fraction),
  !anyDuplicated(rownames(pathway_fraction))
)

package_path <- find.package("curatedMetagenomicData", quiet = TRUE)
if (!nzchar(package_path)) {
  stop("curatedMetagenomicData must be installed to recover sampleMetadata.", call. = FALSE)
}
metadata_env <- new.env(parent = emptyenv())
invisible(lazyLoad(file.path(package_path, "data", "Rdata"), envir = metadata_env))
if (!exists("sampleMetadata", envir = metadata_env, inherits = FALSE)) {
  stop("sampleMetadata was not found in curatedMetagenomicData.", call. = FALSE)
}
metadata_all <- get("sampleMetadata", envir = metadata_env, inherits = FALSE)
metadata_fields <- c(
  "study_name", "sample_id", "subject_id", "body_site", "study_condition",
  "disease", "age", "age_category", "gender", "country", "BMI",
  "number_reads", "NCBI_accession", "PMID", "disease_stage", "fobt"
)
sample_metadata <- metadata_all[
  match(colnames(species_percent), metadata_all$sample_id),
  metadata_fields,
  drop = FALSE
]
stopifnot(
  identical(sample_metadata$sample_id, colnames(species_percent)),
  all(sample_metadata$study_name == "ZellerG_2014"),
  all(sample_metadata$body_site == "stool"),
  length(unique(sample_metadata$subject_id)) == 156L,
  !anyDuplicated(sample_metadata$sample_id),
  setequal(unique(sample_metadata$study_condition), c("control", "adenoma", "CRC")),
  identical(
    as.integer(table(factor(
      sample_metadata$study_condition,
      levels = c("control", "adenoma", "CRC")
    ))),
    c(61L, 42L, 53L)
  )
)
sample_metadata$analysis_group <- factor(
  sample_metadata$study_condition,
  levels = c("control", "adenoma", "CRC"),
  labels = c("Control", "Adenoma", "CRC")
)
sample_metadata$primary_complete_case <-
  sample_metadata$study_condition %in% c("control", "CRC") &
  stats::complete.cases(sample_metadata[, c("age", "gender", "BMI", "number_reads")])
stopifnot(
  sum(sample_metadata$primary_complete_case) == 110L,
  identical(
    as.integer(table(droplevels(sample_metadata$analysis_group[sample_metadata$primary_complete_case]))),
    c(59L, 51L)
  )
)

species_fraction <- species_percent / 100
species_pseudocount <- round(sweep(
  species_fraction,
  2L,
  sample_metadata$number_reads,
  "*"
))
storage.mode(species_pseudocount) <- "double"
pseudocount_rounding_error <- colSums(species_pseudocount) - sample_metadata$number_reads
stopifnot(
  min(species_pseudocount) >= 0,
  max(abs(pseudocount_rounding_error)) <= nrow(species_pseudocount),
  !anyNA(species_pseudocount)
)

species_tsv <- file.path(output_dir, "species-relative-abundance.tsv.gz")
pseudocount_tsv <- file.path(output_dir, "species-pseudocounts.tsv.gz")
pathway_tsv <- file.path(output_dir, "pathway-relative-abundance.tsv.gz")
metadata_tsv <- file.path(output_dir, "sample-metadata.tsv")
contract_tsv <- file.path(output_dir, "analysis-contract.tsv")
manifest_tsv <- file.path(output_dir, "resource-manifest.tsv")
checksum_tsv <- file.path(output_dir, "file-checksums.sha256")
notice_txt <- file.path(output_dir, "data-NOTICE.txt")

write_wide_gz(
  species_fraction,
  species_tsv,
  data.frame(
    Feature = rownames(species_fraction),
    Species = species_label,
    stringsAsFactors = FALSE
  )
)
write_wide_gz(
  species_pseudocount,
  pseudocount_tsv,
  data.frame(
    Feature = rownames(species_pseudocount),
    Species = species_label,
    stringsAsFactors = FALSE
  )
)
write_wide_gz(
  pathway_fraction,
  pathway_tsv,
  data.frame(Pathway = rownames(pathway_fraction), stringsAsFactors = FALSE)
)
sample_metadata$analysis_group <- as.character(sample_metadata$analysis_group)
utils::write.table(
  sample_metadata,
  metadata_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

analysis_contract <- data.frame(
  item = c(
    "seed", "primary_groups", "primary_complete_cases", "primary_covariates",
    "species_min_abundance", "species_min_prevalence",
    "pathway_denominator", "pathway_min_abundance", "pathway_min_prevalence",
    "maaslin3_abundance_null", "maaslin3_prevalence_null",
    "ancombc2_input", "aldex2_input", "aldex2_mc_samples", "fdr_method"
  ),
  value = c(
    "20260724", "Control_vs_CRC", "110",
    "z_age+gender+z_BMI+z_log10_reads", "0.0001", "0.10",
    "ordinary_unstratified_pathways_closed_to_one", "0.00001", "0.20",
    "per-metadatum coefficient median", "zero log-odds coefficient",
    "round(species_relative_fraction*whole_metagenome_number_reads)",
    "same_reconstructed_species_pseudocount", "128", "BH_within_method_feature_space_component"
  ),
  interpretation = c(
    "Fixed before model fitting",
    "Adenoma excluded from the binary primary contrast",
    "59 Control and 51 CRC independent subjects",
    "Read depth is included because detection affects prevalence",
    "Feature must exceed 0.01 percent relative abundance",
    "Threshold applied in the 110-sample primary set",
    "Conditional composition, not total functional explanation",
    "Feature must exceed 0.001 percent of the ordinary-pathway denominator",
    "Threshold applied in the 110-sample primary set",
    "Compositionality-aware primary abundance test",
    "Presence model is not median centered",
    "Sensitivity input; not observed taxon reads",
    "Sensitivity input; not observed taxon reads",
    "Deterministic Monte Carlo design with fixed seed",
    "No shared FDR family across taxa, pathways, or model components"
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  analysis_contract,
  contract_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

resource_manifest <- data.frame(
  Resource = c("Species source", "Pathway source", "Species derived", "Pseudo-count derived", "Pathway derived"),
  Release = c(
    "curatedMetagenomicData 2021-03-31",
    "curatedMetagenomicData 2021-03-31",
    "MetaPhlAn3 species-level fraction",
    "cMD counts convenience transformation",
    "HUMAnN3 ordinary unstratified pathway fraction"
  ),
  URL = c(expected_sources$url, rep("Derived locally from checksum-locked source", 3L)),
  SourceBytes = c(expected_sources$bytes, rep(NA_real_, 3L)),
  SourceMD5 = c(expected_sources$md5, rep(NA_character_, 3L)),
  SourceSHA256 = c(expected_sources$sha256, rep(NA_character_, 3L)),
  OutputFile = c(
    basename(species_path), basename(pathway_path), basename(species_tsv),
    basename(pseudocount_tsv), basename(pathway_tsv)
  ),
  Rows = c(1019L, 22620L, nrow(species_fraction), nrow(species_pseudocount), nrow(pathway_fraction)),
  Samples = rep(ncol(species_fraction), 5L),
  Unit = c(
    "Percent across seven taxonomic ranks",
    "HUMAnN community relative abundance plus stratified rows",
    "Fraction within species-level profile",
    "Rounded reconstructed pseudo-count",
    "Fraction within annotated ordinary pathways"
  ),
  InterpretationBoundary = c(
    "Do not analyze all taxonomic ranks together",
    "Do not combine stratified and community rows",
    "Relative composition; no absolute microbial load",
    "Not observed taxon reads and not absolute abundance",
    "Conditional on annotated ordinary pathways"
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  resource_manifest,
  manifest_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

writeLines(
  c(
    "Article 24 frozen data notice",
    "",
    "Source: curatedMetagenomicData 3.12.0, ZellerG_2014 resources dated 2021-03-31.",
    "Original study: Zeller et al., Molecular Systems Biology 2014, DOI 10.15252/msb.20145645.",
    "The source payloads remain outside Git; this directory contains compact derived tables and checksums.",
    "Species values are MetaPhlAn relative fractions after selecting only species-level rows.",
    "Pathway values are HUMAnN ordinary unstratified pathways reclosed within that annotated denominator.",
    "species-pseudocounts.tsv.gz is round(relative fraction * whole-metagenome read depth).",
    "The pseudo-count table is a method-sensitivity input, not observed taxon reads or absolute abundance.",
    "All 156 ZellerG_2014 profiles are retained; the binary primary model uses 110 complete-case Control/CRC subjects."
  ),
  con = notice_txt,
  useBytes = TRUE
)

derived_files <- c(
  species_tsv, pseudocount_tsv, pathway_tsv, metadata_tsv,
  contract_tsv, manifest_tsv, notice_txt
)
checksum_lines <- paste(
  vapply(derived_files, sha256_file, character(1L)),
  basename(derived_files)
)
writeLines(checksum_lines, con = checksum_tsv, useBytes = TRUE)

cat("Prepared Article 24 frozen inputs\n")
cat("  samples:", ncol(species_fraction), "\n")
cat("  independent subjects:", length(unique(sample_metadata$subject_id)), "\n")
cat("  primary complete cases:", sum(sample_metadata$primary_complete_case), "\n")
cat("  species rows:", nrow(species_fraction), "\n")
cat("  ordinary pathway rows:", nrow(pathway_fraction), "\n")
cat("  pseudo-count rounding error range:", paste(range(pseudocount_rounding_error), collapse = " to "), "\n")
