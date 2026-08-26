#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999, digits = 17)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260728L
set.seed(primary_seed)

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
required <- c("source-dir", "output-dir", "notice")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c("curatedMetagenomicData", "digest")
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

source_dir <- normalizePath(args[["source-dir"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
notice_path <- normalizePath(args[["notice"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(notice_path), recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x, con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}

cohort_contract <- data.frame(
  Study = c(
    "ZellerG_2014", "FengQ_2015", "YuJ_2015", "VogtmannE_2016",
    "WirbelJ_2018", "ThomasAM_2018a", "ThomasAM_2018b",
    "ThomasAM_2019_c"
  ),
  Cohort = c("FR", "AT", "CN", "US", "DE", "IT-A", "IT-B", "JP"),
  Role = c(
    rep("Discovery / meta-analysis", 5L),
    rep("Independent validation", 3L)
  ),
  ExpectedControl = c(61L, 61L, 54L, 52L, 65L, 24L, 28L, 40L),
  ExpectedCRC = c(53L, 46L, 74L, 52L, 60L, 29L, 32L, 40L),
  ExpectedProfileColumns = c(156L, 154L, 128L, 110L, 125L, 80L, 60L, 80L),
  ExpectedSpeciesRows = c(661L, 615L, 584L, 548L, 547L, 486L, 512L, 526L),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expected_source_hashes <- c(
  ZellerG_2014 = "d8e0f3fd00b2339b1aa929197ca0869c43990ff885a04fc675e70d4aff5604b2",
  FengQ_2015 = "7d33b4813f36ee3aa49e5ade0598bb6283c724a44591f316d49d9c349c871304",
  YuJ_2015 = "7c5b4fc3bac6e83371105e147b40c41aa4399f005619da4ed2a92fa6ef777272",
  VogtmannE_2016 = "ba624fa7a51c002909b54c8783073e861f3c512c4295cc82571bc6cce9dfe8da",
  WirbelJ_2018 = "8b7cdc0f4e6ba7546970bef6ae8b62d798488aecb9d0b6af627fb0bd95117fdd",
  ThomasAM_2018a = "f353bf89c9fa749b9308f27bb3aaedad17a1d48615670b7564c689ac888daf48",
  ThomasAM_2018b = "6bb1dd1669c824a25e33de719612e54e597cf83a811b8f2173cb099e778f07f8",
  ThomasAM_2019_c = "1c32d598e1d050cfa4fba8f13d2a644d8e8146b2aba1dde45ebb30d9e015883e"
)

profile_path <- function(study) {
  file.path(
    source_dir,
    paste0("2021-03-31.", study, ".relative_abundance.rda")
  )
}
profile_paths <- stats::setNames(
  vapply(cohort_contract$Study, profile_path, character(1L)),
  cohort_contract$Study
)
missing_profiles <- names(profile_paths)[!file.exists(profile_paths)]
if (length(missing_profiles) > 0L) {
  stop("Missing source profiles: ", paste(missing_profiles, collapse = ", "), call. = FALSE)
}
observed_source_hashes <- vapply(profile_paths, sha256_file, character(1L))
if (!identical(observed_source_hashes[names(expected_source_hashes)], expected_source_hashes)) {
  bad <- names(expected_source_hashes)[
    observed_source_hashes[names(expected_source_hashes)] != expected_source_hashes
  ]
  stop("Source checksum mismatch: ", paste(bad, collapse = ", "), call. = FALSE)
}

metadata_all <- curatedMetagenomicData::sampleMetadata
stopifnot(
  is.data.frame(metadata_all),
  nrow(metadata_all) == 22588L,
  ncol(metadata_all) == 141L,
  identical(as.character(utils::packageVersion("curatedMetagenomicData")), "3.12.0")
)

load_profile <- function(path) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  stopifnot(length(loaded) == 1L)
  object <- env[[loaded[[1L]]]]
  stopifnot(is.matrix(object), is.numeric(object))
  object
}

profile_list <- list()
metadata_list <- list()
species_union <- character()
for (i in seq_len(nrow(cohort_contract))) {
  study <- cohort_contract$Study[[i]]
  profile <- load_profile(profile_paths[[study]])
  species_row <- grepl("\\|s__", rownames(profile)) &
    !grepl("\\|t__", rownames(profile))
  species_profile <- profile[species_row, , drop = FALSE]

  metadata <- metadata_all[
    metadata_all$study_name == study &
      metadata_all$body_site == "stool" &
      metadata_all$study_condition %in% c("control", "CRC"),
    , drop = FALSE
  ]
  metadata <- metadata[match(intersect(colnames(profile), metadata$sample_id), metadata$sample_id), , drop = FALSE]
  stopifnot(
    ncol(profile) == cohort_contract$ExpectedProfileColumns[[i]],
    nrow(species_profile) == cohort_contract$ExpectedSpeciesRows[[i]],
    sum(metadata$study_condition == "control") == cohort_contract$ExpectedControl[[i]],
    sum(metadata$study_condition == "CRC") == cohort_contract$ExpectedCRC[[i]],
    nrow(metadata) == cohort_contract$ExpectedControl[[i]] + cohort_contract$ExpectedCRC[[i]],
    !anyDuplicated(metadata$sample_id),
    !anyDuplicated(metadata$subject_id),
    all(metadata$sample_id %in% colnames(species_profile))
  )

  species_profile <- species_profile[, metadata$sample_id, drop = FALSE] / 100
  stopifnot(
    all(is.finite(species_profile)),
    min(species_profile) >= 0,
    max(species_profile) <= 1,
    max(abs(colSums(species_profile) - 1)) < 2e-6
  )
  profile_list[[study]] <- species_profile
  metadata_list[[study]] <- metadata
  species_union <- union(species_union, rownames(species_profile))
}
stopifnot(length(species_union) == 897L)

extract_species_label <- function(feature) {
  pieces <- strsplit(feature, "\\|", fixed = FALSE)[[1L]]
  species <- sub("^s__", "", pieces[grepl("^s__", pieces)][[1L]])
  if (!nzchar(species) || species == "unclassified") {
    genus_pieces <- pieces[grepl("^g__", pieces)]
    genus <- if (length(genus_pieces) > 0L) {
      sub("^g__", "", genus_pieces[[length(genus_pieces)]])
    } else {
      "Taxon"
    }
    species <- paste(genus, "unclassified")
  }
  gsub("_", " ", species, fixed = TRUE)
}
species_labels <- vapply(species_union, extract_species_label, character(1L))
stopifnot(!anyDuplicated(species_union), !anyDuplicated(species_labels))

metadata <- do.call(rbind, lapply(seq_len(nrow(cohort_contract)), function(i) {
  study <- cohort_contract$Study[[i]]
  x <- metadata_list[[study]]
  x$Cohort <- cohort_contract$Cohort[[i]]
  x$StudyRole <- cohort_contract$Role[[i]]
  x$Outcome <- ifelse(x$study_condition == "CRC", "CRC", "Control")
  x$StudySampleKey <- paste(x$study_name, x$sample_id, sep = "::")
  x
}))
rownames(metadata) <- NULL
metadata$Cohort <- factor(metadata$Cohort, levels = cohort_contract$Cohort)
metadata$Outcome <- factor(metadata$Outcome, levels = c("Control", "CRC"))
stopifnot(
  nrow(metadata) == 771L,
  !anyDuplicated(metadata$StudySampleKey),
  !anyDuplicated(metadata$subject_id),
  identical(as.integer(table(metadata$Outcome)), c(385L, 386L))
)

sample_ids <- metadata$sample_id
combined <- matrix(
  0, nrow = length(species_union), ncol = length(sample_ids),
  dimnames = list(species_union, sample_ids)
)
for (study in cohort_contract$Study) {
  block <- profile_list[[study]]
  combined[rownames(block), colnames(block)] <- block
}
stopifnot(
  all(is.finite(combined)),
  min(combined) >= 0,
  max(combined) <= 1,
  max(abs(colSums(combined) - 1)) < 2e-6
)

species <- data.frame(
  Feature = species_union,
  Species = species_labels,
  combined,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cohort_summary <- do.call(rbind, lapply(seq_len(nrow(cohort_contract)), function(i) {
  study <- cohort_contract$Study[[i]]
  rows <- metadata$study_name == study
  age <- suppressWarnings(as.numeric(metadata$age[rows]))
  data.frame(
    Cohort = cohort_contract$Cohort[[i]],
    Study = study,
    Role = cohort_contract$Role[[i]],
    Country = paste(sort(unique(metadata$country[rows])), collapse = ";"),
    Samples = sum(rows),
    Controls = sum(metadata$Outcome[rows] == "Control"),
    CRC = sum(metadata$Outcome[rows] == "CRC"),
    IndependentSubjects = length(unique(metadata$subject_id[rows])),
    AgeAvailable = sum(is.finite(age)),
    MedianAge = if (any(is.finite(age))) stats::median(age, na.rm = TRUE) else NA_real_,
    GenderAvailable = sum(!is.na(metadata$gender[rows]) & nzchar(metadata$gender[rows])),
    SpeciesAtSource = nrow(profile_list[[study]]),
    SourceProfileSHA256 = observed_source_hashes[[study]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))

cohort_prevalence <- vapply(cohort_contract$Study, function(study) {
  ids <- metadata$sample_id[metadata$study_name == study]
  rowMeans(combined[, ids, drop = FALSE] > 0)
}, numeric(nrow(combined)))
colnames(cohort_prevalence) <- cohort_contract$Cohort
feature_audit <- data.frame(
  Feature = species_union,
  Species = species_labels,
  Samples = ncol(combined),
  NonzeroSamples = rowSums(combined > 0),
  NonzeroPrevalence = rowMeans(combined > 0),
  CohortsWithAnyDetection = rowSums(cohort_prevalence > 0),
  CohortsAtOrAbove5PctPrevalence = rowSums(cohort_prevalence >= 0.05),
  MaximumRelativeAbundance = apply(combined, 1L, max),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
feature_audit <- cbind(
  feature_audit,
  stats::setNames(as.data.frame(cohort_prevalence), paste0("Prevalence_", cohort_contract$Cohort))
)
feature_audit <- feature_audit[order(
  -feature_audit$CohortsAtOrAbove5PctPrevalence,
  -feature_audit$NonzeroPrevalence,
  -feature_audit$MaximumRelativeAbundance,
  feature_audit$Species
), , drop = FALSE]

analysis_contract <- data.frame(
  Item = c(
    "seed", "resource_release", "cohorts", "samples", "control_samples",
    "crc_samples", "raw_species_features", "outcome", "positive_class",
    "outer_validation", "outer_test_unit", "inner_validation",
    "inner_validation_unit", "model", "lambda_grid", "lambda_selection",
    "cohort_class_weights", "training_prevalence_filter",
    "training_abundance_filter", "training_cohort_support_filter",
    "log10_pseudocount", "standardization", "primary_auc_interval",
    "performance_meta_analysis", "hierarchical_bootstrap",
    "single_study_transfer", "batch_correction_policy"
  ),
  Value = c(
    as.character(primary_seed), "2021-03-31", "8", "771", "385", "386",
    "897", "CRC vs Control; adenoma and post-surgery profiles excluded", "CRC",
    "8 outer leave-one-dataset-out iterations; N-1 cohorts train and one complete cohort tests",
    "study_name", "leave-one-training-cohort-out within each outer training set",
    "study_name", "glmnet L1-penalized binomial logistic regression",
    "41 fixed values from 1 to 0.0001", "one-standard-error rule on mean inner-cohort AUROC",
    "each training cohort contributes equal total weight; CRC and Control each contribute half within cohort",
    "nonzero prevalence >= 0.10 in the current training pool",
    "maximum relative abundance >= 0.0001 in the current training pool",
    "prevalence >= 0.05 in at least half of current training cohorts",
    "0.000001", "glmnet training-only feature centering and scaling",
    "DeLong 95% CI within each untouched outer cohort",
    "random-effects REML on logit AUROC; report tau-squared and I-squared",
    "2000 cohort-and-class-stratified bootstrap replicates for macro AUROC",
    "secondary author-faithful 8-by-8 single-study train/test matrix",
    "no joint train-test ComBat; any transductive batch normalization must be labeled sensitivity-only"
  ),
  Interpretation = c(
    "Fixed before splitting, tuning, fitting and bootstrap",
    "All profiles use one cMD release and one MetaPhlAn3 lineage",
    "Five discovery/meta-analysis and three independent-validation cohorts",
    "One stool profile per independent subject",
    "Class count before any model fitting",
    "Class count before any model fitting",
    "Union of species rows across eight profiles; absent rows are structural zero under the shared database",
    "Non-CRC disease states are not relabeled as controls",
    "All probabilities and AUROC values target CRC",
    "The outer cohort never participates in filtering, transformation, tuning or calibration",
    "The generalization unit is an entire study, not a random subject fold",
    "Hyperparameters are selected for transfer to an unseen training cohort",
    "No sample-level random split substitutes for study-level selection",
    "Sparse coefficients permit signature-stability auditing",
    "Grid is prespecified rather than inferred from outer-test performance",
    "Choose the largest eligible lambda to favor a sparse transferable signature",
    "Large cohorts and majority classes cannot dominate pooled training",
    "Fit on training data only and applied unchanged to validation/test data",
    "Fit on training data only and applied unchanged to validation/test data",
    "Requires a taxon to recur across training cohorts without consulting labels",
    "Fixed transformation parameter",
    "Held-out cohorts use training means and scales through the fitted glmnet object",
    "Intervals are cohort-specific external-validation uncertainty",
    "Heterogeneity is treated as evidence, not averaged away",
    "Resamples cohorts and subjects rather than treating cohorts as fixed duplicates",
    "Reproduces the original transfer question but is not the primary N-1 model",
    "Combining the outer test cohort in batch correction would leak its distribution"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

resource_manifest <- do.call(rbind, lapply(seq_len(nrow(cohort_contract)), function(i) {
  study <- cohort_contract$Study[[i]]
  filename <- basename(profile_paths[[study]])
  data.frame(
    Resource = paste0("cMD species profile: ", study),
    Release = paste0("2021-03-31.", study, ".relative_abundance"),
    Source = paste0(
      "https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/",
      "curatedMetagenomicData/2021-03-31/", study, "/", filename
    ),
    SHA256 = observed_source_hashes[[study]],
    Bytes = file.info(profile_paths[[study]])$size,
    UnitAtSource = "percent at each taxonomic rank",
    FrozenUnit = "species relative-abundance fraction",
    InterpretationBoundary = "Uniform cMD reprocessing; not the historical author profile",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
resource_manifest <- rbind(
  resource_manifest,
  data.frame(
    Resource = c("cMD sample metadata", "Anchor publication"),
    Release = c("curatedMetagenomicData 3.12.0", "Nature Medicine 25:679-689 (2019)"),
    Source = c(
      "curatedMetagenomicData::sampleMetadata",
      "https://doi.org/10.1038/s41591-019-0406-6"
    ),
    SHA256 = c(NA_character_, NA_character_),
    Bytes = c(NA_real_, NA_real_),
    UnitAtSource = c("one row per metagenomic profile", NA_character_),
    FrozenUnit = c("one CRC/control row per independent subject", NA_character_),
    InterpretationBoundary = c(
      "Metadata completeness and condition labels differ by cohort",
      "Original analysis used historical profiles and SIAMCAT; numeric equality is not claimed"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)

metadata_columns <- c(
  "sample_id", "subject_id", "StudySampleKey", "study_name", "Cohort",
  "StudyRole", "Outcome", "study_condition", "body_site", "country", "age",
  "age_category", "gender", "BMI", "sequencing_platform", "DNA_extraction_kit",
  "PMID", "number_reads", "number_bases", "NCBI_accession"
)
metadata_columns <- metadata_columns[metadata_columns %in% names(metadata)]
metadata_out_table <- metadata[, metadata_columns, drop = FALSE]
metadata_out_table$Cohort <- as.character(metadata_out_table$Cohort)
metadata_out_table$Outcome <- as.character(metadata_out_table$Outcome)

species_out <- file.path(output_dir, "species-relative-abundance.tsv.gz")
metadata_out <- file.path(output_dir, "sample-metadata.tsv")
cohort_out <- file.path(output_dir, "cohort-summary.tsv")
feature_out <- file.path(output_dir, "feature-universe-audit.tsv.gz")
contract_out <- file.path(output_dir, "analysis-contract.tsv")
manifest_out <- file.path(output_dir, "resource-manifest.tsv")
checksum_out <- file.path(output_dir, "file-checksums.sha256")

write_tsv_gz(species, species_out)
write_tsv(metadata_out_table, metadata_out)
write_tsv(cohort_summary, cohort_out)
write_tsv_gz(feature_audit, feature_out)
write_tsv(analysis_contract, contract_out)
write_tsv(resource_manifest, manifest_out)

writeLines(
  c(
    "Article 28 frozen data notice",
    "",
    "Source: eight stool shotgun-metagenomic CRC cohorts distributed by curatedMetagenomicData 3.12.0.",
    "Resource release: 2021-03-31; profiler lineage: MetaPhlAn 3 / CHOCOPhlAn 201901.",
    "Included cohorts: FR, AT, CN, US, DE, IT-A, IT-B and JP.",
    "The binary task contains 771 independent subjects: 385 Control and 386 CRC.",
    "Adenoma, carcinoma-surgery-history and unlabeled profiles are excluded rather than relabeled.",
    "The union table contains 897 species and relative-abundance fractions summing to one per sample.",
    "Every outer iteration trains on seven complete cohorts and tests once on the eighth complete cohort.",
    "Feature filtering, standardization and lambda selection are repeated using training cohorts only.",
    "The original paper used historical profiles and SIAMCAT; exact numeric reproduction is not claimed.",
    "Raw ExperimentHub RDA files remain outside Git; checksum-locked compact tables are stored in data/small."
  ),
  con = notice_path,
  useBytes = TRUE
)

payloads <- c(
  species_out, metadata_out, cohort_out,
  feature_out, contract_out, manifest_out
)
writeLines(
  paste(vapply(payloads, sha256_file, character(1L)), basename(payloads)),
  checksum_out,
  useBytes = TRUE
)

cat("Prepared Article 28 frozen inputs\n")
cat("  cohorts:", nrow(cohort_summary), "\n")
cat("  samples:", nrow(metadata_out_table), "\n")
cat("  class counts: Control 385 / CRC 386\n")
cat("  union species:", nrow(species), "\n")
cat("  checksum payloads:", length(payloads), "\n")
