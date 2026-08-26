#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999, digits = 17)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260727L
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
required <- c("source-species", "source-metadata", "output-dir", "notice")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required.", call. = FALSE)
}

species_path <- normalizePath(args[["source-species"]], mustWork = TRUE)
metadata_path <- normalizePath(args[["source-metadata"]], mustWork = TRUE)
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

expected_source_hashes <- c(
  species = "b16e688a255d2b9c9aeb8114c7b56b90c9e1ab5d070957c6edc57e7cddc86844",
  metadata = "e748bcf0ed3806d27cc42837ffab2f26fb0e0a09a950183f6d0f2a6241ecc309"
)
observed_source_hashes <- c(
  species = sha256_file(species_path),
  metadata = sha256_file(metadata_path)
)
if (!identical(observed_source_hashes, expected_source_hashes)) {
  bad <- names(expected_source_hashes)[observed_source_hashes != expected_source_hashes]
  stop("Source checksum mismatch: ", paste(bad, collapse = ", "), call. = FALSE)
}

species_all <- utils::read.delim(
  gzfile(species_path), check.names = FALSE,
  quote = "", comment.char = "", stringsAsFactors = FALSE
)
metadata_all <- utils::read.delim(
  metadata_path, check.names = FALSE,
  quote = "", comment.char = "", stringsAsFactors = FALSE
)
stopifnot(
  nrow(species_all) == 661L,
  ncol(species_all) == 158L,
  identical(names(species_all)[seq_len(2L)], c("Feature", "Species")),
  nrow(metadata_all) == 156L,
  identical(names(species_all)[-(1:2)], metadata_all$sample_id),
  !anyDuplicated(metadata_all$sample_id),
  !anyDuplicated(metadata_all$subject_id),
  all(metadata_all$study_name == "ZellerG_2014"),
  all(metadata_all$body_site == "stool")
)

keep <- metadata_all$study_condition %in% c("control", "CRC")
metadata <- metadata_all[keep, , drop = FALSE]
metadata$Outcome <- factor(
  metadata$study_condition,
  levels = c("control", "CRC"),
  labels = c("Control", "CRC")
)
metadata$StudySampleKey <- paste(metadata$study_name, metadata$sample_id, sep = "::")
stopifnot(
  nrow(metadata) == 114L,
  length(unique(metadata$subject_id)) == 114L,
  !anyDuplicated(metadata$StudySampleKey),
  identical(as.integer(table(metadata$Outcome)), c(61L, 53L))
)

sample_columns <- metadata$sample_id
species <- species_all[, c("Feature", "Species", sample_columns), drop = FALSE]
abundance <- as.matrix(species[, sample_columns, drop = FALSE])
storage.mode(abundance) <- "double"
stopifnot(
  all(is.finite(abundance)),
  min(abundance) >= 0,
  max(abundance) <= 1,
  max(abs(colSums(abundance) - 1)) < 2e-6,
  !anyDuplicated(species$Feature),
  !anyDuplicated(species$Species)
)

feature_audit <- data.frame(
  Feature = species$Feature,
  Species = species$Species,
  Samples = ncol(abundance),
  NonzeroSamples = rowSums(abundance > 0),
  NonzeroPrevalence = rowMeans(abundance > 0),
  SamplesAtOrAbove0.01Pct = rowSums(abundance >= 1e-4),
  PrevalenceAtOrAbove0.01Pct = rowMeans(abundance >= 1e-4),
  MaximumRelativeAbundance = apply(abundance, 1L, max),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
feature_audit <- feature_audit[order(
  -feature_audit$NonzeroPrevalence,
  -feature_audit$MaximumRelativeAbundance,
  feature_audit$Species
), , drop = FALSE]

analysis_contract <- data.frame(
  Item = c(
    "seed", "study", "outcome", "positive_class", "samples",
    "control_samples", "crc_samples", "raw_species_features",
    "outer_resamples", "outer_folds", "inner_folds",
    "fold_unit", "training_prevalence_filter", "training_abundance_filter",
    "log10_pseudocount", "algorithms", "tuning_metric",
    "sample_score", "auc_interval", "paired_difference_bootstrap",
    "feature_importance", "leakage_audit_permutations"
  ),
  Value = c(
    as.character(primary_seed), "ZellerG_2014", "CRC vs Control", "CRC", "114",
    "61", "53", "661", "5", "5", "4", "subject_id",
    "nonzero prevalence >= 0.10 within each training fold",
    "maximum relative abundance >= 0.0001 within each training fold",
    "0.000001", "ranger probability forest | xgboost binary logistic",
    "inner-fold AUROC", "mean of five out-of-fold probabilities per subject",
    "DeLong 95% CI on subject-level aggregated OOF scores", "2000 subject bootstraps",
    "outer-test permutation delta AUROC", "50"
  ),
  Interpretation = c(
    "Fixed before every fold, fit, bootstrap and permutation",
    "One study and one processing lineage; external validation is Article 28",
    "Adenoma profiles are excluded rather than relabeled",
    "All reported probabilities target CRC",
    "Every included sample is an independent subject",
    "Class count before any cross-validation split",
    "Class count before any cross-validation split",
    "No global label-informed feature selection",
    "Each subject receives one held-out prediction per repeat",
    "Stratified outer error-estimation folds",
    "Stratified inner hyperparameter-selection folds",
    "Related samples must never cross folds",
    "Fit on training data only and applied unchanged to held-out data",
    "Fit on training data only and applied unchanged to held-out data",
    "Log transform parameter is fixed a priori",
    "Both models are prespecified and both are reported",
    "Outer test folds never select hyperparameters",
    "Reduces repeated-CV predictions to one score per independent subject",
    "Conditional discrimination interval; not external-validation uncertainty",
    "Paired model comparison resamples subjects, not individual predictions",
    "Importance is measured on held-out folds, not by impurity decrease",
    "Safe and deliberately leaky branches use permuted subject labels"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

resource_manifest <- data.frame(
  Resource = c(
    "Species table", "Sample metadata", "Underlying cMD resource",
    "Original publication", "Original Figure 1"
  ),
  Release = c(
    "Article 24 checksum-locked frozen table",
    "Article 24 checksum-locked frozen table",
    "2021-03-31.ZellerG_2014.relative_abundance",
    "Molecular Systems Biology 10:766 (2014)",
    "CC BY 4.0"
  ),
  Source = c(
    "data/small/24-differential-abundance/species-relative-abundance.tsv.gz",
    "data/small/24-differential-abundance/sample-metadata.tsv",
    "curatedMetagenomicData 3.12.0 / MetaPhlAn 3 / CHOCOPhlAn 201901",
    "https://doi.org/10.15252/msb.20145645",
    "figures/27-zeller-fig1-original.png"
  ),
  SHA256 = c(
    observed_source_hashes[["species"]],
    observed_source_hashes[["metadata"]],
    "d8e0f3fd00b2339b1aa929197ca0869c43990ff885a04fc675e70d4aff5604b2",
    NA_character_,
    "6f0dbe5ca4ad7e9bc853fd6568efca093e30f8f46f02a0f090c282b16059ac43"
  ),
  Unit = c(
    "species relative-abundance fraction", "sample/subject metadata",
    "percent at source", NA_character_, NA_character_
  ),
  InterpretationBoundary = c(
    "Relative composition; no absolute microbial load",
    "Only CRC and control are included in the binary task",
    "Uniform cMD reprocessing is not the historical author profile",
    "Original classifier used nested LASSO; current audit compares tree ensembles",
    "Unmodified anchor figure; model scores are not recreated pixel by pixel"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

species_out <- file.path(output_dir, "species-relative-abundance.tsv.gz")
metadata_out <- file.path(output_dir, "sample-metadata.tsv")
feature_out <- file.path(output_dir, "feature-universe-audit.tsv")
contract_out <- file.path(output_dir, "analysis-contract.tsv")
manifest_out <- file.path(output_dir, "resource-manifest.tsv")
checksum_out <- file.path(output_dir, "file-checksums.sha256")

write_tsv_gz(species, species_out)
metadata$Outcome <- as.character(metadata$Outcome)
write_tsv(metadata, metadata_out)
write_tsv(feature_audit, feature_out)
write_tsv(analysis_contract, contract_out)
write_tsv(resource_manifest, manifest_out)

writeLines(
  c(
    "Article 27 frozen data notice",
    "",
    "Source study: ZellerG_2014, Molecular Systems Biology 2014, DOI 10.15252/msb.20145645.",
    "Profile source: curatedMetagenomicData 3.12.0 resource dated 2021-03-31.",
    "Profiler lineage: MetaPhlAn 3 with CHOCOPhlAn 201901 as recorded by cMD3.",
    "The table contains all 661 species rows and 114 independent CRC/control subjects.",
    "Adenoma profiles are excluded from this binary task; they are not relabeled as controls.",
    "Values are species relative-abundance fractions and do not measure absolute microbial load.",
    "Feature filtering, log transformation and hyperparameter tuning are repeated inside training folds.",
    "The original study used nested LASSO; this chapter audits random forest and gradient boosting.",
    "Article 28, not this chapter, performs multi-cohort external validation."
  ),
  con = notice_path,
  useBytes = TRUE
)

payloads <- c(species_out, metadata_out, feature_out, contract_out, manifest_out)
checksum_lines <- paste(
  vapply(payloads, sha256_file, character(1L)),
  basename(payloads)
)
writeLines(checksum_lines, checksum_out, useBytes = TRUE)

cat("Prepared Article 27 frozen inputs\n")
cat("  samples:", nrow(metadata), "\n")
cat("  class counts:", paste(names(table(metadata$Outcome)), table(metadata$Outcome), collapse = " / "), "\n")
cat("  species features:", nrow(species), "\n")
cat("  checksum payloads:", length(payloads), "\n")
