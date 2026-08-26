#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    paste(
      "Usage: prepare_article07_data.R",
      "<otutab.tsv> <taxonomy.tsv> <metadata.tsv>",
      "<result-dir> <summary.json>"
    ),
    call. = FALSE
  )
}

otutab_path <- normalizePath(args[[1]], mustWork = TRUE)
taxonomy_path <- normalizePath(args[[2]], mustWork = TRUE)
metadata_path <- normalizePath(args[[3]], mustWork = TRUE)
result_dir <- args[[4]]
summary_path <- args[[5]]

expected_hashes <- c(
  otutab = "45a093f93a1e9f83e341788c043174882142c91a27b409ff2b01d2696d83624f",
  taxonomy = "da5da1ec8b56b516056457edfd1758db57044bcf3cf7de2caef1bf5a809fcc0a",
  metadata = "e180f70324eb87cebdfb53b6ede10cf279a15a4e0124e82eda6d79c0d64344f2"
)
observed_hashes <- c(
  otutab = digest::digest(
    file = otutab_path,
    algo = "sha256",
    serialize = FALSE
  ),
  taxonomy = digest::digest(
    file = taxonomy_path,
    algo = "sha256",
    serialize = FALSE
  ),
  metadata = digest::digest(
    file = metadata_path,
    algo = "sha256",
    serialize = FALSE
  )
)
stopifnot(identical(observed_hashes, expected_hashes))

runtime_version <- as.character(utils::packageVersion("decontam"))
stopifnot(identical(runtime_version, "1.24.0"))
set.seed(20260719)

read_keyed_table <- function(path, key_name) {
  x <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stopifnot(
    ncol(x) >= 2L,
    identical(colnames(x)[[1]], key_name),
    !anyDuplicated(x[[1]])
  )
  rownames(x) <- x[[1]]
  x[[1]] <- NULL
  x
}

otutab_df <- read_keyed_table(otutab_path, "FeatureID")
taxonomy <- read_keyed_table(taxonomy_path, "FeatureID")
metadata <- read_keyed_table(metadata_path, "SampleID")

otutab <- as.matrix(otutab_df)
storage.mode(otutab) <- "numeric"
metadata$IsNegative <- as.logical(metadata$IsNegative)
metadata$PlateNumber <- factor(metadata$PlateNumber)

required_metadata <- c(
  "PlateNumber",
  "quant_reading",
  "Sample_or_Control",
  "IsNegative"
)
stopifnot(
  identical(rownames(otutab), rownames(taxonomy)),
  identical(colnames(otutab), rownames(metadata)),
  all(required_metadata %in% colnames(metadata)),
  all(is.finite(otutab)),
  all(otutab >= 0),
  all(otutab == floor(otutab)),
  all(colSums(otutab) > 0),
  all(is.finite(metadata$quant_reading)),
  all(metadata$quant_reading > 0),
  !anyNA(metadata$IsNegative)
)

control_by_plate <- with(
  metadata,
  table(PlateNumber, IsNegative)
)
stopifnot(
  nrow(otutab) == 1951L,
  ncol(otutab) == 569L,
  sum(!metadata$IsNegative) == 539L,
  sum(metadata$IsNegative) == 30L,
  nlevels(metadata$PlateNumber) == 6L,
  all(control_by_plate[, "TRUE"] == 5L),
  all(control_by_plate[, "FALSE"] > 0L)
)

seqtab <- t(otutab)
is_negative <- metadata[rownames(seqtab), "IsNegative"]
dna_concentration <- metadata[rownames(seqtab), "quant_reading"]
plate <- metadata[rownames(seqtab), "PlateNumber"]

contam_frequency <- decontam::isContaminant(
  seqtab,
  method = "frequency",
  conc = dna_concentration,
  threshold = 0.10
)
contam_prevalence <- decontam::isContaminant(
  seqtab,
  method = "prevalence",
  neg = is_negative,
  threshold = 0.10
)
contam_combined <- decontam::isContaminant(
  seqtab,
  method = "combined",
  conc = dna_concentration,
  neg = is_negative,
  threshold = 0.10
)
contam_batch <- decontam::isContaminant(
  seqtab,
  method = "combined",
  conc = dna_concentration,
  neg = is_negative,
  batch = plate,
  batch.combine = "minimum",
  threshold = 0.10
)

stopifnot(
  identical(rownames(contam_combined), rownames(taxonomy)),
  identical(rownames(contam_batch), rownames(taxonomy))
)

thresholds <- c(0.05, 0.10, 0.20)
threshold_results <- lapply(
  thresholds,
  function(threshold) {
    if (identical(threshold, 0.10)) {
      result <- contam_combined
    } else {
      result <- decontam::isContaminant(
        seqtab,
        method = "combined",
        conc = dna_concentration,
        neg = is_negative,
        threshold = threshold
      )
    }
    data.frame(
      Method = "Global combined",
      Threshold = threshold,
      ClassifiedFeatures = sum(result$contaminant),
      stringsAsFactors = FALSE
    )
  }
)
threshold_sensitivity <- do.call(rbind, threshold_results)

negative_prevalence <- colMeans(
  seqtab[is_negative, , drop = FALSE] > 0
)
biological_prevalence <- colMeans(
  seqtab[!is_negative, , drop = FALSE] > 0
)
total_feature_reads <- colSums(seqtab)

classification <- data.frame(
  FeatureID = rownames(taxonomy),
  taxonomy,
  NegativePrevalence = negative_prevalence,
  BiologicalPrevalence = biological_prevalence,
  TotalReads = total_feature_reads,
  FrequencyP = contam_frequency$p,
  PrevalenceP = contam_prevalence$p,
  CombinedFrequencyP = contam_combined$p.freq,
  CombinedPrevalenceP = contam_combined$p.prev,
  CombinedP = contam_combined$p,
  FrequencyContaminant = contam_frequency$contaminant,
  PrevalenceContaminant = contam_prevalence$contaminant,
  GlobalCombinedContaminant = contam_combined$contaminant,
  BatchCombinedContaminant = contam_batch$contaminant,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

classifier_summary <- data.frame(
  Analysis = c(
    "Frequency",
    "Prevalence",
    "Global combined",
    "Plate-aware combined",
    "Shared global and plate-aware"
  ),
  Threshold = c(0.10, 0.10, 0.10, 0.10, 0.10),
  ClassifiedFeatures = c(
    sum(contam_frequency$contaminant),
    sum(contam_prevalence$contaminant),
    sum(contam_combined$contaminant),
    sum(contam_batch$contaminant),
    sum(
      contam_combined$contaminant &
        contam_batch$contaminant
    )
  ),
  BatchStrategy = c(
    "Pooled",
    "Pooled",
    "Pooled",
    "PlateNumber; minimum p across batches",
    "Intersection"
  ),
  stringsAsFactors = FALSE
)

sample_class <- ifelse(
  metadata$IsNegative,
  "Negative control",
  "Biological sample"
)
library_size <- colSums(otutab)
contaminant_ids <- rownames(contam_combined)[
  contam_combined$contaminant
]
contaminant_reads <- colSums(
  otutab[contaminant_ids, , drop = FALSE]
)
sample_burden <- data.frame(
  SampleID = colnames(otutab),
  SampleClass = sample_class,
  Plate = as.character(metadata$PlateNumber),
  Subject = metadata$Subject,
  Habitat = metadata$Habitat,
  DNAConcentration = metadata$quant_reading,
  LibrarySize = library_size,
  ContaminantReads = contaminant_reads,
  ContaminantFraction = contaminant_reads / library_size,
  stringsAsFactors = FALSE
)

summarize_vector <- function(x) {
  c(
    minimum = min(x),
    median = stats::median(x),
    mean = mean(x),
    maximum = max(x)
  )
}

library_summary <- do.call(
  rbind,
  lapply(
    split(sample_burden$LibrarySize, sample_burden$SampleClass),
    summarize_vector
  )
)
burden_summary <- do.call(
  rbind,
  lapply(
    split(
      sample_burden$ContaminantFraction,
      sample_burden$SampleClass
    ),
    summarize_vector
  )
)
burden_correlation <- stats::cor.test(
  sample_burden$ContaminantFraction[
    sample_burden$SampleClass == "Biological sample"
  ],
  sample_burden$DNAConcentration[
    sample_burden$SampleClass == "Biological sample"
  ],
  method = "spearman",
  exact = FALSE
)

biological_ids <- rownames(metadata)[!metadata$IsNegative]
kept_features <- setdiff(rownames(otutab), contaminant_ids)
otutab_filtered <- otutab[
  kept_features,
  biological_ids,
  drop = FALSE
]
taxonomy_filtered <- taxonomy[
  kept_features,
  ,
  drop = FALSE
]
metadata_filtered <- metadata[
  biological_ids,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(otutab_filtered),
    rownames(taxonomy_filtered)
  ),
  identical(
    colnames(otutab_filtered),
    rownames(metadata_filtered)
  ),
  all(colSums(otutab_filtered) > 0)
)

biological_reads_before <- sum(
  otutab[, biological_ids, drop = FALSE]
)
biological_reads_after <- sum(otutab_filtered)
biological_reads_removed_fraction <- (
  1 - biological_reads_after / biological_reads_before
)

write_keyed_table <- function(x, key_name, path) {
  out <- data.frame(
    setNames(list(rownames(x)), key_name),
    x,
    check.names = FALSE
  )
  utils::write.table(
    out,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  classification,
  file.path(result_dir, "contaminant-classification.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
utils::write.table(
  classifier_summary,
  file.path(result_dir, "classifier-summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
utils::write.table(
  threshold_sensitivity,
  file.path(result_dir, "threshold-sensitivity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
utils::write.table(
  sample_burden,
  file.path(result_dir, "sample-contamination-burden.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
write_keyed_table(
  otutab_filtered,
  "FeatureID",
  file.path(result_dir, "otutab.filtered.tsv")
)
write_keyed_table(
  taxonomy_filtered,
  "FeatureID",
  file.path(result_dir, "taxonomy.filtered.tsv")
)
write_keyed_table(
  metadata_filtered,
  "SampleID",
  file.path(result_dir, "metadata.filtered.tsv")
)

summary <- list(
  source = "decontam::MUClite.rds exported triad",
  assay = "16S rRNA V4 ASV table used to demonstrate a metagenomics-capable classifier",
  decontam_runtime_version = runtime_version,
  source_otutab_sha256 = observed_hashes[["otutab"]],
  source_taxonomy_sha256 = observed_hashes[["taxonomy"]],
  source_metadata_sha256 = observed_hashes[["metadata"]],
  features = unname(nrow(otutab)),
  samples = unname(ncol(otutab)),
  biological_samples = unname(sum(!metadata$IsNegative)),
  negative_controls = unname(sum(metadata$IsNegative)),
  plates = unname(nlevels(metadata$PlateNumber)),
  controls_per_plate = unname(
    unique(as.integer(control_by_plate[, "TRUE"]))
  ),
  minimum_library_size = unname(min(library_size)),
  median_library_size = unname(stats::median(library_size)),
  maximum_library_size = unname(max(library_size)),
  biological_median_library_size = unname(
    library_summary["Biological sample", "median"]
  ),
  negative_control_median_library_size = unname(
    library_summary["Negative control", "median"]
  ),
  frequency_classified_features = unname(
    sum(contam_frequency$contaminant)
  ),
  prevalence_classified_features = unname(
    sum(contam_prevalence$contaminant)
  ),
  global_combined_classified_features = unname(
    sum(contam_combined$contaminant)
  ),
  plate_aware_classified_features = unname(
    sum(contam_batch$contaminant)
  ),
  global_plate_shared_features = unname(
    sum(
      contam_combined$contaminant &
        contam_batch$contaminant
    )
  ),
  biological_median_contaminant_fraction = unname(
    burden_summary["Biological sample", "median"]
  ),
  biological_mean_contaminant_fraction = unname(
    burden_summary["Biological sample", "mean"]
  ),
  biological_maximum_contaminant_fraction = unname(
    burden_summary["Biological sample", "maximum"]
  ),
  negative_control_median_contaminant_fraction = unname(
    burden_summary["Negative control", "median"]
  ),
  biological_burden_concentration_spearman_rho = unname(
    burden_correlation$estimate
  ),
  biological_burden_concentration_spearman_p = unname(
    burden_correlation$p.value
  ),
  retained_features = unname(nrow(otutab_filtered)),
  retained_biological_samples = unname(ncol(otutab_filtered)),
  biological_reads_removed_fraction = unname(
    biological_reads_removed_fraction
  ),
  random_seed = 20260719
)
jsonlite::write_json(
  summary,
  summary_path,
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 15
)

message(
  sprintf(
    paste(
      "Article 07: %d features x %d samples;",
      "%d global combined candidates;",
      "%d retained features."
    ),
    nrow(otutab),
    ncol(otutab),
    sum(contam_combined$contaminant),
    nrow(otutab_filtered)
  )
)
