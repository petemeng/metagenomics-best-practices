#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260725L
bootstrap_replicates <- 1000L
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
required <- c("project-root", "input-dir", "output-dir", "figure-dir")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c("ggplot2", "patchwork", "scales", "jsonlite", "digest")
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
input_dir <- normalizePath(args[["input-dir"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

font_cache_dir <- file.path(tempdir(), "article25-fontconfig-cache")
dir.create(font_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

validation_log <- file.path(output_dir, "validation.log")
writeLines(
  c(
    "Article 25 composition/core-microbiome validation",
    paste0("StartedUTC\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Seed\t", primary_seed),
    paste0("BootstrapReplicates\t", bootstrap_replicates)
  ),
  validation_log
)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), "\t", paste0(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = validation_log, append = TRUE)
  invisible(line)
}

checks <- data.frame(
  Category = character(), CheckID = character(), Status = character(),
  Detail = character(), stringsAsFactors = FALSE
)
add_check <- function(category, check_id, passed, detail) {
  checks <<- rbind(
    checks,
    data.frame(
      Category = category,
      CheckID = check_id,
      Status = if (isTRUE(passed)) "PASS" else "FAIL",
      Detail = paste(detail, collapse = "; "),
      stringsAsFactors = FALSE
    )
  )
  invisible(passed)
}
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

verify_checksum_manifest <- function(directory) {
  manifest_path <- file.path(directory, "file-checksums.sha256")
  lines <- readLines(manifest_path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  expected_files <- character()
  for (line in lines) {
    pieces <- strsplit(line, "[[:space:]]+", perl = TRUE)[[1L]]
    expected <- pieces[[1L]]
    relative <- paste(pieces[-1L], collapse = " ")
    path <- file.path(directory, relative)
    observed <- if (file.exists(path)) sha256_file(path) else "missing"
    add_check("Frozen input", paste0("sha256-", relative), identical(observed, expected), observed)
    expected_files <- c(expected_files, relative)
  }
  payloads <- sort(setdiff(basename(list.files(directory, full.names = TRUE)), "file-checksums.sha256"))
  add_check(
    "Frozen input", "checksum-manifest-complete",
    identical(payloads, sort(expected_files)),
    paste0("payloads=", length(payloads), "; entries=", length(expected_files))
  )
  length(lines)
}

read_profile <- function(path) {
  tab <- utils::read.delim(
    gzfile(path), check.names = FALSE, quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
  expected_metadata <- c("FeatureID", "Rank", "Label", "Lineage", "Phylum", "CoreEligible")
  stopifnot(identical(names(tab)[seq_along(expected_metadata)], expected_metadata))
  sample_ids <- names(tab)[-(seq_along(expected_metadata))]
  abundance <- data.matrix(tab[, sample_ids, drop = FALSE])
  rownames(abundance) <- tab$FeatureID
  list(
    abundance = abundance,
    features = tab[, expected_metadata, drop = FALSE],
    sample_ids = sample_ids
  )
}

wilson_interval <- function(successes, total, confidence = 0.95) {
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  p <- successes / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half_width <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  lower <- pmax(0, center - half_width)
  upper <- pmin(1, center + half_width)
  lower[successes == 0] <- 0
  upper[successes == total] <- 1
  cbind(Lower = lower, Upper = upper)
}

checksum_entries <- verify_checksum_manifest(input_dir)
notice_path <- file.path(project_root, "data", "small", "25-data-NOTICE.txt")
notice <- readLines(notice_path, warn = FALSE)
add_check("Frozen input", "notice-present", length(notice) >= 14L, notice_path)
add_check(
  "Frozen input", "notice-unit-boundary",
  any(grepl("primary inferential unit is one subject x selected HMP Figure 3 habitat", notice, fixed = TRUE)),
  "subject x selected HMP Figure 3 habitat"
)
add_check(
  "Frozen input", "notice-zero-boundary",
  any(grepl("does not prove biological absence", notice, fixed = TRUE)),
  "zero is a workflow-specific non-detection"
)
add_check(
  "Frozen input", "notice-reprocessing-boundary",
  any(grepl("not the original 2012 paper's exact taxonomic table", notice, fixed = TRUE)),
  "paper/profile version boundary"
)

metadata <- utils::read.delim(
  file.path(input_dir, "sample-metadata.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
selection_audit <- utils::read.delim(
  file.path(input_dir, "sample-selection-audit.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
closure_audit <- utils::read.delim(
  file.path(input_dir, "rank-closure-audit.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
analysis_contract <- utils::read.delim(
  file.path(input_dir, "analysis-contract.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
resource_manifest <- utils::read.delim(
  file.path(input_dir, "resource-manifest.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
profiles <- list(
  Phylum = read_profile(file.path(input_dir, "phylum-relative-abundance.tsv.gz")),
  Genus = read_profile(file.path(input_dir, "genus-relative-abundance.tsv.gz")),
  Species = read_profile(file.path(input_dir, "species-relative-abundance.tsv.gz"))
)

habitat_levels <- c(
  "Anterior nares", "Retroauricular crease", "Buccal mucosa",
  "Tongue dorsum", "Supragingival plaque", "Stool", "Posterior fornix"
)
metadata$Habitat <- factor(metadata$Habitat, levels = habitat_levels)
representative_ids <- metadata$SampleID[metadata$Representative]

add_check("Contract", "contract-row-count", nrow(analysis_contract) == 20L, nrow(analysis_contract))
add_check("Contract", "resource-row-count", nrow(resource_manifest) == 4L, nrow(resource_manifest))
contract_value <- function(item) analysis_contract$Value[match(item, analysis_contract$Item)]
add_check("Contract", "seed-locked", identical(contract_value("seed"), as.character(primary_seed)), contract_value("seed"))
add_check("Contract", "detection-threshold-locked", identical(contract_value("primary_detection_threshold"), "0.0001 fraction (0.01%)"), contract_value("primary_detection_threshold"))
add_check("Contract", "prevalence-threshold-locked", identical(contract_value("primary_prevalence_threshold"), "0.80"), contract_value("primary_prevalence_threshold"))
add_check("Contract", "bootstrap-locked", identical(contract_value("bootstrap_replicates"), as.character(bootstrap_replicates)), contract_value("bootstrap_replicates"))
add_check("Contract", "composition-not-pooled", grepl("never pool reads", contract_value("composition_summary"), fixed = TRUE), contract_value("composition_summary"))
add_check("Contract", "habitat-order-locked", identical(contract_value("habitat_order"), paste(habitat_levels, collapse = " | ")), contract_value("habitat_order"))
add_check("Contract", "bootstrap-unit-locked", identical(contract_value("bootstrap_unit"), "subject within habitat"), contract_value("bootstrap_unit"))

add_check("Input shape", "metadata-rows", nrow(metadata) == 748L, nrow(metadata))
add_check("Input shape", "independent-subjects", length(unique(metadata$SubjectID)) == 103L, length(unique(metadata$SubjectID)))
add_check("Input shape", "representative-units", length(representative_ids) == 490L, length(representative_ids))
representative_counts <- table(factor(metadata$Habitat[metadata$Representative], levels = habitat_levels))
add_check("Input shape", "representative-habitat-counts", identical(as.integer(representative_counts), c(74L, 17L, 84L, 90L, 88L, 95L, 42L)), paste(representative_counts, collapse = "/"))
add_check("Selection", "one-representative-per-subject-habitat", !anyDuplicated(paste(metadata$SubjectID[metadata$Representative], metadata$Habitat[metadata$Representative])), "490 unique subject-habitat keys")
add_check("Selection", "representatives-in-primary-habitats", all(metadata$PrimaryHabitatIncluded[metadata$Representative]) && !anyNA(metadata$Habitat[metadata$Representative]), "all representatives are in the seven declared habitats")
add_check("Selection", "selection-audit-aligned", identical(selection_audit$SampleID, metadata$SampleID) && identical(selection_audit$Representative, metadata$Representative), "selection/metadata")
add_check("Selection", "representatives-rank-one", all(selection_audit$CandidateRank[selection_audit$Representative] == 1L), range(selection_audit$CandidateRank[selection_audit$Representative]))

expected_features <- c(Phylum = 18L, Genus = 232L, Species = 749L)
for (rank_name in names(profiles)) {
  profile <- profiles[[rank_name]]
  add_check("Input shape", paste0(tolower(rank_name), "-shape"), identical(dim(profile$abundance), c(expected_features[[rank_name]], 748L)), paste(dim(profile$abundance), collapse = "x"))
  add_check("Alignment", paste0(tolower(rank_name), "-sample-order"), identical(profile$sample_ids, metadata$SampleID), "profile/metadata")
  add_check("Values", paste0(tolower(rank_name), "-finite-nonnegative"), all(is.finite(profile$abundance)) && min(profile$abundance) >= 0, range(profile$abundance))
  add_check("Values", paste0(tolower(rank_name), "-closure"), max(abs(colSums(profile$abundance) - 1)) < 1e-12, max(abs(colSums(profile$abundance) - 1)))
  add_check("Taxonomy", paste0(tolower(rank_name), "-rank-label"), all(profile$features$Rank == rank_name), unique(profile$features$Rank))
  add_check("Taxonomy", paste0(tolower(rank_name), "-unique-feature-id"), !anyDuplicated(profile$features$FeatureID), nrow(profile$features))
}
add_check(
  "Taxonomy", "species-core-eligibility-declared",
  is.logical(profiles$Species$features$CoreEligible) && all(profiles$Species$features$CoreEligible),
  paste0("eligible=", sum(profiles$Species$features$CoreEligible), "; excluded=", sum(!profiles$Species$features$CoreEligible))
)
add_check("Closure audit", "rank-closure-rows", identical(closure_audit$Rank, names(profiles)), paste(closure_audit$Rank, collapse = "/"))
add_check("Closure audit", "raw-closure-near-100", max(abs(closure_audit$RawSumPercentMedian - 100)) < 1e-12, closure_audit$RawSumPercentMedian)

# Composition summaries ---------------------------------------------------
log_msg("Computing rank-wise composition summaries")
rank_composition_rows <- list()
row_cursor <- 1L
for (rank_name in names(profiles)) {
  profile <- profiles[[rank_name]]
  rep_matrix <- profile$abundance[, representative_ids, drop = FALSE]
  for (habitat in habitat_levels) {
    habitat_ids <- metadata$SampleID[metadata$Representative & metadata$Habitat == habitat]
    habitat_matrix <- rep_matrix[, habitat_ids, drop = FALSE]
    rank_composition_rows[[row_cursor]] <- data.frame(
      Rank = rank_name,
      Habitat = habitat,
      Subjects = ncol(habitat_matrix),
      FeatureID = profile$features$FeatureID,
      Label = profile$features$Label,
      Phylum = profile$features$Phylum,
      MeanRelativeAbundance = rowMeans(habitat_matrix),
      MedianRelativeAbundance = apply(habitat_matrix, 1L, stats::median),
      DetectionPrevalence = rowMeans(habitat_matrix > 0),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    row_cursor <- row_cursor + 1L
  }
}
rank_composition_summary <- do.call(rbind, rank_composition_rows)
composition_sums <- aggregate(
  MeanRelativeAbundance ~ Rank + Habitat,
  rank_composition_summary,
  sum
)
add_check("Composition", "mean-composition-closure", max(abs(composition_sums$MeanRelativeAbundance - 1)) < 1e-12, max(abs(composition_sums$MeanRelativeAbundance - 1)))

# Prevalence, abundance and uncertainty ----------------------------------
primary_detection <- 1e-4
primary_prevalence <- 0.80
species <- profiles$Species
eligible <- species$features$CoreEligible

summarize_carriage <- function(sample_ids, habitat, unit_set) {
  x <- species$abundance[, sample_ids, drop = FALSE]
  detected <- x >= primary_detection
  successes <- rowSums(detected)
  interval <- wilson_interval(successes, ncol(x))
  positive_median <- vapply(
    seq_len(nrow(x)),
    function(i) {
      values <- x[i, detected[i, ], drop = TRUE]
      if (length(values) == 0L) 0 else stats::median(values)
    },
    numeric(1L)
  )
  data.frame(
    UnitSet = unit_set,
    Habitat = habitat,
    Units = ncol(x),
    FeatureID = species$features$FeatureID,
    Label = species$features$Label,
    Phylum = species$features$Phylum,
    CoreEligible = eligible,
    Detected = successes,
    Prevalence = successes / ncol(x),
    WilsonLower = interval[, "Lower"],
    WilsonUpper = interval[, "Upper"],
    MeanRelativeAbundance = rowMeans(x),
    MedianWhenDetected = positive_median,
    OperationalCore = eligible & successes / ncol(x) >= primary_prevalence,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

carriage_rows <- list()
cursor <- 1L
for (habitat in habitat_levels) {
  all_ids <- metadata$SampleID[metadata$PrimaryHabitatIncluded & metadata$Habitat == habitat]
  rep_ids <- metadata$SampleID[metadata$Representative & metadata$Habitat == habitat]
  carriage_rows[[cursor]] <- summarize_carriage(all_ids, habitat, "All samples")
  cursor <- cursor + 1L
  carriage_rows[[cursor]] <- summarize_carriage(rep_ids, habitat, "Subject-habitat representatives")
  cursor <- cursor + 1L
}
carriage_summary <- do.call(rbind, carriage_rows)
primary_core <- carriage_summary[carriage_summary$UnitSet == "Subject-habitat representatives", , drop = FALSE]

add_check("Carriage", "carriage-row-count", nrow(carriage_summary) == 2L * length(habitat_levels) * nrow(species$abundance), nrow(carriage_summary))
add_check("Carriage", "wilson-bounds", all(primary_core$WilsonLower >= 0 & primary_core$WilsonUpper <= 1 & primary_core$WilsonLower <= primary_core$Prevalence & primary_core$WilsonUpper >= primary_core$Prevalence), range(c(primary_core$WilsonLower, primary_core$WilsonUpper)))
add_check("Carriage", "nonzero-median-nonnegative", all(is.finite(primary_core$MedianWhenDetected)) && min(primary_core$MedianWhenDetected) >= 0, range(primary_core$MedianWhenDetected))

# Subject bootstrap preserves shared resampling weights across taxa.
log_msg("Running ", bootstrap_replicates, " subject bootstrap replicates per habitat")
bootstrap_rows <- list()
for (habitat_index in seq_along(habitat_levels)) {
  habitat <- habitat_levels[[habitat_index]]
  habitat_ids <- metadata$SampleID[metadata$Representative & metadata$Habitat == habitat]
  x <- species$abundance[eligible, habitat_ids, drop = FALSE]
  detected <- x >= primary_detection
  n_units <- ncol(detected)
  bootstrap_weights <- stats::rmultinom(
    bootstrap_replicates,
    size = n_units,
    prob = rep(1 / n_units, n_units)
  )
  bootstrap_counts <- detected %*% bootstrap_weights
  inclusion_frequency <- rowMeans(bootstrap_counts / n_units >= primary_prevalence)
  bootstrap_rows[[habitat_index]] <- data.frame(
    Habitat = habitat,
    Units = n_units,
    FeatureID = species$features$FeatureID[eligible],
    Label = species$features$Label[eligible],
    BootstrapReplicates = bootstrap_replicates,
    CoreInclusionFrequency = inclusion_frequency,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
bootstrap_stability <- do.call(rbind, bootstrap_rows)
primary_core$BootstrapInclusionFrequency <- bootstrap_stability$CoreInclusionFrequency[
  match(paste(primary_core$Habitat, primary_core$FeatureID), paste(bootstrap_stability$Habitat, bootstrap_stability$FeatureID))
]
primary_core$BootstrapInclusionFrequency[!primary_core$CoreEligible] <- NA_real_
primary_core$StableCore <- primary_core$OperationalCore & primary_core$BootstrapInclusionFrequency >= 0.80
add_check("Bootstrap", "bootstrap-row-count", nrow(bootstrap_stability) == length(habitat_levels) * sum(eligible), nrow(bootstrap_stability))
add_check("Bootstrap", "bootstrap-frequency-range", all(bootstrap_stability$CoreInclusionFrequency >= 0 & bootstrap_stability$CoreInclusionFrequency <= 1), range(bootstrap_stability$CoreInclusionFrequency))

# Threshold sensitivity grid ---------------------------------------------
detection_thresholds <- c(0, 1e-5, 1e-4, 1e-3)
prevalence_thresholds <- c(0.50, 0.80, 0.90, 1.00)
sensitivity_rows <- list()
cursor <- 1L
for (habitat in habitat_levels) {
  habitat_ids <- metadata$SampleID[metadata$Representative & metadata$Habitat == habitat]
  x <- species$abundance[eligible, habitat_ids, drop = FALSE]
  for (detection_threshold in detection_thresholds) {
    prevalence <- rowMeans(if (detection_threshold == 0) x > 0 else x >= detection_threshold)
    for (prevalence_threshold in prevalence_thresholds) {
      core_keep <- prevalence >= prevalence_threshold
      core_mass <- if (any(core_keep)) colSums(x[core_keep, , drop = FALSE]) else rep(0, ncol(x))
      sensitivity_rows[[cursor]] <- data.frame(
        Habitat = habitat,
        Units = ncol(x),
        DetectionThreshold = detection_threshold,
        PrevalenceThreshold = prevalence_threshold,
        CoreSpecies = sum(core_keep),
        MeanCoreMass = mean(core_mass),
        MedianCoreMass = stats::median(core_mass),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      cursor <- cursor + 1L
    }
  }
}
core_sensitivity <- do.call(rbind, sensitivity_rows)
add_check("Sensitivity", "grid-row-count", nrow(core_sensitivity) == length(habitat_levels) * length(detection_thresholds) * length(prevalence_thresholds), nrow(core_sensitivity))
add_check("Sensitivity", "grid-complete", !anyDuplicated(paste(core_sensitivity$Habitat, core_sensitivity$DetectionThreshold, core_sensitivity$PrevalenceThreshold)), "7 x 4 x 4")
add_check("Sensitivity", "core-mass-range", all(core_sensitivity$MeanCoreMass >= 0 & core_sensitivity$MeanCoreMass <= 1), range(core_sensitivity$MeanCoreMass))

# Core overlap and sample-vs-subject audit -------------------------------
core_species <- primary_core[primary_core$OperationalCore, , drop = FALSE]
core_overlap <- if (nrow(core_species) > 0L) {
  keys <- split(core_species$Habitat, core_species$FeatureID)
  out <- do.call(rbind, lapply(names(keys), function(feature_id) {
    habitats <- habitat_levels[habitat_levels %in% unique(keys[[feature_id]])]
    feature_index <- match(feature_id, species$features$FeatureID)
    data.frame(
      FeatureID = feature_id,
      Label = species$features$Label[[feature_index]],
      Phylum = species$features$Phylum[[feature_index]],
      CoreHabitats = length(habitats),
      Habitats = paste(habitats, collapse = " | "),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  out[order(-out$CoreHabitats, out$Label), , drop = FALSE]
} else {
  data.frame(
    FeatureID = character(), Label = character(), Phylum = character(),
    CoreHabitats = integer(), Habitats = character(), stringsAsFactors = FALSE
  )
}

unit_audit <- merge(
  carriage_summary[carriage_summary$UnitSet == "All samples", c("Habitat", "FeatureID", "Prevalence", "OperationalCore")],
  carriage_summary[carriage_summary$UnitSet == "Subject-habitat representatives", c("Habitat", "FeatureID", "Prevalence", "OperationalCore")],
  by = c("Habitat", "FeatureID"), suffixes = c("AllSamples", "Representatives"), sort = FALSE
)
unit_audit$AbsolutePrevalenceDifference <- abs(unit_audit$PrevalenceAllSamples - unit_audit$PrevalenceRepresentatives)
unit_audit$CoreClassificationChanged <- unit_audit$OperationalCoreAllSamples != unit_audit$OperationalCoreRepresentatives
add_check("Audit", "repeated-sample-impact-present", max(unit_audit$AbsolutePrevalenceDifference) > 0.01, max(unit_audit$AbsolutePrevalenceDifference))

# Figure data -------------------------------------------------------------
top_n_by_rank <- c(Phylum = 6L, Genus = 7L, Species = 7L)
composition_plot_rows <- list()
cursor <- 1L
top_taxa_audit <- list()
for (rank_name in names(profiles)) {
  profile <- profiles[[rank_name]]
  global_mean <- rowMeans(profile$abundance[, representative_ids, drop = FALSE])
  top_index <- head(order(global_mean, decreasing = TRUE), top_n_by_rank[[rank_name]])
  top_ids <- profile$features$FeatureID[top_index]
  top_taxa_audit[[rank_name]] <- data.frame(
    Rank = rank_name,
    FeatureID = top_ids,
    Label = profile$features$Label[top_index],
    GlobalMeanRelativeAbundance = global_mean[top_index],
    DisplayOrder = seq_along(top_index),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rank_summary <- rank_composition_summary[rank_composition_summary$Rank == rank_name, , drop = FALSE]
  for (habitat in habitat_levels) {
    habitat_summary <- rank_summary[rank_summary$Habitat == habitat, , drop = FALSE]
    selected <- habitat_summary[match(top_ids, habitat_summary$FeatureID), , drop = FALSE]
    composition_plot_rows[[cursor]] <- data.frame(
      Rank = rank_name,
      Habitat = habitat,
      Taxon = paste0(rank_name, " · ", selected$Label),
      MeanRelativeAbundance = selected$MeanRelativeAbundance,
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L
    composition_plot_rows[[cursor]] <- data.frame(
      Rank = rank_name,
      Habitat = habitat,
      Taxon = "Other",
      MeanRelativeAbundance = max(0, 1 - sum(selected$MeanRelativeAbundance)),
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L
  }
}
top_taxa_audit <- do.call(rbind, top_taxa_audit)
composition_plot_data <- do.call(rbind, composition_plot_rows)
composition_plot_data$Rank <- factor(composition_plot_data$Rank, levels = names(profiles))
composition_plot_data$Habitat <- factor(composition_plot_data$Habitat, levels = habitat_levels)

stool_ids <- metadata$SampleID[metadata$Representative & metadata$Habitat == "Stool"]
genus <- profiles$Genus
stool_genus <- genus$abundance[, stool_ids, drop = FALSE]
top_stool_index <- head(order(rowMeans(stool_genus), decreasing = TRUE), 12L)
top_stool_ids <- genus$features$FeatureID[top_stool_index]
dominant_id <- top_stool_ids[[1L]]
dominant_label <- genus$features$Label[match(dominant_id, genus$features$FeatureID)]
stool_order <- stool_ids[order(stool_genus[dominant_id, ], decreasing = TRUE, stool_ids)]
stool_plot_rows <- lapply(seq_along(top_stool_index), function(i) {
  feature_index <- top_stool_index[[i]]
  data.frame(
    SampleID = stool_order,
    SubjectOrder = seq_along(stool_order),
    Genus = genus$features$Label[[feature_index]],
    RelativeAbundance = stool_genus[feature_index, stool_order],
    stringsAsFactors = FALSE
  )
})
stool_plot_data <- do.call(rbind, stool_plot_rows)
top_stool_mass <- colSums(stool_genus[top_stool_index, stool_order, drop = FALSE])
stool_plot_data <- rbind(
  stool_plot_data,
  data.frame(
    SampleID = stool_order,
    SubjectOrder = seq_along(stool_order),
    Genus = "Other",
    RelativeAbundance = pmax(0, 1 - top_stool_mass),
    stringsAsFactors = FALSE
  )
)

selection_score <- primary_core$Prevalence * sqrt(primary_core$MeanRelativeAbundance)
primary_core$SelectionScore <- selection_score
candidate_ids <- unique(unlist(lapply(habitat_levels, function(habitat) {
  habitat_data <- primary_core[primary_core$Habitat == habitat & primary_core$CoreEligible, , drop = FALSE]
  habitat_data <- habitat_data[order(-habitat_data$OperationalCore, -habitat_data$SelectionScore, habitat_data$Label), , drop = FALSE]
  head(habitat_data$FeatureID, 6L)
})))
if (length(candidate_ids) > 28L) candidate_ids <- candidate_ids[seq_len(28L)]
selected_species <- species$features$Label[match(candidate_ids, species$features$FeatureID)]
paper_plot_data <- carriage_summary[carriage_summary$FeatureID %in% candidate_ids, , drop = FALSE]
species_order_score <- tapply(primary_core$SelectionScore, primary_core$FeatureID, max)
ordered_ids <- names(sort(species_order_score[candidate_ids], decreasing = FALSE))
paper_plot_data$Label <- factor(paper_plot_data$Label, levels = species$features$Label[match(ordered_ids, species$features$FeatureID)])
paper_plot_data$Habitat <- factor(paper_plot_data$Habitat, levels = habitat_levels)
paper_plot_data$UnitSet <- factor(paper_plot_data$UnitSet, levels = c("All samples", "Subject-habitat representatives"))

primary_core_counts <- aggregate(
  OperationalCore ~ Habitat,
  primary_core,
  sum
)
stable_core_counts <- aggregate(
  StableCore ~ Habitat,
  primary_core,
  sum
)
primary_core_counts <- merge(primary_core_counts, stable_core_counts, by = "Habitat", sort = FALSE)
primary_core_counts <- primary_core_counts[match(habitat_levels, primary_core_counts$Habitat), , drop = FALSE]

cross_habitat_ids <- core_overlap$FeatureID[core_overlap$CoreHabitats >= 2L]
top_core_ids_by_habitat <- unique(unlist(lapply(habitat_levels, function(habitat) {
  habitat_core <- primary_core[
    primary_core$Habitat == habitat & primary_core$OperationalCore,
    , drop = FALSE
  ]
  habitat_core <- habitat_core[
    order(-habitat_core$SelectionScore, habitat_core$Label),
    , drop = FALSE
  ]
  head(habitat_core$FeatureID, 5L)
})))
core_display_ids <- unique(c(cross_habitat_ids, top_core_ids_by_habitat))
core_display_audit <- core_overlap[
  match(core_display_ids, core_overlap$FeatureID),
  , drop = FALSE
]
core_display_audit$DisplayReason <- ifelse(
  core_display_audit$CoreHabitats >= 2L,
  "core in at least two habitats",
  "among five highest prevalence-abundance scores in a core habitat"
)
core_membership_plot_data <- primary_core[primary_core$FeatureID %in% core_display_ids, , drop = FALSE]
if (nrow(core_membership_plot_data) > 0L) {
  core_order <- core_overlap$FeatureID
  core_membership_plot_data$Label <- factor(
    core_membership_plot_data$Label,
    levels = rev(species$features$Label[match(core_order, species$features$FeatureID)])
  )
  core_membership_plot_data$Habitat <- factor(core_membership_plot_data$Habitat, levels = habitat_levels)
}
core_sensitivity$Habitat <- factor(core_sensitivity$Habitat, levels = habitat_levels)
core_sensitivity$DetectionLabel <- factor(
  core_sensitivity$DetectionThreshold,
  levels = detection_thresholds,
  labels = c("Any > 0", "0.001%", "0.01%", "0.1%")
)
core_sensitivity$PrevalenceLabel <- factor(
  core_sensitivity$PrevalenceThreshold,
  levels = rev(prevalence_thresholds),
  labels = scales::percent(rev(prevalence_thresholds), accuracy = 1)
)

# Publication graphics ----------------------------------------------------
theme_pub <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", colour = "#222222"),
      plot.subtitle = ggplot2::element_text(colour = "#444444"),
      axis.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = "#BDBDBD"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      plot.caption = ggplot2::element_text(colour = "#555555", hjust = 0)
    )
}
save_pub <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 350, bg = "white")
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height, units = "in", dpi = 350, compression = "lzw", bg = "white")
}

composition_taxa <- unique(composition_plot_data$Taxon)
composition_taxa <- c(setdiff(composition_taxa, "Other"), "Other")
composition_palette <- c(
  stats::setNames(grDevices::hcl.colors(length(composition_taxa) - 1L, palette = "Dark 3"), composition_taxa[-length(composition_taxa)]),
  Other = "#D9D9D9"
)
composition_plot_data$Taxon <- factor(composition_plot_data$Taxon, levels = composition_taxa)
composition_figure <- ggplot2::ggplot(
  composition_plot_data,
  ggplot2::aes(x = Habitat, y = MeanRelativeAbundance, fill = Taxon)
) +
  ggplot2::geom_col(width = 0.76, colour = "white", linewidth = 0.18) +
  ggplot2::facet_wrap(~ Rank, nrow = 1) +
  ggplot2::scale_fill_manual(values = composition_palette, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = c(0, 0)) +
  ggplot2::labs(
    title = "Community composition changes with taxonomic rank",
    subtitle = "Mean of 490 subject-habitat compositions; each rank is reclosed separately",
    x = NULL, y = "Mean relative abundance", fill = NULL,
    caption = "Top taxa are selected globally within each rank. ‘Other’ is calculated after sample-wise closure; reads are never pooled across subjects."
  ) +
  theme_pub(9.3) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 32, hjust = 1), legend.key.height = grid::unit(0.35, "cm")) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 4, byrow = TRUE))
save_pub(composition_figure, "25-multirank-mean-composition", 14.2, 8.8)

stool_categories <- c(genus$features$Label[top_stool_index], "Other")
stool_palette <- c(
  stats::setNames(grDevices::hcl.colors(length(stool_categories) - 1L, palette = "Dark 3"), stool_categories[-length(stool_categories)]),
  Other = "#D9D9D9"
)
stool_plot_data$Genus <- factor(stool_plot_data$Genus, levels = stool_categories)
stool_figure <- ggplot2::ggplot(
  stool_plot_data,
  ggplot2::aes(x = SubjectOrder, y = RelativeAbundance, fill = Genus)
) +
  ggplot2::geom_col(width = 1, colour = NA) +
  ggplot2::scale_fill_manual(values = stool_palette, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  ggplot2::labs(
    title = "A group mean can hide individualized stool profiles",
    subtitle = paste0("One depth-maximized stool sample per subject (n = 95), ordered by ", dominant_label, " abundance"),
    x = "Subjects", y = "Genus relative abundance", fill = NULL,
    caption = "The twelve genera with the largest mean abundance are shown; all remaining genera are combined as Other."
  ) +
  theme_pub(9.6) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
save_pub(stool_figure, "25-individual-stool-composition", 12.6, 6.8)

paper_figure <- ggplot2::ggplot(
  paper_plot_data,
  ggplot2::aes(x = Habitat, y = Label, size = MedianWhenDetected * 100, fill = Prevalence)
) +
  ggplot2::geom_point(shape = 21, colour = "#333333", stroke = 0.35) +
  ggplot2::facet_wrap(~ UnitSet, ncol = 1) +
  ggplot2::scale_fill_viridis_c(option = "C", limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  ggplot2::scale_size_continuous(range = c(1.2, 9), trans = "sqrt", breaks = c(0.01, 0.1, 1, 10), name = "Median when detected (%)") +
  ggplot2::labs(
    title = "Species carriage depends on the observational unit",
    subtitle = "Detection requires relative abundance ≥0.01%; colour is prevalence and area is non-zero median abundance",
    x = NULL, y = NULL, fill = "Prevalence",
    caption = "The upper panel counts repeated samples; the lower panel counts one representative per subject and habitat. Zero is a profile-specific non-detection."
  ) +
  theme_pub(9.1) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), panel.grid.major.x = ggplot2::element_line(colour = "#E8E8E8"))
save_pub(paper_figure, "25-prevalence-abundance", 11.8, 13.8)

if (nrow(core_membership_plot_data) > 0L) {
  membership_figure <- ggplot2::ggplot(
    core_membership_plot_data,
    ggplot2::aes(x = Habitat, y = Label, fill = Prevalence, alpha = OperationalCore)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(OperationalCore, scales::percent(Prevalence, accuracy = 1), "")), size = 2.5) +
    ggplot2::scale_fill_viridis_c(option = "C", limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.22, `TRUE` = 1), guide = "none") +
    ggplot2::labs(
      title = "A  Selected operational-core membership",
      subtitle = "All multi-habitat members plus five highest-scoring members per habitat",
      x = NULL, y = NULL, fill = "Prevalence"
    ) +
    theme_pub(8.8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
} else {
  membership_figure <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No species met the primary core rule") +
    ggplot2::theme_void() +
    ggplot2::labs(title = "A  Primary operational core")
}

sensitivity_figure <- ggplot2::ggplot(
  core_sensitivity,
  ggplot2::aes(x = DetectionLabel, y = PrevalenceLabel, fill = CoreSpecies)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
  ggplot2::geom_text(ggplot2::aes(label = CoreSpecies), size = 2.7) +
  ggplot2::facet_wrap(~ Habitat, ncol = 3) +
  ggplot2::scale_fill_viridis_c(option = "D") +
  ggplot2::labs(
    title = "B  Core membership is threshold-sensitive",
    subtitle = "Numbers are eligible species retained by each definition",
    x = "Detection threshold", y = "Prevalence threshold", fill = "Core species"
  ) +
  theme_pub(8.8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

core_figure <- membership_figure | sensitivity_figure
core_figure <- core_figure + patchwork::plot_layout(widths = c(1.2, 1)) +
  patchwork::plot_annotation(
    title = "Core microbiome is an operational definition, not a fixed biological list",
    caption = "Primary membership requires species ≥0.01% in ≥80% of subject-habitat units. Panel A is a display subset; all 92 members, Wilson intervals and 1,000-bootstrap stability are retained in the result tables."
  )
save_pub(core_figure, "25-core-membership-sensitivity", 18.0, 12.4)

# Result tables and final validation -------------------------------------
write_tsv(rank_composition_summary, file.path(output_dir, "rank-composition-summary.tsv"))
write_tsv_gz(carriage_summary, file.path(output_dir, "carriage-summary.tsv.gz"))
write_tsv(primary_core, file.path(output_dir, "primary-core-membership.tsv"))
write_tsv(bootstrap_stability, file.path(output_dir, "bootstrap-core-stability.tsv"))
write_tsv(core_sensitivity, file.path(output_dir, "core-threshold-sensitivity.tsv"))
write_tsv(core_overlap, file.path(output_dir, "core-overlap.tsv"))
write_tsv(core_display_audit, file.path(output_dir, "core-display-selection.tsv"))
write_tsv(unit_audit, file.path(output_dir, "sample-vs-subject-habitat-audit.tsv"))
write_tsv(top_taxa_audit, file.path(output_dir, "top-taxa-audit.tsv"))

package_versions <- data.frame(
  Package = c("R", packages),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    vapply(packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1L))
  ),
  stringsAsFactors = FALSE
)
write_tsv(package_versions, file.path(output_dir, "package-versions.tsv"))

primary_counts_vector <- stats::setNames(primary_core_counts$OperationalCore, primary_core_counts$Habitat)
stable_counts_vector <- stats::setNames(primary_core_counts$StableCore, primary_core_counts$Habitat)
add_check("Primary core", "core-counts-nonnegative", all(primary_counts_vector >= 0), primary_counts_vector)
add_check("Primary core", "stable-subset", all(stable_counts_vector <= primary_counts_vector), paste(stable_counts_vector, primary_counts_vector, sep = "/"))
add_check("Primary core", "core-overlap-aligned", nrow(core_overlap) == length(unique(core_species$FeatureID)), paste(nrow(core_overlap), length(unique(core_species$FeatureID)), sep = "/"))
add_check("Primary core", "display-includes-cross-habitat-members", all(cross_habitat_ids %in% core_display_audit$FeatureID), paste(length(cross_habitat_ids), nrow(core_display_audit), sep = "/"))
add_check("Primary core", "display-subset-readable", nrow(core_display_audit) <= 55L, nrow(core_display_audit))
add_check("Primary core", "retroauricular-small-n-visible", representative_counts[["Retroauricular crease"]] == 17L, representative_counts[["Retroauricular crease"]])
add_check(
  "Primary core", "primary-core-counts-locked",
  identical(as.integer(primary_counts_vector), c(2L, 5L, 21L, 50L, 34L, 9L, 0L)),
  paste(primary_counts_vector, collapse = "/")
)
add_check(
  "Primary core", "stable-core-counts-locked",
  identical(as.integer(stable_counts_vector), c(2L, 4L, 20L, 45L, 30L, 6L, 0L)),
  paste(stable_counts_vector, collapse = "/")
)
add_check("Primary core", "cross-habitat-core-count-locked", nrow(core_overlap) == 92L, nrow(core_overlap))
add_check("Audit", "unit-core-classification-changes-locked", sum(unit_audit$CoreClassificationChanged) == 9L, sum(unit_audit$CoreClassificationChanged))
add_check(
  "Audit", "unit-prevalence-shift-locked",
  abs(max(unit_audit$AbsolutePrevalenceDifference) - 0.0871459694989107) < 1e-12,
  max(unit_audit$AbsolutePrevalenceDifference)
)

result_files <- c(
  "rank-composition-summary.tsv", "carriage-summary.tsv.gz",
  "primary-core-membership.tsv", "bootstrap-core-stability.tsv",
  "core-threshold-sensitivity.tsv", "core-overlap.tsv", "core-display-selection.tsv",
  "sample-vs-subject-habitat-audit.tsv", "top-taxa-audit.tsv",
  "package-versions.tsv", "validation.log"
)
for (file_name in result_files) {
  path <- file.path(output_dir, file_name)
  add_check("Outputs", paste0("result-", file_name), file.exists(path) && file.info(path)$size > 0, if (file.exists(path)) file.info(path)$size else "missing")
}
figure_stems <- c(
  "25-multirank-mean-composition", "25-individual-stool-composition",
  "25-prevalence-abundance", "25-core-membership-sensitivity"
)
for (stem in figure_stems) {
  for (extension in c("pdf", "png", "tiff")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    add_check("Figures", paste0(stem, "-", extension), file.exists(path) && file.info(path)$size > 1000, if (file.exists(path)) file.info(path)$size else "missing")
  }
}
original_figure_path <- file.path(project_root, "figures", "25-hmp-fig3-original.jpg")
add_check("Paper anchor", "original-figure-present", file.exists(original_figure_path) && file.info(original_figure_path)$size > 50000, if (file.exists(original_figure_path)) file.info(original_figure_path)$size else "missing")
add_check("Paper anchor", "original-figure-sha256", file.exists(original_figure_path) && identical(sha256_file(original_figure_path), "5afd411dee3ad8604263da4a3069f429942b6605674ce7cc8a3d0a161dca66c2"), if (file.exists(original_figure_path)) sha256_file(original_figure_path) else "missing")

write_tsv(checks, file.path(output_dir, "validation-audit.tsv"))
checks_failed <- sum(checks$Status == "FAIL")
summary <- list(
  status = if (checks_failed == 0L) "passed" else "failed",
  seed = primary_seed,
  source_samples = nrow(metadata),
  independent_subjects = length(unique(metadata$SubjectID)),
  representative_subject_habitat_units = length(representative_ids),
  representative_counts = as.list(as.integer(representative_counts)),
  representative_count_labels = names(representative_counts),
  phylum_features = nrow(profiles$Phylum$abundance),
  genus_features = nrow(profiles$Genus$abundance),
  species_features = nrow(profiles$Species$abundance),
  core_eligible_species = sum(eligible),
  primary_detection_threshold = primary_detection,
  primary_prevalence_threshold = primary_prevalence,
  primary_core_counts = as.list(as.integer(primary_counts_vector)),
  primary_core_count_labels = names(primary_counts_vector),
  stable_core_counts = as.list(as.integer(stable_counts_vector)),
  stable_core_count_labels = names(stable_counts_vector),
  cross_habitat_core_species = nrow(core_overlap),
  multi_habitat_core_species = length(cross_habitat_ids),
  core_display_species = nrow(core_display_audit),
  max_core_habitats = if (nrow(core_overlap) == 0L) 0L else max(core_overlap$CoreHabitats),
  max_sample_vs_subject_prevalence_difference = max(unit_audit$AbsolutePrevalenceDifference),
  core_classification_changes = sum(unit_audit$CoreClassificationChanged),
  sensitivity_branches = nrow(core_sensitivity),
  bootstrap_replicates = bootstrap_replicates,
  checksum_entries = checksum_entries,
  checks_total = nrow(checks),
  checks_passed = sum(checks$Status == "PASS"),
  checks_failed = checks_failed,
  r_version = package_versions$Version[package_versions$Package == "R"],
  ggplot2_version = package_versions$Version[package_versions$Package == "ggplot2"]
)
jsonlite::write_json(
  summary,
  file.path(output_dir, "validation-summary.json"),
  pretty = TRUE, auto_unbox = TRUE, na = "null"
)
log_msg("Validation complete: status=", summary$status, "; checks=", summary$checks_passed, "/", summary$checks_total)

if (checks_failed > 0L) {
  failed <- checks[checks$Status == "FAIL", , drop = FALSE]
  stop(
    "Article 25 validation failed: ",
    paste(paste0(failed$CheckID, " [", failed$Detail, "]"), collapse = "; "),
    call. = FALSE
  )
}
