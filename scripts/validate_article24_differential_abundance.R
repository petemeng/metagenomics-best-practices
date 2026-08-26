#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260724L
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

packages <- c(
  "maaslin3", "ANCOMBC", "ALDEx2", "ggplot2", "patchwork",
  "scales", "jsonlite", "digest"
)
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
font_cache_dir <- file.path(tempdir(), "article24-fontconfig-cache")
dir.create(font_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

validation_log <- file.path(output_dir, "validation.log")
writeLines(
  c(
    "Article 24 differential-abundance validation",
    paste0("StartedUTC\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Seed\t", primary_seed)
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
z_score <- function(x) as.numeric(scale(as.numeric(x)))

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
    add_check(
      "Frozen input", paste0("sha256-", relative),
      identical(observed, expected), observed
    )
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

read_species_profile <- function(path) {
  tab <- utils::read.delim(
    gzfile(path), check.names = FALSE, quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
  stopifnot(identical(names(tab)[1:2], c("Feature", "Species")))
  sample_ids <- names(tab)[-(1:2)]
  abundance <- data.matrix(tab[, sample_ids, drop = FALSE])
  rownames(abundance) <- sprintf("SP%04d", seq_len(nrow(tab)))
  list(
    abundance = abundance,
    map = data.frame(
      FeatureID = rownames(abundance), Feature = tab$Feature,
      Species = tab$Species, stringsAsFactors = FALSE
    )
  )
}

read_pathway_profile <- function(path) {
  tab <- utils::read.delim(
    gzfile(path), check.names = FALSE, quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
  stopifnot(identical(names(tab)[1L], "Pathway"))
  sample_ids <- names(tab)[-1L]
  abundance <- data.matrix(tab[, sample_ids, drop = FALSE])
  rownames(abundance) <- sprintf("PW%04d", seq_len(nrow(tab)))
  list(
    abundance = abundance,
    map = data.frame(
      FeatureID = rownames(abundance), Feature = tab$Pathway,
      Pathway = tab$Pathway, stringsAsFactors = FALSE
    )
  )
}

checksum_entries <- verify_checksum_manifest(input_dir)
analysis_contract <- utils::read.delim(
  file.path(input_dir, "analysis-contract.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
resource_manifest <- utils::read.delim(
  file.path(input_dir, "resource-manifest.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
notice <- readLines(file.path(input_dir, "data-NOTICE.txt"), warn = FALSE)

contract_value <- function(item) {
  analysis_contract$value[match(item, analysis_contract$item)]
}
add_check("Contract", "contract-row-count", nrow(analysis_contract) == 15L, nrow(analysis_contract))
add_check("Contract", "resource-row-count", nrow(resource_manifest) == 5L, nrow(resource_manifest))
add_check("Contract", "seed-locked", identical(contract_value("seed"), as.character(primary_seed)), contract_value("seed"))
add_check("Contract", "primary-n-locked", identical(contract_value("primary_complete_cases"), "110"), contract_value("primary_complete_cases"))
add_check("Contract", "aldex-mc-locked", identical(contract_value("aldex2_mc_samples"), "128"), contract_value("aldex2_mc_samples"))
add_check(
  "Contract", "pseudocount-boundary-in-notice",
  any(grepl("not observed taxon reads", notice, fixed = TRUE)),
  "Pseudo-count boundary is explicit"
)
add_check(
  "Contract", "pathway-denominator-in-notice",
  any(grepl("ordinary unstratified pathways", notice, fixed = TRUE)),
  "Pathway closure boundary is explicit"
)

species_input <- read_species_profile(file.path(input_dir, "species-relative-abundance.tsv.gz"))
pseudocount_input <- read_species_profile(file.path(input_dir, "species-pseudocounts.tsv.gz"))
pathway_input <- read_pathway_profile(file.path(input_dir, "pathway-relative-abundance.tsv.gz"))
metadata_all <- utils::read.delim(
  file.path(input_dir, "sample-metadata.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)

sample_ids <- colnames(species_input$abundance)
add_check("Input shape", "species-shape", identical(dim(species_input$abundance), c(661L, 156L)), paste(dim(species_input$abundance), collapse = "x"))
add_check("Input shape", "pseudocount-shape", identical(dim(pseudocount_input$abundance), c(661L, 156L)), paste(dim(pseudocount_input$abundance), collapse = "x"))
add_check("Input shape", "pathway-shape", identical(dim(pathway_input$abundance), c(493L, 156L)), paste(dim(pathway_input$abundance), collapse = "x"))
add_check("Input shape", "metadata-rows", nrow(metadata_all) == 156L, nrow(metadata_all))
add_check(
  "Alignment", "sample-order",
  identical(sample_ids, colnames(pseudocount_input$abundance)) &&
    identical(sample_ids, colnames(pathway_input$abundance)) &&
    identical(sample_ids, metadata_all$sample_id),
  "species/pseudocount/pathway/metadata"
)
add_check(
  "Alignment", "species-feature-order",
  identical(species_input$map$Feature, pseudocount_input$map$Feature) &&
    identical(species_input$map$Species, pseudocount_input$map$Species),
  "relative abundance/pseudo-count"
)
add_check("Values", "species-finite-nonnegative", all(is.finite(species_input$abundance)) && min(species_input$abundance) >= 0, range(species_input$abundance))
add_check("Values", "pathway-finite-nonnegative", all(is.finite(pathway_input$abundance)) && min(pathway_input$abundance) >= 0, range(pathway_input$abundance))
add_check("Values", "pseudocount-integer-nonnegative", all(is.finite(pseudocount_input$abundance)) && min(pseudocount_input$abundance) >= 0 && max(abs(pseudocount_input$abundance - round(pseudocount_input$abundance))) == 0, range(pseudocount_input$abundance))
add_check("Values", "species-closure", max(abs(colSums(species_input$abundance) - 1)) < 2e-6, max(abs(colSums(species_input$abundance) - 1)))
add_check("Values", "pathway-closure", max(abs(colSums(pathway_input$abundance) - 1)) < 2e-7, max(abs(colSums(pathway_input$abundance) - 1)))
reconstructed <- round(sweep(
  species_input$abundance, 2L, metadata_all$number_reads, "*"
))
reconstruction_error <- max(abs(reconstructed - pseudocount_input$abundance))
add_check("Values", "pseudocount-reconstruction", reconstruction_error == 0, reconstruction_error)

primary_keep <- metadata_all$primary_complete_case
primary_metadata_raw <- metadata_all[primary_keep, , drop = FALSE]
primary_ids <- primary_metadata_raw$sample_id
add_check("Design", "primary-sample-count", length(primary_ids) == 110L, length(primary_ids))
add_check("Design", "independent-subject-count", length(unique(primary_metadata_raw$subject_id)) == 110L, length(unique(primary_metadata_raw$subject_id)))
add_check(
  "Design", "primary-group-counts",
  identical(
    as.integer(table(factor(primary_metadata_raw$analysis_group, levels = c("Control", "CRC")))),
    c(59L, 51L)
  ),
  paste(table(primary_metadata_raw$analysis_group), collapse = "/")
)
add_check(
  "Design", "primary-covariates-complete",
  all(stats::complete.cases(primary_metadata_raw[, c("age", "gender", "BMI", "number_reads")])),
  "age/gender/BMI/number_reads"
)

make_metadata <- function(raw, include_bmi = TRUE, progression = FALSE) {
  out <- data.frame(
    disease = factor(raw$analysis_group, levels = c("Control", "CRC")),
    z_age = z_score(raw$age),
    gender = factor(raw$gender, levels = c("female", "male")),
    z_log10_reads = z_score(log10(raw$number_reads)),
    row.names = raw$sample_id,
    check.names = FALSE
  )
  if (include_bmi) out$z_BMI <- z_score(raw$BMI)
  if (progression) {
    out$disease <- NULL
    out$progression_score <- unname(c(Control = 0, Adenoma = 1, CRC = 2)[raw$analysis_group])
  }
  out
}

primary_metadata <- make_metadata(primary_metadata_raw, include_bmi = TRUE)
primary_metadata <- primary_metadata[, c("disease", "z_age", "gender", "z_BMI", "z_log10_reads")]
primary_model_matrix <- stats::model.matrix(
  ~ disease + z_age + gender + z_BMI + z_log10_reads,
  data = primary_metadata
)
add_check(
  "Design", "primary-model-full-rank",
  qr(primary_model_matrix)$rank == ncol(primary_model_matrix),
  paste0("rank=", qr(primary_model_matrix)$rank, "; columns=", ncol(primary_model_matrix))
)
add_check("Design", "disease-reference-control", levels(primary_metadata$disease)[1L] == "Control", levels(primary_metadata$disease))
add_check("Design", "gender-reference-female", levels(primary_metadata$gender)[1L] == "female", levels(primary_metadata$gender))

species_primary_all <- species_input$abundance[, primary_ids, drop = FALSE]
pseudocount_primary_all <- pseudocount_input$abundance[, primary_ids, drop = FALSE]
pathway_primary_all <- pathway_input$abundance[, primary_ids, drop = FALSE]

species_abundance_threshold <- 0.0001
species_prevalence_threshold <- 0.10
pathway_abundance_threshold <- 0.00001
pathway_prevalence_threshold <- 0.20

species_threshold_prevalence <- rowMeans(species_primary_all > species_abundance_threshold)
species_detected_prevalence <- rowMeans(species_primary_all > 0)
species_keep <- species_threshold_prevalence >= species_prevalence_threshold
pathway_threshold_prevalence <- rowMeans(pathway_primary_all > pathway_abundance_threshold)
pathway_detected_prevalence <- rowMeans(pathway_primary_all > 0)
pathway_keep <- pathway_threshold_prevalence >= pathway_prevalence_threshold

species_group_counts <- t(vapply(
  seq_len(nrow(species_primary_all)),
  function(i) {
    present <- species_primary_all[i, ] > 0
    c(
      ControlPresent = sum(present[primary_metadata$disease == "Control"]),
      ControlAbsent = sum(!present[primary_metadata$disease == "Control"]),
      CRCPresent = sum(present[primary_metadata$disease == "CRC"]),
      CRCAbsent = sum(!present[primary_metadata$disease == "CRC"])
    )
  },
  numeric(4L)
))
species_prevalence_estimable_all <- apply(species_group_counts, 1L, function(x) all(x >= 10L))
species_abundance_estimable_all <- species_group_counts[, "ControlPresent"] >= 10L & species_group_counts[, "CRCPresent"] >= 10L

pathway_group_counts <- t(vapply(
  seq_len(nrow(pathway_primary_all)),
  function(i) {
    present <- pathway_primary_all[i, ] > 0
    c(
      ControlPresent = sum(present[primary_metadata$disease == "Control"]),
      ControlAbsent = sum(!present[primary_metadata$disease == "Control"]),
      CRCPresent = sum(present[primary_metadata$disease == "CRC"]),
      CRCAbsent = sum(!present[primary_metadata$disease == "CRC"])
    )
  },
  numeric(4L)
))
pathway_prevalence_estimable_all <- apply(pathway_group_counts, 1L, function(x) all(x >= 10L))
pathway_abundance_estimable_all <- pathway_group_counts[, "ControlPresent"] >= 10L & pathway_group_counts[, "CRCPresent"] >= 10L

add_check("Feature universe", "species-feature-count", sum(species_keep) == 212L, sum(species_keep))
add_check("Feature universe", "pathway-feature-count", sum(pathway_keep) == 394L, sum(pathway_keep))
add_check("Feature universe", "species-prevalence-estimable", sum(species_prevalence_estimable_all[species_keep]) == 144L, sum(species_prevalence_estimable_all[species_keep]))
add_check("Feature universe", "species-filter-predeclared", identical(species_abundance_threshold, 0.0001) && identical(species_prevalence_threshold, 0.10), paste(species_abundance_threshold, species_prevalence_threshold))
add_check("Feature universe", "pathway-filter-predeclared", identical(pathway_abundance_threshold, 0.00001) && identical(pathway_prevalence_threshold, 0.20), paste(pathway_abundance_threshold, pathway_prevalence_threshold))

feature_filter_audit <- rbind(
  data.frame(
    FeatureSpace = "Species",
    FeatureID = species_input$map$FeatureID,
    Feature = species_input$map$Feature,
    Label = species_input$map$Species,
    AbundanceThreshold = species_abundance_threshold,
    RequiredPrevalence = species_prevalence_threshold,
    ObservedThresholdPrevalence = species_threshold_prevalence,
    ObservedDetectionPrevalence = species_detected_prevalence,
    MaxRelativeAbundance = apply(species_primary_all, 1L, max),
    PrimaryUniverse = species_keep,
    AbundanceEstimable = species_abundance_estimable_all,
    PrevalenceEstimable = species_prevalence_estimable_all,
    species_group_counts,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  data.frame(
    FeatureSpace = "Pathway",
    FeatureID = pathway_input$map$FeatureID,
    Feature = pathway_input$map$Feature,
    Label = pathway_input$map$Pathway,
    AbundanceThreshold = pathway_abundance_threshold,
    RequiredPrevalence = pathway_prevalence_threshold,
    ObservedThresholdPrevalence = pathway_threshold_prevalence,
    ObservedDetectionPrevalence = pathway_detected_prevalence,
    MaxRelativeAbundance = apply(pathway_primary_all, 1L, max),
    PrimaryUniverse = pathway_keep,
    AbundanceEstimable = pathway_abundance_estimable_all,
    PrevalenceEstimable = pathway_prevalence_estimable_all,
    pathway_group_counts,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)

design_balance <- do.call(rbind, lapply(c("Control", "CRC"), function(group) {
  d <- primary_metadata_raw[primary_metadata_raw$analysis_group == group, , drop = FALSE]
  data.frame(
    Group = group,
    Subjects = nrow(d),
    Female = sum(d$gender == "female"),
    Male = sum(d$gender == "male"),
    AgeMean = mean(d$age), AgeSD = stats::sd(d$age),
    BMIMean = mean(d$BMI), BMISD = stats::sd(d$BMI),
    ReadsMedian = stats::median(d$number_reads),
    ReadsIQR = stats::IQR(d$number_reads),
    stringsAsFactors = FALSE
  )
}))

data_lineage <- data.frame(
  Stage = c(
    "Publisher taxonomic payload", "Publisher functional payload",
    "Frozen species profile", "Frozen pathway profile", "Reconstructed pseudo-count"
  ),
  Source = c(
    "curatedMetagenomicData 2021-03-31 ZellerG_2014 relative_abundance",
    "curatedMetagenomicData 2021-03-31 ZellerG_2014 pathway_abundance",
    "MetaPhlAn3 species rows only",
    "Unstratified ordinary MetaCyc rows reclosed to one",
    "round(species relative fraction x whole-metagenome reads)"
  ),
  SourceRows = c(1019L, 22620L, 661L, 493L, 661L),
  Samples = 156L,
  PrimarySubjects = c(NA_integer_, NA_integer_, 110L, 110L, 110L),
  Unit = c(
    "Percent repeated across seven taxonomic ranks",
    "Community plus stratified HUMAnN rows",
    "Relative fraction within species profile",
    "Conditional relative fraction among ordinary annotated pathways",
    "Rounded reconstructed pseudo-count"
  ),
  InferentialRole = c("Source", "Source", "Primary", "Primary", "Sensitivity only"),
  AbsoluteAbundanceEvidence = FALSE,
  stringsAsFactors = FALSE
)

log_msg("Inputs and predeclared feature universes validated")

run_maaslin <- function(
  feature_matrix, metadata, formula, interest_metadata, label,
  median_abundance = TRUE, evaluate_only = NULL
) {
  stopifnot(identical(colnames(feature_matrix), rownames(metadata)))
  output <- tempfile(pattern = paste0("maaslin3-", label, "-"))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  set.seed(primary_seed)
  log_msg(
    "MaAsLin3 start: ", label, "; samples=", nrow(metadata),
    "; features=", nrow(feature_matrix), "; evaluate_only=",
    if (is.null(evaluate_only)) "both" else evaluate_only
  )
  fit <- maaslin3::maaslin3(
    input_data = as.data.frame(t(feature_matrix), check.names = FALSE),
    input_metadata = metadata,
    output = output,
    formula = formula,
    min_abundance = 0,
    min_prevalence = 0,
    max_prevalence = 1.01,
    zero_threshold = 0,
    min_variance = 0,
    max_significance = 0.05,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    median_comparison_abundance = median_abundance,
    median_comparison_prevalence = FALSE,
    subtract_median = FALSE,
    warn_prevalence = is.null(evaluate_only),
    augment = TRUE,
    evaluate_only = evaluate_only,
    plot_summary_plot = FALSE,
    plot_associations = FALSE,
    cores = 1,
    save_models = FALSE,
    save_plots_rds = FALSE,
    verbosity = "ERROR"
  )
  collect <- function(component, component_label) {
    if (is.null(component) || is.null(component$results)) return(NULL)
    x <- as.data.frame(component$results, stringsAsFactors = FALSE)
    x <- x[x$metadata == interest_metadata, , drop = FALSE]
    if (nrow(x) == 0L) return(NULL)
    x$Component <- component_label
    x
  }
  result <- rbind(
    collect(fit$fit_data_abundance, "Abundance"),
    collect(fit$fit_data_prevalence, "Prevalence")
  )
  if (is.null(result) || nrow(result) == 0L) {
    stop("No MaAsLin3 rows found for ", interest_metadata, " in ", label, call. = FALSE)
  }
  error_free <- is.na(result$error) | !nzchar(trimws(result$error))
  p_for_fdr <- ifelse(error_free & is.finite(result$pval_individual), result$pval_individual, NA_real_)
  result$QDiseaseBH <- NA_real_
  for (component_label in unique(result$Component)) {
    idx <- result$Component == component_label
    result$QDiseaseBH[idx] <- stats::p.adjust(p_for_fdr[idx], method = "BH")
  }
  joint_p <- tapply(
    result$pval_joint, result$feature,
    function(x) {
      x <- x[is.finite(x)]
      if (length(x)) x[[1L]] else NA_real_
    }
  )
  joint_q <- stats::p.adjust(joint_p, method = "BH")
  result$QJointDiseaseBH <- unname(joint_q[result$feature])
  result$ErrorFree <- error_free
  out <- data.frame(
    FeatureID = result$feature,
    Component = result$Component,
    Contrast = result$name,
    Level = result$value,
    Coefficient = result$coef,
    NullHypothesis = result$null_hypothesis,
    StdError = result$stderr,
    PValue = result$pval_individual,
    QDiseaseBH = result$QDiseaseBH,
    QPackageWide = result$qval_individual,
    PJoint = result$pval_joint,
    QJointDiseaseBH = result$QJointDiseaseBH,
    QJointPackageWide = result$qval_joint,
    ModelError = result$error,
    ErrorFree = result$ErrorFree,
    N = result$N,
    NNonzero = result$N_not_zero,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rm(fit, result)
  invisible(gc())
  log_msg("MaAsLin3 complete: ", label, "; result rows=", nrow(out))
  out
}

decorate_maaslin <- function(result, feature_map, feature_space, abundance_estimable, prevalence_estimable) {
  idx <- match(result$FeatureID, feature_map$FeatureID)
  stopifnot(!anyNA(idx))
  label_column <- if (feature_space == "Species") "Species" else "Pathway"
  result$FeatureSpace <- feature_space
  result$Feature <- feature_map$Feature[idx]
  result$Label <- feature_map[[label_column]][idx]
  result$Estimable <- ifelse(
    result$Component == "Abundance",
    abundance_estimable[result$FeatureID],
    prevalence_estimable[result$FeatureID]
  )
  result$Reportable <- result$ErrorFree & result$Estimable & is.finite(result$PValue)
  result$EffectScale <- ifelse(
    result$Component == "Abundance",
    "log2 relative-abundance coefficient; p-value tests the per-metadatum median null",
    "presence log-odds coefficient; p-value tests zero"
  )
  result$FoldChangeOrOddsRatio <- ifelse(
    result$Component == "Abundance",
    2^result$Coefficient,
    exp(result$Coefficient)
  )
  result$InputScale <- "Relative fraction"
  result[, c(
    "FeatureSpace", "FeatureID", "Feature", "Label", "Component", "Contrast",
    "Level", "Coefficient", "NullHypothesis", "StdError", "FoldChangeOrOddsRatio",
    "PValue", "QDiseaseBH", "QPackageWide", "PJoint", "QJointDiseaseBH",
    "QJointPackageWide", "N", "NNonzero", "Estimable", "ErrorFree",
    "Reportable", "ModelError", "EffectScale", "InputScale"
  )]
}

species_primary_matrix <- species_primary_all[species_keep, , drop = FALSE]
pathway_primary_matrix <- pathway_primary_all[pathway_keep, , drop = FALSE]
pseudocount_primary_matrix <- pseudocount_primary_all[species_keep, , drop = FALSE]

species_abundance_estimable <- setNames(species_abundance_estimable_all[species_keep], rownames(species_primary_matrix))
species_prevalence_estimable <- setNames(species_prevalence_estimable_all[species_keep], rownames(species_primary_matrix))
pathway_abundance_estimable <- setNames(pathway_abundance_estimable_all[pathway_keep], rownames(pathway_primary_matrix))
pathway_prevalence_estimable <- setNames(pathway_prevalence_estimable_all[pathway_keep], rownames(pathway_primary_matrix))

maaslin_formula <- "~ disease + z_age + gender + z_BMI + z_log10_reads"
maaslin_species_raw <- run_maaslin(
  species_primary_matrix, primary_metadata, maaslin_formula,
  "disease", "species-primary", median_abundance = TRUE
)
maaslin_pathway_raw <- run_maaslin(
  pathway_primary_matrix, primary_metadata, maaslin_formula,
  "disease", "pathway-primary", median_abundance = TRUE
)
maaslin_species <- decorate_maaslin(
  maaslin_species_raw, species_input$map, "Species",
  species_abundance_estimable, species_prevalence_estimable
)
maaslin_pathways <- decorate_maaslin(
  maaslin_pathway_raw, pathway_input$map, "Pathway",
  pathway_abundance_estimable, pathway_prevalence_estimable
)

add_check("MaAsLin3", "species-result-rows", nrow(maaslin_species) == 424L, nrow(maaslin_species))
add_check("MaAsLin3", "pathway-result-rows", nrow(maaslin_pathways) == 788L, nrow(maaslin_pathways))
add_check("MaAsLin3", "species-two-components", identical(sort(unique(maaslin_species$Component)), c("Abundance", "Prevalence")), unique(maaslin_species$Component))
add_check("MaAsLin3", "pathway-two-components", identical(sort(unique(maaslin_pathways$Component)), c("Abundance", "Prevalence")), unique(maaslin_pathways$Component))
abundance_null <- maaslin_species$NullHypothesis[maaslin_species$Component == "Abundance"]
abundance_error <- maaslin_species$ModelError[maaslin_species$Component == "Abundance"]
abundance_has_error <- !is.na(abundance_error) & nzchar(trimws(abundance_error))
add_check(
  "MaAsLin3", "abundance-null-median",
  all(abs(abundance_null[is.finite(abundance_null)]) > 0) &&
    all(is.na(abundance_null) == abundance_has_error),
  c(unique(abundance_null[is.finite(abundance_null)]), paste0("model errors=", sum(abundance_has_error)))
)
add_check("MaAsLin3", "prevalence-null-zero", all(maaslin_species$NullHypothesis[maaslin_species$Component == "Prevalence"] == 0), unique(maaslin_species$NullHypothesis[maaslin_species$Component == "Prevalence"]))
add_check("MaAsLin3", "species-disease-bh-range", all(maaslin_species$QDiseaseBH[is.finite(maaslin_species$QDiseaseBH)] >= 0 & maaslin_species$QDiseaseBH[is.finite(maaslin_species$QDiseaseBH)] <= 1), range(maaslin_species$QDiseaseBH, na.rm = TRUE))
add_check("MaAsLin3", "pathway-disease-bh-range", all(maaslin_pathways$QDiseaseBH[is.finite(maaslin_pathways$QDiseaseBH)] >= 0 & maaslin_pathways$QDiseaseBH[is.finite(maaslin_pathways$QDiseaseBH)] <= 1), range(maaslin_pathways$QDiseaseBH, na.rm = TRUE))

# Predeclared MaAsLin3 sensitivity branches --------------------------------
maaslin_median_off <- run_maaslin(
  species_primary_matrix, primary_metadata, maaslin_formula,
  "disease", "species-median-off", median_abundance = FALSE,
  evaluate_only = "abundance"
)

species_keep20 <- species_threshold_prevalence >= 0.20
maaslin_prevalence20 <- run_maaslin(
  species_primary_all[species_keep20, , drop = FALSE],
  primary_metadata, maaslin_formula,
  "disease", "species-prevalence20", median_abundance = TRUE,
  evaluate_only = "abundance"
)

no_bmi_keep <- metadata_all$analysis_group %in% c("Control", "CRC") &
  stats::complete.cases(metadata_all[, c("age", "gender", "number_reads")])
no_bmi_raw <- metadata_all[no_bmi_keep, , drop = FALSE]
no_bmi_metadata <- make_metadata(no_bmi_raw, include_bmi = FALSE)
no_bmi_metadata <- no_bmi_metadata[, c("disease", "z_age", "gender", "z_log10_reads")]
maaslin_no_bmi <- run_maaslin(
  species_input$abundance[species_keep, no_bmi_raw$sample_id, drop = FALSE],
  no_bmi_metadata, "~ disease + z_age + gender + z_log10_reads",
  "disease", "species-no-bmi", median_abundance = TRUE,
  evaluate_only = "abundance"
)

progression_keep <- stats::complete.cases(metadata_all[, c("age", "gender", "BMI", "number_reads")])
progression_raw <- metadata_all[progression_keep, , drop = FALSE]
progression_metadata <- make_metadata(progression_raw, include_bmi = TRUE, progression = TRUE)
progression_metadata <- progression_metadata[, c("progression_score", "z_age", "gender", "z_BMI", "z_log10_reads")]
maaslin_progression <- run_maaslin(
  species_input$abundance[species_keep, progression_raw$sample_id, drop = FALSE],
  progression_metadata,
  "~ progression_score + z_age + gender + z_BMI + z_log10_reads",
  "progression_score", "species-ordered-progression",
  median_abundance = TRUE, evaluate_only = "abundance"
)

add_check("Sensitivity", "median-off-null-zero", all(maaslin_median_off$NullHypothesis == 0), unique(maaslin_median_off$NullHypothesis))
add_check("Sensitivity", "prevalence20-features-reduced", nrow(maaslin_prevalence20) < 212L && nrow(maaslin_prevalence20) > 0L, nrow(maaslin_prevalence20))
add_check("Sensitivity", "no-bmi-samples", nrow(no_bmi_metadata) == 114L, nrow(no_bmi_metadata))
add_check("Sensitivity", "progression-samples", nrow(progression_metadata) == 151L, nrow(progression_metadata))
add_check("Sensitivity", "progression-score-order", identical(sort(unique(progression_metadata$progression_score)), c(0, 1, 2)), unique(progression_metadata$progression_score))

run_ancombc2 <- function(count_matrix, metadata, label) {
  stopifnot(identical(colnames(count_matrix), rownames(metadata)))
  set.seed(primary_seed)
  log_msg("ANCOM-BC2 start: ", label, "; samples=", nrow(metadata), "; features=", nrow(count_matrix))
  fit <- ANCOMBC::ancombc2(
    data = count_matrix,
    taxa_are_rows = TRUE,
    aggregate_data = count_matrix,
    meta_data = metadata,
    fix_formula = "disease + z_age + gender + z_BMI",
    rand_formula = NULL,
    p_adj_method = "BH",
    pseudo = 0,
    pseudo_sens = TRUE,
    prv_cut = 0,
    lib_cut = 0,
    s0_perc = 0.05,
    group = NULL,
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = 1,
    verbose = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE
  )
  res <- fit$res
  out <- data.frame(
    FeatureID = res$taxon,
    Coefficient = res$lfc_diseaseCRC,
    StdError = res$se_diseaseCRC,
    Statistic = res$W_diseaseCRC,
    PValue = res$p_diseaseCRC,
    QPackageBH = res$q_diseaseCRC,
    PassedPseudoSensitivity = res$passed_ss_diseaseCRC,
    SignificantPackage = res$diff_diseaseCRC,
    SignificantRobust = res$diff_robust_diseaseCRC,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$QDiseaseBH <- stats::p.adjust(out$PValue, method = "BH")
  out$EffectScale <- "natural-log ANCOM-BC2 coefficient"
  out$InputScale <- "Reconstructed pseudo-count; not observed taxon reads"
  rm(fit)
  invisible(gc())
  log_msg("ANCOM-BC2 complete: ", label, "; result rows=", nrow(out))
  out
}

ancom_metadata <- primary_metadata[, c("disease", "z_age", "gender", "z_BMI")]
ancombc2_species <- run_ancombc2(pseudocount_primary_matrix, ancom_metadata, "primary-universe")
ancom_idx <- match(ancombc2_species$FeatureID, species_input$map$FeatureID)
ancombc2_species$Feature <- species_input$map$Feature[ancom_idx]
ancombc2_species$Species <- species_input$map$Species[ancom_idx]
ancombc2_species <- ancombc2_species[, c(
  "FeatureID", "Feature", "Species", "Coefficient", "StdError", "Statistic",
  "PValue", "QDiseaseBH", "QPackageBH", "PassedPseudoSensitivity",
  "SignificantPackage", "SignificantRobust", "EffectScale", "InputScale"
)]
add_check("ANCOM-BC2", "result-feature-count", nrow(ancombc2_species) == 212L, nrow(ancombc2_species))
add_check("ANCOM-BC2", "disease-bh-recomputed", max(abs(ancombc2_species$QDiseaseBH - ancombc2_species$QPackageBH), na.rm = TRUE) < 1e-12, max(abs(ancombc2_species$QDiseaseBH - ancombc2_species$QPackageBH), na.rm = TRUE))
add_check("ANCOM-BC2", "pseudo-sensitivity-recorded", all(!is.na(ancombc2_species$PassedPseudoSensitivity)), table(ancombc2_species$PassedPseudoSensitivity))

no_filter_keep <- rowSums(pseudocount_primary_all) > 0
ancombc2_no_filter <- run_ancombc2(
  pseudocount_primary_all[no_filter_keep, , drop = FALSE],
  ancom_metadata, "no-predeclared-feature-filter"
)
add_check("Sensitivity", "ancom-no-filter-feature-count", nrow(ancombc2_no_filter) == 634L, nrow(ancombc2_no_filter))

run_aldex2 <- function(
  count_matrix, metadata, denominator, label,
  denominator_label = if (is.character(denominator)) denominator else paste0("explicit feature vector (n=", length(denominator), ")")
) {
  stopifnot(identical(colnames(count_matrix), rownames(metadata)))
  count_matrix <- round(count_matrix)
  storage.mode(count_matrix) <- "integer"
  model <- stats::model.matrix(
    ~ disease + z_age + gender + z_BMI,
    data = metadata
  )
  stopifnot("diseaseCRC" %in% colnames(model))
  set.seed(primary_seed)
  log_msg(
    "ALDEx2 start: ", label, "; samples=", nrow(metadata),
    "; features=", nrow(count_matrix), "; MC=128; denominator=", denominator_label
  )
  clr <- ALDEx2::aldex.clr(
    count_matrix, model, mc.samples = 128,
    denom = denominator, verbose = FALSE, useMC = FALSE
  )
  glm <- ALDEx2::aldex.glm(clr, verbose = FALSE, fdr.method = "BH")
  out <- data.frame(
    FeatureID = rownames(glm),
    Coefficient = glm[, "diseaseCRC:Est"],
    StdError = glm[, "diseaseCRC:SE"],
    Statistic = glm[, "diseaseCRC:t.val"],
    PValue = glm[, "diseaseCRC:pval"],
    QPackageBH = glm[, "diseaseCRC:pval.padj"],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$QDiseaseBH <- stats::p.adjust(out$PValue, method = "BH")
  out$Denominator <- denominator_label
  out$MonteCarloSamples <- 128L
  out$EffectScale <- "expected CLR coefficient (log2-ratio scale)"
  out$InputScale <- "Reconstructed pseudo-count; not observed taxon reads"
  rm(clr, glm)
  invisible(gc())
  log_msg("ALDEx2 complete: ", label, "; result rows=", nrow(out))
  out
}

aldex_metadata <- ancom_metadata
aldex2_species <- run_aldex2(
  pseudocount_primary_matrix, aldex_metadata, "all", "all-features-denominator"
)
aldex_idx <- match(aldex2_species$FeatureID, species_input$map$FeatureID)
aldex2_species$Feature <- species_input$map$Feature[aldex_idx]
aldex2_species$Species <- species_input$map$Species[aldex_idx]
aldex2_species <- aldex2_species[, c(
  "FeatureID", "Feature", "Species", "Coefficient", "StdError", "Statistic",
  "PValue", "QDiseaseBH", "QPackageBH", "Denominator", "MonteCarloSamples",
  "EffectScale", "InputScale"
)]
add_check("ALDEx2", "result-feature-count", nrow(aldex2_species) == 212L, nrow(aldex2_species))
add_check("ALDEx2", "mc-samples", all(aldex2_species$MonteCarloSamples == 128L), unique(aldex2_species$MonteCarloSamples))
add_check(
  "ALDEx2", "package-q-preserved",
  all(aldex2_species$QPackageBH[is.finite(aldex2_species$QPackageBH)] >= 0 &
        aldex2_species$QPackageBH[is.finite(aldex2_species$QPackageBH)] <= 1) &&
    any(abs(aldex2_species$QDiseaseBH - aldex2_species$QPackageBH) > 1e-12, na.rm = TRUE),
  paste0(
    "maximum difference from disease-only BH=",
    max(abs(aldex2_species$QDiseaseBH - aldex2_species$QPackageBH), na.rm = TRUE)
  )
)

compute_iqlr_denominator <- function(reads, groups) {
  stopifnot(length(groups) == ncol(reads))
  work <- reads
  if (min(work) == 0) work <- work + 0.5
  clr_by_sample <- t(apply(work, 2L, function(x) log2(x) - mean(log2(x))))
  typical_variance <- function(x) {
    variance <- apply(x, 2L, stats::var)
    quartiles <- stats::quantile(variance, probs = c(0.25, 0.75), na.rm = TRUE)
    which(variance > quartiles[[1L]] & variance < quartiles[[2L]])
  }
  sets <- lapply(levels(groups), function(level) {
    typical_variance(clr_by_sample[groups == level, , drop = FALSE])
  })
  sets[[length(sets) + 1L]] <- typical_variance(clr_by_sample)
  Reduce(intersect, sets)
}

iqlr_denominator <- compute_iqlr_denominator(
  pseudocount_primary_matrix, aldex_metadata$disease
)
add_check(
  "Sensitivity", "aldex-iqlr-denominator-size",
  length(iqlr_denominator) > 5L && length(iqlr_denominator) < nrow(pseudocount_primary_matrix),
  length(iqlr_denominator)
)
aldex2_iqlr <- run_aldex2(
  pseudocount_primary_matrix, aldex_metadata, iqlr_denominator,
  "iqlr-explicit-denominator",
  paste0("IQLR explicit feature vector (n=", length(iqlr_denominator), ")")
)

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) return(NA_real_)
  unname(stats::cor(x[keep], y[keep], method = method))
}
direction_agreement <- function(x, y) {
  keep <- is.finite(x) & is.finite(y) & x != 0 & y != 0
  if (!any(keep)) return(NA_real_)
  mean(sign(x[keep]) == sign(y[keep]))
}

primary_abundance <- maaslin_species[maaslin_species$Component == "Abundance", , drop = FALSE]
match_primary <- function(x) match(x$FeatureID, primary_abundance$FeatureID)
make_maaslin_sensitivity <- function(method_result, variant, units, formula, primary = FALSE, notes) {
  idx <- match_primary(method_result)
  data.frame(
    Method = "MaAsLin3",
    Variant = variant,
    Units = units,
    Features = nrow(method_result),
    Formula = formula,
    Seed = primary_seed,
    Primary = primary,
    SignificantFDR05 = sum(method_result$QDiseaseBH < 0.05, na.rm = TRUE),
    SharedWithPrimary = sum(!is.na(idx)),
    EffectSpearmanWithPrimary = safe_cor(method_result$Coefficient, primary_abundance$Coefficient[idx]),
    DirectionAgreementWithPrimary = direction_agreement(method_result$Coefficient, primary_abundance$Coefficient[idx]),
    EffectScale = "log2 relative-abundance coefficient",
    InputScale = "Relative fraction",
    Notes = notes,
    stringsAsFactors = FALSE
  )
}

sensitivity_audit <- rbind(
  make_maaslin_sensitivity(
    primary_abundance, "Primary median comparison", 110L, maaslin_formula,
    TRUE, "Per-metadatum median null"
  ),
  make_maaslin_sensitivity(
    maaslin_median_off, "Median comparison off", 110L, maaslin_formula,
    FALSE, "Zero null on the relative scale"
  ),
  make_maaslin_sensitivity(
    maaslin_prevalence20, "Species threshold prevalence >=20%", 110L,
    maaslin_formula, FALSE, "Feature-universe sensitivity"
  ),
  make_maaslin_sensitivity(
    maaslin_no_bmi, "Remove BMI; complete cases", 114L,
    "~ disease + z_age + gender + z_log10_reads", FALSE,
    "Covariate and sample-set sensitivity"
  ),
  make_maaslin_sensitivity(
    maaslin_progression, "Ordered Control-Adenoma-CRC score", 151L,
    "~ progression_score + z_age + gender + z_BMI + z_log10_reads", FALSE,
    "Linear ordinal-trend sensitivity; score fixed to 0/1/2"
  )
)

ancom_primary_idx <- match(ancombc2_no_filter$FeatureID, ancombc2_species$FeatureID)
sensitivity_audit <- rbind(
  sensitivity_audit,
  data.frame(
    Method = "ANCOM-BC2",
    Variant = c("Primary predeclared universe", "No predeclared feature filter"),
    Units = 110L,
    Features = c(nrow(ancombc2_species), nrow(ancombc2_no_filter)),
    Formula = "disease + z_age + gender + z_BMI",
    Seed = primary_seed,
    Primary = c(FALSE, FALSE),
    SignificantFDR05 = c(
      sum(ancombc2_species$SignificantRobust, na.rm = TRUE),
      sum(ancombc2_no_filter$SignificantRobust, na.rm = TRUE)
    ),
    SharedWithPrimary = c(nrow(ancombc2_species), sum(!is.na(ancom_primary_idx))),
    EffectSpearmanWithPrimary = c(
      1,
      safe_cor(
        ancombc2_no_filter$Coefficient,
        ancombc2_species$Coefficient[ancom_primary_idx]
      )
    ),
    DirectionAgreementWithPrimary = c(
      1,
      direction_agreement(
        ancombc2_no_filter$Coefficient,
        ancombc2_species$Coefficient[ancom_primary_idx]
      )
    ),
    EffectScale = "natural-log ANCOM-BC2 coefficient",
    InputScale = "Reconstructed pseudo-count",
    Notes = c(
      "pseudo_sens=TRUE; only robust hits are counted",
      "634 non-all-zero species; method sensitivity only"
    ),
    stringsAsFactors = FALSE
  )
)

iqlr_idx <- match(aldex2_iqlr$FeatureID, aldex2_species$FeatureID)
sensitivity_audit <- rbind(
  sensitivity_audit,
  data.frame(
    Method = "ALDEx2",
    Variant = c("All-features denominator", "IQLR denominator"),
    Units = 110L,
    Features = 212L,
    Formula = "disease + z_age + gender + z_BMI",
    Seed = primary_seed,
    Primary = c(FALSE, FALSE),
    SignificantFDR05 = c(
      sum(aldex2_species$QDiseaseBH < 0.05, na.rm = TRUE),
      sum(aldex2_iqlr$QDiseaseBH < 0.05, na.rm = TRUE)
    ),
    SharedWithPrimary = 212L,
    EffectSpearmanWithPrimary = c(
      1,
      safe_cor(aldex2_iqlr$Coefficient, aldex2_species$Coefficient[iqlr_idx])
    ),
    DirectionAgreementWithPrimary = c(
      1,
      direction_agreement(aldex2_iqlr$Coefficient, aldex2_species$Coefficient[iqlr_idx])
    ),
    EffectScale = "expected CLR coefficient (log2-ratio scale)",
    InputScale = "Reconstructed pseudo-count",
    Notes = c(
      "128 Monte Carlo instances; all-features denominator",
      paste0(
        "128 Monte Carlo instances; explicit IQLR denominator with ",
        length(iqlr_denominator),
        " features because ALDEx2 1.42.0 cannot derive iqlr directly from a model matrix"
      )
    ),
    stringsAsFactors = FALSE
  )
)

add_check("Sensitivity", "branch-count", nrow(sensitivity_audit) == 9L, nrow(sensitivity_audit))
add_check("Sensitivity", "all-runs-fixed-seed", all(sensitivity_audit$Seed == primary_seed), unique(sensitivity_audit$Seed))
add_check("Sensitivity", "single-maaslin-primary", sum(sensitivity_audit$Primary) == 1L, sum(sensitivity_audit$Primary))

# Cross-method evidence keeps incomparable scales separate -----------------
cross_method_evidence <- data.frame(
  FeatureID = primary_abundance$FeatureID,
  Feature = primary_abundance$Feature,
  Species = primary_abundance$Label,
  MaAsLin3CoefficientLog2Relative = primary_abundance$Coefficient,
  MaAsLin3NullMedian = primary_abundance$NullHypothesis,
  MaAsLin3QDiseaseBH = primary_abundance$QDiseaseBH,
  MaAsLin3Estimable = primary_abundance$Estimable,
  MaAsLin3ErrorFree = primary_abundance$ErrorFree,
  stringsAsFactors = FALSE
)
an_idx <- match(cross_method_evidence$FeatureID, ancombc2_species$FeatureID)
al_idx <- match(cross_method_evidence$FeatureID, aldex2_species$FeatureID)
stopifnot(!anyNA(an_idx), !anyNA(al_idx))
cross_method_evidence$ANCOMBC2CoefficientNaturalLog <- ancombc2_species$Coefficient[an_idx]
cross_method_evidence$ANCOMBC2QDiseaseBH <- ancombc2_species$QDiseaseBH[an_idx]
cross_method_evidence$ANCOMBC2PassedPseudoSensitivity <- ancombc2_species$PassedPseudoSensitivity[an_idx]
cross_method_evidence$ANCOMBC2RobustHit <- ancombc2_species$SignificantRobust[an_idx]
cross_method_evidence$ALDEx2CoefficientCLRLog2 <- aldex2_species$Coefficient[al_idx]
cross_method_evidence$ALDEx2QDiseaseBH <- aldex2_species$QDiseaseBH[al_idx]
cross_method_evidence$MaAsLin3Direction <- sign(cross_method_evidence$MaAsLin3CoefficientLog2Relative)
cross_method_evidence$ANCOMBC2Direction <- sign(cross_method_evidence$ANCOMBC2CoefficientNaturalLog)
cross_method_evidence$ALDEx2Direction <- sign(cross_method_evidence$ALDEx2CoefficientCLRLog2)
cross_method_evidence$AllThreeDirectionConcordant <- with(
  cross_method_evidence,
  MaAsLin3Direction == ANCOMBC2Direction & ANCOMBC2Direction == ALDEx2Direction
)
cross_method_evidence$MaAsLin3Hit <- with(
  cross_method_evidence,
  MaAsLin3Estimable & MaAsLin3ErrorFree & MaAsLin3QDiseaseBH < 0.05
)
cross_method_evidence$ALDEx2Hit <- cross_method_evidence$ALDEx2QDiseaseBH < 0.05
cross_method_evidence$SignificantMethodCount <- rowSums(cbind(
  cross_method_evidence$MaAsLin3Hit,
  cross_method_evidence$ANCOMBC2RobustHit,
  cross_method_evidence$ALDEx2Hit
), na.rm = TRUE)
cross_method_evidence$EvidenceLabel <- ifelse(
  cross_method_evidence$AllThreeDirectionConcordant,
  "Direction concordant across three methods",
  "Direction differs across methods"
)
cross_method_evidence$IndependentReplication <- FALSE
cross_method_evidence$InputBoundary <- "ANCOM-BC2 and ALDEx2 use reconstructed pseudo-counts"

add_check("Cross-method", "feature-count", nrow(cross_method_evidence) == 212L, nrow(cross_method_evidence))
add_check("Cross-method", "no-independent-replication-claim", all(!cross_method_evidence$IndependentReplication), unique(cross_method_evidence$IndependentReplication))
add_check("Cross-method", "three-effect-scales-retained", all(c(
  "MaAsLin3CoefficientLog2Relative", "ANCOMBC2CoefficientNaturalLog",
  "ALDEx2CoefficientCLRLog2"
) %in% names(cross_method_evidence)), "Three non-comparable coefficients")

log_msg("Primary and sensitivity models complete")

# Original-paper marker audit ----------------------------------------------
marker_species <- c(
  "Fusobacterium nucleatum",
  "Peptostreptococcus stomatis",
  "Porphyromonas asaccharolytica",
  "Eubacterium rectale"
)
marker_map_idx <- match(marker_species, species_input$map$Species)
add_check("Paper anchor", "four-marker-identities", !anyNA(marker_map_idx), marker_species)
marker_ids <- species_input$map$FeatureID[marker_map_idx]
marker_long <- do.call(rbind, lapply(seq_along(marker_species), function(i) {
  data.frame(
    FeatureID = marker_ids[[i]],
    Species = marker_species[[i]],
    SampleID = sample_ids,
    Group = factor(metadata_all$analysis_group, levels = c("Control", "Adenoma", "CRC")),
    RelativeAbundance = species_input$abundance[marker_ids[[i]], ],
    RelativeAbundancePercent = 100 * species_input$abundance[marker_ids[[i]], ],
    Present = species_input$abundance[marker_ids[[i]], ] > 0,
    stringsAsFactors = FALSE
  )
}))
marker_long$Species <- factor(marker_long$Species, levels = marker_species)
marker_long$Group <- factor(marker_long$Group, levels = c("Control", "Adenoma", "CRC"))

original_marker_audit <- do.call(rbind, lapply(marker_species, function(species_name) {
  do.call(rbind, lapply(c("Control", "Adenoma", "CRC"), function(group_name) {
    x <- marker_long[
      as.character(marker_long$Species) == species_name &
        as.character(marker_long$Group) == group_name,
      , drop = FALSE
    ]
    data.frame(
      Species = species_name,
      Group = group_name,
      Subjects = nrow(x),
      Present = sum(x$Present),
      Absent = sum(!x$Present),
      Prevalence = mean(x$Present),
      MedianRelativeAbundancePercent = stats::median(x$RelativeAbundancePercent),
      Q1RelativeAbundancePercent = unname(stats::quantile(x$RelativeAbundancePercent, 0.25)),
      Q3RelativeAbundancePercent = unname(stats::quantile(x$RelativeAbundancePercent, 0.75)),
      MeanRelativeAbundancePercent = mean(x$RelativeAbundancePercent),
      PaperAnchor = "Zeller et al. 2014 Figure 1A/E species context; harmonized species resolution",
      stringsAsFactors = FALSE
    )
  }))
}))
original_marker_audit$Species <- factor(
  original_marker_audit$Species,
  levels = marker_species
)
original_marker_audit$Group <- factor(
  original_marker_audit$Group,
  levels = c("Control", "Adenoma", "CRC")
)
marker_primary_idx <- match(original_marker_audit$Species, species_input$map$Species)
marker_cross_idx <- match(original_marker_audit$Species, cross_method_evidence$Species)
original_marker_audit$PrimarySpeciesUniverse <- species_keep[marker_primary_idx]
original_marker_audit$PrimaryPrevalenceEstimable <- species_prevalence_estimable_all[marker_primary_idx]
original_marker_audit$MaAsLin3AbundanceQ <- cross_method_evidence$MaAsLin3QDiseaseBH[marker_cross_idx]
original_marker_audit$ANCOMBC2RobustHit <- cross_method_evidence$ANCOMBC2RobustHit[marker_cross_idx]
original_marker_audit$ALDEx2Q <- cross_method_evidence$ALDEx2QDiseaseBH[marker_cross_idx]
original_marker_audit$CurrentAnalysisRole <- "Descriptive paper anchor; F. nucleatum subspecies are collapsed in the harmonized profile"

add_check("Paper anchor", "marker-audit-rows", nrow(original_marker_audit) == 12L, nrow(original_marker_audit))
add_check(
  "Paper anchor", "marker-group-counts",
  identical(
    as.integer(tapply(original_marker_audit$Subjects, original_marker_audit$Group, unique)[c("Control", "Adenoma", "CRC")]),
    c(61L, 42L, 53L)
  ),
  paste(tapply(original_marker_audit$Subjects, original_marker_audit$Group, unique), collapse = "/")
)

interpretation_boundaries <- data.frame(
  Topic = c(
    "Inference unit", "Primary contrast", "MaAsLin3 abundance",
    "MaAsLin3 prevalence", "Joint association", "Species FDR",
    "Pathway denominator", "ANCOM-BC2 input", "ALDEx2 input",
    "Cross-method overlap", "Original paper anchor", "Causality"
  ),
  AuthorizedInterpretation = c(
    "One independent subject",
    "CRC versus Control among 110 complete cases",
    "A log2 relative-abundance coefficient tested against the disease-coefficient median",
    "A presence log-odds coefficient reported only when each group has at least 10 present and 10 absent subjects",
    "Evidence that at least one estimable component departs from its declared null",
    "BH correction within the 212-species disease contrast and within each model component",
    "Composition conditional on 493 annotated ordinary unstratified pathways",
    "Pseudo-count sensitivity using round(relative fraction x whole-metagenome reads)",
    "CLR sensitivity using the same reconstructed pseudo-count convention",
    "Direction concordance among methods with distinct assumptions and effect scales",
    "Conceptual redraw from harmonized MetaPhlAn3 profiles; not a point-for-point reconstruction of the original pipeline",
    "Adjusted cross-sectional association"
  ),
  ProhibitedClaim = c(
    "Repeated samples or technical replicates as independent subjects",
    "Adenoma silently treated as Control",
    "An unqualified absolute-abundance fold change",
    "A separated or sparse logistic coefficient presented as stable evidence",
    "Both abundance and prevalence must be significant",
    "One FDR family combining species, pathways, abundance and prevalence",
    "The fraction of all microbial metabolic potential",
    "Observed taxon reads or absolute microbial load",
    "Raw counts measured directly by the sequencer for each taxon",
    "Independent replication, majority vote or ground truth",
    "Exact numerical replication of Zeller et al. KEGG/CAZy results",
    "CRC causes the microbial association or vice versa"
  ),
  stringsAsFactors = FALSE
)

# Reproducibility records --------------------------------------------------
core_versions <- c(
  R = paste(R.version$major, R.version$minor, sep = "."),
  maaslin3 = as.character(utils::packageVersion("maaslin3")),
  ANCOMBC = as.character(utils::packageVersion("ANCOMBC")),
  ALDEx2 = as.character(utils::packageVersion("ALDEx2")),
  ggplot2 = as.character(utils::packageVersion("ggplot2")),
  patchwork = as.character(utils::packageVersion("patchwork")),
  digest = as.character(utils::packageVersion("digest")),
  jsonlite = as.character(utils::packageVersion("jsonlite"))
)
required_versions <- c(
  R = "4.5.2", maaslin3 = "1.5.3", ANCOMBC = "2.12.0", ALDEx2 = "1.42.0",
  ggplot2 = NA_character_, patchwork = NA_character_, digest = NA_character_,
  jsonlite = NA_character_
)
maaslin_archive <- if ("maaslin-archive" %in% names(args)) {
  normalizePath(args[["maaslin-archive"]], mustWork = TRUE)
} else {
  file.path(
    project_root, "data", "raw", "article24",
    "maaslin3-3a194ece449ef249354df394b58bfe3e6f951ca3.tar.gz"
  )
}
maaslin_archive_sha <- if (file.exists(maaslin_archive)) sha256_file(maaslin_archive) else NA_character_
package_versions <- data.frame(
  Component = names(core_versions),
  Version = unname(core_versions),
  RequiredVersion = unname(required_versions[names(core_versions)]),
  SourceLock = c(
    "Conda explicit lock",
    "Git commit 3a194ece449ef249354df394b58bfe3e6f951ca3",
    "Bioconductor 3.22 / conda package",
    "Bioconductor 3.22 source release",
    "Conda explicit lock", "Conda plus R package lock",
    "Conda explicit lock", "Conda explicit lock"
  ),
  SourceArchiveSHA256 = c(
    NA_character_, maaslin_archive_sha, rep(NA_character_, 6L)
  ),
  stringsAsFactors = FALSE
)
add_check("Environment", "r-version", identical(core_versions[["R"]], "4.5.2"), core_versions[["R"]])
add_check("Environment", "maaslin3-version", identical(core_versions[["maaslin3"]], "1.5.3"), core_versions[["maaslin3"]])
add_check("Environment", "ancombc-version", identical(core_versions[["ANCOMBC"]], "2.12.0"), core_versions[["ANCOMBC"]])
add_check("Environment", "aldex2-version", identical(core_versions[["ALDEx2"]], "1.42.0"), core_versions[["ALDEx2"]])
add_check(
  "Environment", "maaslin3-source-sha256",
  identical(maaslin_archive_sha, "358f35e094e03026c8d694d731c0d1e1e6c9a15dec3c466649b9db9329ca1f07"),
  maaslin_archive_sha
)

# Publication graphics ----------------------------------------------------
pal_group <- c(Control = "#0072B2", Adenoma = "#E69F00", CRC = "#D55E00")
pal_direction <- c(
  Depleted = "#0072B2",
  Neutral = "#BDBDBD",
  Enriched = "#D55E00",
  `Not estimable` = "#F2F2F2"
)
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
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  tiff_path <- file.path(figure_dir, paste0(stem, ".tiff"))
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = 350, bg = "white")
  ggplot2::ggsave(tiff_path, plot = plot, width = width, height = height, units = "in", dpi = 350, compression = "lzw", bg = "white")
  invisible(c(pdf_path, png_path, tiff_path))
}

marker_abundance_plot <- ggplot2::ggplot(
  marker_long,
  ggplot2::aes(x = Group, y = RelativeAbundancePercent, colour = Group, fill = Group)
) +
  ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.22, linewidth = 0.45) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.14, height = 0, seed = primary_seed),
    size = 1.15, alpha = 0.62
  ) +
  ggplot2::facet_wrap(~ Species, scales = "free_y", ncol = 2) +
  ggplot2::scale_colour_manual(values = pal_group, drop = FALSE) +
  ggplot2::scale_fill_manual(values = pal_group, drop = FALSE) +
  ggplot2::scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.0001),
    labels = scales::label_number(accuracy = 0.001),
    n.breaks = 4
  ) +
  ggplot2::labs(
    title = "A  Relative abundance",
    x = NULL, y = "Relative abundance (%) — pseudo-log10",
    colour = NULL, fill = NULL
  ) +
  theme_pub(10) +
  ggplot2::theme(legend.position = "none")

marker_prevalence_plot <- ggplot2::ggplot(
  original_marker_audit,
  ggplot2::aes(x = Group, y = Prevalence, fill = Group)
) +
  ggplot2::geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(Prevalence, accuracy = 1)),
    vjust = -0.25, size = 2.7, colour = "#222222"
  ) +
  ggplot2::facet_wrap(~ Species, ncol = 2) +
  ggplot2::scale_fill_manual(values = pal_group, drop = FALSE) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.08), breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1), expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "B  Detection prevalence",
    x = NULL, y = "Subjects with detected species", fill = NULL
  ) +
  theme_pub(10) +
  ggplot2::theme(legend.position = "none")

marker_figure <- marker_abundance_plot / marker_prevalence_plot +
  patchwork::plot_annotation(
    title = "Selected Zeller Figure 1 species in the harmonized stool cohort",
    subtitle = "Descriptive profiles use all 156 independent subjects; the adjusted primary contrast uses 110 complete cases",
    caption = "MetaPhlAn3 collapses the two original F. nucleatum subspecies. E. rectale represents a depleted Figure 1A signature; zeros are retained."
  )
save_pub(marker_figure, "24-zeller-marker-redraw", 11.5, 10.5)

select_top_results <- function(x, n = 16L) {
  candidate <- x[x$Reportable & is.finite(x$Coefficient) & is.finite(x$StdError), , drop = FALSE]
  candidate$SortQ <- ifelse(is.finite(candidate$QDiseaseBH), candidate$QDiseaseBH, 1)
  candidate <- candidate[order(candidate$SortQ, -abs(candidate$Coefficient - candidate$NullHypothesis)), , drop = FALSE]
  utils::head(candidate, n)
}

forest_component <- function(x, component_label, x_label) {
  x <- select_top_results(x[x$Component == component_label, , drop = FALSE], 16L)
  x$Display <- factor(x$Label, levels = rev(x$Label))
  x$Lower <- x$Coefficient - 1.96 * x$StdError
  x$Upper <- x$Coefficient + 1.96 * x$StdError
  x$FDR <- ifelse(x$QDiseaseBH < 0.05, "FDR < 0.05", "FDR ≥ 0.05")
  null_value <- stats::median(x$NullHypothesis, na.rm = TRUE)
  ggplot2::ggplot(x, ggplot2::aes(y = Display, x = Coefficient, colour = FDR)) +
    ggplot2::geom_vline(xintercept = null_value, linetype = 2, colour = "#666666", linewidth = 0.45) +
    ggplot2::geom_segment(ggplot2::aes(x = Lower, xend = Upper, yend = Display), linewidth = 0.55) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_colour_manual(
      values = c("FDR < 0.05" = "#D55E00", "FDR ≥ 0.05" = "#777777"),
      limits = c("FDR < 0.05", "FDR ≥ 0.05"),
      drop = FALSE
    ) +
    ggplot2::labs(title = component_label, x = x_label, y = NULL, colour = NULL) +
    theme_pub(9.3)
}

abundance_component_figure <- forest_component(
  maaslin_species, "Abundance", "CRC coefficient (log2 relative abundance)"
)
prevalence_component_figure <- forest_component(
  maaslin_species, "Prevalence", "CRC coefficient (presence log odds)"
) + ggplot2::guides(colour = "none")
two_part_figure <- abundance_component_figure | prevalence_component_figure
two_part_figure <- two_part_figure +
  patchwork::plot_annotation(
  title = "MaAsLin3 separates abundance from prevalence evidence",
  subtitle = "Disease-specific BH correction within each 212-species component",
  caption = "Dashed lines show the declared null: the disease-coefficient median for abundance and zero for prevalence. Only estimable rows are displayed."
  ) & ggplot2::theme(legend.position = "bottom")
save_pub(two_part_figure, "24-maaslin3-two-part", 13.2, 8.2)

pathway_abundance_results <- maaslin_pathways[maaslin_pathways$Component == "Abundance", , drop = FALSE]
top_pathways <- select_top_results(pathway_abundance_results, 20L)
top_pathways$Display <- factor(top_pathways$Label, levels = rev(top_pathways$Label))
top_pathways$Lower <- top_pathways$Coefficient - 1.96 * top_pathways$StdError
top_pathways$Upper <- top_pathways$Coefficient + 1.96 * top_pathways$StdError
top_pathways$FDR <- ifelse(top_pathways$QDiseaseBH < 0.05, "FDR < 0.05", "FDR ≥ 0.05")
pathway_null <- stats::median(top_pathways$NullHypothesis, na.rm = TRUE)
functional_figure <- ggplot2::ggplot(
  top_pathways,
  ggplot2::aes(y = Display, x = Coefficient, colour = FDR)
) +
  ggplot2::geom_vline(xintercept = pathway_null, linetype = 2, colour = "#666666", linewidth = 0.45) +
  ggplot2::geom_segment(ggplot2::aes(x = Lower, xend = Upper, yend = Display), linewidth = 0.58) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_colour_manual(values = c("FDR < 0.05" = "#D55E00", "FDR ≥ 0.05" = "#777777")) +
  ggplot2::labs(
    title = "Top MetaCyc pathway abundance associations",
    subtitle = "CRC versus Control, adjusted for age, gender, BMI and log10 read depth",
    x = "CRC coefficient (log2 relative pathway abundance)", y = NULL, colour = NULL,
    caption = "Pathways are reclosed within 493 ordinary unstratified annotated rows. The dashed line is the pathway-coefficient median null."
  ) +
  theme_pub(9.6)
save_pub(functional_figure, "24-functional-associations", 11.3, 8.8)

cross_method_evidence$MinimumQ <- pmin(
  cross_method_evidence$MaAsLin3QDiseaseBH,
  cross_method_evidence$ANCOMBC2QDiseaseBH,
  cross_method_evidence$ALDEx2QDiseaseBH,
  na.rm = TRUE
)
cross_method_evidence$AnyMethodHit <- with(
  cross_method_evidence,
  MaAsLin3Hit | ANCOMBC2RobustHit | ALDEx2Hit
)
concordance_selected <- cross_method_evidence[
  order(-cross_method_evidence$SignificantMethodCount, cross_method_evidence$MinimumQ),
  , drop = FALSE
]
concordance_selected <- utils::head(concordance_selected, 24L)
concordance_long <- rbind(
  data.frame(
    Species = concordance_selected$Species, Method = "MaAsLin3",
    Coefficient = concordance_selected$MaAsLin3CoefficientLog2Relative,
    QValue = concordance_selected$MaAsLin3QDiseaseBH,
    Hit = concordance_selected$MaAsLin3Hit,
    EffectScale = "log2 relative coefficient", stringsAsFactors = FALSE
  ),
  data.frame(
    Species = concordance_selected$Species, Method = "ANCOM-BC2",
    Coefficient = concordance_selected$ANCOMBC2CoefficientNaturalLog,
    QValue = concordance_selected$ANCOMBC2QDiseaseBH,
    Hit = concordance_selected$ANCOMBC2RobustHit,
    EffectScale = "natural-log coefficient", stringsAsFactors = FALSE
  ),
  data.frame(
    Species = concordance_selected$Species, Method = "ALDEx2",
    Coefficient = concordance_selected$ALDEx2CoefficientCLRLog2,
    QValue = concordance_selected$ALDEx2QDiseaseBH,
    Hit = concordance_selected$ALDEx2Hit,
    EffectScale = "CLR log2 coefficient", stringsAsFactors = FALSE
  )
)
concordance_long$Direction <- factor(
  ifelse(
    !is.finite(concordance_long$Coefficient),
    "Not estimable",
    ifelse(concordance_long$Coefficient > 0, "Enriched", ifelse(concordance_long$Coefficient < 0, "Depleted", "Neutral"))
  ),
  levels = c("Depleted", "Neutral", "Enriched", "Not estimable")
)
concordance_long$Evidence <- -log10(pmax(concordance_long$QValue, 1e-12))
concordance_long$Evidence[!is.finite(concordance_long$Evidence)] <- 0
concordance_long$Method <- factor(concordance_long$Method, levels = c("MaAsLin3", "ANCOM-BC2", "ALDEx2"))
concordance_long$Species <- factor(
  concordance_long$Species,
  levels = rev(concordance_selected$Species)
)
concordance_figure <- ggplot2::ggplot(
  concordance_long,
  ggplot2::aes(x = Method, y = Species, fill = Direction, size = Evidence, shape = Hit)
) +
  ggplot2::geom_point(colour = "#333333", stroke = 0.45, alpha = 0.92) +
  ggplot2::scale_fill_manual(values = pal_direction, drop = TRUE) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 22), labels = c(`FALSE` = "No", `TRUE` = "Yes")) +
  ggplot2::scale_size_continuous(range = c(1.8, 6.2), name = expression(-log[10](q))) +
  ggplot2::labs(
    title = "Cross-method direction concordance",
    subtitle = "Point size represents within-method disease FDR; square points meet each method's reporting rule",
    x = NULL, y = NULL, fill = "Direction", shape = "Reported hit",
    caption = paste0(
      "Coefficient magnitudes are not compared across methods.\n",
      "Pale grey denotes a coefficient that was not estimable.\n",
      "ANCOM-BC2 and ALDEx2 use reconstructed pseudo-counts, so concordance is not independent replication."
    )
  ) +
  ggplot2::guides(
    size = ggplot2::guide_legend(order = 1, nrow = 1),
    shape = ggplot2::guide_legend(order = 2, nrow = 1),
    fill = ggplot2::guide_legend(order = 3, nrow = 1)
  ) +
  theme_pub(9.6) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(colour = "#E6E6E6"),
    legend.box = "vertical"
  )
save_pub(concordance_figure, "24-method-concordance", 10.5, 9.5)

# Write result tables ------------------------------------------------------
maaslin_species$ModelError <- gsub("[\r\n\t]+", " ", maaslin_species$ModelError)
maaslin_pathways$ModelError <- gsub("[\r\n\t]+", " ", maaslin_pathways$ModelError)
write_tsv(data_lineage, file.path(output_dir, "data-lineage.tsv"))
write_tsv(design_balance, file.path(output_dir, "design-balance.tsv"))
write_tsv(feature_filter_audit, file.path(output_dir, "feature-filter-audit.tsv"))
write_tsv(original_marker_audit, file.path(output_dir, "original-marker-audit.tsv"))
write_tsv_gz(maaslin_species, file.path(output_dir, "maaslin3-species.tsv.gz"))
write_tsv_gz(maaslin_pathways, file.path(output_dir, "maaslin3-pathways.tsv.gz"))
write_tsv_gz(ancombc2_species, file.path(output_dir, "ancombc2-species.tsv.gz"))
write_tsv_gz(aldex2_species, file.path(output_dir, "aldex2-species.tsv.gz"))
write_tsv(cross_method_evidence, file.path(output_dir, "cross-method-evidence.tsv"))
write_tsv(sensitivity_audit, file.path(output_dir, "sensitivity-audit.tsv"))
write_tsv(interpretation_boundaries, file.path(output_dir, "interpretation-boundaries.tsv"))
write_tsv(package_versions, file.path(output_dir, "package-versions.tsv"))

expected_result_files <- c(
  "data-lineage.tsv", "design-balance.tsv", "feature-filter-audit.tsv",
  "original-marker-audit.tsv", "maaslin3-species.tsv.gz", "maaslin3-pathways.tsv.gz",
  "ancombc2-species.tsv.gz", "aldex2-species.tsv.gz", "cross-method-evidence.tsv",
  "sensitivity-audit.tsv", "interpretation-boundaries.tsv", "package-versions.tsv",
  "validation.log"
)
for (file_name in expected_result_files) {
  path <- file.path(output_dir, file_name)
  add_check("Outputs", paste0("result-", file_name), file.exists(path) && file.info(path)$size > 0, if (file.exists(path)) file.info(path)$size else "missing")
}
for (stem in c(
  "24-zeller-marker-redraw", "24-maaslin3-two-part",
  "24-functional-associations", "24-method-concordance"
)) {
  for (extension in c("pdf", "png", "tiff")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    add_check("Figures", paste0(stem, "-", extension), file.exists(path) && file.info(path)$size > 1000, if (file.exists(path)) file.info(path)$size else "missing")
  }
}
paper_figure_path <- file.path(project_root, "figures", "24-zeller-fig1-original.png")
add_check(
  "Paper anchor", "original-figure-present",
  file.exists(paper_figure_path) && file.info(paper_figure_path)$size > 10000,
  if (file.exists(paper_figure_path)) file.info(paper_figure_path)$size else "missing"
)

write_tsv(checks, file.path(output_dir, "validation-audit.tsv"))
checks_failed <- sum(checks$Status == "FAIL")
summary <- list(
  status = if (checks_failed == 0L) "passed" else "failed",
  seed = primary_seed,
  primary_samples = length(primary_ids),
  control_samples = sum(primary_metadata$disease == "Control"),
  crc_samples = sum(primary_metadata$disease == "CRC"),
  independent_subjects = length(unique(primary_metadata_raw$subject_id)),
  species_features = nrow(species_primary_matrix),
  pathway_features = nrow(pathway_primary_matrix),
  species_prevalence_estimable = sum(species_prevalence_estimable),
  species_abundance_estimable = sum(species_abundance_estimable),
  maaslin3_species_abundance_hits = sum(
    maaslin_species$Component == "Abundance" & maaslin_species$Reportable &
      maaslin_species$QDiseaseBH < 0.05,
    na.rm = TRUE
  ),
  maaslin3_species_prevalence_hits = sum(
    maaslin_species$Component == "Prevalence" & maaslin_species$Reportable &
      maaslin_species$QDiseaseBH < 0.05,
    na.rm = TRUE
  ),
  maaslin3_pathway_abundance_hits = sum(
    maaslin_pathways$Component == "Abundance" & maaslin_pathways$Reportable &
      maaslin_pathways$QDiseaseBH < 0.05,
    na.rm = TRUE
  ),
  ancombc2_robust_hits = sum(ancombc2_species$SignificantRobust, na.rm = TRUE),
  aldex2_hits = sum(aldex2_species$QDiseaseBH < 0.05, na.rm = TRUE),
  all_three_direction_concordant = sum(cross_method_evidence$AllThreeDirectionConcordant, na.rm = TRUE),
  all_three_reported_hits = sum(cross_method_evidence$SignificantMethodCount == 3L, na.rm = TRUE),
  sensitivity_branches = nrow(sensitivity_audit),
  checksum_entries = checksum_entries,
  checks_total = nrow(checks),
  checks_passed = sum(checks$Status == "PASS"),
  checks_failed = checks_failed,
  r_version = core_versions[["R"]],
  maaslin3_version = core_versions[["maaslin3"]],
  ancombc_version = core_versions[["ANCOMBC"]],
  aldex2_version = core_versions[["ALDEx2"]],
  maaslin3_commit = "3a194ece449ef249354df394b58bfe3e6f951ca3",
  maaslin3_source_sha256 = maaslin_archive_sha
)
jsonlite::write_json(
  summary,
  file.path(output_dir, "validation-summary.json"),
  pretty = TRUE, auto_unbox = TRUE, na = "null"
)
log_msg(
  "Validation complete: status=", summary$status,
  "; checks=", summary$checks_passed, "/", summary$checks_total
)

if (checks_failed > 0L) {
  failed <- checks[checks$Status == "FAIL", , drop = FALSE]
  stop(
    "Article 24 validation failed: ",
    paste(paste0(failed$CheckID, " [", failed$Detail, "]"), collapse = "; "),
    call. = FALSE
  )
}
