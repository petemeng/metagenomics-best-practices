#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260723L
primary_nperm <- 9999L
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
  "vegan", "ape", "ggplot2", "patchwork", "scales",
  "jsonlite", "digest", "magick"
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
close_rows <- function(x) {
  totals <- rowSums(x)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Every unit must have a positive finite denominator.", call. = FALSE)
  }
  x / totals
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
    add_check(
      "Frozen input", paste0("sha256-", relative),
      identical(observed, expected), observed
    )
    expected_files <- c(expected_files, relative)
  }
  payloads <- basename(list.files(directory, full.names = TRUE))
  payloads <- sort(setdiff(payloads, "file-checksums.sha256"))
  add_check(
    "Frozen input", "checksum-manifest-complete",
    identical(payloads, sort(expected_files)),
    paste0("payloads=", length(payloads), "; entries=", length(expected_files))
  )
  length(lines)
}

read_wide_profile <- function(path) {
  tab <- utils::read.delim(
    gzfile(path), check.names = FALSE, quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
  stopifnot(identical(names(tab)[1:2], c("MAG", "Taxonomy")))
  unit_columns <- names(tab)[-(1:2)]
  abundance <- t(data.matrix(tab[, unit_columns, drop = FALSE]))
  rownames(abundance) <- unit_columns
  colnames(abundance) <- tab$MAG
  list(
    abundance = abundance,
    taxonomy = tab[, c("MAG", "Taxonomy"), drop = FALSE]
  )
}

checksum_entries <- verify_checksum_manifest(input_dir)
notice_path <- file.path(project_root, "data", "small", "23-data-NOTICE.txt")
notice <- readLines(notice_path, warn = FALSE)
add_check("Frozen input", "notice-present", length(notice) > 10L, notice_path)
add_check(
  "Frozen input", "notice-inference-unit",
  any(grepl("inferential unit is one hot spring", notice, fixed = TRUE)),
  "Hot spring is declared as the inferential unit"
)
add_check(
  "Frozen input", "notice-catalog-boundary",
  any(grepl("not whole-community effects", notice, fixed = TRUE)),
  "Recovered-MAG catalog boundary is declared"
)

resource_manifest <- utils::read.delim(
  file.path(input_dir, "resource-manifest.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
analysis_contract <- utils::read.delim(
  file.path(input_dir, "analysis-contract.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
regime_contract <- utils::read.delim(
  file.path(input_dir, "temperature-regime-contract.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
spring_metadata <- utils::read.delim(
  file.path(input_dir, "spring-metadata.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
sample_ledger <- utils::read.delim(
  file.path(input_dir, "sample-to-spring-ledger.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
equal_input <- read_wide_profile(
  file.path(input_dir, "spring-mag-relative-abundance.tsv.gz")
)
weighted_input <- read_wide_profile(
  file.path(input_dir, "spring-mag-recruitment-weighted.tsv.gz")
)

add_check("Contract", "resource-manifest-rows", nrow(resource_manifest) == 9L, nrow(resource_manifest))
add_check("Contract", "analysis-contract-rows", nrow(analysis_contract) == 20L, nrow(analysis_contract))
add_check(
  "Contract", "inference-unit-locked",
  identical(
    analysis_contract$PrimaryValue[analysis_contract$Parameter == "Inference unit"],
    "Hot spring (n=56)"
  ),
  analysis_contract$PrimaryValue[analysis_contract$Parameter == "Inference unit"]
)
add_check(
  "Contract", "distance-locked",
  identical(
    analysis_contract$PrimaryValue[analysis_contract$Parameter == "Primary distance"],
    "Bray-Curtis"
  ),
  analysis_contract$PrimaryValue[analysis_contract$Parameter == "Primary distance"]
)
add_check(
  "Contract", "marginal-test-locked",
  identical(
    analysis_contract$PrimaryValue[analysis_contract$Parameter == "Primary sums of squares"],
    "Marginal (by=margin)"
  ),
  analysis_contract$PrimaryValue[analysis_contract$Parameter == "Primary sums of squares"]
)
add_check(
  "Contract", "seed-locked",
  identical(
    analysis_contract$PrimaryValue[analysis_contract$Parameter == "Seed"],
    as.character(primary_seed)
  ),
  analysis_contract$PrimaryValue[analysis_contract$Parameter == "Seed"]
)
add_check(
  "Contract", "regime-contract-counts",
  identical(as.integer(regime_contract$Springs), c(3L, 16L, 30L, 4L, 3L)),
  paste(regime_contract$Springs, collapse = "/")
)
add_check(
  "Contract", "regime-rule-predeclared",
  all(regime_contract$RuleLockedBeforeInference),
  paste(regime_contract$RuleLockedBeforeInference, collapse = ",")
)

spring_equal <- close_rows(equal_input$abundance)
spring_weighted <- close_rows(weighted_input$abundance)
add_check("Input shape", "equal-shape", identical(dim(spring_equal), c(56L, 780L)), paste(dim(spring_equal), collapse = "x"))
add_check("Input shape", "weighted-shape", identical(dim(spring_weighted), c(56L, 780L)), paste(dim(spring_weighted), collapse = "x"))
add_check("Input shape", "spring-metadata-rows", nrow(spring_metadata) == 56L, nrow(spring_metadata))
add_check("Input shape", "sample-ledger-rows", nrow(sample_ledger) == 500L, nrow(sample_ledger))
add_check("Input shape", "mag-features", ncol(spring_equal) == 780L, ncol(spring_equal))
add_check(
  "Alignment", "spring-order",
  identical(rownames(spring_equal), spring_metadata$Spring) &&
    identical(rownames(spring_weighted), spring_metadata$Spring),
  "equal/weighted/metadata"
)
add_check(
  "Alignment", "mag-order",
  identical(colnames(spring_equal), colnames(spring_weighted)) &&
    identical(colnames(spring_equal), equal_input$taxonomy$MAG),
  "equal/weighted/taxonomy"
)
add_check("Values", "equal-finite-nonnegative", all(is.finite(spring_equal)) && min(spring_equal) >= 0, range(spring_equal))
add_check("Values", "weighted-finite-nonnegative", all(is.finite(spring_weighted)) && min(spring_weighted) >= 0, range(spring_weighted))
add_check("Values", "equal-closure", max(abs(rowSums(spring_equal) - 1)) < 1e-12, max(abs(rowSums(spring_equal) - 1)))
add_check("Values", "weighted-closure", max(abs(rowSums(spring_weighted) - 1)) < 1e-12, max(abs(rowSums(spring_weighted) - 1)))

# Independently rederive both frozen spring matrices from the Article 22 sample table.
source_dir <- file.path(project_root, "data", "small", "22-diversity-inputs")
source_mag_path <- file.path(source_dir, "mag-relative-abundance.tsv.gz")
source_metadata_path <- file.path(source_dir, "hot-spring-sample-metadata.tsv")
source_recruitment_path <- file.path(source_dir, "mag-recruitment.tsv")
source_expected <- c(
  mag = "8447575babb47ed660c40f8b87a4fb9bebaf0f4b870c02dd573a4f496b74e58b",
  metadata = "999fd1f521712651d903163b4d238adbe2ad256cbd2c74a5174b62e2a608cc23",
  recruitment = "2f5858c591e0e2ca5c8a026edfa97c8de43d94504af0f9426c83db0408a89c3a"
)
source_observed <- c(
  mag = sha256_file(source_mag_path),
  metadata = sha256_file(source_metadata_path),
  recruitment = sha256_file(source_recruitment_path)
)
add_check(
  "Lineage", "article22-source-checksums",
  identical(unname(source_observed), unname(source_expected)),
  paste(source_observed, collapse = ",")
)
source_mag <- read_wide_profile(source_mag_path)$abundance
source_metadata <- utils::read.delim(
  source_metadata_path, check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
source_recruitment <- utils::read.delim(
  source_recruitment_path, check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
source_metadata <- source_metadata[match(rownames(source_mag), source_metadata$sample), , drop = FALSE]
source_recruitment <- source_recruitment[match(rownames(source_mag), source_recruitment$sample_id), , drop = FALSE]
add_check("Lineage", "source-sample-count", nrow(source_mag) == 500L, nrow(source_mag))
add_check("Lineage", "source-spring-count", length(unique(source_metadata$hotspring)) == 56L, length(unique(source_metadata$hotspring)))
add_check("Lineage", "source-feature-count", ncol(source_mag) == 780L, ncol(source_mag))
add_check(
  "Lineage", "source-catalog-relative-sums",
  max(abs(rowSums(source_mag) - 1)) < 2e-5,
  max(abs(rowSums(source_mag) - 1))
)
add_check(
  "Lineage", "source-sample-alignment",
  identical(rownames(source_mag), source_metadata$sample) &&
    identical(rownames(source_mag), source_recruitment$sample_id),
  "sample composition/metadata/recruitment"
)

spring_indices <- split(
  seq_len(nrow(source_metadata)),
  factor(source_metadata$hotspring, levels = rownames(spring_equal))
)
recomputed_equal <- t(vapply(
  spring_indices,
  function(idx) colMeans(source_mag[idx, , drop = FALSE]),
  numeric(ncol(source_mag))
))
recomputed_equal <- close_rows(recomputed_equal)
catalog_read_mass <- source_metadata$filtered_Nsequences * source_recruitment$total_hit_rate
recomputed_weighted <- t(vapply(
  spring_indices,
  function(idx) {
    weights <- catalog_read_mass[idx] / sum(catalog_read_mass[idx])
    colSums(source_mag[idx, , drop = FALSE] * weights)
  },
  numeric(ncol(source_mag))
))
recomputed_weighted <- close_rows(recomputed_weighted)
colnames(recomputed_equal) <- colnames(source_mag)
colnames(recomputed_weighted) <- colnames(source_mag)

equal_max_error <- max(abs(recomputed_equal - spring_equal))
weighted_max_error <- max(abs(recomputed_weighted - spring_weighted))
add_check("Aggregation", "equal-recalculation", equal_max_error < 1e-12, equal_max_error)
add_check("Aggregation", "weighted-recalculation", weighted_max_error < 1e-12, weighted_max_error)
add_check(
  "Aggregation", "sample-ledger-counts",
  identical(
    as.integer(table(factor(sample_ledger$Spring, levels = spring_metadata$Spring))),
    spring_metadata$Samples
  ),
  paste(range(spring_metadata$Samples), collapse = "..")
)
equal_weight_sums <- tapply(sample_ledger$EqualSampleWeight, sample_ledger$Spring, sum)
catalog_weight_sums <- tapply(sample_ledger$CatalogReadWeight, sample_ledger$Spring, sum)
add_check("Aggregation", "equal-weights-close", max(abs(equal_weight_sums - 1)) < 1e-12, max(abs(equal_weight_sums - 1)))
add_check("Aggregation", "catalog-weights-close", max(abs(catalog_weight_sums - 1)) < 1e-12, max(abs(catalog_weight_sums - 1)))
add_check("Aggregation", "sample-count-range", identical(range(spring_metadata$Samples), c(1L, 33L)), paste(range(spring_metadata$Samples), collapse = ".."))

regime_levels <- c("Below 30", "30–<50", "50–<70", "70–<80", "80–100")
spring_metadata$TemperatureRegime <- factor(
  spring_metadata$TemperatureRegime,
  levels = regime_levels,
  ordered = FALSE
)
spring_metadata$BroadRegion <- factor(spring_metadata$BroadRegion)
rownames(spring_metadata) <- spring_metadata$Spring
add_check("Design", "temperature-complete", !anyNA(spring_metadata$MedianTemperatureC), sum(is.na(spring_metadata$MedianTemperatureC)))
add_check("Design", "ph-complete", !anyNA(spring_metadata$MedianPH), sum(is.na(spring_metadata$MedianPH)))
add_check("Design", "conductivity-complete-cases", sum(!is.na(spring_metadata$MedianConductivity)) == 54L, sum(!is.na(spring_metadata$MedianConductivity)))
add_check(
  "Design", "temperature-group-counts",
  identical(as.integer(table(spring_metadata$TemperatureRegime)), c(3L, 16L, 30L, 4L, 3L)),
  paste(table(spring_metadata$TemperatureRegime), collapse = "/")
)
add_check("Design", "broad-regions", nlevels(spring_metadata$BroadRegion) == 6L, levels(spring_metadata$BroadRegion))
model_matrix <- stats::model.matrix(
  ~ BroadRegion + scale(MedianPH) + TemperatureRegime,
  data = spring_metadata
)
add_check(
  "Design", "primary-model-full-rank",
  qr(model_matrix)$rank == ncol(model_matrix),
  paste0("rank=", qr(model_matrix)$rank, "; columns=", ncol(model_matrix))
)

design_balance <- as.data.frame(
  table(spring_metadata$BroadRegion, spring_metadata$TemperatureRegime),
  stringsAsFactors = FALSE
)
names(design_balance) <- c("BroadRegion", "TemperatureRegime", "Springs")
design_balance$BroadRegion <- as.character(design_balance$BroadRegion)
design_balance$TemperatureRegime <- as.character(design_balance$TemperatureRegime)

spring_aggregation_audit <- data.frame(
  Spring = spring_metadata$Spring,
  BroadRegion = spring_metadata$BroadRegion,
  TemperatureRegime = spring_metadata$TemperatureRegime,
  Samples = spring_metadata$Samples,
  EqualWeightSum = as.numeric(equal_weight_sums[spring_metadata$Spring]),
  CatalogReadWeightSum = as.numeric(catalog_weight_sums[spring_metadata$Spring]),
  EqualRecalculationMaxAbsError = apply(
    abs(recomputed_equal - spring_equal), 1L, max
  ),
  WeightedRecalculationMaxAbsError = apply(
    abs(recomputed_weighted - spring_weighted), 1L, max
  ),
  MeanRecruitment = spring_metadata$MeanRecruitment,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

data_lineage <- data.frame(
  Stage = c(
    "Publisher sample table", "Frozen sample table", "Primary inference table",
    "Sensitivity inference table"
  ),
  Source = c(
    "Figshare 30284068 v2 file 61153444",
    "Article 22 checksum-locked MAG profile",
    "Equal-sample spring aggregation",
    "Estimated catalog-read-weighted spring aggregation"
  ),
  Samples = c(500L, 500L, NA_integer_, NA_integer_),
  InferenceUnits = c(NA_integer_, NA_integer_, 56L, 56L),
  Features = 780L,
  Denominator = "Recovered 780-MAG catalog closure",
  InferentialRole = c(
    "Source", "Source", "Primary", "Sensitivity only"
  ),
  stringsAsFactors = FALSE
)

# Restricted permutation matrices ------------------------------------------
make_block_permutations <- function(block, nset, seed) {
  block <- factor(block)
  identity <- seq_along(block)
  out <- matrix(NA_integer_, nrow = nset, ncol = length(block))
  seen <- new.env(hash = TRUE, parent = emptyenv())
  set.seed(seed)
  i <- 0L
  attempts <- 0L
  while (i < nset) {
    attempts <- attempts + 1L
    candidate <- identity
    for (level in levels(block)) {
      idx <- which(block == level)
      candidate[idx] <- sample(idx, length(idx), replace = FALSE)
    }
    if (identical(candidate, identity)) next
    key <- paste(candidate, collapse = ",")
    if (exists(key, envir = seen, inherits = FALSE)) next
    assign(key, TRUE, envir = seen)
    i <- i + 1L
    out[i, ] <- candidate
  }
  list(matrix = out, attempts = attempts)
}

make_free_permutations <- function(n, nset, seed) {
  identity <- seq_len(n)
  out <- matrix(NA_integer_, nrow = nset, ncol = n)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  set.seed(seed)
  i <- 0L
  while (i < nset) {
    candidate <- sample.int(n)
    if (identical(candidate, identity)) next
    key <- paste(candidate, collapse = ",")
    if (exists(key, envir = seen, inherits = FALSE)) next
    assign(key, TRUE, envir = seen)
    i <- i + 1L
    out[i, ] <- candidate
  }
  out
}

primary_permutation <- make_block_permutations(
  spring_metadata$BroadRegion, primary_nperm, primary_seed
)
primary_permutation_matrix <- primary_permutation$matrix
identity_rows <- sum(apply(
  primary_permutation_matrix, 1L,
  function(x) identical(as.integer(x), seq_len(nrow(spring_metadata)))
))
legal_rows <- vapply(
  seq_len(nrow(primary_permutation_matrix)),
  function(i) all(
    spring_metadata$BroadRegion[primary_permutation_matrix[i, ]] ==
      spring_metadata$BroadRegion
  ),
  logical(1L)
)
region_counts <- as.integer(table(spring_metadata$BroadRegion))
primary_possible <- exp(sum(lgamma(region_counts + 1)))
minimum_primary_p <- 1 / (primary_nperm + 1)
add_check("Permutation", "primary-count", nrow(primary_permutation_matrix) == 9999L, nrow(primary_permutation_matrix))
add_check("Permutation", "primary-width", ncol(primary_permutation_matrix) == 56L, ncol(primary_permutation_matrix))
add_check("Permutation", "primary-unique", nrow(unique(primary_permutation_matrix)) == 9999L, nrow(unique(primary_permutation_matrix)))
add_check("Permutation", "primary-no-identity", identity_rows == 0L, identity_rows)
add_check("Permutation", "primary-all-legal", all(legal_rows), sum(legal_rows))
add_check("Permutation", "primary-minimum-p", identical(minimum_primary_p, 0.0001), minimum_primary_p)

permutation_space_audit <- data.frame(
  Scheme = "Primary region-restricted unit permutations",
  Contrast = "Omnibus",
  Seed = primary_seed,
  Units = nrow(spring_metadata),
  Blocks = nlevels(spring_metadata$BroadRegion),
  SharedBlocks = NA_integer_,
  PossibleIncludingIdentity = primary_possible,
  NonIdentityAvailable = primary_possible - 1,
  PermutationsUsed = primary_nperm,
  Exact = FALSE,
  UniqueRows = nrow(unique(primary_permutation_matrix)),
  IdentityRows = identity_rows,
  LegalRows = sum(legal_rows),
  MinimumP = minimum_primary_p,
  MatrixSHA256 = digest::digest(
    primary_permutation_matrix, algo = "sha256", serialize = TRUE
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Primary distances, ordinations, and location tests -----------------------
primary_distance <- vegan::vegdist(spring_equal, method = "bray")
primary_distance_matrix <- as.matrix(primary_distance)
add_check("Distance", "primary-finite", all(is.finite(primary_distance_matrix)), range(primary_distance_matrix))
add_check("Distance", "primary-nonnegative", min(primary_distance_matrix) >= 0, min(primary_distance_matrix))
add_check("Distance", "primary-symmetric", max(abs(primary_distance_matrix - t(primary_distance_matrix))) < 1e-12, max(abs(primary_distance_matrix - t(primary_distance_matrix))))
add_check("Distance", "primary-diagonal-zero", max(abs(diag(primary_distance_matrix))) < 1e-12, max(abs(diag(primary_distance_matrix))))
add_check("Distance", "primary-order", identical(rownames(primary_distance_matrix), spring_metadata$Spring), "distance/metadata")

set.seed(primary_seed)
primary_marginal_fit <- vegan::adonis2(
  primary_distance ~ BroadRegion + scale(MedianPH) + TemperatureRegime,
  data = spring_metadata,
  permutations = primary_permutation_matrix,
  by = "margin",
  parallel = 1
)
set.seed(primary_seed)
sequential_registered_fit <- vegan::adonis2(
  primary_distance ~ BroadRegion + scale(MedianPH) + TemperatureRegime,
  data = spring_metadata,
  permutations = primary_permutation_matrix,
  by = "terms",
  parallel = 1
)
set.seed(primary_seed)
sequential_swapped_fit <- vegan::adonis2(
  primary_distance ~ BroadRegion + TemperatureRegime + scale(MedianPH),
  data = spring_metadata,
  permutations = primary_permutation_matrix,
  by = "terms",
  parallel = 1
)

fit_rows <- function(fit, test, formula_label, primary_term = "TemperatureRegime") {
  terms <- setdiff(rownames(fit), c("Residual", "Total"))
  residual_ss <- fit["Residual", "SumOfSqs"]
  data.frame(
    Test = test,
    Formula = formula_label,
    Term = terms,
    Df = fit[terms, "Df"],
    SumOfSquares = fit[terms, "SumOfSqs"],
    R2Total = fit[terms, "R2"],
    R2Partial = fit[terms, "SumOfSqs"] /
      (fit[terms, "SumOfSqs"] + residual_ss),
    PseudoF = fit[terms, "F"],
    PValue = fit[terms, "Pr(>F)"],
    Permutations = primary_nperm,
    MinimumP = minimum_primary_p,
    PrimaryInference = test == "Marginal primary" & terms == primary_term,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

permanova_audit <- rbind(
  fit_rows(
    primary_marginal_fit,
    "Marginal primary",
    "BroadRegion + scale(MedianPH) + TemperatureRegime"
  ),
  fit_rows(
    sequential_registered_fit,
    "Sequential registered order",
    "BroadRegion -> scale(MedianPH) -> TemperatureRegime"
  ),
  fit_rows(
    sequential_swapped_fit,
    "Sequential swapped order",
    "BroadRegion -> TemperatureRegime -> scale(MedianPH)"
  )
)
primary_temperature <- permanova_audit[
  permanova_audit$PrimaryInference, , drop = FALSE
]
add_check("PERMANOVA", "primary-target-single-row", nrow(primary_temperature) == 1L, nrow(primary_temperature))
add_check("PERMANOVA", "primary-p-grid", abs(primary_temperature$PValue * 10000 - round(primary_temperature$PValue * 10000)) < 1e-8, primary_temperature$PValue)
primary_effect_values <- unlist(
  primary_temperature[, c("R2Total", "R2Partial", "PseudoF", "PValue")],
  use.names = FALSE
)
add_check(
  "PERMANOVA", "primary-effect-finite",
  all(is.finite(primary_effect_values)), primary_effect_values
)
add_check("PERMANOVA", "primary-r2-bounds", primary_temperature$R2Total > 0 && primary_temperature$R2Total < 1 && primary_temperature$R2Partial > 0 && primary_temperature$R2Partial < 1, c(primary_temperature$R2Total, primary_temperature$R2Partial))

pcoa_fit <- ape::pcoa(primary_distance, correction = "lingoes")
pcoa_coordinates <- if (!is.null(pcoa_fit$vectors.cor)) {
  pcoa_fit$vectors.cor
} else {
  pcoa_fit$vectors
}
pcoa_relative <- if ("Rel_corr_eig" %in% names(pcoa_fit$values)) {
  pcoa_fit$values$Rel_corr_eig
} else {
  pcoa_fit$values$Relative_eig
}
pcoa_raw_eigen <- pcoa_fit$values$Eigenvalues
pcoa_negative_count <- sum(pcoa_raw_eigen < -1e-10)
pcoa_negative_mass <- sum(abs(pcoa_raw_eigen[pcoa_raw_eigen < -1e-10]))
add_check("Ordination", "pcoa-finite", all(is.finite(pcoa_coordinates[, 1:2])), range(pcoa_coordinates[, 1:2]))
add_check(
  "Ordination", "pcoa-lingoes-policy",
  grepl("Lingoes", pcoa_fit$note, ignore.case = TRUE) ||
    grepl("no negative eigenvalues", pcoa_fit$note, ignore.case = TRUE),
  pcoa_fit$note
)
add_check("Ordination", "pcoa-negative-audited", pcoa_negative_count >= 0L && is.finite(pcoa_negative_mass), paste(pcoa_negative_count, pcoa_negative_mass))

set.seed(primary_seed)
cap_fit <- vegan::dbrda(
  primary_distance ~ TemperatureRegime + scale(MedianPH) + Condition(BroadRegion),
  data = spring_metadata,
  add = "lingoes"
)
cap_r2 <- vegan::RsquareAdj(cap_fit)
set.seed(primary_seed)
cap_test <- stats::anova(
  cap_fit,
  permutations = primary_permutation_matrix,
  by = "margin",
  parallel = 1
)
cap_sites <- vegan::scores(
  cap_fit, display = "sites", choices = 1:2, scaling = 1
)
cap_biplot <- vegan::scores(
  cap_fit, display = "bp", choices = 1:2, scaling = 1
)
cap_eigen <- cap_fit$CCA$eig
cap_axis_constrained <- cap_eigen[1:2] / sum(cap_eigen)
cap_axis_total <- cap_eigen[1:2] / cap_fit$tot.chi
add_check("Ordination", "cap-finite", all(is.finite(cap_sites)), range(cap_sites))
add_check("Ordination", "cap-adjusted-r2-finite", is.finite(cap_r2$adj.r.squared), cap_r2$adj.r.squared)
add_check("Ordination", "cap-r2-bounds", cap_r2$r.squared > 0 && cap_r2$r.squared < 1, cap_r2$r.squared)
add_check("Ordination", "cap-ph-vector", "scale(MedianPH)" %in% rownames(cap_biplot), rownames(cap_biplot))

cap_terms <- setdiff(rownames(cap_test), "Residual")
cap_residual_ss <- cap_test["Residual", "SumOfSqs"]
cap_rows <- data.frame(
  Test = "Partial CAP marginal",
  Formula = "TemperatureRegime + scale(MedianPH) + Condition(BroadRegion)",
  Term = cap_terms,
  Df = cap_test[cap_terms, "Df"],
  SumOfSquares = cap_test[cap_terms, "SumOfSqs"],
  R2Total = cap_test[cap_terms, "SumOfSqs"] / cap_fit$tot.chi,
  R2Partial = cap_test[cap_terms, "SumOfSqs"] /
    (cap_test[cap_terms, "SumOfSqs"] + cap_residual_ss),
  PseudoF = cap_test[cap_terms, "F"],
  PValue = cap_test[cap_terms, "Pr(>F)"],
  Permutations = primary_nperm,
  MinimumP = minimum_primary_p,
  PrimaryInference = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
permanova_audit <- rbind(permanova_audit, cap_rows)

pcoa_rows <- data.frame(
  Method = "PCoA (Lingoes policy)",
  Spring = rownames(pcoa_coordinates),
  Axis1 = pcoa_coordinates[, 1L],
  Axis2 = pcoa_coordinates[, 2L],
  Axis1VarianceTotal = pcoa_relative[[1L]],
  Axis2VarianceTotal = pcoa_relative[[2L]],
  Axis1VarianceWithinConstrained = NA_real_,
  Axis2VarianceWithinConstrained = NA_real_,
  NegativeEigenvaluesRaw = pcoa_negative_count,
  NegativeEigenvalueMassRaw = pcoa_negative_mass,
  Correction = if (pcoa_negative_count > 0L) {
    "Lingoes applied"
  } else {
    "Lingoes requested; not needed"
  },
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cap_rows_scores <- data.frame(
  Method = "Partial CAP/dbRDA",
  Spring = rownames(cap_sites),
  Axis1 = cap_sites[, 1L],
  Axis2 = cap_sites[, 2L],
  Axis1VarianceTotal = cap_axis_total[[1L]],
  Axis2VarianceTotal = cap_axis_total[[2L]],
  Axis1VarianceWithinConstrained = cap_axis_constrained[[1L]],
  Axis2VarianceWithinConstrained = cap_axis_constrained[[2L]],
  NegativeEigenvaluesRaw = pcoa_negative_count,
  NegativeEigenvalueMassRaw = pcoa_negative_mass,
  Correction = "Lingoes policy in dbrda",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
ordination_scores <- rbind(pcoa_rows, cap_rows_scores)
ordination_scores$BroadRegion <- spring_metadata$BroadRegion[
  match(ordination_scores$Spring, spring_metadata$Spring)
]
ordination_scores$TemperatureRegime <- spring_metadata$TemperatureRegime[
  match(ordination_scores$Spring, spring_metadata$Spring)
]
ordination_scores$MedianTemperatureC <- spring_metadata$MedianTemperatureC[
  match(ordination_scores$Spring, spring_metadata$Spring)
]
ordination_scores$MedianPH <- spring_metadata$MedianPH[
  match(ordination_scores$Spring, spring_metadata$Spring)
]
add_check("Ordination", "score-row-count", nrow(ordination_scores) == 112L, nrow(ordination_scores))
add_check("Ordination", "score-methods", identical(sort(unique(ordination_scores$Method)), sort(c("PCoA (Lingoes policy)", "Partial CAP/dbRDA"))), unique(ordination_scores$Method))

# Multivariate dispersion ---------------------------------------------------
dispersion_model <- vegan::betadisper(
  primary_distance,
  spring_metadata$TemperatureRegime,
  type = "median",
  bias.adjust = TRUE,
  add = "lingoes"
)
set.seed(primary_seed)
dispersion_test <- vegan::permutest(
  dispersion_model,
  permutations = primary_permutation_matrix,
  parallel = 1
)
dispersion_distances <- dispersion_model$distances
dispersion_global <- data.frame(
  RowType = "Global test",
  TemperatureRegime = "All groups",
  N = length(dispersion_distances),
  MeanDistance = mean(dispersion_distances),
  MedianDistance = stats::median(dispersion_distances),
  Df = dispersion_test$tab["Groups", "Df"],
  SumOfSquares = dispersion_test$tab["Groups", "Sum Sq"],
  PseudoF = dispersion_test$tab["Groups", "F"],
  PValue = dispersion_test$tab["Groups", "Pr(>F)"],
  Permutations = primary_nperm,
  MinimumP = minimum_primary_p,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
dispersion_groups <- do.call(
  rbind,
  lapply(regime_levels, function(group) {
    values <- dispersion_distances[spring_metadata$TemperatureRegime == group]
    data.frame(
      RowType = "Group summary",
      TemperatureRegime = group,
      N = length(values),
      MeanDistance = mean(values),
      MedianDistance = stats::median(values),
      Df = NA_real_, SumOfSquares = NA_real_, PseudoF = NA_real_,
      PValue = NA_real_, Permutations = NA_integer_, MinimumP = NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)
dispersion_audit <- rbind(dispersion_global, dispersion_groups)
add_check("PERMDISP", "distance-count", length(dispersion_distances) == 56L, length(dispersion_distances))
add_check("PERMDISP", "distances-finite", all(is.finite(dispersion_distances)) && min(dispersion_distances) >= 0, range(dispersion_distances))
add_check("PERMDISP", "permutation-count", dispersion_global$Permutations == 9999L, dispersion_global$Permutations)
add_check("PERMDISP", "p-grid", abs(dispersion_global$PValue * 10000 - round(dispersion_global$PValue * 10000)) < 1e-8, dispersion_global$PValue)

# Exact/finite pairwise temperature-regime contrasts ----------------------
safe_combinations <- function(indices, size) {
  if (size == 0L) return(list(integer()))
  if (size == length(indices)) return(list(indices))
  utils::combn(indices, size, simplify = FALSE)
}

make_label_allocation_permutations <- function(meta, nset, seed) {
  stopifnot(is.factor(meta$Group), nlevels(meta$Group) == 2L)
  region_indices <- split(seq_len(nrow(meta)), meta$BroadRegion)
  combinations <- lapply(region_indices, function(indices) {
    n_a <- sum(meta$Group[indices] == levels(meta$Group)[[1L]])
    safe_combinations(indices, n_a)
  })
  radices <- vapply(combinations, length, integer(1L))
  possible <- prod(as.numeric(radices))
  nonidentity <- possible - 1
  used <- min(nset, nonidentity)
  if (used == 0L) {
    return(list(
      matrix = matrix(integer(), nrow = 0L, ncol = nrow(meta)),
      possible = possible,
      radices = radices,
      exact = TRUE
    ))
  }
  codes <- if (nonidentity <= nset) {
    seq_len(nonidentity)
  } else {
    set.seed(seed)
    sample.int(nonidentity, used, replace = FALSE)
  }
  permutations <- matrix(
    rep(seq_len(nrow(meta)), each = used),
    nrow = used,
    byrow = FALSE
  )
  for (row_i in seq_along(codes)) {
    code <- codes[[row_i]]
    for (region_i in seq_along(region_indices)) {
      state <- code %% radices[[region_i]]
      code <- code %/% radices[[region_i]]
      destination <- region_indices[[region_i]]
      source_a <- combinations[[region_i]][[state + 1L]]
      source_b <- setdiff(destination, source_a)
      destination_a <- destination[
        meta$Group[destination] == levels(meta$Group)[[1L]]
      ]
      destination_b <- destination[
        meta$Group[destination] == levels(meta$Group)[[2L]]
      ]
      permutations[row_i, destination_a] <- sort(source_a)
      permutations[row_i, destination_b] <- sort(source_b)
    }
  }
  list(
    matrix = permutations,
    possible = possible,
    radices = radices,
    exact = nonidentity <= nset
  )
}

pairwise_rows <- list()
pairwise_permutation_rows <- list()
contrast_index <- 0L
for (i in seq_len(length(regime_levels) - 1L)) {
  for (j in (i + 1L):length(regime_levels)) {
    contrast_index <- contrast_index + 1L
    group_a <- regime_levels[[i]]
    group_b <- regime_levels[[j]]
    contrast <- paste(group_a, "vs", group_b)
    keep <- spring_metadata$TemperatureRegime %in% c(group_a, group_b)
    pair_meta <- droplevels(spring_metadata[keep, , drop = FALSE])
    pair_meta$Group <- factor(
      as.character(pair_meta$TemperatureRegime),
      levels = c(group_a, group_b)
    )
    pair_meta <- pair_meta[
      order(pair_meta$BroadRegion, pair_meta$Group, pair_meta$Spring),
      , drop = FALSE
    ]
    pair_meta$BroadRegion <- droplevels(pair_meta$BroadRegion)
    shared_regions <- sum(vapply(
      split(as.character(pair_meta$Group), pair_meta$BroadRegion),
      function(x) length(unique(x)) == 2L,
      logical(1L)
    ))
    allocation <- make_label_allocation_permutations(
      pair_meta, primary_nperm, primary_seed + contrast_index
    )
    estimable <- shared_regions >= 2L
    ph_a <- pair_meta$MedianPH[pair_meta$Group == group_a]
    ph_b <- pair_meta$MedianPH[pair_meta$Group == group_b]
    pooled_sd <- stats::sd(pair_meta$MedianPH)
    ph_standardized_difference <- if (is.finite(pooled_sd) && pooled_sd > 0) {
      (mean(ph_a) - mean(ph_b)) / pooled_sd
    } else {
      NA_real_
    }
    result <- data.frame(
      Contrast = contrast,
      GroupA = group_a,
      GroupB = group_b,
      Units = nrow(pair_meta),
      Regions = nlevels(pair_meta$BroadRegion),
      SharedRegions = shared_regions,
      PossibleIncludingIdentity = allocation$possible,
      NonIdentityAvailable = allocation$possible - 1,
      PermutationsUsed = if (estimable) nrow(allocation$matrix) else 0L,
      MinimumP = if (estimable) 1 / allocation$possible else NA_real_,
      Exact = if (estimable) allocation$exact else NA,
      Status = if (estimable) "Estimable" else "Not estimable",
      Reason = if (estimable) {
        "At least two regions contain both groups"
      } else {
        "Fewer than two regions contain both groups"
      },
      Df = NA_real_, SumOfSquares = NA_real_, R2Total = NA_real_,
      R2Partial = NA_real_, PseudoF = NA_real_, PValue = NA_real_,
      PAdjustedHolm = NA_real_, RejectHolm05 = NA,
      PHMeanA = mean(ph_a), PHMeanB = mean(ph_b),
      PHStandardizedDifference = ph_standardized_difference,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (estimable) {
      pair_distance <- stats::as.dist(
        primary_distance_matrix[pair_meta$Spring, pair_meta$Spring, drop = FALSE]
      )
      set.seed(primary_seed + contrast_index)
      pair_fit <- vegan::adonis2(
        pair_distance ~ BroadRegion + Group,
        data = pair_meta,
        permutations = allocation$matrix,
        by = "terms",
        parallel = 1
      )
      pair_ss <- pair_fit["Group", "SumOfSqs"]
      pair_residual_ss <- pair_fit["Residual", "SumOfSqs"]
      result$Df <- pair_fit["Group", "Df"]
      result$SumOfSquares <- pair_ss
      result$R2Total <- pair_fit["Group", "R2"]
      result$R2Partial <- pair_ss / (pair_ss + pair_residual_ss)
      result$PseudoF <- pair_fit["Group", "F"]
      result$PValue <- pair_fit["Group", "Pr(>F)"]
    }
    pairwise_rows[[contrast_index]] <- result
    allocation_matrix <- allocation$matrix
    legal_pair_rows <- if (nrow(allocation_matrix) > 0L) {
      vapply(
        seq_len(nrow(allocation_matrix)),
        function(k) all(
          pair_meta$BroadRegion[allocation_matrix[k, ]] ==
            pair_meta$BroadRegion
        ),
        logical(1L)
      )
    } else {
      logical()
    }
    pairwise_permutation_rows[[contrast_index]] <- data.frame(
      Scheme = if (estimable) {
        if (allocation$exact) "Exact label allocation" else "Monte Carlo label allocation"
      } else {
        "Not tested: insufficient region overlap"
      },
      Contrast = contrast,
      Seed = if (estimable && !allocation$exact) primary_seed + contrast_index else NA_integer_,
      Units = nrow(pair_meta),
      Blocks = nlevels(pair_meta$BroadRegion),
      SharedBlocks = shared_regions,
      PossibleIncludingIdentity = allocation$possible,
      NonIdentityAvailable = allocation$possible - 1,
      PermutationsUsed = if (estimable) nrow(allocation_matrix) else 0L,
      Exact = if (estimable) allocation$exact else NA,
      UniqueRows = if (estimable) nrow(unique(allocation_matrix)) else 0L,
      IdentityRows = if (estimable && nrow(allocation_matrix) > 0L) {
        sum(apply(
          allocation_matrix, 1L,
          function(x) identical(as.integer(x), seq_len(nrow(pair_meta)))
        ))
      } else 0L,
      LegalRows = if (estimable) sum(legal_pair_rows) else 0L,
      MinimumP = if (estimable) 1 / allocation$possible else NA_real_,
      MatrixSHA256 = if (estimable) {
        digest::digest(allocation_matrix, algo = "sha256", serialize = TRUE)
      } else NA_character_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}
pairwise_permanova <- do.call(rbind, pairwise_rows)
estimable_index <- pairwise_permanova$Status == "Estimable"
pairwise_permanova$PAdjustedHolm[estimable_index] <- stats::p.adjust(
  pairwise_permanova$PValue[estimable_index], method = "holm"
)
pairwise_permanova$RejectHolm05[estimable_index] <-
  pairwise_permanova$PAdjustedHolm[estimable_index] < 0.05
pairwise_permutation_audit <- do.call(rbind, pairwise_permutation_rows)
permutation_space_audit <- rbind(
  permutation_space_audit,
  pairwise_permutation_audit
)
add_check("Pairwise", "contrast-count", nrow(pairwise_permanova) == 10L, nrow(pairwise_permanova))
add_check("Pairwise", "estimable-count", sum(estimable_index) == 6L, sum(estimable_index))
add_check("Pairwise", "not-estimable-count", sum(!estimable_index) == 4L, sum(!estimable_index))
add_check(
  "Pairwise", "possible-spaces",
  identical(
    pairwise_permanova$PossibleIncludingIdentity,
    c(96, 224, 8, 1, 41126400, 80, 2, 240, 10, 1)
  ),
  paste(pairwise_permanova$PossibleIncludingIdentity, collapse = ",")
)
add_check(
  "Pairwise", "used-permutations",
  identical(
    pairwise_permanova$PermutationsUsed,
    c(95L, 223L, 7L, 0L, 9999L, 79L, 0L, 239L, 0L, 0L)
  ),
  paste(pairwise_permanova$PermutationsUsed, collapse = ",")
)
add_check(
  "Pairwise", "holm-family-six",
  sum(!is.na(pairwise_permanova$PAdjustedHolm)) == 6L,
  sum(!is.na(pairwise_permanova$PAdjustedHolm))
)
add_check(
  "Pairwise", "holm-monotonic",
  all(
    pairwise_permanova$PAdjustedHolm[estimable_index] + 1e-15 >=
      pairwise_permanova$PValue[estimable_index]
  ),
  paste(signif(pairwise_permanova$PAdjustedHolm[estimable_index], 5), collapse = ",")
)
exact_index <- estimable_index & pairwise_permanova$Exact
add_check(
  "Pairwise", "exact-p-grid",
  all(abs(
    pairwise_permanova$PValue[exact_index] *
      pairwise_permanova$PossibleIncludingIdentity[exact_index] -
      round(pairwise_permanova$PValue[exact_index] *
        pairwise_permanova$PossibleIncludingIdentity[exact_index])
  ) < 1e-8),
  paste(signif(pairwise_permanova$PValue[exact_index], 6), collapse = ",")
)
add_check(
  "Pairwise", "monte-carlo-p-grid",
  all(abs(
    pairwise_permanova$PValue[estimable_index & !pairwise_permanova$Exact] * 10000 -
      round(pairwise_permanova$PValue[estimable_index & !pairwise_permanova$Exact] * 10000)
  ) < 1e-8),
  pairwise_permanova$PValue[estimable_index & !pairwise_permanova$Exact]
)

# Predeclared sensitivity analyses -----------------------------------------
run_location <- function(distance, metadata, permutations, model) {
  fit <- switch(
    model,
    primary = vegan::adonis2(
      distance ~ BroadRegion + scale(MedianPH) + TemperatureRegime,
      data = metadata, permutations = permutations, by = "margin", parallel = 1
    ),
    no_ph = vegan::adonis2(
      distance ~ BroadRegion + TemperatureRegime,
      data = metadata, permutations = permutations, by = "margin", parallel = 1
    ),
    conductivity = vegan::adonis2(
      distance ~ BroadRegion + scale(MedianPH) +
        scale(MedianConductivity) + TemperatureRegime,
      data = metadata, permutations = permutations, by = "margin", parallel = 1
    ),
    stop("Unknown model branch: ", model, call. = FALSE)
  )
  ss <- fit["TemperatureRegime", "SumOfSqs"]
  residual_ss <- fit["Residual", "SumOfSqs"]
  c(
    Df = fit["TemperatureRegime", "Df"],
    SumOfSquares = ss,
    R2Total = fit["TemperatureRegime", "R2"],
    R2Partial = ss / (ss + residual_ss),
    PseudoF = fit["TemperatureRegime", "F"],
    PValue = fit["TemperatureRegime", "Pr(>F)"]
  )
}

make_sensitivity_row <- function(
  variant, distance_name, features, units, formula_label,
  permutation_scheme, permutations, stats, distance_rho,
  primary_inference = FALSE, interpretation
) {
  data.frame(
    Variant = variant,
    Distance = distance_name,
    Features = features,
    Units = units,
    Formula = formula_label,
    PermutationScheme = permutation_scheme,
    Permutations = permutations,
    Df = stats[["Df"]],
    SumOfSquares = stats[["SumOfSquares"]],
    R2Total = stats[["R2Total"]],
    R2Partial = stats[["R2Partial"]],
    PseudoF = stats[["PseudoF"]],
    PValue = stats[["PValue"]],
    DistanceSpearmanRho = distance_rho,
    PrimaryInference = primary_inference,
    Interpretation = interpretation,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

primary_stats <- c(
  Df = primary_temperature$Df,
  SumOfSquares = primary_temperature$SumOfSquares,
  R2Total = primary_temperature$R2Total,
  R2Partial = primary_temperature$R2Partial,
  PseudoF = primary_temperature$PseudoF,
  PValue = primary_temperature$PValue
)
sensitivity_rows <- list(
  make_sensitivity_row(
    "Primary", "Bray-Curtis", 780L, 56L,
    "BroadRegion + scale(MedianPH) + TemperatureRegime",
    "Within BroadRegion", primary_nperm, primary_stats, 1,
    TRUE, "Primary marginal location test"
  )
)

jaccard_distance <- vegan::vegdist(
  spring_equal > 0, method = "jaccard", binary = TRUE
)
jaccard_stats <- run_location(
  jaccard_distance, spring_metadata, primary_permutation_matrix, "primary"
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- make_sensitivity_row(
  "Binary membership", "Binary Jaccard", 780L, 56L,
  "BroadRegion + scale(MedianPH) + TemperatureRegime",
  "Within BroadRegion", primary_nperm, jaccard_stats,
  stats::cor(as.vector(primary_distance), as.vector(jaccard_distance), method = "spearman"),
  FALSE, "Presence/absence geometry; sensitivity only"
)

prevalence_keep <- colSums(spring_equal > 0) >= 6L
prevalence_matrix <- close_rows(spring_equal[, prevalence_keep, drop = FALSE])
prevalence_distance <- vegan::vegdist(prevalence_matrix, method = "bray")
prevalence_stats <- run_location(
  prevalence_distance, spring_metadata, primary_permutation_matrix, "primary"
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- make_sensitivity_row(
  "MAG prevalence >=10%", "Bray-Curtis", sum(prevalence_keep), 56L,
  "BroadRegion + scale(MedianPH) + TemperatureRegime",
  "Within BroadRegion", primary_nperm, prevalence_stats,
  stats::cor(as.vector(primary_distance), as.vector(prevalence_distance), method = "spearman"),
  FALSE, "Feature-universe sensitivity"
)

weighted_distance <- vegan::vegdist(spring_weighted, method = "bray")
weighted_stats <- run_location(
  weighted_distance, spring_metadata, primary_permutation_matrix, "primary"
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- make_sensitivity_row(
  "Catalog-read-weighted aggregation", "Bray-Curtis", 780L, 56L,
  "BroadRegion + scale(MedianPH) + TemperatureRegime",
  "Within BroadRegion", primary_nperm, weighted_stats,
  stats::cor(as.vector(primary_distance), as.vector(weighted_distance), method = "spearman"),
  FALSE, "Sampling-weight sensitivity; not the primary unit weighting"
)

no_ph_stats <- run_location(
  primary_distance, spring_metadata, primary_permutation_matrix, "no_ph"
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- make_sensitivity_row(
  "Remove pH covariate", "Bray-Curtis", 780L, 56L,
  "BroadRegion + TemperatureRegime",
  "Within BroadRegion", primary_nperm, no_ph_stats, 1,
  FALSE, "Covariate sensitivity"
)

free_permutation_matrix <- make_free_permutations(
  nrow(spring_metadata), primary_nperm, primary_seed + 100L
)
unrestricted_stats <- run_location(
  primary_distance, spring_metadata, free_permutation_matrix, "primary"
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- make_sensitivity_row(
  "Unrestricted permutations", "Bray-Curtis", 780L, 56L,
  "BroadRegion + scale(MedianPH) + TemperatureRegime",
  "Unrestricted (diagnostic only)", primary_nperm, unrestricted_stats, 1,
  FALSE, "Invalid as primary because geography is ignored"
)

conductivity_keep <- !is.na(spring_metadata$MedianConductivity)
conductivity_meta <- droplevels(spring_metadata[conductivity_keep, , drop = FALSE])
conductivity_distance <- stats::as.dist(
  primary_distance_matrix[
    conductivity_meta$Spring, conductivity_meta$Spring, drop = FALSE
  ]
)
conductivity_permutation_matrix <- make_block_permutations(
  conductivity_meta$BroadRegion, primary_nperm, primary_seed + 200L
)$matrix
conductivity_model_matrix <- stats::model.matrix(
  ~ BroadRegion + scale(MedianPH) + scale(MedianConductivity) +
    TemperatureRegime,
  data = conductivity_meta
)
add_check(
  "Sensitivity", "conductivity-model-full-rank",
  qr(conductivity_model_matrix)$rank == ncol(conductivity_model_matrix),
  paste0("rank=", qr(conductivity_model_matrix)$rank, "; columns=", ncol(conductivity_model_matrix))
)
conductivity_stats <- run_location(
  conductivity_distance, conductivity_meta,
  conductivity_permutation_matrix, "conductivity"
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- make_sensitivity_row(
  "Add conductivity (complete cases)", "Bray-Curtis", 780L, 54L,
  "BroadRegion + scale(MedianPH) + scale(MedianConductivity) + TemperatureRegime",
  "Within BroadRegion", primary_nperm, conductivity_stats,
  1, FALSE, "54-spring complete-case covariate sensitivity"
)

registered_temp <- permanova_audit[
  permanova_audit$Test == "Sequential registered order" &
    permanova_audit$Term == "TemperatureRegime", , drop = FALSE
]
swapped_temp <- permanova_audit[
  permanova_audit$Test == "Sequential swapped order" &
    permanova_audit$Term == "TemperatureRegime", , drop = FALSE
]
for (branch in list(
  list(label = "Sequential registered order", row = registered_temp),
  list(label = "Sequential swapped order", row = swapped_temp)
)) {
  row <- branch$row
  branch_stats <- c(
    Df = row$Df, SumOfSquares = row$SumOfSquares,
    R2Total = row$R2Total, R2Partial = row$R2Partial,
    PseudoF = row$PseudoF, PValue = row$PValue
  )
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- make_sensitivity_row(
    branch$label, "Bray-Curtis", 780L, 56L, row$Formula,
    "Within BroadRegion", primary_nperm, branch_stats, 1,
    FALSE, "Sequential sums of squares diagnostic; not primary inference"
  )
}
sensitivity_audit <- do.call(rbind, sensitivity_rows)
add_check("Sensitivity", "branch-count", nrow(sensitivity_audit) == 9L, nrow(sensitivity_audit))
add_check("Sensitivity", "single-primary", sum(sensitivity_audit$PrimaryInference) == 1L, sum(sensitivity_audit$PrimaryInference))
add_check("Sensitivity", "all-effects-finite", all(is.finite(sensitivity_audit$R2Total)) && all(is.finite(sensitivity_audit$PValue)), range(sensitivity_audit$R2Total))
add_check("Sensitivity", "prevalence-features-reduced", sum(prevalence_keep) < 780L && sum(prevalence_keep) > 0L, sum(prevalence_keep))
add_check("Sensitivity", "free-permutations-unique", nrow(unique(free_permutation_matrix)) == 9999L, nrow(unique(free_permutation_matrix)))
add_check("Sensitivity", "free-permutations-not-primary", !sensitivity_audit$PrimaryInference[sensitivity_audit$Variant == "Unrestricted permutations"], "Diagnostic only")

interpretation_boundaries <- data.frame(
  Topic = c(
    "Inference unit", "PCoA", "Partial CAP/dbRDA", "PERMANOVA",
    "PERMDISP", "Pairwise contrasts", "Permutation p-value",
    "pH and conductivity", "MAG denominator", "Sampling design", "Causality"
  ),
  AuthorizedInterpretation = c(
    "One equally weighted hot spring",
    "Unconstrained Bray-Curtis geometry after Lingoes correction",
    "Variation constrained by temperature regime and pH after conditioning broad region",
    "Marginal temperature-regime association with centroid location after declared adjustment",
    "Differences in distance to group spatial medians",
    "Only six contrasts with overlap in at least two broad regions",
    "Discrete resolution set by the observed plus legal non-identity allocations",
    "Measured environmental covariates; conductivity uses 54 complete springs",
    "Composition within 780 recovered MAGs; recruitment remains separate",
    "Purposeful local-diversity sampling with 1-33 samples per spring",
    "Observed association requiring independent mechanistic or longitudinal evidence"
  ),
  ProhibitedClaim = c(
    "Five hundred independent biological replicates",
    "Significant separation by visual inspection",
    "The first axes explain all community variation",
    "Pure centroid shift when dispersion also differs",
    "Difference in centroid location",
    "A non-significant result for a non-estimable contrast",
    "More decimal precision than the permutation grid supports",
    "Randomized or manipulation variables",
    "Whole-community abundance or recovered-read percentage",
    "Spatially balanced probability survey",
    "Temperature causes the complete community change"
  ),
  stringsAsFactors = FALSE
)

# Publication graphics -----------------------------------------------------
pal_regime <- c(
  "Below 30" = "#0072B2",
  "30–<50" = "#56B4E9",
  "50–<70" = "#009E73",
  "70–<80" = "#E69F00",
  "80–100" = "#D55E00"
)
pal_region <- c(
  "ALK" = "#332288", "ALV" = "#88CCEE", "BLR" = "#44AA99",
  "GBD" = "#117733", "SIE" = "#DDCC77", "YEL" = "#CC6677"
)
shape_region <- c(
  "ALK" = 15, "ALV" = 16, "BLR" = 17,
  "GBD" = 18, "SIE" = 3, "YEL" = 8
)
theme_pub <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.25),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey45"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.05), lineheight = 0.98
      ),
      plot.subtitle = ggplot2::element_text(
        size = ggplot2::rel(0.88), lineheight = 0.98
      )
    )
}
save_pub <- function(plot, stem, width_mm, height_mm) {
  base <- file.path(figure_dir, stem)
  ggplot2::ggsave(
    paste0(base, ".pdf"), plot, width = width_mm, height = height_mm,
    units = "mm", device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".png"), plot, width = width_mm, height = height_mm,
    units = "mm", dpi = 350, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".tiff"), plot, width = width_mm, height = height_mm,
    units = "mm", dpi = 350, compression = "lzw", bg = "white"
  )
}
format_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2), sprintf("%.4f", x))
  )
}

# Figure 1: design balance, unequal local sampling, and legal pairwise spaces.
design_plot <- design_balance
design_plot$TemperatureRegime <- factor(
  design_plot$TemperatureRegime, levels = regime_levels
)
design_plot$BroadRegion <- factor(
  design_plot$BroadRegion, levels = levels(spring_metadata$BroadRegion)
)
p_design_balance <- ggplot2::ggplot(
  design_plot,
  ggplot2::aes(TemperatureRegime, Springs, fill = BroadRegion)
) +
  ggplot2::geom_col(color = "white", linewidth = 0.25) +
  ggplot2::scale_fill_manual(values = pal_region, drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = seq(0, 30, 5), expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::labs(
    title = "A. Geographic balance",
    x = "Median spring temperature (C)",
    y = "Hot springs",
    fill = "Broad region"
  ) +
  theme_pub(9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

sampling_plot <- spring_metadata
p_sampling <- ggplot2::ggplot(
  sampling_plot,
  ggplot2::aes(
    MedianTemperatureC, Samples,
    color = TemperatureRegime, shape = BroadRegion
  )
) +
  ggplot2::geom_hline(
    yintercept = stats::median(sampling_plot$Samples),
    linetype = "dashed", color = "grey45", linewidth = 0.45
  ) +
  ggplot2::geom_point(size = 2.5, alpha = 0.88) +
  ggplot2::scale_color_manual(values = pal_regime, drop = FALSE) +
  ggplot2::scale_shape_manual(values = shape_region, drop = FALSE) +
  ggplot2::scale_y_continuous(
    trans = scales::log1p_trans(), breaks = c(1, 2, 4, 8, 16, 32)
  ) +
  ggplot2::labs(
    title = "B. Unequal local sampling",
    subtitle = "Dashed line: median of four samples per spring",
    x = "Median temperature (C)",
    y = "Samples per spring (log1p scale)",
    color = "Temperature regime",
    shape = "Broad region"
  ) +
  theme_pub(9)

pair_space_plot <- pairwise_permanova
pair_space_plot$Contrast <- factor(
  pair_space_plot$Contrast,
  levels = rev(pair_space_plot$Contrast)
)
pair_space_plot$Decision <- ifelse(
  pair_space_plot$Status == "Estimable", "Estimable", "Not estimable"
)
pair_space_plot$PermutationLabel <- ifelse(
  pair_space_plot$Status != "Estimable",
  "Not estimable",
  ifelse(
    pair_space_plot$Exact,
    paste0("Exact: ", scales::comma(pair_space_plot$PermutationsUsed)),
    paste0("MC: ", scales::comma(pair_space_plot$PermutationsUsed))
  )
)
pair_space_plot$LabelHjust <- ifelse(
  log10(pair_space_plot$PossibleIncludingIdentity) > 7,
  1.05, -0.08
)
p_pair_space <- ggplot2::ggplot(
  pair_space_plot,
  ggplot2::aes(
    log10(PossibleIncludingIdentity), Contrast,
    color = Decision
  )
) +
  ggplot2::geom_point(ggplot2::aes(size = pmax(SharedRegions, 0)), alpha = 0.9) +
  ggplot2::geom_text(
    ggplot2::aes(label = PermutationLabel, hjust = LabelHjust),
    size = 2.45, color = "black"
  ) +
  ggplot2::scale_color_manual(
    values = c("Estimable" = "#009E73", "Not estimable" = "#999999")
  ) +
  ggplot2::scale_size_continuous(range = c(2.0, 4.5), breaks = 0:4) +
  ggplot2::scale_x_continuous(
    limits = c(-0.1, 8.9),
    breaks = 0:8,
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::labs(
    title = "C. Discrete pairwise precision",
    x = "log10 legal allocations",
    y = NULL,
    color = NULL,
    size = "Shared regions"
  ) +
  theme_pub(8.6) +
  ggplot2::theme(legend.position = "bottom")

p_design_balance <- p_design_balance + ggplot2::theme(legend.position = "none")
p_pair_space <- p_pair_space + ggplot2::theme(legend.position = "none")
p_design <- p_design_balance + p_sampling + p_pair_space +
  patchwork::plot_layout(widths = c(0.82, 0.92, 1.35)) +
  patchwork::plot_annotation(
    title = "The permutation design begins with the sampling design",
    subtitle = "Primary inference uses 56 equally weighted hot springs; pairwise tests require overlap in at least two regions"
  )

# Figure 2: unconstrained PCoA versus region-conditioned partial CAP/dbRDA.
pcoa_plot_data <- ordination_scores[
  ordination_scores$Method == "PCoA (Lingoes policy)", , drop = FALSE
]
cap_plot_data <- ordination_scores[
  ordination_scores$Method == "Partial CAP/dbRDA", , drop = FALSE
]
centroid_segments <- function(data) {
  centroids <- stats::aggregate(
    cbind(Axis1, Axis2) ~ TemperatureRegime,
    data = data,
    FUN = mean
  )
  names(centroids)[2:3] <- c("Centroid1", "Centroid2")
  merge(data, centroids, by = "TemperatureRegime", sort = FALSE)
}
pcoa_segments <- centroid_segments(pcoa_plot_data)
cap_segments <- centroid_segments(cap_plot_data)
pcoa_axis1 <- 100 * unique(pcoa_plot_data$Axis1VarianceTotal)[[1L]]
pcoa_axis2 <- 100 * unique(pcoa_plot_data$Axis2VarianceTotal)[[1L]]
cap_axis1 <- 100 * unique(cap_plot_data$Axis1VarianceWithinConstrained)[[1L]]
cap_axis2 <- 100 * unique(cap_plot_data$Axis2VarianceWithinConstrained)[[1L]]

p_pcoa <- ggplot2::ggplot(
  pcoa_plot_data,
  ggplot2::aes(Axis1, Axis2, color = TemperatureRegime, shape = BroadRegion)
) +
  ggplot2::geom_segment(
    data = pcoa_segments,
    ggplot2::aes(xend = Centroid1, yend = Centroid2, color = TemperatureRegime),
    linewidth = 0.35, alpha = 0.18, show.legend = FALSE
  ) +
  ggplot2::geom_point(size = 2.4, alpha = 0.9) +
  ggplot2::scale_color_manual(values = pal_regime, drop = FALSE) +
  ggplot2::scale_shape_manual(values = shape_region, drop = FALSE) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "A. Unconstrained PCoA",
    subtitle = sprintf(
      "Raw negative axes: %d | Lingoes %s",
      pcoa_negative_count,
      if (pcoa_negative_count > 0L) "applied" else "not needed"
    ),
    x = sprintf("PCoA1 (%.1f%% inertia)", pcoa_axis1),
    y = sprintf("PCoA2 (%.1f%% inertia)", pcoa_axis2),
    color = "Temperature regime",
    shape = "Broad region"
  ) +
  theme_pub(9)

ph_arrow <- cap_biplot["scale(MedianPH)", 1:2]
arrow_multiplier <- 1.35
p_cap <- ggplot2::ggplot(
  cap_plot_data,
  ggplot2::aes(Axis1, Axis2, color = TemperatureRegime, shape = BroadRegion)
) +
  ggplot2::geom_segment(
    data = cap_segments,
    ggplot2::aes(xend = Centroid1, yend = Centroid2, color = TemperatureRegime),
    linewidth = 0.35, alpha = 0.18, show.legend = FALSE
  ) +
  ggplot2::geom_point(size = 2.4, alpha = 0.9) +
  ggplot2::annotate(
    "segment", x = 0, y = 0,
    xend = arrow_multiplier * ph_arrow[[1L]],
    yend = arrow_multiplier * ph_arrow[[2L]],
    arrow = grid::arrow(length = grid::unit(2.3, "mm")),
    color = "grey20", linewidth = 0.55
  ) +
  ggplot2::annotate(
    "text",
    x = arrow_multiplier * ph_arrow[[1L]],
    y = arrow_multiplier * ph_arrow[[2L]],
    label = "Median pH", hjust = 0, vjust = -0.25, size = 3
  ) +
  ggplot2::scale_color_manual(values = pal_regime, drop = FALSE) +
  ggplot2::scale_shape_manual(values = shape_region, drop = FALSE) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "B. Partial CAP/dbRDA",
    subtitle = sprintf(
      "Region conditioned | adjusted R2 = %.3f",
      cap_r2$adj.r.squared
    ),
    x = sprintf("dbRDA1 (%.1f%% constrained inertia)", cap_axis1),
    y = sprintf("dbRDA2 (%.1f%% constrained inertia)", cap_axis2),
    color = "Temperature regime",
    shape = "Broad region"
  ) +
  theme_pub(9)

p_ordination <- p_pcoa + p_cap +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "PCoA and CAP answer different geometric questions",
    subtitle = "Centroid spokes are shown; sparse groups have no confidence ellipses"
  ) &
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

# Figure 3: marginal location effects alongside multivariate dispersion.
location_plot <- permanova_audit[
  permanova_audit$Test == "Marginal primary" &
    permanova_audit$Term %in% c("scale(MedianPH)", "TemperatureRegime"),
  , drop = FALSE
]
location_plot$TermLabel <- c(
  "Median pH", "Temperature regime"
)[match(location_plot$Term, c("scale(MedianPH)", "TemperatureRegime"))]
location_long <- rbind(
  data.frame(
    TermLabel = location_plot$TermLabel,
    Effect = "Total R2", Value = location_plot$R2Total,
    PValue = location_plot$PValue
  ),
  data.frame(
    TermLabel = location_plot$TermLabel,
    Effect = "Partial R2", Value = location_plot$R2Partial,
    PValue = location_plot$PValue
  )
)
location_long$TermLabel <- factor(
  location_long$TermLabel,
  levels = c("Median pH", "Temperature regime")
)
p_location <- ggplot2::ggplot(
  location_long,
  ggplot2::aes(TermLabel, Value, fill = Effect)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.72),
    width = 0.64, color = "white"
  ) +
  ggplot2::geom_text(
    data = location_plot,
    ggplot2::aes(
      x = TermLabel, y = R2Partial + 0.025,
      label = paste0("p = ", format_p(PValue))
    ),
    inherit.aes = FALSE, size = 3
  ) +
  ggplot2::scale_fill_manual(
    values = c("Total R2" = "#56B4E9", "Partial R2" = "#0072B2")
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, max(location_long$Value) + 0.08),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "A. Marginal location effects",
    subtitle = "Temperature adjusted for region and pH",
    x = NULL, y = "Explained dissimilarity", fill = NULL
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "bottom")

dispersion_plot_data <- spring_metadata
dispersion_plot_data$DistanceToMedian <- dispersion_distances[
  match(dispersion_plot_data$Spring, names(dispersion_distances))
]
if (anyNA(dispersion_plot_data$DistanceToMedian)) {
  dispersion_plot_data$DistanceToMedian <- dispersion_distances
}
p_dispersion <- ggplot2::ggplot(
  dispersion_plot_data,
  ggplot2::aes(TemperatureRegime, DistanceToMedian, fill = TemperatureRegime)
) +
  ggplot2::geom_boxplot(
    width = 0.58, outlier.shape = NA, alpha = 0.55, linewidth = 0.4
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = TemperatureRegime),
    width = 0.10, height = 0, size = 1.55, alpha = 0.75,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = pal_regime, guide = "none") +
  ggplot2::scale_color_manual(values = pal_regime, guide = "none") +
  ggplot2::labs(
    title = "B. Multivariate dispersion",
    subtitle = sprintf(
      "Bias adjusted | F = %.2f | p = %s",
      dispersion_global$PseudoF, format_p(dispersion_global$PValue)
    ),
    x = "Median spring temperature (C)",
    y = "Distance to group spatial median"
  ) +
  theme_pub(9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

dispersion_message <- if (dispersion_global$PValue < 0.05) {
  "Location and dispersion both vary; PERMANOVA is not an isolated centroid-shift result"
} else {
  "No dispersion difference was detected at alpha = 0.05; power remains group-size dependent"
}
p_location_dispersion <- p_location + p_dispersion +
  patchwork::plot_layout(widths = c(0.88, 1.12)) +
  patchwork::plot_annotation(
    title = "Location and within-group spread must be read together",
    subtitle = dispersion_message
  )

# Figure 4: estimability-aware pairwise results and all preregistered audits.
pair_result_plot <- pairwise_permanova
pair_result_plot$Decision <- ifelse(
  pair_result_plot$Status != "Estimable",
  "Not estimable",
  ifelse(pair_result_plot$RejectHolm05, "Holm < 0.05", "Holm >= 0.05")
)
pair_result_plot$PlotR2 <- ifelse(
  pair_result_plot$Status == "Estimable", pair_result_plot$R2Partial, 0
)
pair_result_plot$Label <- ifelse(
  pair_result_plot$Status != "Estimable",
  "Not estimable",
  paste0("Holm p = ", format_p(pair_result_plot$PAdjustedHolm))
)
pair_result_plot$Contrast <- factor(
  pair_result_plot$Contrast,
  levels = rev(pair_result_plot$Contrast)
)
p_pair_result <- ggplot2::ggplot(
  pair_result_plot,
  ggplot2::aes(PlotR2, Contrast, color = Decision, shape = Decision)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = PlotR2, yend = Contrast),
    color = "grey82", linewidth = 0.65
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label),
    hjust = -0.08, size = 2.75, color = "black"
  ) +
  ggplot2::scale_color_manual(values = c(
    "Holm < 0.05" = "#009E73",
    "Holm >= 0.05" = "#D55E00",
    "Not estimable" = "#999999"
  )) +
  ggplot2::scale_shape_manual(values = c(
    "Holm < 0.05" = 16,
    "Holm >= 0.05" = 17,
    "Not estimable" = 4
  )) +
  ggplot2::scale_x_continuous(
    limits = c(0, max(pair_result_plot$PlotR2, na.rm = TRUE) + 0.28),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "A. Pairwise results respect estimability",
    subtitle = "Holm adjustment covers the six tested contrasts only",
    x = "Region-adjusted partial R2",
    y = NULL, color = NULL, shape = NULL
  ) +
  theme_pub(8.8) +
  ggplot2::theme(legend.position = "bottom")

sensitivity_plot <- sensitivity_audit
sensitivity_plot$Role <- ifelse(
  sensitivity_plot$PrimaryInference,
  "Primary",
  ifelse(
    sensitivity_plot$Variant == "Unrestricted permutations",
    "Invalid diagnostic",
    "Sensitivity"
  )
)
sensitivity_plot$Variant <- factor(
  sensitivity_plot$Variant,
  levels = rev(sensitivity_plot$Variant)
)
sensitivity_plot$Label <- paste0("p = ", format_p(sensitivity_plot$PValue))
p_sensitivity <- ggplot2::ggplot(
  sensitivity_plot,
  ggplot2::aes(R2Total, Variant, color = Role, shape = Role)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = R2Total, yend = Variant),
    color = "grey82", linewidth = 0.65
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label),
    hjust = -0.08, size = 2.7, color = "black"
  ) +
  ggplot2::scale_color_manual(values = c(
    "Primary" = "#0072B2",
    "Sensitivity" = "#009E73",
    "Invalid diagnostic" = "#999999"
  )) +
  ggplot2::scale_shape_manual(values = c(
    "Primary" = 16,
    "Sensitivity" = 17,
    "Invalid diagnostic" = 4
  )) +
  ggplot2::scale_x_continuous(
    limits = c(0, max(sensitivity_plot$R2Total) + 0.065),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "B. The primary result is audited, not replaced",
    subtitle = "Unrestricted permutations are diagnostic-only",
    x = "Temperature-regime total R2",
    y = NULL, color = NULL, shape = NULL
  ) +
  theme_pub(8.8) +
  ggplot2::theme(legend.position = "bottom")

p_pairwise_sensitivity <- p_pair_result + p_sensitivity +
  patchwork::plot_layout(widths = c(1.05, 0.95)) +
  patchwork::plot_annotation(
    title = "Effect size, p-value resolution, and sensitivity form one result",
    subtitle = "Not estimable is a design conclusion, not a non-significant result"
  )

figure_specs <- data.frame(
  Stem = c(
    "23-design-permutation-space", "23-pcoa-cap",
    "23-permanova-dispersion", "23-pairwise-sensitivity"
  ),
  WidthMM = c(225, 210, 200, 225),
  HeightMM = c(132, 126, 116, 138),
  stringsAsFactors = FALSE
)
plots <- list(
  p_design, p_ordination, p_location_dispersion, p_pairwise_sensitivity
)
for (i in seq_len(nrow(figure_specs))) {
  save_pub(
    plots[[i]], figure_specs$Stem[[i]],
    figure_specs$WidthMM[[i]], figure_specs$HeightMM[[i]]
  )
}

visible_labels <- c(
  "The permutation design begins with the sampling design",
  "Temperature groups are geographically sparse",
  "Local sampling intensity is unequal",
  "Pairwise precision is discrete",
  "PCoA and CAP answer different geometric questions",
  "Unconstrained PCoA", "Partial CAP/dbRDA", "Median pH",
  "Centroid location and within-group spread must be read together",
  "Marginal PERMANOVA location effects", "PERMDISP around spatial medians",
  "Pairwise results respect estimability",
  "The primary result is audited, not replaced",
  regime_levels, levels(spring_metadata$BroadRegion),
  pairwise_permanova$Contrast, sensitivity_audit$Variant
)
add_check(
  "Graphics", "visible-labels-no-cjk",
  !any(grepl("[\u4e00-\u9fff]", visible_labels, perl = TRUE)),
  "All visible plot labels are English"
)
for (i in seq_len(nrow(figure_specs))) {
  stem <- figure_specs$Stem[[i]]
  width_px <- round(figure_specs$WidthMM[[i]] / 25.4 * 350)
  height_px <- round(figure_specs$HeightMM[[i]] / 25.4 * 350)
  for (extension in c("pdf", "png", "tiff")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    add_check(
      "Graphics", paste0(stem, "-", extension, "-exists"),
      file.exists(path) && file.info(path)$size > 10000,
      if (file.exists(path)) file.info(path)$size else "missing"
    )
  }
  for (extension in c("png", "tiff")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    info <- magick::image_info(magick::image_read(path))
    add_check(
      "Graphics", paste0(stem, "-", extension, "-350dpi-pixels"),
      abs(info$width[[1L]] - width_px) <= 2L &&
        abs(info$height[[1L]] - height_px) <= 2L,
      paste0(
        info$width[[1L]], "x", info$height[[1L]],
        "; expected ", width_px, "x", height_px
      )
    )
  }
}

# Machine-readable outputs -------------------------------------------------
write_tsv(data_lineage, file.path(output_dir, "data-lineage.tsv"))
write_tsv(
  spring_aggregation_audit,
  file.path(output_dir, "spring-aggregation-audit.tsv")
)
write_tsv(design_balance, file.path(output_dir, "design-balance.tsv"))
write_tsv(
  permutation_space_audit,
  file.path(output_dir, "permutation-space-audit.tsv")
)
write_tsv(ordination_scores, file.path(output_dir, "ordination-scores.tsv"))
write_tsv(permanova_audit, file.path(output_dir, "permanova-audit.tsv"))
write_tsv(dispersion_audit, file.path(output_dir, "dispersion-audit.tsv"))
write_tsv(
  pairwise_permanova,
  file.path(output_dir, "pairwise-permanova.tsv")
)
write_tsv(sensitivity_audit, file.path(output_dir, "sensitivity-audit.tsv"))
write_tsv(
  interpretation_boundaries,
  file.path(output_dir, "interpretation-boundaries.tsv")
)

table_outputs <- c(
  "data-lineage.tsv", "spring-aggregation-audit.tsv", "design-balance.tsv",
  "permutation-space-audit.tsv", "ordination-scores.tsv",
  "permanova-audit.tsv", "dispersion-audit.tsv", "pairwise-permanova.tsv",
  "sensitivity-audit.tsv", "interpretation-boundaries.tsv"
)
for (relative in table_outputs) {
  path <- file.path(output_dir, relative)
  add_check(
    "Output", paste0(relative, "-exists"),
    file.exists(path) && file.info(path)$size > 50,
    if (file.exists(path)) file.info(path)$size else "missing"
  )
}

write_tsv(checks, file.path(output_dir, "validation-audit.tsv"))
checks_failed <- sum(checks$Status == "FAIL")
checks_passed <- sum(checks$Status == "PASS")
summary <- list(
  status = if (checks_failed == 0L) "passed" else "failed",
  source_samples = nrow(source_mag),
  inference_units = nrow(spring_metadata),
  mag_features = ncol(spring_equal),
  temperature_groups = nlevels(spring_metadata$TemperatureRegime),
  temperature_group_counts = as.list(
    stats::setNames(
      as.integer(table(spring_metadata$TemperatureRegime)),
      levels(spring_metadata$TemperatureRegime)
    )
  ),
  broad_regions = nlevels(spring_metadata$BroadRegion),
  primary_distance = "Bray-Curtis",
  primary_formula = "BroadRegion + scale(MedianPH) + TemperatureRegime",
  primary_sums_of_squares = "marginal",
  primary_permutations = primary_nperm,
  primary_minimum_p = minimum_primary_p,
  primary_permutation_sha256 = permutation_space_audit$MatrixSHA256[[1L]],
  primary_temperature_r2_total = unname(primary_temperature$R2Total),
  primary_temperature_r2_partial = unname(primary_temperature$R2Partial),
  primary_temperature_pseudo_f = unname(primary_temperature$PseudoF),
  primary_temperature_pvalue = unname(primary_temperature$PValue),
  pcoa_raw_negative_eigenvalues = pcoa_negative_count,
  pcoa_raw_negative_eigenvalue_mass = pcoa_negative_mass,
  cap_r2 = unname(cap_r2$r.squared),
  cap_adjusted_r2 = unname(cap_r2$adj.r.squared),
  dispersion_pseudo_f = unname(dispersion_global$PseudoF),
  dispersion_pvalue = unname(dispersion_global$PValue),
  pairwise_estimable = sum(estimable_index),
  pairwise_not_estimable = sum(!estimable_index),
  pairwise_holm_rejections = sum(
    pairwise_permanova$RejectHolm05, na.rm = TRUE
  ),
  sensitivity_branches = nrow(sensitivity_audit),
  prevalence10_features = sum(prevalence_keep),
  conductivity_complete_cases = sum(conductivity_keep),
  checksum_entries = checksum_entries,
  checks_passed = checks_passed,
  checks_failed = checks_failed,
  random_seed = primary_seed,
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  package_versions = as.list(vapply(
    packages,
    function(pkg) as.character(utils::packageVersion(pkg)),
    character(1L)
  ))
)
jsonlite::write_json(
  summary,
  file.path(output_dir, "validation-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 10
)

log_lines <- c(
  "Article 23 ordination/PERMANOVA/dispersion validation",
  paste0("Status: ", summary$status),
  paste0("Random seed: ", summary$random_seed),
  paste0("Frozen checksum entries: ", summary$checksum_entries),
  paste0(
    "Source samples / inference springs / MAGs: ",
    summary$source_samples, " / ", summary$inference_units, " / ",
    summary$mag_features
  ),
  paste0(
    "Temperature groups: ",
    paste(unlist(summary$temperature_group_counts), collapse = "/")
  ),
  paste0(
    "Primary temperature regime: R2 total=",
    sprintf("%.6f", summary$primary_temperature_r2_total),
    ", R2 partial=", sprintf("%.6f", summary$primary_temperature_r2_partial),
    ", pseudo-F=", sprintf("%.6f", summary$primary_temperature_pseudo_f),
    ", p=", format_p(summary$primary_temperature_pvalue)
  ),
  paste0(
    "PERMDISP: pseudo-F=", sprintf("%.6f", summary$dispersion_pseudo_f),
    ", p=", format_p(summary$dispersion_pvalue)
  ),
  paste0(
    "Partial CAP R2 / adjusted R2: ",
    sprintf("%.6f", summary$cap_r2), " / ",
    sprintf("%.6f", summary$cap_adjusted_r2)
  ),
  paste0(
    "Pairwise estimable / not estimable / Holm rejections: ",
    summary$pairwise_estimable, " / ", summary$pairwise_not_estimable,
    " / ", summary$pairwise_holm_rejections
  ),
  paste0(
    "Sensitivity branches / prevalence>=10% MAGs / conductivity cases: ",
    summary$sensitivity_branches, " / ", summary$prevalence10_features,
    " / ", summary$conductivity_complete_cases
  ),
  paste0("Checks passed/failed: ", checks_passed, "/", checks_failed)
)
writeLines(log_lines, file.path(output_dir, "validation.log"), useBytes = TRUE)

cat(paste(log_lines, collapse = "\n"), "\n", sep = "")
if (checks_failed > 0L) {
  failed <- checks[checks$Status == "FAIL", , drop = FALSE]
  print(failed, row.names = FALSE)
  stop("Article 23 validation failed.", call. = FALSE)
}
