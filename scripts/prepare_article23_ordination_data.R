#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, digits = 17)
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
required <- c("source-dir", "output-dir", "notice")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Missing R package: digest", call. = FALSE)
}

source_dir <- normalizePath(args[["source-dir"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
notice_path <- normalizePath(args[["notice"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(notice_path), recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

source_paths <- c(
  mag = file.path(source_dir, "mag-relative-abundance.tsv.gz"),
  metadata = file.path(source_dir, "hot-spring-sample-metadata.tsv"),
  recruitment = file.path(source_dir, "mag-recruitment.tsv")
)
if (!all(file.exists(source_paths))) {
  stop("Article 22 frozen MAG inputs are incomplete.", call. = FALSE)
}
expected_sha256 <- c(
  mag = "8447575babb47ed660c40f8b87a4fb9bebaf0f4b870c02dd573a4f496b74e58b",
  metadata = "999fd1f521712651d903163b4d238adbe2ad256cbd2c74a5174b62e2a608cc23",
  recruitment = "2f5858c591e0e2ca5c8a026edfa97c8de43d94504af0f9426c83db0408a89c3a"
)
observed_sha256 <- vapply(source_paths, sha256_file, character(1L))
if (!identical(unname(observed_sha256), unname(expected_sha256))) {
  stop("Article 22 MAG input checksum mismatch.", call. = FALSE)
}

mag_table <- utils::read.delim(
  source_paths[["mag"]], check.names = FALSE, quote = "", comment.char = ""
)
metadata <- utils::read.delim(
  source_paths[["metadata"]], check.names = FALSE, quote = "", comment.char = ""
)
recruitment <- utils::read.delim(
  source_paths[["recruitment"]], check.names = FALSE, quote = "", comment.char = ""
)

stopifnot(
  identical(names(mag_table)[1:2], c("MAG", "Taxonomy")),
  nrow(mag_table) == 780L,
  ncol(mag_table) == 502L,
  nrow(metadata) == 500L,
  nrow(recruitment) == 500L,
  !anyDuplicated(mag_table$MAG),
  !anyDuplicated(metadata$sample),
  !anyDuplicated(recruitment$sample_id)
)

sample_ids <- names(mag_table)[-(1:2)]
metadata <- metadata[match(sample_ids, metadata$sample), , drop = FALSE]
recruitment <- recruitment[match(sample_ids, recruitment$sample_id), , drop = FALSE]
stopifnot(
  identical(metadata$sample, sample_ids),
  identical(recruitment$sample_id, sample_ids),
  all(is.finite(recruitment$total_hit_rate)),
  max(abs(
    metadata$metagenome2MAG_read_recruitment_rate - recruitment$total_hit_rate
  )) < 1e-12
)

sample_matrix <- as.matrix(mag_table[, -(1:2), drop = FALSE])
storage.mode(sample_matrix) <- "double"
rownames(sample_matrix) <- mag_table$MAG
stopifnot(
  !anyNA(sample_matrix),
  min(sample_matrix) >= 0,
  max(abs(colSums(sample_matrix) - 1)) < 2e-5,
  length(unique(metadata$hotspring)) == 56L
)

spring_ids <- sort(unique(metadata$hotspring))
spring_indices <- split(seq_len(nrow(metadata)), factor(metadata$hotspring, levels = spring_ids))

equal_matrix <- vapply(
  spring_indices,
  function(idx) rowMeans(sample_matrix[, idx, drop = FALSE]),
  numeric(nrow(sample_matrix))
)
rownames(equal_matrix) <- rownames(sample_matrix)
colnames(equal_matrix) <- spring_ids
equal_matrix <- sweep(equal_matrix, 2L, colSums(equal_matrix), "/")

catalog_read_mass <- metadata$filtered_Nsequences * recruitment$total_hit_rate
stopifnot(all(is.finite(catalog_read_mass)), min(catalog_read_mass) > 0)
recruitment_weighted_matrix <- vapply(
  spring_indices,
  function(idx) {
    weights <- catalog_read_mass[idx] / sum(catalog_read_mass[idx])
    as.numeric(sample_matrix[, idx, drop = FALSE] %*% weights)
  },
  numeric(nrow(sample_matrix))
)
rownames(recruitment_weighted_matrix) <- rownames(sample_matrix)
colnames(recruitment_weighted_matrix) <- spring_ids
recruitment_weighted_matrix <- sweep(
  recruitment_weighted_matrix, 2L, colSums(recruitment_weighted_matrix), "/"
)
stopifnot(
  max(abs(colSums(equal_matrix) - 1)) < 1e-12,
  max(abs(colSums(recruitment_weighted_matrix) - 1)) < 1e-12
)

one_value <- function(x, label) {
  observed <- unique(x[!is.na(x) & nzchar(as.character(x))])
  if (length(observed) != 1L) {
    stop(label, " is not constant within a hot spring.", call. = FALSE)
  }
  observed[[1L]]
}
one_value_or <- function(x, fallback, label) {
  observed <- unique(x[!is.na(x) & nzchar(as.character(x))])
  if (length(observed) == 0L) return(fallback)
  if (length(observed) != 1L) {
    stop(label, " is not constant within a hot spring.", call. = FALSE)
  }
  observed[[1L]]
}
safe_median <- function(x) {
  if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
}

spring_rows <- lapply(spring_ids, function(spring_id) {
  idx <- spring_indices[[spring_id]]
  years <- sort(unique(metadata$year[idx]))
  data.frame(
    Spring = spring_id,
    HotSpringName = one_value_or(
      metadata$hotspring_common_name[idx], spring_id, "Hot-spring name"
    ),
    BroadRegion = one_value(metadata$broad_region_short[idx], "Broad region"),
    BroadRegionName = one_value(metadata$region[idx], "Broad region name"),
    LatitudeMedian = safe_median(metadata$latitude[idx]),
    LongitudeMedian = safe_median(metadata$longitude[idx]),
    Samples = length(idx),
    SamplingYears = paste(years, collapse = ","),
    YearsObserved = length(years),
    MedianPH = safe_median(metadata$pH[idx]),
    PHSamples = sum(!is.na(metadata$pH[idx])),
    MedianTemperatureC = safe_median(metadata$temperature[idx]),
    TemperatureSamples = sum(!is.na(metadata$temperature[idx])),
    MedianConductivity = safe_median(metadata$conductivity[idx]),
    ConductivitySamples = sum(!is.na(metadata$conductivity[idx])),
    MeanRecruitment = mean(recruitment$total_hit_rate[idx]),
    MedianRecruitment = stats::median(recruitment$total_hit_rate[idx]),
    TotalFilteredSequences = sum(metadata$filtered_Nsequences[idx]),
    EstimatedCatalogReads = sum(catalog_read_mass[idx]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
spring_metadata <- do.call(rbind, spring_rows)

regime_levels <- c("Below 30", "30–<50", "50–<70", "70–<80", "80–100")
spring_metadata$TemperatureRegime <- cut(
  spring_metadata$MedianTemperatureC,
  breaks = c(-Inf, 30, 50, 70, 80, Inf),
  right = FALSE,
  labels = regime_levels,
  ordered_result = TRUE
)
stopifnot(
  !anyNA(spring_metadata$MedianPH),
  !anyNA(spring_metadata$MedianTemperatureC),
  identical(as.integer(table(spring_metadata$TemperatureRegime)), c(3L, 16L, 30L, 4L, 3L)),
  sum(is.na(spring_metadata$MedianConductivity)) == 2L
)

sample_count_by_spring <- ave(
  seq_len(nrow(metadata)), metadata$hotspring,
  FUN = length
)
equal_weight <- 1 / sample_count_by_spring
catalog_weight <- ave(
  catalog_read_mass, metadata$hotspring,
  FUN = function(x) x / sum(x)
)
sample_ledger <- data.frame(
  Sample = metadata$sample,
  Spring = metadata$hotspring,
  BroadRegion = metadata$broad_region_short,
  TemperatureC = metadata$temperature,
  PH = metadata$pH,
  Conductivity = metadata$conductivity,
  FilteredSequences = metadata$filtered_Nsequences,
  Recruitment = recruitment$total_hit_rate,
  EstimatedCatalogReads = catalog_read_mass,
  EqualSampleWeight = equal_weight,
  CatalogReadWeight = catalog_weight,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(
  max(abs(tapply(sample_ledger$EqualSampleWeight, sample_ledger$Spring, sum) - 1)) < 1e-12,
  max(abs(tapply(sample_ledger$CatalogReadWeight, sample_ledger$Spring, sum) - 1)) < 1e-12
)

write_wide_gz <- function(x, path) {
  con <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(con), add = TRUE)
  out <- data.frame(
    MAG = mag_table$MAG,
    Taxonomy = mag_table$Taxonomy,
    x,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    out, file = con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
  invisible(path)
}

equal_path <- file.path(output_dir, "spring-mag-relative-abundance.tsv.gz")
weighted_path <- file.path(output_dir, "spring-mag-recruitment-weighted.tsv.gz")
metadata_path <- file.path(output_dir, "spring-metadata.tsv")
ledger_path <- file.path(output_dir, "sample-to-spring-ledger.tsv")
analysis_contract_path <- file.path(output_dir, "analysis-contract.tsv")
regime_contract_path <- file.path(output_dir, "temperature-regime-contract.tsv")
resource_manifest_path <- file.path(output_dir, "resource-manifest.tsv")

write_wide_gz(equal_matrix, equal_path)
write_wide_gz(recruitment_weighted_matrix, weighted_path)
utils::write.table(
  spring_metadata, metadata_path, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)
utils::write.table(
  sample_ledger, ledger_path, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

analysis_contract <- data.frame(
  Parameter = c(
    "Inference unit", "Primary aggregation", "Primary distance",
    "Primary formula order", "Primary sums of squares", "Primary target term",
    "Permutation blocks", "Primary permutations", "Seed", "PCoA correction",
    "Partial CAP formula", "Dispersion test", "Pairwise estimability",
    "Pairwise permutations", "Pairwise multiplicity", "Alpha",
    "Prevalence sensitivity", "Aggregation sensitivity",
    "Distance sensitivity", "Covariate sensitivity"
  ),
  PrimaryValue = c(
    "Hot spring (n=56)",
    "Equal mean of sample-relative MAG profiles, then closure",
    "Bray-Curtis",
    "BroadRegion + scale(MedianPH) + TemperatureRegime",
    "Marginal (by=margin)",
    "TemperatureRegime",
    "Permute springs within BroadRegion only",
    "9999 unique non-identity unit permutations",
    "20260723",
    "Lingoes",
    "TemperatureRegime + scale(MedianPH) + Condition(BroadRegion)",
    "Spatial median; bias.adjust=TRUE; same restricted matrix",
    "At least two BroadRegions contain both contrast groups",
    "All non-identity label allocations if <9999; otherwise 9999 unique",
    "Holm across six estimable contrasts",
    "0.05",
    "MAG presence in at least 10% of springs (>=6/56)",
    "Estimated catalog-read weights = filtered sequences x recruitment",
    "Binary Jaccard",
    "Remove pH; add conductivity in 54-spring complete cases"
  ),
  Role = c(
    rep("Primary", 16L),
    rep("Sensitivity", 4L)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
utils::write.table(
  analysis_contract, analysis_contract_path, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

regime_contract <- data.frame(
  Order = seq_along(regime_levels),
  Label = regime_levels,
  LowerInclusiveC = c(NA, 30, 50, 70, 80),
  UpperExclusiveC = c(30, 50, 70, 80, NA),
  Springs = as.integer(table(spring_metadata$TemperatureRegime)),
  RuleLockedBeforeInference = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
utils::write.table(
  regime_contract, regime_contract_path, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

derived_paths <- c(
  equal_path, weighted_path, metadata_path, ledger_path,
  analysis_contract_path, regime_contract_path
)
resource_manifest <- data.frame(
  Resource = c(
    "Sample-level MAG composition", "Sample metadata", "Sample recruitment",
    "Equal-sample spring composition", "Catalog-read-weighted spring composition",
    "Spring metadata", "Sample-to-spring ledger", "Analysis contract",
    "Temperature-regime contract"
  ),
  Source = c(
    rep("Article 22 frozen Figshare 30284068 v2 payload", 3L),
    rep("Deterministic Article 23 derivation", 6L)
  ),
  File = c(basename(source_paths), basename(derived_paths)),
  Bytes = c(
    vapply(source_paths, function(path) as.numeric(file.info(path)$size), numeric(1L)),
    vapply(derived_paths, function(path) as.numeric(file.info(path)$size), numeric(1L))
  ),
  SHA256 = c(observed_sha256, vapply(derived_paths, sha256_file, character(1L))),
  Rows = c(780L, 500L, 500L, 780L, 780L, 56L, 500L, nrow(analysis_contract), 5L),
  Units = c(500L, 500L, 500L, 56L, 56L, 56L, 500L, 56L, 56L),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
utils::write.table(
  resource_manifest, resource_manifest_path, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

notice <- c(
  "Article 23 data notice",
  "",
  "- Source: Korchagina et al. 2026, Scientific Data.",
  "  DOI: 10.1038/s41597-026-07139-w.",
  "- Publisher payload: Figshare record 30284068 v2, files 61153444,",
  "  61153429, and 61153471, released under CC BY 4.0.",
  "- Article 23 reuses the checksum-locked Article 22 frozen tables; routine QA",
  "  does not access Figshare, FASTQ, BAM, MAG FASTA, or the network.",
  "- The inferential unit is one hot spring (n=56). The primary profile is the",
  "  equal mean of sample-relative 780-MAG profiles within each spring, reclosed",
  "  to one. This prevents 1-33 local samples from silently becoming unequal",
  "  inferential weights.",
  "- The catalog-read-weighted profile is sensitivity-only. Its sample weight is",
  "  filtered read count multiplied by the reported MAG recruitment fraction.",
  "- MAG abundances are closed within the recovered-genome catalog. Mean sample",
  "  recruitment is approximately 13%; results are not whole-community effects.",
  "- Sampling intentionally maximized local diversity and was not a random,",
  "  spatially balanced survey. Associations must not be described as causal."
)
writeLines(notice, notice_path, useBytes = TRUE)

payloads <- sort(setdiff(
  list.files(output_dir, recursive = FALSE, full.names = TRUE),
  file.path(output_dir, "file-checksums.sha256")
))
checksum_lines <- vapply(
  payloads,
  function(path) paste(sha256_file(path), basename(path)),
  character(1L)
)
writeLines(
  checksum_lines,
  file.path(output_dir, "file-checksums.sha256"),
  useBytes = TRUE
)

cat("Prepared Article 23 ordination inputs\n")
cat("Source samples: ", nrow(metadata), "\n", sep = "")
cat("Inference units: ", nrow(spring_metadata), "\n", sep = "")
cat("MAG features: ", nrow(equal_matrix), "\n", sep = "")
cat(
  "Temperature groups: ",
  paste(as.integer(table(spring_metadata$TemperatureRegime)), collapse = "/"),
  "\n", sep = ""
)
cat("Conductivity complete cases: ", sum(!is.na(spring_metadata$MedianConductivity)), "\n", sep = "")
