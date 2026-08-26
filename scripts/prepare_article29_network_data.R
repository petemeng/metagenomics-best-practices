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
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}

source_paths <- c(
  composition = file.path(source_dir, "spring-mag-relative-abundance.tsv.gz"),
  metadata = file.path(source_dir, "spring-metadata.tsv")
)
if (!all(file.exists(source_paths))) {
  stop("Article 23 spring-level inputs are incomplete.", call. = FALSE)
}
expected_source_sha256 <- c(
  composition = "b37278a7640c76fca357ad16105d835f745693e9799b4e163b6171cbb6df63a8",
  metadata = "943b7d213d0fa6eb4153e9812989f403d9063ceaf468b1137e956aee0ad9cee6"
)
observed_source_sha256 <- vapply(source_paths, sha256_file, character(1L))
if (!identical(unname(observed_source_sha256), unname(expected_source_sha256))) {
  stop("Article 23 spring-level input checksum mismatch.", call. = FALSE)
}

composition <- utils::read.delim(
  gzfile(source_paths[["composition"]]),
  check.names = FALSE, quote = "", comment.char = "",
  stringsAsFactors = FALSE
)
metadata <- utils::read.delim(
  source_paths[["metadata"]],
  check.names = FALSE, quote = "", comment.char = "",
  stringsAsFactors = FALSE
)
stopifnot(
  identical(names(composition)[1:2], c("MAG", "Taxonomy")),
  nrow(composition) == 780L,
  ncol(composition) == 58L,
  nrow(metadata) == 56L,
  !anyDuplicated(composition$MAG),
  !anyDuplicated(metadata$Spring)
)

spring_ids <- names(composition)[-(1:2)]
metadata <- metadata[match(spring_ids, metadata$Spring), , drop = FALSE]
abundance <- t(as.matrix(composition[, spring_ids, drop = FALSE]))
storage.mode(abundance) <- "double"
rownames(abundance) <- spring_ids
colnames(abundance) <- composition$MAG
stopifnot(
  identical(rownames(abundance), metadata$Spring),
  all(is.finite(abundance)),
  min(abundance) >= 0,
  max(abs(rowSums(abundance) - 1)) < 1e-12,
  !anyNA(metadata$BroadRegion),
  !anyNA(metadata$MedianTemperatureC),
  !anyNA(metadata$MedianPH),
  length(unique(metadata$BroadRegion)) == 6L
)

prevalence <- colMeans(abundance > 0)
mean_abundance <- colMeans(abundance)
maximum_abundance <- apply(abundance, 2L, max)
filter_audit <- data.frame(
  MAG = composition$MAG,
  Taxonomy = composition$Taxonomy,
  Springs = nrow(abundance),
  NonzeroSprings = colSums(abundance > 0),
  Prevalence = prevalence,
  MeanRelativeAbundance = mean_abundance,
  MaximumRelativeAbundance = maximum_abundance,
  RelaxedFilter = prevalence >= 0.60 & mean_abundance >= 0.001,
  PrimaryFilter = prevalence >= 0.70 & mean_abundance >= 0.001,
  StrictFilter = prevalence >= 0.70 & mean_abundance >= 0.002,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(
  sum(filter_audit$RelaxedFilter) == 183L,
  sum(filter_audit$PrimaryFilter) == 93L,
  sum(filter_audit$StrictFilter) == 63L
)

analysis_contract <- data.frame(
  Item = c(
    "seed", "source_record", "inference_unit", "raw_features",
    "primary_feature_filter", "relaxed_feature_filter",
    "strict_feature_filter", "subcomposition_policy", "zero_policy",
    "transformation", "environment_adjustment", "residual_standardization",
    "primary_estimator", "regularization_path", "primary_selection",
    "stars_subsample_ratio", "edge_weight", "edge_bootstrap",
    "bootstrap_consensus", "high_stability_edge", "module_detection",
    "zi_pi_definition", "topology_candidate_gate", "topology_null",
    "deletion_robustness", "sensitivity_branches", "interpretation_boundary"
  ),
  Value = c(
    "20260729",
    "Figshare 30284068 v2; Scientific Data 2026; DOI 10.1038/s41597-026-07139-w",
    "one hot spring; n=56; every spring has equal inferential weight",
    "780 recovered MAGs in the catalog-closed spring profile",
    "prevalence >=0.70 and mean catalog-relative abundance >=0.001; 93 MAGs",
    "prevalence >=0.60 and mean catalog-relative abundance >=0.001; 183 MAGs",
    "prevalence >=0.70 and mean catalog-relative abundance >=0.002; 63 MAGs",
    "subset selected MAGs and reclose each spring to one before zero replacement",
    "add fixed fraction 0.000001 to every selected component before CLR",
    "natural-log centered log ratio (CLR)",
    "feature-wise OLS residuals from BroadRegion + z(MedianTemperatureC) + z(MedianPH)",
    "center and scale every residual feature; residual design rank=8 and df=48",
    "huge 1.5.1 graphical lasso sparse inverse covariance",
    "30 lambdas; lambda.min.ratio=0.05",
    "StARS instability threshold=0.05; 100 subsamples",
    "0.80 because n<=144 in huge 1.5.1",
    "signed partial correlation from the selected precision matrix",
    "1000 BroadRegion-stratified spring bootstraps; fixed nodes and primary lambda; refit adjustment",
    "selection frequency >=0.70 is the consensus audit threshold",
    "selection frequency >=0.80",
    "Louvain on the unweighted primary non-isolate graph; seed fixed",
    "classic unweighted within-module degree Zi and participation Pi; hub Zi>2.5; connector Pi>0.62",
    "stable degree >=2 and either Zi>2.5, Pi>0.62, or primary degree in the top decile",
    "1000 degree-preserving rewires; empirical upper-tail P for modularity and transitivity",
    "adaptive highest-degree removal versus 1000 random deletion orders; LCC normalized to initial nodes",
    "pseudocount 1e-7/1e-5; no environment adjustment; MB-OR; relaxed/strict filters; StARS 0.025/0.10",
    "undirected conditional associations in catalog-relative cross-sectional data; not direct, causal, or validated ecological interactions"
  ),
  Role = c(
    rep("Primary", 17L),
    "Uncertainty", "Audit", "Audit", "Topology", "Topology", "Audit",
    "Audit", "Audit", "Sensitivity", "Boundary"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(nrow(analysis_contract) == 27L, length(analysis_contract$Role) == 27L)

composition_out <- file.path(output_dir, "spring-mag-relative-abundance.tsv.gz")
metadata_out <- file.path(output_dir, "spring-metadata.tsv")
filter_out <- file.path(output_dir, "feature-filter-audit.tsv")
contract_out <- file.path(output_dir, "analysis-contract.tsv")
manifest_out <- file.path(output_dir, "resource-manifest.tsv")
checksum_out <- file.path(output_dir, "file-checksums.sha256")

if (!file.copy(source_paths[["composition"]], composition_out, overwrite = TRUE)) {
  stop("Could not copy the frozen spring composition.", call. = FALSE)
}
if (!file.copy(source_paths[["metadata"]], metadata_out, overwrite = TRUE)) {
  stop("Could not copy the frozen spring metadata.", call. = FALSE)
}
write_tsv(filter_audit, filter_out)
write_tsv(analysis_contract, contract_out)

resource_manifest <- data.frame(
  Resource = c(
    "Equal-sample spring MAG composition", "Spring metadata",
    "Feature-filter audit", "Analysis contract", "Anchor publication",
    "Network-method publication"
  ),
  Release = c(
    "Figshare 30284068 v2 / Article 23 deterministic aggregation",
    "Figshare 30284068 v2 / Article 23 deterministic aggregation",
    "Article 29 deterministic derivation", "Article 29 preregistered contract",
    "Scientific Data 13 (2026)", "PLOS Computational Biology 11 (2015)"
  ),
  Source = c(
    "https://doi.org/10.6084/m9.figshare.30284068.v2",
    "https://doi.org/10.6084/m9.figshare.30284068.v2",
    "spring-level composition", "analysis specification",
    "https://doi.org/10.1038/s41597-026-07139-w",
    "https://doi.org/10.1371/journal.pcbi.1004226"
  ),
  File = c(
    basename(composition_out), basename(metadata_out), basename(filter_out),
    basename(contract_out), NA_character_, NA_character_
  ),
  Bytes = c(
    as.numeric(file.info(composition_out)$size),
    as.numeric(file.info(metadata_out)$size),
    as.numeric(file.info(filter_out)$size),
    as.numeric(file.info(contract_out)$size),
    NA_real_, NA_real_
  ),
  SHA256 = c(
    sha256_file(composition_out), sha256_file(metadata_out),
    sha256_file(filter_out), sha256_file(contract_out),
    NA_character_, NA_character_
  ),
  Rows = c(780L, 56L, 780L, nrow(analysis_contract), NA_integer_, NA_integer_),
  Units = c(56L, 56L, 56L, 56L, NA_integer_, NA_integer_),
  InterpretationBoundary = c(
    "closed within the 780 recovered-MAG catalog; not whole-community abundance",
    "one row per hot spring; median environment and equal spring weight",
    "filter is outcome-free but still changes the estimand",
    "must be changed before rerunning, not after inspecting network appeal",
    "sampling was designed to capture local diversity, not probability sampling",
    "method motivates CLR, sparse graphical models, and StARS; no numeric replication is claimed"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(resource_manifest, manifest_out)

writeLines(
  c(
    "Article 29 frozen data notice",
    "",
    "Source: Korchagina et al. 2026, Scientific Data; DOI 10.1038/s41597-026-07139-w.",
    "Publisher payload: Figshare record 30284068 v2 under CC BY 4.0.",
    "The 500 local metagenomes were aggregated into 56 equal-weight hot-spring profiles before this chapter.",
    "The abundance table contains 780 recovered MAGs and is closed within that catalog.",
    "Mean sample recruitment into the catalog is approximately 13%, so values are not whole-community fractions.",
    "The primary network keeps 93 MAGs with prevalence >=70% and mean catalog-relative abundance >=0.1%.",
    "The analysis is cross-sectional and observational. Edges are conditional associations, not validated interactions.",
    "Topology-priority nodes and Zi-Pi roles are hypothesis generators, not ecological keystone proof.",
    "Routine analysis and rendering use checksum-locked compact tables and do not access the network."
  ),
  con = notice_path,
  useBytes = TRUE
)

payloads <- c(composition_out, metadata_out, filter_out, contract_out, manifest_out)
writeLines(
  paste(vapply(payloads, sha256_file, character(1L)), basename(payloads)),
  checksum_out,
  useBytes = TRUE
)

cat("Prepared Article 29 network inputs\n")
cat("  inference units:", nrow(abundance), "hot springs\n")
cat("  raw MAG features:", ncol(abundance), "\n")
cat("  relaxed / primary / strict:",
    sum(filter_audit$RelaxedFilter), "/",
    sum(filter_audit$PrimaryFilter), "/",
    sum(filter_audit$StrictFilter), "\n")
cat("  checksum payloads:", length(payloads), "\n")
