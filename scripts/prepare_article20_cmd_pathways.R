#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
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
required <- c("abundance-rda", "coverage-rda", "output-dir")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

if (!requireNamespace("curatedMetagenomicData", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE)) {
  stop("curatedMetagenomicData and digest are required.", call. = FALSE)
}

abundance_rda <- normalizePath(args[["abundance-rda"]], mustWork = TRUE)
coverage_rda <- normalizePath(args[["coverage-rda"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected <- c(
  abundance = "ead7c78c075fec92a7d641b731594e068b2ba2a47479151d081c338f615af121",
  coverage = "73a1b77b70f88e9028e8707ba3e99b93f0ff99cd91401a3966c4f7a31dbfc3a1"
)
sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
observed <- c(
  abundance = sha256_file(abundance_rda),
  coverage = sha256_file(coverage_rda)
)
stopifnot(identical(unname(observed), unname(expected)))

load_one <- function(path, object_name) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  stopifnot(object_name %in% loaded)
  get(object_name, envir = env, inherits = FALSE)
}

abundance_name <- "2021-10-14.AsnicarF_2017.pathway_abundance"
coverage_name <- "2021-10-14.AsnicarF_2017.pathway_coverage"
abundance <- load_one(abundance_rda, abundance_name)
coverage <- load_one(coverage_rda, coverage_name)

stopifnot(
  is.matrix(abundance),
  is.matrix(coverage),
  identical(dim(abundance), c(11173L, 24L)),
  identical(dim(coverage), c(11173L, 24L)),
  setequal(rownames(abundance), rownames(coverage)),
  identical(colnames(abundance), colnames(coverage)),
  !anyNA(abundance),
  !anyNA(coverage),
  min(abundance) >= 0,
  min(coverage) >= 0,
  max(coverage) <= 1
)
coverage <- coverage[rownames(abundance), colnames(abundance), drop = FALSE]
stopifnot(identical(rownames(abundance), rownames(coverage)))

metadata_all <- curatedMetagenomicData::sampleMetadata
metadata <- metadata_all[
  match(colnames(abundance), metadata_all$sample_id),
  c(
    "study_name", "sample_id", "subject_id", "body_site", "family_role",
    "visit_number", "days_from_first_collection", "NCBI_accession", "PMID"
  ),
  drop = FALSE
]
stopifnot(
  identical(metadata$sample_id, colnames(abundance)),
  all(metadata$study_name == "AsnicarF_2017"),
  length(unique(metadata$subject_id)) == 15L
)

write_matrix <- function(x, path, feature_name) {
  con <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(con), add = TRUE)
  out <- data.frame(feature = rownames(x), x, check.names = FALSE)
  names(out)[[1L]] <- feature_name
  utils::write.table(
    out, file = con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}

abundance_tsv <- file.path(output_dir, "pathway-abundance.tsv.gz")
coverage_tsv <- file.path(output_dir, "pathway-coverage.tsv.gz")
metadata_tsv <- file.path(output_dir, "sample-metadata.tsv")
manifest_tsv <- file.path(output_dir, "resource-manifest.tsv")

write_matrix(abundance, abundance_tsv, "Pathway")
write_matrix(coverage, coverage_tsv, "Pathway")
utils::write.table(
  metadata, metadata_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

manifest <- data.frame(
  resource = c("pathway_abundance", "pathway_coverage"),
  experimenthub_id = c("EH7089", "EH7090"),
  selected_title = c(abundance_name, coverage_name),
  source_url = c(
    "https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/curatedMetagenomicData/2021-10-14/AsnicarF_2017/2021-10-14.AsnicarF_2017.pathway_abundance.rda",
    "https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/curatedMetagenomicData/2021-10-14/AsnicarF_2017/2021-10-14.AsnicarF_2017.pathway_coverage.rda"
  ),
  raw_file = c(
    "AsnicarF_2017.pathway_abundance.rda",
    "AsnicarF_2017.pathway_coverage.rda"
  ),
  raw_bytes = c(file.info(abundance_rda)$size, file.info(coverage_rda)$size),
  raw_sha256 = unname(observed),
  derived_file = c("pathway-abundance.tsv.gz", "pathway-coverage.tsv.gz"),
  derived_bytes = c(file.info(abundance_tsv)$size, file.info(coverage_tsv)$size),
  derived_sha256 = c(sha256_file(abundance_tsv), sha256_file(coverage_tsv)),
  rows = c(nrow(abundance), nrow(coverage)),
  samples = c(ncol(abundance), ncol(coverage)),
  value_min = c(min(abundance), min(coverage)),
  value_max = c(max(abundance), max(coverage)),
  package_version = as.character(utils::packageVersion("curatedMetagenomicData")),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
utils::write.table(
  manifest, manifest_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

notice <- c(
  "Article 20 data notice",
  "",
  "The two .rda files are the curatedMetagenomicData 3.12.0 ExperimentHub resources",
  "EH7089 (AsnicarF_2017 pathway_abundance) and EH7090 (pathway_coverage).",
  "They contain HUMAnN3 functional profiles for 24 longitudinal samples from 15 subjects.",
  "The compressed TSV files are lossless tabular exports generated by",
  "scripts/prepare_article20_cmd_pathways.R; sample-metadata.tsv is aligned by sample_id.",
  "",
  "Primary study: Asnicar et al. 2017, mSystems 2:e00164-16.",
  "DOI: 10.1128/mSystems.00164-16; PMID: 28144631.",
  "Resource: curatedMetagenomicData; DOI: 10.1038/nmeth.4468.",
  "No biological group-comparison test is authorized by this tutorial branch."
)
writeLines(notice, file.path(output_dir, "data-NOTICE.txt"), useBytes = TRUE)

file.copy(abundance_rda, file.path(output_dir, "AsnicarF_2017.pathway_abundance.rda"), overwrite = TRUE)
file.copy(coverage_rda, file.path(output_dir, "AsnicarF_2017.pathway_coverage.rda"), overwrite = TRUE)

payloads <- sort(setdiff(
  list.files(output_dir, recursive = FALSE, full.names = TRUE),
  file.path(output_dir, "file-checksums.sha256")
))
checksum_lines <- vapply(
  payloads,
  function(path) paste(sha256_file(path), basename(path)),
  character(1L)
)
writeLines(checksum_lines, file.path(output_dir, "file-checksums.sha256"), useBytes = TRUE)

cat("Prepared Article 20 curatedMetagenomicData pathway resources\n")
cat("Rows: ", nrow(abundance), "\n", sep = "")
cat("Samples: ", ncol(abundance), "\n", sep = "")
cat("Subjects: ", length(unique(metadata$subject_id)), "\n", sep = "")
