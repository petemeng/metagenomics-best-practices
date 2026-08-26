#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "UTC")

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c(
  "project-root", "cache-dir", "snapshot", "manifest", "log"
)
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0L) {
  stop(
    "Missing required arguments: ",
    paste(paste0("--", missing_args), collapse = ", "),
    call. = FALSE
  )
}

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
cache_dir <- normalizePath(args[["cache-dir"]], mustWork = FALSE)
snapshot_path <- normalizePath(args[["snapshot"]], mustWork = FALSE)
manifest_path <- normalizePath(args[["manifest"]], mustWork = FALSE)
log_path <- normalizePath(args[["log"]], mustWork = FALSE)

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snapshot_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "curatedMetagenomicData", "ExperimentHub", "BiocFileCache",
  "SummarizedExperiment", "TreeSummarizedExperiment", "S4Vectors", "digest"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

sink(log_path, split = TRUE)
on.exit({
  while (sink.number() > 0L) {
    sink()
  }
}, add = TRUE)

cat("Article 12 curatedMetagenomicData one-time retrieval\n")
cat("UTC start: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
cat("Project root: <PROJECT_ROOT>\n")
cat("Cache: <PROJECT_ROOT>/.cache/R/ExperimentHub\n")
cat("R: ", paste(R.version$major, R.version$minor, sep = "."), "\n", sep = "")
cat(
  "Bioconductor: ",
  as.character(BiocManager::version()),
  "\n",
  sep = ""
)
cat(
  "curatedMetagenomicData: ",
  as.character(utils::packageVersion("curatedMetagenomicData")),
  "\n",
  sep = ""
)

stopifnot(
  identical(paste(R.version$major, R.version$minor, sep = "."), "4.4.1"),
  identical(as.character(BiocManager::version()), "3.19"),
  identical(
    as.character(utils::packageVersion("curatedMetagenomicData")),
    "3.12.0"
  )
)

resource_pattern <- "AsnicarF_2017.relative_abundance"
selected_title <- "2021-10-14.AsnicarF_2017.relative_abundance"

ExperimentHub::setExperimentHubOption("CACHE", cache_dir)
ExperimentHub::setExperimentHubOption("LOCAL", FALSE)

dryrun_titles <- suppressMessages(
  curatedMetagenomicData::curatedMetagenomicData(
    resource_pattern,
    dryrun = TRUE
  )
)
stopifnot(selected_title %in% dryrun_titles)
cat("Dry-run candidates:\n")
cat(paste0("  ", dryrun_titles, "\n"))

hub <- ExperimentHub::ExperimentHub(localHub = FALSE)
hub_hit <- AnnotationHub::query(hub, selected_title)
hub_hit <- hub_hit[
  ExperimentHub::package(hub_hit) %in% "curatedMetagenomicData"
]
stopifnot(length(hub_hit) == 1L)

experimenthub_id <- names(hub_hit)[[1L]]
hub_metadata <- as.data.frame(S4Vectors::mcols(hub_hit))
cat("Selected ExperimentHub ID: ", experimenthub_id, "\n", sep = "")

online_list <- suppressMessages(
  curatedMetagenomicData::curatedMetagenomicData(
    resource_pattern,
    dryrun = FALSE,
    counts = FALSE,
    rownames = "long"
  )
)
stopifnot(
  length(online_list) == 1L,
  identical(names(online_list), selected_title)
)
online_object <- online_list[[1L]]

ExperimentHub::setExperimentHubOption("LOCAL", TRUE)
offline_list <- suppressMessages(
  curatedMetagenomicData::curatedMetagenomicData(
    resource_pattern,
    dryrun = FALSE,
    counts = FALSE,
    rownames = "long"
  )
)
stopifnot(
  length(offline_list) == 1L,
  identical(names(offline_list), selected_title)
)
offline_object <- offline_list[[1L]]
offline_replay_identical <- identical(online_object, offline_object)
stopifnot(offline_replay_identical)

stopifnot(
  methods::is(offline_object, "TreeSummarizedExperiment"),
  length(SummarizedExperiment::assays(offline_object)) == 1L,
  identical(
    names(SummarizedExperiment::assays(offline_object)),
    "relative_abundance"
  )
)

saveRDS(offline_object, snapshot_path, version = 3L, compress = "xz")

bfc <- BiocFileCache::BiocFileCache(cache_dir, ask = FALSE)
cache_info <- as.data.frame(BiocFileCache::bfcinfo(bfc))
cache_rows <- cache_info[
  startsWith(cache_info[["rname"]], paste0(experimenthub_id, " : ")),
  ,
  drop = FALSE
]
stopifnot(nrow(cache_rows) == 1L)
resource_cache_path <- cache_rows[["rpath"]][[1L]]
stopifnot(file.exists(resource_cache_path))

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

assay_matrix <- SummarizedExperiment::assay(
  offline_object,
  "relative_abundance"
)
column_sums <- colSums(assay_matrix)
metadata <- as.data.frame(SummarizedExperiment::colData(offline_object))

metadata_field <- function(name) {
  if (name %in% names(hub_metadata)) {
    as.character(hub_metadata[[name]][[1L]])
  } else {
    NA_character_
  }
}

manifest <- data.frame(
  resource_pattern = resource_pattern,
  selected_title = selected_title,
  experimenthub_id = experimenthub_id,
  source_url = metadata_field("sourceurl"),
  rdatapath = metadata_field("rdatapath"),
  rdataclass = metadata_field("rdataclass"),
  resource_cache_bytes = file.info(resource_cache_path)$size,
  resource_cache_sha256 = sha256_file(resource_cache_path),
  snapshot_file = file.path(
    "data",
    "small",
    basename(snapshot_path)
  ),
  snapshot_bytes = file.info(snapshot_path)$size,
  snapshot_sha256 = sha256_file(snapshot_path),
  object_class = class(offline_object)[[1L]],
  assay_name = "relative_abundance",
  features = nrow(offline_object),
  samples = ncol(offline_object),
  metadata_columns = ncol(metadata),
  minimum_column_sum = min(column_sums),
  median_column_sum = stats::median(column_sums),
  maximum_column_sum = max(column_sums),
  offline_replay_identical = offline_replay_identical,
  retrieved_at_utc = format(
    Sys.time(),
    tz = "UTC",
    format = "%Y-%m-%dT%H:%M:%SZ"
  ),
  package_version = as.character(
    utils::packageVersion("curatedMetagenomicData")
  ),
  bioconductor_version = as.character(BiocManager::version()),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  check.names = FALSE
)

utils::write.table(
  manifest,
  file = manifest_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

cat("Object class: ", class(offline_object)[[1L]], "\n", sep = "")
cat(
  "Object dimensions: ",
  nrow(offline_object),
  " features x ",
  ncol(offline_object),
  " samples\n",
  sep = ""
)
cat(
  "Relative-abundance column sums: ",
  paste(signif(range(column_sums), 8L), collapse = " to "),
  "\n",
  sep = ""
)
cat(
  "Raw cache SHA-256: ",
  manifest$resource_cache_sha256,
  "\n",
  sep = ""
)
cat(
  "Snapshot SHA-256: ",
  manifest$snapshot_sha256,
  "\n",
  sep = ""
)
cat("Offline replay identical: TRUE\n")
cat("UTC end: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
