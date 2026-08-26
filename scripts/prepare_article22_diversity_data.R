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
required <- c(
  "species-rds", "gene-rda", "mag-biom", "mag-metadata",
  "mag-recruitment", "output-dir", "notice"
)
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c(
  "Matrix", "TreeSummarizedExperiment", "SummarizedExperiment",
  "curatedMetagenomicData", "jsonlite", "digest"
)
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

paths <- lapply(args[required[1:5]], normalizePath, mustWork = TRUE)
names(paths) <- required[1:5]
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
notice_path <- normalizePath(args[["notice"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(notice_path), recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

md5_file <- function(path) {
  digest::digest(file = path, algo = "md5", serialize = FALSE)
}

expected_sources <- data.frame(
  key = c("species-rds", "gene-rda", "mag-biom", "mag-metadata", "mag-recruitment"),
  bytes = c(199904L, 55739524L, 2017689L, 133160L, 9059L),
  sha256 = c(
    "2952774730bff2af9e13c9c40058320aed524dc2dc7408b44dc6e697c06564b2",
    "f171e511667a5d265c529d731b81407af22da3be27db2ba8587cbc6bc159257a",
    "d3188968399136dfe595bacf9b2fda01568284a03945550f7ae0153eb4a7b131",
    "8dc11ab08a4c53e800038f3ab87958b085267417aa7ced06ce0b78a9985f1c37",
    "832f801eab0f6fbfcc859ad53787751856cf3dd9dc3343c179ae315bb37c40a8"
  ),
  md5 = c(
    NA_character_,
    "49caa4e88cbd51e7f5700cbcf4590e55",
    "0c1762f31a5c4473dec78571a7b74287",
    "402095b04e2c5518cbec462f41528d4f",
    "51288407e9927f23a25b210edcda7b47"
  ),
  stringsAsFactors = FALSE
)

source_audit <- expected_sources
source_audit$observed_bytes <- vapply(
  source_audit$key, function(key) as.numeric(file.info(paths[[key]])$size), numeric(1L)
)
source_audit$observed_sha256 <- vapply(
  source_audit$key, function(key) sha256_file(paths[[key]]), character(1L)
)
source_audit$observed_md5 <- vapply(
  source_audit$key, function(key) md5_file(paths[[key]]), character(1L)
)
stopifnot(
  identical(source_audit$observed_bytes, as.numeric(source_audit$bytes)),
  identical(source_audit$observed_sha256, source_audit$sha256),
  all(is.na(source_audit$md5) | source_audit$observed_md5 == source_audit$md5)
)

write_wide_gz <- function(x, path, annotations, chunk_size = 25000L) {
  stopifnot(nrow(x) == nrow(annotations), ncol(x) > 0L)
  con <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(con), add = TRUE)
  writeLines(
    paste(c(names(annotations), colnames(x)), collapse = "\t"),
    con = con,
    useBytes = TRUE
  )
  starts <- seq.int(1L, nrow(x), by = chunk_size)
  for (first in starts) {
    last <- min(nrow(x), first + chunk_size - 1L)
    block <- data.frame(
      annotations[first:last, , drop = FALSE],
      as.matrix(x[first:last, , drop = FALSE]),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    utils::write.table(
      block,
      file = con,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      na = ""
    )
  }
  invisible(path)
}

# Human species profiles: preserve the original MetaPhlAn percentage scale.
species_object <- readRDS(paths[["species-rds"]])
stopifnot(
  inherits(species_object, "TreeSummarizedExperiment"),
  identical(dim(species_object), c(298L, 24L)),
  identical(SummarizedExperiment::assayNames(species_object), "relative_abundance")
)
species_matrix <- SummarizedExperiment::assay(
  species_object, "relative_abundance", withDimnames = TRUE
)
stopifnot(
  all(grepl("|s__", rownames(species_matrix), fixed = TRUE)),
  !anyNA(species_matrix),
  min(species_matrix) >= 0,
  min(colSums(species_matrix)) > 90,
  max(colSums(species_matrix)) < 101
)
species_label <- sub("^.*\\|s__", "", rownames(species_matrix))
species_label <- gsub("_", " ", species_label, fixed = TRUE)

# Human gene-family profiles: community-level rows only, no special rows.
suppressPackageStartupMessages(library(Matrix))
gene_env <- new.env(parent = emptyenv())
gene_object_name <- "2021-10-14.AsnicarF_2017.gene_families"
loaded <- load(paths[["gene-rda"]], envir = gene_env)
stopifnot(gene_object_name %in% loaded)
gene_matrix_raw <- get(gene_object_name, envir = gene_env, inherits = FALSE)
stopifnot(
  inherits(gene_matrix_raw, "dgCMatrix"),
  identical(dim(gene_matrix_raw), c(2704846L, 24L)),
  identical(colnames(gene_matrix_raw), colnames(species_matrix)),
  !anyNA(gene_matrix_raw@x),
  min(gene_matrix_raw@x) > 0
)
special_features <- c("UNMAPPED", "UNGROUPED", "UNINTEGRATED")
gene_unstratified <- !grepl("|", rownames(gene_matrix_raw), fixed = TRUE)
gene_ordinary <- !rownames(gene_matrix_raw) %in% special_features
gene_matrix <- gene_matrix_raw[gene_unstratified & gene_ordinary, , drop = FALSE]
gene_prevalence_count <- Matrix::rowSums(gene_matrix > 0)
threshold_counts <- c(`0.10` = 3L, `0.20` = 5L, `0.50` = 12L)
gene_keep10 <- gene_prevalence_count >= threshold_counts[["0.10"]]
gene_matrix10 <- gene_matrix[gene_keep10, , drop = FALSE]
stopifnot(
  identical(dim(gene_matrix10), c(415581L, 24L)),
  Matrix::nnzero(gene_matrix10) == 1963233L
)

# Align richer sample metadata from the same curatedMetagenomicData release.
metadata_all <- curatedMetagenomicData::sampleMetadata
human_fields <- c(
  "study_name", "sample_id", "subject_id", "family_role", "visit_number",
  "days_from_first_collection", "age_category", "gender", "NCBI_accession",
  "PMID", "number_reads"
)
human_metadata <- metadata_all[
  match(colnames(species_matrix), metadata_all$sample_id),
  human_fields,
  drop = FALSE
]
stopifnot(
  identical(human_metadata$sample_id, colnames(species_matrix)),
  all(human_metadata$study_name == "AsnicarF_2017"),
  length(unique(human_metadata$subject_id)) == 15L
)

# Environmental MAG catalog: parse BIOM 1.0 sparse coordinates without BIOM I/O drift.
mag_biom <- jsonlite::fromJSON(paths[["mag-biom"]], simplifyVector = FALSE)
stopifnot(
  identical(mag_biom$format, "Biological Observation Matrix 1.0.0"),
  identical(mag_biom$matrix_type, "sparse"),
  identical(as.integer(unlist(mag_biom$shape)), c(780L, 500L))
)
triplet <- do.call(
  rbind,
  lapply(mag_biom$data, function(entry) unlist(entry, use.names = FALSE))
)
mag_ids <- vapply(mag_biom$rows, `[[`, character(1L), "id")
mag_sample_ids <- vapply(mag_biom$columns, `[[`, character(1L), "id")
mag_taxonomy <- vapply(
  mag_biom$rows,
  function(entry) paste(unlist(entry$metadata$taxonomy), collapse = ";"),
  character(1L)
)
mag_matrix <- Matrix::sparseMatrix(
  i = as.integer(triplet[, 1L]) + 1L,
  j = as.integer(triplet[, 2L]) + 1L,
  x = as.numeric(triplet[, 3L]),
  dims = c(780L, 500L),
  dimnames = list(mag_ids, mag_sample_ids)
)
stopifnot(
  Matrix::nnzero(mag_matrix) == 86510L,
  min(Matrix::colSums(mag_matrix)) > 0.99998,
  max(Matrix::colSums(mag_matrix)) < 1.00002
)

mag_metadata_raw <- utils::read.delim(
  paths[["mag-metadata"]], check.names = FALSE, quote = "", comment.char = ""
)
mag_metadata <- mag_metadata_raw[
  match(mag_sample_ids, mag_metadata_raw$sample),
  c(
    "sample", "hotspring", "hotspring_common_name", "region",
    "broad_region_short", "latitude", "longitude", "date", "year", "pH",
    "temperature", "conductivity", "temperature_regime", "type",
    "sequencing_platform", "metagenome_SRR", "metagenome_SAMN",
    "raw_Nsequences", "filtered_Nsequences", "metagenome2MAG_read_recruitment_rate"
  ),
  drop = FALSE
]
mag_recruitment <- utils::read.delim(
  paths[["mag-recruitment"]], check.names = FALSE, quote = "", comment.char = "#"
)
mag_recruitment <- mag_recruitment[
  match(mag_sample_ids, mag_recruitment$metagenome), , drop = FALSE
]
names(mag_recruitment) <- c("sample_id", "total_hit_rate")
stopifnot(
  identical(mag_metadata$sample, mag_sample_ids),
  identical(mag_recruitment$sample_id, mag_sample_ids),
  length(unique(mag_metadata$hotspring)) == 56L,
  max(abs(
    mag_metadata$metagenome2MAG_read_recruitment_rate - mag_recruitment$total_hit_rate
  )) < 1e-12,
  abs(mean(mag_recruitment$total_hit_rate) - 0.12983) < 1e-4
)

# Frozen, network-free routine inputs.
species_tsv <- file.path(output_dir, "species-relative-abundance.tsv.gz")
gene_tsv <- file.path(output_dir, "gene-family-prevalence10.tsv.gz")
human_metadata_tsv <- file.path(output_dir, "human-sample-metadata.tsv")
feature_audit_tsv <- file.path(output_dir, "human-feature-audit.tsv")
mag_tsv <- file.path(output_dir, "mag-relative-abundance.tsv.gz")
mag_metadata_tsv <- file.path(output_dir, "hot-spring-sample-metadata.tsv")
mag_recruitment_tsv <- file.path(output_dir, "mag-recruitment.tsv")
filter_contract_tsv <- file.path(output_dir, "filter-contract.tsv")
resource_manifest_tsv <- file.path(output_dir, "resource-manifest.tsv")

write_wide_gz(
  species_matrix,
  species_tsv,
  data.frame(
    Feature = rownames(species_matrix),
    Species = species_label,
    stringsAsFactors = FALSE
  ),
  chunk_size = 1000L
)
write_wide_gz(
  gene_matrix10,
  gene_tsv,
  data.frame(GeneFamily = rownames(gene_matrix10), stringsAsFactors = FALSE),
  chunk_size = 20000L
)
write_wide_gz(
  mag_matrix,
  mag_tsv,
  data.frame(MAG = mag_ids, Taxonomy = mag_taxonomy, stringsAsFactors = FALSE),
  chunk_size = 1000L
)
utils::write.table(
  human_metadata, human_metadata_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)
utils::write.table(
  mag_metadata, mag_metadata_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)
utils::write.table(
  mag_recruitment, mag_recruitment_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

summarize_stage <- function(feature_space, stage, x, prevalence_threshold = NA_real_) {
  sums <- Matrix::colSums(x)
  data.frame(
    FeatureSpace = feature_space,
    Stage = stage,
    PrevalenceThreshold = prevalence_threshold,
    Features = nrow(x),
    Samples = ncol(x),
    NonzeroCells = Matrix::nnzero(x),
    ColumnSumMin = min(sums),
    ColumnSumMedian = stats::median(sums),
    ColumnSumMax = max(sums),
    stringsAsFactors = FALSE
  )
}
feature_audit <- do.call(
  rbind,
  c(
    list(
      summarize_stage("Species", "MetaPhlAn species percentages", species_matrix),
      summarize_stage("Gene family", "Raw community plus stratified", gene_matrix_raw),
      summarize_stage("Gene family", "Unstratified ordinary", gene_matrix)
    ),
    lapply(names(threshold_counts), function(label) {
      count <- threshold_counts[[label]]
      summarize_stage(
        "Gene family",
        paste0("Prevalence >= ", label),
        gene_matrix[gene_prevalence_count >= count, , drop = FALSE],
        as.numeric(label)
      )
    })
  )
)
utils::write.table(
  feature_audit, feature_audit_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

filter_contract <- data.frame(
  FeatureSpace = c("Species", "Gene family", "MAG catalog"),
  Source = c(
    "curatedMetagenomicData EH7091",
    "curatedMetagenomicData EH7086",
    "Figshare 30284068 v2 file 61153444"
  ),
  NativeScale = c(
    "MetaPhlAn relative abundance (%)",
    "HUMAnN gene-family relative abundance (fraction)",
    "Genome-coverage relative abundance (catalog closure)"
  ),
  FeatureUniverse = c(
    "298 species-level MetaPhlAn features",
    "Unstratified ordinary UniRef90 families",
    "780 recovered MAGs"
  ),
  PrimaryFilter = c("None", "Sample prevalence >= 20% (5/24)", "None"),
  Sensitivity = c(
    "Aitchison pseudocount 1e-06 versus 1e-05",
    "Prevalence 10%, 20%, 50%; pseudocount 1e-08 versus 1e-07",
    "Reads recruited to catalog reported separately"
  ),
  Closure = c(
    "Within observed species per sample",
    "After each prevalence filter per sample",
    "Publisher-provided catalog closure; rechecked"
  ),
  Alpha = "Hill q=0, q=1, q=2",
  Beta = c(
    "Bray-Curtis; binary Jaccard; Aitchison sensitivity",
    "Bray-Curtis; binary Jaccard; Aitchison sensitivity",
    "Bray-Curtis; binary Jaccard"
  ),
  Inference = "Descriptive only; no group test or permutation",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
utils::write.table(
  filter_contract, filter_contract_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

derived_files <- c(
  species_tsv, gene_tsv, human_metadata_tsv, feature_audit_tsv,
  mag_tsv, mag_metadata_tsv, mag_recruitment_tsv, filter_contract_tsv
)
resource_manifest <- data.frame(
  Resource = c(
    "Human species", "Human gene families", "Hot-spring MAG abundance",
    "Hot-spring metadata", "Hot-spring recruitment"
  ),
  RepositoryID = c("EH7091", "EH7086", "61153444", "61153429", "61153471"),
  Release = c(
    "curatedMetagenomicData 2021-10-14",
    "curatedMetagenomicData 2021-10-14",
    rep("Figshare 30284068 v2", 3L)
  ),
  SourceURL = c(
    "https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/curatedMetagenomicData/2021-10-14/AsnicarF_2017/2021-10-14.AsnicarF_2017.relative_abundance.rda",
    "https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/curatedMetagenomicData/2021-10-14/AsnicarF_2017/2021-10-14.AsnicarF_2017.gene_families.rda",
    "https://ndownloader.figshare.com/files/61153444",
    "https://ndownloader.figshare.com/files/61153429",
    "https://ndownloader.figshare.com/files/61153471"
  ),
  RawBytes = source_audit$observed_bytes,
  RawMD5 = source_audit$observed_md5,
  RawSHA256 = source_audit$observed_sha256,
  DerivedFile = c(
    basename(species_tsv), basename(gene_tsv), basename(mag_tsv),
    basename(mag_metadata_tsv), basename(mag_recruitment_tsv)
  ),
  DerivedBytes = vapply(
    c(species_tsv, gene_tsv, mag_tsv, mag_metadata_tsv, mag_recruitment_tsv),
    function(path) as.numeric(file.info(path)$size), numeric(1L)
  ),
  DerivedSHA256 = vapply(
    c(species_tsv, gene_tsv, mag_tsv, mag_metadata_tsv, mag_recruitment_tsv),
    sha256_file,
    character(1L)
  ),
  Rows = c(298L, nrow(gene_matrix10), 780L, nrow(mag_metadata), nrow(mag_recruitment)),
  Samples = c(24L, 24L, 500L, 500L, 500L),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
utils::write.table(
  resource_manifest, resource_manifest_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = ""
)

notice <- c(
  "Article 22 data notice",
  "",
  "Human branch:",
  "- Species profiles are the checksum-locked curatedMetagenomicData 3.12.0",
  "  AsnicarF_2017.relative_abundance resource (ExperimentHub EH7091).",
  "- Gene families are the matching 2021-10-14 AsnicarF_2017.gene_families",
  "  resource (ExperimentHub EH7086). Community-level UniRef90 rows were retained;",
  "  taxon-stratified rows and UNMAPPED/UNGROUPED/UNINTEGRATED were excluded.",
  "- The frozen gene-family table retains every feature present in at least 3/24",
  "  profiles (10% threshold); routine analysis derives the 20% primary and 50%",
  "  sensitivity branches from that table.",
  "- Primary study: Asnicar et al. 2017, mSystems 2:e00164-16.",
  "  DOI: 10.1128/mSystems.00164-16; PMID: 28144631.",
  "",
  "MAG branch:",
  "- The BIOM table, sample metadata, and read-recruitment table are publisher",
  "  files 61153444, 61153429, and 61153471 from Figshare record 30284068 v2",
  "  under CC BY 4.0.",
  "- Primary study: Korchagina et al. 2026, Scientific Data.",
  "  DOI: 10.1038/s41597-026-07139-w.",
  "- MAG abundance is closed within the 780 recovered-genome catalog. Total read",
  "  recruitment (mean 0.12983) remains separate and is not a whole-community sum.",
  "",
  "The human and hot-spring branches are independent datasets. They must not be",
  "treated as one matched cohort or entered into a shared biological group test.",
  "Routine QA is network-free and does not read the original RDA or BIOM files."
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

cat("Prepared Article 22 diversity inputs\n")
cat("Human profiles: ", ncol(species_matrix), "\n", sep = "")
cat("Human subjects: ", length(unique(human_metadata$subject_id)), "\n", sep = "")
cat("Species features: ", nrow(species_matrix), "\n", sep = "")
cat("Gene families at prevalence >=10%: ", nrow(gene_matrix10), "\n", sep = "")
cat("Gene families at prevalence >=20%: ", sum(gene_prevalence_count >= 5L), "\n", sep = "")
cat("Gene families at prevalence >=50%: ", sum(gene_prevalence_count >= 12L), "\n", sep = "")
cat("MAG samples: ", ncol(mag_matrix), "\n", sep = "")
cat("Hot springs: ", length(unique(mag_metadata$hotspring)), "\n", sep = "")
cat("MAG features: ", nrow(mag_matrix), "\n", sep = "")
cat("Mean read recruitment: ", sprintf("%.5f", mean(mag_recruitment$total_hit_rate)), "\n", sep = "")
