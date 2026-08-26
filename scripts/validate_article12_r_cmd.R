#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "UTC")
set.seed(20260720)

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
required_args <- c("project-root", "output-dir", "figure-dir")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0L) {
  stop(
    "Missing required arguments: ",
    paste(paste0("--", missing_args), collapse = ", "),
    call. = FALSE
  )
}

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_packages <- c(
  "jsonlite", "digest", "readr", "dplyr", "tidyr", "tibble", "ggplot2",
  "scales", "BiocManager", "SummarizedExperiment",
  "TreeSummarizedExperiment", "sessioninfo"
)
missing_bootstrap <- bootstrap_packages[
  !vapply(bootstrap_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bootstrap) > 0L) {
  stop(
    "Missing validator package(s): ",
    paste(missing_bootstrap, collapse = ", "),
    call. = FALSE
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) y else x
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

write_tsv <- function(x, filename) {
  readr::write_tsv(x, file.path(output_dir, filename), na = "")
}

sanitize_text <- function(x) {
  x <- gsub(project_root, "<PROJECT_ROOT>", x, fixed = TRUE)
  x <- gsub(output_dir, "<OUTPUT_DIR>", x, fixed = TRUE)
  x <- gsub(path.expand("~"), "<HOME>", x, fixed = TRUE)
  x
}

save_publication_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 350,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(stem, ".tiff")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 350,
    compression = "lzw",
    bg = "white"
  )
}

theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2",
        color = "#666666",
        linewidth = 0.3
      )
    )
}

small_dir <- file.path(project_root, "data", "small")
lock_path <- file.path(project_root, "env", "renv.lock")
package_contract_path <- file.path(small_dir, "12-package-contract.tsv")
resource_manifest_path <- file.path(
  small_dir,
  "12-cmd-resource-manifest.tsv"
)
snapshot_path <- file.path(
  small_dir,
  "12-cmd-asnicarf-2017-relative-abundance.rds"
)
retrieval_log_path <- file.path(small_dir, "12-resource-retrieval.log")
notice_path <- file.path(small_dir, "12-data-NOTICE.txt")

required_inputs <- c(
  lock_path,
  package_contract_path,
  resource_manifest_path,
  snapshot_path,
  retrieval_log_path,
  notice_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Article 12 input(s): ",
    paste(basename(missing_inputs), collapse = ", "),
    call. = FALSE
  )
}

package_contract <- readr::read_tsv(
  package_contract_path,
  show_col_types = FALSE,
  progress = FALSE
)
resource_manifest <- readr::read_tsv(
  resource_manifest_path,
  show_col_types = FALSE,
  progress = FALSE
)
lock <- jsonlite::read_json(lock_path, simplifyVector = FALSE)

if (nrow(resource_manifest) != 1L) {
  stop("Article 12 resource manifest must have exactly one row.", call. = FALSE)
}

library_scope <- function(package) {
  path <- normalizePath(find.package(package), mustWork = TRUE)
  if (startsWith(path, project_root)) {
    "Project library"
  } else if (startsWith(path, "/usr/") || startsWith(path, "/usr/local/")) {
    "System / site library"
  } else {
    "User library"
  }
}

observed_package <- lapply(package_contract$package, function(package) {
  loaded <- requireNamespace(package, quietly = TRUE)
  if (!loaded) {
    return(list(
      loaded = FALSE,
      version = NA_character_,
      repository = NA_character_,
      scope = NA_character_
    ))
  }
  description <- utils::packageDescription(package)
  list(
    loaded = TRUE,
    version = as.character(description$Version),
    repository = as.character(description$Repository %||% "unknown"),
    scope = library_scope(package)
  )
})

lock_record <- lapply(package_contract$package, function(package) {
  record <- lock$Packages[[package]]
  if (is.null(record)) {
    return(list(
      version = NA_character_,
      source = NA_character_,
      repository = NA_character_
    ))
  }
  list(
    version = as.character(record$Version %||% NA_character_),
    source = as.character(record$Source %||% NA_character_),
    repository = as.character(record$Repository %||% NA_character_)
  )
})

package_audit <- tibble::tibble(
  package = package_contract$package,
  role = package_contract$role,
  expected_version = package_contract$expected_version,
  observed_version = vapply(
    observed_package,
    `[[`,
    character(1),
    "version"
  ),
  expected_repository = package_contract$repository,
  observed_repository = vapply(
    observed_package,
    `[[`,
    character(1),
    "repository"
  ),
  library_scope = vapply(observed_package, `[[`, character(1), "scope"),
  namespace_loaded = vapply(observed_package, `[[`, logical(1), "loaded"),
  lock_version = vapply(lock_record, `[[`, character(1), "version"),
  lock_source = vapply(lock_record, `[[`, character(1), "source"),
  lock_repository = vapply(lock_record, `[[`, character(1), "repository")
) |>
  dplyr::mutate(
    version_status = ifelse(
      namespace_loaded & observed_version == expected_version,
      "PASS",
      "FAIL"
    ),
    repository_status = ifelse(
      observed_repository == expected_repository,
      "PASS",
      "FAIL"
    ),
    lock_status = ifelse(
      lock_version == expected_version &
        !is.na(lock_source) &
        nzchar(lock_source),
      "PASS",
      "FAIL"
    ),
    status = ifelse(
      version_status == "PASS" &
        repository_status == "PASS" &
        lock_status == "PASS",
      "PASS",
      "FAIL"
    )
  )

lock_audit <- package_audit |>
  dplyr::select(
    package,
    expected_version,
    lock_version,
    lock_source,
    lock_repository,
    lock_status
  )

object <- readRDS(snapshot_path)
assay_name <- names(SummarizedExperiment::assays(object))[[1L]]
assay_matrix <- SummarizedExperiment::assay(object, assay_name)
row_metadata <- as.data.frame(SummarizedExperiment::rowData(object))
sample_metadata <- as.data.frame(SummarizedExperiment::colData(object))
tree <- TreeSummarizedExperiment::rowTree(object)
phyloseq_object <- suppressWarnings(
  mia::makePhyloseqFromTreeSummarizedExperiment(
    object,
    assay.type = assay_name
  )
)

column_sums <- colSums(assay_matrix)
phyloseq_column_sums <- phyloseq::sample_sums(phyloseq_object)
missing_percent <- pmax(0, 100 - column_sums)
sample_metadata_audit <- tibble::tibble(
  sample_id = colnames(object),
  subject_id = as.character(sample_metadata$subject_id),
  body_site = as.character(sample_metadata$body_site),
  study_condition = as.character(sample_metadata$study_condition),
  disease = as.character(sample_metadata$disease),
  number_reads = as.numeric(sample_metadata$number_reads),
  assay_sum_percent = as.numeric(column_sums),
  unrepresented_percent = as.numeric(missing_percent),
  finite_nonnegative = vapply(
    seq_len(ncol(assay_matrix)),
    function(i) {
      all(is.finite(assay_matrix[, i])) &&
        all(assay_matrix[, i] >= 0)
    },
    logical(1)
  )
) |>
  dplyr::mutate(
    status = ifelse(
      finite_nonnegative &
        number_reads > 0 &
        assay_sum_percent >= 90 &
        assay_sum_percent <= 100.1,
      "PASS",
      "FAIL"
    )
  )

check_rows <- list()
add_check <- function(category, id, observed, expected, passed) {
  check_rows[[length(check_rows) + 1L]] <<- data.frame(
    category = category,
    check_id = id,
    observed = paste(observed, collapse = ","),
    expected = paste(expected, collapse = ","),
    status = ifelse(isTRUE(passed), "PASS", "FAIL"),
    check.names = FALSE
  )
}

r_version <- paste(R.version$major, R.version$minor, sep = ".")
bioconductor_version <- as.character(BiocManager::version())
lock_sha256 <- sha256_file(lock_path)
snapshot_sha256 <- sha256_file(snapshot_path)
retrieval_log <- paste(readLines(retrieval_log_path, warn = FALSE), collapse = "\n")
notice <- paste(readLines(notice_path, warn = FALSE), collapse = "\n")
manifest <- resource_manifest[1L, , drop = FALSE]

add_check("runtime", "r-version", r_version, "4.4.1", r_version == "4.4.1")
add_check(
  "runtime",
  "r-architecture",
  paste0(.Machine$sizeof.pointer * 8L, "-bit"),
  "64-bit",
  .Machine$sizeof.pointer == 8L
)
add_check(
  "runtime",
  "bioconductor-version",
  bioconductor_version,
  "3.19",
  bioconductor_version == "3.19"
)
add_check(
  "runtime",
  "lock-r-version",
  lock$R$Version,
  "4.4.1",
  identical(lock$R$Version, "4.4.1")
)
add_check(
  "runtime",
  "lock-bioconductor-version",
  lock$Bioconductor$Version,
  "3.19",
  identical(lock$Bioconductor$Version, "3.19")
)
add_check(
  "runtime",
  "lock-record-count",
  length(lock$Packages),
  185L,
  length(lock$Packages) == 185L
)
add_check(
  "runtime",
  "lock-sha256",
  lock_sha256,
  "cf448ad154eb7412d7c069cbb9ea5c5bcef8c48a047a775409d9982f513d540e",
  identical(
    lock_sha256,
    "cf448ad154eb7412d7c069cbb9ea5c5bcef8c48a047a775409d9982f513d540e"
  )
)
add_check(
  "runtime",
  "package-contract-rows",
  nrow(package_contract),
  17L,
  nrow(package_contract) == 17L
)

for (i in seq_len(nrow(package_audit))) {
  package <- package_audit$package[[i]]
  add_check(
    "package",
    paste0(package, "-namespace"),
    package_audit$namespace_loaded[[i]],
    TRUE,
    package_audit$namespace_loaded[[i]]
  )
  add_check(
    "package",
    paste0(package, "-version"),
    package_audit$observed_version[[i]],
    package_audit$expected_version[[i]],
    package_audit$version_status[[i]] == "PASS"
  )
  add_check(
    "package",
    paste0(package, "-repository"),
    package_audit$observed_repository[[i]],
    package_audit$expected_repository[[i]],
    package_audit$repository_status[[i]] == "PASS"
  )
  add_check(
    "package",
    paste0(package, "-lock"),
    package_audit$lock_version[[i]],
    package_audit$expected_version[[i]],
    package_audit$lock_status[[i]] == "PASS"
  )
}

add_check(
  "resource",
  "resource-pattern",
  manifest$resource_pattern,
  "AsnicarF_2017.relative_abundance",
  manifest$resource_pattern == "AsnicarF_2017.relative_abundance"
)
add_check(
  "resource",
  "selected-title",
  manifest$selected_title,
  "2021-10-14.AsnicarF_2017.relative_abundance",
  manifest$selected_title ==
    "2021-10-14.AsnicarF_2017.relative_abundance"
)
add_check(
  "resource",
  "experimenthub-id",
  manifest$experimenthub_id,
  "EH7091",
  manifest$experimenthub_id == "EH7091"
)
add_check(
  "resource",
  "raw-cache-bytes",
  manifest$resource_cache_bytes,
  18187,
  manifest$resource_cache_bytes == 18187
)
add_check(
  "resource",
  "raw-cache-sha256",
  manifest$resource_cache_sha256,
  "ad631532fbbab39dfc3550a676a58310887d317e947392dcae9e0d4e4d69db27",
  manifest$resource_cache_sha256 ==
    "ad631532fbbab39dfc3550a676a58310887d317e947392dcae9e0d4e4d69db27"
)
add_check(
  "resource",
  "snapshot-bytes",
  file.info(snapshot_path)$size,
  manifest$snapshot_bytes,
  file.info(snapshot_path)$size == manifest$snapshot_bytes
)
add_check(
  "resource",
  "snapshot-sha256",
  snapshot_sha256,
  manifest$snapshot_sha256,
  snapshot_sha256 == manifest$snapshot_sha256
)
add_check(
  "resource",
  "manifest-snapshot-sha256",
  manifest$snapshot_sha256,
  "2952774730bff2af9e13c9c40058320aed524dc2dc7408b44dc6e697c06564b2",
  manifest$snapshot_sha256 ==
    "2952774730bff2af9e13c9c40058320aed524dc2dc7408b44dc6e697c06564b2"
)
add_check(
  "resource",
  "offline-replay",
  manifest$offline_replay_identical,
  TRUE,
  isTRUE(manifest$offline_replay_identical)
)
add_check(
  "resource",
  "manifest-package-version",
  manifest$package_version,
  "3.12.0",
  manifest$package_version == "3.12.0"
)
add_check(
  "resource",
  "manifest-bioconductor-version",
  manifest$bioconductor_version,
  "3.19",
  as.character(manifest$bioconductor_version) == "3.19"
)
add_check(
  "resource",
  "retrieval-log-id",
  grepl("Selected ExperimentHub ID: EH7091", retrieval_log, fixed = TRUE),
  TRUE,
  grepl("Selected ExperimentHub ID: EH7091", retrieval_log, fixed = TRUE)
)
add_check(
  "resource",
  "retrieval-log-offline",
  grepl("Offline replay identical: TRUE", retrieval_log, fixed = TRUE),
  TRUE,
  grepl("Offline replay identical: TRUE", retrieval_log, fixed = TRUE)
)
add_check(
  "resource",
  "package-source-sha256-notice",
  grepl(
    "611df3b405002fa5c221c90f12c7e3b0b4e8463c5b1a7b6b6db11f5376d1e1c2",
    notice,
    fixed = TRUE
  ),
  TRUE,
  grepl(
    "611df3b405002fa5c221c90f12c7e3b0b4e8463c5b1a7b6b6db11f5376d1e1c2",
    notice,
    fixed = TRUE
  )
)

add_check(
  "object",
  "object-class",
  class(object)[[1L]],
  "TreeSummarizedExperiment",
  class(object)[[1L]] == "TreeSummarizedExperiment"
)
add_check("object", "object-features", nrow(object), 298L, nrow(object) == 298L)
add_check("object", "object-samples", ncol(object), 24L, ncol(object) == 24L)
add_check(
  "object",
  "assay-count",
  length(SummarizedExperiment::assays(object)),
  1L,
  length(SummarizedExperiment::assays(object)) == 1L
)
add_check(
  "object",
  "assay-name",
  assay_name,
  "relative_abundance",
  assay_name == "relative_abundance"
)
add_check(
  "object",
  "row-metadata-columns",
  ncol(row_metadata),
  7L,
  ncol(row_metadata) == 7L
)
add_check(
  "object",
  "sample-metadata-columns",
  ncol(sample_metadata),
  22L,
  ncol(sample_metadata) == 22L
)
add_check(
  "object",
  "feature-ids-unique",
  anyDuplicated(rownames(object)),
  0L,
  anyDuplicated(rownames(object)) == 0L
)
add_check(
  "object",
  "sample-ids-unique",
  anyDuplicated(colnames(object)),
  0L,
  anyDuplicated(colnames(object)) == 0L
)
add_check(
  "object",
  "assay-feature-alignment",
  identical(rownames(assay_matrix), rownames(object)),
  TRUE,
  identical(rownames(assay_matrix), rownames(object))
)
add_check(
  "object",
  "assay-sample-alignment",
  identical(colnames(assay_matrix), colnames(object)),
  TRUE,
  identical(colnames(assay_matrix), colnames(object))
)
add_check(
  "object",
  "row-metadata-alignment",
  identical(rownames(row_metadata), rownames(object)),
  TRUE,
  identical(rownames(row_metadata), rownames(object))
)
add_check(
  "object",
  "sample-metadata-alignment",
  identical(rownames(sample_metadata), colnames(object)),
  TRUE,
  identical(rownames(sample_metadata), colnames(object))
)
add_check(
  "object",
  "tree-tip-alignment",
  sum(rownames(object) %in% tree$tip.label),
  nrow(object),
  all(rownames(object) %in% tree$tip.label)
)
add_check(
  "object",
  "assay-numeric",
  is.numeric(assay_matrix),
  TRUE,
  is.numeric(assay_matrix)
)
add_check(
  "object",
  "assay-finite",
  sum(!is.finite(assay_matrix)),
  0L,
  all(is.finite(assay_matrix))
)
add_check(
  "object",
  "assay-nonnegative",
  min(assay_matrix),
  ">=0",
  min(assay_matrix) >= 0
)
add_check(
  "object",
  "assay-percent-upper-bound",
  max(assay_matrix),
  "<=100",
  max(assay_matrix) <= 100.0001
)
add_check(
  "object",
  "assay-is-not-integer-counts",
  sum(abs(assay_matrix - round(assay_matrix)) > 1e-8),
  ">0",
  any(abs(assay_matrix - round(assay_matrix)) > 1e-8)
)
add_check(
  "object",
  "minimum-column-sum",
  signif(min(column_sums), 8L),
  94.65166,
  isTRUE(all.equal(min(column_sums), 94.65166, tolerance = 1e-8))
)
add_check(
  "object",
  "median-column-sum",
  signif(stats::median(column_sums), 8L),
  99.999995,
  isTRUE(
    all.equal(
      stats::median(column_sums),
      99.999995,
      tolerance = 1e-8
    )
  )
)
add_check(
  "object",
  "maximum-column-sum",
  signif(max(column_sums), 8L),
  100.00003,
  isTRUE(all.equal(max(column_sums), 100.00003, tolerance = 1e-8))
)
add_check(
  "object",
  "study-name",
  paste(unique(sample_metadata$study_name), collapse = ","),
  "AsnicarF_2017",
  identical(unique(sample_metadata$study_name), "AsnicarF_2017")
)
add_check(
  "object",
  "subject-count",
  length(unique(sample_metadata$subject_id)),
  15L,
  length(unique(sample_metadata$subject_id)) == 15L
)
add_check(
  "object",
  "pmid",
  paste(unique(sample_metadata$PMID), collapse = ","),
  "28144631",
  identical(as.character(unique(sample_metadata$PMID)), "28144631")
)
add_check(
  "object",
  "accession-count",
  length(unique(sample_metadata$NCBI_accession)),
  24L,
  length(unique(sample_metadata$NCBI_accession)) == 24L
)
add_check(
  "object",
  "sequencing-depth-positive",
  min(sample_metadata$number_reads),
  ">0",
  all(sample_metadata$number_reads > 0)
)
add_check(
  "object",
  "sample-audit-pass",
  sum(sample_metadata_audit$status == "PASS"),
  24L,
  all(sample_metadata_audit$status == "PASS")
)
add_check(
  "object",
  "phyloseq-feature-parity",
  phyloseq::ntaxa(phyloseq_object),
  nrow(object),
  phyloseq::ntaxa(phyloseq_object) == nrow(object)
)
add_check(
  "object",
  "phyloseq-sample-parity",
  phyloseq::nsamples(phyloseq_object),
  ncol(object),
  phyloseq::nsamples(phyloseq_object) == ncol(object)
)
add_check(
  "object",
  "phyloseq-library-parity",
  max(abs(phyloseq_column_sums[colnames(object)] - column_sums)),
  0,
  isTRUE(
    all.equal(
      as.numeric(phyloseq_column_sums[colnames(object)]),
      as.numeric(column_sums),
      tolerance = 0
    )
  )
)

validation_log <- dplyr::bind_rows(check_rows)
checks_total <- nrow(validation_log)
checks_passed <- sum(validation_log$status == "PASS")
checks_failed <- sum(validation_log$status == "FAIL")
package_checks_failed <- sum(
  validation_log$category == "package" &
    validation_log$status == "FAIL"
)

resource_audit <- validation_log |>
  dplyr::filter(category == "resource")
object_audit <- validation_log |>
  dplyr::filter(category == "object")

write_tsv(package_audit, "package-audit.tsv")
write_tsv(lock_audit, "lockfile-audit.tsv")
write_tsv(resource_audit, "resource-audit.tsv")
write_tsv(object_audit, "object-audit.tsv")
write_tsv(sample_metadata_audit, "sample-metadata-audit.tsv")
write_tsv(validation_log, "validation.log")

session_text <- capture.output(sessioninfo::session_info())
writeLines(
  sanitize_text(session_text),
  file.path(output_dir, "r-session-info.txt"),
  useBytes = TRUE
)

pal <- c(
  "Fixed runtime" = "#0072B2",
  "Online once" = "#D55E00",
  "Frozen offline" = "#009E73"
)

boundary_nodes <- tibble::tribble(
  ~x, ~title, ~evidence, ~scope,
  1, "R runtime", "R 4.4.1\nBioconductor 3.19", "Fixed runtime",
  2, "Package lock", "17 core packages\n185 lock records", "Fixed runtime",
  3, "Hub registry", "ExperimentHub\nEH7091", "Online once",
  4, "Resource cache", "18,187 bytes\nSHA-256 locked", "Online once",
  5, "TSE snapshot", "298 features\n24 samples", "Frozen offline",
  6, "Routine QA", "Local RDS only\nNo network", "Frozen offline"
)

boundary_plot <- ggplot2::ggplot(boundary_nodes, ggplot2::aes(x = x, y = 1)) +
  ggplot2::geom_segment(
    data = data.frame(x = 1:5, xend = 2:6, y = 1, yend = 1),
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.7,
    color = "#777777",
    arrow = grid::arrow(length = grid::unit(0.12, "inches"))
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = scope),
    size = 11,
    shape = 21,
    fill = "white",
    stroke = 2
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = title),
    vjust = -2.2,
    fontface = "bold",
    size = 3.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = evidence),
    vjust = 2.0,
    lineheight = 0.95,
    size = 3.1,
    color = "#333333"
  ) +
  ggplot2::scale_color_manual(values = pal) +
  ggplot2::scale_x_continuous(limits = c(0.7, 6.3)) +
  ggplot2::coord_cartesian(ylim = c(0.72, 1.28), clip = "off") +
  ggplot2::labs(
    title = "Separate installation, retrieval, and offline evidence",
    subtitle = "Only the Hub registry and first resource cache require network access",
    color = "Evidence scope",
    caption = "AsnicarF_2017.relative_abundance | ExperimentHub EH7091"
  ) +
  ggplot2::theme_void(base_size = 11, base_family = "sans") +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(size = 10.5),
    plot.caption = ggplot2::element_text(color = "#555555"),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(18, 20, 18, 20)
  )

save_publication_plot(
  boundary_plot,
  "12-r-data-access-boundaries",
  width = 10.2,
  height = 4.0
)

package_tiles <- tidyr::crossing(
  package = package_audit$package,
  check = c("Namespace", "Exact version", "Repository", "Lock record")
) |>
  dplyr::left_join(
    package_audit |>
      dplyr::select(
        package,
        namespace_loaded,
        version_status,
        repository_status,
        lock_status
      ),
    by = "package"
  ) |>
  dplyr::mutate(
    status = dplyr::case_when(
      check == "Namespace" & namespace_loaded ~ "PASS",
      check == "Exact version" & version_status == "PASS" ~ "PASS",
      check == "Repository" & repository_status == "PASS" ~ "PASS",
      check == "Lock record" & lock_status == "PASS" ~ "PASS",
      TRUE ~ "FAIL"
    ),
    package = factor(
      package,
      levels = rev(package_audit$package)
    ),
    check = factor(
      check,
      levels = c("Namespace", "Exact version", "Repository", "Lock record")
    )
  )

version_labels <- package_audit |>
  dplyr::mutate(
    package = factor(package, levels = rev(package_audit$package)),
    label = paste0("v", observed_version)
  )

package_plot <- ggplot2::ggplot(
  package_tiles,
  ggplot2::aes(x = check, y = package, fill = status)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = status),
    color = "white",
    fontface = "bold",
    size = 2.5
  ) +
  ggplot2::geom_text(
    data = version_labels,
    ggplot2::aes(x = 4.65, y = package, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.8,
    color = "#333333"
  ) +
  ggplot2::scale_fill_manual(values = c(PASS = "#009E73", FAIL = "#D55E00")) +
  ggplot2::coord_cartesian(xlim = c(0.5, 5.55), clip = "off") +
  ggplot2::labs(
    title = "The R ecosystem is a versioned contract",
    subtitle = "Every required namespace, repository, version, and lock record passes",
    x = NULL,
    y = NULL,
    caption = "R 4.4.1 | Bioconductor 3.19 | renv.lock: 185 records"
  ) +
  theme_pub(10.5) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8.5),
    plot.margin = ggplot2::margin(8, 65, 8, 8)
  )

save_publication_plot(
  package_plot,
  "12-package-role-contract",
  width = 8.6,
  height = 8.0
)

feature_means <- rowMeans(assay_matrix)
top_index <- head(order(feature_means, decreasing = TRUE), 15L)
species <- as.character(row_metadata$species[top_index])
species[is.na(species) | !nzchar(species)] <- rownames(object)[top_index][
  is.na(species) | !nzchar(species)
]
species <- make.unique(species)

heatmap_data <- as.data.frame(
  assay_matrix[top_index, , drop = FALSE],
  check.names = FALSE
)
heatmap_data$feature <- species
heatmap_long <- tidyr::pivot_longer(
  heatmap_data,
  cols = -feature,
  names_to = "sample_id",
  values_to = "relative_abundance"
) |>
  dplyr::left_join(
    tibble::tibble(
      sample_id = colnames(object),
      subject_id = sample_metadata$subject_id,
      body_site = sample_metadata$body_site
    ),
    by = "sample_id"
  ) |>
  dplyr::mutate(
    sample_id = factor(
      sample_id,
      levels = colnames(object)[
        order(sample_metadata$body_site, sample_metadata$subject_id)
      ]
    ),
    feature = factor(feature, levels = rev(species)),
    log_percent = log10(relative_abundance + 0.01)
  )

heatmap_plot <- ggplot2::ggplot(
  heatmap_long,
  ggplot2::aes(x = sample_id, y = feature, fill = log_percent)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.15) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(body_site),
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_fill_gradientn(
    colors = c("#F7FBFF", "#9ECAE1", "#3182BD", "#08306B"),
    name = "log10(percent + 0.01)"
  ) +
  ggplot2::labs(
    title = "A real TreeSummarizedExperiment survives the offline handoff",
    subtitle = "Top 15 taxa by mean abundance across the complete 24-sample object",
    x = "Sample ID",
    y = "Species",
    caption = paste0(
      "298 features × 24 samples | Sample sums: ",
      format(min(column_sums), digits = 4),
      "–",
      format(max(column_sums), digits = 4),
      "%"
    )
  ) +
  theme_pub(10.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 65,
      hjust = 1,
      vjust = 1,
      size = 6.8
    ),
    axis.text.y = ggplot2::element_text(size = 8.5),
    legend.position = "right"
  )

save_publication_plot(
  heatmap_plot,
  "12-cmd-object-contract",
  width = 10.0,
  height = 6.8
)

summary <- list(
  status = if (checks_failed == 0L) "passed" else "failed",
  validation_scope = "offline-sha256-locked-object-audit",
  qa_network_access = FALSE,
  r_version = r_version,
  r_architecture = paste0(.Machine$sizeof.pointer * 8L, "-bit"),
  bioconductor_version = bioconductor_version,
  curated_metagenomic_data_version = package_audit$observed_version[
    package_audit$package == "curatedMetagenomicData"
  ][[1L]],
  lock_records = length(lock$Packages),
  lock_sha256 = lock_sha256,
  package_contract_rows = nrow(package_contract),
  package_checks_failed = package_checks_failed,
  resource_pattern = manifest$resource_pattern[[1L]],
  selected_resource_title = manifest$selected_title[[1L]],
  experimenthub_id = manifest$experimenthub_id[[1L]],
  resource_cache_bytes = manifest$resource_cache_bytes[[1L]],
  resource_cache_sha256 = manifest$resource_cache_sha256[[1L]],
  snapshot_bytes = file.info(snapshot_path)$size,
  snapshot_sha256 = snapshot_sha256,
  retrieval_offline_replay_identical = isTRUE(
    manifest$offline_replay_identical[[1L]]
  ),
  object_class = class(object)[[1L]],
  object_features = nrow(object),
  object_samples = ncol(object),
  object_assay = assay_name,
  object_row_metadata_columns = ncol(row_metadata),
  object_sample_metadata_columns = ncol(sample_metadata),
  object_subjects = length(unique(sample_metadata$subject_id)),
  phyloseq_features = phyloseq::ntaxa(phyloseq_object),
  phyloseq_samples = phyloseq::nsamples(phyloseq_object),
  phyloseq_maximum_library_delta = max(
    abs(phyloseq_column_sums[colnames(object)] - column_sums)
  ),
  minimum_relative_abundance_sum = min(column_sums),
  median_relative_abundance_sum = stats::median(column_sums),
  maximum_relative_abundance_sum = max(column_sums),
  checks_total = checks_total,
  checks_passed = checks_passed,
  checks_failed = checks_failed,
  seed = 20260720
)

jsonlite::write_json(
  summary,
  file.path(output_dir, "environment-summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 12
)

if (checks_failed > 0L) {
  failed <- validation_log |>
    dplyr::filter(status == "FAIL")
  stop(
    "Article 12 validation failed: ",
    paste(failed$check_id, collapse = ", "),
    call. = FALSE
  )
}

cat(
  "Article 12 validation passed: ",
  checks_passed,
  "/",
  checks_total,
  " checks; object ",
  nrow(object),
  " × ",
  ncol(object),
  "; ExperimentHub ",
  manifest$experimenthub_id[[1L]],
  "; no network used.\n",
  sep = ""
)
