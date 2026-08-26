args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: prepare_article04_data.R <meta_all.tsv> <LakeLanier.npo> ",
    "<library_sizes.tsv> <coverage_curve.tsv> <summary.json>"
  )
}

metadata_path <- args[[1]]
npo_path <- args[[2]]
library_path <- args[[3]]
curve_path <- args[[4]]
summary_path <- args[[5]]

metadata <- utils::read.delim(
  metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)

required_columns <- c("Sample_ID", "Study", "Group", "Library_Size")
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(metadata$Sample_ID)) {
  stop("Sample_ID values must be unique")
}

discovery_studies <- c("FR-CRC", "AT-CRC", "CN-CRC", "US-CRC", "DE-CRC")
validation_studies <- c("IT-CRC", "IT-CRC-2", "JP-CRC")
study_order <- c(discovery_studies, validation_studies)
group_order <- c("CTR", "CRC")

crc <- metadata[
  metadata$Study %in% study_order & metadata$Group %in% group_order,
  required_columns,
  drop = FALSE
]
if (anyNA(crc$Library_Size) || any(crc$Library_Size <= 0)) {
  stop("Library_Size must be present and positive for all CRC analysis samples")
}
crc$Role <- ifelse(
  crc$Study %in% discovery_studies,
  "Discovery / meta-analysis",
  "Independent validation"
)
crc$Group <- ifelse(crc$Group == "CTR", "Control", "CRC")
crc$Study <- factor(crc$Study, levels = study_order)
crc$Group <- factor(crc$Group, levels = c("Control", "CRC"))
crc <- crc[order(crc$Study, crc$Group, crc$Sample_ID), , drop = FALSE]
crc$Study <- as.character(crc$Study)
crc$Group <- as.character(crc$Group)
rownames(crc) <- NULL

read_npo_metadata <- function(path) {
  headers <- grep("^# @", readLines(path, warn = FALSE), value = TRUE)
  keys <- sub(":.*$", "", sub("^# @", "", headers))
  values <- sub("^.*: ", "", headers)
  stats::setNames(values, keys)
}

npo_metadata <- read_npo_metadata(npo_path)
required_headers <- c("version", "L", "R", "overlap")
missing_headers <- setdiff(required_headers, names(npo_metadata))
if (length(missing_headers) > 0L) {
  stop("Missing NPO headers: ", paste(missing_headers, collapse = ", "))
}

npo <- utils::read.table(
  npo_path,
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  col.names = c(
    "Read_subsample",
    "Mean_redundancy",
    "SD_redundancy",
    "Q1_redundancy",
    "Median_redundancy",
    "Q3_redundancy"
  )
)
npo <- npo[order(npo$Read_subsample), , drop = FALSE]

mean_read_length <- as.numeric(npo_metadata[["L"]])
total_reads <- as.numeric(npo_metadata[["R"]])
overlap_pct <- as.numeric(npo_metadata[["overlap"]])
if (
  !is.finite(mean_read_length) ||
  !is.finite(total_reads) ||
  !is.finite(overlap_pct)
) {
  stop("NPO L, R and overlap headers must be numeric")
}

# Nonpareil converts alignment-kernel redundancy to abundance-weighted
# coverage with the overlap correction reported by Rodriguez-R &
# Konstantinidis (2014) and implemented in Nonpareil 3.5.5.
coverage_exponent <- 1 - exp(2.23e-2 * overlap_pct - 3.5698)
redundancy_columns <- c(
  "Mean_redundancy",
  "Q1_redundancy",
  "Median_redundancy",
  "Q3_redundancy"
)
coverage_columns <- c(
  "Mean_coverage",
  "Q1_coverage",
  "Median_coverage",
  "Q3_coverage"
)
for (i in seq_along(redundancy_columns)) {
  npo[[coverage_columns[[i]]]] <-
    npo[[redundancy_columns[[i]]]]^coverage_exponent
}

final_coverage <- utils::tail(npo$Mean_coverage, 1L)
max_observed_reads <- max(npo$Read_subsample)
effort_shape <- final_coverage^0.27
npo$Adjusted_read_subsample <- ifelse(
  npo$Read_subsample == 0,
  0,
  exp(
    log(max_observed_reads) +
      effort_shape * (
        log(npo$Read_subsample) - log(max_observed_reads)
      )
  )
)
npo$Adjusted_effort_bp <-
  npo$Adjusted_read_subsample *
  mean_read_length *
  total_reads /
  max(npo$Adjusted_read_subsample)

curve_columns <- c(
  "Read_subsample",
  "Adjusted_effort_bp",
  "Mean_redundancy",
  "Q1_redundancy",
  "Median_redundancy",
  "Q3_redundancy",
  "Mean_coverage",
  "Q1_coverage",
  "Median_coverage",
  "Q3_coverage"
)
curve <- npo[, curve_columns, drop = FALSE]

summary_values <- list(
  crc_samples = nrow(crc),
  crc_studies = length(unique(crc$Study)),
  crc_median_library_reads = stats::median(crc$Library_Size),
  crc_minimum_library_reads = min(crc$Library_Size),
  crc_maximum_library_reads = max(crc$Library_Size),
  lake_curve_points = nrow(curve),
  lake_total_reads = total_reads,
  lake_mean_read_length = mean_read_length,
  lake_sequencing_effort_bp = total_reads * mean_read_length,
  lake_overlap_pct = overlap_pct,
  lake_final_redundancy = utils::tail(npo$Mean_redundancy, 1L),
  lake_final_coverage = final_coverage,
  coverage_exponent = coverage_exponent,
  npo_generator_version = unname(npo_metadata[["version"]])
)

stopifnot(
  summary_values$crc_samples == 768L,
  summary_values$crc_studies == 8L,
  summary_values$crc_median_library_reads == 24019758.5,
  summary_values$crc_minimum_library_reads == 1932775,
  summary_values$crc_maximum_library_reads == 90924580,
  summary_values$lake_curve_points == 47L,
  summary_values$lake_total_reads == 13551190,
  abs(summary_values$lake_mean_read_length - 86.158) < 1e-8,
  abs(summary_values$lake_sequencing_effort_bp - 1167543428.02) < 0.01,
  abs(summary_values$lake_final_redundancy - 0.65906) < 1e-8,
  round(100 * summary_values$lake_final_coverage) == 68L,
  summary_values$npo_generator_version == "2.40"
)

dir.create(dirname(library_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(curve_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  crc,
  library_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
utils::write.table(
  curve,
  curve_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

json_scalar <- function(value) {
  if (length(value) != 1L || is.na(value)) {
    stop("Summary JSON values must be scalar and non-missing")
  }
  if (is.character(value)) {
    escaped <- gsub("\\\\", "\\\\\\\\", value)
    escaped <- gsub('"', '\\\\"', escaped, fixed = TRUE)
    return(sprintf('"%s"', escaped))
  }
  if (is.logical(value)) {
    return(tolower(as.character(value)))
  }
  format(value, scientific = FALSE, digits = 15, trim = TRUE)
}

json_lines <- vapply(
  names(summary_values),
  function(name) {
    sprintf('  "%s": %s', name, json_scalar(summary_values[[name]]))
  },
  character(1)
)
json <- paste0(
  "{\n",
  paste(json_lines, collapse = ",\n"),
  "\n}\n"
)
writeLines(json, summary_path, useBytes = TRUE)
cat(json)
