args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5L) {
  stop(
    paste(
      "Usage: Rscript scripts/prepare_article06_data.R",
      "<marotz_source.tsv> <longhi_source.tsv>",
      "<read_budget.tsv> <saponin_tradeoff.tsv> <summary.json>"
    )
  )
}

marotz_path <- args[[1]]
longhi_path <- args[[2]]
read_budget_path <- args[[3]]
saponin_output_path <- args[[4]]
summary_path <- args[[5]]

required_files <- c(marotz_path, longhi_path)
if (!all(file.exists(required_files))) {
  stop("One or more Article 06 source tables are missing")
}

marotz <- read.delim(
  marotz_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)
longhi <- read.delim(
  longhi_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_marotz <- c(
  "method_code", "method_label", "treatment_stage", "mechanism",
  "host_mean_pct", "reported_plusminus_pct",
  "reagent_cost_usd_2018", "significant_vs_untreated",
  "participants", "replicates_per_participant", "source_location"
)
required_longhi <- c(
  "saponin_pct", "treatment_label", "sequenced_reads",
  "high_quality_reads", "nonhuman_reads", "host_filtered_pct",
  "gram_positive_pct", "gram_negative_pct", "source_location"
)

stopifnot(
  identical(names(marotz), required_marotz),
  identical(names(longhi), required_longhi),
  nrow(marotz) == 6L,
  nrow(longhi) == 7L,
  identical(marotz$method_code, c("raw", "fil", "neb", "mol", "qia", "lypma")),
  identical(longhi$saponin_pct, c(0, 0.0125, 0.05, 0.1, 0.5, 1.5, 2)),
  all(marotz$participants == 8L),
  all(marotz$replicates_per_participant == 3L),
  all(marotz$host_mean_pct >= 0 & marotz$host_mean_pct <= 100),
  all(longhi$host_filtered_pct >= 0 & longhi$host_filtered_pct <= 100),
  all(abs(
    longhi$gram_positive_pct + longhi$gram_negative_pct - 100
  ) < 1e-8),
  all(longhi$nonhuman_reads <= longhi$high_quality_reads),
  all(longhi$high_quality_reads <= longhi$sequenced_reads)
)

expected_host <- c(89.29, 89.69, 90.83, 62.88, 29.17, 8.53)
expected_longhi_host <- c(53.15, 36.38, 5.73, 2.58, 2.16, 2.73, 3.21)
expected_longhi_gram_negative <- c(
  82.75, 42.25, 40.25, 36.95, 34.98, 22.40, 24.25
)

stopifnot(
  identical(marotz$host_mean_pct, expected_host),
  identical(longhi$host_filtered_pct, expected_longhi_host),
  identical(longhi$gram_negative_pct, expected_longhi_gram_negative)
)

target_microbial_reads <- 1e7
marotz$microbial_pct <- 100 - marotz$host_mean_pct
marotz$microbial_yield_fold_vs_untreated <-
  marotz$microbial_pct / marotz$microbial_pct[marotz$method_code == "raw"]
marotz$total_reads_for_10m_microbial <-
  target_microbial_reads / (marotz$microbial_pct / 100)
marotz$total_reads_for_10m_microbial_millions <-
  marotz$total_reads_for_10m_microbial / 1e6

untreated_host <- longhi$host_filtered_pct[longhi$saponin_pct == 0]
untreated_gram_negative <-
  longhi$gram_negative_pct[longhi$saponin_pct == 0]

longhi$host_reduction_percentage_points <-
  untreated_host - longhi$host_filtered_pct
longhi$gram_negative_retention_vs_untreated_pct <-
  100 * longhi$gram_negative_pct / untreated_gram_negative
longhi$observed_nonhuman_fraction_pct <-
  100 * longhi$nonhuman_reads / longhi$high_quality_reads

stopifnot(
  max(abs(
    longhi$observed_nonhuman_fraction_pct -
      (100 - longhi$host_filtered_pct)
  )) < 0.01
)

dir.create(
  dirname(read_budget_path),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  dirname(saponin_output_path),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  dirname(summary_path),
  recursive = TRUE,
  showWarnings = FALSE
)

write.table(
  marotz,
  read_budget_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)
write.table(
  longhi,
  saponin_output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

raw_row <- marotz[marotz$method_code == "raw", , drop = FALSE]
lypma_row <- marotz[marotz$method_code == "lypma", , drop = FALSE]
lowest_host_row <- longhi[
  which.min(longhi$host_filtered_pct),
  ,
  drop = FALSE
]

summary <- list(
  article = 6L,
  seed = 20260719L,
  marotz_methods = nrow(marotz),
  marotz_participants = unique(marotz$participants),
  marotz_replicates_per_method = unique(
    marotz$replicates_per_participant
  ),
  untreated_host_percent = raw_row$host_mean_pct,
  lypma_host_percent = lypma_row$host_mean_pct,
  lypma_host_reduction_percentage_points =
    raw_row$host_mean_pct - lypma_row$host_mean_pct,
  lypma_microbial_yield_fold =
    lypma_row$microbial_yield_fold_vs_untreated,
  untreated_reads_for_10m_microbial =
    raw_row$total_reads_for_10m_microbial,
  lypma_reads_for_10m_microbial =
    lypma_row$total_reads_for_10m_microbial,
  longhi_saponin_conditions = nrow(longhi),
  longhi_specimens = 1L,
  longhi_untreated_host_percent = untreated_host,
  longhi_untreated_gram_negative_percent =
    untreated_gram_negative,
  longhi_lowest_host_saponin_percent =
    lowest_host_row$saponin_pct,
  longhi_lowest_host_percent =
    lowest_host_row$host_filtered_pct,
  longhi_gram_negative_at_lowest_host_percent =
    lowest_host_row$gram_negative_pct,
  longhi_gram_negative_retention_at_lowest_host_percent =
    lowest_host_row$gram_negative_retention_vs_untreated_pct
)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required")
}

jsonlite::write_json(
  summary,
  summary_path,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

message(
  sprintf(
    paste(
      "Prepared Article 06: %d Marotz methods;",
      "%d Longhi conditions; lyPMA microbial-yield fold %.3f"
    ),
    nrow(marotz),
    nrow(longhi),
    lypma_row$microbial_yield_fold_vs_untreated
  )
)
