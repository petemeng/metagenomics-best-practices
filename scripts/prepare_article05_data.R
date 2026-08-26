args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "Usage: prepare_article05_data.R <sample_data.rda> ",
    "<metaphlan2_profiles.rda> <mock_composition.rda> ",
    "<mock_profiles.tsv> <bias_summary.tsv> <syndna_benchmark.tsv> ",
    "<summary.json>"
  )
}

sample_data_path <- args[[1]]
profiles_path <- args[[2]]
mock_path <- args[[3]]
mock_profiles_path <- args[[4]]
bias_summary_path <- args[[5]]
syndna_benchmark_path <- args[[6]]
summary_path <- args[[7]]

load_named_object <- function(path, object_name) {
  object_env <- new.env(parent = baseenv())
  loaded <- load(path, envir = object_env)
  if (!(object_name %in% loaded)) {
    stop("Expected object ", object_name, " was not found in ", path)
  }
  get(object_name, envir = object_env, inherits = FALSE)
}

geometric_mean <- function(x) {
  if (any(!is.finite(x)) || any(x <= 0)) {
    stop("Geometric means require finite positive values")
  }
  exp(mean(log(x)))
}

center_composition <- function(x) {
  x / geometric_mean(x)
}

sample_data <- as.data.frame(
  load_named_object(
    sample_data_path,
    "costea2017_sample_data"
  ),
  stringsAsFactors = FALSE
)
profiles <- as.data.frame(
  load_named_object(
    profiles_path,
    "costea2017_metaphlan2_profiles"
  ),
  stringsAsFactors = FALSE
)
mock_raw <- as.data.frame(
  load_named_object(
    mock_path,
    "costea2017_mock_composition"
  ),
  stringsAsFactors = FALSE
)

required_sample_columns <- c(
  "Sample",
  "Protocol",
  "Individual",
  "Run_accession"
)
missing_sample_columns <- setdiff(
  required_sample_columns,
  names(sample_data)
)
if (length(missing_sample_columns) > 0L) {
  stop(
    "Missing sample-data columns: ",
    paste(missing_sample_columns, collapse = ", ")
  )
}
sample_data <- sample_data[
  !(sample_data$Sample %in% c("QA", "QB")),
  required_sample_columns,
  drop = FALSE
]
if (
  nrow(sample_data) != 27L ||
  !identical(
    as.integer(table(factor(
      sample_data$Protocol,
      levels = c("H", "Q", "W")
    ))),
    c(9L, 9L, 9L)
  )
) {
  stop("Expected 27 Costea samples, with 9 samples per protocol")
}

facs_column <- "bacterial cells in spike in Mix"
if (
  !all(c("Taxon", facs_column) %in% names(mock_raw)) ||
  anyNA(mock_raw[[facs_column]]) ||
  any(mock_raw[[facs_column]] <= 0)
) {
  stop("Mock composition is missing positive FACS cell measurements")
}
mock_truth <- stats::aggregate(
  mock_raw[[facs_column]],
  list(Taxon = mock_raw$Taxon),
  mean
)
names(mock_truth)[[2]] <- "FACS_cells"
mock_truth$Actual_mock_fraction <-
  mock_truth$FACS_cells /
  sum(mock_truth$FACS_cells)
mock_truth <- mock_truth[
  order(-mock_truth$FACS_cells, mock_truth$Taxon),
  ,
  drop = FALSE
]
rownames(mock_truth) <- NULL

if (!("Clade" %in% names(profiles))) {
  stop("MetaPhlAn2 profiles must contain a Clade column")
}
missing_runs <- setdiff(
  sample_data$Run_accession,
  names(profiles)
)
if (length(missing_runs) > 0L) {
  stop(
    "MetaPhlAn2 profiles are missing runs: ",
    paste(missing_runs, collapse = ", ")
  )
}

species_rows <-
  grepl("^k__Bacteria(?:\\||$)", profiles$Clade) &
  grepl("\\|s__[^|]+$", profiles$Clade)
species_taxa <- sub(
  ".*\\|s__",
  "",
  profiles$Clade[species_rows]
)
if (anyDuplicated(species_taxa)) {
  stop("Species-level MetaPhlAn2 identifiers must be unique")
}

profile_matrix <- as.matrix(
  profiles[
    species_rows,
    sample_data$Run_accession,
    drop = FALSE
  ]
)
storage.mode(profile_matrix) <- "double"
rownames(profile_matrix) <- species_taxa
colnames(profile_matrix) <- sample_data$Sample[
  match(
    colnames(profile_matrix),
    sample_data$Run_accession
  )
]
if (
  anyNA(profile_matrix) ||
  any(profile_matrix < 0) ||
  any(colSums(profile_matrix) <= 0)
) {
  stop("Species-level MetaPhlAn2 values must be finite and non-negative")
}
profile_matrix <- sweep(
  profile_matrix,
  2,
  colSums(profile_matrix),
  "/"
)

# These ambiguous profiles were excluded in the archived McLaren analysis
# before the mock subcomposition was evaluated.
excluded_taxa <- c(
  "Salmonella_unclassified",
  "Escherichia_coli",
  "Escherichia_unclassified"
)
profile_matrix <- profile_matrix[
  !(rownames(profile_matrix) %in% excluded_taxa),
  ,
  drop = FALSE
]
profile_matrix <- sweep(
  profile_matrix,
  2,
  colSums(profile_matrix),
  "/"
)

missing_mock_taxa <- setdiff(
  mock_truth$Taxon,
  rownames(profile_matrix)
)
if (length(missing_mock_taxa) > 0L) {
  stop(
    "MetaPhlAn2 profiles are missing mock taxa: ",
    paste(missing_mock_taxa, collapse = ", ")
  )
}
mock_matrix <- profile_matrix[
  mock_truth$Taxon,
  ,
  drop = FALSE
]
if (any(mock_matrix <= 0)) {
  stop("Every mock taxon must be detected in every Costea sample")
}
observed_mock <- sweep(
  mock_matrix,
  2,
  colSums(mock_matrix),
  "/"
)
actual_mock <- mock_truth$Actual_mock_fraction[
  match(
    rownames(observed_mock),
    mock_truth$Taxon
  )
]
error_matrix <- sweep(
  observed_mock,
  1,
  actual_mock,
  "/"
)
relative_efficiency_observations <- apply(
  error_matrix,
  2,
  center_composition
)

mock_profiles <- expand.grid(
  Taxon = rownames(observed_mock),
  Sample = colnames(observed_mock),
  stringsAsFactors = FALSE
)
mock_profiles$Protocol <- sample_data$Protocol[
  match(mock_profiles$Sample, sample_data$Sample)
]
mock_profiles$Specimen <- sample_data$Individual[
  match(mock_profiles$Sample, sample_data$Sample)
]
mock_profiles$FACS_cells <- mock_truth$FACS_cells[
  match(mock_profiles$Taxon, mock_truth$Taxon)
]
mock_profiles$Actual_mock_fraction <- actual_mock[
  match(mock_profiles$Taxon, rownames(observed_mock))
]
mock_profiles$Observed_mock_fraction <- as.vector(observed_mock)
mock_profiles$Relative_efficiency_observation <-
  as.vector(relative_efficiency_observations)
protocol_order <- c("H", "Q", "W")
specimen_order <- c(as.character(1:8), "M")
mock_profiles <- mock_profiles[
  order(
    match(mock_profiles$Protocol, protocol_order),
    match(mock_profiles$Specimen, specimen_order),
    match(mock_profiles$Taxon, mock_truth$Taxon)
  ),
  ,
  drop = FALSE
]
mock_profiles <- mock_profiles[
  ,
  c(
    "Sample",
    "Protocol",
    "Specimen",
    "Taxon",
    "FACS_cells",
    "Actual_mock_fraction",
    "Observed_mock_fraction",
    "Relative_efficiency_observation"
  ),
  drop = FALSE
]
rownames(mock_profiles) <- NULL

bias_parts <- lapply(
  protocol_order,
  function(protocol) {
    protocol_samples <- sample_data$Sample[
      sample_data$Protocol == protocol
    ]
    bias <- apply(
      relative_efficiency_observations[
        ,
        protocol_samples,
        drop = FALSE
      ],
      1,
      geometric_mean
    )
    bias <- center_composition(bias)
    pairwise_bias <- utils::combn(
      bias,
      2,
      function(values) {
        max(values) / min(values)
      }
    )
    data.frame(
      Protocol = protocol,
      Taxon = names(bias),
      Relative_efficiency = unname(bias),
      Log2_relative_efficiency = log2(unname(bias)),
      Maximum_pairwise_bias = max(pairwise_bias),
      Average_pairwise_bias = geometric_mean(pairwise_bias),
      stringsAsFactors = FALSE
    )
  }
)
bias_summary <- do.call(rbind, bias_parts)
bias_summary <- bias_summary[
  order(
    match(bias_summary$Protocol, protocol_order),
    match(bias_summary$Taxon, mock_truth$Taxon)
  ),
  ,
  drop = FALSE
]
rownames(bias_summary) <- NULL

syndna <- utils::read.delim(
  syndna_benchmark_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_syndna_columns <- c(
  "Taxon",
  "Expected_cell_percent",
  "SynDNA_pool2_percent",
  "Relative_read_percent"
)
missing_syndna_columns <- setdiff(
  required_syndna_columns,
  names(syndna)
)
if (length(missing_syndna_columns) > 0L) {
  stop(
    "Missing synDNA benchmark columns: ",
    paste(missing_syndna_columns, collapse = ", ")
  )
}
if (
  nrow(syndna) != 8L ||
  anyNA(syndna[, required_syndna_columns[-1], drop = FALSE]) ||
  any(
    as.matrix(
      syndna[
        ,
        required_syndna_columns[-1],
        drop = FALSE
      ]
    ) < 0
  )
) {
  stop("Expected eight non-negative synDNA mock benchmark rows")
}

syndna_correlation <- stats::cor(
  syndna$Expected_cell_percent,
  syndna$SynDNA_pool2_percent,
  method = "pearson"
)
relative_correlation <- stats::cor(
  syndna$Expected_cell_percent,
  syndna$Relative_read_percent,
  method = "pearson"
)
syndna_rmse <- sqrt(mean(
  (
    syndna$SynDNA_pool2_percent -
    syndna$Expected_cell_percent
  )^2
))
relative_rmse <- sqrt(mean(
  (
    syndna$Relative_read_percent -
    syndna$Expected_cell_percent
  )^2
))

protocol_metrics <- unique(
  bias_summary[
    ,
    c(
      "Protocol",
      "Maximum_pairwise_bias",
      "Average_pairwise_bias"
    ),
    drop = FALSE
  ]
)
protocol_metric <- function(protocol, column) {
  protocol_metrics[
    protocol_metrics$Protocol == protocol,
    column
  ][[1]]
}

summary_values <- list(
  costea_samples = ncol(observed_mock),
  costea_samples_per_protocol = 9,
  costea_protocols = length(protocol_order),
  costea_mock_taxa = nrow(observed_mock),
  costea_true_cell_span_fold =
    max(mock_truth$FACS_cells) /
    min(mock_truth$FACS_cells),
  protocol_h_maximum_pairwise_bias =
    protocol_metric("H", "Maximum_pairwise_bias"),
  protocol_q_maximum_pairwise_bias =
    protocol_metric("Q", "Maximum_pairwise_bias"),
  protocol_w_maximum_pairwise_bias =
    protocol_metric("W", "Maximum_pairwise_bias"),
  protocol_h_lowest_efficiency_taxon =
    bias_summary$Taxon[
      bias_summary$Protocol == "H"
    ][which.min(
      bias_summary$Relative_efficiency[
        bias_summary$Protocol == "H"
      ]
    )],
  protocol_h_highest_efficiency_taxon =
    bias_summary$Taxon[
      bias_summary$Protocol == "H"
    ][which.max(
      bias_summary$Relative_efficiency[
        bias_summary$Protocol == "H"
      ]
    )],
  syndna_mock_taxa = nrow(syndna),
  syndna_pool2_expected_pearson = syndna_correlation,
  relative_reads_expected_pearson = relative_correlation,
  syndna_pool2_rmse_percentage_points = syndna_rmse,
  relative_reads_rmse_percentage_points = relative_rmse
)

stopifnot(
  summary_values$costea_samples == 27L,
  summary_values$costea_samples_per_protocol == 9L,
  summary_values$costea_protocols == 3L,
  summary_values$costea_mock_taxa == 10L,
  summary_values$costea_true_cell_span_fold > 440,
  summary_values$costea_true_cell_span_fold < 450,
  abs(
    summary_values$protocol_h_maximum_pairwise_bias -
    2751.00178
  ) < 0.01,
  abs(
    summary_values$protocol_q_maximum_pairwise_bias -
    55.82398
  ) < 0.01,
  abs(
    summary_values$protocol_w_maximum_pairwise_bias -
    74.00882
  ) < 0.01,
  summary_values$protocol_h_lowest_efficiency_taxon ==
    "Lactobacillus_plantarum",
  summary_values$protocol_h_highest_efficiency_taxon ==
    "Fusobacterium_nucleatum",
  summary_values$syndna_mock_taxa == 8L,
  summary_values$syndna_pool2_expected_pearson > 0.99,
  summary_values$relative_reads_expected_pearson < -0.74,
  summary_values$syndna_pool2_rmse_percentage_points < 1,
  summary_values$relative_reads_rmse_percentage_points > 12
)

dir.create(
  dirname(mock_profiles_path),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  dirname(bias_summary_path),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  dirname(summary_path),
  recursive = TRUE,
  showWarnings = FALSE
)
utils::write.table(
  mock_profiles,
  mock_profiles_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
utils::write.table(
  bias_summary,
  bias_summary_path,
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
    sprintf(
      '  "%s": %s',
      name,
      json_scalar(summary_values[[name]])
    )
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
