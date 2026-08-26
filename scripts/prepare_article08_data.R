#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(
    paste(
      "Usage: prepare_article08_data.R",
      "<read-prefix-metrics.tsv> <fastq-anatomy.tsv>",
      "<prefix-source-summary.json> <platform-benchmark.tsv>",
      "<ena-fastq-sources.tsv> <native-format-contract.tsv>",
      "<result-dir> <summary.json>"
    ),
    call. = FALSE
  )
}

metrics_path <- normalizePath(args[[1]], mustWork = TRUE)
anatomy_path <- normalizePath(args[[2]], mustWork = TRUE)
source_summary_path <- normalizePath(args[[3]], mustWork = TRUE)
benchmark_path <- normalizePath(args[[4]], mustWork = TRUE)
sources_path <- normalizePath(args[[5]], mustWork = TRUE)
native_contract_path <- normalizePath(args[[6]], mustWork = TRUE)
result_dir <- args[[7]]
summary_path <- args[[8]]

input_paths <- c(
  metrics = metrics_path,
  anatomy = anatomy_path,
  source_summary = source_summary_path,
  benchmark = benchmark_path,
  sources = sources_path,
  native_contract = native_contract_path
)
expected_hashes <- c(
  metrics = "6c426e0c89a95b1c1cc7ec631629e3635a7fd83c61c90d4086801b28495b5a61",
  anatomy = "591c7a7387232da0aa48c0cba7a3e1bd8d8a8b9b92ef6a0ed1f58e427199b62b",
  source_summary = "f7bc59e303f26f723aacb2ff8a1d065411c1dcface225ce6d16969baefa66d2d",
  benchmark = "4f5db378ec40b0e7ca16579aa3aadec770bba4bdd097c3bd8fa9ee826e7e7bf1",
  sources = "96638ae0d16ce5953760b596887c9f8f443c1dc33e0f82bef2e6ac040e69d962",
  native_contract = "26ddd05a1f56df7d13a5769090241b87ed264399b0a9044076ce0467fa351af8"
)
observed_hashes <- vapply(
  input_paths,
  function(path) {
    digest::digest(
      file = path,
      algo = "sha256",
      serialize = FALSE
    )
  },
  character(1)
)
stopifnot(identical(observed_hashes, expected_hashes))
set.seed(20260719)

read_tsv <- function(path) {
  utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

metrics <- read_tsv(metrics_path)
anatomy <- read_tsv(anatomy_path)
benchmark <- read_tsv(benchmark_path)
sources <- read_tsv(sources_path)
native_contract <- read_tsv(native_contract_path)
source_summary <- jsonlite::read_json(
  source_summary_path,
  simplifyVector = TRUE
)

platform_order <- c("Illumina", "ONT", "PacBio")
run_contract <- c(
  Illumina = "ERR9765746",
  ONT = "ERR9765780",
  PacBio = "ERR9765783"
)
expected_prefix_means <- c(
  Illumina = 149.8392,
  ONT = 4322.7572,
  PacBio = 9691.6746
)
expected_prefix_medians <- c(
  Illumina = 150,
  ONT = 4033.5,
  PacBio = 9204
)
expected_prefix_maxima <- c(
  Illumina = 150,
  ONT = 24698,
  PacBio = 25821
)

stopifnot(
  nrow(metrics) == 15000L,
  nrow(anatomy) == 4L,
  nrow(benchmark) == 3L,
  nrow(sources) == 4L,
  nrow(native_contract) == 3L,
  identical(sort(unique(metrics$PlatformKey)), sort(platform_order)),
  identical(sort(unique(benchmark$PlatformKey)), sort(platform_order)),
  identical(sort(unique(native_contract$PlatformKey)), sort(platform_order)),
  identical(unique(metrics$SampleAccession), "SAMEA14435832"),
  identical(unique(benchmark$SampleAccession), "SAMEA14435832"),
  identical(unique(sources$SampleAccession), "SAMEA14435832"),
  identical(source_summary$source_project, "PRJEB52977"),
  identical(source_summary$sample_accession, "SAMEA14435832"),
  isFALSE(source_summary$raw_fastq_stored),
  isTRUE(source_summary$illumina_prefix_mates_synchronized),
  source_summary$metrics_rows == 15000L,
  source_summary$anatomy_rows == 4L,
  identical(
    source_summary$metrics_sha256,
    expected_hashes[["metrics"]]
  ),
  identical(
    source_summary$anatomy_sha256,
    expected_hashes[["anatomy"]]
  ),
  all(table(metrics$PlatformKey) == 5000L),
  !anyDuplicated(
    paste(
      metrics$PlatformKey,
      metrics$ReadIndex,
      sep = "::"
    )
  ),
  all(metrics$ReadLength > 0),
  all(metrics$NCount >= 0),
  all(metrics$MinimumBaseQ >= 0),
  all(metrics$MaximumBaseQ <= 93),
  all(metrics$ExpectedAccuracyPercent >= 0),
  all(metrics$ExpectedAccuracyPercent <= 100),
  all(nchar(anatomy$SequencePrefix60) == 60L),
  all(nchar(anatomy$QualityPrefix60) == 60L),
  all(anatomy$Separator == "+"),
  all(vapply(
    strsplit(anatomy$DecodedPhredFirst12, ",", fixed = TRUE),
    length,
    integer(1)
  ) == 12L),
  all(sources$ENABytes > 0),
  all(nchar(sources$ENAReportedMD5) == 32L),
  all(sources$PrefixRecords == 5000L),
  identical(
    setNames(
      benchmark$RunAccession,
      benchmark$PlatformKey
    )[platform_order],
    run_contract
  ),
  identical(
    setNames(
      metrics$RunAccession[
        match(platform_order, metrics$PlatformKey)
      ],
      platform_order
    ),
    run_contract
  )
)

for (platform in platform_order) {
  platform_metrics <- metrics[
    metrics$PlatformKey == platform,
    ,
    drop = FALSE
  ]
  stopifnot(
    isTRUE(all.equal(
      mean(platform_metrics$ReadLength),
      expected_prefix_means[[platform]],
      tolerance = 1e-10
    )),
    isTRUE(all.equal(
      stats::median(platform_metrics$ReadLength),
      expected_prefix_medians[[platform]],
      tolerance = 1e-10
    )),
    max(platform_metrics$ReadLength) ==
      expected_prefix_maxima[[platform]]
  )
}

stopifnot(
  benchmark$ENAReadCount[
    benchmark$PlatformKey == "Illumina"
  ] == 41195050,
  benchmark$ENAReadCount[
    benchmark$PlatformKey == "ONT"
  ] == 696944,
  benchmark$ENAReadCount[
    benchmark$PlatformKey == "PacBio"
  ] == 524805,
  benchmark$MeanReadLengthBp[
    benchmark$PlatformKey == "Illumina"
  ] == 149,
  benchmark$MeanReadLengthBp[
    benchmark$PlatformKey == "ONT"
  ] == 4408.41,
  benchmark$MeanReadLengthBp[
    benchmark$PlatformKey == "PacBio"
  ] == 10289.7,
  benchmark$MeanMappedIdentityPercent[
    benchmark$PlatformKey == "Illumina"
  ] == 99.45,
  benchmark$MeanMappedIdentityPercent[
    benchmark$PlatformKey == "ONT"
  ] == 89.08,
  benchmark$MeanMappedIdentityPercent[
    benchmark$PlatformKey == "PacBio"
  ] == 99.72,
  benchmark$AssemblyN50Bp[
    benchmark$PlatformKey == "Illumina"
  ] == 13707,
  benchmark$AssemblyN50Bp[
    benchmark$PlatformKey == "ONT"
  ] == 759940,
  benchmark$AssemblyN50Bp[
    benchmark$PlatformKey == "PacBio"
  ] == 2013697,
  benchmark$RecoveredFullGenomes[
    benchmark$PlatformKey == "Illumina"
  ] == 12,
  benchmark$RecoveredFullGenomes[
    benchmark$PlatformKey == "ONT"
  ] == 22,
  benchmark$RecoveredFullGenomes[
    benchmark$PlatformKey == "PacBio"
  ] == 36
)

read_n50 <- function(lengths) {
  ordered <- sort(lengths, decreasing = TRUE)
  ordered[which(cumsum(ordered) >= sum(ordered) / 2)[1]]
}

prefix_summary <- do.call(
  rbind,
  lapply(
    platform_order,
    function(platform) {
      x <- metrics[
        metrics$PlatformKey == platform,
        ,
        drop = FALSE
      ]
      data.frame(
        PlatformKey = platform,
        DisplayLabel = unique(x$PlatformLabel),
        RunAccession = unique(x$RunAccession),
        PrefixRecords = nrow(x),
        MinimumReadLengthBp = min(x$ReadLength),
        MedianReadLengthBp = stats::median(x$ReadLength),
        MeanReadLengthBp = mean(x$ReadLength),
        P95ReadLengthBp = unname(
          stats::quantile(x$ReadLength, 0.95, type = 7)
        ),
        MaximumReadLengthBp = max(x$ReadLength),
        ReadN50Bp = read_n50(x$ReadLength),
        MedianErrorProbabilityMeanQ = stats::median(
          x$ErrorProbabilityMeanQ
        ),
        MeanGCPercent = mean(x$GCPercent),
        FractionAtLeast1kb = mean(x$ReadLength >= 1000),
        FractionAtLeast5kb = mean(x$ReadLength >= 5000),
        FractionAtLeast10kb = mean(x$ReadLength >= 10000),
        FractionAtLeast20kb = mean(x$ReadLength >= 20000),
        stringsAsFactors = FALSE
      )
    }
  )
)

span_thresholds <- c(
  50,
  100,
  150,
  200,
  300,
  500,
  750,
  1000,
  1500,
  2000,
  3000,
  5000,
  7500,
  10000,
  15000,
  20000,
  25000
)
span_survival <- do.call(
  rbind,
  lapply(
    platform_order,
    function(platform) {
      lengths <- metrics$ReadLength[
        metrics$PlatformKey == platform
      ]
      data.frame(
        PlatformKey = platform,
        DisplayLabel = unique(
          metrics$PlatformLabel[
            metrics$PlatformKey == platform
          ]
        ),
        ThresholdBp = span_thresholds,
        FractionSpanning = vapply(
          span_thresholds,
          function(threshold) {
            mean(lengths >= threshold)
          },
          numeric(1)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

anatomy_audit <- data.frame(
  PlatformKey = anatomy$PlatformKey,
  PlatformLabel = anatomy$PlatformLabel,
  RunAccession = anatomy$RunAccession,
  Mate = anatomy$Mate,
  HeaderStartsAt = startsWith(anatomy$Header, "@"),
  SeparatorStartsPlus = startsWith(anatomy$Separator, "+"),
  DisplaySequenceCharacters = nchar(anatomy$SequencePrefix60),
  DisplayQualityCharacters = nchar(anatomy$QualityPrefix60),
  FullReadLength = anatomy$FullReadLength,
  FullSequenceSHA256 = anatomy$FullSequenceSHA256,
  FullQualitySHA256 = anatomy$FullQualitySHA256,
  stringsAsFactors = FALSE
)

illumina_r1_header <- anatomy$Header[
  anatomy$PlatformKey == "Illumina" &
    anatomy$Mate == "R1"
]
illumina_r2_header <- anatomy$Header[
  anatomy$PlatformKey == "Illumina" &
    anatomy$Mate == "R2"
]
stopifnot(
  length(illumina_r1_header) == 1L,
  length(illumina_r2_header) == 1L,
  sub("/1$", "", illumina_r1_header) ==
    sub("/2$", "", illumina_r2_header)
)

illumina_length <- benchmark$MeanReadLengthBp[
  benchmark$PlatformKey == "Illumina"
]
illumina_n50 <- benchmark$AssemblyN50Bp[
  benchmark$PlatformKey == "Illumina"
]
benchmark_audit <- benchmark
benchmark_audit$MeanReadLengthFoldVsIllumina <- (
  benchmark$MeanReadLengthBp / illumina_length
)
benchmark_audit$AssemblyN50FoldVsIllumina <- (
  benchmark$AssemblyN50Bp / illumina_n50
)
benchmark_audit$AccessionSource <- (
  "ENA API checked 2026-07-19"
)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

write_tsv(
  prefix_summary,
  file.path(result_dir, "read-prefix-summary.tsv")
)
write_tsv(
  span_survival,
  file.path(result_dir, "span-survival.tsv")
)
write_tsv(
  anatomy_audit,
  file.path(result_dir, "fastq-anatomy-audit.tsv")
)
write_tsv(
  benchmark_audit,
  file.path(result_dir, "platform-benchmark-audit.tsv")
)
write_tsv(
  native_contract,
  file.path(result_dir, "native-format-audit.tsv")
)

value_for <- function(table, platform, column) {
  unname(table[table$PlatformKey == platform, column][[1]])
}

summary <- list(
  source_project = "PRJEB52977",
  source_sample = "SAMEA14435832",
  source_mock = "MOCK1; 71-strain synthetic metagenome",
  platforms = 3,
  source_fastq_files = 4,
  prefix_records_per_platform = 5000,
  prefix_metric_rows = nrow(metrics),
  anatomy_records = nrow(anatomy),
  raw_fastq_stored = FALSE,
  illumina_mates_synchronized = TRUE,
  illumina_run = run_contract[["Illumina"]],
  ont_run = run_contract[["ONT"]],
  pacbio_run = run_contract[["PacBio"]],
  illumina_full_mean_read_length_bp = value_for(
    benchmark,
    "Illumina",
    "MeanReadLengthBp"
  ),
  ont_full_mean_read_length_bp = value_for(
    benchmark,
    "ONT",
    "MeanReadLengthBp"
  ),
  pacbio_full_mean_read_length_bp = value_for(
    benchmark,
    "PacBio",
    "MeanReadLengthBp"
  ),
  illumina_full_mapped_identity_percent = value_for(
    benchmark,
    "Illumina",
    "MeanMappedIdentityPercent"
  ),
  ont_full_mapped_identity_percent = value_for(
    benchmark,
    "ONT",
    "MeanMappedIdentityPercent"
  ),
  pacbio_full_mapped_identity_percent = value_for(
    benchmark,
    "PacBio",
    "MeanMappedIdentityPercent"
  ),
  illumina_assembly_n50_bp = value_for(
    benchmark,
    "Illumina",
    "AssemblyN50Bp"
  ),
  ont_assembly_n50_bp = value_for(
    benchmark,
    "ONT",
    "AssemblyN50Bp"
  ),
  pacbio_assembly_n50_bp = value_for(
    benchmark,
    "PacBio",
    "AssemblyN50Bp"
  ),
  illumina_recovered_full_genomes = value_for(
    benchmark,
    "Illumina",
    "RecoveredFullGenomes"
  ),
  ont_recovered_full_genomes = value_for(
    benchmark,
    "ONT",
    "RecoveredFullGenomes"
  ),
  pacbio_recovered_full_genomes = value_for(
    benchmark,
    "PacBio",
    "RecoveredFullGenomes"
  ),
  ont_prefix_fraction_at_least_5kb = value_for(
    prefix_summary,
    "ONT",
    "FractionAtLeast5kb"
  ),
  pacbio_prefix_fraction_at_least_5kb = value_for(
    prefix_summary,
    "PacBio",
    "FractionAtLeast5kb"
  ),
  ont_prefix_fraction_at_least_10kb = value_for(
    prefix_summary,
    "ONT",
    "FractionAtLeast10kb"
  ),
  pacbio_prefix_fraction_at_least_10kb = value_for(
    prefix_summary,
    "PacBio",
    "FractionAtLeast10kb"
  ),
  metrics_sha256 = observed_hashes[["metrics"]],
  anatomy_sha256 = observed_hashes[["anatomy"]],
  selection_bias_note = paste(
    "Deterministic file prefixes support format and span auditing;",
    "they are not random full-run samples or current platform specifications."
  ),
  random_seed = 20260719
)
jsonlite::write_json(
  summary,
  summary_path,
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 15
)

message(
  sprintf(
    paste(
      "Article 08: %d prefix reads across %d platforms;",
      "assembly N50 %s/%s/%s bp."
    ),
    nrow(metrics),
    length(platform_order),
    format(summary$illumina_assembly_n50_bp, big.mark = ","),
    format(summary$ont_assembly_n50_bp, big.mark = ","),
    format(summary$pacbio_assembly_n50_bp, big.mark = ",")
  )
)
