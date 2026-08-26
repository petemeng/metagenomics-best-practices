#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(png)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop("Missing argument: ", flag)
  args[[index + 1L]]
}
cache_dir <- normalizePath(value_after("--cache-dir"), mustWork = TRUE)
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

SEED <- 70001L
PLOT_SEED <- 20260770L
set.seed(SEED)
Sys.setenv(TZ = "UTC")

write_tsv <- function(x, name) {
  connection <- if (grepl("\\.gz$", name)) {
    gzfile(file.path(output_dir, name), "wt")
  } else {
    file.path(output_dir, name)
  }
  on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

manifest <- fromJSON(file.path(cache_dir, "download-manifest.json"), simplifyVector = FALSE)
stopifnot(identical(manifest$article, 70L))
source_rdata <- file.path(cache_dir, "twosamplemr-vig-perform-mr.RData")
loaded <- load(source_rdata)
stopifnot(
  all(c("bmi_exp_dat", "chd_out_dat") %in% loaded),
  nrow(bmi_exp_dat) == 79L,
  nrow(chd_out_dat) == 79L,
  length(unique(bmi_exp_dat$SNP)) == 79L,
  length(unique(chd_out_dat$SNP)) == 79L,
  setequal(bmi_exp_dat$SNP, chd_out_dat$SNP),
  unique(bmi_exp_dat$id.exposure) == "ieu-a-2",
  unique(chd_out_dat$id.outcome) == "ieu-a-7",
  all(bmi_exp_dat$pval.exposure < 5e-8),
  all(chd_out_dat$samplesize.outcome == 184305L)
)

write_tsv(bmi_exp_dat, "bmi-instruments-raw.tsv.gz")
write_tsv(chd_out_dat, "chd-associations-raw.tsv.gz")
file.copy(
  source_rdata,
  file.path(output_dir, "twosamplemr-vig-perform-mr.RData"),
  overwrite = TRUE,
  copy.date = TRUE
)
file.copy(
  file.path(cache_dir, "download-manifest.json"),
  file.path(output_dir, "source-manifest.json"),
  overwrite = TRUE,
  copy.date = TRUE
)

attrition <- data.frame(
  StageOrder = 1:4,
  Stage = c(
    "Official vignette BMI instruments",
    "CHD associations available",
    "Unique SNP overlap",
    "Genome-wide significant exposure instruments"
  ),
  SNPs = c(
    nrow(bmi_exp_dat),
    nrow(chd_out_dat),
    length(intersect(bmi_exp_dat$SNP, chd_out_dat$SNP)),
    sum(bmi_exp_dat$pval.exposure < 5e-8)
  )
)
write_tsv(attrition, "input-attrition.tsv")

metrics <- list(
  article = 70,
  analysis_seed = SEED,
  plot_seed = PLOT_SEED,
  exposure_id = "ieu-a-2",
  outcome_id = "ieu-a-7",
  exposure_instruments = nrow(bmi_exp_dat),
  outcome_associations = nrow(chd_out_dat),
  snp_overlap = length(intersect(bmi_exp_dat$SNP, chd_out_dat$SNP)),
  exposure_sample_size_min = min(bmi_exp_dat$samplesize.exposure),
  exposure_sample_size_max = max(bmi_exp_dat$samplesize.exposure),
  exposure_study_maximum_n = 339224,
  outcome_sample_size = unique(chd_out_dat$samplesize.outcome),
  outcome_cases = 60801,
  outcome_controls = 123504,
  exposure_pvalue_max = max(bmi_exp_dat$pval.exposure),
  two_sample_mr_version = "0.7.9",
  two_sample_mr_commit = "3d119f20d6fc164b0c7f710f5590fee9580f2c7b"
)
writeLines(
  toJSON(metrics, pretty = TRUE, auto_unbox = TRUE, digits = 16),
  file.path(output_dir, "analysis-metrics.json")
)

contract <- list(
  article = 70,
  design = "two-sample Mendelian randomization with summary associations",
  exposure = "body mass index, original GWAS standardized scale (OpenGWAS ieu-a-2)",
  outcome = "coronary heart disease log odds (OpenGWAS ieu-a-7; 60,801 cases and 123,504 controls)",
  instruments = paste(
    "79 genome-wide significant variants from the official TwoSampleMR vignette cache;",
    "the cache does not preserve the LD reference-panel release or clumping ledger"
  ),
  harmonisation = "TwoSampleMR harmonise_data action=2; effect-allele frequency used for resolvable palindromes",
  primary_estimator = "multiplicative random-effects inverse-variance weighted MR",
  robustness_estimators = c("weighted median", "MR-Egger", "simple mode", "weighted mode"),
  diagnostics = c(
    "per-SNP F statistic and I2GX",
    "Cochran heterogeneity",
    "MR-Egger intercept",
    "MR-PRESSO global/outlier/distortion tests",
    "radial IVW",
    "leave-one-out",
    "Steiger directionality across CHD prevalence assumptions"
  ),
  mr_presso_null_distributions = 2000,
  steiger_prevalence_grid = c(0.03, 0.06, 0.10, 0.20),
  analysis_seed = SEED,
  plot_seed = PLOT_SEED,
  interpretation_limit = paste(
    "This BMI-to-CHD example teaches the statistical workflow; a microbiome MR",
    "claim additionally requires an external host-genetic microbiome GWAS with",
    "trait definition, ancestry, LD reference, sample overlap, and replication metadata."
  )
)
writeLines(
  toJSON(contract, pretty = TRUE, auto_unbox = TRUE, digits = 16),
  file.path(output_dir, "methods-contract.json")
)

# Render and crop the published Hemani et al. Figure 1 from the open eLife PDF.
temporary <- tempfile("article70-anchor-")
dir.create(temporary)
prefix <- file.path(temporary, "page")
status <- system2(
  "pdftoppm",
  args = c(
    "-f", "4", "-l", "4", "-singlefile", "-png", "-r", "240",
    file.path(cache_dir, "hemani-mrbase-paper.pdf"), prefix
  ),
  stdout = TRUE,
  stderr = TRUE
)
exit_status <- attr(status, "status")
if (!is.null(exit_status) && exit_status != 0L) stop(paste(status, collapse = "\n"))
page <- readPNG(paste0(prefix, ".png"))
height <- dim(page)[1]
width <- dim(page)[2]
rows <- seq.int(round(height * 0.067), round(height * 0.835))
columns <- seq.int(round(width * 0.060), round(width * 0.825))
anchor <- page[rows, columns, , drop = FALSE]
writePNG(anchor, file.path(output_dir, "hemani-figure1-original.png"))
unlink(temporary, recursive = TRUE)

versions <- data.frame(
  Package = c("R", "jsonlite", "png"),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("jsonlite")),
    as.character(packageVersion("png"))
  )
)
write_tsv(versions, "software-versions-preparation.tsv")

cat(toJSON(metrics, pretty = TRUE, auto_unbox = TRUE), "\n")
