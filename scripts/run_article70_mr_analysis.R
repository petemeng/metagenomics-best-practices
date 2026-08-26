#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MRPRESSO)
  library(RadialMR)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop("Missing argument: ", flag)
  args[[index + 1L]]
}
input_dir <- normalizePath(value_after("--input-dir"), mustWork = TRUE)
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

SEED <- 70001L
PLOT_SEED <- 20260770L
PRESSO_DISTRIBUTIONS <- 2000L
OUTCOME_CASES <- 60801L
OUTCOME_CONTROLS <- 123504L
PREVALENCE_GRID <- c(0.03, 0.06, 0.10, 0.20)
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

exposure <- read.delim(
  gzfile(file.path(input_dir, "bmi-instruments-raw.tsv.gz")),
  check.names = FALSE
)
outcome <- read.delim(
  gzfile(file.path(input_dir, "chd-associations-raw.tsv.gz")),
  check.names = FALSE
)
stopifnot(
  nrow(exposure) == 79L,
  nrow(outcome) == 79L,
  setequal(exposure$SNP, outcome$SNP),
  all(exposure$pval.exposure < 5e-8),
  packageVersion("TwoSampleMR") == "0.7.9",
  packageVersion("ieugwasr") == "1.1.0"
)

harmonised <- harmonise_data(exposure, outcome, action = 2)
harmonisation_log <- attr(harmonised, "log")
stopifnot(nrow(harmonised) == 79L, sum(harmonised$mr_keep) == 79L)

raw_outcome_index <- match(harmonised$SNP, outcome$SNP)
harmonised$beta.outcome.raw <- outcome$beta.outcome[raw_outcome_index]
harmonised$effect_allele.outcome.raw <- outcome$effect_allele.outcome[raw_outcome_index]
harmonised$other_allele.outcome.raw <- outcome$other_allele.outcome[raw_outcome_index]
harmonised$OutcomeBetaFlipped <- abs(
  harmonised$beta.outcome + harmonised$beta.outcome.raw
) < 1e-12 & abs(harmonised$beta.outcome - harmonised$beta.outcome.raw) > 1e-12
harmonised$FStatistic <- (harmonised$beta.exposure / harmonised$se.exposure)^2
harmonised$r.exposure <- get_r_from_bsen(
  harmonised$beta.exposure,
  harmonised$se.exposure,
  harmonised$samplesize.exposure
)
harmonised$R2Exposure <- harmonised$r.exposure^2
orientation <- ifelse(harmonised$beta.exposure < 0, -1, 1)
harmonised$BetaExposureOriented <- harmonised$beta.exposure * orientation
harmonised$BetaOutcomeOriented <- harmonised$beta.outcome * orientation
write_tsv(harmonised, "harmonised-instruments.tsv.gz")
saveRDS(harmonised, file.path(output_dir, "harmonised-instruments.rds"), compress = "xz")

harmonisation_audit <- data.frame(
  Quantity = c(
    "Exposure instruments",
    "Outcome associations",
    "Matched SNPs",
    "Retained for MR",
    "Palindromic SNPs",
    "Ambiguous palindromic SNPs",
    "Incompatible alleles",
    "Outcome beta sign flips",
    "Proxy variants"
  ),
  Count = c(
    nrow(exposure),
    nrow(outcome),
    nrow(harmonised),
    sum(harmonised$mr_keep),
    sum(harmonised$palindromic),
    sum(harmonised$ambiguous),
    sum(harmonised$remove),
    sum(harmonised$OutcomeBetaFlipped),
    harmonisation_log$proxy_variants
  )
)
write_tsv(harmonisation_audit, "harmonisation-audit.tsv")

i2gx <- Isq(harmonised$beta.exposure, harmonised$se.exposure)
strength <- data.frame(
  Quantity = c(
    "Instrument count",
    "Minimum F",
    "Median F",
    "Mean F",
    "Maximum F",
    "Instruments with F < 10",
    "I2GX",
    "Sum of approximate exposure R2",
    "Approximate exposure variance explained (%)"
  ),
  Value = c(
    nrow(harmonised),
    min(harmonised$FStatistic),
    median(harmonised$FStatistic),
    mean(harmonised$FStatistic),
    max(harmonised$FStatistic),
    sum(harmonised$FStatistic < 10),
    i2gx,
    sum(harmonised$R2Exposure),
    100 * sum(harmonised$R2Exposure)
  )
)
write_tsv(strength, "instrument-strength-summary.tsv")

method_list <- c(
  "mr_ivw",
  "mr_weighted_median",
  "mr_egger_regression",
  "mr_simple_mode",
  "mr_weighted_mode"
)
set.seed(SEED + 100L)
mr_results <- mr(harmonised, method_list = method_list)
mr_results <- generate_odds_ratios(mr_results)
mr_results$EffectScale <- "CHD odds ratio per genetically predicted 1-SD higher BMI"
write_tsv(mr_results, "mr-estimates.tsv")

heterogeneity <- mr_heterogeneity(
  harmonised,
  method_list = c("mr_egger_regression", "mr_ivw")
)
heterogeneity$I2Percent <- pmax(0, 100 * (heterogeneity$Q - heterogeneity$Q_df) / heterogeneity$Q)
write_tsv(heterogeneity, "mr-heterogeneity.tsv")
egger_intercept <- mr_pleiotropy_test(harmonised)
write_tsv(egger_intercept, "egger-intercept.tsv")

set.seed(SEED + 200L)
single_snp <- mr_singlesnp(
  harmonised,
  all_method = c("mr_ivw", "mr_egger_regression")
)
single_snp$CILower <- single_snp$b - 1.96 * single_snp$se
single_snp$CIUpper <- single_snp$b + 1.96 * single_snp$se
single_snp$OR <- exp(single_snp$b)
single_snp$ORLower <- exp(single_snp$CILower)
single_snp$ORUpper <- exp(single_snp$CIUpper)
write_tsv(single_snp, "single-snp-estimates.tsv.gz")

leave_one_out <- mr_leaveoneout(harmonised)
leave_one_out$CILower <- leave_one_out$b - 1.96 * leave_one_out$se
leave_one_out$CIUpper <- leave_one_out$b + 1.96 * leave_one_out$se
leave_one_out$OR <- exp(leave_one_out$b)
leave_one_out$ORLower <- exp(leave_one_out$CILower)
leave_one_out$ORUpper <- exp(leave_one_out$CIUpper)
write_tsv(leave_one_out, "leave-one-out.tsv")

set.seed(SEED + 1000L)
presso_list <- run_mr_presso(
  harmonised,
  NbDistribution = PRESSO_DISTRIBUTIONS,
  SignifThreshold = 0.05
)
presso <- presso_list[[1L]]
presso_main <- presso[["Main MR results"]]
presso_estimates <- data.frame(
  Analysis = presso_main[["MR Analysis"]],
  Estimate = presso_main[["Causal Estimate"]],
  StandardError = presso_main[["Sd"]],
  TStatistic = presso_main[["T-stat"]],
  PValue = presso_main[["P-value"]]
)
presso_estimates$CILower <- presso_estimates$Estimate - 1.96 * presso_estimates$StandardError
presso_estimates$CIUpper <- presso_estimates$Estimate + 1.96 * presso_estimates$StandardError
presso_estimates$OddsRatio <- exp(presso_estimates$Estimate)
presso_estimates$OddsRatioLower <- exp(presso_estimates$CILower)
presso_estimates$OddsRatioUpper <- exp(presso_estimates$CIUpper)
write_tsv(presso_estimates, "mr-presso-estimates.tsv")

presso_results <- presso[["MR-PRESSO results"]]
presso_indices <- as.integer(presso_results[["Distortion Test"]][["Outliers Indices"]])
stopifnot(length(presso_indices) >= 1L, all(presso_indices >= 1L & presso_indices <= nrow(harmonised)))
presso_outlier_table <- presso_results[["Outlier Test"]]
presso_outliers <- data.frame(
  RowIndex = seq_len(nrow(harmonised)),
  SNP = harmonised$SNP,
  RSSObserved = presso_outlier_table[["RSSobs"]],
  PValueText = as.character(presso_outlier_table[["Pvalue"]]),
  DistortionOutlier = seq_len(nrow(harmonised)) %in% presso_indices
)
write_tsv(presso_outliers, "mr-presso-outliers.tsv")
presso_tests <- data.frame(
  Test = c("Global", "Distortion"),
  Statistic = c(
    presso_results[["Global Test"]][["RSSobs"]],
    unname(presso_results[["Distortion Test"]][["Distortion Coefficient"]])
  ),
  PValueText = c(
    as.character(presso_results[["Global Test"]][["Pvalue"]]),
    format(presso_results[["Distortion Test"]][["Pvalue"]], digits = 16)
  ),
  OutlierCount = c(length(presso_indices), length(presso_indices)),
  OutlierSNPs = paste(harmonised$SNP[presso_indices], collapse = ";")
)
write_tsv(presso_tests, "mr-presso-tests.tsv")
saveRDS(presso_list, file.path(output_dir, "mr-presso-object.rds"), compress = "xz")

radial_data <- format_radial(
  harmonised$beta.exposure,
  harmonised$beta.outcome,
  harmonised$se.exposure,
  harmonised$se.outcome,
  RSID = harmonised$SNP
)
radial_fit <- ivw_radial(
  radial_data,
  alpha = 0.05,
  weights = 3,
  tol = 1e-4,
  summary = FALSE
)
write_tsv(radial_fit$outliers, "radial-ivw-outliers.tsv")
radial_estimates <- data.frame(
  Model = rownames(radial_fit$coef),
  radial_fit$coef,
  row.names = NULL,
  check.names = FALSE
)
write_tsv(radial_estimates, "radial-ivw-estimates.tsv")

set.seed(SEED + 1100L)
outlier_exclusion <- mr(
  harmonised[-presso_indices, , drop = FALSE],
  method_list = method_list
)
outlier_exclusion <- generate_odds_ratios(outlier_exclusion)
outlier_exclusion$RemovedSNPs <- paste(harmonised$SNP[presso_indices], collapse = ";")
write_tsv(outlier_exclusion, "presso-outlier-exclusion-estimates.tsv")

steiger_rows <- list()
for (index in seq_along(PREVALENCE_GRID)) {
  prevalence <- PREVALENCE_GRID[[index]]
  steiger_data <- harmonised
  steiger_data$r.outcome <- get_r_from_lor(
    steiger_data$beta.outcome,
    steiger_data$eaf.outcome,
    OUTCOME_CASES,
    OUTCOME_CONTROLS,
    prevalence
  )
  result <- directionality_test(steiger_data)
  steiger_rows[[index]] <- data.frame(
    AssumedCHDPrevalence = prevalence,
    ExposureR2 = result$snp_r2.exposure,
    OutcomeLiabilityR2 = result$snp_r2.outcome,
    CorrectCausalDirection = result$correct_causal_direction,
    SteigerPValue = result$steiger_pval
  )
}
steiger <- do.call(rbind, steiger_rows)
write_tsv(steiger, "steiger-directionality.tsv")

steiger_primary <- harmonised
steiger_primary$r.outcome <- get_r_from_lor(
  steiger_primary$beta.outcome,
  steiger_primary$eaf.outcome,
  OUTCOME_CASES,
  OUTCOME_CONTROLS,
  0.06
)
steiger_per_snp <- data.frame(
  SNP = steiger_primary$SNP,
  ExposureR2 = steiger_primary$r.exposure^2,
  OutcomeLiabilityR2 = steiger_primary$r.outcome^2,
  ExposureExplainsMore = steiger_primary$r.exposure^2 > steiger_primary$r.outcome^2
)
write_tsv(steiger_per_snp, "steiger-per-snp.tsv")

design_audit <- data.frame(
  Item = c(
    "Genome-wide exposure threshold",
    "Effect-allele harmonisation",
    "Palindromic variants",
    "LD clumping reference and release",
    "Exposure scale in cache",
    "Exposure and outcome ancestry",
    "Participant overlap",
    "Outcome case/control counts",
    "Outcome population prevalence",
    "External replication"
  ),
  Status = c(
    "Verified in frozen associations",
    "Verified with harmonisation ledger",
    "Four retained; none ambiguous",
    "Not encoded in vignette cache",
    "Recovered from source publication; not encoded in cache",
    "Not encoded in cache; audit source GWAS",
    "Not encoded in cache",
    "Recovered from source publication",
    "Not encoded; sensitivity grid used",
    "Not part of this teaching example"
  ),
  Consequence = c(
    "All 79 SNPs have P < 5e-8",
    "No outcome beta sign flips were required",
    "Allele frequencies support action=2 retention",
    "Do not claim clumping was independently reconstructed",
    "Interpret estimate per original standardized BMI scale",
    "Transportability and LD comparability remain assumptions",
    "Bias direction depends on overlap and instrument strength",
    "60,801 cases and 123,504 controls",
    "Steiger liability R2 depends on prevalence",
    "Do not treat this single analysis as final causal proof"
  )
)
write_tsv(design_audit, "design-audit.tsv")

ivw <- mr_results[mr_results$method == "Inverse variance weighted", , drop = FALSE]
weighted_median <- mr_results[mr_results$method == "Weighted median", , drop = FALSE]
egger <- mr_results[mr_results$method == "MR Egger", , drop = FALSE]
ivw_heterogeneity <- heterogeneity[
  heterogeneity$method == "Inverse variance weighted", , drop = FALSE
]
metrics <- list(
  article = 70,
  analysis_seed = SEED,
  plot_seed = PLOT_SEED,
  instruments = nrow(harmonised),
  retained_instruments = sum(harmonised$mr_keep),
  palindromic_instruments = sum(harmonised$palindromic),
  ambiguous_instruments = sum(harmonised$ambiguous),
  outcome_beta_flips = sum(harmonised$OutcomeBetaFlipped),
  minimum_f = min(harmonised$FStatistic),
  mean_f = mean(harmonised$FStatistic),
  i2gx = i2gx,
  exposure_r2_sum = sum(harmonised$R2Exposure),
  ivw_beta = ivw$b,
  ivw_se = ivw$se,
  ivw_p = ivw$pval,
  ivw_or = ivw$or,
  ivw_or_ci_lower = ivw$or_lci95,
  ivw_or_ci_upper = ivw$or_uci95,
  weighted_median_beta = weighted_median$b,
  weighted_median_p = weighted_median$pval,
  egger_beta = egger$b,
  egger_p = egger$pval,
  ivw_q = ivw_heterogeneity$Q,
  ivw_q_df = ivw_heterogeneity$Q_df,
  ivw_q_p = ivw_heterogeneity$Q_pval,
  egger_intercept = egger_intercept$egger_intercept,
  egger_intercept_p = egger_intercept$pval,
  mr_presso_distributions = PRESSO_DISTRIBUTIONS,
  mr_presso_global_p = as.character(presso_results[["Global Test"]][["Pvalue"]]),
  mr_presso_outlier_count = length(presso_indices),
  mr_presso_outliers = harmonised$SNP[presso_indices],
  mr_presso_corrected_beta = presso_estimates$Estimate[presso_estimates$Analysis == "Outlier-corrected"],
  mr_presso_distortion_p = presso_results[["Distortion Test"]][["Pvalue"]],
  radial_nominal_outliers = nrow(radial_fit$outliers),
  leave_one_out_beta_min = min(leave_one_out$b[leave_one_out$SNP != "All"]),
  leave_one_out_beta_max = max(leave_one_out$b[leave_one_out$SNP != "All"]),
  steiger_all_prevalences_forward = all(steiger$CorrectCausalDirection),
  steiger_per_snp_forward = sum(steiger_per_snp$ExposureExplainsMore)
)
writeLines(
  toJSON(metrics, pretty = TRUE, auto_unbox = TRUE, digits = 16),
  file.path(output_dir, "model-metrics.json")
)

versions <- data.frame(
  Package = c(
    "R", "TwoSampleMR", "ieugwasr", "MRPRESSO", "RadialMR",
    "MendelianRandomization", "jsonlite"
  ),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("TwoSampleMR")),
    as.character(packageVersion("ieugwasr")),
    as.character(packageVersion("MRPRESSO")),
    as.character(packageVersion("RadialMR")),
    as.character(packageVersion("MendelianRandomization")),
    as.character(packageVersion("jsonlite"))
  )
)
write_tsv(versions, "software-versions-analysis.tsv")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "r-session-info.txt"))
saveRDS(
  list(
    harmonised = harmonised,
    mr = mr_results,
    heterogeneity = heterogeneity,
    egger_intercept = egger_intercept,
    leave_one_out = leave_one_out,
    steiger = steiger
  ),
  file.path(output_dir, "mr-analysis-object.rds"),
  compress = "xz"
)

cat(toJSON(metrics, pretty = TRUE, auto_unbox = TRUE), "\n")
