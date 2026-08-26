#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(statmod)
  library(mediation)
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

SEED <- 69001L
PRIMARY_BOOTSTRAP <- 5000L
SENSITIVITY_BOOTSTRAP <- 2000L
QUADRATURE_NODES <- 20L
set.seed(SEED)
Sys.setenv(TZ = "UTC")

write_tsv <- function(x, name) {
  connection <- if (grepl("\\.gz$", name)) gzfile(file.path(output_dir, name), "wt") else file.path(output_dir, name)
  on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

cohort <- read.delim(file.path(input_dir, "mediation-cohort.tsv"), check.names = FALSE)
cohort$QualitySensitivityPass <- tolower(as.character(cohort$QualitySensitivityPass)) == "true"
cohort$A <- as.integer(cohort$ExposureSufficient)
cohort$Y <- as.integer(cohort$OutcomeResponder)
cohort$M <- cohort$FaecalibacteriumLog2
cohort$MRumi <- cohort$RuminococcaceaeLog2
cohort$BMIz <- (cohort$BMI - mean(cohort$BMI)) / sd(cohort$BMI)
cohort$Mucosal <- as.integer(cohort$PrimarySubtype == "Mucosal_or_acral")
cohort$StageM1D <- as.integer(cohort$AdvancedSubstage == "Stage_M1D")
cohort$LDHHigh <- as.integer(cohort$LDH == "Yes")
cohort$PrimarySubtype <- factor(cohort$PrimarySubtype, levels = c("Cutaneous_or_unknown", "Mucosal_or_acral"))
cohort$AdvancedSubstage <- factor(cohort$AdvancedSubstage, levels = c("Stage_M1C", "Stage_M1D"))
cohort$LDH <- factor(cohort$LDH, levels = c("No", "Yes"))
stopifnot(nrow(cohort) == 94L, sum(cohort$A) == 23L, sum(cohort$Y) == 60L)

mediator_fit <- lm(
  M ~ A + PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort
)
outcome_fit <- glm(
  Y ~ A + M + PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort,
  family = binomial()
)
total_fit <- glm(
  Y ~ A + PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort,
  family = binomial()
)
interaction_fit <- glm(
  Y ~ A * M + PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort,
  family = binomial()
)
continuous_mediator_fit <- lm(
  M ~ I(FiberGrams / 5) + PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort
)
continuous_outcome_fit <- glm(
  Y ~ I(FiberGrams / 5) + M + PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort,
  family = binomial()
)

tidy_lm <- function(fit, model, scale) {
  result <- summary(fit)$coefficients
  data.frame(
    Model = model,
    Scale = scale,
    Term = rownames(result),
    Estimate = result[, "Estimate"],
    StandardError = result[, "Std. Error"],
    CILower = result[, "Estimate"] - 1.96 * result[, "Std. Error"],
    CIUpper = result[, "Estimate"] + 1.96 * result[, "Std. Error"],
    TestStatistic = result[, "t value"],
    PValue = result[, "Pr(>|t|)"],
    OddsRatio = NA_real_,
    OddsRatioLower = NA_real_,
    OddsRatioUpper = NA_real_,
    stringsAsFactors = FALSE
  )
}

tidy_glm <- function(fit, model, scale) {
  result <- summary(fit)$coefficients
  data.frame(
    Model = model,
    Scale = scale,
    Term = rownames(result),
    Estimate = result[, "Estimate"],
    StandardError = result[, "Std. Error"],
    CILower = result[, "Estimate"] - 1.96 * result[, "Std. Error"],
    CIUpper = result[, "Estimate"] + 1.96 * result[, "Std. Error"],
    TestStatistic = result[, "z value"],
    PValue = result[, "Pr(>|z|)"],
    OddsRatio = exp(result[, "Estimate"]),
    OddsRatioLower = exp(result[, "Estimate"] - 1.96 * result[, "Std. Error"]),
    OddsRatioUpper = exp(result[, "Estimate"] + 1.96 * result[, "Std. Error"]),
    stringsAsFactors = FALSE
  )
}

path_estimates <- rbind(
  tidy_lm(mediator_fit, "Mediator model", "Mediator log2 PPM difference"),
  tidy_glm(outcome_fit, "Outcome model", "Log odds / odds ratio"),
  tidy_glm(total_fit, "Total-association model", "Log odds / odds ratio"),
  tidy_glm(interaction_fit, "Outcome interaction sensitivity", "Log odds / odds ratio"),
  tidy_lm(continuous_mediator_fit, "Continuous-fiber mediator model", "Mediator log2 PPM difference per 5 g/day"),
  tidy_glm(continuous_outcome_fit, "Continuous-fiber outcome model", "Log odds / odds ratio per 5 g/day")
)
write_tsv(path_estimates, "path-model-estimates.tsv")

propensity_fit <- glm(
  A ~ PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort,
  family = binomial()
)
cohort$Propensity <- predict(propensity_fit, type = "response")
cohort$IPTW <- ifelse(cohort$A == 1, 1 / cohort$Propensity, 1 / (1 - cohort$Propensity))
overlap <- cohort[, c("SampleID", "A", "FiberCategory", "Propensity", "IPTW", "BMIz", "Mucosal", "StageM1D", "LDHHigh")]
write_tsv(overlap, "exposure-overlap.tsv")
overlap_summary <- do.call(
  rbind,
  lapply(0:1, function(exposure) {
    values <- cohort$Propensity[cohort$A == exposure]
    data.frame(
      Exposure = exposure,
      Label = ifelse(exposure == 1, "Sufficient", "Insufficient"),
      N = length(values),
      Minimum = min(values),
      Q05 = quantile(values, 0.05),
      Median = median(values),
      Q95 = quantile(values, 0.95),
      Maximum = max(values),
      MaximumIPTW = max(cohort$IPTW[cohort$A == exposure])
    )
  })
)
write_tsv(overlap_summary, "exposure-overlap-summary.tsv")

weighted_mean <- function(x, w) sum(x * w) / sum(w)
weighted_var <- function(x, w) sum(w * (x - weighted_mean(x, w))^2) / sum(w)
smd_binary <- function(x, a, w = rep(1, length(x))) {
  mean1 <- weighted_mean(x[a == 1], w[a == 1])
  mean0 <- weighted_mean(x[a == 0], w[a == 0])
  variance <- (weighted_var(x[a == 1], w[a == 1]) + weighted_var(x[a == 0], w[a == 0])) / 2
  (mean1 - mean0) / sqrt(variance)
}
balance <- do.call(
  rbind,
  lapply(c("BMIz", "Mucosal", "StageM1D", "LDHHigh"), function(variable) {
    data.frame(
      Covariate = variable,
      SMDUnweighted = smd_binary(cohort[[variable]], cohort$A),
      SMDIPTW = smd_binary(cohort[[variable]], cohort$A, cohort$IPTW)
    )
  })
)
write_tsv(balance, "exposure-balance.tsv")

quadrature <- gauss.quad.prob(QUADRATURE_NODES, dist = "normal")
confounder_columns <- c("Mucosal", "StageM1D", "LDHHigh", "BMIz")

fast_gformula <- function(data, mediator_column = "M", interaction = FALSE, target = data, individual = FALSE) {
  confounders <- as.matrix(data[, confounder_columns, drop = FALSE])
  mediator_design <- cbind(1, data$A, confounders)
  mediator_model <- lm.fit(mediator_design, data[[mediator_column]])
  if (mediator_model$df.residual <= 0 || any(!is.finite(mediator_model$coefficients))) stop("Invalid mediator model")
  mediator_sigma <- sqrt(sum(mediator_model$residuals^2) / mediator_model$df.residual)
  if (!is.finite(mediator_sigma) || mediator_sigma <= 0) stop("Invalid mediator residual scale")

  if (interaction) {
    outcome_design <- cbind(1, data$A, data[[mediator_column]], data$A * data[[mediator_column]], confounders)
  } else {
    outcome_design <- cbind(1, data$A, data[[mediator_column]], confounders)
  }
  outcome_model <- suppressWarnings(glm.fit(outcome_design, data$Y, family = binomial()))
  if (any(!is.finite(outcome_model$coefficients)) || max(abs(outcome_model$coefficients)) > 30) stop("Unstable outcome model")

  target_confounders <- as.matrix(target[, confounder_columns, drop = FALSE])
  scenario <- function(outcome_exposure, mediator_exposure) {
    mediator_mean <- mediator_model$coefficients[[1L]] +
      mediator_model$coefficients[[2L]] * mediator_exposure +
      as.vector(target_confounders %*% mediator_model$coefficients[3:6])
    if (interaction) {
      base <- outcome_model$coefficients[[1L]] +
        outcome_model$coefficients[[2L]] * outcome_exposure +
        as.vector(target_confounders %*% outcome_model$coefficients[5:8])
      mediator_coefficient <- outcome_model$coefficients[[3L]] + outcome_model$coefficients[[4L]] * outcome_exposure
    } else {
      base <- outcome_model$coefficients[[1L]] +
        outcome_model$coefficients[[2L]] * outcome_exposure +
        as.vector(target_confounders %*% outcome_model$coefficients[4:7])
      mediator_coefficient <- outcome_model$coefficients[[3L]]
    }
    eta <- outer(base + mediator_coefficient * mediator_mean, rep(1, QUADRATURE_NODES)) +
      outer(rep(1, nrow(target)), mediator_coefficient * mediator_sigma * quadrature$nodes)
    probabilities <- as.vector(plogis(eta) %*% quadrature$weights)
    if (individual) probabilities else mean(probabilities)
  }

  p00 <- scenario(0, 0)
  p10 <- scenario(1, 0)
  p11 <- scenario(1, 1)
  p01 <- scenario(0, 1)
  if (individual) {
    return(data.frame(P00 = p00, P10 = p10, P11 = p11, P01 = p01, Direct = p10 - p00, Indirect = p11 - p10, Total = p11 - p00))
  }
  total <- p11 - p00
  c(
    P00 = p00,
    P10 = p10,
    P11 = p11,
    P01 = p01,
    Direct = p10 - p00,
    Indirect = p11 - p10,
    Total = total,
    ProportionMediated = ifelse(abs(total) < 1e-8, NA_real_, (p11 - p10) / total)
  )
}

bootstrap_effects <- function(data, mediator_column, interaction, iterations, seed, variant) {
  set.seed(seed)
  result <- matrix(NA_real_, nrow = iterations, ncol = 8)
  colnames(result) <- c("P00", "P10", "P11", "P01", "Direct", "Indirect", "Total", "ProportionMediated")
  for (iteration in seq_len(iterations)) {
    index <- sample(seq_len(nrow(data)), nrow(data), replace = TRUE)
    result[iteration, ] <- tryCatch(
      fast_gformula(data[index, , drop = FALSE], mediator_column, interaction),
      error = function(error) rep(NA_real_, 8)
    )
  }
  data.frame(Variant = variant, Iteration = seq_len(iterations), result, check.names = FALSE)
}

primary_point <- fast_gformula(cohort, "M", FALSE)
interaction_point <- fast_gformula(cohort, "M", TRUE)
quality_cohort <- cohort[cohort$QualitySensitivityPass, , drop = FALSE]
quality_point <- fast_gformula(quality_cohort, "M", FALSE)
rumi_point <- fast_gformula(cohort, "MRumi", FALSE)

primary_bootstrap <- bootstrap_effects(cohort, "M", FALSE, PRIMARY_BOOTSTRAP, SEED + 1000L, "Primary")
interaction_bootstrap <- bootstrap_effects(cohort, "M", TRUE, SENSITIVITY_BOOTSTRAP, SEED + 2000L, "Exposure-mediator interaction")
quality_bootstrap <- bootstrap_effects(quality_cohort, "M", FALSE, SENSITIVITY_BOOTSTRAP, SEED + 3000L, "Sequencing-QC subset")
rumi_bootstrap <- bootstrap_effects(cohort, "MRumi", FALSE, SENSITIVITY_BOOTSTRAP, SEED + 4000L, "Ruminococcaceae mediator")
write_tsv(primary_bootstrap, "primary-gformula-bootstrap.tsv.gz")
write_tsv(rbind(interaction_bootstrap, quality_bootstrap, rumi_bootstrap), "sensitivity-gformula-bootstrap.tsv.gz")

point_map <- list(
  Primary = primary_point,
  `Exposure-mediator interaction` = interaction_point,
  `Sequencing-QC subset` = quality_point,
  `Ruminococcaceae mediator` = rumi_point
)
bootstrap_map <- list(
  Primary = primary_bootstrap,
  `Exposure-mediator interaction` = interaction_bootstrap,
  `Sequencing-QC subset` = quality_bootstrap,
  `Ruminococcaceae mediator` = rumi_bootstrap
)
summary_rows <- list()
cursor <- 1L
for (variant in names(point_map)) {
  draws <- bootstrap_map[[variant]]
  for (effect in c("P00", "P10", "P11", "P01", "Direct", "Indirect", "Total", "ProportionMediated")) {
    values <- draws[[effect]][is.finite(draws[[effect]])]
    estimate <- point_map[[variant]][[effect]]
    summary_rows[[cursor]] <- data.frame(
      Variant = variant,
      Effect = effect,
      Estimate = estimate,
      CILower = unname(quantile(values, 0.025)),
      CIUpper = unname(quantile(values, 0.975)),
      BootstrapP = if (effect %in% c("Direct", "Indirect", "Total")) min(1, 2 * min(mean(values <= 0), mean(values >= 0))) else NA_real_,
      ValidBootstrap = length(values),
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L
  }
}
gformula_summary <- do.call(rbind, summary_rows)
write_tsv(gformula_summary, "gformula-effect-summary.tsv")

individual <- fast_gformula(cohort, "M", FALSE, individual = TRUE)
individual <- cbind(cohort[, c("SampleID", "A", "Y", "FiberCategory", "M")], individual)
write_tsv(individual, "individual-standardized-risks.tsv")

probit_mediator_fit <- lm(
  M ~ A + PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort
)
probit_outcome_fit <- glm(
  Y ~ A + M + PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort,
  family = binomial(link = "probit")
)
set.seed(SEED + 8000L)
mediate_fit <- mediate(
  probit_mediator_fit,
  probit_outcome_fit,
  treat = "A",
  mediator = "M",
  control.value = 0,
  treat.value = 1,
  sims = 5000,
  boot = FALSE
)
mediate_summary <- data.frame(
  Effect = c("ACME control", "ACME treated", "ACME average", "ADE average", "Total effect", "Proportion mediated average"),
  Estimate = c(mediate_fit$d0, mediate_fit$d1, mediate_fit$d.avg, mediate_fit$z.avg, mediate_fit$tau.coef, mediate_fit$n.avg),
  CILower = c(mediate_fit$d0.ci[1], mediate_fit$d1.ci[1], mediate_fit$d.avg.ci[1], mediate_fit$z.avg.ci[1], mediate_fit$tau.ci[1], mediate_fit$n.avg.ci[1]),
  CIUpper = c(mediate_fit$d0.ci[2], mediate_fit$d1.ci[2], mediate_fit$d.avg.ci[2], mediate_fit$z.avg.ci[2], mediate_fit$tau.ci[2], mediate_fit$n.avg.ci[2]),
  PValue = c(mediate_fit$d0.p, mediate_fit$d1.p, mediate_fit$d.avg.p, mediate_fit$z.avg.p, mediate_fit$tau.p, mediate_fit$n.avg.p)
)
write_tsv(mediate_summary, "probit-mediate-summary.tsv")

set.seed(SEED + 9000L)
sensitivity <- medsens(mediate_fit, rho.by = 0.05, effect.type = "indirect", sims = 1000)
sensitivity_rho <- data.frame(
  Rho = sensitivity$rho,
  ACMEControl = sensitivity$d0,
  ACMETreated = sensitivity$d1,
  ACMEAverage = (sensitivity$d0 + sensitivity$d1) / 2,
  LowerAverage = (sensitivity$lower.d0 + sensitivity$lower.d1) / 2,
  UpperAverage = (sensitivity$upper.d0 + sensitivity$upper.d1) / 2
)
write_tsv(sensitivity_rho, "residual-correlation-sensitivity.tsv")
sensitivity_audit <- data.frame(
  Quantity = c("Residual correlation where ACME crosses zero (control)", "Residual correlation where ACME crosses zero (treated)", "R2-star product threshold (control)", "R2-star product threshold (treated)", "R2-tilde product threshold (control)", "R2-tilde product threshold (treated)"),
  Value = c(sensitivity$err.cr.d[1], sensitivity$err.cr.d[2], sensitivity$R2star.d.thresh[1], sensitivity$R2star.d.thresh[2], sensitivity$R2tilde.d.thresh[1], sensitivity$R2tilde.d.thresh[2])
)
write_tsv(sensitivity_audit, "unmeasured-confounding-sensitivity-audit.tsv")

saveRDS(mediator_fit, file.path(output_dir, "mediator-model.rds"))
saveRDS(outcome_fit, file.path(output_dir, "outcome-model.rds"))
saveRDS(mediate_fit, file.path(output_dir, "probit-mediate-object.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "r-session-info.txt"))
versions <- data.frame(
  Package = c("R", "mediation", "statmod", "sandwich", "MASS", "Matrix", "mvtnorm"),
  Version = c(
    as.character(getRversion()),
    as.character(packageVersion("mediation")),
    as.character(packageVersion("statmod")),
    as.character(packageVersion("sandwich")),
    as.character(packageVersion("MASS")),
    as.character(packageVersion("Matrix")),
    as.character(packageVersion("mvtnorm"))
  )
)
write_tsv(versions, "software-versions-r.tsv")

mediator_a <- summary(mediator_fit)$coefficients["A", ]
outcome_m <- summary(outcome_fit)$coefficients["M", ]
outcome_a <- summary(outcome_fit)$coefficients["A", ]
total_a <- summary(total_fit)$coefficients["A", ]
primary_indirect <- subset(gformula_summary, Variant == "Primary" & Effect == "Indirect")
primary_direct <- subset(gformula_summary, Variant == "Primary" & Effect == "Direct")
primary_total <- subset(gformula_summary, Variant == "Primary" & Effect == "Total")
model_metrics <- list(
  article = 69L,
  analysis_seed = SEED,
  primary_bootstrap = PRIMARY_BOOTSTRAP,
  sensitivity_bootstrap = SENSITIVITY_BOOTSTRAP,
  quadrature_nodes = QUADRATURE_NODES,
  mediator_exposure_beta = unname(mediator_a[["Estimate"]]),
  mediator_exposure_p = unname(mediator_a[["Pr(>|t|)"]]),
  outcome_mediator_or = exp(unname(outcome_m[["Estimate"]])),
  outcome_mediator_p = unname(outcome_m[["Pr(>|z|)"]]),
  outcome_direct_exposure_or = exp(unname(outcome_a[["Estimate"]])),
  outcome_direct_exposure_p = unname(outcome_a[["Pr(>|z|)"]]),
  total_exposure_or = exp(unname(total_a[["Estimate"]])),
  total_exposure_p = unname(total_a[["Pr(>|z|)"]]),
  primary_direct_rd = primary_direct$Estimate,
  primary_direct_ci_lower = primary_direct$CILower,
  primary_direct_ci_upper = primary_direct$CIUpper,
  primary_indirect_rd = primary_indirect$Estimate,
  primary_indirect_ci_lower = primary_indirect$CILower,
  primary_indirect_ci_upper = primary_indirect$CIUpper,
  primary_total_rd = primary_total$Estimate,
  primary_total_ci_lower = primary_total$CILower,
  primary_total_ci_upper = primary_total$CIUpper,
  primary_proportion_mediated = primary_point[["ProportionMediated"]],
  propensity_min = min(cohort$Propensity),
  propensity_max = max(cohort$Propensity),
  maximum_iptw = max(cohort$IPTW),
  acme_zero_rho = sensitivity$err.cr.d[1]
)
jsonlite::write_json(model_metrics, file.path(output_dir, "model-metrics.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
cat(jsonlite::toJSON(model_metrics, pretty = TRUE, auto_unbox = TRUE, digits = 8), "\n")
