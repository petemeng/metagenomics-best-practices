#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
  library(splines)
  library(timeROC)
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

SEED <- 68001L
BOOTSTRAP_ITERATIONS <- 1000L
CV_REPEATS <- 20L
CV_FOLDS <- 5L
ROC_TIMES <- c(180, 365, 548)
set.seed(SEED)
Sys.setenv(TZ = "UTC")

write_tsv <- function(x, name) {
  write.table(
    x,
    file.path(output_dir, name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

near_time_survival <- function(fit, time) {
  summary(fit, times = time, extend = TRUE)$surv[[1L]]
}

prepare_factors <- function(data) {
  data$PrimarySubtype <- factor(
    data$PrimarySubtype,
    levels = c("Cutaneous_or_unknown", "Mucosal_or_acral")
  )
  data$AdvancedSubstage <- factor(
    data$AdvancedSubstage,
    levels = c("Stage_M1C", "Stage_M1D")
  )
  data$LDH <- factor(data$LDH, levels = c("No", "Yes"))
  data$Treatment <- factor(data$Treatment, levels = c("all ICB", "anti-PD1"))
  data
}

cohort <- read.delim(file.path(input_dir, "survival-cohort.tsv"), check.names = FALSE)
cohort <- prepare_factors(cohort)
cohort$QualitySensitivityPass <- tolower(as.character(cohort$QualitySensitivityPass)) == "true"
stopifnot(nrow(cohort) == 110L, sum(cohort$Event) == 61L, !anyDuplicated(cohort$SampleID))
cohort$BMIz <- as.numeric(scale(cohort$BMI))
cohort$FaecalibacteriumZ <- as.numeric(scale(cohort$FaecalibacteriumLog2))

survival_formula <- Surv(PFS_days, Event) ~ FaecalibacteriumLog2
clinical_formula <- Surv(PFS_days, Event) ~ PrimarySubtype + AdvancedSubstage + LDH + BMIz
adjusted_formula <- update(clinical_formula, . ~ . + FaecalibacteriumLog2)

fit_unadjusted <- coxph(survival_formula, data = cohort, x = TRUE, y = TRUE)
fit_clinical <- coxph(clinical_formula, data = cohort, x = TRUE, y = TRUE)
fit_adjusted <- coxph(adjusted_formula, data = cohort, x = TRUE, y = TRUE)
quality <- subset(cohort, QualitySensitivityPass)
fit_quality <- coxph(adjusted_formula, data = quality, x = TRUE, y = TRUE)
anti_pd1 <- subset(cohort, Treatment == "anti-PD1")
fit_anti_pd1 <- coxph(adjusted_formula, data = anti_pd1, x = TRUE, y = TRUE)

tidy_cox <- function(fit, model_name, cohort_name) {
  estimate <- summary(fit)$coefficients
  ci <- summary(fit)$conf.int
  data.frame(
    Model = model_name,
    Cohort = cohort_name,
    N = fit$n,
    Events = fit$nevent,
    Term = rownames(estimate),
    LogHazard = estimate[, "coef"],
    StandardError = estimate[, "se(coef)"],
    HazardRatio = ci[, "exp(coef)"],
    CILower = ci[, "lower .95"],
    CIUpper = ci[, "upper .95"],
    WaldZ = estimate[, "z"],
    PValue = estimate[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

estimates <- do.call(
  rbind,
  list(
    tidy_cox(fit_unadjusted, "Unadjusted", "All complete PFS"),
    tidy_cox(fit_adjusted, "Adjusted primary", "All complete PFS"),
    tidy_cox(fit_quality, "Sequencing-QC sensitivity", "Quality sensitivity subset"),
    tidy_cox(fit_anti_pd1, "Anti-PD1 sensitivity", "Anti-PD1 monotherapy")
  )
)
write_tsv(estimates, "cox-model-estimates.tsv")

lrt <- anova(fit_clinical, fit_adjusted, test = "LRT")
model_audit <- data.frame(
  Model = c("Unadjusted", "Clinical only", "Adjusted primary", "Sequencing-QC sensitivity", "Anti-PD1 sensitivity"),
  N = c(fit_unadjusted$n, fit_clinical$n, fit_adjusted$n, fit_quality$n, fit_anti_pd1$n),
  Events = c(fit_unadjusted$nevent, fit_clinical$nevent, fit_adjusted$nevent, fit_quality$nevent, fit_anti_pd1$nevent),
  Concordance = c(
    summary(fit_unadjusted)$concordance[[1L]],
    summary(fit_clinical)$concordance[[1L]],
    summary(fit_adjusted)$concordance[[1L]],
    summary(fit_quality)$concordance[[1L]],
    summary(fit_anti_pd1)$concordance[[1L]]
  ),
  LogLikelihood = c(
    fit_unadjusted$loglik[[2L]], fit_clinical$loglik[[2L]], fit_adjusted$loglik[[2L]],
    fit_quality$loglik[[2L]], fit_anti_pd1$loglik[[2L]]
  ),
  AIC = c(AIC(fit_unadjusted), AIC(fit_clinical), AIC(fit_adjusted), AIC(fit_quality), AIC(fit_anti_pd1)),
  stringsAsFactors = FALSE
)
write_tsv(model_audit, "cox-model-audit.tsv")
incremental <- data.frame(
  Comparison = "Clinical plus Faecalibacterium versus clinical only",
  DegreesFreedom = lrt[["Df"]][[2L]],
  LikelihoodRatio = lrt[["Chisq"]][[2L]],
  PValue = lrt[["Pr(>|Chi|)"]][[2L]],
  ApparentDeltaConcordance = summary(fit_adjusted)$concordance[[1L]] - summary(fit_clinical)$concordance[[1L]]
)
write_tsv(incremental, "incremental-model-test.tsv")

ph <- cox.zph(fit_adjusted, transform = "km")
ph_rho <- c(
  vapply(seq_len(ncol(ph$y)), function(index) cor(ph$x, ph$y[, index]), numeric(1L)),
  NA_real_
)
ph_table <- data.frame(
  Term = rownames(ph$table),
  Rho = ph_rho,
  ChiSquare = ph$table[, "chisq"],
  DegreesFreedom = ph$table[, "df"],
  PValue = ph$table[, "p"],
  stringsAsFactors = FALSE
)
write_tsv(ph_table, "proportional-hazards-tests.tsv")
ph_residual <- data.frame(
  EventTimeDays = ph$time,
  TransformedTime = ph$x,
  ScaledSchoenfeld = ph$y[, "FaecalibacteriumLog2"],
  stringsAsFactors = FALSE
)
write_tsv(ph_residual, "faecalibacterium-schoenfeld.tsv")

fit_spline <- coxph(
  Surv(PFS_days, Event) ~ ns(FaecalibacteriumLog2, df = 3) +
    PrimarySubtype + AdvancedSubstage + LDH + BMIz,
  data = cohort,
  x = TRUE,
  y = TRUE
)
spline_lrt <- anova(fit_adjusted, fit_spline, test = "LRT")
spline_audit <- data.frame(
  Comparison = "Three-df natural spline versus linear Faecalibacterium effect",
  DegreesFreedom = spline_lrt[["Df"]][[2L]],
  LikelihoodRatio = spline_lrt[["Chisq"]][[2L]],
  PValue = spline_lrt[["Pr(>|Chi|)"]][[2L]],
  stringsAsFactors = FALSE
)
write_tsv(spline_audit, "nonlinearity-test.tsv")

grid_values <- seq(
  quantile(cohort$FaecalibacteriumLog2, 0.05),
  quantile(cohort$FaecalibacteriumLog2, 0.95),
  length.out = 121
)
reference_value <- median(cohort$FaecalibacteriumLog2)
prediction_grid <- data.frame(
  FaecalibacteriumLog2 = grid_values,
  PrimarySubtype = factor("Cutaneous_or_unknown", levels = levels(cohort$PrimarySubtype)),
  AdvancedSubstage = factor("Stage_M1C", levels = levels(cohort$AdvancedSubstage)),
  LDH = factor("No", levels = levels(cohort$LDH)),
  BMIz = 0
)
prediction_reference <- prediction_grid[1, , drop = FALSE]
prediction_reference$FaecalibacteriumLog2 <- reference_value
x_grid <- model.matrix(fit_spline, data = prediction_grid)
x_reference <- model.matrix(fit_spline, data = prediction_reference)
contrast <- sweep(x_grid, 2, x_reference[1, ], "-")
log_hr <- as.vector(contrast %*% coef(fit_spline))
se_log_hr <- sqrt(rowSums((contrast %*% vcov(fit_spline)) * contrast))
spline_prediction <- data.frame(
  FaecalibacteriumLog2 = grid_values,
  FaecalibacteriumPPM = pmax(0, 2^grid_values - 25),
  ReferenceLog2 = reference_value,
  HazardRatio = exp(log_hr),
  CILower = exp(log_hr - 1.96 * se_log_hr),
  CIUpper = exp(log_hr + 1.96 * se_log_hr)
)
write_tsv(spline_prediction, "spline-effect-curve.tsv")

leave_one_out <- lapply(seq_len(nrow(cohort)), function(index) {
  fit <- coxph(adjusted_formula, data = cohort[-index, , drop = FALSE])
  coefficient <- summary(fit)$coefficients["FaecalibacteriumLog2", ]
  data.frame(
    OmittedSampleID = cohort$SampleID[[index]],
    OmittedEvent = cohort$Event[[index]],
    OmittedFaecalibacteriumPPM = cohort$FaecalibacteriumPPM[[index]],
    OmittedQualitySensitivityPass = cohort$QualitySensitivityPass[[index]],
    HazardRatio = exp(coefficient[["coef"]]),
    CILower = exp(coefficient[["coef"]] - 1.96 * coefficient[["se(coef)"]]),
    CIUpper = exp(coefficient[["coef"]] + 1.96 * coefficient[["se(coef)"]]),
    PValue = coefficient[["Pr(>|z|)"]],
    stringsAsFactors = FALSE
  )
})
leave_one_out <- do.call(rbind, leave_one_out)
write_tsv(leave_one_out, "leave-one-out-influence.tsv")

stratified_split <- function(data, fraction = 0.70, seed = SEED) {
  set.seed(seed)
  stratum <- interaction(data$Event, data$AdvancedSubstage, drop = TRUE)
  train <- rep(FALSE, nrow(data))
  for (level in levels(stratum)) {
    indices <- which(stratum == level)
    indices <- sample(indices, length(indices), replace = FALSE)
    number <- max(1L, min(length(indices) - 1L, round(length(indices) * fraction)))
    train[indices[seq_len(number)]] <- TRUE
  }
  train
}

train_flag <- stratified_split(cohort)
cohort$Split <- ifelse(train_flag, "Training", "Held-out test")
split_ledger <- cohort[
  , c("SampleID", "Split", "PFS_days", "Event", "AdvancedSubstage", "FaecalibacteriumLog2", "FaecalibacteriumPPM")
]
write_tsv(split_ledger, "cutoff-split-ledger.tsv")

cutoff_search <- function(data, label) {
  lower <- unname(quantile(data$FaecalibacteriumLog2, 0.20))
  upper <- unname(quantile(data$FaecalibacteriumLog2, 0.80))
  candidates <- sort(unique(data$FaecalibacteriumLog2[data$FaecalibacteriumLog2 >= lower & data$FaecalibacteriumLog2 <= upper]))
  minimum_group <- max(10L, ceiling(0.15 * nrow(data)))
  rows <- lapply(candidates, function(cutoff) {
    group <- ifelse(data$FaecalibacteriumLog2 >= cutoff, "High", "Low")
    counts <- table(group)
    if (length(counts) != 2L || min(counts) < minimum_group) return(NULL)
    fit <- survdiff(Surv(PFS_days, Event) ~ group, data = data)
    p_value <- pchisq(fit$chisq, df = 1, lower.tail = FALSE)
    data.frame(
      SearchData = label,
      CutoffLog2 = cutoff,
      CutoffPPM = pmax(0, 2^cutoff - 25),
      LowN = unname(counts[["Low"]]),
      HighN = unname(counts[["High"]]),
      ChiSquare = unname(fit$chisq),
      PValue = p_value
    )
  })
  result <- do.call(rbind, rows)
  result <- result[order(result$PValue, result$CutoffLog2), , drop = FALSE]
  result$Selected <- seq_len(nrow(result)) == 1L
  result
}

training <- cohort[train_flag, , drop = FALSE]
test <- cohort[!train_flag, , drop = FALSE]
search_training <- cutoff_search(training, "Training only")
search_full <- cutoff_search(cohort, "Full data (leaky)")
write_tsv(search_training[order(search_training$CutoffLog2), ], "cutoff-search-training.tsv")
write_tsv(search_full[order(search_full$CutoffLog2), ], "cutoff-search-full-leaky.tsv")
cutoff_training <- search_training$CutoffLog2[[1L]]
cutoff_full <- search_full$CutoffLog2[[1L]]

evaluate_cutoff <- function(data, cutoff, evaluation, cutoff_source) {
  group <- factor(ifelse(data$FaecalibacteriumLog2 >= cutoff, "High", "Low"), levels = c("Low", "High"))
  fit <- coxph(Surv(PFS_days, Event) ~ group, data = data)
  logrank <- survdiff(Surv(PFS_days, Event) ~ group, data = data)
  coefficient <- summary(fit)$coefficients["groupHigh", ]
  counts <- table(group)
  events <- tapply(data$Event, group, sum)
  data.frame(
    EvaluationData = evaluation,
    CutoffSource = cutoff_source,
    CutoffLog2 = cutoff,
    CutoffPPM = pmax(0, 2^cutoff - 25),
    LowN = unname(counts[["Low"]]),
    LowEvents = unname(events[["Low"]]),
    HighN = unname(counts[["High"]]),
    HighEvents = unname(events[["High"]]),
    HazardRatioHighVsLow = exp(coefficient[["coef"]]),
    CILower = exp(coefficient[["coef"]] - 1.96 * coefficient[["se(coef)"]]),
    CIUpper = exp(coefficient[["coef"]] + 1.96 * coefficient[["se(coef)"]]),
    CoxPValue = coefficient[["Pr(>|z|)"]],
    LogRankPValue = pchisq(logrank$chisq, 1, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

cutoff_audit <- rbind(
  evaluate_cutoff(training, cutoff_training, "Training", "Training"),
  evaluate_cutoff(test, cutoff_training, "Held-out test", "Training"),
  evaluate_cutoff(cohort, cutoff_full, "Full data", "Full data (leaky)")
)
write_tsv(cutoff_audit, "cutoff-evaluation.tsv")

survfit_frame <- function(data, cutoff, evaluation, cutoff_source) {
  data$CutoffGroup <- factor(
    ifelse(data$FaecalibacteriumLog2 >= cutoff, "High", "Low"),
    levels = c("Low", "High")
  )
  fit <- survfit(Surv(PFS_days, Event) ~ CutoffGroup, data = data, conf.type = "log-log")
  strata_names <- rep(names(fit$strata), fit$strata)
  result <- data.frame(
    EvaluationData = evaluation,
    CutoffSource = cutoff_source,
    Group = sub("CutoffGroup=", "", strata_names, fixed = TRUE),
    TimeDays = fit$time,
    TimeMonths = fit$time / 30.4375,
    Survival = fit$surv,
    CILower = fit$lower,
    CIUpper = fit$upper,
    AtRisk = fit$n.risk,
    Events = fit$n.event,
    Censored = fit$n.censor,
    stringsAsFactors = FALSE
  )
  starts <- data.frame(
    EvaluationData = evaluation,
    CutoffSource = cutoff_source,
    Group = levels(data$CutoffGroup),
    TimeDays = 0,
    TimeMonths = 0,
    Survival = 1,
    CILower = 1,
    CIUpper = 1,
    AtRisk = as.integer(table(data$CutoffGroup)[levels(data$CutoffGroup)]),
    Events = 0,
    Censored = 0,
    stringsAsFactors = FALSE
  )
  rbind(starts, result)
}

km_curves <- rbind(
  survfit_frame(test, cutoff_training, "Held-out test", "Training"),
  survfit_frame(cohort, cutoff_full, "Full data", "Full data (leaky)")
)
write_tsv(km_curves, "cutoff-km-curves.tsv")

risk_table <- function(data, cutoff, evaluation, cutoff_source) {
  data$CutoffGroup <- factor(ifelse(data$FaecalibacteriumLog2 >= cutoff, "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(PFS_days, Event) ~ CutoffGroup, data = data)
  times <- c(0, 180, 365, 548, 730, 1095)
  summary_fit <- summary(fit, times = times, extend = TRUE)
  data.frame(
    EvaluationData = evaluation,
    CutoffSource = cutoff_source,
    Group = sub("CutoffGroup=", "", summary_fit$strata, fixed = TRUE),
    TimeDays = summary_fit$time,
    TimeMonths = summary_fit$time / 30.4375,
    AtRisk = summary_fit$n.risk,
    Survival = summary_fit$surv,
    stringsAsFactors = FALSE
  )
}
write_tsv(
  rbind(
    risk_table(test, cutoff_training, "Held-out test", "Training"),
    risk_table(cohort, cutoff_full, "Full data", "Full data (leaky)")
  ),
  "cutoff-risk-table.tsv"
)

assign_folds <- function(data, repeat_index) {
  set.seed(SEED + repeat_index)
  strata <- interaction(data$Event, data$AdvancedSubstage, drop = TRUE)
  folds <- integer(nrow(data))
  for (level in levels(strata)) {
    indices <- which(strata == level)
    indices <- sample(indices, length(indices), replace = FALSE)
    assigned <- rep(seq_len(CV_FOLDS), length.out = length(indices))
    folds[indices] <- assigned
  }
  folds
}

baseline_hazard_at <- function(fit, times) {
  hazard <- basehaz(fit, centered = FALSE)
  vapply(times, function(time) {
    eligible <- hazard$time <= time
    if (!any(eligible)) 0 else tail(hazard$hazard[eligible], 1L)
  }, numeric(1L))
}

cv_rows <- list()
fold_rows <- list()
cursor <- 1L
for (repeat_index in seq_len(CV_REPEATS)) {
  folds <- assign_folds(cohort, repeat_index)
  for (fold in seq_len(CV_FOLDS)) {
    train_index <- which(folds != fold)
    test_index <- which(folds == fold)
    train_data <- cohort[train_index, , drop = FALSE]
    test_data <- cohort[test_index, , drop = FALSE]
    bmi_mean <- mean(train_data$BMI)
    bmi_sd <- sd(train_data$BMI)
    faec_mean <- mean(train_data$FaecalibacteriumLog2)
    faec_sd <- sd(train_data$FaecalibacteriumLog2)
    train_data$BMIcv <- (train_data$BMI - bmi_mean) / bmi_sd
    test_data$BMIcv <- (test_data$BMI - bmi_mean) / bmi_sd
    train_data$FaecCV <- (train_data$FaecalibacteriumLog2 - faec_mean) / faec_sd
    test_data$FaecCV <- (test_data$FaecalibacteriumLog2 - faec_mean) / faec_sd
    clinical_fit <- coxph(
      Surv(PFS_days, Event) ~ PrimarySubtype + AdvancedSubstage + LDH + BMIcv,
      data = train_data,
      x = TRUE,
      y = TRUE
    )
    microbiome_fit <- coxph(
      Surv(PFS_days, Event) ~ PrimarySubtype + AdvancedSubstage + LDH + BMIcv + FaecCV,
      data = train_data,
      x = TRUE,
      y = TRUE
    )
    clinical_lp <- as.numeric(predict(clinical_fit, newdata = test_data, type = "lp", reference = "zero"))
    microbiome_lp <- as.numeric(predict(microbiome_fit, newdata = test_data, type = "lp", reference = "zero"))
    clinical_hazard <- baseline_hazard_at(clinical_fit, ROC_TIMES)
    microbiome_hazard <- baseline_hazard_at(microbiome_fit, ROC_TIMES)
    clinical_risk <- sapply(clinical_hazard, function(hazard) 1 - exp(-hazard * exp(clinical_lp)))
    microbiome_risk <- sapply(microbiome_hazard, function(hazard) 1 - exp(-hazard * exp(microbiome_lp)))
    if (length(test_index) == 1L) {
      clinical_risk <- matrix(clinical_risk, nrow = 1L)
      microbiome_risk <- matrix(microbiome_risk, nrow = 1L)
    }
    cv_rows[[cursor]] <- data.frame(
      Repeat = repeat_index,
      Fold = fold,
      SampleID = test_data$SampleID,
      ClinicalLP = clinical_lp,
      MicrobiomeLP = microbiome_lp,
      ClinicalRisk180 = clinical_risk[, 1L],
      ClinicalRisk365 = clinical_risk[, 2L],
      ClinicalRisk548 = clinical_risk[, 3L],
      MicrobiomeRisk180 = microbiome_risk[, 1L],
      MicrobiomeRisk365 = microbiome_risk[, 2L],
      MicrobiomeRisk548 = microbiome_risk[, 3L],
      stringsAsFactors = FALSE
    )
    fold_rows[[cursor]] <- data.frame(
      Repeat = repeat_index,
      Fold = fold,
      TrainingN = length(train_index),
      TrainingEvents = sum(cohort$Event[train_index]),
      TestN = length(test_index),
      TestEvents = sum(cohort$Event[test_index]),
      TrainingSampleHash = paste(sort(cohort$SampleID[train_index]), collapse = ";"),
      TestSampleHash = paste(sort(cohort$SampleID[test_index]), collapse = ";"),
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L
  }
}
cv_raw <- do.call(rbind, cv_rows)
fold_audit <- do.call(rbind, fold_rows)
write_tsv(cv_raw, "repeated-cv-predictions-long.tsv")
write_tsv(fold_audit, "repeated-cv-fold-audit.tsv")

mean_columns <- setdiff(names(cv_raw), c("Repeat", "Fold", "SampleID"))
cv_prediction <- aggregate(cv_raw[, mean_columns], list(SampleID = cv_raw$SampleID), mean)
cv_prediction <- merge(
  cohort[, c("SampleID", "PFS_days", "Event", "AdvancedSubstage")],
  cv_prediction,
  by = "SampleID",
  sort = FALSE
)
cv_prediction <- cv_prediction[match(cohort$SampleID, cv_prediction$SampleID), ]
stopifnot(all(cv_prediction$SampleID == cohort$SampleID))
write_tsv(cv_prediction, "repeated-cv-predictions.tsv")

harrell_c <- function(time, event, marker) {
  concordance(Surv(time, event) ~ marker, reverse = TRUE)$concordance
}

auc_at <- function(time, event, marker, horizon) {
  result <- tryCatch(
    timeROC(
      T = time,
      delta = event,
      marker = marker,
      cause = 1,
      weighting = "marginal",
      times = horizon,
      iid = FALSE
    ),
    error = function(error) NULL
  )
  if (is.null(result)) return(NA_real_)
  tail(result$AUC, 1L)
}

metric_point <- list(
  CindexClinical = harrell_c(cv_prediction$PFS_days, cv_prediction$Event, cv_prediction$ClinicalLP),
  CindexMicrobiome = harrell_c(cv_prediction$PFS_days, cv_prediction$Event, cv_prediction$MicrobiomeLP)
)
for (time_index in seq_along(ROC_TIMES)) {
  horizon <- ROC_TIMES[[time_index]]
  metric_point[[paste0("AUC", horizon, "Clinical")]] <- auc_at(
    cv_prediction$PFS_days, cv_prediction$Event,
    cv_prediction[[paste0("ClinicalRisk", horizon)]], horizon
  )
  metric_point[[paste0("AUC", horizon, "Microbiome")]] <- auc_at(
    cv_prediction$PFS_days, cv_prediction$Event,
    cv_prediction[[paste0("MicrobiomeRisk", horizon)]], horizon
  )
}

set.seed(SEED + 7000L)
bootstrap_metrics <- vector("list", BOOTSTRAP_ITERATIONS)
for (iteration in seq_len(BOOTSTRAP_ITERATIONS)) {
  indices <- sample(seq_len(nrow(cv_prediction)), nrow(cv_prediction), replace = TRUE)
  sample_data <- cv_prediction[indices, , drop = FALSE]
  row <- data.frame(
    Iteration = iteration,
    CindexClinical = harrell_c(sample_data$PFS_days, sample_data$Event, sample_data$ClinicalLP),
    CindexMicrobiome = harrell_c(sample_data$PFS_days, sample_data$Event, sample_data$MicrobiomeLP)
  )
  for (horizon in ROC_TIMES) {
    row[[paste0("AUC", horizon, "Clinical")]] <- auc_at(
      sample_data$PFS_days, sample_data$Event,
      sample_data[[paste0("ClinicalRisk", horizon)]], horizon
    )
    row[[paste0("AUC", horizon, "Microbiome")]] <- auc_at(
      sample_data$PFS_days, sample_data$Event,
      sample_data[[paste0("MicrobiomeRisk", horizon)]], horizon
    )
  }
  bootstrap_metrics[[iteration]] <- row
}
bootstrap_metrics <- do.call(rbind, bootstrap_metrics)
for (metric in c("Cindex", paste0("AUC", ROC_TIMES))) {
  bootstrap_metrics[[paste0(metric, "Delta")]] <-
    bootstrap_metrics[[paste0(metric, "Microbiome")]] - bootstrap_metrics[[paste0(metric, "Clinical")]]
}
write_tsv(bootstrap_metrics, "performance-bootstrap.tsv")

summary_rows <- list()
cursor <- 1L
summarize_metric <- function(values, estimate) {
  values <- values[is.finite(values)]
  c(Estimate = unname(estimate), CILower = unname(quantile(values, 0.025)), CIUpper = unname(quantile(values, 0.975)), ValidBootstrap = length(values))
}
for (metric in c("Cindex", paste0("AUC", ROC_TIMES))) {
  horizon <- if (metric == "Cindex") NA_integer_ else as.integer(sub("AUC", "", metric))
  clinical_name <- paste0(metric, "Clinical")
  microbiome_name <- paste0(metric, "Microbiome")
  clinical_point <- metric_point[[clinical_name]]
  microbiome_point <- metric_point[[microbiome_name]]
  for (model in c("Clinical", "Clinical + Faecalibacterium", "Increment")) {
    if (model == "Clinical") {
      values <- bootstrap_metrics[[clinical_name]]
      estimate <- clinical_point
    } else if (model == "Clinical + Faecalibacterium") {
      values <- bootstrap_metrics[[microbiome_name]]
      estimate <- microbiome_point
    } else {
      values <- bootstrap_metrics[[paste0(metric, "Delta")]]
      estimate <- microbiome_point - clinical_point
    }
    interval <- summarize_metric(values, estimate)
    summary_rows[[cursor]] <- data.frame(
      Metric = ifelse(metric == "Cindex", "Harrell C-index", "Cumulative/dynamic AUC"),
      HorizonDays = horizon,
      Model = model,
      Estimate = interval[["Estimate"]],
      CILower = interval[["CILower"]],
      CIUpper = interval[["CIUpper"]],
      ValidBootstrap = interval[["ValidBootstrap"]],
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L
  }
}
performance <- do.call(rbind, summary_rows)
write_tsv(performance, "cv-performance-summary.tsv")

roc_curve <- function(model, marker) {
  fit <- timeROC(
    T = cv_prediction$PFS_days,
    delta = cv_prediction$Event,
    marker = marker,
    cause = 1,
    weighting = "marginal",
    times = 365,
    iid = FALSE
  )
  column <- ncol(fit$TP)
  data.frame(
    Model = model,
    HorizonDays = 365,
    FalsePositiveRate = fit$FP[, column],
    TruePositiveRate = fit$TP[, column],
    AUC = tail(fit$AUC, 1L),
    stringsAsFactors = FALSE
  )
}
roc_curves <- rbind(
  roc_curve("Clinical", cv_prediction$ClinicalRisk365),
  roc_curve("Clinical + Faecalibacterium", cv_prediction$MicrobiomeRisk365)
)
write_tsv(roc_curves, "time-dependent-roc-365d.tsv")

calibration_rows <- list()
cursor <- 1L
for (model in c("Clinical", "Clinical + Faecalibacterium")) {
  risk_column <- if (model == "Clinical") "ClinicalRisk365" else "MicrobiomeRisk365"
  breaks <- unique(quantile(cv_prediction[[risk_column]], seq(0, 1, length.out = 6), type = 2))
  if (length(breaks) != 6L) stop("Calibration quintile boundaries are not unique")
  groups <- cut(cv_prediction[[risk_column]], breaks = breaks, include.lowest = TRUE, labels = FALSE)
  for (group in sort(unique(groups))) {
    indices <- which(groups == group)
    fit <- survfit(Surv(PFS_days, Event) ~ 1, data = cv_prediction[indices, , drop = FALSE], conf.type = "log-log")
    result <- summary(fit, times = 365, extend = TRUE)
    calibration_rows[[cursor]] <- data.frame(
      Model = model,
      RiskQuintile = group,
      N = length(indices),
      Events = sum(cv_prediction$Event[indices]),
      MeanPredictedRisk = mean(cv_prediction[[risk_column]][indices]),
      ObservedRiskKM = 1 - result$surv[[1L]],
      ObservedRiskLower = 1 - result$upper[[1L]],
      ObservedRiskUpper = 1 - result$lower[[1L]],
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L
  }
}
calibration <- do.call(rbind, calibration_rows)
write_tsv(calibration, "calibration-365d.tsv")

saveRDS(fit_adjusted, file.path(output_dir, "adjusted-cox-model.rds"))
saveRDS(fit_spline, file.path(output_dir, "spline-cox-model.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "r-session-info.txt"))
versions <- data.frame(
  Package = c("R", "survival", "timeROC"),
  Version = c(as.character(getRversion()), as.character(packageVersion("survival")), as.character(packageVersion("timeROC"))),
  stringsAsFactors = FALSE
)
write_tsv(versions, "software-versions-r.tsv")

faec_unadjusted <- subset(estimates, Model == "Unadjusted" & Term == "FaecalibacteriumLog2")
faec_adjusted <- subset(estimates, Model == "Adjusted primary" & Term == "FaecalibacteriumLog2")
faec_quality <- subset(estimates, Model == "Sequencing-QC sensitivity" & Term == "FaecalibacteriumLog2")
test_cutoff <- subset(cutoff_audit, EvaluationData == "Held-out test")
model_metrics <- list(
  article = 68L,
  analysis_seed = SEED,
  bootstrap_iterations = BOOTSTRAP_ITERATIONS,
  cv_repeats = CV_REPEATS,
  cv_folds = CV_FOLDS,
  train_samples = sum(train_flag),
  test_samples = sum(!train_flag),
  train_events = sum(cohort$Event[train_flag]),
  test_events = sum(cohort$Event[!train_flag]),
  unadjusted_faecalibacterium_hr_per_doubling = faec_unadjusted$HazardRatio,
  unadjusted_faecalibacterium_p = faec_unadjusted$PValue,
  adjusted_faecalibacterium_hr_per_doubling = faec_adjusted$HazardRatio,
  adjusted_faecalibacterium_p = faec_adjusted$PValue,
  quality_sensitivity_faecalibacterium_hr_per_doubling = faec_quality$HazardRatio,
  quality_sensitivity_faecalibacterium_p = faec_quality$PValue,
  global_ph_p = ph_table$PValue[ph_table$Term == "GLOBAL"],
  nonlinearity_p = spline_audit$PValue,
  training_cutoff_log2 = cutoff_training,
  training_cutoff_ppm = pmax(0, 2^cutoff_training - 25),
  leaky_full_cutoff_log2 = cutoff_full,
  heldout_cutoff_hr = test_cutoff$HazardRatioHighVsLow,
  heldout_cutoff_p = test_cutoff$LogRankPValue,
  cv_cindex_clinical = metric_point$CindexClinical,
  cv_cindex_microbiome = metric_point$CindexMicrobiome,
  cv_auc365_clinical = metric_point$AUC365Clinical,
  cv_auc365_microbiome = metric_point$AUC365Microbiome
)
jsonlite::write_json(
  model_metrics,
  file.path(output_dir, "model-metrics.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
cat(jsonlite::toJSON(model_metrics, pretty = TRUE, auto_unbox = TRUE, digits = 8), "\n")
