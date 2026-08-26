#!/usr/bin/env Rscript

# Fit the predeclared subject-aware mixed models for Article 67.
Sys.setenv(TZ = "UTC")
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  where <- match(flag, args)
  if (is.na(where) || where == length(args)) stop("Missing argument: ", flag)
  args[[where + 1]]
}

input_dir <- normalizePath(value_after("--input-dir"), mustWork = TRUE)
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(67001)

write_tsv <- function(x, name) {
  write.table(
    x,
    file = file.path(output_dir, name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

sample <- read.delim(file.path(input_dir, "sample-ledger.tsv"), check.names = FALSE)
feature_audit <- read.delim(
  file.path(input_dir, "species-feature-audit.tsv"),
  check.names = FALSE
)
clr <- read.delim(
  gzfile(file.path(input_dir, "selected-species-clr.tsv.gz")),
  check.names = FALSE
)

stopifnot(nrow(sample) == 1523L, length(unique(sample$SubjectID)) == 107L)
stopifnot(identical(sample$SampleID, clr$SampleID))
sample$Diagnosis <- factor(sample$Diagnosis, levels = c("Control", "CD", "UC"))
sample$Antibiotics <- factor(sample$Antibiotics, levels = c("No", "Yes"))
sample$SubjectID <- factor(sample$SubjectID)

fixed_terms <- paste(
  "DysbiosisScore ~ Diagnosis * WeekYearCentered +",
  "Antibiotics + Log10ReadsCentered"
)
primary_formula <- as.formula(paste(fixed_terms, "+ (1 + WeekYearCentered | SubjectID)"))
fit_control <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 200000),
  check.conv.singular = "ignore"
)
primary <- lmerTest::lmer(
  primary_formula,
  data = sample,
  REML = TRUE,
  control = fit_control
)

coefficient_matrix <- coef(summary(primary))
fixed_effects <- data.frame(
  Term = rownames(coefficient_matrix),
  Estimate = coefficient_matrix[, "Estimate"],
  SE = coefficient_matrix[, "Std. Error"],
  DF = coefficient_matrix[, "df"],
  TStatistic = coefficient_matrix[, "t value"],
  PValue = coefficient_matrix[, "Pr(>|t|)"],
  row.names = NULL,
  check.names = FALSE
)
fixed_effects$CILower <- fixed_effects$Estimate - 1.96 * fixed_effects$SE
fixed_effects$CIUpper <- fixed_effects$Estimate + 1.96 * fixed_effects$SE
write_tsv(fixed_effects, "primary-fixed-effects.tsv")

anova_primary <- as.data.frame(anova(primary, type = 3, ddf = "Kenward-Roger"))
anova_primary$Term <- rownames(anova_primary)
rownames(anova_primary) <- NULL
anova_primary <- anova_primary[, c("Term", setdiff(names(anova_primary), "Term"))]
write_tsv(anova_primary, "primary-type3-anova.tsv")

random_effects <- as.data.frame(VarCorr(primary))
write_tsv(random_effects, "primary-random-effects.tsv")

ml_intercept <- lmer(
  as.formula(paste(fixed_terms, "+ (1 | SubjectID)")),
  data = sample,
  REML = FALSE,
  control = fit_control
)
ml_slope <- lmer(
  primary_formula,
  data = sample,
  REML = FALSE,
  control = fit_control
)
comparison <- as.data.frame(anova(ml_intercept, ml_slope, refit = FALSE))
comparison$Model <- rownames(comparison)
rownames(comparison) <- NULL
comparison <- comparison[, c("Model", setdiff(names(comparison), "Model"))]
write_tsv(comparison, "random-slope-model-comparison.tsv")

newdata <- expand.grid(
  Diagnosis = factor(c("Control", "CD", "UC"), levels = levels(sample$Diagnosis)),
  Week = seq(0, 52, by = 2),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
newdata$Diagnosis <- factor(newdata$Diagnosis, levels = levels(sample$Diagnosis))
newdata$WeekYearCentered <- (newdata$Week - 26) / 52
newdata$Antibiotics <- factor("No", levels = levels(sample$Antibiotics))
newdata$Log10ReadsCentered <- 0
fixed_formula <- formula(lme4::nobars(formula(primary)))
design <- model.matrix(delete.response(terms(fixed_formula)), newdata)
beta <- fixef(primary)
design <- design[, names(beta), drop = FALSE]
newdata$Predicted <- as.numeric(design %*% beta)
prediction_se <- sqrt(diag(design %*% vcov(primary) %*% t(design)))
newdata$CILower <- newdata$Predicted - 1.96 * prediction_se
newdata$CIUpper <- newdata$Predicted + 1.96 * prediction_se
write_tsv(newdata, "primary-marginal-predictions.tsv")

diagnostics <- data.frame(
  SampleID = sample$SampleID,
  SubjectID = as.character(sample$SubjectID),
  Diagnosis = as.character(sample$Diagnosis),
  Week = sample$Week,
  Fitted = fitted(primary),
  Residual = residuals(primary),
  PearsonResidual = residuals(primary, type = "pearson"),
  stringsAsFactors = FALSE
)
write_tsv(diagnostics, "primary-model-diagnostics.tsv")

gradient <- primary@optinfo$derivs$gradient
messages <- primary@optinfo$conv$lme4$messages
model_audit <- data.frame(
  Model = c("Primary REML random intercept+slope", "ML random intercept", "ML random intercept+slope"),
  Observations = c(nobs(primary), nobs(ml_intercept), nobs(ml_slope)),
  Subjects = length(unique(sample$SubjectID)),
  Singular = c(isSingular(primary), isSingular(ml_intercept), isSingular(ml_slope)),
  MaxAbsGradient = c(
    if (is.null(gradient)) NA_real_ else max(abs(gradient)),
    if (is.null(ml_intercept@optinfo$derivs$gradient)) NA_real_ else max(abs(ml_intercept@optinfo$derivs$gradient)),
    if (is.null(ml_slope@optinfo$derivs$gradient)) NA_real_ else max(abs(ml_slope@optinfo$derivs$gradient))
  ),
  ConvergenceMessage = c(
    if (is.null(messages)) "none" else paste(messages, collapse = "; "),
    if (is.null(ml_intercept@optinfo$conv$lme4$messages)) "none" else paste(ml_intercept@optinfo$conv$lme4$messages, collapse = "; "),
    if (is.null(ml_slope@optinfo$conv$lme4$messages)) "none" else paste(ml_slope@optinfo$conv$lme4$messages, collapse = "; ")
  ),
  stringsAsFactors = FALSE
)
write_tsv(model_audit, "primary-model-audit.tsv")
saveRDS(primary, file.path(output_dir, "primary-dysbiosis-lmm.rds"))

selected_flag <- tolower(as.character(feature_audit$SelectedForMixedModel)) == "true"
selected_ids <- feature_audit$FeatureID[selected_flag]
stopifnot(length(selected_ids) == 63L)
selected_terms <- c(
  "DiagnosisCD",
  "DiagnosisUC",
  "AntibioticsYes",
  "DiagnosisCD:WeekYearCentered",
  "DiagnosisUC:WeekYearCentered"
)
secondary_records <- list()
secondary_audit <- list()
for (feature_id in selected_ids) {
  frame <- sample
  frame$FeatureCLR <- clr[[feature_id]]
  warning_text <- character()
  model <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        FeatureCLR ~ Diagnosis * WeekYearCentered + Antibiotics +
          Log10ReadsCentered + (1 | SubjectID),
        data = frame,
        REML = TRUE,
        control = fit_control
      ),
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  species <- feature_audit$Species[match(feature_id, feature_audit$FeatureID)]
  if (inherits(model, "error")) {
    secondary_audit[[length(secondary_audit) + 1L]] <- data.frame(
      FeatureID = feature_id,
      Species = species,
      Success = FALSE,
      Singular = NA,
      MaxAbsGradient = NA,
      Message = conditionMessage(model),
      stringsAsFactors = FALSE
    )
    next
  }
  matrix <- coef(summary(model))
  available <- intersect(selected_terms, rownames(matrix))
  for (term in available) {
    secondary_records[[length(secondary_records) + 1L]] <- data.frame(
      FeatureID = feature_id,
      Species = species,
      Term = term,
      Estimate = matrix[term, "Estimate"],
      SE = matrix[term, "Std. Error"],
      DF = matrix[term, "df"],
      TStatistic = matrix[term, "t value"],
      PValue = matrix[term, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }
  feature_gradient <- model@optinfo$derivs$gradient
  model_message <- model@optinfo$conv$lme4$messages
  secondary_audit[[length(secondary_audit) + 1L]] <- data.frame(
    FeatureID = feature_id,
    Species = species,
    Success = TRUE,
    Singular = isSingular(model),
    MaxAbsGradient = if (is.null(feature_gradient)) NA_real_ else max(abs(feature_gradient)),
    Message = paste(c(warning_text, model_message), collapse = "; "),
    stringsAsFactors = FALSE
  )
}

secondary <- do.call(rbind, secondary_records)
secondary$PValue <- as.numeric(secondary$PValue)
secondary$QValueWithinTerm <- NA_real_
for (term in unique(secondary$Term)) {
  term_rows <- which(secondary$Term == term)
  secondary$QValueWithinTerm[term_rows] <- p.adjust(
    secondary$PValue[term_rows], method = "BH"
  )
}
secondary$QValueGlobal <- p.adjust(secondary$PValue, method = "BH")
secondary$CILower <- secondary$Estimate - 1.96 * secondary$SE
secondary$CIUpper <- secondary$Estimate + 1.96 * secondary$SE
secondary <- secondary[order(secondary$QValueGlobal, secondary$PValue), ]
write_tsv(secondary, "species-mixed-model-results.tsv")
write_tsv(do.call(rbind, secondary_audit), "species-model-audit.tsv")

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "r-session-info.txt"))
versions <- data.frame(
  Package = c("R", "Matrix", "lme4", "lmerTest", "pbkrtest"),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("Matrix")),
    as.character(packageVersion("lme4")),
    as.character(packageVersion("lmerTest")),
    as.character(packageVersion("pbkrtest"))
  ),
  stringsAsFactors = FALSE
)
write_tsv(versions, "software-versions-r.tsv")
cat("primary rows", nrow(sample), "secondary models", length(selected_ids), "\n")
