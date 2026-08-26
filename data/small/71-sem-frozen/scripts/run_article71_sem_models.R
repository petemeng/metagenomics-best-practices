#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(piecewiseSEM)
  library(sandwich)
  library(lmtest)
  library(car)
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

SEED <- 71001L
PLOT_SEED <- 20260771L
BOOTSTRAP <- 5000L
SENSITIVITY_BOOTSTRAP <- 2000L
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

primary <- read.delim(
  file.path(input_dir, "sem-primary-cohort.tsv"),
  check.names = FALSE
)
validation <- read.delim(
  file.path(input_dir, "sem-validation-cohort.tsv"),
  check.names = FALSE
)
stopifnot(
  nrow(primary) == 90L,
  sum(primary$Antibiotic) == 13L,
  nrow(validation) == 38L,
  sum(validation$Antibiotic) == 0L,
  packageVersion("piecewiseSEM") == "2.3.0.1"
)

covariates <- c(
  "CD", "UC", "AgeZ", "Immunosuppressant", "Mesalamine", "Steroids"
)
formula_mediator <- as.formula(
  paste("ShannonZ ~ Antibiotic +", paste(covariates, collapse = " + "))
)
formula_outcome <- as.formula(
  paste(
    "LogCalprotectinZ ~ Antibiotic + ShannonZ +",
    paste(covariates, collapse = " + ")
  )
)
formula_total <- as.formula(
  paste("LogCalprotectinZ ~ Antibiotic +", paste(covariates, collapse = " + "))
)
formula_constrained <- as.formula(
  paste("LogCalprotectinZ ~ ShannonZ +", paste(covariates, collapse = " + "))
)
formula_reverse_m <- as.formula(
  paste(
    "ShannonZ ~ Antibiotic + LogCalprotectinZ +",
    paste(covariates, collapse = " + ")
  )
)

mediator_fit <- lm(formula_mediator, data = primary)
outcome_fit <- lm(formula_outcome, data = primary)
total_fit <- lm(formula_total, data = primary)
constrained_outcome_fit <- lm(formula_constrained, data = primary)
reverse_mediator_fit <- lm(formula_reverse_m, data = primary)

primary_sem <- psem(mediator_fit, outcome_fit, data = primary)
constrained_sem <- psem(
  mediator_fit,
  constrained_outcome_fit,
  data = primary
)
reverse_sem <- psem(total_fit, reverse_mediator_fit, data = primary)

primary_aic <- AIC(primary_sem)
constrained_aic <- AIC(constrained_sem)
reverse_aic <- AIC(reverse_sem)
constrained_dsep <- dSep(constrained_sem)
constrained_fisher <- fisherC(constrained_dsep)

model_fit <- data.frame(
  Model = c(
    "Primary partial-path model",
    "Constrained microbiome-only path",
    "Reverse cross-sectional orientation"
  ),
  Direction = c(
    "Antibiotic -> Shannon -> calprotectin; antibiotic -> calprotectin",
    "Antibiotic -> Shannon -> calprotectin",
    "Antibiotic -> calprotectin -> Shannon; antibiotic -> Shannon"
  ),
  AIC = c(primary_aic$AIC, constrained_aic$AIC, reverse_aic$AIC),
  Parameters = c(primary_aic$K, constrained_aic$K, reverse_aic$K),
  N = c(primary_aic$n, constrained_aic$n, reverse_aic$n),
  IndependenceClaims = c(0L, nrow(constrained_dsep), 0L),
  FisherC = c(NA_real_, constrained_fisher$Fisher.C, NA_real_),
  FisherDF = c(0L, constrained_fisher$df, 0L),
  FisherP = c(NA_real_, constrained_fisher$P.Value, NA_real_),
  Saturated = c(TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
model_fit$DeltaAIC <- model_fit$AIC - min(model_fit$AIC)
write_tsv(model_fit, "sem-fit-comparison.tsv")

dsep_table <- data.frame(
  Model = "Constrained microbiome-only path",
  IndependenceClaim = constrained_dsep$Independ.Claim,
  TestType = constrained_dsep$Test.Type,
  DF = constrained_dsep$DF,
  Estimate = constrained_dsep$Crit.Value,
  PValue = constrained_dsep$P.Value
)
write_tsv(dsep_table, "directed-separation-claims.tsv")

hc3_table <- function(model, model_name) {
  tested <- coeftest(model, vcov. = vcovHC(model, type = "HC3"))
  data.frame(
    Model = model_name,
    Term = rownames(tested),
    Estimate = tested[, 1L],
    RobustSE = tested[, 2L],
    TStatistic = tested[, 3L],
    PValue = tested[, 4L],
    CILower = tested[, 1L] - 1.96 * tested[, 2L],
    CIUpper = tested[, 1L] + 1.96 * tested[, 2L],
    stringsAsFactors = FALSE
  )
}
local_coefficients <- rbind(
  hc3_table(mediator_fit, "Microbiome node"),
  hc3_table(outcome_fit, "Phenotype node"),
  hc3_table(total_fit, "Total-association node"),
  hc3_table(constrained_outcome_fit, "Constrained phenotype node"),
  hc3_table(reverse_mediator_fit, "Reverse microbiome node")
)
write_tsv(local_coefficients, "local-path-coefficients-hc3.tsv")

matrix_contract <- function(model) {
  list(
    x = model.matrix(model),
    y = model.response(model.frame(model)),
    terms = colnames(model.matrix(model))
  )
}
fit_matrix <- function(contract, index) {
  fit <- .lm.fit(
    x = contract$x[index, , drop = FALSE],
    y = contract$y[index]
  )
  if (fit$rank != ncol(contract$x) || any(!is.finite(fit$coefficients))) {
    return(rep(NA_real_, ncol(contract$x)))
  }
  names(fit$coefficients) <- contract$terms
  fit$coefficients
}

bootstrap_paths <- function(
  mediator_model,
  outcome_model,
  total_model,
  repetitions,
  seed
) {
  set.seed(seed)
  mediator_contract <- matrix_contract(mediator_model)
  outcome_contract <- matrix_contract(outcome_model)
  total_contract <- matrix_contract(total_model)
  n <- nrow(mediator_contract$x)
  result <- matrix(
    NA_real_,
    nrow = repetitions,
    ncol = 7L,
    dimnames = list(
      NULL,
      c("A", "B", "Direct", "Indirect", "Total", "PathSum", "IdentityError")
    )
  )
  for (iteration in seq_len(repetitions)) {
    index <- sample.int(n, size = n, replace = TRUE)
    mediator_coef <- fit_matrix(mediator_contract, index)
    outcome_coef <- fit_matrix(outcome_contract, index)
    total_coef <- fit_matrix(total_contract, index)
    if (
      any(!is.finite(mediator_coef)) ||
      any(!is.finite(outcome_coef)) ||
      any(!is.finite(total_coef))
    ) {
      next
    }
    a <- mediator_coef[["Antibiotic"]]
    mediator_term <- setdiff(
      intersect(c("ShannonZ", "LogFaecalibacteriumZ"), names(outcome_coef)),
      character()
    )
    if (length(mediator_term) != 1L) stop("Cannot identify mediator term")
    b <- outcome_coef[[mediator_term]]
    direct <- outcome_coef[["Antibiotic"]]
    indirect <- a * b
    total <- total_coef[["Antibiotic"]]
    path_sum <- direct + indirect
    result[iteration, ] <- c(
      a, b, direct, indirect, total, path_sum, total - path_sum
    )
  }
  data.frame(
    Iteration = seq_len(repetitions),
    as.data.frame(result),
    check.names = FALSE
  )
}

primary_bootstrap <- bootstrap_paths(
  mediator_fit,
  outcome_fit,
  total_fit,
  BOOTSTRAP,
  SEED + 1000L
)
write_tsv(primary_bootstrap, "sem-path-bootstrap.tsv.gz")

bootstrap_p <- function(values) {
  values <- values[is.finite(values)]
  min(1, 2 * min(mean(values <= 0), mean(values >= 0)))
}
effect_summary <- function(
  draws,
  points,
  model_name,
  bootstrap_repetitions
) {
  effects <- c("A", "B", "Direct", "Indirect", "Total")
  rows <- lapply(effects, function(effect) {
    values <- draws[[effect]]
    valid <- values[is.finite(values)]
    data.frame(
      Model = model_name,
      Effect = effect,
      Estimate = unname(points[[effect]]),
      CILower = unname(quantile(valid, 0.025)),
      CIUpper = unname(quantile(valid, 0.975)),
      BootstrapP = bootstrap_p(valid),
      ValidBootstrap = length(valid),
      RequestedBootstrap = bootstrap_repetitions
    )
  })
  do.call(rbind, rows)
}

primary_points <- c(
  A = coef(mediator_fit)[["Antibiotic"]],
  B = coef(outcome_fit)[["ShannonZ"]],
  Direct = coef(outcome_fit)[["Antibiotic"]]
)
primary_points[["Indirect"]] <- primary_points[["A"]] * primary_points[["B"]]
primary_points[["Total"]] <- coef(total_fit)[["Antibiotic"]]
stopifnot(
  abs(
    primary_points[["Total"]] -
      primary_points[["Direct"]] -
      primary_points[["Indirect"]]
  ) < 1e-10
)
primary_effects <- effect_summary(
  primary_bootstrap,
  primary_points,
  "Primary Shannon path",
  BOOTSTRAP
)

formula_faec_m <- as.formula(
  paste(
    "LogFaecalibacteriumZ ~ Antibiotic +",
    paste(covariates, collapse = " + ")
  )
)
formula_faec_y <- as.formula(
  paste(
    "LogCalprotectinZ ~ Antibiotic + LogFaecalibacteriumZ +",
    paste(covariates, collapse = " + ")
  )
)
faec_mediator_fit <- lm(formula_faec_m, data = primary)
faec_outcome_fit <- lm(formula_faec_y, data = primary)
faec_bootstrap <- bootstrap_paths(
  faec_mediator_fit,
  faec_outcome_fit,
  total_fit,
  SENSITIVITY_BOOTSTRAP,
  SEED + 2000L
)
write_tsv(faec_bootstrap, "faecalibacterium-path-bootstrap.tsv.gz")
faec_points <- c(
  A = coef(faec_mediator_fit)[["Antibiotic"]],
  B = coef(faec_outcome_fit)[["LogFaecalibacteriumZ"]],
  Direct = coef(faec_outcome_fit)[["Antibiotic"]]
)
faec_points[["Indirect"]] <- faec_points[["A"]] * faec_points[["B"]]
faec_points[["Total"]] <- coef(total_fit)[["Antibiotic"]]
faec_effects <- effect_summary(
  faec_bootstrap,
  faec_points,
  "Faecalibacterium sensitivity",
  SENSITIVITY_BOOTSTRAP
)
path_effects <- rbind(primary_effects, faec_effects)
write_tsv(path_effects, "path-effect-summary.tsv")

leave_one_out <- lapply(seq_len(nrow(primary)), function(omitted) {
  data <- primary[-omitted, , drop = FALSE]
  m <- lm(formula_mediator, data = data)
  y <- lm(formula_outcome, data = data)
  total <- lm(formula_total, data = data)
  a <- coef(m)[["Antibiotic"]]
  b <- coef(y)[["ShannonZ"]]
  direct <- coef(y)[["Antibiotic"]]
  indirect <- a * b
  total_effect <- coef(total)[["Antibiotic"]]
  data.frame(
    OmittedSample = primary$Sample[omitted],
    OmittedDiagnosis = primary$Diagnosis[omitted],
    OmittedAntibiotic = primary$Antibiotic[omitted],
    A = a,
    B = b,
    Direct = direct,
    Indirect = indirect,
    Total = total_effect,
    IdentityError = total_effect - direct - indirect
  )
})
leave_one_out <- do.call(rbind, leave_one_out)
write_tsv(leave_one_out, "leave-one-out-paths.tsv")

transport_formula <- formula_constrained
transport_rows <- lapply(
  list(
    "PRISM same-variable model" = primary,
    "Validation same-variable model" = validation
  ),
  function(data) {
    fit <- lm(transport_formula, data = data)
    row <- hc3_table(fit, "transport")
    row <- row[row$Term == "ShannonZ", , drop = FALSE]
    data.frame(
      N = nrow(data),
      AntibioticExposed = sum(data$Antibiotic),
      Estimate = row$Estimate,
      RobustSE = row$RobustSE,
      CILower = row$CILower,
      CIUpper = row$CIUpper,
      PValue = row$PValue,
      R2 = summary(fit)$r.squared
    )
  }
)
transport <- do.call(rbind, transport_rows)
transport$CohortModel <- rownames(transport)
rownames(transport) <- NULL
transport <- transport[
  ,
  c(
    "CohortModel", "N", "AntibioticExposed", "Estimate",
    "RobustSE", "CILower", "CIUpper", "PValue", "R2"
  )
]
write_tsv(transport, "outcome-path-transport.tsv")

diagnosis_exposure <- as.data.frame.matrix(
  table(primary$Diagnosis, primary$Antibiotic)
)
names(diagnosis_exposure) <- c("Unexposed", "Exposed")
diagnosis_exposure$Diagnosis <- rownames(diagnosis_exposure)
rownames(diagnosis_exposure) <- NULL
diagnosis_exposure <- diagnosis_exposure[
  match(c("Control", "CD", "UC"), diagnosis_exposure$Diagnosis),
  c("Diagnosis", "Unexposed", "Exposed")
]
diagnosis_exposure$ExposureFraction <- with(
  diagnosis_exposure,
  Exposed / (Exposed + Unexposed)
)
write_tsv(diagnosis_exposure, "antibiotic-overlap-by-diagnosis.tsv")

propensity_warnings <- character()
propensity_fit <- withCallingHandlers(
  glm(
    Antibiotic ~ CD + UC + AgeZ + Immunosuppressant + Mesalamine + Steroids,
    data = primary,
    family = binomial()
  ),
  warning = function(warning) {
    propensity_warnings <<- c(propensity_warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  }
)
propensity <- predict(propensity_fit, type = "response")
illustrative_weight <- ifelse(
  primary$Antibiotic == 1,
  1 / propensity,
  1 / (1 - propensity)
)
propensity_audit <- data.frame(
  Quantity = c(
    "GLM converged",
    "Iterations",
    "Maximum absolute coefficient",
    "Control antibiotic exposed",
    "PRISM antibiotic exposed",
    "Validation antibiotic exposed",
    "Minimum fitted propensity",
    "Maximum fitted propensity",
    "Maximum illustrative unstabilized weight",
    "Warning count",
    "Weights used for inference"
  ),
  Value = c(
    propensity_fit$converged,
    propensity_fit$iter,
    max(abs(coef(propensity_fit))),
    sum(primary$Diagnosis == "Control" & primary$Antibiotic == 1),
    sum(primary$Antibiotic),
    sum(validation$Antibiotic),
    min(propensity),
    max(propensity),
    max(illustrative_weight),
    length(propensity_warnings),
    FALSE
  ),
  Interpretation = c(
    "Convergence does not repair structural absence in controls",
    "Binomial fitting iterations",
    "Large diagnosis coefficients indicate separation",
    "Structural positivity warning",
    "Exposure path is estimated from 13 subjects",
    "External exposure path is not estimable",
    "Diagnostic only",
    "Diagnostic only",
    "Not used because the propensity model is separated/unstable",
    paste(propensity_warnings, collapse = " | "),
    "Primary path models are covariate-adjusted regressions"
  )
)
write_tsv(propensity_audit, "propensity-positivity-audit.tsv")

diagnose_model <- function(model, model_name) {
  bp <- bptest(model)
  cook <- cooks.distance(model)
  data.frame(
    Model = model_name,
    N = nobs(model),
    R2 = summary(model)$r.squared,
    AdjustedR2 = summary(model)$adj.r.squared,
    ConditionNumber = kappa(model.matrix(model), exact = TRUE),
    BreuschPaganStatistic = unname(bp$statistic),
    BreuschPaganP = bp$p.value,
    MaximumCookDistance = max(cook),
    CookAbove4OverN = sum(cook > 4 / nobs(model)),
    MaximumAbsoluteResidual = max(abs(residuals(model)))
  )
}
diagnostics <- rbind(
  diagnose_model(mediator_fit, "Microbiome node"),
  diagnose_model(outcome_fit, "Phenotype node"),
  diagnose_model(total_fit, "Total-association node")
)
write_tsv(diagnostics, "local-model-diagnostics.tsv")

vif_rows <- lapply(
  list(
    "Microbiome node" = mediator_fit,
    "Phenotype node" = outcome_fit
  ),
  function(model) {
    values <- car::vif(model)
    data.frame(Term = names(values), VIF = unname(values))
  }
)
vif_table <- do.call(rbind, vif_rows)
vif_table$Model <- rep(names(vif_rows), vapply(vif_rows, nrow, integer(1L)))
rownames(vif_table) <- NULL
vif_table <- vif_table[, c("Model", "Term", "VIF")]
write_tsv(vif_table, "variance-inflation.tsv")

software <- data.frame(
  Package = c(
    "R", "piecewiseSEM", "sandwich", "lmtest", "car", "jsonlite"
  ),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("piecewiseSEM")),
    as.character(packageVersion("sandwich")),
    as.character(packageVersion("lmtest")),
    as.character(packageVersion("car")),
    as.character(packageVersion("jsonlite"))
  )
)
write_tsv(software, "software-versions-r.tsv")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "r-session-info.txt"))
saveRDS(
  list(
    mediator_fit = mediator_fit,
    outcome_fit = outcome_fit,
    total_fit = total_fit,
    constrained_outcome_fit = constrained_outcome_fit,
    reverse_mediator_fit = reverse_mediator_fit,
    faec_mediator_fit = faec_mediator_fit,
    faec_outcome_fit = faec_outcome_fit
  ),
  file.path(output_dir, "sem-model-objects.rds"),
  compress = "xz"
)

model_metrics <- list(
  article = 71L,
  analysis_seed = SEED,
  plot_seed = PLOT_SEED,
  primary_subjects = nrow(primary),
  primary_antibiotic_exposed = sum(primary$Antibiotic),
  bootstrap = BOOTSTRAP,
  bootstrap_valid = sum(complete.cases(primary_bootstrap)),
  shannon_a = unname(primary_points[["A"]]),
  shannon_b = unname(primary_points[["B"]]),
  shannon_direct = unname(primary_points[["Direct"]]),
  shannon_indirect = unname(primary_points[["Indirect"]]),
  shannon_total = unname(primary_points[["Total"]]),
  shannon_indirect_ci_lower = primary_effects$CILower[primary_effects$Effect == "Indirect"],
  shannon_indirect_ci_upper = primary_effects$CIUpper[primary_effects$Effect == "Indirect"],
  shannon_indirect_p = primary_effects$BootstrapP[primary_effects$Effect == "Indirect"],
  mediator_r2 = summary(mediator_fit)$r.squared,
  outcome_r2 = summary(outcome_fit)$r.squared,
  constrained_fisher_c = constrained_fisher$Fisher.C,
  constrained_fisher_df = constrained_fisher$df,
  constrained_fisher_p = constrained_fisher$P.Value,
  primary_aic = primary_aic$AIC,
  constrained_aic = constrained_aic$AIC,
  reverse_aic = reverse_aic$AIC,
  propensity_converged = propensity_fit$converged,
  propensity_max_abs_coefficient = max(abs(coef(propensity_fit))),
  validation_antibiotic_exposed = sum(validation$Antibiotic),
  leave_one_out_indirect_min = min(leave_one_out$Indirect),
  leave_one_out_indirect_max = max(leave_one_out$Indirect)
)
write_json(
  model_metrics,
  file.path(output_dir, "model-metrics.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
print(toJSON(model_metrics, pretty = TRUE, auto_unbox = TRUE, digits = 6))
