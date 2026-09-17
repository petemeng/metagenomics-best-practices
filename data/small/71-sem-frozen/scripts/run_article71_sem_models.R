#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sandwich)
  library(lmtest)
  library(car)
  library(jsonlite)
  suppressWarnings(library(mediation))
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
POWER_REPETITIONS <- 1000L
POSITIVE_CONTROL_BOOTSTRAP <- 2000L
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
  packageVersion("mediation") == "4.5.1"
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

# For this two-equation Gaussian SEM, the graph AIC is the sum of the local
# model AIC values. Writing the calculation explicitly makes the notebook
# reproducible without hiding the local equations inside a package object.
graph_aic <- function(...) sum(vapply(list(...), AIC, numeric(1L)))
graph_parameters <- function(...) {
  sum(vapply(list(...), function(model) attr(logLik(model), "df"), numeric(1L)))
}
primary_aic <- graph_aic(mediator_fit, outcome_fit)
constrained_aic <- graph_aic(mediator_fit, constrained_outcome_fit)
reverse_aic <- graph_aic(total_fit, reverse_mediator_fit)

# The constrained graph omits Antibiotic -> calprotectin. Its single basis-set
# claim is tested by the conventional (not HC3) coefficient test used by
# directed-separation/Fisher-C calculations.
dsep_test <- summary(outcome_fit)$coefficients["Antibiotic", , drop = FALSE]
dsep_p <- unname(dsep_test[1L, "Pr(>|t|)"])
constrained_fisher_c <- -2 * log(dsep_p)
constrained_fisher_df <- 2L
constrained_fisher_p <- pchisq(
  constrained_fisher_c,
  df = constrained_fisher_df,
  lower.tail = FALSE
)

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
  AIC = c(primary_aic, constrained_aic, reverse_aic),
  Parameters = c(
    graph_parameters(mediator_fit, outcome_fit),
    graph_parameters(mediator_fit, constrained_outcome_fit),
    graph_parameters(total_fit, reverse_mediator_fit)
  ),
  N = nrow(primary),
  IndependenceClaims = c(0L, 1L, 0L),
  FisherC = c(NA_real_, constrained_fisher_c, NA_real_),
  FisherDF = c(0L, constrained_fisher_df, 0L),
  FisherP = c(NA_real_, constrained_fisher_p, NA_real_),
  Saturated = c(TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
model_fit$DeltaAIC <- model_fit$AIC - min(model_fit$AIC)
write_tsv(model_fit, "sem-fit-comparison.tsv")

dsep_table <- data.frame(
  Model = "Constrained microbiome-only path",
  IndependenceClaim = "LogCalprotectinZ ~ Antibiotic + ...",
  TestType = "coef",
  DF = df.residual(outcome_fit),
  Estimate = unname(dsep_test[1L, "t value"]),
  PValue = dsep_p
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

# ---------------------------------------------------------------------------
# Design simulation: how much information is needed for the observed b path?
# ---------------------------------------------------------------------------

draw_exposure_with_overlap <- function(diagnosis, probability) {
  exposure <- integer(length(diagnosis))
  for (level in unique(diagnosis)) {
    index <- which(diagnosis == level)
    values <- rbinom(length(index), size = 1L, prob = probability)
    if (length(index) >= 2L && sum(values) == 0L) {
      values[sample.int(length(index), 1L)] <- 1L
    }
    if (length(index) >= 2L && sum(values) == length(index)) {
      values[sample.int(length(index), 1L)] <- 0L
    }
    exposure[index] <- values
  }
  exposure
}

fast_hc3_p <- function(model, term) {
  x <- model.matrix(model)
  y <- model.response(model.frame(model))
  fit <- .lm.fit(x, y)
  if (fit$rank != ncol(x)) return(NA_real_)
  inverse <- chol2inv(chol(crossprod(x)))
  residual <- as.vector(y - x %*% fit$coefficients)
  leverage <- rowSums((x %*% inverse) * x)
  omega <- (residual / pmax(1 - leverage, 1e-10))^2
  meat <- crossprod(x, x * omega)
  covariance <- inverse %*% meat %*% inverse
  index <- match(term, colnames(x))
  statistic <- fit$coefficients[index] / sqrt(covariance[index, index])
  2 * pt(-abs(statistic), df = nrow(x) - ncol(x))
}

simulate_once <- function(sample_size, exposure_probability) {
  index <- sample.int(nrow(primary), sample_size, replace = TRUE)
  simulated <- primary[index, covariates, drop = FALSE]
  diagnosis <- ifelse(simulated$CD == 1L, "CD", ifelse(simulated$UC == 1L, "UC", "Control"))
  simulated$Antibiotic <- draw_exposure_with_overlap(
    diagnosis,
    exposure_probability
  )
  mediator_x <- model.matrix(
    as.formula(paste("~ Antibiotic +", paste(covariates, collapse = " + "))),
    data = simulated
  )
  simulated$ShannonZ <- as.vector(mediator_x %*% coef(mediator_fit)) +
    rnorm(sample_size, sd = sigma(mediator_fit))
  outcome_x <- model.matrix(
    as.formula(
      paste(
        "~ Antibiotic + ShannonZ +",
        paste(covariates, collapse = " + ")
      )
    ),
    data = simulated
  )
  simulated$LogCalprotectinZ <- as.vector(outcome_x %*% coef(outcome_fit)) +
    rnorm(sample_size, sd = sigma(outcome_fit))
  mediator <- lm(formula_mediator, data = simulated)
  outcome <- lm(formula_outcome, data = simulated)
  c(
    A = fast_hc3_p(mediator, "Antibiotic"),
    B = fast_hc3_p(outcome, "ShannonZ")
  )
}

binomial_interval <- function(successes, repetitions) {
  unname(binom.test(successes, repetitions)$conf.int)
}

set.seed(SEED + 2000L)
power_sizes <- c(90L, 150L, 250L, 350L, 450L, 500L, 550L, 700L)
power_rows <- lapply(power_sizes, function(sample_size) {
  p_values <- replicate(
    POWER_REPETITIONS,
    simulate_once(sample_size, mean(primary$Antibiotic))
  )
  detected_b <- is.finite(p_values["B", ]) & p_values["B", ] < 0.05
  detected_joint <- detected_b & is.finite(p_values["A", ]) & p_values["A", ] < 0.05
  do.call(
    rbind,
    lapply(
      list("B path (HC3)" = detected_b, "Joint a and b (HC3)" = detected_joint),
      function(detected) {
        successes <- sum(detected)
        interval <- binomial_interval(successes, POWER_REPETITIONS)
        data.frame(
          SampleSize = sample_size,
          Target = deparse(substitute(detected)),
          ExposureFraction = mean(primary$Antibiotic),
          Repetitions = POWER_REPETITIONS,
          Detections = successes,
          Power = successes / POWER_REPETITIONS,
          CILower = interval[1L],
          CIUpper = interval[2L]
        )
      }
    )
  )
})
power_table <- do.call(rbind, power_rows)
# Replace deparse-generated labels with the stable list names.
power_table$Target <- rep(
  c("B path (HC3)", "Joint a and b (HC3)"),
  times = length(power_sizes)
)
rownames(power_table) <- NULL
write_tsv(power_table, "path-power-simulation.tsv")

# ---------------------------------------------------------------------------
# Positive control: known temporal DAG, overlap in every diagnosis stratum.
# ---------------------------------------------------------------------------

set.seed(SEED + 3000L)
positive_n <- 500L
positive_index <- sample.int(nrow(primary), positive_n, replace = TRUE)
positive <- primary[positive_index, covariates, drop = FALSE]
positive$Diagnosis <- ifelse(
  positive$CD == 1L,
  "CD",
  ifelse(positive$UC == 1L, "UC", "Control")
)
positive$Antibiotic <- draw_exposure_with_overlap(positive$Diagnosis, 0.30)
positive_m_coefficients <- coef(mediator_fit)
positive_m_coefficients[["Antibiotic"]] <- -0.60
positive_m_x <- model.matrix(
  as.formula(paste("~ Antibiotic +", paste(covariates, collapse = " + "))),
  data = positive
)
positive$ShannonZ <- as.vector(positive_m_x %*% positive_m_coefficients) +
  rnorm(positive_n, sd = sigma(mediator_fit))
positive_y_coefficients <- coef(outcome_fit)
positive_y_coefficients[["Antibiotic"]] <- -0.15
positive_y_coefficients[["ShannonZ"]] <- -0.30
positive_y_x <- model.matrix(
  as.formula(
    paste(
      "~ Antibiotic + ShannonZ +",
      paste(covariates, collapse = " + ")
    )
  ),
  data = positive
)
positive$LogCalprotectinZ <- as.vector(positive_y_x %*% positive_y_coefficients) +
  rnorm(positive_n, sd = sigma(outcome_fit))

positive_m_fit <- lm(formula_mediator, data = positive)
positive_y_fit <- lm(formula_outcome, data = positive)
positive_total_fit <- lm(formula_total, data = positive)
positive_bootstrap <- bootstrap_paths(
  positive_m_fit,
  positive_y_fit,
  positive_total_fit,
  POSITIVE_CONTROL_BOOTSTRAP,
  SEED + 4000L
)
positive_points <- c(
  A = coef(positive_m_fit)[["Antibiotic"]],
  B = coef(positive_y_fit)[["ShannonZ"]],
  Direct = coef(positive_y_fit)[["Antibiotic"]]
)
positive_points[["Indirect"]] <- positive_points[["A"]] * positive_points[["B"]]
positive_points[["Total"]] <- coef(positive_total_fit)[["Antibiotic"]]
positive_effects <- effect_summary(
  positive_bootstrap,
  positive_points,
  "Simulated positive control",
  POSITIVE_CONTROL_BOOTSTRAP
)
positive_truth <- c(
  A = -0.60,
  B = -0.30,
  Direct = -0.15,
  Indirect = 0.18,
  Total = 0.03
)
positive_effects$TrueValue <- positive_truth[positive_effects$Effect]
write_tsv(positive_effects, "positive-control-paths.tsv")

positive_overlap <- as.data.frame.matrix(
  table(positive$Diagnosis, positive$Antibiotic)
)
names(positive_overlap) <- c("Unexposed", "Exposed")
positive_overlap$Diagnosis <- rownames(positive_overlap)
rownames(positive_overlap) <- NULL
positive_overlap <- positive_overlap[
  match(c("Control", "CD", "UC"), positive_overlap$Diagnosis),
  c("Diagnosis", "Unexposed", "Exposed")
]
write_tsv(positive_overlap, "positive-control-overlap.tsv")

positive_reverse_fit <- lm(formula_reverse_m, data = positive)
positive_forward_aic <- graph_aic(positive_m_fit, positive_y_fit)
positive_reverse_aic <- graph_aic(positive_total_fit, positive_reverse_fit)
positive_audit <- data.frame(
  Criterion = c(
    "Temporal order",
    "Exposure overlap",
    "A path interval",
    "B path interval",
    "Indirect interval",
    "Forward versus reverse AIC"
  ),
  Result = c(
    "Known from generator",
    "Present in every diagnosis stratum",
    ifelse(
      positive_effects$CIUpper[positive_effects$Effect == "A"] < 0,
      "Excludes zero",
      "Crosses zero"
    ),
    ifelse(
      positive_effects$CIUpper[positive_effects$Effect == "B"] < 0,
      "Excludes zero",
      "Crosses zero"
    ),
    ifelse(
      positive_effects$CILower[positive_effects$Effect == "Indirect"] > 0,
      "Excludes zero",
      "Crosses zero"
    ),
    "Identical"
  ),
  Evidence = c(
    "A generated before M; M generated before Y",
    paste0(
      "Control ",
      positive_overlap$Exposed[positive_overlap$Diagnosis == "Control"],
      " exposed"
    ),
    "2,000 subject-level bootstrap refits",
    "2,000 subject-level bootstrap refits",
    "2,000 subject-level bootstrap refits",
    sprintf("%.6f versus %.6f", positive_forward_aic, positive_reverse_aic)
  )
)
write_tsv(positive_audit, "positive-control-audit.tsv")

# ---------------------------------------------------------------------------
# Residual-correlation sensitivity under a hypothetical causal interpretation.
# ---------------------------------------------------------------------------

set.seed(SEED + 5000L)
mediation_fit <- mediate(
  mediator_fit,
  outcome_fit,
  treat = "Antibiotic",
  mediator = "ShannonZ",
  control.value = 0,
  treat.value = 1,
  sims = 2000,
  boot = FALSE
)
mediation_sensitivity <- medsens(
  mediation_fit,
  rho.by = 0.05,
  effect.type = "indirect"
)
sensitivity_curve <- data.frame(
  Rho = mediation_sensitivity$rho,
  Indirect = mediation_sensitivity$d0,
  CILower = mediation_sensitivity$lower.d0,
  CIUpper = mediation_sensitivity$upper.d0
)
write_tsv(sensitivity_curve, "mediation-rho-sensitivity.tsv")
rho_zero_row <- sensitivity_curve[which.min(abs(sensitivity_curve$Rho)), ]
sensitivity_summary <- data.frame(
  Quantity = c(
    "Residual rho where point indirect effect is zero",
    "R2-star product threshold",
    "R2-tilde product threshold",
    "Indirect interval at rho=0 crosses zero"
  ),
  Value = c(
    mediation_sensitivity$err.cr.d,
    mediation_sensitivity$R2star.d.thresh,
    mediation_sensitivity$R2tilde.d.thresh,
    with(rho_zero_row, CILower <= 0 & CIUpper >= 0)
  )
)
write_tsv(sensitivity_summary, "mediation-rho-summary.tsv")

software <- data.frame(
  Package = c(
    "R", "mediation", "sandwich", "lmtest", "car", "jsonlite"
  ),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("mediation")),
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
    faec_outcome_fit = faec_outcome_fit,
    positive_mediator_fit = positive_m_fit,
    positive_outcome_fit = positive_y_fit,
    positive_total_fit = positive_total_fit,
    mediation_fit = mediation_fit,
    mediation_sensitivity = mediation_sensitivity
  ),
  file.path(output_dir, "sem-model-objects.rds"),
  compress = "xz"
)

power_candidates <- power_table$SampleSize[
  power_table$Target == "B path (HC3)" & power_table$Power >= 0.80
]
first_tested_n_with_b_power_80 <- if (length(power_candidates)) {
  min(power_candidates)
} else {
  NA_integer_
}

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
  constrained_fisher_c = constrained_fisher_c,
  constrained_fisher_df = constrained_fisher_df,
  constrained_fisher_p = constrained_fisher_p,
  primary_aic = primary_aic,
  constrained_aic = constrained_aic,
  reverse_aic = reverse_aic,
  propensity_converged = propensity_fit$converged,
  propensity_max_abs_coefficient = max(abs(coef(propensity_fit))),
  validation_antibiotic_exposed = sum(validation$Antibiotic),
  leave_one_out_indirect_min = min(leave_one_out$Indirect),
  leave_one_out_indirect_max = max(leave_one_out$Indirect),
  power_repetitions = POWER_REPETITIONS,
  b_power_n90 = power_table$Power[
    power_table$SampleSize == 90L & power_table$Target == "B path (HC3)"
  ],
  first_tested_n_with_b_power_80 = first_tested_n_with_b_power_80,
  positive_control_n = positive_n,
  positive_control_indirect = positive_points[["Indirect"]],
  positive_control_indirect_ci_lower = positive_effects$CILower[
    positive_effects$Effect == "Indirect"
  ],
  positive_control_indirect_ci_upper = positive_effects$CIUpper[
    positive_effects$Effect == "Indirect"
  ],
  positive_forward_aic = positive_forward_aic,
  positive_reverse_aic = positive_reverse_aic,
  medsens_rho_zero = mediation_sensitivity$err.cr.d,
  medsens_r2star_product_zero = mediation_sensitivity$R2star.d.thresh,
  medsens_r2tilde_product_zero = mediation_sensitivity$R2tilde.d.thresh
)
write_json(
  model_metrics,
  file.path(output_dir, "model-metrics.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
print(toJSON(model_metrics, pretty = TRUE, auto_unbox = TRUE, digits = 6))
