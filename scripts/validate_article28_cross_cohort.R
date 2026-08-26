#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999, digits = 17)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260728L
set.seed(primary_seed)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
    }
    out[[substring(args[[i]], 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("project-root", "input-dir", "output-dir", "figure-dir", "chapter")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c(
  "glmnet", "metafor", "pROC", "ggplot2", "patchwork",
  "scales", "jsonlite", "digest"
)
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
input_dir <- normalizePath(args[["input-dir"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
chapter_path <- normalizePath(args[["chapter"]], mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

font_cache_dir <- file.path(tempdir(), "article28-fontconfig-cache")
dir.create(font_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

validation_log <- file.path(output_dir, "validation.log")
writeLines(
  c(
    "Article 28 cross-cohort validation",
    paste0("StartedUTC\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Seed\t", primary_seed)
  ),
  validation_log
)
log_msg <- function(...) {
  line <- paste0(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    "\t", paste0(..., collapse = "")
  )
  cat(line, "\n")
  cat(line, "\n", file = validation_log, append = TRUE)
  invisible(line)
}

checks <- data.frame(
  Category = character(), CheckID = character(), Status = character(),
  Detail = character(), stringsAsFactors = FALSE
)
add_check <- function(category, check_id, passed, detail) {
  checks <<- rbind(
    checks,
    data.frame(
      Category = category,
      CheckID = check_id,
      Status = if (isTRUE(passed)) "PASS" else "FAIL",
      Detail = paste(detail, collapse = "; "),
      stringsAsFactors = FALSE
    )
  )
  invisible(passed)
}
sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x, con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
read_tsv <- function(name) {
  utils::read.delim(
    file.path(input_dir, name), check.names = FALSE,
    quote = "", comment.char = "", stringsAsFactors = FALSE
  )
}
read_tsv_gz <- function(name) {
  utils::read.delim(
    gzfile(file.path(input_dir, name)), check.names = FALSE,
    quote = "", comment.char = "", stringsAsFactors = FALSE
  )
}

verify_checksum_manifest <- function(directory) {
  manifest_path <- file.path(directory, "file-checksums.sha256")
  lines <- readLines(manifest_path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  expected_files <- character()
  for (line in lines) {
    pieces <- strsplit(line, "[[:space:]]+", perl = TRUE)[[1L]]
    expected <- pieces[[1L]]
    relative <- paste(pieces[-1L], collapse = " ")
    path <- file.path(directory, relative)
    observed <- if (file.exists(path)) sha256_file(path) else "missing"
    add_check(
      "Frozen input", paste0("sha256-", relative),
      identical(observed, expected), observed
    )
    expected_files <- c(expected_files, relative)
  }
  payloads <- sort(setdiff(
    basename(list.files(directory, full.names = TRUE)),
    "file-checksums.sha256"
  ))
  add_check(
    "Frozen input", "checksum-manifest-complete",
    identical(payloads, sort(expected_files)),
    paste0("payloads=", length(payloads), "; entries=", length(expected_files))
  )
  length(lines)
}

auc_score <- function(y, score) {
  y <- as.integer(y)
  keep <- is.finite(score) & y %in% c(0L, 1L)
  y <- y[keep]
  score <- score[keep]
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

average_precision <- function(y, score) {
  y <- as.integer(y)
  keep <- is.finite(score) & y %in% c(0L, 1L)
  y <- y[keep]
  score <- score[keep]
  ord <- order(-score, seq_along(score))
  y <- y[ord]
  score <- score[ord]
  threshold_ends <- which(!duplicated(score, fromLast = TRUE))
  tp <- cumsum(y == 1L)[threshold_ends]
  fp <- cumsum(y == 0L)[threshold_ends]
  recall <- tp / sum(y == 1L)
  precision <- tp / (tp + fp)
  sum(diff(c(0, recall)) * precision)
}

roc_curve <- function(y, score, cohort) {
  y <- as.integer(y)
  ord <- order(-score, seq_along(score))
  y <- y[ord]
  score <- score[ord]
  threshold_ends <- which(!duplicated(score, fromLast = TRUE))
  tp <- cumsum(y == 1L)[threshold_ends]
  fp <- cumsum(y == 0L)[threshold_ends]
  data.frame(
    Cohort = cohort,
    Threshold = c(Inf, score[threshold_ends], -Inf),
    FPR = c(0, fp / sum(y == 0L), 1),
    TPR = c(0, tp / sum(y == 1L), 1),
    stringsAsFactors = FALSE
  )
}

make_stratified_folds <- function(y, k, seed) {
  y <- factor(y)
  fold <- integer(length(y))
  set.seed(seed)
  for (level in levels(y)) {
    idx <- which(y == level)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

fit_preprocessor <- function(x, cohort) {
  cohort <- as.character(cohort)
  prevalence <- colMeans(x > 0)
  maximum <- apply(x, 2L, max)
  cohort_levels <- sort(unique(cohort))
  cohort_prevalence <- vapply(
    cohort_levels,
    function(z) colMeans(x[cohort == z, , drop = FALSE] > 0),
    numeric(ncol(x))
  )
  if (is.null(dim(cohort_prevalence))) {
    cohort_prevalence <- matrix(
      cohort_prevalence, ncol = 1L,
      dimnames = list(colnames(x), cohort_levels)
    )
  }
  support_required <- ceiling(length(cohort_levels) / 2)
  support <- rowSums(cohort_prevalence >= 0.05)
  variability <- apply(x, 2L, stats::sd) > 0
  keep <- prevalence >= 0.10 & maximum >= 1e-4 &
    support >= support_required & variability
  if (sum(keep) < 2L) {
    stop("Training data retained fewer than two features.", call. = FALSE)
  }
  list(
    features = colnames(x)[keep],
    pseudocount = 1e-6,
    prevalence = prevalence[keep],
    support = support[keep],
    cohorts = cohort_levels,
    support_required = support_required
  )
}

apply_preprocessor <- function(preprocessor, x) {
  out <- log10(
    x[, preprocessor$features, drop = FALSE] + preprocessor$pseudocount
  )
  storage.mode(out) <- "double"
  out
}

cohort_class_weights <- function(cohort, y) {
  cohort <- as.character(cohort)
  y <- as.character(y)
  cohort_levels <- sort(unique(cohort))
  weights <- numeric(length(y))
  for (z in cohort_levels) {
    for (class in c("Control", "CRC")) {
      idx <- which(cohort == z & y == class)
      if (length(idx) == 0L) {
        stop("Every training cohort must contain both outcome classes.", call. = FALSE)
      }
      weights[idx] <- 1 / (2 * length(cohort_levels) * length(idx))
    }
  }
  stopifnot(abs(sum(weights) - 1) < 1e-12)
  weights
}

lambda_grid <- exp(seq(log(1), log(1e-4), length.out = 41L))

fit_glmnet_path <- function(x, y, cohort, lambda = lambda_grid) {
  glmnet::glmnet(
    x = x,
    y = as.integer(y == "CRC"),
    family = "binomial",
    alpha = 1,
    lambda = lambda,
    weights = cohort_class_weights(cohort, y),
    standardize = TRUE,
    intercept = TRUE,
    thresh = 1e-8,
    maxit = 100000L
  )
}

predict_path <- function(model, x, lambda = lambda_grid) {
  out <- as.matrix(stats::predict(
    model, newx = x, type = "response", s = lambda
  ))
  if (ncol(out) != length(lambda)) {
    stop("Prediction path has an unexpected lambda dimension.", call. = FALSE)
  }
  out
}

select_lambda_one_se <- function(fold_auc, lambda = lambda_grid) {
  stopifnot(ncol(fold_auc) == length(lambda))
  mean_auc <- colMeans(fold_auc, na.rm = TRUE)
  valid_n <- colSums(is.finite(fold_auc))
  se_auc <- apply(fold_auc, 2L, stats::sd, na.rm = TRUE) / sqrt(valid_n)
  best <- which.max(mean_auc)
  threshold <- mean_auc[[best]] - se_auc[[best]]
  eligible <- is.finite(mean_auc) & mean_auc >= threshold
  selected <- which(eligible)[which.max(lambda[eligible])]
  list(
    selected_index = selected,
    best_index = best,
    selected_lambda = lambda[[selected]],
    best_lambda = lambda[[best]],
    threshold = threshold,
    mean_auc = mean_auc,
    se_auc = se_auc,
    valid_n = valid_n
  )
}

hedges_g <- function(values, outcome) {
  cases <- values[outcome == "CRC"]
  controls <- values[outcome == "Control"]
  if (
    length(cases) < 3L || length(controls) < 3L ||
      !is.finite(stats::sd(cases)) || !is.finite(stats::sd(controls)) ||
      stats::sd(cases) == 0 && stats::sd(controls) == 0
  ) {
    return(c(Effect = NA_real_, Variance = NA_real_))
  }
  effect <- tryCatch(
    metafor::escalc(
      measure = "SMD",
      m1i = mean(cases), sd1i = stats::sd(cases), n1i = length(cases),
      m2i = mean(controls), sd2i = stats::sd(controls), n2i = length(controls)
    ),
    error = function(e) NULL
  )
  if (is.null(effect)) return(c(Effect = NA_real_, Variance = NA_real_))
  c(Effect = as.numeric(effect$yi), Variance = as.numeric(effect$vi))
}

checksum_entries <- verify_checksum_manifest(input_dir)
notice_path <- file.path(project_root, "data", "small", "28-data-NOTICE.txt")
notice <- readLines(notice_path, warn = FALSE)
add_check("Frozen input", "notice-present", length(notice) >= 11L, notice_path)
for (token in c(
  "771 independent subjects", "385 Control and 386 CRC",
  "Every outer iteration trains on seven complete cohorts",
  "Feature filtering, standardization and lambda selection",
  "exact numeric reproduction is not claimed"
)) {
  add_check(
    "Frozen input", paste0("notice-", gsub("[^a-z0-9]+", "-", tolower(token))),
    any(grepl(token, notice, fixed = TRUE)), token
  )
}

species_table <- read_tsv_gz("species-relative-abundance.tsv.gz")
metadata <- read_tsv("sample-metadata.tsv")
cohort_summary <- read_tsv("cohort-summary.tsv")
feature_audit <- read_tsv_gz("feature-universe-audit.tsv.gz")
contract <- read_tsv("analysis-contract.tsv")
resource_manifest <- read_tsv("resource-manifest.tsv")

sample_ids <- metadata$sample_id
x <- t(as.matrix(species_table[, sample_ids, drop = FALSE]))
storage.mode(x) <- "double"
colnames(x) <- species_table$Feature
rownames(x) <- sample_ids
y <- factor(metadata$Outcome, levels = c("Control", "CRC"))
cohort_levels <- as.character(cohort_summary$Cohort)
cohort <- factor(metadata$Cohort, levels = cohort_levels)
role <- stats::setNames(cohort_summary$Role, cohort_summary$Cohort)
species_labels <- stats::setNames(species_table$Species, species_table$Feature)

add_check("Data", "samples", nrow(x) == 771L, nrow(x))
add_check("Data", "features", ncol(x) == 897L, ncol(x))
add_check("Data", "cohorts", length(cohort_levels) == 8L, paste(cohort_levels, collapse = "/"))
add_check("Data", "sample-alignment", identical(rownames(x), metadata$sample_id), "matrix vs metadata")
add_check("Data", "independent-subjects", !anyDuplicated(metadata$subject_id), length(unique(metadata$subject_id)))
add_check("Data", "study-qualified-keys", !anyDuplicated(metadata$StudySampleKey), length(unique(metadata$StudySampleKey)))
add_check("Data", "stool-only", all(metadata$body_site == "stool"), paste(unique(metadata$body_site), collapse = "/"))
add_check("Data", "class-counts", identical(as.integer(table(y)), c(385L, 386L)), paste(table(y), collapse = "/"))
add_check("Data", "non-target-excluded", all(metadata$study_condition %in% c("control", "CRC")), paste(unique(metadata$study_condition), collapse = "/"))
add_check("Data", "finite", all(is.finite(x)), sum(!is.finite(x)))
add_check("Data", "nonnegative", min(x) >= 0, min(x))
add_check("Data", "fraction-unit", max(x) <= 1, max(x))
add_check("Data", "closure", max(abs(rowSums(x) - 1)) < 2e-6, range(rowSums(x)))
add_check("Data", "feature-audit-rows", nrow(feature_audit) == 897L, nrow(feature_audit))
add_check("Data", "resource-manifest-rows", nrow(resource_manifest) == 10L, nrow(resource_manifest))
add_check("Data", "cohort-summary-rows", nrow(cohort_summary) == 8L, nrow(cohort_summary))
add_check("Data", "cohort-sample-counts", identical(as.integer(table(cohort)), cohort_summary$Samples), paste(table(cohort), collapse = "/"))
add_check("Data", "cohort-control-counts", identical(as.integer(table(cohort[y == "Control"])), cohort_summary$Controls), paste(table(cohort[y == "Control"]), collapse = "/"))
add_check("Data", "cohort-crc-counts", identical(as.integer(table(cohort[y == "CRC"])), cohort_summary$CRC), paste(table(cohort[y == "CRC"]), collapse = "/"))

contract_value <- function(item) contract$Value[match(item, contract$Item)]
expected_contract <- c(
  seed = "20260728", resource_release = "2021-03-31", cohorts = "8",
  samples = "771", control_samples = "385", crc_samples = "386",
  raw_species_features = "897", positive_class = "CRC",
  outer_test_unit = "study_name", inner_validation_unit = "study_name",
  log10_pseudocount = "0.000001"
)
for (item in names(expected_contract)) {
  observed <- contract_value(item)
  add_check(
    "Contract", paste0("contract-", gsub("_", "-", item)),
    identical(observed, expected_contract[[item]]), observed
  )
}
add_check("Contract", "true-n-minus-one", grepl("N-1 cohorts train", contract_value("outer_validation"), fixed = TRUE), contract_value("outer_validation"))
add_check("Contract", "inner-cohort-holdout", grepl("leave-one-training-cohort-out", contract_value("inner_validation"), fixed = TRUE), contract_value("inner_validation"))
add_check("Contract", "no-joint-batch-correction", grepl("no joint train-test ComBat", contract_value("batch_correction_policy"), fixed = TRUE), contract_value("batch_correction_policy"))

log_msg("starting discovery-only species meta-analysis")
discovery_cohorts <- cohort_levels[grepl("Discovery", role[cohort_levels], fixed = TRUE)]
validation_cohorts <- setdiff(cohort_levels, discovery_cohorts)
discovery_idx <- cohort %in% discovery_cohorts
discovery_support <- vapply(
  discovery_cohorts,
  function(z) colMeans(x[cohort == z, , drop = FALSE] > 0),
  numeric(ncol(x))
)
meta_features <- colnames(x)[
  colMeans(x[discovery_idx, , drop = FALSE] > 0) >= 0.10 &
    rowSums(discovery_support >= 0.05) >= 3L &
    apply(x[discovery_idx, , drop = FALSE], 2L, max) >= 1e-4
]

association_effects <- vector("list", length(meta_features) * length(cohort_levels))
effect_index <- 0L
for (feature in meta_features) {
  transformed <- log10(x[, feature] + 1e-6)
  for (z in cohort_levels) {
    idx <- cohort == z
    effect <- hedges_g(transformed[idx], y[idx])
    effect_index <- effect_index + 1L
    association_effects[[effect_index]] <- data.frame(
      Feature = feature,
      Species = unname(species_labels[[feature]]),
      Cohort = z,
      Role = unname(role[[z]]),
      Controls = sum(y[idx] == "Control"),
      CRC = sum(y[idx] == "CRC"),
      HedgesG = effect[["Effect"]],
      Variance = effect[["Variance"]],
      StandardError = sqrt(effect[["Variance"]]),
      stringsAsFactors = FALSE
    )
  }
}
association_effects <- do.call(rbind, association_effects)

association_meta_rows <- vector("list", length(meta_features))
for (i in seq_along(meta_features)) {
  feature <- meta_features[[i]]
  z <- association_effects[
    association_effects$Feature == feature &
      association_effects$Cohort %in% discovery_cohorts &
      is.finite(association_effects$HedgesG) &
      is.finite(association_effects$Variance) &
      association_effects$Variance > 0,
    , drop = FALSE
  ]
  if (nrow(z) >= 3L) {
    fit <- metafor::rma.uni(
      yi = z$HedgesG, vi = z$Variance, method = "REML"
    )
    estimate <- as.numeric(fit$b)
    direction_fraction <- mean(sign(z$HedgesG) == sign(estimate))
    association_meta_rows[[i]] <- data.frame(
      Feature = feature,
      Species = unname(species_labels[[feature]]),
      DiscoveryCohorts = nrow(z),
      PooledHedgesG = estimate,
      Lower95 = fit$ci.lb,
      Upper95 = fit$ci.ub,
      PValue = fit$pval,
      TauSquared = fit$tau2,
      I2Percent = fit$I2,
      HeterogeneityP = fit$QEp,
      DiscoveryDirectionFraction = direction_fraction,
      stringsAsFactors = FALSE
    )
  } else {
    association_meta_rows[[i]] <- data.frame(
      Feature = feature,
      Species = unname(species_labels[[feature]]),
      DiscoveryCohorts = nrow(z),
      PooledHedgesG = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_,
      PValue = NA_real_, TauSquared = NA_real_, I2Percent = NA_real_,
      HeterogeneityP = NA_real_, DiscoveryDirectionFraction = NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
association_meta <- do.call(rbind, association_meta_rows)
association_meta$QValue <- stats::p.adjust(association_meta$PValue, method = "BH")
association_meta <- association_meta[order(
  association_meta$QValue, -abs(association_meta$PooledHedgesG),
  association_meta$Species
), , drop = FALSE]

validation_concordance <- do.call(
  rbind,
  lapply(seq_len(nrow(association_meta)), function(i) {
    row <- association_meta[i, , drop = FALSE]
    z <- association_effects[
      association_effects$Feature == row$Feature &
        association_effects$Cohort %in% validation_cohorts,
      , drop = FALSE
    ]
    valid <- is.finite(z$HedgesG)
    data.frame(
      Feature = row$Feature,
      Species = row$Species,
      DiscoveryQValue = row$QValue,
      PooledHedgesG = row$PooledHedgesG,
      ValidationCohorts = sum(valid),
      ConcordantValidationCohorts = sum(
        sign(z$HedgesG[valid]) == sign(row$PooledHedgesG)
      ),
      ValidationDirectionFraction = if (sum(valid) > 0L) {
        mean(sign(z$HedgesG[valid]) == sign(row$PooledHedgesG))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
)
meta_signature <- merge(
  association_meta, validation_concordance,
  by = c("Feature", "Species", "PooledHedgesG"), sort = FALSE
)
meta_signature <- meta_signature[order(
  meta_signature$QValue,
  -meta_signature$ValidationDirectionFraction,
  -abs(meta_signature$PooledHedgesG)
), , drop = FALSE]
log_msg(
  "meta-analysis eligible species=", length(meta_features),
  "; q<0.05=", sum(meta_signature$QValue < 0.05, na.rm = TRUE)
)

log_msg("starting true N-1 outer validation")
lodo_predictions <- list()
lodo_performance <- list()
lodo_tuning <- list()
lodo_preprocessing <- list()
lodo_coefficients <- list()
prediction_index <- performance_index <- tuning_index <- preprocessing_index <- coefficient_index <- 0L

for (outer_index in seq_along(cohort_levels)) {
  test_cohort <- cohort_levels[[outer_index]]
  test <- cohort == test_cohort
  train <- !test
  training_cohorts <- setdiff(cohort_levels, test_cohort)
  inner_auc <- matrix(
    NA_real_, nrow = length(training_cohorts), ncol = length(lambda_grid),
    dimnames = list(training_cohorts, sprintf("%.10g", lambda_grid))
  )

  for (inner_index in seq_along(training_cohorts)) {
    validation_cohort <- training_cohorts[[inner_index]]
    inner_validation <- train & cohort == validation_cohort
    inner_train <- train & cohort != validation_cohort
    preprocessor <- fit_preprocessor(x[inner_train, , drop = FALSE], cohort[inner_train])
    x_inner_train <- apply_preprocessor(preprocessor, x[inner_train, , drop = FALSE])
    x_inner_validation <- apply_preprocessor(preprocessor, x[inner_validation, , drop = FALSE])
    model <- fit_glmnet_path(
      x_inner_train, y[inner_train], cohort[inner_train], lambda_grid
    )
    probability <- predict_path(model, x_inner_validation, lambda_grid)
    inner_auc[inner_index, ] <- vapply(
      seq_along(lambda_grid),
      function(j) auc_score(y[inner_validation] == "CRC", probability[, j]),
      numeric(1L)
    )
    for (j in seq_along(lambda_grid)) {
      tuning_index <- tuning_index + 1L
      lodo_tuning[[tuning_index]] <- data.frame(
        OuterTestCohort = test_cohort,
        InnerValidationCohort = validation_cohort,
        Lambda = lambda_grid[[j]],
        AUROC = inner_auc[inner_index, j],
        InnerTrainingCohorts = paste(
          setdiff(training_cohorts, validation_cohort), collapse = ";"
        ),
        FeaturesKept = length(preprocessor$features),
        stringsAsFactors = FALSE
      )
    }
  }

  selected <- select_lambda_one_se(inner_auc, lambda_grid)
  preprocessor <- fit_preprocessor(x[train, , drop = FALSE], cohort[train])
  x_train <- apply_preprocessor(preprocessor, x[train, , drop = FALSE])
  x_test <- apply_preprocessor(preprocessor, x[test, , drop = FALSE])
  model <- fit_glmnet_path(
    x_train, y[train], cohort[train], selected$selected_lambda
  )
  probability <- as.numeric(stats::predict(
    model, newx = x_test, type = "response", s = selected$selected_lambda
  ))

  truth <- factor(y[test], levels = c("Control", "CRC"))
  binary <- as.integer(truth == "CRC")
  roc_object <- pROC::roc(
    response = truth, predictor = probability,
    levels = c("Control", "CRC"), direction = "<", quiet = TRUE
  )
  interval <- as.numeric(pROC::ci.auc(
    roc_object, conf.level = 0.95, method = "delong"
  ))
  auc_value <- auc_score(binary, probability)
  auc_variance <- as.numeric(pROC::var(roc_object, method = "delong"))
  predicted <- probability >= 0.5

  prediction_index <- prediction_index + 1L
  lodo_predictions[[prediction_index]] <- data.frame(
    SampleID = sample_ids[test],
    SubjectID = metadata$subject_id[test],
    Study = metadata$study_name[test],
    Cohort = test_cohort,
    Outcome = as.character(truth),
    ProbabilityCRC = probability,
    OuterTrainingCohorts = paste(training_cohorts, collapse = ";"),
    stringsAsFactors = FALSE
  )
  performance_index <- performance_index + 1L
  lodo_performance[[performance_index]] <- data.frame(
    Cohort = test_cohort,
    Study = unique(metadata$study_name[test]),
    Role = unname(role[[test_cohort]]),
    Samples = sum(test),
    Controls = sum(truth == "Control"),
    CRC = sum(truth == "CRC"),
    TrainingCohorts = length(training_cohorts),
    TrainingSamples = sum(train),
    FeaturesKept = length(preprocessor$features),
    SelectedLambda = selected$selected_lambda,
    BestMeanInnerAUROC = selected$mean_auc[[selected$best_index]],
    SelectedMeanInnerAUROC = selected$mean_auc[[selected$selected_index]],
    AUROC = auc_value,
    AUROCLower95 = interval[[1L]],
    AUROCUpper95 = interval[[3L]],
    AUROCVariance = auc_variance,
    AUPRC = average_precision(binary, probability),
    Brier = mean((binary - probability)^2),
    SensitivityAt05 = sum(predicted & binary == 1L) / sum(binary == 1L),
    SpecificityAt05 = sum(!predicted & binary == 0L) / sum(binary == 0L),
    stringsAsFactors = FALSE
  )
  preprocessing_index <- preprocessing_index + 1L
  lodo_preprocessing[[preprocessing_index]] <- data.frame(
    OuterTestCohort = test_cohort,
    OuterTestSamples = sum(test),
    TrainingCohorts = paste(training_cohorts, collapse = ";"),
    TrainingSamples = sum(train),
    FeaturesAvailable = ncol(x),
    FeaturesKept = length(preprocessor$features),
    CohortSupportRequired = preprocessor$support_required,
    Pseudocount = preprocessor$pseudocount,
    SelectedLambda = selected$selected_lambda,
    FittedOnOuterTrainingOnly = TRUE,
    stringsAsFactors = FALSE
  )

  coefficient_matrix <- as.matrix(stats::coef(
    model, s = selected$selected_lambda
  ))
  coefficient_matrix <- coefficient_matrix[setdiff(
    rownames(coefficient_matrix), "(Intercept)"
  ), , drop = FALSE]
  for (feature in rownames(coefficient_matrix)) {
    coefficient_index <- coefficient_index + 1L
    coefficient <- as.numeric(coefficient_matrix[feature, 1L])
    lodo_coefficients[[coefficient_index]] <- data.frame(
      OuterTestCohort = test_cohort,
      Feature = feature,
      Species = unname(species_labels[[feature]]),
      Coefficient = coefficient,
      Selected = coefficient != 0,
      stringsAsFactors = FALSE
    )
  }
  log_msg(
    "LODO ", outer_index, "/8: test=", test_cohort,
    "; AUROC=", sprintf("%.3f", auc_value),
    "; lambda=", sprintf("%.5g", selected$selected_lambda),
    "; features=", length(preprocessor$features)
  )
}

lodo_predictions <- do.call(rbind, lodo_predictions)
lodo_performance <- do.call(rbind, lodo_performance)
lodo_tuning <- do.call(rbind, lodo_tuning)
lodo_preprocessing <- do.call(rbind, lodo_preprocessing)
lodo_coefficients <- do.call(rbind, lodo_coefficients)
lodo_performance$Cohort <- factor(lodo_performance$Cohort, levels = cohort_levels)
lodo_performance <- lodo_performance[order(lodo_performance$Cohort), , drop = FALSE]
lodo_performance$Cohort <- as.character(lodo_performance$Cohort)

tuning_summary <- do.call(
  rbind,
  lapply(split(lodo_tuning, lodo_tuning$OuterTestCohort), function(z) {
    by_lambda <- split(z$AUROC, z$Lambda)
    means <- vapply(by_lambda, mean, numeric(1L))
    ses <- vapply(by_lambda, function(v) stats::sd(v) / sqrt(length(v)), numeric(1L))
    lambdas <- as.numeric(names(by_lambda))
    order_idx <- order(lambdas, decreasing = TRUE)
    data.frame(
      OuterTestCohort = z$OuterTestCohort[[1L]],
      Lambda = lambdas[order_idx],
      MeanInnerAUROC = means[order_idx],
      SEInnerAUROC = ses[order_idx],
      Selected = lambdas[order_idx] == lodo_performance$SelectedLambda[
        match(z$OuterTestCohort[[1L]], lodo_performance$Cohort)
      ],
      stringsAsFactors = FALSE
    )
  })
)

logit_auc <- stats::qlogis(lodo_performance$AUROC)
se_logit_auc <- sqrt(lodo_performance$AUROCVariance) /
  (lodo_performance$AUROC * (1 - lodo_performance$AUROC))
performance_meta_fit <- metafor::rma.uni(
  yi = logit_auc, sei = se_logit_auc, method = "REML"
)
performance_meta <- data.frame(
  Cohorts = nrow(lodo_performance),
  Model = "Random-effects REML on logit AUROC",
  PooledAUROC = stats::plogis(as.numeric(performance_meta_fit$b)),
  Lower95 = stats::plogis(performance_meta_fit$ci.lb),
  Upper95 = stats::plogis(performance_meta_fit$ci.ub),
  TauSquaredLogitScale = performance_meta_fit$tau2,
  I2Percent = performance_meta_fit$I2,
  HeterogeneityP = performance_meta_fit$QEp,
  stringsAsFactors = FALSE
)

set.seed(primary_seed + 800000L)
bootstrap_replicates <- 2000L
bootstrap_macro_auc <- numeric(bootstrap_replicates)
for (b in seq_len(bootstrap_replicates)) {
  sampled_cohorts <- sample(
    cohort_levels, length(cohort_levels), replace = TRUE
  )
  cohort_aucs <- numeric(length(sampled_cohorts))
  for (j in seq_along(sampled_cohorts)) {
    z <- lodo_predictions[lodo_predictions$Cohort == sampled_cohorts[[j]], , drop = FALSE]
    control_rows <- which(z$Outcome == "Control")
    crc_rows <- which(z$Outcome == "CRC")
    sampled_rows <- c(
      sample(control_rows, length(control_rows), replace = TRUE),
      sample(crc_rows, length(crc_rows), replace = TRUE)
    )
    cohort_aucs[[j]] <- auc_score(
      z$Outcome[sampled_rows] == "CRC", z$ProbabilityCRC[sampled_rows]
    )
  }
  bootstrap_macro_auc[[b]] <- mean(cohort_aucs)
}
bootstrap_summary <- data.frame(
  Estimand = "Unweighted macro mean of eight cohort AUROCs",
  Estimate = mean(lodo_performance$AUROC),
  Replicates = bootstrap_replicates,
  Lower95 = unname(stats::quantile(bootstrap_macro_auc, 0.025, type = 6)),
  Upper95 = unname(stats::quantile(bootstrap_macro_auc, 0.975, type = 6)),
  Median = stats::median(bootstrap_macro_auc),
  stringsAsFactors = FALSE
)
bootstrap_table <- data.frame(
  Replicate = seq_len(bootstrap_replicates),
  MacroAUROC = bootstrap_macro_auc,
  stringsAsFactors = FALSE
)

coefficient_stability <- do.call(
  rbind,
  lapply(split(lodo_coefficients, lodo_coefficients$Feature), function(z) {
    selected <- z$Coefficient != 0
    nonzero <- z$Coefficient[selected]
    data.frame(
      Feature = z$Feature[[1L]],
      Species = z$Species[[1L]],
      ModelsAvailable = nrow(z),
      ModelsSelected = sum(selected),
      SelectionFraction = sum(selected) / length(cohort_levels),
      PositiveModels = sum(nonzero > 0),
      NegativeModels = sum(nonzero < 0),
      SignConsistency = if (length(nonzero) > 0L) {
        max(sum(nonzero > 0), sum(nonzero < 0)) / length(nonzero)
      } else {
        NA_real_
      },
      MedianAbsoluteCoefficient = if (length(nonzero) > 0L) {
        stats::median(abs(nonzero))
      } else {
        0
      },
      stringsAsFactors = FALSE
    )
  })
)
coefficient_stability <- coefficient_stability[order(
  -coefficient_stability$ModelsSelected,
  -coefficient_stability$SignConsistency,
  -coefficient_stability$MedianAbsoluteCoefficient,
  coefficient_stability$Species
), , drop = FALSE]

roc_data <- do.call(
  rbind,
  lapply(cohort_levels, function(z) {
    rows <- lodo_predictions$Cohort == z
    roc_curve(
      lodo_predictions$Outcome[rows] == "CRC",
      lodo_predictions$ProbabilityCRC[rows], z
    )
  })
)

log_msg("starting secondary 8-by-8 single-study transfer")
transfer_rows <- list()
transfer_tuning_rows <- list()
transfer_index <- transfer_tuning_index <- 0L
for (train_index in seq_along(cohort_levels)) {
  train_cohort <- cohort_levels[[train_index]]
  train_rows <- which(cohort == train_cohort)
  folds <- make_stratified_folds(
    y[train_rows], 5L, primary_seed + 900000L + train_index * 1000L
  )
  fold_auc <- matrix(
    NA_real_, nrow = 5L, ncol = length(lambda_grid),
    dimnames = list(seq_len(5L), sprintf("%.10g", lambda_grid))
  )
  oof_path <- matrix(
    NA_real_, nrow = length(train_rows), ncol = length(lambda_grid)
  )
  for (fold_id in seq_len(5L)) {
    inner_train_rows <- train_rows[folds != fold_id]
    inner_test_rows <- train_rows[folds == fold_id]
    preprocessor <- fit_preprocessor(
      x[inner_train_rows, , drop = FALSE], cohort[inner_train_rows]
    )
    model <- fit_glmnet_path(
      apply_preprocessor(preprocessor, x[inner_train_rows, , drop = FALSE]),
      y[inner_train_rows], cohort[inner_train_rows], lambda_grid
    )
    probability <- predict_path(
      model,
      apply_preprocessor(preprocessor, x[inner_test_rows, , drop = FALSE]),
      lambda_grid
    )
    oof_path[folds == fold_id, ] <- probability
    fold_auc[fold_id, ] <- vapply(
      seq_along(lambda_grid),
      function(j) auc_score(y[inner_test_rows] == "CRC", probability[, j]),
      numeric(1L)
    )
    for (j in seq_along(lambda_grid)) {
      transfer_tuning_index <- transfer_tuning_index + 1L
      transfer_tuning_rows[[transfer_tuning_index]] <- data.frame(
        TrainingCohort = train_cohort,
        Fold = fold_id,
        Lambda = lambda_grid[[j]],
        AUROC = fold_auc[fold_id, j],
        FeaturesKept = length(preprocessor$features),
        stringsAsFactors = FALSE
      )
    }
  }
  selected <- select_lambda_one_se(fold_auc, lambda_grid)
  full_preprocessor <- fit_preprocessor(
    x[train_rows, , drop = FALSE], cohort[train_rows]
  )
  full_model <- fit_glmnet_path(
    apply_preprocessor(full_preprocessor, x[train_rows, , drop = FALSE]),
    y[train_rows], cohort[train_rows], selected$selected_lambda
  )
  for (test_cohort in cohort_levels) {
    test_rows <- which(cohort == test_cohort)
    if (test_cohort == train_cohort) {
      probability <- oof_path[, selected$selected_index]
      truth <- y[train_rows]
      validation_type <- "Within-cohort 5-fold CV"
    } else {
      probability <- as.numeric(stats::predict(
        full_model,
        newx = apply_preprocessor(
          full_preprocessor, x[test_rows, , drop = FALSE]
        ),
        type = "response", s = selected$selected_lambda
      ))
      truth <- y[test_rows]
      validation_type <- "External single-study transfer"
    }
    transfer_index <- transfer_index + 1L
    transfer_rows[[transfer_index]] <- data.frame(
      TrainingCohort = train_cohort,
      TestCohort = test_cohort,
      ValidationType = validation_type,
      TrainingSamples = length(train_rows),
      TestSamples = length(truth),
      FeaturesKept = length(full_preprocessor$features),
      SelectedLambda = selected$selected_lambda,
      AUROC = auc_score(truth == "CRC", probability),
      stringsAsFactors = FALSE
    )
  }
  log_msg("single-study transfer source ", train_index, "/8: ", train_cohort)
}
single_transfer <- do.call(rbind, transfer_rows)
single_transfer_tuning <- do.call(rbind, transfer_tuning_rows)

validation_gap <- do.call(
  rbind,
  lapply(cohort_levels, function(z) {
    within <- single_transfer$AUROC[
      single_transfer$TrainingCohort == z & single_transfer$TestCohort == z
    ]
    external <- single_transfer$AUROC[
      single_transfer$TrainingCohort != z & single_transfer$TestCohort == z
    ]
    lodo <- lodo_performance$AUROC[lodo_performance$Cohort == z]
    data.frame(
      Cohort = z,
      WithinCohortCV = within,
      MeanSingleStudyTransfer = mean(external),
      MinimumSingleStudyTransfer = min(external),
      MaximumSingleStudyTransfer = max(external),
      NMinusOneLODO = lodo,
      LODOGainOverMeanSingleStudy = lodo - mean(external),
      stringsAsFactors = FALSE
    )
  })
)

add_check("Association meta-analysis", "discovery-cohorts", length(discovery_cohorts) == 5L, paste(discovery_cohorts, collapse = "/"))
add_check("Association meta-analysis", "validation-cohorts", length(validation_cohorts) == 3L, paste(validation_cohorts, collapse = "/"))
add_check("Association meta-analysis", "eligible-features", length(meta_features) >= 50L, length(meta_features))
add_check("Association meta-analysis", "effect-rows", nrow(association_effects) == length(meta_features) * 8L, nrow(association_effects))
add_check("Association meta-analysis", "meta-rows", nrow(association_meta) == length(meta_features), nrow(association_meta))
add_check("Association meta-analysis", "finite-models", sum(is.finite(association_meta$PooledHedgesG)) >= 50L, sum(is.finite(association_meta$PooledHedgesG)))
add_check("Association meta-analysis", "bh-range", all(association_meta$QValue[is.finite(association_meta$QValue)] >= 0 & association_meta$QValue[is.finite(association_meta$QValue)] <= 1), range(association_meta$QValue, na.rm = TRUE))
add_check("Association meta-analysis", "signature-detected", sum(meta_signature$QValue < 0.05, na.rm = TRUE) >= 1L, sum(meta_signature$QValue < 0.05, na.rm = TRUE))

add_check("LODO", "outer-models", nrow(lodo_performance) == 8L, nrow(lodo_performance))
add_check("LODO", "one-test-per-cohort", identical(lodo_performance$Cohort, cohort_levels), paste(lodo_performance$Cohort, collapse = "/"))
add_check("LODO", "predictions", nrow(lodo_predictions) == 771L, nrow(lodo_predictions))
add_check("LODO", "one-prediction-per-sample", !anyDuplicated(lodo_predictions$SampleID), sum(duplicated(lodo_predictions$SampleID)))
add_check("LODO", "prediction-range", all(lodo_predictions$ProbabilityCRC >= 0 & lodo_predictions$ProbabilityCRC <= 1), range(lodo_predictions$ProbabilityCRC))
add_check("LODO", "test-counts", identical(lodo_performance$Samples, cohort_summary$Samples), paste(lodo_performance$Samples, collapse = "/"))
add_check("LODO", "seven-training-cohorts", all(lodo_performance$TrainingCohorts == 7L), paste(lodo_performance$TrainingCohorts, collapse = "/"))
add_check("LODO", "training-sample-complement", all(lodo_performance$TrainingSamples + lodo_performance$Samples == 771L), paste(lodo_performance$TrainingSamples, collapse = "/"))
add_check("LODO", "finite-auroc", all(is.finite(lodo_performance$AUROC)), paste(round(lodo_performance$AUROC, 3), collapse = "/"))
add_check("LODO", "auc-range", all(lodo_performance$AUROC >= 0 & lodo_performance$AUROC <= 1), range(lodo_performance$AUROC))
add_check("LODO", "ci-contains-estimate", all(lodo_performance$AUROCLower95 <= lodo_performance$AUROC & lodo_performance$AUROC <= lodo_performance$AUROCUpper95), "DeLong intervals")
add_check("LODO", "inner-audit-rows", nrow(lodo_tuning) == 8L * 7L * 41L, nrow(lodo_tuning))
add_check("LODO", "inner-summary-rows", nrow(tuning_summary) == 8L * 41L, nrow(tuning_summary))
add_check("LODO", "selected-one-per-outer", all(table(tuning_summary$OuterTestCohort[tuning_summary$Selected]) == 1L), paste(table(tuning_summary$OuterTestCohort[tuning_summary$Selected]), collapse = "/"))
add_check("LODO", "selected-lambda-on-grid", all(vapply(lodo_performance$SelectedLambda, function(z) any(abs(lambda_grid - z) < 1e-12), logical(1L))), paste(signif(lodo_performance$SelectedLambda, 4), collapse = "/"))
add_check("LODO", "features-retained", all(lodo_performance$FeaturesKept >= 10L & lodo_performance$FeaturesKept <= 897L), paste(lodo_performance$FeaturesKept, collapse = "/"))
add_check("LODO", "preprocessing-audit", nrow(lodo_preprocessing) == 8L && all(lodo_preprocessing$FittedOnOuterTrainingOnly), paste(lodo_preprocessing$FittedOnOuterTrainingOnly, collapse = "/"))
add_check("LODO", "outer-test-excluded", all(vapply(seq_len(nrow(lodo_preprocessing)), function(i) !lodo_preprocessing$OuterTestCohort[[i]] %in% strsplit(lodo_preprocessing$TrainingCohorts[[i]], ";", fixed = TRUE)[[1L]], logical(1L))), "test cohort absent from training list")
add_check("LODO", "coefficient-models", length(unique(lodo_coefficients$OuterTestCohort)) == 8L, length(unique(lodo_coefficients$OuterTestCohort)))
add_check("LODO", "nonzero-coefficients", sum(lodo_coefficients$Selected) >= 8L, sum(lodo_coefficients$Selected))

add_check("Performance synthesis", "meta-row", nrow(performance_meta) == 1L, nrow(performance_meta))
add_check("Performance synthesis", "pooled-auc-range", performance_meta$PooledAUROC > 0.5 && performance_meta$PooledAUROC < 1, performance_meta$PooledAUROC)
add_check("Performance synthesis", "meta-ci", performance_meta$Lower95 < performance_meta$PooledAUROC && performance_meta$PooledAUROC < performance_meta$Upper95, paste(performance_meta$Lower95, performance_meta$Upper95, sep = "/"))
add_check("Performance synthesis", "heterogeneity", performance_meta$I2Percent >= 0 && performance_meta$I2Percent <= 100, performance_meta$I2Percent)
add_check("Performance synthesis", "bootstrap-replicates", nrow(bootstrap_table) == 2000L, nrow(bootstrap_table))
add_check("Performance synthesis", "bootstrap-finite", all(is.finite(bootstrap_table$MacroAUROC)), range(bootstrap_table$MacroAUROC))
add_check("Performance synthesis", "bootstrap-ci", bootstrap_summary$Lower95 < bootstrap_summary$Estimate && bootstrap_summary$Estimate < bootstrap_summary$Upper95, paste(bootstrap_summary$Lower95, bootstrap_summary$Upper95, sep = "/"))

add_check("Single-study transfer", "matrix-cells", nrow(single_transfer) == 64L, nrow(single_transfer))
add_check("Single-study transfer", "diagonal", sum(single_transfer$TrainingCohort == single_transfer$TestCohort) == 8L, sum(single_transfer$TrainingCohort == single_transfer$TestCohort))
add_check("Single-study transfer", "external", sum(single_transfer$TrainingCohort != single_transfer$TestCohort) == 56L, sum(single_transfer$TrainingCohort != single_transfer$TestCohort))
add_check("Single-study transfer", "finite-auc", all(is.finite(single_transfer$AUROC)), range(single_transfer$AUROC))
add_check("Single-study transfer", "tuning-rows", nrow(single_transfer_tuning) == 8L * 5L * 41L, nrow(single_transfer_tuning))
add_check("Single-study transfer", "validation-gap-rows", nrow(validation_gap) == 8L, nrow(validation_gap))

pal_pub <- c(
  "FR" = "#0072B2", "AT" = "#D55E00", "CN" = "#009E73",
  "US" = "#CC79A7", "DE" = "#E69F00", "IT-A" = "#56B4E9",
  "IT-B" = "#A6761D", "JP" = "#666666"
)
theme_pub <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.7), color = NA
      ),
      plot.title.position = "plot"
    )
}
label_fixed <- function(digits = 2L) {
  function(x) formatC(x, format = "f", digits = digits)
}
save_pub <- function(plot, stem, width, height, units = "mm", dpi = 350L) {
  base <- file.path(figure_dir, stem)
  ggplot2::ggsave(
    paste0(base, ".pdf"), plot, width = width, height = height,
    units = units, device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    paste0(base, ".png"), plot, width = width, height = height,
    units = units, dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".tiff"), plot, width = width, height = height,
    units = units, dpi = dpi, compression = "lzw", bg = "white"
  )
  invisible(plot)
}

signature_display <- meta_signature[
  is.finite(meta_signature$QValue), , drop = FALSE
]
signature_display <- head(signature_display, 15L)
effect_display <- association_effects[
  association_effects$Feature %in% signature_display$Feature,
  c("Feature", "Species", "Cohort", "HedgesG"), drop = FALSE
]
pooled_display <- signature_display[, c("Feature", "Species", "PooledHedgesG")]
names(pooled_display)[[3L]] <- "HedgesG"
pooled_display$Cohort <- "Discovery meta"
effect_display <- rbind(effect_display, pooled_display[, names(effect_display)])
effect_display$Cohort <- factor(
  effect_display$Cohort,
  levels = c(discovery_cohorts, "Discovery meta", validation_cohorts)
)
effect_display$Species <- factor(
  effect_display$Species, levels = rev(signature_display$Species)
)
effect_limit <- max(1, stats::quantile(abs(effect_display$HedgesG), 0.95, na.rm = TRUE))
p_meta_signature <- ggplot2::ggplot(
  effect_display,
  ggplot2::aes(Cohort, Species, fill = pmax(-effect_limit, pmin(effect_limit, HedgesG)))
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 5.5, linewidth = 0.5, color = "grey25") +
  ggplot2::geom_vline(xintercept = 6.5, linewidth = 0.5, color = "grey25") +
  ggplot2::scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00",
    midpoint = 0, limits = c(-effect_limit, effect_limit),
    oob = scales::squish
  ) +
  ggplot2::labs(
    title = "Discovery associations face three untouched cohorts",
    subtitle = "Top 15 discovery random-effects results; positive values indicate CRC enrichment",
    x = NULL, y = NULL, fill = "Hedges' g",
    caption = "FR–DE: discovery · IT-A, IT-B and JP: direction-only external audit"
  ) +
  theme_pub(9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
save_pub(p_meta_signature, "28-meta-signature", width = 210, height = 145)

single_transfer$TrainingCohort <- factor(
  single_transfer$TrainingCohort, levels = rev(cohort_levels)
)
single_transfer$TestCohort <- factor(
  single_transfer$TestCohort, levels = cohort_levels
)
p_transfer <- ggplot2::ggplot(
  single_transfer,
  ggplot2::aes(TestCohort, TrainingCohort, fill = AUROC)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", AUROC)), size = 3
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#D9D9D9", mid = "#F0E442", high = "#0072B2",
    midpoint = 0.65, limits = c(0.3, 1), oob = scales::squish,
    breaks = seq(0.4, 1, 0.1), labels = label_fixed(1),
    guide = ggplot2::guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(42, "mm"), barheight = grid::unit(4, "mm")
    )
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "A model can travel badly even when its diagonal looks strong",
    subtitle = "Diagonal: within-study 5-fold CV · Off-diagonal: untouched study transfer",
    x = "Test cohort", y = "Training cohort", fill = "AUROC",
    caption = "This 8 × 8 matrix is descriptive; the primary estimate uses seven-study training."
  ) +
  theme_pub(10) +
  ggplot2::theme(legend.position = "bottom")
save_pub(p_transfer, "28-single-study-transfer", width = 175, height = 155)

forest_rows <- rbind(
  data.frame(
    Cohort = lodo_performance$Cohort,
    Samples = lodo_performance$Samples,
    AUROC = lodo_performance$AUROC,
    Lower95 = lodo_performance$AUROCLower95,
    Upper95 = lodo_performance$AUROCUpper95,
    Type = "Cohort-specific",
    stringsAsFactors = FALSE
  ),
  data.frame(
    Cohort = "Random-effects pooled",
    Samples = sum(lodo_performance$Samples),
    AUROC = performance_meta$PooledAUROC,
    Lower95 = performance_meta$Lower95,
    Upper95 = performance_meta$Upper95,
    Type = "Pooled",
    stringsAsFactors = FALSE
  )
)
forest_rows$Cohort <- factor(
  forest_rows$Cohort,
  levels = rev(c(cohort_levels, "Random-effects pooled"))
)
p_forest <- ggplot2::ggplot(
  forest_rows,
  ggplot2::aes(AUROC, Cohort, color = Type)
) +
  ggplot2::geom_vline(xintercept = 0.5, linetype = 2, color = "grey55") +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = Lower95, xmax = Upper95),
    height = 0, linewidth = 0.6
  ) +
  ggplot2::geom_point(ggplot2::aes(size = Samples), alpha = 0.95) +
  ggplot2::scale_color_manual(
    values = c("Cohort-specific" = "#0072B2", "Pooled" = "#D55E00"),
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0.4, 1, 0.1), labels = label_fixed(1)
  ) +
  ggplot2::coord_cartesian(xlim = c(0.35, 1)) +
  ggplot2::labs(
    title = "Every point is an untouched cohort",
    subtitle = paste0(
      "Seven-cohort training; random-effects I² = ",
      sprintf("%.1f%%", performance_meta$I2Percent)
    ),
    x = "External-validation AUROC (95% CI)", y = NULL,
    size = "Subjects",
    caption = "Confidence intervals are cohort-specific DeLong intervals."
  ) +
  theme_pub(10) +
  ggplot2::theme(legend.position = "bottom")
save_pub(p_forest, "28-lodo-forest", width = 175, height = 132)

roc_data$Cohort <- factor(roc_data$Cohort, levels = cohort_levels)
roc_labels <- stats::setNames(
  sprintf(
    "%s · AUC %.2f", lodo_performance$Cohort, lodo_performance$AUROC
  ),
  lodo_performance$Cohort
)
p_roc <- ggplot2::ggplot(
  roc_data, ggplot2::aes(FPR, TPR, color = Cohort)
) +
  ggplot2::geom_abline(
    slope = 1, intercept = 0, linetype = 2, color = "grey65"
  ) +
  ggplot2::geom_step(linewidth = 0.75, direction = "hv") +
  ggplot2::facet_wrap(
    ~Cohort, ncol = 4,
    labeller = ggplot2::as_labeller(roc_labels)
  ) +
  ggplot2::scale_color_manual(values = pal_pub, guide = "none") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    title = "One prediction per subject, one complete cohort per panel",
    subtitle = "No outer-test sample entered filtering, standardization or lambda selection",
    x = "False-positive rate", y = "True-positive rate"
  ) +
  theme_pub(8.5) +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
save_pub(p_roc, "28-lodo-roc", width = 210, height = 115)

gap_long <- rbind(
  data.frame(
    Cohort = validation_gap$Cohort,
    Validation = "Within-study CV",
    AUROC = validation_gap$WithinCohortCV
  ),
  data.frame(
    Cohort = validation_gap$Cohort,
    Validation = "Mean one-study transfer",
    AUROC = validation_gap$MeanSingleStudyTransfer
  ),
  data.frame(
    Cohort = validation_gap$Cohort,
    Validation = "Seven-study LODO",
    AUROC = validation_gap$NMinusOneLODO
  )
)
gap_long$Cohort <- factor(gap_long$Cohort, levels = cohort_levels)
gap_long$Validation <- factor(
  gap_long$Validation,
  levels = c("Within-study CV", "Mean one-study transfer", "Seven-study LODO")
)
p_gap <- ggplot2::ggplot(
  gap_long,
  ggplot2::aes(Cohort, AUROC, group = Cohort, color = Validation)
) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = 2, color = "grey60") +
  ggplot2::geom_line(ggplot2::aes(group = Cohort), color = "grey75", linewidth = 0.5) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::scale_color_manual(values = c(
    "Within-study CV" = "#999999",
    "Mean one-study transfer" = "#E69F00",
    "Seven-study LODO" = "#0072B2"
  )) +
  ggplot2::scale_y_continuous(
    breaks = seq(0.4, 1, 0.1), labels = label_fixed(1)
  ) +
  ggplot2::coord_cartesian(ylim = c(0.35, 1)) +
  ggplot2::labs(
    title = "Validation design changes the answer",
    x = NULL, y = "AUROC", color = NULL
  ) +
  theme_pub(8.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

stable_features <- head(coefficient_stability$Feature, 15L)
coefficient_display <- expand.grid(
  Feature = stable_features,
  OuterTestCohort = cohort_levels,
  stringsAsFactors = FALSE
)
coefficient_display <- merge(
  coefficient_display,
  lodo_coefficients[, c("OuterTestCohort", "Feature", "Coefficient")],
  by = c("OuterTestCohort", "Feature"), all.x = TRUE, sort = FALSE
)
coefficient_display$Coefficient[is.na(coefficient_display$Coefficient)] <- 0
max_by_model <- tapply(
  abs(coefficient_display$Coefficient),
  coefficient_display$OuterTestCohort, max
)
coefficient_display$ScaledCoefficient <- coefficient_display$Coefficient /
  pmax(max_by_model[coefficient_display$OuterTestCohort], 1e-12)
coefficient_display$Species <- species_labels[coefficient_display$Feature]
coefficient_display$Species <- factor(
  coefficient_display$Species,
  levels = rev(coefficient_stability$Species[match(
    stable_features, coefficient_stability$Feature
  )])
)
coefficient_display$OuterTestCohort <- factor(
  coefficient_display$OuterTestCohort, levels = cohort_levels
)
p_stability <- ggplot2::ggplot(
  coefficient_display,
  ggplot2::aes(OuterTestCohort, Species, fill = ScaledCoefficient)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00",
    midpoint = 0, limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = "Sparse signatures are not identical",
    subtitle = "Top 15 by selection frequency; columns name the held-out cohort",
    x = "Outer test cohort", y = NULL, fill = "Scaled\ncoefficient"
  ) +
  theme_pub(8.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

p_stability_combined <- p_gap + p_stability +
  patchwork::plot_layout(widths = c(0.78, 1.35)) +
  patchwork::plot_annotation(
    title = "Generalization improves with cohort diversity, but the selected taxa still move",
    caption = "Coefficient colors are normalized within each outer model; magnitude is not comparable across columns.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.caption = ggplot2::element_text(color = "grey35", size = 8.5)
    )
  )
save_pub(p_stability_combined, "28-validation-stability", width = 240, height = 145)

figure_stems <- c(
  "28-meta-signature", "28-single-study-transfer", "28-lodo-forest",
  "28-lodo-roc", "28-validation-stability"
)
figure_audit <- do.call(
  rbind,
  lapply(figure_stems, function(stem) {
    do.call(
      rbind,
      lapply(c("pdf", "png", "tiff"), function(extension) {
        path <- file.path(figure_dir, paste0(stem, ".", extension))
        data.frame(
          Figure = stem,
          Format = extension,
          Exists = file.exists(path),
          Bytes = if (file.exists(path)) file.info(path)$size else NA_real_,
          SHA256 = if (file.exists(path)) sha256_file(path) else NA_character_,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
add_check("Figures", "figure-files", nrow(figure_audit) == 15L && all(figure_audit$Exists), paste0(sum(figure_audit$Exists), "/15"))
add_check("Figures", "figure-nonempty", all(figure_audit$Bytes > 10000), paste(range(figure_audit$Bytes), collapse = "/"))
anchor_path <- file.path(figure_dir, "28-wirbel-fig3-original.png")
anchor_hash <- if (file.exists(anchor_path)) sha256_file(anchor_path) else "missing"
add_check("Figures", "anchor-present", file.exists(anchor_path) && file.info(anchor_path)$size > 100000, anchor_path)
add_check(
  "Figures", "anchor-sha256",
  identical(anchor_hash, "e430cb0afefaed0b1c33333350758fac74c1d165cfa484ef7984c15ed747f43f"),
  anchor_hash
)

chapter_text <- paste(readLines(chapter_path, warn = FALSE), collapse = "\n")
for (token in c(
  "draft: false", "eval: true", "freeze: auto", "expected_images: 6",
  "真正的 N−1", "七个完整队列", "one-standard-error",
  "Random-effects", "ComBat", "MetaPhlAn 3", "CHOCOPhlAn 201901",
  "e430cb0afefaed0b1c33333350758fac74c1d165cfa484ef7984c15ed747f43f",
  "不属于本仓库 CC BY/MIT 授权内容",
  "## 审计与升级", "## 出版级美化", "## 常见坑",
  "## 这段 Methods 怎么写", "## 换成你自己的数据怎么做", "## 参考"
)) {
  add_check(
    "Chapter", paste0("chapter-", gsub("[^a-z0-9]+", "-", tolower(token))),
    grepl(token, chapter_text, fixed = TRUE), token
  )
}
add_check("Chapter", "no-source-theme-dependency", !grepl("source(\"R/theme_pub.R\")", chapter_text, fixed = TRUE), "inline plotting functions")
add_check("Chapter", "no-draft-meta-copy", !grepl("本篇可独立跑通", chapter_text, fixed = TRUE), "reader-facing copy")

write_tsv_gz(association_effects, file.path(output_dir, "association-cohort-effects.tsv.gz"))
write_tsv(association_meta, file.path(output_dir, "association-meta-analysis.tsv"))
write_tsv(meta_signature, file.path(output_dir, "meta-signature-validation.tsv"))
write_tsv_gz(lodo_predictions, file.path(output_dir, "lodo-predictions.tsv.gz"))
write_tsv(lodo_performance, file.path(output_dir, "lodo-performance.tsv"))
write_tsv_gz(lodo_tuning, file.path(output_dir, "lodo-inner-tuning.tsv.gz"))
write_tsv(tuning_summary, file.path(output_dir, "lodo-tuning-summary.tsv"))
write_tsv(lodo_preprocessing, file.path(output_dir, "lodo-preprocessing-audit.tsv"))
write_tsv_gz(lodo_coefficients, file.path(output_dir, "lodo-coefficients.tsv.gz"))
write_tsv(coefficient_stability, file.path(output_dir, "lodo-coefficient-stability.tsv"))
write_tsv(performance_meta, file.path(output_dir, "performance-meta-analysis.tsv"))
write_tsv(bootstrap_summary, file.path(output_dir, "hierarchical-bootstrap-summary.tsv"))
write_tsv_gz(bootstrap_table, file.path(output_dir, "hierarchical-bootstrap-replicates.tsv.gz"))
write_tsv(roc_data, file.path(output_dir, "lodo-roc-curves.tsv"))
write_tsv(single_transfer, file.path(output_dir, "single-study-transfer.tsv"))
write_tsv_gz(single_transfer_tuning, file.path(output_dir, "single-study-tuning.tsv.gz"))
write_tsv(validation_gap, file.path(output_dir, "validation-gap.tsv"))
write_tsv(figure_audit, file.path(output_dir, "figure-audit.tsv"))
write_tsv(checks, file.path(output_dir, "validation-checks.tsv"))

failures <- checks[checks$Status != "PASS", , drop = FALSE]
summary <- list(
  status = if (nrow(failures) == 0L) "passed" else "failed",
  article = 28L,
  seed = primary_seed,
  cohorts = length(cohort_levels),
  samples = nrow(x),
  controls = sum(y == "Control"),
  crc = sum(y == "CRC"),
  raw_species_features = ncol(x),
  discovery_meta_features = length(meta_features),
  discovery_q_lt_005 = sum(meta_signature$QValue < 0.05, na.rm = TRUE),
  lodo_models = nrow(lodo_performance),
  lodo_macro_auroc = mean(lodo_performance$AUROC),
  lodo_macro_bootstrap_lower95 = bootstrap_summary$Lower95,
  lodo_macro_bootstrap_upper95 = bootstrap_summary$Upper95,
  pooled_auroc = performance_meta$PooledAUROC,
  pooled_auroc_lower95 = performance_meta$Lower95,
  pooled_auroc_upper95 = performance_meta$Upper95,
  performance_i2_percent = performance_meta$I2Percent,
  stable_features_selected_in_all_models = sum(coefficient_stability$ModelsSelected == 8L),
  transfer_matrix_cells = nrow(single_transfer),
  checksum_entries = checksum_entries,
  checks = nrow(checks),
  passed = sum(checks$Status == "PASS"),
  failed = nrow(failures),
  generated_figure_files = nrow(figure_audit),
  package_versions = list(
    glmnet = as.character(utils::packageVersion("glmnet")),
    metafor = as.character(utils::packageVersion("metafor")),
    pROC = as.character(utils::packageVersion("pROC")),
    R = paste(R.version$major, R.version$minor, sep = ".")
  ),
  generated_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
jsonlite::write_json(
  summary, file.path(output_dir, "validation-summary.json"),
  auto_unbox = TRUE, pretty = TRUE
)

log_msg(
  "checks=", nrow(checks), "; passed=", sum(checks$Status == "PASS"),
  "; failed=", nrow(failures)
)
if (nrow(failures) > 0L) {
  print(failures, row.names = FALSE)
  stop("Article 28 validation failed.", call. = FALSE)
}
log_msg("Article 28 validation passed.")
cat("Article 28 validation passed:", nrow(checks), "checks.\n")
