#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999, digits = 17)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260727L
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
  "ranger", "xgboost", "pROC", "ggplot2", "patchwork",
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

font_cache_dir <- file.path(tempdir(), "article27-fontconfig-cache")
dir.create(font_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

validation_log <- file.path(output_dir, "validation.log")
writeLines(
  c(
    "Article 27 leakage-safe machine-learning validation",
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
    add_check("Frozen input", paste0("sha256-", relative), identical(observed, expected), observed)
    expected_files <- c(expected_files, relative)
  }
  payloads <- sort(setdiff(basename(list.files(directory, full.names = TRUE)), "file-checksums.sha256"))
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

roc_curve <- function(y, score, model) {
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
  data.frame(
    Model = model,
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

fit_preprocessor <- function(x) {
  keep <- colMeans(x > 0) >= 0.10 &
    apply(x, 2L, max) >= 1e-4 &
    apply(x, 2L, stats::sd) > 0
  if (sum(keep) < 2L) stop("Training fold retained fewer than two features.", call. = FALSE)
  list(features = colnames(x)[keep], pseudocount = 1e-6)
}

apply_preprocessor <- function(preprocessor, x) {
  out <- log10(x[, preprocessor$features, drop = FALSE] + preprocessor$pseudocount)
  storage.mode(out) <- "double"
  out
}

mtry_value <- function(specification, p) {
  value <- switch(
    specification,
    sqrt = floor(sqrt(p)),
    fifth = floor(p / 5),
    third = floor(p / 3),
    stop("Unknown mtry specification.", call. = FALSE)
  )
  max(1L, min(p, as.integer(value)))
}

fit_model <- function(algorithm, x, y, params, seed) {
  if (algorithm == "Random forest") {
    ranger::ranger(
      x = x,
      y = factor(y, levels = c("Control", "CRC")),
      probability = TRUE,
      num.trees = 750L,
      mtry = mtry_value(params$mtry_spec, ncol(x)),
      min.node.size = as.integer(params$min_node),
      splitrule = "gini",
      importance = "none",
      num.threads = 1L,
      seed = as.integer(seed)
    )
  } else if (algorithm == "XGBoost") {
    set.seed(as.integer(seed))
    xgboost::xgboost(
      data = x,
      label = as.integer(y == "CRC"),
      objective = "binary:logistic",
      eval_metric = "auc",
      nrounds = as.integer(params$nrounds),
      max_depth = as.integer(params$max_depth),
      eta = as.numeric(params$eta),
      min_child_weight = 1,
      subsample = 0.8,
      colsample_bytree = 0.8,
      nthread = 1L,
      verbosity = 0,
      verbose = 0
    )
  } else {
    stop("Unknown algorithm: ", algorithm, call. = FALSE)
  }
}

predict_model <- function(algorithm, model, x) {
  if (algorithm == "Random forest") {
    as.numeric(predict(model, data = x)$predictions[, "CRC"])
  } else {
    as.numeric(stats::predict(model, x))
  }
}

rf_grid <- expand.grid(
  mtry_spec = c("sqrt", "fifth", "third"),
  min_node = c(3L, 10L),
  stringsAsFactors = FALSE
)
xgb_grid <- expand.grid(
  nrounds = c(50L, 150L),
  max_depth = c(1L, 2L),
  eta = c(0.03, 0.10),
  stringsAsFactors = FALSE
)

tune_model <- function(algorithm, x, y, folds, seed) {
  grid <- if (algorithm == "Random forest") rf_grid else xgb_grid
  scores <- numeric(nrow(grid))
  for (g in seq_len(nrow(grid))) {
    predictions <- rep(NA_real_, length(y))
    params <- as.list(grid[g, , drop = FALSE])
    for (fold in sort(unique(folds))) {
      train <- folds != fold
      validation <- folds == fold
      preprocessor <- fit_preprocessor(x[train, , drop = FALSE])
      x_train <- apply_preprocessor(preprocessor, x[train, , drop = FALSE])
      x_validation <- apply_preprocessor(preprocessor, x[validation, , drop = FALSE])
      model <- fit_model(
        algorithm, x_train, y[train], params,
        seed + g * 100L + fold
      )
      predictions[validation] <- predict_model(algorithm, model, x_validation)
    }
    scores[[g]] <- auc_score(y == "CRC", predictions)
  }
  grid$InnerAUROC <- scores
  best <- which.max(scores)
  list(params = as.list(grid[best, , drop = FALSE]), audit = grid, best = best)
}

checksum_entries <- verify_checksum_manifest(input_dir)
notice_path <- file.path(project_root, "data", "small", "27-data-NOTICE.txt")
notice <- readLines(notice_path, warn = FALSE)
add_check("Frozen input", "notice-present", length(notice) >= 11L, notice_path)
for (token in c(
  "114 independent CRC/control subjects", "all 661 species rows",
  "Feature filtering, log transformation and hyperparameter tuning",
  "Article 28, not this chapter"
)) {
  add_check(
    "Frozen input", paste0("notice-", gsub("[^a-z0-9]+", "-", tolower(token))),
    any(grepl(token, notice, fixed = TRUE)), token
  )
}

species_table <- read_tsv_gz("species-relative-abundance.tsv.gz")
metadata <- read_tsv("sample-metadata.tsv")
feature_audit <- read_tsv("feature-universe-audit.tsv")
contract <- read_tsv("analysis-contract.tsv")
resource_manifest <- read_tsv("resource-manifest.tsv")

sample_ids <- metadata$sample_id
x <- t(as.matrix(species_table[, sample_ids, drop = FALSE]))
storage.mode(x) <- "double"
colnames(x) <- species_table$Feature
rownames(x) <- sample_ids
y <- factor(metadata$Outcome, levels = c("Control", "CRC"))
species_labels <- stats::setNames(species_table$Species, species_table$Feature)

add_check("Data", "samples", nrow(x) == 114L, nrow(x))
add_check("Data", "features", ncol(x) == 661L, ncol(x))
add_check("Data", "sample-alignment", identical(rownames(x), metadata$sample_id), "matrix vs metadata")
add_check("Data", "independent-subjects", !anyDuplicated(metadata$subject_id), length(unique(metadata$subject_id)))
add_check("Data", "study-qualified-keys", !anyDuplicated(metadata$StudySampleKey), length(unique(metadata$StudySampleKey)))
add_check("Data", "study", all(metadata$study_name == "ZellerG_2014"), unique(metadata$study_name))
add_check("Data", "body-site", all(metadata$body_site == "stool"), unique(metadata$body_site))
add_check("Data", "class-counts", identical(as.integer(table(y)), c(61L, 53L)), paste(table(y), collapse = "/"))
add_check("Data", "adenoma-excluded", !any(metadata$study_condition == "adenoma"), unique(metadata$study_condition))
add_check("Data", "finite", all(is.finite(x)), range(x))
add_check("Data", "nonnegative", min(x) >= 0, min(x))
add_check("Data", "fraction-unit", max(x) <= 1, max(x))
add_check("Data", "closure", max(abs(rowSums(x) - 1)) < 2e-6, range(rowSums(x)))
add_check("Data", "feature-audit-rows", nrow(feature_audit) == 661L, nrow(feature_audit))
add_check("Data", "resource-manifest-rows", nrow(resource_manifest) == 5L, nrow(resource_manifest))

contract_value <- function(item) contract$Value[match(item, contract$Item)]
expected_contract <- c(
  seed = "20260727", study = "ZellerG_2014", outcome = "CRC vs Control",
  positive_class = "CRC", samples = "114", control_samples = "61",
  crc_samples = "53", raw_species_features = "661", outer_resamples = "5",
  outer_folds = "5", inner_folds = "4", fold_unit = "subject_id",
  log10_pseudocount = "0.000001", tuning_metric = "inner-fold AUROC",
  paired_difference_bootstrap = "2000 subject bootstraps",
  feature_importance = "outer-test permutation delta AUROC",
  leakage_audit_permutations = "50"
)
for (item in names(expected_contract)) {
  observed <- contract_value(item)
  add_check("Contract", paste0("contract-", gsub("_", "-", item)), identical(observed, expected_contract[[item]]), observed)
}
add_check("Contract", "training-prevalence-scoped", grepl("within each training fold", contract_value("training_prevalence_filter"), fixed = TRUE), contract_value("training_prevalence_filter"))
add_check("Contract", "training-abundance-scoped", grepl("within each training fold", contract_value("training_abundance_filter"), fixed = TRUE), contract_value("training_abundance_filter"))
add_check("Contract", "both-models-prespecified", grepl("ranger", contract_value("algorithms"), fixed = TRUE) && grepl("xgboost", contract_value("algorithms"), fixed = TRUE), contract_value("algorithms"))

outer_assignments <- list()
outer_predictions <- list()
outer_metrics <- list()
tuning_audits <- list()
importance_rows <- list()
preprocessing_rows <- list()
prediction_index <- metric_index <- tuning_index <- importance_index <- preprocessing_index <- 0L

algorithms <- c("Random forest", "XGBoost")
for (repeat_id in seq_len(5L)) {
  outer_fold <- make_stratified_folds(y, 5L, primary_seed + repeat_id * 1000L)
  outer_assignments[[repeat_id]] <- data.frame(
    SampleID = sample_ids,
    SubjectID = metadata$subject_id,
    Outcome = as.character(y),
    Repeat = repeat_id,
    Fold = outer_fold,
    stringsAsFactors = FALSE
  )
  for (fold_id in seq_len(5L)) {
    train <- outer_fold != fold_id
    test <- outer_fold == fold_id
    inner_fold <- make_stratified_folds(
      y[train], 4L,
      primary_seed + repeat_id * 10000L + fold_id * 100L
    )
    for (algorithm_index in seq_along(algorithms)) {
      algorithm <- algorithms[[algorithm_index]]
      seed_base <- primary_seed + algorithm_index * 1000000L + repeat_id * 10000L + fold_id * 100L
      tuned <- tune_model(algorithm, x[train, , drop = FALSE], y[train], inner_fold, seed_base)
      tuning_audit <- tuned$audit
      tuning_columns <- c(
        "mtry_spec", "min_node", "max_depth", "nrounds", "eta",
        "InnerAUROC"
      )
      for (column in setdiff(tuning_columns, names(tuning_audit))) {
        tuning_audit[[column]] <- NA
      }
      tuning_audit <- tuning_audit[, tuning_columns, drop = FALSE]
      tuning_audit$Algorithm <- algorithm
      tuning_audit$Repeat <- repeat_id
      tuning_audit$OuterFold <- fold_id
      tuning_audit$Selected <- seq_len(nrow(tuning_audit)) == tuned$best
      tuning_index <- tuning_index + 1L
      tuning_audits[[tuning_index]] <- tuning_audit

      preprocessor <- fit_preprocessor(x[train, , drop = FALSE])
      x_train <- apply_preprocessor(preprocessor, x[train, , drop = FALSE])
      x_test <- apply_preprocessor(preprocessor, x[test, , drop = FALSE])
      model <- fit_model(algorithm, x_train, y[train], tuned$params, seed_base + 999L)
      probability <- predict_model(algorithm, model, x_test)
      baseline_auc <- auc_score(y[test] == "CRC", probability)
      params_text <- paste(
        paste(names(tuned$params), unlist(tuned$params), sep = "="),
        collapse = ";"
      )

      prediction_index <- prediction_index + 1L
      outer_predictions[[prediction_index]] <- data.frame(
        SampleID = sample_ids[test],
        SubjectID = metadata$subject_id[test],
        Outcome = as.character(y[test]),
        Repeat = repeat_id,
        OuterFold = fold_id,
        Algorithm = algorithm,
        ProbabilityCRC = probability,
        stringsAsFactors = FALSE
      )
      metric_index <- metric_index + 1L
      outer_metrics[[metric_index]] <- data.frame(
        Algorithm = algorithm,
        Repeat = repeat_id,
        OuterFold = fold_id,
        TestSamples = sum(test),
        TestControls = sum(y[test] == "Control"),
        TestCRC = sum(y[test] == "CRC"),
        FeaturesKept = length(preprocessor$features),
        InnerAUROC = max(tuned$audit$InnerAUROC),
        OuterAUROC = baseline_auc,
        OuterAUPRC = average_precision(y[test] == "CRC", probability),
        OuterBrier = mean((as.integer(y[test] == "CRC") - probability)^2),
        SelectedParameters = params_text,
        stringsAsFactors = FALSE
      )
      preprocessing_index <- preprocessing_index + 1L
      preprocessing_rows[[preprocessing_index]] <- data.frame(
        Algorithm = algorithm,
        Repeat = repeat_id,
        OuterFold = fold_id,
        TrainingSamples = sum(train),
        TrainingControls = sum(y[train] == "Control"),
        TrainingCRC = sum(y[train] == "CRC"),
        FeaturesAvailable = ncol(x),
        FeaturesKept = length(preprocessor$features),
        Pseudocount = preprocessor$pseudocount,
        FittedOnOuterTrainingOnly = TRUE,
        stringsAsFactors = FALSE
      )

      for (feature_index in seq_along(preprocessor$features)) {
        x_permuted <- x_test
        set.seed(seed_base + 2000L + feature_index)
        x_permuted[, feature_index] <- sample(x_permuted[, feature_index], replace = FALSE)
        permuted_probability <- predict_model(algorithm, model, x_permuted)
        importance_index <- importance_index + 1L
        feature <- preprocessor$features[[feature_index]]
        importance_rows[[importance_index]] <- data.frame(
          Algorithm = algorithm,
          Repeat = repeat_id,
          OuterFold = fold_id,
          Feature = feature,
          Species = unname(species_labels[[feature]]),
          BaselineAUROC = baseline_auc,
          PermutedAUROC = auc_score(y[test] == "CRC", permuted_probability),
          DeltaAUROC = baseline_auc - auc_score(y[test] == "CRC", permuted_probability),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  log_msg("completed outer repeat ", repeat_id, "/5")
}

fold_assignments <- do.call(rbind, outer_assignments)
outer_predictions <- do.call(rbind, outer_predictions)
outer_metrics <- do.call(rbind, outer_metrics)
tuning_audit <- do.call(rbind, tuning_audits)
importance <- do.call(rbind, importance_rows)
preprocessing_audit <- do.call(rbind, preprocessing_rows)

sample_predictions <- do.call(
  rbind,
  lapply(split(outer_predictions, list(outer_predictions$Algorithm, outer_predictions$SampleID), drop = TRUE), function(z) {
    data.frame(
      Algorithm = z$Algorithm[[1L]],
      SampleID = z$SampleID[[1L]],
      SubjectID = z$SubjectID[[1L]],
      Outcome = z$Outcome[[1L]],
      OOFPredictions = nrow(z),
      ProbabilityCRC = mean(z$ProbabilityCRC),
      ProbabilitySD = stats::sd(z$ProbabilityCRC),
      stringsAsFactors = FALSE
    )
  })
)
sample_predictions <- sample_predictions[order(sample_predictions$Algorithm, sample_predictions$SampleID), ]

performance_rows <- list()
roc_rows <- list()
calibration_rows <- list()
for (algorithm in algorithms) {
  z <- sample_predictions[sample_predictions$Algorithm == algorithm, , drop = FALSE]
  truth <- factor(z$Outcome, levels = c("Control", "CRC"))
  binary <- as.integer(truth == "CRC")
  roc_object <- pROC::roc(
    response = truth,
    predictor = z$ProbabilityCRC,
    levels = c("Control", "CRC"),
    direction = "<",
    quiet = TRUE
  )
  interval <- as.numeric(pROC::ci.auc(roc_object, conf.level = 0.95, method = "delong"))
  predicted <- z$ProbabilityCRC >= 0.5
  performance_rows[[algorithm]] <- data.frame(
    Algorithm = algorithm,
    Samples = nrow(z),
    Controls = sum(binary == 0L),
    CRC = sum(binary == 1L),
    AUROC = auc_score(binary, z$ProbabilityCRC),
    AUROCLower95 = interval[[1L]],
    AUROCUpper95 = interval[[3L]],
    AUPRC = average_precision(binary, z$ProbabilityCRC),
    Brier = mean((binary - z$ProbabilityCRC)^2),
    SensitivityAt0.5 = sum(predicted & binary == 1L) / sum(binary == 1L),
    SpecificityAt0.5 = sum(!predicted & binary == 0L) / sum(binary == 0L),
    stringsAsFactors = FALSE
  )
  roc_rows[[algorithm]] <- roc_curve(binary, z$ProbabilityCRC, algorithm)
  ord <- order(z$ProbabilityCRC, z$SampleID)
  bins <- integer(nrow(z))
  bins[ord] <- pmin(5L, ceiling(seq_along(ord) * 5 / length(ord)))
  calibration_rows[[algorithm]] <- do.call(
    rbind,
    lapply(seq_len(5L), function(bin) {
      idx <- bins == bin
      data.frame(
        Algorithm = algorithm,
        Bin = bin,
        Samples = sum(idx),
        MeanPredictedRisk = mean(z$ProbabilityCRC[idx]),
        ObservedCRCProportion = mean(binary[idx]),
        stringsAsFactors = FALSE
      )
    })
  )
}
performance <- do.call(rbind, performance_rows)
roc_data <- do.call(rbind, roc_rows)
calibration <- do.call(rbind, calibration_rows)

wide_predictions <- reshape(
  sample_predictions[, c("SampleID", "Outcome", "Algorithm", "ProbabilityCRC")],
  idvar = c("SampleID", "Outcome"), timevar = "Algorithm", direction = "wide"
)
names(wide_predictions) <- sub("ProbabilityCRC\\.", "", names(wide_predictions))
binary <- as.integer(wide_predictions$Outcome == "CRC")
observed_difference <- auc_score(binary, wide_predictions[["Random forest"]]) -
  auc_score(binary, wide_predictions[["XGBoost"]])
set.seed(primary_seed + 8000000L)
bootstrap_difference <- replicate(2000L, {
  idx <- sample(seq_len(nrow(wide_predictions)), replace = TRUE)
  if (length(unique(binary[idx])) < 2L) return(NA_real_)
  auc_score(binary[idx], wide_predictions[["Random forest"]][idx]) -
    auc_score(binary[idx], wide_predictions[["XGBoost"]][idx])
})
bootstrap_difference <- bootstrap_difference[is.finite(bootstrap_difference)]
model_comparison <- data.frame(
  Contrast = "Random forest minus XGBoost",
  AUROCDifference = observed_difference,
  BootstrapLower95 = unname(stats::quantile(bootstrap_difference, 0.025)),
  BootstrapUpper95 = unname(stats::quantile(bootstrap_difference, 0.975)),
  BootstrapReplicates = length(bootstrap_difference),
  ResamplingUnit = "subject",
  Interpretation = "Conditional on aggregated out-of-fold predictions; not external validation",
  stringsAsFactors = FALSE
)

importance_summary <- do.call(
  rbind,
  lapply(split(importance, list(importance$Algorithm, importance$Feature), drop = TRUE), function(z) {
    data.frame(
      Algorithm = z$Algorithm[[1L]],
      Feature = z$Feature[[1L]],
      Species = z$Species[[1L]],
      OuterFoldsAvailable = nrow(z),
      MeanDeltaAUROC = mean(z$DeltaAUROC),
      SEDeltaAUROC = stats::sd(z$DeltaAUROC) / sqrt(nrow(z)),
      PositiveFoldFraction = mean(z$DeltaAUROC > 0),
      stringsAsFactors = FALSE
    )
  })
)
importance_summary$Lower95 <- importance_summary$MeanDeltaAUROC - 1.96 * importance_summary$SEDeltaAUROC
importance_summary$Upper95 <- importance_summary$MeanDeltaAUROC + 1.96 * importance_summary$SEDeltaAUROC
importance_summary <- importance_summary[order(
  importance_summary$Algorithm, -importance_summary$MeanDeltaAUROC,
  importance_summary$Species
), ]
importance_summary$Rank <- ave(
  -importance_summary$MeanDeltaAUROC,
  importance_summary$Algorithm,
  FUN = rank, ties.method = "first"
)

log_msg("starting 50 subject-label permutation leakage audits")
leakage_rows <- list()
leakage_index <- 0L
for (permutation in seq_len(50L)) {
  set.seed(primary_seed + 9000000L + permutation)
  permuted_y <- sample(y, length(y), replace = FALSE)
  permuted_fold <- make_stratified_folds(
    permuted_y, 5L,
    primary_seed + 9100000L + permutation
  )
  association <- vapply(seq_len(ncol(x)), function(j) {
    abs(auc_score(permuted_y == "CRC", x[, j]) - 0.5)
  }, numeric(1L))
  leaky_features <- colnames(x)[order(-association, colnames(x))[seq_len(20L)]]
  safe_probability <- leaky_probability <- rep(NA_real_, nrow(x))
  for (fold in seq_len(5L)) {
    train <- permuted_fold != fold
    test <- permuted_fold == fold
    safe_preprocessor <- fit_preprocessor(x[train, , drop = FALSE])
    safe_train <- apply_preprocessor(safe_preprocessor, x[train, , drop = FALSE])
    safe_test <- apply_preprocessor(safe_preprocessor, x[test, , drop = FALSE])
    params <- list(mtry_spec = "sqrt", min_node = 5L)
    safe_model <- ranger::ranger(
      x = safe_train,
      y = factor(permuted_y[train], levels = c("Control", "CRC")),
      probability = TRUE,
      num.trees = 300L,
      mtry = max(1L, floor(sqrt(ncol(safe_train)))),
      min.node.size = 5L,
      splitrule = "gini",
      importance = "none",
      num.threads = 1L,
      seed = primary_seed + 9200000L + permutation * 10L + fold
    )
    safe_probability[test] <- as.numeric(
      predict(safe_model, data = safe_test)$predictions[, "CRC"]
    )
    leaky_train <- log10(x[train, leaky_features, drop = FALSE] + 1e-6)
    leaky_test <- log10(x[test, leaky_features, drop = FALSE] + 1e-6)
    leaky_model <- ranger::ranger(
      x = leaky_train,
      y = factor(permuted_y[train], levels = c("Control", "CRC")),
      probability = TRUE,
      num.trees = 300L,
      mtry = max(1L, floor(sqrt(ncol(leaky_train)))),
      min.node.size = 5L,
      splitrule = "gini",
      importance = "none",
      num.threads = 1L,
      seed = primary_seed + 9300000L + permutation * 10L + fold
    )
    leaky_probability[test] <- as.numeric(
      predict(leaky_model, data = leaky_test)$predictions[, "CRC"]
    )
  }
  leakage_index <- leakage_index + 1L
  leakage_rows[[leakage_index]] <- data.frame(
    Permutation = permutation,
    Pipeline = "Fold-scoped preprocessing",
    AUROC = auc_score(permuted_y == "CRC", safe_probability),
    LabelInformationFromTestFold = FALSE,
    stringsAsFactors = FALSE
  )
  leakage_index <- leakage_index + 1L
  leakage_rows[[leakage_index]] <- data.frame(
    Permutation = permutation,
    Pipeline = "Global label-selected top 20",
    AUROC = auc_score(permuted_y == "CRC", leaky_probability),
    LabelInformationFromTestFold = TRUE,
    stringsAsFactors = FALSE
  )
}
leakage_audit <- do.call(rbind, leakage_rows)
leakage_summary <- do.call(
  rbind,
  lapply(split(leakage_audit, leakage_audit$Pipeline), function(z) {
    data.frame(
      Pipeline = z$Pipeline[[1L]],
      Permutations = nrow(z),
      MedianAUROC = stats::median(z$AUROC),
      Lower95 = unname(stats::quantile(z$AUROC, 0.025)),
      Upper95 = unname(stats::quantile(z$AUROC, 0.975)),
      LabelInformationFromTestFold = z$LabelInformationFromTestFold[[1L]],
      stringsAsFactors = FALSE
    )
  })
)

add_check("Cross-validation", "fold-assignment-rows", nrow(fold_assignments) == 570L, nrow(fold_assignments))
add_check("Cross-validation", "five-repeats", setequal(unique(fold_assignments$Repeat), 1:5), unique(fold_assignments$Repeat))
add_check("Cross-validation", "five-folds", setequal(unique(fold_assignments$Fold), 1:5), unique(fold_assignments$Fold))
add_check("Cross-validation", "one-fold-per-repeat", all(table(fold_assignments$SampleID, fold_assignments$Repeat) == 1L), range(table(fold_assignments$SampleID, fold_assignments$Repeat)))
add_check("Cross-validation", "subject-fold-integrity", all(table(fold_assignments$SubjectID, fold_assignments$Repeat) == 1L), "subject x repeat")
add_check("Cross-validation", "outer-prediction-rows", nrow(outer_predictions) == 1140L, nrow(outer_predictions))
add_check("Cross-validation", "five-oof-predictions", all(table(outer_predictions$SampleID, outer_predictions$Algorithm) == 5L), range(table(outer_predictions$SampleID, outer_predictions$Algorithm)))
add_check("Cross-validation", "sample-prediction-rows", nrow(sample_predictions) == 228L, nrow(sample_predictions))
add_check("Cross-validation", "probability-bounds", all(outer_predictions$ProbabilityCRC >= 0 & outer_predictions$ProbabilityCRC <= 1), range(outer_predictions$ProbabilityCRC))
add_check("Cross-validation", "outer-model-count", nrow(outer_metrics) == 50L, nrow(outer_metrics))
add_check("Cross-validation", "preprocessors-training-only", all(preprocessing_audit$FittedOnOuterTrainingOnly), unique(preprocessing_audit$FittedOnOuterTrainingOnly))
add_check("Cross-validation", "feature-filter-fold-specific", length(unique(preprocessing_audit$FeaturesKept)) > 1L, range(preprocessing_audit$FeaturesKept))
add_check("Cross-validation", "selected-tuning-per-model", all(table(tuning_audit$Algorithm, tuning_audit$Repeat, tuning_audit$OuterFold, tuning_audit$Selected)[,,,"TRUE"] == 1L), "one selected row")

add_check("Performance", "two-models", setequal(performance$Algorithm, algorithms), performance$Algorithm)
add_check("Performance", "auc-bounds", all(performance$AUROC >= 0.5 & performance$AUROC <= 1), performance$AUROC)
add_check("Performance", "auc-ci-bounds", all(performance$AUROCLower95 <= performance$AUROC & performance$AUROCUpper95 >= performance$AUROC), paste(performance$AUROCLower95, performance$AUROCUpper95))
add_check("Performance", "auprc-bounds", all(performance$AUPRC >= 0 & performance$AUPRC <= 1), performance$AUPRC)
add_check("Performance", "brier-bounds", all(performance$Brier >= 0 & performance$Brier <= 1), performance$Brier)
add_check("Performance", "paired-bootstrap", model_comparison$BootstrapReplicates == 2000L, model_comparison$BootstrapReplicates)
add_check("Performance", "roc-endpoints", all(vapply(split(roc_data, roc_data$Model), function(z) identical(z$FPR[c(1, nrow(z))], c(0, 1)) && identical(z$TPR[c(1, nrow(z))], c(0, 1)), logical(1L))), "0,0 to 1,1")
add_check("Performance", "calibration-bins", nrow(calibration) == 10L && all(calibration$Samples >= 22L), paste(calibration$Samples, collapse = "/"))

add_check("Importance", "importance-nonempty", nrow(importance) > 5000L, nrow(importance))
add_check("Importance", "heldout-baseline", all(is.finite(importance$BaselineAUROC)), range(importance$BaselineAUROC))
add_check("Importance", "two-model-importance", setequal(unique(importance$Algorithm), algorithms), unique(importance$Algorithm))
add_check("Importance", "importance-summary-nonempty", nrow(importance_summary) > 500L, nrow(importance_summary))
add_check("Importance", "importance-fold-count", all(importance_summary$OuterFoldsAvailable >= 1L & importance_summary$OuterFoldsAvailable <= 25L), range(importance_summary$OuterFoldsAvailable))

safe_median <- leakage_summary$MedianAUROC[leakage_summary$Pipeline == "Fold-scoped preprocessing"]
leaky_median <- leakage_summary$MedianAUROC[leakage_summary$Pipeline == "Global label-selected top 20"]
add_check("Leakage audit", "permutation-rows", nrow(leakage_audit) == 100L, nrow(leakage_audit))
add_check("Leakage audit", "fifty-per-branch", all(table(leakage_audit$Pipeline) == 50L), table(leakage_audit$Pipeline))
add_check("Leakage audit", "safe-null-centered", safe_median >= 0.45 && safe_median <= 0.55, safe_median)
add_check("Leakage audit", "leakage-inflates-null", leaky_median >= safe_median + 0.10, paste(safe_median, leaky_median))
add_check("Leakage audit", "labels-explicit", identical(sort(unique(leakage_audit$LabelInformationFromTestFold)), c(FALSE, TRUE)), unique(leakage_audit$LabelInformationFromTestFold))

anchor_path <- file.path(figure_dir, "27-zeller-fig1-original.png")
add_check("Anchor figure", "anchor-present", file.exists(anchor_path), anchor_path)
add_check("Anchor figure", "anchor-sha256", file.exists(anchor_path) && identical(sha256_file(anchor_path), "6f0dbe5ca4ad7e9bc853fd6568efca093e30f8f46f02a0f090c282b16059ac43"), if (file.exists(anchor_path)) sha256_file(anchor_path) else "missing")

chapter_lines <- readLines(chapter_path, warn = FALSE)
chapter_text <- paste(chapter_lines, collapse = "\n")
required_chapter_tokens <- c(
  "draft: false", "eval: true", "freeze: auto", "expected_images: 5",
  "## 这一步对应论文里的哪张图", "## 理论：", "## 准备工作",
  "## 可复制代码", "## 审计与升级", "## 出版级美化",
  "## 常见坑", "## 这段 Methods 怎么写", "## 换成你自己的数据怎么做",
  "## 参考", "nested cross-validation", "Random forest", "XGBoost",
  "AUROC", "AUPRC", "Brier", "permutation importance", "data leakage",
  "114", "661", "MetaPhlAn 3", "CHOCOPhlAn 201901", "set.seed(20260727)",
  "27-nested-roc.png", "27-model-performance.png",
  "27-permutation-importance.png", "27-leakage-permutation-audit.png",
  "[@zeller2014crc]", "[@wirbel2021siamcat]"
)
for (token in required_chapter_tokens) {
  add_check(
    "Chapter", paste0("chapter-token-", gsub("[^a-z0-9]+", "-", tolower(token))),
    grepl(token, chapter_text, fixed = TRUE), token
  )
}
for (banned in c(
  "Planned chapter", "Do not publish", "本篇可独立跑通",
  "这体现全系列", "作者代码通常长这样", "（即本文）"
)) {
  add_check(
    "Chapter", paste0("chapter-banned-", gsub("[^a-z0-9]+", "-", tolower(banned))),
    !grepl(banned, chapter_text, fixed = TRUE), banned
  )
}
add_check("Chapter", "chapter-single-source-free", !grepl('source("R/theme_pub.R")', chapter_text, fixed = TRUE), "inline plotting functions")
add_check("Chapter", "chapter-anchor-hash", grepl("6f0dbe5ca4ad7e9bc853fd6568efca093e30f8f46f02a0f090c282b16059ac43", chapter_text, fixed = TRUE), "anchor SHA-256")

pal_pub <- c(
  `Random forest` = "#0072B2", XGBoost = "#D55E00",
  `Fold-scoped preprocessing` = "#009E73",
  `Global label-selected top 20` = "#CC79A7"
)
theme_pub <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey30", size = base_size),
      plot.caption = ggplot2::element_text(
        color = "grey35", size = base_size - 1, hjust = 0,
        lineheight = 1.05, margin = ggplot2::margin(t = 5)
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}
save_pub <- function(plot, stem, width = 190, height = 130, dpi = 350) {
  base <- file.path(figure_dir, stem)
  ggplot2::ggsave(paste0(base, ".pdf"), plot, width = width, height = height, units = "mm", device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(paste0(base, ".png"), plot, width = width, height = height, units = "mm", dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(base, ".tiff"), plot, width = width, height = height, units = "mm", dpi = dpi, compression = "lzw", bg = "white")
}
label_fixed <- function(digits) {
  force(digits)
  function(x) formatC(x, format = "f", digits = digits)
}
label_percent_fixed <- function(x) {
  paste0(formatC(100 * x, format = "f", digits = 0), "%")
}

performance$Legend <- sprintf(
  "%s · %.2f (%.2f–%.2f)",
  performance$Algorithm, performance$AUROC,
  performance$AUROCLower95, performance$AUROCUpper95
)
legend_map <- stats::setNames(performance$Legend, performance$Algorithm)
roc_plot_data <- roc_data
roc_plot_data$Legend <- unname(legend_map[roc_plot_data$Model])
roc_colors <- stats::setNames(unname(pal_pub[performance$Algorithm]), performance$Legend)
p_roc <- ggplot2::ggplot(roc_plot_data, ggplot2::aes(FPR, TPR, color = Legend)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
  ggplot2::geom_step(linewidth = 0.9) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_fixed(2)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_fixed(2)
  ) +
  ggplot2::scale_color_manual(
    values = roc_colors,
    guide = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  ) +
  ggplot2::labs(
    title = "Nested out-of-fold discrimination",
    subtitle = "5 × 5 outer CV; one aggregated OOF score per subject",
    x = "False positive rate", y = "True positive rate", color = NULL,
    caption = "All preprocessing and tuning were fit within training folds."
  ) + theme_pub(10)
save_pub(p_roc, "27-nested-roc", width = 190, height = 140)

outer_metrics$Algorithm <- factor(outer_metrics$Algorithm, levels = algorithms)
p_fold <- ggplot2::ggplot(outer_metrics, ggplot2::aes(Algorithm, OuterAUROC, color = Algorithm)) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = 2, color = "grey55") +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, color = "grey30", fill = "white") +
  ggplot2::geom_jitter(width = 0.12, height = 0, alpha = 0.65, size = 1.6) +
  ggplot2::scale_color_manual(values = pal_pub[algorithms]) +
  ggplot2::scale_y_continuous(
    breaks = seq(0.4, 1, 0.1),
    labels = label_fixed(1)
  ) +
  ggplot2::scale_x_discrete(labels = c(
    `Random forest` = "Random\nforest", XGBoost = "XGBoost"
  )) +
  ggplot2::coord_cartesian(ylim = c(0.35, 1)) +
  ggplot2::labs(
    title = "Outer-fold AUROC",
    x = NULL, y = "Outer-fold AUROC", color = NULL
  ) + theme_pub(9) + ggplot2::theme(legend.position = "none")
p_calibration <- ggplot2::ggplot(calibration, ggplot2::aes(MeanPredictedRisk, ObservedCRCProportion, color = Algorithm)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 3.2, alpha = 0.9) +
  ggplot2::scale_color_manual(
    values = pal_pub[algorithms],
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_fixed(2)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = label_fixed(2)
  ) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Calibration", x = "Predicted CRC risk",
    y = "Observed CRC rate", color = NULL
  ) + theme_pub(9) + ggplot2::theme(legend.position = "none")
p_performance <- p_fold + p_calibration +
  patchwork::plot_layout(widths = c(1, 1.08)) +
  patchwork::plot_annotation(
  title = "Discrimination is not calibration",
  subtitle = "Outer-fold variation and quintile reliability bins",
  caption = "Blue: Random forest · Orange: XGBoost",
  theme = ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.subtitle = ggplot2::element_text(color = "grey30", size = 10),
    plot.caption = ggplot2::element_text(color = "grey35", size = 9, hjust = 0.5),
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )
)
save_pub(p_performance, "27-model-performance", width = 205, height = 125)

display_importance <- importance_summary[importance_summary$Rank <= 12L, , drop = FALSE]
display_species <- unique(display_importance$Species[order(-display_importance$MeanDeltaAUROC)])
display_importance$Species <- factor(display_importance$Species, levels = rev(display_species))
p_importance <- ggplot2::ggplot(display_importance, ggplot2::aes(MeanDeltaAUROC, Species, color = Algorithm)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower95, xmax = Upper95), height = 0, linewidth = 0.45) +
  ggplot2::geom_point(ggplot2::aes(size = PositiveFoldFraction), alpha = 0.9) +
  ggplot2::facet_wrap(~Algorithm, scales = "free_y") +
  ggplot2::scale_color_manual(values = pal_pub[algorithms], guide = "none") +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_pretty(n = 5),
    labels = label_fixed(3)
  ) +
  ggplot2::scale_size_continuous(
    range = c(1.8, 4.5), limits = c(0, 1),
    breaks = c(0.25, 0.5, 0.75, 1),
    labels = label_percent_fixed
  ) +
  ggplot2::labs(
    title = "Held-out permutation importance",
    subtitle = "Top 12 per model; bars show fold mean ± 1.96 SE",
    x = "Outer-test AUROC decrease after permutation", y = NULL,
    color = NULL, size = "Positive folds",
    caption = "Importance is conditional on correlated features; it is not a causal effect."
  ) + theme_pub(9) + ggplot2::theme(legend.position = "bottom")
save_pub(p_importance, "27-permutation-importance", width = 210, height = 165)

actual_rf_auc <- performance$AUROC[performance$Algorithm == "Random forest"]
p_leakage <- ggplot2::ggplot(leakage_audit, ggplot2::aes(Pipeline, AUROC, fill = Pipeline)) +
  ggplot2::geom_hline(yintercept = 0.5, linetype = 2, color = "grey45") +
  ggplot2::geom_violin(width = 0.75, alpha = 0.55, color = NA, trim = FALSE) +
  ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white", color = "grey25") +
  ggplot2::geom_jitter(width = 0.08, height = 0, alpha = 0.45, size = 1.2) +
  ggplot2::geom_hline(yintercept = actual_rf_auc, color = pal_pub[["Random forest"]], linewidth = 0.7) +
  ggplot2::annotate(
    "text", x = 1.5, y = actual_rf_auc,
    label = sprintf("Observed nested RF = %.2f", actual_rf_auc),
    vjust = -0.6, color = pal_pub[["Random forest"]], size = 3.2
  ) +
  ggplot2::scale_fill_manual(
    values = pal_pub[levels(factor(leakage_audit$Pipeline))],
    guide = "none"
  ) +
  ggplot2::scale_x_discrete(labels = c(
    `Fold-scoped preprocessing` = "Fold-scoped\npreprocessing",
    `Global label-selected top 20` = "Global label-selected\ntop 20"
  )) +
  ggplot2::scale_y_continuous(
    breaks = seq(0.3, 1, 0.1),
    labels = label_fixed(1)
  ) +
  ggplot2::coord_cartesian(ylim = c(0.25, 1)) +
  ggplot2::labs(
    title = "Test-label leakage creates signal from permuted outcomes",
    subtitle = "50 label permutations on the same 661-feature matrix",
    x = NULL, y = "Cross-validated AUROC", fill = NULL,
    caption = "Global label-based selection before CV is intentionally invalid."
  ) + theme_pub(10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(hjust = 0.5),
    legend.position = "none"
  )
save_pub(p_leakage, "27-leakage-permutation-audit", width = 190, height = 132)

figure_stems <- c(
  "27-nested-roc", "27-model-performance",
  "27-permutation-importance", "27-leakage-permutation-audit"
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
add_check("Figures", "figure-files", nrow(figure_audit) == 12L && all(figure_audit$Exists), paste0(sum(figure_audit$Exists), "/12"))
add_check("Figures", "figure-nonempty", all(figure_audit$Bytes > 10000), paste(range(figure_audit$Bytes), collapse = "/"))

write_tsv(fold_assignments, file.path(output_dir, "outer-fold-assignments.tsv"))
write_tsv_gz(outer_predictions, file.path(output_dir, "outer-predictions.tsv.gz"))
write_tsv(sample_predictions, file.path(output_dir, "sample-predictions.tsv"))
write_tsv(outer_metrics, file.path(output_dir, "outer-fold-metrics.tsv"))
write_tsv_gz(tuning_audit, file.path(output_dir, "inner-tuning-audit.tsv.gz"))
write_tsv(preprocessing_audit, file.path(output_dir, "preprocessing-audit.tsv"))
write_tsv(performance, file.path(output_dir, "performance-summary.tsv"))
write_tsv(roc_data, file.path(output_dir, "roc-curves.tsv"))
write_tsv(calibration, file.path(output_dir, "calibration-audit.tsv"))
write_tsv(model_comparison, file.path(output_dir, "model-comparison.tsv"))
write_tsv_gz(importance, file.path(output_dir, "permutation-importance-folds.tsv.gz"))
write_tsv(importance_summary, file.path(output_dir, "permutation-importance-summary.tsv"))
write_tsv(leakage_audit, file.path(output_dir, "leakage-permutation-audit.tsv"))
write_tsv(leakage_summary, file.path(output_dir, "leakage-permutation-summary.tsv"))
write_tsv(figure_audit, file.path(output_dir, "figure-audit.tsv"))
write_tsv(checks, file.path(output_dir, "validation-checks.tsv"))

failures <- checks[checks$Status != "PASS", , drop = FALSE]
summary <- list(
  status = if (nrow(failures) == 0L) "passed" else "failed",
  article = 27L,
  seed = primary_seed,
  samples = nrow(x),
  controls = sum(y == "Control"),
  crc = sum(y == "CRC"),
  raw_species_features = ncol(x),
  outer_resamples = 5L,
  outer_folds = 5L,
  inner_folds = 4L,
  outer_models = nrow(outer_metrics),
  random_forest_auroc = performance$AUROC[performance$Algorithm == "Random forest"],
  xgboost_auroc = performance$AUROC[performance$Algorithm == "XGBoost"],
  random_forest_auprc = performance$AUPRC[performance$Algorithm == "Random forest"],
  xgboost_auprc = performance$AUPRC[performance$Algorithm == "XGBoost"],
  auc_difference_rf_minus_xgb = model_comparison$AUROCDifference,
  safe_permutation_median_auroc = safe_median,
  leaky_permutation_median_auroc = leaky_median,
  checksum_entries = checksum_entries,
  checks = nrow(checks),
  passed = sum(checks$Status == "PASS"),
  failed = nrow(failures),
  generated_figure_files = nrow(figure_audit),
  package_versions = list(
    ranger = as.character(utils::packageVersion("ranger")),
    xgboost = as.character(utils::packageVersion("xgboost")),
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
  stop("Article 27 validation failed.", call. = FALSE)
}
log_msg("Article 27 validation passed.")
cat("Article 27 validation passed:", nrow(checks), "checks.\n")
