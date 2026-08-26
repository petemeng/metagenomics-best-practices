#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mixOmics)
  library(BiocParallel)
})

options(stringsAsFactors = FALSE, scipen = 999)

parse_args <- function(args) {
  result <- list(input_dir = ".tutorial_runs/article63",
                 output_dir = ".tutorial_runs/article63/diablo",
                 folds = 5L, repeats = 5L, bootstraps = 100L,
                 permutations = 200L, metric_bootstraps = 2000L,
                 seed = 63001L, force_tuning = FALSE)
  index <- 1L
  while (index <= length(args)) {
    key <- sub("^--", "", args[[index]])
    if (!key %in% names(result) || index == length(args)) {
      stop("Unknown or incomplete argument: ", args[[index]])
    }
    value <- args[[index + 1L]]
    if (key %in% c("folds", "repeats", "bootstraps", "permutations",
                   "metric_bootstraps", "seed")) value <- as.integer(value)
    if (key == "force_tuning") value <- tolower(value) %in% c("1", "true", "yes")
    result[[key]] <- value
    index <- index + 2L
  }
  result
}

read_matrix <- function(path) {
  frame <- read.delim(path, check.names = FALSE)
  if (!"Sample" %in% names(frame)) stop("Missing Sample column in ", path)
  sample_id <- frame$Sample
  values <- as.matrix(frame[, setdiff(names(frame), "Sample"), drop = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- sample_id
  values
}

write_tsv <- function(frame, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "", eol = "\n")
}

prediction_class <- function(prediction, component) {
  vote <- prediction$WeightedVote$centroids.dist
  column <- paste0("comp", component)
  if (is.null(dim(vote)) || !column %in% colnames(vote)) {
    stop("Weighted-vote prediction does not contain ", column)
  }
  as.character(vote[, column])
}

classification_metrics <- function(truth, predicted, levels) {
  truth <- factor(truth, levels = levels)
  predicted <- factor(predicted, levels = levels)
  confusion <- table(Truth = truth, Predicted = predicted)
  per_class <- lapply(levels, function(level) {
    tp <- confusion[level, level]
    fn <- sum(confusion[level, ]) - tp
    fp <- sum(confusion[, level]) - tp
    tn <- sum(confusion) - tp - fn - fp
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    f1 <- if (is.finite(precision + recall) && (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else 0
    data.frame(Class = level, Recall = recall, Precision = precision,
               Specificity = specificity, F1 = f1)
  })
  per_class <- do.call(rbind, per_class)
  overall <- data.frame(
    Accuracy = sum(diag(confusion)) / sum(confusion),
    BalancedAccuracy = mean(per_class$Recall),
    MacroF1 = mean(per_class$F1)
  )
  list(overall = overall, per_class = per_class, confusion = confusion)
}

metric_ci <- function(truth, predicted, levels, n_boot, seed) {
  set.seed(seed)
  indices <- split(seq_along(truth), factor(truth, levels = levels))
  draws <- matrix(NA_real_, nrow = n_boot, ncol = 3L,
                  dimnames = list(NULL, c("Accuracy", "BalancedAccuracy", "MacroF1")))
  for (iteration in seq_len(n_boot)) {
    sampled <- unlist(lapply(indices, function(index) {
      sample(index, length(index), replace = TRUE)
    }), use.names = FALSE)
    current <- classification_metrics(truth[sampled], predicted[sampled], levels)$overall
    draws[iteration, ] <- unlist(current[1, ], use.names = FALSE)
  }
  data.frame(
    Metric = colnames(draws),
    Low = apply(draws, 2L, quantile, probs = 0.025, names = FALSE, type = 7),
    High = apply(draws, 2L, quantile, probs = 0.975, names = FALSE, type = 7)
  )
}

sample_stratified <- function(y) {
  unlist(lapply(levels(y), function(level) {
    index <- which(y == level)
    sample(index, length(index), replace = TRUE)
  }), use.names = FALSE)
}

fit_diablo <- function(x, y, keep_x, design, seed) {
  set.seed(seed)
  invisible(capture.output(
    fit <- mixOmics::block.splsda(
      X = x, Y = y, ncomp = 2L, keepX = keep_x,
      design = design, scale = TRUE, near.zero.var = FALSE,
      max.iter = 100L
    )
  ))
  fit
}

selected_features <- function(fit, blocks, replicate_id) {
  rows <- list()
  position <- 1L
  for (block in blocks) {
    for (component in seq_len(2L)) {
      selected <- mixOmics::selectVar(fit, block = block, comp = component)[[block]]
      rows[[position]] <- data.frame(
        Replicate = replicate_id,
        Block = block,
        Component = component,
        FeatureID = selected$name,
        Loading = selected$value$value.var
      )
      position <- position + 1L
    }
  }
  do.call(rbind, rows)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- normalizePath(args$input_dir, mustWork = TRUE)
output_dir <- args$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(file.path(input_dir, "sample-metadata.tsv"),
                       check.names = FALSE)
features <- read.delim(file.path(input_dir, "feature-audit.tsv"),
                       check.names = FALSE)
microbiome <- read_matrix(file.path(input_dir, "microbiome-clr.tsv.gz"))
metabolome <- read_matrix(file.path(input_dir, "metabolome-log1p.tsv.gz"))
if (!identical(metadata$Sample, rownames(microbiome)) ||
    !identical(metadata$Sample, rownames(metabolome))) {
  stop("Metadata and matrices are not identically sample aligned")
}
if (anyDuplicated(metadata$Subject)) stop("DIABLO requires independent subjects")

class_levels <- c("Control", "CD", "UC")
training <- metadata$Cohort == "PRISM"
validation <- metadata$Cohort == "Validation"
x_train <- list(
  microbiome = microbiome[training, , drop = FALSE],
  metabolome = metabolome[training, , drop = FALSE]
)
x_validation <- list(
  microbiome = microbiome[validation, , drop = FALSE],
  metabolome = metabolome[validation, , drop = FALSE]
)
y_train <- factor(metadata$Study.Group[training], levels = class_levels)
y_validation <- factor(metadata$Study.Group[validation], levels = class_levels)
if (anyNA(y_train) || anyNA(y_validation)) stop("Unexpected diagnostic group")
if (nrow(x_train$microbiome) != 155L || nrow(x_validation$microbiome) != 65L) {
  stop("Unexpected discovery/external-validation sample counts")
}

blocks <- names(x_train)
design <- matrix(c(0, 0.1, 0.1, 0), nrow = 2L, byrow = TRUE,
                 dimnames = list(blocks, blocks))
test_keep_x <- list(microbiome = c(5L, 10L, 20L),
                    metabolome = c(5L, 10L, 20L))
tuning_path <- file.path(output_dir, "tuning-object.rds")
if (file.exists(tuning_path) && !args$force_tuning) {
  tuning <- readRDS(tuning_path)
} else {
  set.seed(args$seed + 100L)
  tuning <- mixOmics::tune.block.splsda(
    X = x_train, Y = y_train, ncomp = 2L,
    test.keepX = test_keep_x,
    validation = "Mfold", folds = args$folds,
    dist = "centroids.dist", measure = "BER", weighted = TRUE,
    nrepeat = args$repeats, design = design, scale = TRUE,
    near.zero.var = FALSE, progressBar = FALSE,
    BPPARAM = BiocParallel::SerialParam(), seed = args$seed + 100L
  )
  saveRDS(tuning, tuning_path, version = 3)
}

keep_x <- lapply(blocks, function(block) {
  chosen <- as.integer(tuning$choice.keepX[[block]])
  if (!length(chosen)) stop("Tuning did not return keepX for ", block)
  if (length(chosen) < 2L) chosen <- c(chosen, rep(tail(chosen, 1L), 2L - length(chosen)))
  chosen[seq_len(2L)]
})
names(keep_x) <- blocks

final_fit <- fit_diablo(x_train, y_train, keep_x, design, args$seed + 200L)
saveRDS(final_fit, file.path(output_dir, "final-model.rds"), version = 3)
external_prediction <- predict(final_fit, newdata = x_validation)
predicted <- prediction_class(external_prediction, 2L)
metrics <- classification_metrics(y_validation, predicted, class_levels)
confidence <- metric_ci(as.character(y_validation), predicted, class_levels,
                        args$metric_bootstraps, args$seed + 300L)
overall <- merge(
  data.frame(Metric = names(metrics$overall),
             Estimate = as.numeric(metrics$overall[1, ])),
  confidence, by = "Metric", sort = FALSE
)

tuning_summary <- do.call(rbind, lapply(blocks, function(block) {
  data.frame(Block = block, Component = seq_len(2L), KeepX = keep_x[[block]])
}))
tuning_summary$DesignWeight <- 0.1
tuning_summary$Folds <- args$folds
tuning_summary$Repeats <- args$repeats
tuning_summary$Measure <- "Balanced error rate"
tuning_summary$Distance <- "Centroids"
write_tsv(tuning_summary, file.path(output_dir, "tuning-summary.tsv"))

confusion <- as.data.frame(metrics$confusion)
names(confusion) <- c("Truth", "Predicted", "Samples")
predictions <- data.frame(
  Sample = metadata$Sample[validation],
  Truth = as.character(y_validation),
  Predicted = predicted,
  Correct = as.character(y_validation) == predicted
)
write_tsv(overall, file.path(output_dir, "external-metrics.tsv"))
write_tsv(metrics$per_class, file.path(output_dir, "external-class-metrics.tsv"))
write_tsv(confusion, file.path(output_dir, "external-confusion.tsv"))
write_tsv(predictions, file.path(output_dir, "external-predictions.tsv"))

latent_rows <- list()
for (cohort in c("PRISM", "Validation")) {
  if (cohort == "PRISM") {
    variates <- final_fit$variates[blocks]
    current_metadata <- metadata[training, , drop = FALSE]
  } else {
    variates <- external_prediction$variates[blocks]
    current_metadata <- metadata[validation, , drop = FALSE]
  }
  for (component in seq_len(2L)) {
    latent_rows[[length(latent_rows) + 1L]] <- data.frame(
      Sample = current_metadata$Sample,
      Cohort = cohort,
      Diagnosis = current_metadata$Study.Group,
      Component = component,
      MicrobiomeScore = variates$microbiome[, component],
      MetabolomeScore = variates$metabolome[, component]
    )
  }
}
latent <- do.call(rbind, latent_rows)
latent_correlations <- do.call(rbind, lapply(split(latent, list(latent$Cohort,
                                                                latent$Component)),
                                             function(frame) {
  data.frame(Cohort = frame$Cohort[[1]], Component = frame$Component[[1]],
             SpearmanRho = cor(frame$MicrobiomeScore, frame$MetabolomeScore,
                               method = "spearman"), Samples = nrow(frame))
}))
write_tsv(latent, file.path(output_dir, "latent-scores.tsv"))
write_tsv(latent_correlations, file.path(output_dir, "latent-correlations.tsv"))

final_selection <- selected_features(final_fit, blocks, "Final")
feature_labels <- features[, c("FeatureID", "DisplayName", "Modality",
                               "Phylum", "HMDB", "ChemicalClass")]
final_selection <- merge(final_selection, feature_labels, by = "FeatureID",
                         all.x = TRUE, sort = FALSE)
write_tsv(final_selection, file.path(output_dir, "final-selected-features.tsv"))

# Stratified PRISM bootstrap: quantify how often each feature survives a refit.
bootstrap_checkpoint <- file.path(output_dir, "bootstrap-checkpoint.rds")
if (file.exists(bootstrap_checkpoint)) {
  bootstrap_state <- readRDS(bootstrap_checkpoint)
} else {
  bootstrap_state <- list(completed = 0L, selected = list())
}
if (bootstrap_state$completed < args$bootstraps) {
  for (iteration in seq.int(bootstrap_state$completed + 1L, args$bootstraps)) {
    set.seed(args$seed + 10000L + iteration)
    sampled <- sample_stratified(y_train)
    current_x <- lapply(x_train, function(block) block[sampled, , drop = FALSE])
    current_y <- factor(y_train[sampled], levels = class_levels)
    unique_names <- sprintf("Bootstrap%03d_%03d", iteration, seq_along(sampled))
    current_x <- lapply(current_x, function(block) {
      rownames(block) <- unique_names
      block
    })
    names(current_y) <- unique_names
    current_fit <- fit_diablo(current_x, current_y, keep_x, design,
                              args$seed + 20000L + iteration)
    bootstrap_state$selected[[iteration]] <- selected_features(
      current_fit, blocks, iteration
    )[, c("Replicate", "Block", "Component", "FeatureID")]
    bootstrap_state$completed <- iteration
    if (iteration %% 5L == 0L || iteration == args$bootstraps) {
      saveRDS(bootstrap_state, bootstrap_checkpoint, version = 3)
    }
  }
}
bootstrap_selected <- do.call(rbind, bootstrap_state$selected[seq_len(args$bootstraps)])
stability <- aggregate(Replicate ~ Block + Component + FeatureID,
                       bootstrap_selected, function(value) length(unique(value)))
names(stability)[names(stability) == "Replicate"] <- "SelectedBootstraps"
stability$Bootstraps <- args$bootstraps
stability$SelectionFrequency <- stability$SelectedBootstraps / args$bootstraps
stability <- merge(stability, feature_labels, by = "FeatureID", all.x = TRUE,
                   sort = FALSE)
stability <- stability[order(stability$Block, stability$Component,
                             -stability$SelectionFrequency, stability$FeatureID), ]
write_tsv(stability, file.path(output_dir, "bootstrap-feature-stability.tsv"))

# Fixed-complexity label-permutation controls: the external cohort never changes
# model selection, and these null fits are not used to revise the final model.
permutation_checkpoint <- file.path(output_dir, "permutation-checkpoint.rds")
if (file.exists(permutation_checkpoint)) {
  null_accuracy <- readRDS(permutation_checkpoint)
  if (length(null_accuracy) != args$permutations) null_accuracy <- rep(NA_real_, args$permutations)
} else {
  null_accuracy <- rep(NA_real_, args$permutations)
}
for (iteration in which(is.na(null_accuracy))) {
  set.seed(args$seed + 30000L + iteration)
  permuted_y <- factor(sample(as.character(y_train), replace = FALSE),
                       levels = class_levels)
  null_fit <- fit_diablo(x_train, permuted_y, keep_x, design,
                         args$seed + 40000L + iteration)
  null_prediction <- prediction_class(predict(null_fit, newdata = x_validation), 2L)
  null_accuracy[[iteration]] <- classification_metrics(
    y_validation, null_prediction, class_levels
  )$overall$BalancedAccuracy
  if (iteration %% 5L == 0L || iteration == args$permutations) {
    saveRDS(null_accuracy, permutation_checkpoint, version = 3)
  }
}
observed_ba <- metrics$overall$BalancedAccuracy
null_table <- data.frame(Iteration = seq_len(args$permutations),
                         BalancedAccuracy = null_accuracy)
null_summary <- data.frame(
  ObservedBalancedAccuracy = observed_ba,
  NullMean = mean(null_accuracy),
  NullSD = sd(null_accuracy),
  NullQ025 = unname(quantile(null_accuracy, 0.025)),
  NullQ975 = unname(quantile(null_accuracy, 0.975)),
  EmpiricalP = (1 + sum(null_accuracy >= observed_ba)) / (1 + args$permutations),
  Permutations = args$permutations
)
write_tsv(null_table, file.path(output_dir, "label-permutation-null.tsv"))
write_tsv(null_summary, file.path(output_dir, "label-permutation-summary.tsv"))

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session-info.txt"), useBytes = TRUE)
print(tuning_summary, row.names = FALSE)
print(overall, row.names = FALSE)
print(null_summary, row.names = FALSE)
