#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(permute)
})

options(stringsAsFactors = FALSE, scipen = 999)

parse_args <- function(args) {
  result <- list(input_dir = ".tutorial_runs/article63",
                 output_dir = ".tutorial_runs/article63/global",
                 permutations = 9999L,
                 bootstraps = 2000L,
                 seed = 63001L)
  index <- 1L
  while (index <= length(args)) {
    key <- sub("^--", "", args[[index]])
    if (!key %in% names(result) || index == length(args)) {
      stop("Unknown or incomplete argument: ", args[[index]])
    }
    value <- args[[index + 1L]]
    result[[key]] <- if (key %in% c("permutations", "bootstraps", "seed")) {
      as.integer(value)
    } else {
      value
    }
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

procrustes_r <- function(x, y) {
  fit <- vegan::procrustes(x, y, symmetric = TRUE)
  sqrt(max(0, 1 - fit$ss))
}

bootstrap_ci <- function(x, y, n_boot, seed) {
  set.seed(seed)
  estimates <- numeric(n_boot)
  n <- nrow(x)
  for (iteration in seq_len(n_boot)) {
    sampled <- sample.int(n, n, replace = TRUE)
    estimates[[iteration]] <- procrustes_r(x[sampled, , drop = FALSE],
                                            y[sampled, , drop = FALSE])
  }
  unname(quantile(estimates, c(0.025, 0.975), names = FALSE, type = 7))
}

permutation_control <- function(groups, n_perm, restricted) {
  control <- permute::how(nperm = n_perm)
  if (restricted) permute::setBlocks(control) <- factor(groups)
  control
}

run_global_test <- function(cohort, metadata, microbe_scores, metabolite_scores,
                            microbe_values, metabolite_values, n_perm, n_boot,
                            seed) {
  rows <- list()
  boot_ci <- bootstrap_ci(microbe_scores, metabolite_scores, n_boot, seed + 1L)
  for (restricted in c(FALSE, TRUE)) {
    control <- permutation_control(metadata$Study.Group, n_perm, restricted)
    set.seed(seed + ifelse(restricted, 20L, 10L))
    protest_fit <- vegan::protest(
      microbe_scores, metabolite_scores,
      permutations = control, symmetric = TRUE
    )
    set.seed(seed + ifelse(restricted, 40L, 30L))
    mantel_fit <- vegan::mantel(
      dist(microbe_values), dist(metabolite_values),
      method = "spearman", permutations = control
    )
    rows[[length(rows) + 1L]] <- data.frame(
      Cohort = cohort,
      Restriction = ifelse(restricted, "Within diagnosis", "Unrestricted"),
      Samples = nrow(metadata),
      ProcrustesR = as.numeric(protest_fit$t0),
      ProcrustesP = as.numeric(protest_fit$signif),
      ProcrustesBootLow = boot_ci[[1]],
      ProcrustesBootHigh = boot_ci[[2]],
      MantelRho = as.numeric(mantel_fit$statistic),
      MantelP = as.numeric(mantel_fit$signif),
      Permutations = n_perm,
      Bootstraps = n_boot
    )
  }
  do.call(rbind, rows)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- normalizePath(args$input_dir, mustWork = TRUE)
output_dir <- args$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(file.path(input_dir, "sample-metadata.tsv"),
                       check.names = FALSE)
if (anyDuplicated(metadata$Sample) || anyDuplicated(metadata$Subject)) {
  stop("Global concordance requires one row per independent subject")
}
microbiome <- read_matrix(file.path(input_dir, "microbiome-clr.tsv.gz"))
metabolome <- read_matrix(file.path(input_dir, "metabolome-log1p.tsv.gz"))
if (!identical(metadata$Sample, rownames(microbiome)) ||
    !identical(metadata$Sample, rownames(metabolome))) {
  stop("Metadata and matrices are not identically sample aligned")
}
if (any(!is.finite(microbiome)) || any(!is.finite(metabolome))) {
  stop("Non-finite values in transformed matrices")
}

training <- metadata$Cohort == "PRISM"
validation <- metadata$Cohort == "Validation"
if (sum(training) != 155L || sum(validation) != 65L) {
  stop("Unexpected discovery/external-validation sample counts")
}

# Both ordination bases and metabolite scaling are estimated in PRISM only.
# The independent Validation cohort is projected without refitting either basis.
set.seed(args$seed)
microbe_pca <- prcomp(microbiome[training, , drop = FALSE],
                     center = TRUE, scale. = FALSE, rank. = 10L)
metabolite_pca <- prcomp(metabolome[training, , drop = FALSE],
                        center = TRUE, scale. = TRUE, rank. = 10L)

microbe_scores <- rbind(
  microbe_pca$x[, seq_len(10L), drop = FALSE],
  predict(microbe_pca, microbiome[validation, , drop = FALSE])[, seq_len(10L), drop = FALSE]
)
metabolite_scores <- rbind(
  metabolite_pca$x[, seq_len(10L), drop = FALSE],
  predict(metabolite_pca, metabolome[validation, , drop = FALSE])[, seq_len(10L), drop = FALSE]
)
# Restore source order after rbind(PRISM, Validation).
microbe_scores <- microbe_scores[metadata$Sample, , drop = FALSE]
metabolite_scores <- metabolite_scores[metadata$Sample, , drop = FALSE]

metabolite_scaled <- sweep(metabolome, 2L, metabolite_pca$center, "-")
metabolite_scaled <- sweep(metabolite_scaled, 2L, metabolite_pca$scale, "/")
microbiome_centered <- sweep(microbiome, 2L, microbe_pca$center, "-")

results <- list()
score_rows <- list()
for (cohort_index in seq_along(c("PRISM", "Validation"))) {
  cohort <- c("PRISM", "Validation")[[cohort_index]]
  selected <- metadata$Cohort == cohort
  current_metadata <- metadata[selected, , drop = FALSE]
  x_scores <- microbe_scores[selected, , drop = FALSE]
  y_scores <- metabolite_scores[selected, , drop = FALSE]
  fit <- vegan::procrustes(x_scores, y_scores, symmetric = TRUE)
  current_scores <- data.frame(
    Sample = current_metadata$Sample,
    Cohort = cohort,
    Diagnosis = current_metadata$Study.Group,
    MicrobiomeAxis1 = fit$X[, 1],
    MicrobiomeAxis2 = fit$X[, 2],
    MetabolomeAxis1 = fit$Yrot[, 1],
    MetabolomeAxis2 = fit$Yrot[, 2],
    Residual = residuals(fit)
  )
  score_rows[[cohort_index]] <- current_scores
  results[[cohort_index]] <- run_global_test(
    cohort = cohort,
    metadata = current_metadata,
    microbe_scores = x_scores,
    metabolite_scores = y_scores,
    microbe_values = microbiome_centered[selected, , drop = FALSE],
    metabolite_values = metabolite_scaled[selected, , drop = FALSE],
    n_perm = args$permutations,
    n_boot = args$bootstraps,
    seed = args$seed + cohort_index * 1000L
  )
}

global <- do.call(rbind, results)
scores <- do.call(rbind, score_rows)
variance <- rbind(
  data.frame(Modality = "Microbiome CLR",
             Axis = seq_len(10L),
             VarianceExplained = microbe_pca$sdev[seq_len(10L)]^2 /
               sum(microbe_pca$sdev^2)),
  data.frame(Modality = "Metabolome log1p autoscaled",
             Axis = seq_len(10L),
             VarianceExplained = metabolite_pca$sdev[seq_len(10L)]^2 /
               sum(metabolite_pca$sdev^2))
)

write_tsv(global, file.path(output_dir, "global-concordance.tsv"))
write_tsv(scores, file.path(output_dir, "procrustes-scores.tsv"))
write_tsv(variance, file.path(output_dir, "ordination-variance.tsv"))
saveRDS(list(microbiome = microbe_pca, metabolome = metabolite_pca),
        file.path(output_dir, "prism-ordination-models.rds"), version = 3)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session-info.txt"), useBytes = TRUE)

print(global, row.names = FALSE)
