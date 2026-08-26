#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999, digits = 17)
Sys.setenv(TZ = "UTC")

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
required <- c("result-dir", "figure-dir")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c("ggplot2", "patchwork", "scales", "jsonlite", "digest")
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

result_dir <- normalizePath(args[["result-dir"]], mustWork = TRUE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

font_cache_dir <- file.path(tempdir(), "article27-figure-fontconfig-cache")
dir.create(font_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

read_tsv <- function(name) {
  path <- file.path(result_dir, name)
  if (!file.exists(path)) stop("Missing frozen result: ", path, call. = FALSE)
  utils::read.delim(
    path, check.names = FALSE, quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

performance <- read_tsv("performance-summary.tsv")
roc_data <- read_tsv("roc-curves.tsv")
outer_metrics <- read_tsv("outer-fold-metrics.tsv")
calibration <- read_tsv("calibration-audit.tsv")
importance_summary <- read_tsv("permutation-importance-summary.tsv")
leakage_audit <- read_tsv("leakage-permutation-audit.tsv")

algorithms <- c("Random forest", "XGBoost")
stopifnot(
  setequal(performance$Algorithm, algorithms),
  setequal(unique(roc_data$Model), algorithms),
  setequal(unique(outer_metrics$Algorithm), algorithms),
  setequal(unique(calibration$Algorithm), algorithms),
  setequal(unique(importance_summary$Algorithm), algorithms),
  nrow(leakage_audit) == 100L
)

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
      plot.title = ggplot2::element_text(
        face = "bold", size = base_size + 2, hjust = 0
      ),
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
  ggplot2::ggsave(
    paste0(base, ".pdf"), plot, width = width, height = height,
    units = "mm", device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".png"), plot, width = width, height = height,
    units = "mm", dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".tiff"), plot, width = width, height = height,
    units = "mm", dpi = dpi, compression = "lzw", bg = "white"
  )
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
roc_colors <- stats::setNames(
  unname(pal_pub[performance$Algorithm]), performance$Legend
)
p_roc <- ggplot2::ggplot(
  roc_plot_data, ggplot2::aes(FPR, TPR, color = Legend)
) +
  ggplot2::geom_abline(
    slope = 1, intercept = 0, linetype = 2, color = "grey55"
  ) +
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
  ) +
  theme_pub(10)
save_pub(p_roc, "27-nested-roc", width = 190, height = 140)

outer_metrics$Algorithm <- factor(outer_metrics$Algorithm, levels = algorithms)
p_fold <- ggplot2::ggplot(
  outer_metrics, ggplot2::aes(Algorithm, OuterAUROC, color = Algorithm)
) +
  ggplot2::geom_hline(
    yintercept = 0.5, linetype = 2, color = "grey55"
  ) +
  ggplot2::geom_boxplot(
    width = 0.55, outlier.shape = NA, color = "grey30", fill = "white"
  ) +
  ggplot2::geom_jitter(
    width = 0.12, height = 0, alpha = 0.65, size = 1.6
  ) +
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
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "none")

p_calibration <- ggplot2::ggplot(
  calibration,
  ggplot2::aes(MeanPredictedRisk, ObservedCRCProportion, color = Algorithm)
) +
  ggplot2::geom_abline(
    slope = 1, intercept = 0, linetype = 2, color = "grey55"
  ) +
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
    title = "Calibration",
    x = "Predicted CRC risk",
    y = "Observed CRC rate", color = NULL
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "none")

p_performance <- p_fold + p_calibration +
  patchwork::plot_layout(widths = c(1, 1.08)) +
  patchwork::plot_annotation(
    title = "Discrimination is not calibration",
    subtitle = "Outer-fold variation and quintile reliability bins",
    caption = "Blue: Random forest · Orange: XGBoost",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(color = "grey30", size = 10),
      plot.caption = ggplot2::element_text(
        color = "grey35", size = 9, hjust = 0.5
      ),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
  )
save_pub(p_performance, "27-model-performance", width = 205, height = 125)

display_importance <- importance_summary[
  importance_summary$Rank <= 12L, , drop = FALSE
]
display_species <- unique(
  display_importance$Species[order(-display_importance$MeanDeltaAUROC)]
)
display_importance$Species <- factor(
  display_importance$Species, levels = rev(display_species)
)
p_importance <- ggplot2::ggplot(
  display_importance,
  ggplot2::aes(MeanDeltaAUROC, Species, color = Algorithm)
) +
  ggplot2::geom_vline(
    xintercept = 0, linetype = 2, color = "grey55"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = Lower95, xmax = Upper95),
    height = 0, linewidth = 0.45
  ) +
  ggplot2::geom_point(
    ggplot2::aes(size = PositiveFoldFraction), alpha = 0.9
  ) +
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
    caption = paste0(
      "Importance is conditional on correlated features; ",
      "it is not a causal effect."
    )
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "bottom")
save_pub(p_importance, "27-permutation-importance", width = 210, height = 165)

actual_rf_auc <- performance$AUROC[
  performance$Algorithm == "Random forest"
]
p_leakage <- ggplot2::ggplot(
  leakage_audit, ggplot2::aes(Pipeline, AUROC, fill = Pipeline)
) +
  ggplot2::geom_hline(
    yintercept = 0.5, linetype = 2, color = "grey45"
  ) +
  ggplot2::geom_violin(
    width = 0.75, alpha = 0.55, color = NA, trim = FALSE
  ) +
  ggplot2::geom_boxplot(
    width = 0.22, outlier.shape = NA, fill = "white", color = "grey25"
  ) +
  ggplot2::geom_jitter(
    width = 0.08, height = 0, alpha = 0.45, size = 1.2
  ) +
  ggplot2::geom_hline(
    yintercept = actual_rf_auc, color = pal_pub[["Random forest"]],
    linewidth = 0.7
  ) +
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
  ) +
  theme_pub(10) +
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

if (nrow(figure_audit) != 12L || !all(figure_audit$Exists)) {
  stop("Article 27 figure rendering did not create all 12 files.", call. = FALSE)
}
if (!all(figure_audit$Bytes > 10000)) {
  stop("One or more Article 27 figure files are unexpectedly small.", call. = FALSE)
}

write_tsv(figure_audit, file.path(result_dir, "figure-audit.tsv"))
summary <- list(
  status = "passed",
  article = 27L,
  source = "frozen validated result tables",
  figure_stems = length(figure_stems),
  figure_files = nrow(figure_audit),
  minimum_bytes = min(figure_audit$Bytes),
  generated_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
jsonlite::write_json(
  summary, file.path(result_dir, "figure-render-summary.json"),
  auto_unbox = TRUE, pretty = TRUE
)

cat("Article 27 figures rendered: 4 stems, 12 files.\n")
