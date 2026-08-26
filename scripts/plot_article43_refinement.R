#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(forcats)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop("Missing argument: ", flag)
  args[[index + 1]]
}
input_dir <- value_after("--input-dir")
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260743)

method_colors <- c("DAS Tool" = "#0072B2", "Binette" = "#D55E00")
stage_colors <- c("Input binner" = "#999999", "Refinement" = "#009E73")
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
      axis.text = element_text(color = "black"), legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA),
      plot.title.position = "plot"
    )
}
save_pub <- function(plot, file_base, width = 190, height = 125, dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = "mm", dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = "mm", dpi = dpi, compression = "lzw", bg = "white")
}

comparison <- read_tsv(file.path(input_dir, "input-vs-refinement-summary.tsv"), show_col_types = FALSE) |>
  mutate(
    Label = recode(
      Method,
      `MetaBAT2-MOCK1-only` = "MetaBAT2\nMOCK1 only",
      `MetaBAT2-multisample` = "MetaBAT2\nmultisample",
      `SemiBin2-self-supervised` = "SemiBin2\nself-supervised",
      `VAMB-taxonomy-free` = "VAMB\ntaxonomy-free",
      `TaxVAMB-Kraken2` = "TaxVAMB\nKraken2",
      `DAS Tool` = "DAS Tool", Binette = "Binette"
    ),
    Label = factor(Label, levels = Label)
  )
p_yield <- ggplot(comparison, aes(x = Label, y = MinimumPassBins, fill = Stage)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  geom_text(aes(label = MinimumPassBins), vjust = -0.35, size = 3.3) +
  scale_fill_manual(values = stage_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "Refinement is evaluated against every input partition",
    subtitle = "Minimum pass: CheckM2 completeness >=50%, contamination <10%, and GUNC pass",
    x = NULL, y = "Minimum-pass bins", fill = NULL
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 8.5), legend.position = "bottom")
save_pub(p_yield, file.path(output_dir, "43-refinement-yield"), width = 220, height = 130)

quality <- read_tsv(file.path(input_dir, "refinement-quality-truth-audit.tsv"), show_col_types = FALSE) |>
  mutate(
    GUNCPass = as.logical(GUNCPass),
    ReferenceFreeMinimumPass = as.logical(ReferenceFreeMinimumPass),
    Status = if_else(ReferenceFreeMinimumPass, "Minimum pass", "Excluded")
  )
p_quality <- ggplot(
  quality,
  aes(x = CheckM2Contamination, y = CheckM2Completeness, color = Method,
      shape = GUNCPass, size = BinBp, alpha = Status)
) +
  annotate("rect", xmin = -Inf, xmax = 10, ymin = 50, ymax = Inf,
           fill = "#009E73", alpha = 0.06) +
  geom_vline(xintercept = 10, linetype = 3, color = "#555555", linewidth = 0.45) +
  geom_hline(yintercept = 50, linetype = 3, color = "#555555", linewidth = 0.45) +
  geom_point(stroke = 0.55) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "GUNC fail", `TRUE` = "GUNC pass")) +
  scale_size_continuous(range = c(1.8, 6), labels = label_number(scale = 1e-6, suffix = " MB")) +
  scale_alpha_manual(values = c("Excluded" = 0.45, "Minimum pass" = 0.88), guide = "none") +
  coord_cartesian(xlim = c(0, max(15, max(quality$CheckM2Contamination, na.rm = TRUE) * 1.08)), ylim = c(0, 102)) +
  labs(
    title = "Completeness, contamination, and chimerism answer different questions",
    subtitle = "Green region marks the CheckM2 minimum; filled points also pass GUNC",
    x = "CheckM2 contamination (%)", y = "CheckM2 completeness (%)",
    color = NULL, shape = NULL, size = "Bin size"
  ) +
  theme_pub() +
  theme(panel.grid.major = element_line(color = "#EEEEEE"), legend.position = "right")
save_pub(p_quality, file.path(output_dir, "43-quality-landscape"), width = 205, height = 145)

provenance <- read_tsv(file.path(input_dir, "refinement-provenance.tsv"), show_col_types = FALSE)
provenance_y_min <- max(0, floor(min(provenance$DominantInputCoveragePct, na.rm = TRUE) / 5) * 5 - 5)
p_provenance <- ggplot(
  provenance,
  aes(x = ContributingInputBins, y = DominantInputCoveragePct, color = Method)
) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.76, size = 2.35) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(provenance_y_min, 102), breaks = pretty_breaks()) +
  labs(
    title = "Refined bins retain an auditable source lineage",
    subtitle = "Dominant overlap by assembled bp; full observed y-range shown",
    x = "Overlapping input bins", y = "Dominant input coverage (%)", color = NULL
  ) +
  theme_pub() + theme(legend.position = "bottom")
save_pub(p_provenance, file.path(output_dir, "43-refinement-provenance"), width = 180, height = 130)

refinement <- read_tsv(file.path(input_dir, "refinement-summary.tsv"), show_col_types = FALSE)
selection <- read_tsv(file.path(input_dir, "final-method-selection.tsv"), show_col_types = FALSE)
p_count <- ggplot(refinement, aes(x = Method, y = MinimumPassBins, fill = Method)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = MinimumPassBins), vjust = -0.35, size = 4) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Primary criterion", subtitle = "Minimum-pass bins", x = NULL, y = "Bins") +
  theme_pub()
p_score <- ggplot(refinement, aes(x = Method, y = PassingScoreSum, fill = Method)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = number(PassingScoreSum, accuracy = 0.1)), vjust = -0.35, size = 4) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Predeclared tie-break", subtitle = "Sum of completeness - 5 x contamination", x = NULL, y = "Score") +
  theme_pub()
p_selection <- (p_count | p_score) +
  plot_annotation(
    title = paste0("Selected refinement: ", selection$SelectedMethod),
    subtitle = "Known mock truth was not used for method selection"
  )
save_pub(p_selection, file.path(output_dir, "43-method-selection"), width = 205, height = 125)
