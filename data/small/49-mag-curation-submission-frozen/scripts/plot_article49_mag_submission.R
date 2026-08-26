#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(forcats)
  library(scales)
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
set.seed(20260749)

pal <- c(
  "Pass" = "#009E73", "Review" = "#E69F00", "Blocked" = "#D55E00",
  "High quality" = "#0072B2", "Medium quality" = "#CC79A7"
)
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EEEEEE"),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
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

catalog <- read_tsv(file.path(input_dir, "catalog-disposition.tsv"), show_col_types = FALSE)
funnel <- tibble::tribble(
  ~Stage, ~MAGs, ~Status,
  "95%-ANI representatives", nrow(catalog), "Pass",
  "Complete Article 44 audit", sum(catalog$Article44CompleteAudit), "Pass",
  "NCBI numeric gate", sum(catalog$Disposition == "TECHNICAL_REVIEW_SET"), "Review",
  "External submission ready", 0, "Blocked"
) |>
  mutate(Stage = factor(Stage, levels = rev(Stage)))
p_funnel <- ggplot(funnel, aes(x = MAGs, y = Stage, fill = Status)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = MAGs), hjust = if_else(funnel$MAGs > 2, 1.2, -0.4),
            color = if_else(funnel$MAGs > 2, "white", "black"), fontface = "bold") +
  scale_fill_manual(values = pal[c("Pass", "Review", "Blocked")]) +
  scale_x_continuous(limits = c(0, 26), breaks = seq(0, 24, 4), expand = c(0, 0)) +
  labs(
    title = "Quality gates narrow a catalog; metadata gates can stop submission",
    subtitle = "Zero is the correct external-ready count until manual review and accessions exist",
    x = "Species-level genome bins", y = NULL, fill = NULL
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_funnel, file.path(output_dir, "49-submission-funnel"), width = 190, height = 120)

curation <- read_tsv(file.path(input_dir, "manual-review-sheet.tsv"), show_col_types = FALSE) |>
  mutate(Isolate = fct_reorder(Isolate, AutomatedOutlierBpPct))
p_flags <- ggplot(curation, aes(x = AutomatedOutlierBpPct, y = Isolate)) +
  geom_segment(aes(x = 0, xend = AutomatedOutlierBpPct, yend = Isolate), color = "#BBBBBB") +
  geom_point(aes(size = AutomatedOutlierContigs, color = GraphBoundaryPct), alpha = 0.9) +
  scale_color_viridis_c(option = "C", end = 0.9, name = "Graph boundary\nread-pair weight (%)") +
  scale_size_area(max_size = 8, breaks = pretty_breaks(4), name = "Flagged contigs") +
  scale_x_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Automated GC/depth flags prioritize—not replace—manual review",
    subtitle = "Flagged sequence is retained; no contig is removed automatically",
    x = "Genome length flagged by robust screen (%)", y = NULL
  ) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_flags, file.path(output_dir, "49-curation-flags"), width = 200, height = 145)

checklist <- read_tsv(file.path(input_dir, "submission-readiness-checklist.tsv"), show_col_types = FALSE) |>
  mutate(
    DisplayStatus = case_when(
      Status == "PASS" ~ "Pass",
      grepl("PASS", Status) ~ "Review",
      TRUE ~ "Blocked"
    ),
    Gate = fct_rev(factor(Gate, levels = Gate)),
    Label = if_else(DisplayStatus == "Pass", "PASS", if_else(DisplayStatus == "Review", "REVIEW", "BLOCKED"))
  )
p_ready <- ggplot(checklist, aes(x = 1, y = Gate, fill = DisplayStatus)) +
  geom_tile(width = 0.9, height = 0.8, color = "white") +
  geom_text(aes(label = Label), color = "white", fontface = "bold", size = 3.4) +
  scale_fill_manual(values = pal[c("Pass", "Review", "Blocked")]) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "Sequence quality alone does not make a submission package",
    subtitle = "Manual signoff, taxonomy coordination, source metadata and accessions remain blocking gates",
    x = NULL, y = NULL, fill = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid = element_blank(), legend.position = "top"
  )
save_pub(p_ready, file.path(output_dir, "49-readiness-matrix"), width = 210, height = 150)

mimag <- read_tsv(file.path(input_dir, "mimag-quality-supplement.tsv"), show_col_types = FALSE)
p_quality <- ggplot(mimag, aes(x = CheckM2Contamination, y = CheckM2Completeness,
                              color = MIMAGQuality, label = Isolate)) +
  annotate("rect", xmin = -Inf, xmax = 5, ymin = 90, ymax = Inf,
           fill = "#009E73", alpha = 0.08) +
  geom_hline(yintercept = 90, linetype = 2, color = "#555555") +
  geom_vline(xintercept = 5, linetype = 2, color = "#555555") +
  geom_point(size = 3.2, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 2.7, max.overlaps = Inf, seed = 20260749,
                           box.padding = 0.35, min.segment.length = 0,
                           show.legend = FALSE) +
  scale_color_manual(values = pal[c("High quality", "Medium quality")]) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_continuous(limits = c(89.5, 100.5), breaks = seq(90, 100, 2)) +
  labs(
    title = "All 12 review MAGs clear the NCBI numeric genome gate",
    subtitle = "MIMAG high quality additionally requires complete rRNA and tRNA evidence",
    x = "CheckM2 contamination (%)", y = "CheckM2 completeness (%)", color = "MIMAG tier"
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_quality, file.path(output_dir, "49-mimag-quality"), width = 195, height = 135)
