#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: plot_article54_virus_discovery.R SUMMARY_DIR FIGURE_DIR")
}
summary_dir <- normalizePath(args[[1]], mustWork = TRUE)
figure_dir <- args[[2]]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260754)

pal_pub <- c(
  "Blue" = "#0072B2",
  "Orange" = "#D55E00",
  "Green" = "#009E73",
  "Sky" = "#56B4E9",
  "Yellow" = "#E69F00",
  "Purple" = "#CC79A7",
  "Gray" = "#999999",
  "Light" = "#E9EEF2",
  "Dark" = "#253238"
)
scale_color_pub <- function(...) scale_color_manual(values = pal_pub, ...)
scale_fill_pub <- function(...) scale_fill_manual(values = pal_pub, ...)
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#ECEFF1", linewidth = 0.3),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA),
      plot.title.position = "plot",
      plot.caption = element_text(color = "#455A64", hjust = 0)
    )
}
save_pub <- function(plot, file_base, width = 205, height = 145,
                     units = "mm", dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = units, device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = units, dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = units, dpi = dpi, compression = "lzw", bg = "white")
  invisible(plot)
}
wrap_text <- function(x, width = 34) {
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), "")
}

evidence <- read_tsv(file.path(summary_dir, "virus-evidence-matrix.tsv"),
                     show_col_types = FALSE)
overlap <- read_tsv(file.path(summary_dir, "discovery-overlap.tsv"),
                    show_col_types = FALSE) %>%
  mutate(
    Pattern = factor(Pattern, levels = c("Both", "geNomad only", "VirSorter2 only", "Neither")),
    FillKey = recode(as.character(Pattern),
                     "Both" = "Green", "geNomad only" = "Blue",
                     "VirSorter2 only" = "Orange", "Neither" = "Gray")
  )

p_overlap <- ggplot(overlap, aes(Pattern, Contigs, fill = FillKey)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.4) +
  geom_text(aes(label = Contigs), vjust = -0.35, fontface = "bold", size = 4) +
  scale_fill_pub(guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "A. Discovery calls overlap, but are not identical",
    x = NULL, y = "Contigs"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

pattern_rank <- c("Both" = 1, "geNomad only" = 2, "VirSorter2 only" = 3, "Neither" = 4)
matrix_data <- evidence %>%
  mutate(
    PatternRank = unname(pattern_rank[DiscoveryPattern]),
    RowOrder = row_number()
  ) %>%
  arrange(PatternRank, desc(as.numeric(geNomadVirusScore)), desc(as.numeric(VirSorter2MaxScoreAll))) %>%
  mutate(ContigIndex = row_number()) %>%
  select(ContigIndex, geNomadDetected, VirSorter2Detected, VirSorter2HighConfidence) %>%
  pivot_longer(-ContigIndex, names_to = "Call", values_to = "Detected") %>%
  mutate(
    Call = recode(Call,
                  "geNomadDetected" = "geNomad final",
                  "VirSorter2Detected" = "VirSorter2 >=0.5",
                  "VirSorter2HighConfidence" = "VirSorter2 high confidence"),
    Call = factor(Call, levels = rev(c("geNomad final", "VirSorter2 >=0.5", "VirSorter2 high confidence"))),
    Detected = tolower(as.character(Detected)) == "true"
  )

p_matrix <- ggplot(matrix_data, aes(ContigIndex, Call, fill = Detected)) +
  geom_tile(color = "white", linewidth = 0.22) +
  scale_fill_manual(values = c(`TRUE` = pal_pub[["Blue"]], `FALSE` = pal_pub[["Light"]]),
                    labels = c(`TRUE` = "Called", `FALSE` = "Not called"), name = NULL) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 46)) +
  labs(
    title = "B. The high-confidence rule changes the retained set",
    subtitle = "VirSorter2 high confidence: score >=0.9, or score >=0.7 with at least one hallmark gene",
    x = "Contigs ordered by discovery pattern", y = NULL
  ) +
  theme_pub() +
  theme(panel.grid = element_blank(), legend.position = "bottom")

p_discovery <- p_overlap / p_matrix +
  plot_annotation(
    title = "Two discovery tools provide complementary evidence",
    caption = "CheckV is intentionally absent: it assesses supplied sequences and is not a virus detector."
  )
save_pub(p_discovery, file.path(figure_dir, "54-discovery-consensus"), height = 175)

quality_counts <- read_tsv(file.path(summary_dir, "checkv-quality-counts.tsv"),
                           show_col_types = FALSE) %>%
  mutate(
    CheckVQuality = factor(CheckVQuality, levels = rev(c("High-quality", "Medium-quality", "Low-quality", "Not-determined"))),
    FillKey = recode(as.character(CheckVQuality),
                     "High-quality" = "Green", "Medium-quality" = "Blue",
                     "Low-quality" = "Yellow", "Not-determined" = "Gray")
  )
p_quality_bar <- ggplot(quality_counts, aes(Contigs, CheckVQuality, fill = FillKey)) +
  geom_col(width = 0.68, color = "white") +
  geom_text(aes(label = Contigs), hjust = -0.25, fontface = "bold", size = 4) +
  scale_fill_pub(guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "A. CheckV quality tiers", x = "Contigs", y = NULL) +
  theme_pub()

quality_points <- evidence %>%
  mutate(
    Completeness = suppressWarnings(as.numeric(CompletenessPct)),
    Contamination = suppressWarnings(as.numeric(ContaminationPct)),
    Provirus = factor(CheckVProvirus, levels = c("No", "Yes")),
    QualityKey = recode(CheckVQuality,
                        "High-quality" = "Green", "Medium-quality" = "Blue",
                        "Low-quality" = "Yellow", "Not-determined" = "Gray"),
    Label = if_else(ContigID == "UHGV-0000346", "UHGV-0000346", "")
  ) %>%
  filter(!is.na(Completeness))
p_quality_scatter <- ggplot(quality_points,
                            aes(Completeness, Contamination, color = QualityKey, shape = Provirus)) +
  annotate("rect", xmin = 90, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = pal_pub[["Green"]], alpha = 0.07) +
  geom_vline(xintercept = c(50, 90), linetype = "dashed", color = "#607D8B", linewidth = 0.45) +
  geom_point(size = 3.1, alpha = 0.86) +
  geom_text(data = filter(quality_points, Label != ""), aes(label = Label),
            nudge_x = -4, nudge_y = 5, hjust = 1, color = pal_pub[["Dark"]], size = 3.2) +
  scale_color_pub(name = "Quality", breaks = c("Green", "Blue", "Yellow", "Gray"),
                  labels = c("High", "Medium", "Low", "Not determined")) +
  scale_shape_manual(values = c("No" = 16, "Yes" = 17), name = "Provirus") +
  coord_cartesian(xlim = c(0, 105), ylim = c(0, 75), clip = "off") +
  labs(
    title = "B. Completeness and host flanks are separate",
    x = "Estimated completeness (%)", y = "Estimated contamination (%)"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical")

p_quality <- p_quality_bar + p_quality_scatter +
  plot_layout(widths = c(0.78, 1.55)) +
  plot_annotation(
    title = "Quality assessment starts after virus discovery",
    caption = "UHGV-0000346 is a complete proviral region inside a contig with 30.29% host-flank contamination.\nTrim boundaries before gene or vOTU analysis."
  )
save_pub(p_quality, file.path(figure_dir, "54-checkv-quality"), width = 235, height = 160)

votu_pairs <- read_tsv(file.path(summary_dir, "votu-pairwise-threshold-audit.tsv"),
                       show_col_types = FALSE) %>%
  mutate(
    Gate = case_when(
      tolower(as.character(SameVOTU)) == "true" ~ "Both gates",
      tolower(as.character(PassANI95)) == "true" ~ "ANI only",
      TRUE ~ "Neither"
    ),
    ColorKey = recode(Gate, "Both gates" = "Green", "ANI only" = "Orange", "Neither" = "Gray"),
    Boundary = tolower(as.character(BoundaryPair)) == "true"
  )
p_votu <- ggplot(votu_pairs, aes(ANIpct, ShorterAlignmentFractionPct)) +
  annotate("rect", xmin = 95, xmax = Inf, ymin = 85, ymax = Inf,
           fill = pal_pub[["Green"]], alpha = 0.10) +
  geom_vline(xintercept = 95, linetype = "dashed", color = pal_pub[["Dark"]]) +
  geom_hline(yintercept = 85, linetype = "dashed", color = pal_pub[["Dark"]]) +
  geom_point(aes(color = ColorKey, shape = Boundary), size = 3.7, alpha = 0.9) +
  geom_text(data = filter(votu_pairs, Boundary),
            aes(label = "95.16% ANI\n84.76% shorter-sequence AF"),
            nudge_x = -0.5, nudge_y = -9, hjust = 1, size = 3.4, color = pal_pub[["Dark"]]) +
  annotate("text", x = 98.0, y = 96, label = "Same-vOTU acceptance region",
           color = pal_pub[["Green"]], fontface = "bold", size = 3.5) +
  scale_color_pub(name = NULL, breaks = c("Green", "Orange", "Gray"),
                  labels = c("Both gates", "ANI gate only", "Neither gate")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  scale_x_continuous(limits = c(79, 101), breaks = c(80, 85, 90, 95, 100)) +
  scale_y_continuous(limits = c(0, 102), breaks = c(0, 25, 50, 75, 85, 100)) +
  labs(
    title = "vOTU assignment requires identity and shared span",
    subtitle = "Ten non-self alignable pairs from the 46-sequence official CheckV fixture",
    x = "Average nucleotide identity (%)",
    y = "Alignment fraction of the shorter sequence (%)",
    caption = "The highlighted UHGV pair passes 95% ANI but misses the 85% alignment-fraction gate.\nAll 46 sequences therefore remain separate vOTUs."
  ) +
  theme_pub() +
  theme(legend.position = "bottom")
save_pub(p_votu, file.path(figure_dir, "54-votu-threshold"), width = 195, height = 150)

ladder <- read_tsv(file.path(summary_dir, "virus-evidence-ladder.tsv"),
                   show_col_types = FALSE) %>%
  mutate(
    y = rev(seq_len(n())),
    RequiredWrapped = wrap_text(RequiredEvidence, 46),
    ForbiddenWrapped = wrap_text(ForbiddenShortcut, 34)
  )
ladder_spine <- tibble(
  x = 0.18, xend = 0.18,
  y = max(ladder$y), yend = min(ladder$y)
)
p_ladder <- ggplot(ladder) +
  geom_segment(data = ladder_spine,
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 1.1, color = pal_pub[["Sky"]],
               arrow = grid::arrow(type = "closed", length = grid::unit(0.16, "inches")),
               inherit.aes = FALSE) +
  geom_point(aes(0.18, y), size = 8, color = "white", fill = pal_pub[["Blue"]], shape = 21, stroke = 1.2) +
  geom_text(aes(0.18, y, label = Order), color = "white", fontface = "bold", size = 3.3) +
  geom_text(aes(0.30, y + 0.17, label = Decision), hjust = 0, fontface = "bold", size = 3.8, color = pal_pub[["Dark"]]) +
  geom_text(aes(0.30, y - 0.13, label = RequiredWrapped), hjust = 0, vjust = 1, size = 3.0, lineheight = 0.95) +
  geom_label(aes(1.05, y, label = ForbiddenWrapped), hjust = 0, size = 2.75,
             fill = "#FFF4E8", color = pal_pub[["Dark"]], label.size = 0.15) +
  annotate("text", x = 1.05, y = max(ladder$y) + 0.65, label = "Do not substitute",
           hjust = 0, fontface = "bold", color = pal_pub[["Orange"]], size = 3.6) +
  coord_cartesian(xlim = c(0.05, 1.82), ylim = c(0.45, max(ladder$y) + 0.9), clip = "off") +
  labs(
    title = "A viral contig becomes publishable through five distinct decisions",
    subtitle = "The weakest unresolved decision limits the final claim"
  ) +
  theme_void(base_size = 12) +
  theme(plot.title.position = "plot", plot.margin = margin(8, 18, 8, 8))
save_pub(p_ladder, file.path(figure_dir, "54-virus-evidence-ladder"), width = 220, height = 160)

library_audit <- read_tsv(file.path(summary_dir, "library-design-bias-audit.tsv"),
                          show_col_types = FALSE) %>%
  mutate(
    Row = rev(seq_len(n())),
    LibraryLabel = wrap_text(Library, 15),
    StepLabel = wrap_text(Step, 27),
    BlindLabel = wrap_text(KnownBlindSpot, 36),
    QuantLabel = wrap_text(QuantitativeUse, 27),
    FillKey = recode(Library,
                     "Total metagenome" = "Blue",
                     "Virus-enriched virome" = "Green",
                     "Amplified virome" = "Orange")
  )
p_library <- ggplot(library_audit) +
  geom_tile(aes(x = 0.15, y = Row, fill = FillKey), width = 0.25, height = 0.74, color = "white") +
  geom_text(aes(0.31, Row, label = LibraryLabel), hjust = 0, fontface = "bold", size = 3.35) +
  geom_text(aes(0.90, Row, label = StepLabel), hjust = 0, size = 3.1, lineheight = 0.95) +
  geom_text(aes(1.55, Row, label = BlindLabel), hjust = 0, size = 3.0, lineheight = 0.95) +
  geom_text(aes(2.42, Row, label = QuantLabel), hjust = 0, size = 3.0, lineheight = 0.95) +
  annotate("text", x = c(0.31, 0.90, 1.55, 2.42), y = max(library_audit$Row) + 0.72,
           label = c("Library", "Defining step", "Main blind spot", "Abundance use"),
           hjust = 0, fontface = "bold", color = pal_pub[["Dark"]], size = 3.5) +
  scale_fill_pub(guide = "none") +
  coord_cartesian(xlim = c(0, 3.20), ylim = c(0.45, max(library_audit$Row) + 1), clip = "off") +
  labs(
    title = "Total metagenomes and enriched viromes answer different sampling questions",
    subtitle = "Filtration, nuclease treatment, and amplification change the observable population",
    caption = "Methodological evidence map; rows are not a quantitative ranking.\nAmplified read depth should not be interpreted as original particle abundance."
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.caption = element_text(color = "#455A64", hjust = 0),
    plot.margin = margin(8, 12, 8, 8)
  )
save_pub(p_library, file.path(figure_dir, "54-library-design-bias"), width = 240, height = 150)

write_tsv(
  tibble(
    Package = c("R", "ggplot2", "dplyr", "readr", "tidyr", "patchwork"),
    Version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      as.character(packageVersion("ggplot2")),
      as.character(packageVersion("dplyr")),
      as.character(packageVersion("readr")),
      as.character(packageVersion("tidyr")),
      as.character(packageVersion("patchwork"))
    )
  ),
  file.path(summary_dir, "plot-software-versions.tsv")
)

message("Article 54 figures written to ", normalizePath(figure_dir))
