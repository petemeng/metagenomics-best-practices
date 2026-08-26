#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(ggdendro)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) stop(paste("Missing", flag))
  args[[position + 1]]
}
work_dir <- normalizePath(value_after("--work-dir"), mustWork = TRUE)
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260752)

country_pal <- c(
  "AUT" = "#CC79A7", "CHN" = "#D55E00", "DEU" = "#E69F00",
  "FJI" = "#009E73", "GBR" = "#56B4E9", "NLD" = "#0072B2",
  "RUS" = "#F0E442", "SWE" = "#6A3D9A", "USA" = "#8C564B"
)
category_pal <- c(
  "Core >=95%" = "#0072B2", "Accessory 5-<95%" = "#E69F00",
  "Rare >0-<5%" = "#CC79A7", "Undetected" = "#999999"
)
status_pal <- c("Retained" = "#0072B2", "Excluded" = "#D55E00")

scale_color_pub <- function(...) scale_color_manual(values = country_pal, ...)
scale_fill_pub <- function(...) scale_fill_manual(values = country_pal, ...)
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.3),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}
save_pub <- function(plot, file_base, width = 200, height = 145,
                     units = "mm", dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = units, device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = units, dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = units, dpi = dpi, compression = "lzw", bg = "white")
  invisible(plot)
}

summary_dir <- file.path(work_dir, "summary")
write_tsv(
  tibble(
    Software = c("R", "ggplot2", "dplyr", "readr", "ggrepel", "patchwork", "ggdendro"),
    Version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      as.character(packageVersion("ggplot2")),
      as.character(packageVersion("dplyr")),
      as.character(packageVersion("readr")),
      as.character(packageVersion("ggrepel")),
      as.character(packageVersion("patchwork")),
      as.character(packageVersion("ggdendro"))
    )
  ),
  file.path(summary_dir, "plot-software-versions.tsv")
)
prevalence <- read_tsv(
  file.path(summary_dir, "gene-family-prevalence.tsv.gz"),
  show_col_types = FALSE
) %>%
  mutate(
    PrimaryCategory = factor(
      PrimaryCategory,
      levels = c("Core >=95%", "Accessory 5-<95%", "Rare >0-<5%", "Undetected")
    )
  )

category_counts <- prevalence %>%
  count(PrimaryCategory, name = "GeneFamilies")
p_category <- ggplot(category_counts, aes(PrimaryCategory, GeneFamilies, fill = PrimaryCategory)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.4) +
  geom_text(aes(label = comma(GeneFamilies)), vjust = -0.35, size = 3.5) +
  scale_fill_manual(values = category_pal, guide = "none", drop = FALSE) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Operational pangenome partitions",
    subtitle = "22 samples retained by the primary plateau filters",
    x = NULL, y = "Gene families"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

count_spectrum <- prevalence %>%
  count(PrimaryCount, PrimaryCategory, name = "GeneFamilies")
p_spectrum <- ggplot(
  count_spectrum,
  aes(PrimaryCount, GeneFamilies, fill = PrimaryCategory)
) +
  geom_col(width = 0.82) +
  scale_fill_manual(values = category_pal, name = "Partition", drop = FALSE) +
  scale_x_continuous(breaks = c(0, 1, 5, 10, 15, 21, 22)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Prevalence is a continuum",
    subtitle = "Core is an explicit >=95% operational cutoff, not a natural boundary",
    x = "Samples with gene family present (n = 22)", y = "Gene families"
  ) +
  theme_pub() +
  theme(legend.position = "bottom")

save_pub(
  p_category / p_spectrum + plot_layout(heights = c(0.9, 1.2)),
  file.path(output_dir, "52-pangenome-prevalence"),
  width = 215, height = 205
)

pcoa <- read_tsv(file.path(summary_dir, "pcoa-jaccard.tsv"), show_col_types = FALSE)
axis1 <- unique(pcoa$PCoA1Pct)
axis2 <- unique(pcoa$PCoA2Pct)
p_pcoa <- ggplot(pcoa, aes(PCoA1, PCoA2, color = Country)) +
  geom_hline(yintercept = 0, color = "#DDDDDD", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "#DDDDDD", linewidth = 0.35) +
  geom_point(size = 3.2, alpha = 0.92) +
  geom_text_repel(
    aes(label = Sample), size = 2.55, max.overlaps = Inf,
    box.padding = 0.25, point.padding = 0.15, min.segment.length = 0,
    seed = 20260752, show.legend = FALSE
  ) +
  scale_color_pub(name = "Country") +
  labs(
    title = "Within-species gene-content structure",
    subtitle = "Jaccard PCoA of binary PanPhlAn gene-family calls; color is metadata only",
    x = sprintf("PCoA 1 (%.1f%%)", axis1),
    y = sprintf("PCoA 2 (%.1f%%)", axis2)
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
save_pub(p_pcoa, file.path(output_dir, "52-gene-content-pcoa"), width = 205, height = 155)

distance_table <- read_tsv(
  file.path(summary_dir, "sample-jaccard-matrix.tsv"),
  show_col_types = FALSE
)
sample_metadata <- read_tsv(
  file.path(summary_dir, "sample-filter-audit.tsv"),
  show_col_types = FALSE
) %>% select(Sample, Country)
distance_matrix <- as.matrix(distance_table[, -1])
rownames(distance_matrix) <- distance_table$Sample
storage.mode(distance_matrix) <- "double"
tree <- hclust(as.dist(distance_matrix), method = "average")
dendro <- dendro_data(tree, type = "rectangle")
tree_labels <- label(dendro) %>%
  rename(Sample = label) %>%
  left_join(sample_metadata, by = "Sample")
p_tree <- ggplot() +
  geom_segment(
    data = segment(dendro),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "#4D4D4D", linewidth = 0.45
  ) +
  geom_text(
    data = tree_labels,
    aes(x = x, y = -0.012, label = Sample, color = Country),
    angle = 55, hjust = 1, vjust = 1, size = 2.65
  ) +
  scale_color_pub(name = "Country") +
  scale_y_continuous(expand = expansion(mult = c(0.24, 0.05))) +
  labs(
    title = "A gene-content dendrogram is not a species phylogeny",
    subtitle = "Average linkage on pairwise Jaccard distances; branch height is gene-content dissimilarity",
    x = NULL, y = "Jaccard dissimilarity"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid = element_blank(), legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
save_pub(p_tree, file.path(output_dir, "52-gene-content-dendrogram"), width = 220, height = 160)

heatmap <- read_tsv(
  file.path(summary_dir, "accessory-feature-heatmap.tsv"),
  show_col_types = FALSE
)
sample_levels <- heatmap %>% distinct(Sample, SampleOrder) %>% arrange(SampleOrder) %>% pull(Sample)
feature_levels <- heatmap %>% distinct(FeatureLabel, FeatureOrder) %>% arrange(desc(FeatureOrder)) %>% pull(FeatureLabel)
heatmap <- heatmap %>%
  mutate(
    Sample = factor(Sample, levels = sample_levels),
    FeatureLabel = factor(FeatureLabel, levels = feature_levels),
    Call = factor(if_else(Present == 1, "Present", "Absent"), levels = c("Absent", "Present"))
  )
p_heatmap <- ggplot(heatmap, aes(Sample, FeatureLabel, fill = Call)) +
  geom_tile(color = "white", linewidth = 0.18) +
  scale_fill_manual(values = c("Absent" = "#F0F0F0", "Present" = "#0072B2"), drop = FALSE) +
  labs(
    title = "Variable annotated accessory gene families",
    subtitle = "Twenty deterministic high-variance features; labels show UniRef90 suffix and first available KO/Pfam/GO term",
    x = "Metagenome sample", y = "Gene family | annotation", fill = "Call"
  ) +
  theme_pub(base_size = 10.5) +
  theme(
    panel.grid = element_blank(), axis.text.x = element_text(angle = 58, hjust = 1, vjust = 1),
    legend.position = "bottom"
  )
save_pub(p_heatmap, file.path(output_dir, "52-accessory-gene-heatmap"), width = 230, height = 175)

filter_audit <- read_tsv(
  file.path(summary_dir, "sample-filter-audit.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    PrimaryStatus = factor(if_else(PrimaryRetained, "Retained", "Excluded"), levels = c("Retained", "Excluded")),
    Warning = factor(if_else(MultiStrainWarning, "Potential multi-strain", "No flag"))
  )
p_plateau <- ggplot(
  filter_audit,
  aes(MedianCoverage, LeftCoverage, color = PrimaryStatus, shape = Warning)
) +
  geom_hline(yintercept = 1.25, linetype = "dashed", color = "#D55E00", linewidth = 0.55) +
  geom_hline(yintercept = 1.70, linetype = "dotted", color = "#009E73", linewidth = 0.55) +
  geom_point(size = 3.1, alpha = 0.9) +
  geom_text_repel(
    data = filter_audit %>% filter(!PrimaryRetained),
    aes(label = Sample), size = 3.1, seed = 20260752,
    box.padding = 0.35, min.segment.length = 0, show.legend = FALSE
  ) +
  annotate("text", x = 94, y = 1.235, label = "Primary left maximum = 1.25", hjust = 1, vjust = 1, size = 3.1, color = "#D55E00") +
  annotate("text", x = 94, y = 1.685, label = "Sensitivity left maximum = 1.70", hjust = 1, vjust = 1, size = 3.1, color = "#009E73") +
  scale_color_manual(values = status_pal, name = "Primary status") +
  scale_shape_manual(values = c("No flag" = 16, "Potential multi-strain" = 17), name = "Plateau warning") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  scale_y_continuous(limits = c(1.02, 1.75), expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Three samples fail only the primary left-side plateau gate",
    subtitle = "All pass median/right gates; sensitivity thresholds retain all 25 samples",
    x = "Median normalized coverage", y = "Left-side plateau ratio"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(
    color = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 2, nrow = 1)
  )
save_pub(p_plateau, file.path(output_dir, "52-plateau-sensitivity"), width = 205, height = 160)
