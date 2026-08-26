#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop("Missing argument: ", flag)
  args[[index + 1]]
}
work_dir <- value_after("--work-dir")
input_dir <- file.path(work_dir, "summary")
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260751)

country_pal <- c(
  "AUT" = "#CC79A7", "CHN" = "#D55E00", "DEU" = "#E69F00",
  "FJI" = "#009E73", "GBR" = "#56B4E9", "NLD" = "#0072B2",
  "RUS" = "#F0E442", "SWE" = "#6A3D9A", "USA" = "#8C564B",
  "Reference" = "#555555"
)
stratum_pal <- c(
  "Same study" = "#0072B2",
  "Same country, different study" = "#009E73",
  "Different country" = "#D55E00"
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
save_pub <- function(plot, file_base, width = 190, height = 135, dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = "mm", dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = "mm", dpi = dpi, compression = "lzw", bg = "white")
}

paths <- read_tsv(file.path(work_dir, "output-paths.tsv"), show_col_types = FALSE)
tips <- read_tsv(file.path(input_dir, "tip-metadata.tsv"), show_col_types = FALSE)
tree <- read.tree(paths$Tree[[1]])
tip_annotation <- tips |>
  transmute(label = Tip, Type, Study, Country)
tree_height <- max(node.depth.edgelength(tree))
p_tree <- ggtree(tree, linewidth = 0.48) %<+% tip_annotation +
  geom_tippoint(aes(color = Country, shape = Type), size = 2.5, alpha = 0.9) +
  geom_tiplab(aes(color = Country), size = 2.55, align = FALSE,
              offset = tree_height * 0.015, show.legend = FALSE) +
  scale_color_manual(values = country_pal, drop = FALSE) +
  scale_shape_manual(values = c("Metagenome sample" = 16, "Reference genome" = 17)) +
  coord_cartesian(xlim = c(0, tree_height * 1.47), clip = "off") +
  labs(
    title = "Marker-gene phylogeny resolves within-species structure",
    subtitle = "25 stool metagenomes from 13 studies; references are triangles",
    x = "Substitutions per site", y = NULL, color = "Country", shape = NULL
  ) +
  theme_tree2() +
  theme(
    text = element_text(size = 11),
    plot.title.position = "plot", axis.text = element_text(color = "black"),
    axis.title.x = element_text(size = 10.5, color = "black", margin = margin(t = 5)),
    legend.position = "right", plot.margin = margin(6, 72, 16, 6)
  )
save_pub(p_tree, file.path(output_dir, "51-strainphlan-tree"), 225, 205)

pairs <- read_tsv(file.path(input_dir, "pairwise-p-distance.tsv"), show_col_types = FALSE)
sample_order <- sort(unique(c(pairs$Sample1, pairs$Sample2)))
distance_matrix <- matrix(0, nrow = length(sample_order), ncol = length(sample_order),
                          dimnames = list(sample_order, sample_order))
for (index in seq_len(nrow(pairs))) {
  first <- pairs$Sample1[[index]]
  second <- pairs$Sample2[[index]]
  distance_matrix[first, second] <- pairs$PDistance[[index]]
  distance_matrix[second, first] <- pairs$PDistance[[index]]
}
ordination <- cmdscale(as.dist(distance_matrix), k = 2, eig = TRUE, add = TRUE)
positive_eigenvalues <- ordination$eig[ordination$eig > 0]
explained <- 100 * ordination$eig[1:2] / sum(positive_eigenvalues)
sample_metadata <- tips |>
  filter(Type == "Metagenome sample") |>
  select(Sample = Tip, Study, Country)
pcoa <- as.data.frame(ordination$points) |>
  tibble::rownames_to_column("Sample") |>
  rename(Axis1 = V1, Axis2 = V2) |>
  left_join(sample_metadata, by = "Sample")
p_pcoa <- ggplot(pcoa, aes(Axis1, Axis2, color = Country)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#BBBBBB") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#BBBBBB") +
  geom_point(size = 3.2, alpha = 0.9) +
  geom_text_repel(
    aes(label = Sample), size = 2.35, seed = 20260751,
    box.padding = 0.24, point.padding = 0.17, min.segment.length = 0,
    force = 2, max.time = 5, max.iter = 100000,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(values = country_pal[names(country_pal) != "Reference"]) +
  labs(
    title = "Pairwise marker divergence provides an ordination view",
    subtitle = "Additive PCoA of gap-aware nucleotide p-distances",
    x = sprintf("PCoA 1 (%.1f%% of positive eigenvalues)", explained[[1]]),
    y = sprintf("PCoA 2 (%.1f%% of positive eigenvalues)", explained[[2]]),
    color = "Country"
  ) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_pcoa, file.path(output_dir, "51-marker-distance-pcoa"), 195, 145)

pairs <- pairs |>
  mutate(
    PairStratum = factor(
      PairStratum,
      levels = c("Same study", "Same country, different study", "Different country")
    )
  )
pair_counts <- pairs |>
  count(PairStratum) |>
  mutate(Label = paste0("n = ", n))
p_pair <- ggplot(pairs, aes(PairStratum, DifferencesPer10kb, fill = PairStratum)) +
  geom_violin(width = 0.82, alpha = 0.28, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.72, linewidth = 0.45) +
  geom_jitter(aes(color = PairStratum), width = 0.11, height = 0,
              size = 1.35, alpha = 0.42, show.legend = FALSE) +
  geom_text(
    data = pair_counts,
    aes(x = PairStratum, y = Inf, label = Label),
    vjust = 1.5, inherit.aes = FALSE, size = 3.1
  ) +
  scale_fill_manual(values = stratum_pal) +
  scale_color_manual(values = stratum_pal) +
  scale_x_discrete(labels = c(
    "Same study" = "Same\nstudy",
    "Same country, different study" = "Same country,\ndifferent study",
    "Different country" = "Different\ncountry"
  )) +
  labs(
    title = "Study and geography strata are descriptive, not causal",
    subtitle = "Country, cohort and laboratory are partly confounded in this tutorial panel",
    x = NULL, y = "Pairwise marker differences per 10,000 comparable sites"
  ) +
  theme_pub() + theme(legend.position = "none")
save_pub(p_pair, file.path(output_dir, "51-pairwise-distance-strata"), 185, 140)

polymorphism <- read_tsv(
  file.path(input_dir, "polymorphism-by-sample.tsv"), show_col_types = FALSE
) |>
  mutate(sample = reorder(sample, percentage_of_polymorphic_sites))
p_poly <- ggplot(
  polymorphism,
  aes(percentage_of_polymorphic_sites, sample, color = Country)
) +
  geom_segment(aes(x = 0, xend = percentage_of_polymorphic_sites,
                   y = sample, yend = sample), color = "#D9D9D9", linewidth = 0.65) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = country_pal[names(country_pal) != "Reference"]) +
  scale_x_continuous(labels = label_number(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Within-sample polymorphism is a separate quality signal",
    subtitle = "High values can reflect mixed strains, mapping ambiguity or true diversity",
    x = "Polymorphic marker sites (%)", y = NULL, color = "Country"
  ) +
  theme_pub(10.5) + theme(panel.grid.major.y = element_blank())
save_pub(p_poly, file.path(output_dir, "51-polymorphism-audit"), 195, 170)
