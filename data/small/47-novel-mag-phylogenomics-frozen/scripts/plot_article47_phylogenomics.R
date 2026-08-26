#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(scales)
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
set.seed(20260747)

domain_pal <- c("Bacteria" = "#D55E00", "Archaea" = "#0072B2")
type_pal <- c("Query MAG" = "#D55E00", "Reference genome" = "#666666")
quality_pal <- c("Quality eligible" = "#009E73", "Below quality gate" = "#E69F00")

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
save_pub <- function(plot, file_base, width = 195, height = 145, dpi = 350) {
  ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
         units = "mm", dpi = dpi, bg = "white")
  ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
         units = "mm", dpi = dpi, compression = "lzw", bg = "white")
}
as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "yes", "1")

paths <- read_tsv(file.path(work_dir, "alignment-paths.tsv"), show_col_types = FALSE)
domain_summary <- read_tsv(
  file.path(input_dir, "domain-tree-summary.tsv"), show_col_types = FALSE
)
tip_ledger <- read_tsv(
  file.path(input_dir, "tree-tip-ledger.tsv"), show_col_types = FALSE
)

plot_tree <- function(domain, file_stem, query_color, width, height) {
  path_row <- paths |> filter(Domain == domain)
  summary_row <- domain_summary |> filter(Domain == domain)
  annotation <- tip_ledger |>
    filter(Domain == domain) |>
    mutate(
      DisplayLabel = if_else(
        Type == "Query MAG", StableID,
        sub("^REF_", "", StableID)
      )
    ) |>
    transmute(label = Tip, DisplayLabel, Type)
  tree <- read.tree(path_row$TreeFile[[1]])
  tree_height <- max(node.depth.edgelength(tree))
  type_colors <- c("Query MAG" = query_color, "Reference genome" = "#666666")
  p <- ggtree(tree, linewidth = 0.42)
  p$data$UFBoot <- ifelse(
    !p$data$isTip & grepl("/", p$data$label),
    suppressWarnings(as.numeric(sub(".*/", "", p$data$label))),
    NA_real_
  )
  p <- p %<+% annotation +
    geom_tippoint(aes(color = Type, shape = Type), size = 2.25, alpha = 0.92) +
    geom_tiplab(
      aes(label = DisplayLabel, color = Type), size = 2.45,
      offset = tree_height * 0.012, show.legend = FALSE
    ) +
    geom_text2(
      aes(
        subset = !isTip & !is.na(UFBoot) & UFBoot < 95,
        label = round(UFBoot)
      ),
      color = "#CC79A7", size = 2.35, hjust = -0.12, vjust = -0.3
    ) +
    scale_color_manual(values = type_colors) +
    scale_shape_manual(values = c("Query MAG" = 16, "Reference genome" = 17)) +
    coord_cartesian(xlim = c(0, tree_height * 1.52), clip = "off") +
    labs(
      title = paste0(domain, " phylogenomics with domain-specific SCGs"),
      subtitle = sprintf(
        "%d query MAGs + %d references; %s aa sites; %s; labels show UFBoot <95",
        as.integer(summary_row$QueryMAGs[[1]]),
        as.integer(summary_row$UniqueReferenceTips[[1]]),
        comma(summary_row$AlignmentSites[[1]]), summary_row$BestFitModel[[1]]
      ),
      x = "Substitutions per site", y = NULL, color = NULL, shape = NULL
    ) +
    theme_tree2() +
    theme(
      text = element_text(size = 10.5), plot.title.position = "plot",
      axis.text = element_text(color = "black"), axis.title = element_text(color = "black"),
      axis.title.x = element_text(color = "black", margin = margin(t = 5)),
      legend.position = "right", plot.margin = margin(6, 76, 16, 6)
    ) +
    xlab("Substitutions per site")
  save_pub(p, file.path(output_dir, file_stem), width, height)
}

plot_tree("Bacteria", "47-bacteria-phylogenomics", "#D55E00", 230, 205)
plot_tree("Archaea", "47-archaea-phylogenomics", "#0072B2", 215, 155)

novelty <- read_tsv(
  file.path(input_dir, "novelty-decision-audit.tsv"), show_col_types = FALSE
) |>
  mutate(
    ANIMarginPct = FastANIANI - SpeciesRadiusPct,
    TreeStatus = if_else(
      as_bool(PhylogenomicTreeIncluded), "Included in tree", "Excluded by SCG filter"
    ),
    QualityStatus = if_else(
      as_bool(QualityEligibleForNoveltyReview), "Quality eligible", "Below quality gate"
    ),
    Label = if_else(
      SGB %in% c("SGB_006", "SGB_015"), SGB, NA_character_
    )
  )
novelty_labels <- novelty |> filter(!is.na(Label))
p_novelty <- ggplot(
  novelty,
  aes(AlignmentFractionPct, ANIMarginPct, color = Domain, shape = TreeStatus)
) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.55, color = "#555555") +
  geom_point(aes(fill = QualityStatus), size = 3.4, stroke = 1.05) +
  geom_text_repel(
    data = novelty_labels,
    aes(AlignmentFractionPct, ANIMarginPct, label = Label),
    inherit.aes = FALSE, seed = 20260747, size = 3.2,
    box.padding = 0.32, min.segment.length = 0, max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = domain_pal) +
  scale_fill_manual(values = quality_pal) +
  scale_shape_manual(values = c("Included in tree" = 21, "Excluded by SCG filter" = 24)) +
  scale_x_continuous(labels = label_number(accuracy = 1), expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = "All 24 SGBs clear their reference-specific species radius",
    subtitle = "Positive ANI margin plus sufficient alignment fraction supports existing species assignments",
    x = "FastANI alignment fraction (%)",
    y = "ANI minus GTDB species radius (percentage points)",
    color = "Domain", fill = "MAG quality", shape = "Phylogeny status"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, color = "#555555"))
  ) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_novelty, file.path(output_dir, "47-novelty-gates"), 205, 145)

genomes <- read_tsv(file.path(work_dir, "genome-ledger.tsv"), show_col_types = FALSE) |>
  select(SGB, Domain, CheckM2Completeness = Completeness)
gtt <- read_tsv(
  file.path(input_dir, "gtotree-genome-audit.tsv"), show_col_types = FALSE
) |>
  filter(Type == "Query MAG") |>
  transmute(
    SGB = AssemblyID, Domain, SCGCompletenessPct,
    UniqueSCGHits, SCGsAfterLengthFilter,
    TreeStatus = if_else(as_bool(InFinalTree), "Included in tree", "Excluded by SCG filter")
  ) |>
  left_join(genomes, by = c("SGB", "Domain")) |>
  mutate(Label = if_else(TreeStatus == "Excluded by SCG filter", SGB, NA_character_))
gtt_labels <- gtt |> filter(!is.na(Label))
p_recovery <- ggplot(
  gtt,
  aes(CheckM2Completeness, SCGCompletenessPct, color = Domain, shape = TreeStatus)
) +
  geom_point(size = 3.4, alpha = 0.9) +
  geom_text_repel(
    data = gtt_labels,
    aes(CheckM2Completeness, SCGCompletenessPct, label = Label),
    inherit.aes = FALSE, seed = 20260747, size = 3.2,
    box.padding = 0.4, min.segment.length = 0, max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = domain_pal) +
  scale_shape_manual(values = c("Included in tree" = 16, "Excluded by SCG filter" = 17)) +
  labs(
    title = "MAG completeness does not guarantee phylogenomic marker recovery",
    subtitle = "SGB_015: 38 unique SCGs, 33 after length filtering; excluded from the tree",
    x = "CheckM2 completeness (%)", y = "GToTree SCG completeness (%)",
    color = "Domain", shape = NULL
  ) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_recovery, file.path(output_dir, "47-marker-recovery-audit"), 195, 140)
