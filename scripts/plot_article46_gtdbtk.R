#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
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
set.seed(20260746)

pal <- c("Bacteria" = "#0072B2", "Archaea" = "#E69F00",
         "Resolved" = "#009E73", "Unresolved" = "#999999")
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#EEEEEE"),
      axis.text = element_text(color = "black"), legend.key = element_blank(),
      strip.background = element_rect(fill = "#EEF3F5", color = NA), plot.title.position = "plot"
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

ranks <- read_tsv(file.path(input_dir, "rank-resolution.tsv"), show_col_types = FALSE) |>
  mutate(
    Domain = sub("^d__", "", Domain),
    Rank = factor(Rank, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  ) |>
  select(Domain, Rank, ResolvedSGBs, UnresolvedSGBs) |>
  pivot_longer(c(ResolvedSGBs, UnresolvedSGBs), names_to = "Status", values_to = "SGBs") |>
  mutate(Status = recode(Status, ResolvedSGBs = "Resolved", UnresolvedSGBs = "Unresolved"))
p_rank <- ggplot(ranks, aes(x = Rank, y = SGBs, fill = Status)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = if_else(SGBs > 0, as.character(SGBs), "")),
            position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 3.2) +
  facet_wrap(~ Domain, scales = "free_y") +
  scale_fill_manual(values = pal[c("Resolved", "Unresolved")]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Taxonomic resolution must be reported rank by rank",
    subtitle = "An empty GTDB species field is not a named novel species",
    x = NULL, y = "Species-level genome bins", fill = NULL
  ) +
  theme_pub() + theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1))
save_pub(p_rank, file.path(output_dir, "46-rank-resolution"), width = 190, height = 125)

phyla <- read_tsv(file.path(input_dir, "phylum-summary.tsv"), show_col_types = FALSE) |>
  mutate(
    Domain = sub("^d__", "", Domain),
    Phylum = sub("^p__", "", Phylum),
    Phylum = if_else(Phylum == "", "Unresolved phylum", Phylum),
    Phylum = fct_reorder(Phylum, SGBs)
  )
p_phylum <- ggplot(phyla, aes(x = SGBs, y = Phylum, fill = Domain)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = SGBs), hjust = -0.2, size = 3.4) +
  scale_fill_manual(values = pal[c("Bacteria", "Archaea")]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16)), breaks = pretty_breaks()) +
  labs(
    title = "GTDB R232 assigns the non-redundant catalog across phyla",
    subtitle = "Counts refer to 95%-ANI representatives, not raw bins",
    x = "Species-level genome bins", y = NULL, fill = "Domain"
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_phylum, file.path(output_dir, "46-phylum-composition"), width = 190, height = 130)

taxonomy <- read_tsv(file.path(input_dir, "taxonomy-summary.tsv"), show_col_types = FALSE) |>
  mutate(
    Domain = sub("^d__", "", Domain),
    SpeciesStatus = if_else(SpeciesAssigned, "Resolved", "Unresolved"),
    Route = case_when(
      grepl("ANI", ClassificationMethod, ignore.case = TRUE) &
        grepl("topology|placement", ClassificationMethod, ignore.case = TRUE) ~ "Topology + ANI",
      grepl("ANI", ClassificationMethod, ignore.case = TRUE) ~ "ANI",
      grepl("topology|placement", ClassificationMethod, ignore.case = TRUE) ~ "Tree topology",
      TRUE ~ "Other"
    )
  )
route <- taxonomy |> count(Route, SpeciesStatus, name = "SGBs")
p_route <- ggplot(route, aes(x = Route, y = SGBs, fill = SpeciesStatus)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = if_else(SGBs > 0, as.character(SGBs), "")),
            position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
  scale_fill_manual(values = pal[c("Resolved", "Unresolved")]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), breaks = pretty_breaks()) +
  labs(
    title = "The ANI prescreen resolves every catalog SGB",
    subtitle = "Marker-tree placement is not invoked for this dataset",
    x = NULL, y = "Species-level genome bins", fill = "Species field"
  ) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_route, file.path(output_dir, "46-classification-route"), width = 180, height = 120)

ani <- read_tsv(file.path(input_dir, "fastani-reference-audit.tsv"), show_col_types = FALSE) |>
  mutate(
    Domain = sub("^d__", "", Domain),
    SpeciesStatus = if_else(SpeciesAssigned, "Resolved", "Unresolved")
  )
if (nrow(ani) > 0) {
  p_ani <- ggplot(ani, aes(x = ANIMarginToReferenceRadiusPctPoints, y = FastANIAFPct,
                          color = Domain, shape = SpeciesStatus)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#555555") +
    geom_hline(yintercept = 50, linetype = 2, color = "#555555") +
    geom_point(size = 3.4, alpha = 0.82) +
    scale_color_manual(values = pal[c("Bacteria", "Archaea")]) +
    scale_shape_manual(values = c("Resolved" = 16, "Unresolved" = 1)) +
    labs(
      title = "ANI must clear a reference-specific radius and alignment fraction",
      subtitle = "Positive x values exceed the GTDB species radius; AF is shown independently",
      x = "ANI minus reference radius (percentage points)",
      y = "Alignment fraction (%)", color = "Domain", shape = "Species field"
    ) +
    theme_pub() + theme(legend.position = "top")
} else {
  p_ani <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No fastANI reference hits", size = 6) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = "ANI screening returned no reference hits", x = NULL, y = NULL) +
    theme_pub() + theme(axis.text = element_blank(), axis.ticks = element_blank())
}
save_pub(p_ani, file.path(output_dir, "46-ani-af-audit"), width = 190, height = 130)
