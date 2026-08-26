args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
  stop(
    paste(
      "Usage: Rscript scripts/render_article06_figures.R",
      "<read_budget.tsv> <saponin_tradeoff.tsv> <figure_dir>"
    )
  )
}

read_budget_path <- args[[1]]
saponin_path <- args[[2]]
figure_dir <- args[[3]]

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package 'scales' is required")
}

library(ggplot2)
library(patchwork)
library(scales)

set.seed(20260719)

read_budget <- read.delim(
  read_budget_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
saponin <- read.delim(
  saponin_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(read_budget) == 6L,
  nrow(saponin) == 7L,
  all(c(
    "method_label", "host_mean_pct", "microbial_pct",
    "total_reads_for_10m_microbial_millions"
  ) %in% names(read_budget)),
  all(c(
    "treatment_label", "host_filtered_pct",
    "gram_negative_pct",
    "gram_negative_retention_vs_untreated_pct"
  ) %in% names(saponin))
)

pal_pub <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#F0E442", "#999999"
)

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        color = "#E8ECF1",
        linewidth = 0.25
      ),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      legend.key = element_blank(),
      legend.background = element_rect(
        fill = alpha("white", 0.75),
        color = NA
      ),
      plot.title.position = "plot"
    )
}

save_pub <- function(
  plot,
  file_base,
  width = 190,
  height = 125,
  units = "mm",
  dpi = 300
) {
  dir.create(
    dirname(file_base),
    recursive = TRUE,
    showWarnings = FALSE
  )
  ggsave(
    paste0(file_base, ".pdf"),
    plot,
    width = width,
    height = height,
    units = units,
    device = grDevices::cairo_pdf
  )
  ggsave(
    paste0(file_base, ".png"),
    plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = "white"
  )
  ggsave(
    paste0(file_base, ".tiff"),
    plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  invisible(plot)
}

method_levels <- rev(read_budget$method_label)
read_budget$method_label <- factor(
  read_budget$method_label,
  levels = method_levels
)

fraction_long <- rbind(
  data.frame(
    method_label = read_budget$method_label,
    fraction = "Host",
    percent = read_budget$host_mean_pct,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method_label = read_budget$method_label,
    fraction = "Microbial",
    percent = read_budget$microbial_pct,
    stringsAsFactors = FALSE
  )
)
fraction_long$fraction <- factor(
  fraction_long$fraction,
  levels = c("Host", "Microbial")
)

p_fraction <- ggplot(
  fraction_long,
  aes(percent, method_label, fill = fraction)
) +
  geom_col(width = 0.68) +
  geom_text(
    data = read_budget,
    aes(
      x = 122,
      y = method_label,
      label = sprintf("%.1f%% host", host_mean_pct)
    ),
    inherit.aes = FALSE,
    hjust = 1,
    color = "#172033",
    fontface = "bold",
    size = 2.7
  ) +
  scale_fill_manual(
    values = c("Host" = "#B7BDC7", "Microbial" = "#0072B2")
  ) +
  scale_x_continuous(
    limits = c(0, 124),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "A. Observed read fractions",
    x = "Quality-filtered reads",
    y = NULL,
    fill = NULL
  ) +
  theme_pub(10.8) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    plot.title = element_text(face = "bold", size = 10)
  )

p_budget <- ggplot(
  read_budget,
  aes(total_reads_for_10m_microbial_millions, method_label)
) +
  geom_col(
    width = 0.68,
    fill = "#009E73",
    alpha = 0.88
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1fM",
        total_reads_for_10m_microbial_millions
      )
    ),
    hjust = -0.12,
    fontface = "bold",
    size = 3.15
  ) +
  scale_x_continuous(
    limits = c(0, 132),
    breaks = c(0, 25, 50, 75, 100, 125),
    labels = function(x) paste0(x, "M"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "B. Reads needed for 10M microbial",
    x = "Expected total reads",
    y = NULL
  ) +
  theme_pub(10.8) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", size = 10)
  )

p_efficiency <- p_fraction + p_budget +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(
    title = "Computational filtering cannot recover sequencing capacity",
    subtitle = paste0(
      "Observed saliva host fractions from eight participants;\n",
      "read budgets are arithmetic projections, not purchase thresholds"
    ),
    caption = paste(
      "Source: Marotz et al. (2018), Figure 2 and Results.",
      "Each participant contributed triplicate aliquots per method."
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 13,
        color = "#172033"
      ),
      plot.subtitle = element_text(
        size = 9.6,
        color = "#3D4757",
        lineheight = 1.05
      ),
      plot.caption = element_text(
        size = 8.5,
        color = "#596273",
        hjust = 0
      )
    )
  )

save_pub(
  p_efficiency,
  file.path(figure_dir, "06-host-depletion-efficiency"),
  width = 216,
  height = 128
)

dose_levels <- saponin$treatment_label
saponin$treatment_label <- factor(
  saponin$treatment_label,
  levels = dose_levels
)

p_host <- ggplot(
  saponin,
  aes(treatment_label, host_filtered_pct, group = 1)
) +
  geom_hline(
    yintercept = saponin$host_filtered_pct[saponin$saponin_pct == 0],
    linetype = 2,
    color = "#7B8492",
    linewidth = 0.45
  ) +
  geom_line(color = "#D55E00", linewidth = 0.85) +
  geom_point(
    shape = 21,
    size = 3.6,
    stroke = 0.8,
    color = "#D55E00",
    fill = "white"
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", host_filtered_pct)),
    vjust = -0.85,
    size = 3.05,
    fontface = "bold"
  ) +
  scale_y_continuous(
    limits = c(0, 63),
    breaks = seq(0, 60, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "A. Host reads fall",
    x = "Saponin (w/v)",
    y = "Reads filtered as human"
  ) +
  theme_pub(10.8) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 11)
  )

p_gram_negative <- ggplot(
  saponin,
  aes(treatment_label, gram_negative_pct, group = 1)
) +
  geom_hline(
    yintercept = saponin$gram_negative_pct[saponin$saponin_pct == 0],
    linetype = 2,
    color = "#7B8492",
    linewidth = 0.45
  ) +
  geom_line(color = "#0072B2", linewidth = 0.85) +
  geom_point(
    shape = 21,
    size = 3.6,
    stroke = 0.8,
    color = "#0072B2",
    fill = "white"
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", gram_negative_pct)),
    vjust = -0.85,
    size = 3.05,
    fontface = "bold"
  ) +
  scale_y_continuous(
    limits = c(0, 97),
    breaks = c(0, 25, 50, 75),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "B. Gram-negative representation also falls",
    x = "Saponin (w/v)",
    y = "Gram-negative fraction"
  ) +
  theme_pub(10.8) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 11)
  )

p_tradeoff <- p_host + p_gram_negative +
  plot_annotation(
    title = "Host depletion efficiency is not composition fidelity",
    subtitle = paste(
      "Descriptive dose series from one sputum specimen;",
      "no population-level optimum can be inferred"
    ),
    caption = paste(
      "Source: Longhi et al. (2024), main Table 2",
      "and Supplementary Table 6."
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        color = "#172033"
      ),
      plot.subtitle = element_text(size = 10.5, color = "#3D4757"),
      plot.caption = element_text(
        size = 8.5,
        color = "#596273",
        hjust = 0
      )
    )
  )

save_pub(
  p_tradeoff,
  file.path(figure_dir, "06-saponin-tradeoff"),
  width = 200,
  height = 122
)

boxes <- data.frame(
  xmin = c(
    3.2, 0.35, 3.2, 6.05,
    0.35, 3.2, 6.05, 2.1, 2.1
  ),
  xmax = c(
    6.8, 3.15, 6.8, 8.85,
    3.15, 6.8, 8.85, 7.9, 7.9
  ),
  ymin = c(
    8.25, 6.55, 6.55, 6.55,
    4.65, 4.65, 4.65, 2.65, 0.65
  ),
  ymax = c(
    9.35, 7.65, 7.65, 7.65,
    5.75, 5.75, 5.75, 3.75, 1.75
  ),
  label = c(
    "Split representative\nsample aliquots",
    "Untreated\nreference",
    "Candidate wet-lab\ndepletion",
    "Matrix blank +\ncellular mock",
    "Host and microbial\nqPCR/ddPCR",
    "Library yield +\nnon-host reads",
    "Taxon / virus\nrecovery audit",
    "Advance only if usable reads rise\nwithout unacceptable target loss",
    "Lock one protocol across groups\nthen filter host reads and audit privacy"
  ),
  type = c(
    "start", "branch", "branch", "control",
    "measure", "measure", "measure", "decision", "finish"
  ),
  text_size = c(
    3.35, 3.35, 3.35, 3.35,
    3.2, 3.2, 3.2, 3.05, 3.05
  ),
  stringsAsFactors = FALSE
)

arrows <- data.frame(
  x = c(5, 5, 5, 7.45, 1.75, 5, 7.45, 1.75, 5, 7.45, 5),
  y = c(8.25, 8.25, 8.25, 6.55, 6.55, 6.55, 6.55, 4.65, 4.65, 4.65, 2.65),
  xend = c(1.75, 5, 7.45, 7.45, 1.75, 5, 7.45, 5, 5, 5, 5),
  yend = c(7.65, 7.65, 7.65, 5.75, 5.75, 5.75, 5.75, 3.75, 3.75, 3.75, 1.75)
)

box_colors <- c(
  "start" = "#172033",
  "branch" = "#0072B2",
  "control" = "#CC79A7",
  "measure" = "#009E73",
  "decision" = "#E69F00",
  "finish" = "#D55E00"
)

p_workflow <- ggplot() +
  geom_segment(
    data = arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "#8B94A3",
    linewidth = 0.6,
    arrow = grid::arrow(
      length = grid::unit(2.1, "mm"),
      type = "closed"
    )
  ) +
  geom_rect(
    data = boxes,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = type
    ),
    color = "white",
    linewidth = 0.8,
    alpha = 0.95
  ) +
  geom_text(
    data = boxes,
    aes(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = label,
      size = text_size
    ),
    color = "white",
    lineheight = 0.93,
    fontface = "bold"
  ) +
  scale_fill_manual(values = box_colors) +
  scale_size_identity() +
  coord_cartesian(
    xlim = c(0.1, 9.1),
    ylim = c(0.35, 9.6),
    clip = "off"
  ) +
  labs(
    title = "A depletion protocol earns its place in a paired pilot",
    subtitle = paste(
      "Efficiency, target retention, controls, and privacy",
      "are separate acceptance gates"
    ),
    caption = paste(
      "Use sample- and endpoint-specific failure rules;",
      "there is no universal host-fraction cutoff."
    )
  ) +
  theme_void(base_size = 11.5, base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 14,
      color = "#172033"
    ),
    plot.subtitle = element_text(size = 10.5, color = "#3D4757"),
    plot.caption = element_text(
      size = 8.5,
      color = "#596273",
      hjust = 0
    ),
    plot.margin = margin(8, 10, 8, 10)
  )

save_pub(
  p_workflow,
  file.path(figure_dir, "06-host-depletion-decision"),
  width = 198,
  height = 145
)

message("Article 06 figures written to ", normalizePath(figure_dir))
