#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999, digits = 17)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260729L
set.seed(primary_seed)

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
required <- c("project-root", "input-dir", "output-dir", "figure-dir", "chapter")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c(
  "huge", "igraph", "ggplot2", "patchwork", "ggrepel",
  "scales", "jsonlite", "digest"
)
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
input_dir <- normalizePath(args[["input-dir"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
chapter_path <- normalizePath(args[["chapter"]], mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

font_cache_dir <- file.path(tempdir(), "article29-fontconfig-cache")
dir.create(font_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

validation_log <- file.path(output_dir, "validation.log")
writeLines(
  c(
    "Article 29 conditional-association network validation",
    paste0("StartedUTC\t", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Seed\t", primary_seed)
  ),
  validation_log
)
log_msg <- function(...) {
  line <- paste0(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    "\t", paste0(..., collapse = "")
  )
  cat(line, "\n")
  cat(line, "\n", file = validation_log, append = TRUE)
  invisible(line)
}

checks <- data.frame(
  Category = character(), CheckID = character(), Status = character(),
  Detail = character(), stringsAsFactors = FALSE
)
add_check <- function(category, check_id, passed, detail) {
  checks <<- rbind(
    checks,
    data.frame(
      Category = category,
      CheckID = check_id,
      Status = if (isTRUE(passed)) "PASS" else "FAIL",
      Detail = paste(detail, collapse = "; "),
      stringsAsFactors = FALSE
    )
  )
  invisible(passed)
}
sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x, con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
read_tsv <- function(name) {
  utils::read.delim(
    file.path(input_dir, name), check.names = FALSE,
    quote = "", comment.char = "", stringsAsFactors = FALSE
  )
}
read_tsv_gz <- function(name) {
  utils::read.delim(
    gzfile(file.path(input_dir, name)), check.names = FALSE,
    quote = "", comment.char = "", stringsAsFactors = FALSE
  )
}

verify_checksum_manifest <- function(directory) {
  manifest_path <- file.path(directory, "file-checksums.sha256")
  lines <- readLines(manifest_path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  expected_files <- character()
  for (line in lines) {
    pieces <- strsplit(line, "[[:space:]]+", perl = TRUE)[[1L]]
    expected <- pieces[[1L]]
    relative <- paste(pieces[-1L], collapse = " ")
    path <- file.path(directory, relative)
    observed <- if (file.exists(path)) sha256_file(path) else "missing"
    add_check(
      "Frozen input", paste0("sha256-", relative),
      identical(observed, expected), observed
    )
    expected_files <- c(expected_files, relative)
  }
  payloads <- sort(setdiff(
    basename(list.files(directory, full.names = TRUE)),
    "file-checksums.sha256"
  ))
  add_check(
    "Frozen input", "checksum-manifest-complete",
    identical(payloads, sort(expected_files)),
    paste0("payloads=", length(payloads), "; entries=", length(expected_files))
  )
  length(lines)
}

checksum_entries <- verify_checksum_manifest(input_dir)
composition <- read_tsv_gz("spring-mag-relative-abundance.tsv.gz")
metadata <- read_tsv("spring-metadata.tsv")
filter_audit <- read_tsv("feature-filter-audit.tsv")
analysis_contract <- read_tsv("analysis-contract.tsv")
resource_manifest <- read_tsv("resource-manifest.tsv")

stopifnot(identical(names(composition)[1:2], c("MAG", "Taxonomy")))
spring_ids <- names(composition)[-(1:2)]
abundance <- t(as.matrix(composition[, spring_ids, drop = FALSE]))
storage.mode(abundance) <- "double"
rownames(abundance) <- spring_ids
colnames(abundance) <- composition$MAG
metadata <- metadata[match(spring_ids, metadata$Spring), , drop = FALSE]
metadata$BroadRegion <- factor(metadata$BroadRegion, levels = sort(unique(metadata$BroadRegion)))

add_check("Data", "composition-dimensions", identical(dim(abundance), c(56L, 780L)), paste(dim(abundance), collapse = "x"))
add_check("Data", "metadata-dimensions", nrow(metadata) == 56L, nrow(metadata))
add_check("Data", "unit-order", identical(rownames(abundance), metadata$Spring), "spring IDs aligned")
add_check("Data", "feature-ids-unique", !anyDuplicated(colnames(abundance)), ncol(abundance))
add_check("Data", "finite-nonnegative", all(is.finite(abundance)) && min(abundance) >= 0, range(abundance))
add_check("Data", "catalog-closure", max(abs(rowSums(abundance) - 1)) < 1e-12, max(abs(rowSums(abundance) - 1)))
add_check("Data", "environment-complete", !anyNA(metadata[, c("BroadRegion", "MedianTemperatureC", "MedianPH")]), "region / temperature / pH")
add_check("Data", "six-broad-regions", nlevels(metadata$BroadRegion) == 6L, levels(metadata$BroadRegion))
add_check("Data", "filter-audit-rows", nrow(filter_audit) == 780L, nrow(filter_audit))
add_check("Data", "contract-rows", nrow(analysis_contract) == 27L, nrow(analysis_contract))
add_check("Data", "resource-manifest", nrow(resource_manifest) == 6L, nrow(resource_manifest))

prevalence <- colMeans(abundance > 0)
mean_abundance <- colMeans(abundance)
relaxed_features <- colnames(abundance)[prevalence >= 0.60 & mean_abundance >= 0.001]
primary_features <- colnames(abundance)[prevalence >= 0.70 & mean_abundance >= 0.001]
strict_features <- colnames(abundance)[prevalence >= 0.70 & mean_abundance >= 0.002]
add_check("Filter", "relaxed-feature-count", length(relaxed_features) == 183L, length(relaxed_features))
add_check("Filter", "primary-feature-count", length(primary_features) == 93L, length(primary_features))
add_check("Filter", "strict-feature-count", length(strict_features) == 63L, length(strict_features))
add_check(
  "Filter", "primary-filter-audit-match",
  identical(sort(primary_features), sort(filter_audit$MAG[filter_audit$PrimaryFilter])),
  "computed against frozen audit"
)

taxonomy_rank <- function(taxonomy, rank_index, fallback = "Unclassified") {
  vapply(strsplit(taxonomy, ";", fixed = TRUE), function(x) {
    if (length(x) >= rank_index && nzchar(x[[rank_index]])) x[[rank_index]] else fallback
  }, character(1L))
}

transform_network_data <- function(x, meta, pseudocount = 1e-6, adjust = TRUE) {
  if (any(!is.finite(x)) || min(x) < 0 || any(rowSums(x) <= 0)) {
    stop("Network input contains invalid abundance rows.", call. = FALSE)
  }
  x <- sweep(x, 1L, rowSums(x), "/")
  logged <- log(x + pseudocount)
  clr <- logged - rowMeans(logged)
  design_rank <- NA_integer_
  residual_df <- nrow(clr) - 1L
  if (isTRUE(adjust)) {
    meta$BroadRegion <- factor(meta$BroadRegion, levels = levels(metadata$BroadRegion))
    design <- stats::model.matrix(
      ~ BroadRegion + scale(MedianTemperatureC) + scale(MedianPH),
      data = meta
    )
    design_qr <- qr(design)
    design_rank <- design_qr$rank
    if (design_rank != ncol(design)) {
      stop("Environmental design matrix is rank deficient.", call. = FALSE)
    }
    transformed <- qr.resid(design_qr, clr)
    residual_df <- nrow(clr) - design_rank
  } else {
    transformed <- clr
  }
  standard_deviation <- apply(transformed, 2L, stats::sd)
  if (any(!is.finite(standard_deviation)) || any(standard_deviation <= 0)) {
    stop("At least one transformed feature has zero variance.", call. = FALSE)
  }
  transformed <- scale(transformed, center = TRUE, scale = TRUE)
  storage.mode(transformed) <- "double"
  attr(transformed, "design_rank") <- design_rank
  attr(transformed, "residual_df") <- residual_df
  transformed
}

fit_stars_network <- function(
  x, meta, pseudocount = 1e-6, adjust = TRUE,
  method = "glasso", stars_threshold = 0.05,
  rep_num = 100L, seed = primary_seed
) {
  transformed <- transform_network_data(x, meta, pseudocount, adjust)
  set.seed(seed)
  path_fit <- huge::huge(
    transformed,
    method = method,
    nlambda = 30L,
    lambda.min.ratio = 0.05,
    scr = FALSE,
    sym = "or",
    verbose = FALSE
  )
  selected <- huge::huge.select(
    path_fit,
    criterion = "stars",
    stars.thresh = stars_threshold,
    stars.subsample.ratio = 0.80,
    rep.num = rep_num,
    verbose = FALSE
  )
  adjacency <- as.matrix(selected$refit != 0)
  adjacency <- adjacency | t(adjacency)
  diag(adjacency) <- FALSE
  dimnames(adjacency) <- list(colnames(transformed), colnames(transformed))
  partial <- matrix(
    NA_real_, ncol(transformed), ncol(transformed),
    dimnames = dimnames(adjacency)
  )
  if (method == "glasso") {
    precision <- as.matrix(selected$opt.icov)
    partial <- -precision / outer(sqrt(diag(precision)), sqrt(diag(precision)))
    diag(partial) <- 0
    dimnames(partial) <- dimnames(adjacency)
  }
  list(
    adjacency = adjacency,
    partial = partial,
    lambda = selected$opt.lambda,
    opt_index = selected$opt.index,
    sparsity = selected$opt.sparsity,
    variability = selected$variability,
    selected_variability = selected$variability[[selected$opt.index]],
    features = colnames(transformed),
    transformed = transformed,
    method = method,
    pseudocount = pseudocount,
    adjust = adjust,
    stars_threshold = stars_threshold,
    rep_num = rep_num
  )
}

log_msg("Fitting primary 93-MAG CLR-glasso-StARS network")
primary_x <- abundance[, primary_features, drop = FALSE]
primary <- fit_stars_network(
  primary_x, metadata,
  pseudocount = 1e-6,
  adjust = TRUE,
  method = "glasso",
  stars_threshold = 0.05,
  rep_num = 100L,
  seed = primary_seed
)
primary_edges_n <- sum(primary$adjacency[upper.tri(primary$adjacency)])
primary_nodes_n <- sum(rowSums(primary$adjacency) > 0)
add_check("Primary model", "design-rank", attr(primary$transformed, "design_rank") == 8L, attr(primary$transformed, "design_rank"))
add_check("Primary model", "residual-df", attr(primary$transformed, "residual_df") == 48L, attr(primary$transformed, "residual_df"))
add_check("Primary model", "selected-lambda-finite", is.finite(primary$lambda) && primary$lambda > 0, primary$lambda)
add_check("Primary model", "selected-index", primary$opt_index >= 1L && primary$opt_index <= 30L, primary$opt_index)
add_check("Primary model", "sparse-edge-count", primary_edges_n > 0L && primary_edges_n < choose(93L, 2L) * 0.10, primary_edges_n)
add_check("Primary model", "nonisolates", primary_nodes_n >= 20L && primary_nodes_n <= 93L, primary_nodes_n)
add_check("Primary model", "adjacency-symmetric", identical(primary$adjacency, t(primary$adjacency)), "logical symmetric matrix")
add_check("Primary model", "partial-finite-on-edges", all(is.finite(primary$partial[primary$adjacency])), range(primary$partial[primary$adjacency]))

log_msg("Running 1000 BroadRegion-stratified fixed-lambda bootstraps")
bootstrap_replicates <- 1000L
p <- length(primary_features)
selection_count <- matrix(0L, p, p, dimnames = dimnames(primary$adjacency))
positive_count <- matrix(0L, p, p, dimnames = dimnames(primary$adjacency))
negative_count <- matrix(0L, p, p, dimnames = dimnames(primary$adjacency))
partial_sum <- matrix(0, p, p, dimnames = dimnames(primary$adjacency))
partial_squared_sum <- matrix(0, p, p, dimnames = dimnames(primary$adjacency))
region_indices <- split(seq_len(nrow(metadata)), metadata$BroadRegion)
set.seed(primary_seed + 1000L)
completed_bootstraps <- 0L
bootstrap_attempts <- 0L
while (completed_bootstraps < bootstrap_replicates) {
  bootstrap_attempts <- bootstrap_attempts + 1L
  if (bootstrap_attempts > bootstrap_replicates + 100L) {
    stop("Too many rank-deficient bootstrap attempts.", call. = FALSE)
  }
  index <- unlist(lapply(
    region_indices,
    function(z) sample(z, length(z), replace = TRUE)
  ), use.names = FALSE)
  transformed <- try(
    transform_network_data(
      primary_x[index, , drop = FALSE],
      metadata[index, , drop = FALSE],
      pseudocount = 1e-6,
      adjust = TRUE
    ),
    silent = TRUE
  )
  if (inherits(transformed, "try-error")) next
  fit <- huge::huge(
    transformed,
    method = "glasso",
    lambda = primary$lambda,
    scr = FALSE,
    verbose = FALSE
  )
  adjacency <- as.matrix(fit$path[[1L]] != 0)
  adjacency <- adjacency | t(adjacency)
  diag(adjacency) <- FALSE
  precision <- as.matrix(fit$icov[[1L]])
  partial <- -precision / outer(sqrt(diag(precision)), sqrt(diag(precision)))
  diag(partial) <- 0
  selected_partial <- ifelse(adjacency, partial, 0)
  selection_count <- selection_count + adjacency
  positive_count <- positive_count + (adjacency & partial > 0)
  negative_count <- negative_count + (adjacency & partial < 0)
  partial_sum <- partial_sum + selected_partial
  partial_squared_sum <- partial_squared_sum + selected_partial^2
  completed_bootstraps <- completed_bootstraps + 1L
}
bootstrap_frequency <- selection_count / bootstrap_replicates
bootstrap_ij <- which(upper.tri(bootstrap_frequency), arr.ind = TRUE)
edge_bootstrap <- data.frame(
  Source = rownames(bootstrap_frequency)[bootstrap_ij[, 1L]],
  Target = colnames(bootstrap_frequency)[bootstrap_ij[, 2L]],
  PrimaryEdge = primary$adjacency[bootstrap_ij],
  PrimaryPartialCorrelation = primary$partial[bootstrap_ij],
  SelectionFrequency = bootstrap_frequency[bootstrap_ij],
  PositiveFrequency = positive_count[bootstrap_ij] / bootstrap_replicates,
  NegativeFrequency = negative_count[bootstrap_ij] / bootstrap_replicates,
  MeanPartialAcrossAllBootstraps = partial_sum[bootstrap_ij] / bootstrap_replicates,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
edge_bootstrap$MeanPartialWhenSelected <- ifelse(
  selection_count[bootstrap_ij] > 0,
  partial_sum[bootstrap_ij] / selection_count[bootstrap_ij],
  NA_real_
)
edge_bootstrap$SDPartialAcrossAllBootstraps <- sqrt(pmax(
  partial_squared_sum[bootstrap_ij] / bootstrap_replicates -
    edge_bootstrap$MeanPartialAcrossAllBootstraps^2,
  0
))
edge_bootstrap$PrimarySignAgreement <- ifelse(
  edge_bootstrap$PrimaryPartialCorrelation > 0,
  ifelse(selection_count[bootstrap_ij] > 0,
         positive_count[bootstrap_ij] / selection_count[bootstrap_ij], NA_real_),
  ifelse(selection_count[bootstrap_ij] > 0,
         negative_count[bootstrap_ij] / selection_count[bootstrap_ij], NA_real_)
)
edge_bootstrap$Consensus70 <- edge_bootstrap$SelectionFrequency >= 0.70
edge_bootstrap$HighStability80 <- edge_bootstrap$SelectionFrequency >= 0.80
edge_bootstrap$EdgeKey <- paste(edge_bootstrap$Source, edge_bootstrap$Target, sep = "||")

primary_edge_table <- edge_bootstrap[edge_bootstrap$PrimaryEdge, , drop = FALSE]
primary_edge_table$Sign <- ifelse(primary_edge_table$PrimaryPartialCorrelation > 0, "Positive", "Negative")
primary_edge_table$AbsolutePartialCorrelation <- abs(primary_edge_table$PrimaryPartialCorrelation)
primary_edge_table <- primary_edge_table[order(
  -primary_edge_table$SelectionFrequency,
  -primary_edge_table$AbsolutePartialCorrelation,
  primary_edge_table$Source,
  primary_edge_table$Target
), ]
rownames(primary_edge_table) <- NULL
add_check("Bootstrap", "replicates", completed_bootstraps == 1000L, completed_bootstraps)
add_check("Bootstrap", "attempts-bounded", bootstrap_attempts <= 1100L, bootstrap_attempts)
add_check("Bootstrap", "primary-edge-rows", nrow(primary_edge_table) == primary_edges_n, nrow(primary_edge_table))
add_check("Bootstrap", "frequency-range", all(edge_bootstrap$SelectionFrequency >= 0 & edge_bootstrap$SelectionFrequency <= 1), range(edge_bootstrap$SelectionFrequency))
add_check("Bootstrap", "consensus-subset", sum(primary_edge_table$Consensus70) > 0L, sum(primary_edge_table$Consensus70))
add_check("Bootstrap", "high-stability-subset", sum(primary_edge_table$HighStability80) > 0L, sum(primary_edge_table$HighStability80))

log_msg("Computing modules, Zi-Pi roles and topology-priority candidates")
graph_all <- igraph::graph_from_adjacency_matrix(
  primary$adjacency, mode = "undirected", diag = FALSE
)
core_names <- igraph::V(graph_all)$name[igraph::degree(graph_all) > 0]
graph_core <- igraph::induced_subgraph(graph_all, vids = core_names)
core_edge_ends <- igraph::as_edgelist(graph_core, names = TRUE)
edge_lookup <- paste(
  pmin(core_edge_ends[, 1L], core_edge_ends[, 2L]),
  pmax(core_edge_ends[, 1L], core_edge_ends[, 2L]),
  sep = "||"
)
primary_lookup <- match(edge_lookup, primary_edge_table$EdgeKey)
stopifnot(!anyNA(primary_lookup))
graph_core <- igraph::set_edge_attr(
  graph_core, "PartialCorrelation",
  value = primary_edge_table$PrimaryPartialCorrelation[primary_lookup]
)
graph_core <- igraph::set_edge_attr(
  graph_core, "BootstrapFrequency",
  value = primary_edge_table$SelectionFrequency[primary_lookup]
)

set.seed(primary_seed + 2000L)
community <- igraph::cluster_louvain(
  graph_core, weights = rep(1, igraph::ecount(graph_core))
)
raw_membership <- igraph::membership(community)
module_members <- split(names(raw_membership), raw_membership)
module_order <- order(
  -vapply(module_members, length, integer(1L)),
  vapply(module_members, function(x) sort(x)[[1L]], character(1L))
)
module_map <- stats::setNames(
  paste0("M", seq_along(module_order)),
  names(module_members)[module_order]
)
membership <- module_map[as.character(raw_membership)]
names(membership) <- names(raw_membership)

degree_value <- igraph::degree(graph_core)
strength_value <- igraph::strength(
  graph_core,
  weights = abs(igraph::edge_attr(graph_core, "PartialCorrelation"))
)
stable_degree70 <- rowSums(bootstrap_frequency[core_names, , drop = FALSE] >= 0.70)
stable_degree80 <- rowSums(bootstrap_frequency[core_names, , drop = FALSE] >= 0.80)
expected_degree <- rowSums(bootstrap_frequency[core_names, , drop = FALSE])
mean_incident_stability <- vapply(core_names, function(node) {
  edge_ids <- igraph::incident(graph_core, node, mode = "all")
  mean(igraph::edge_attr(graph_core, "BootstrapFrequency", index = edge_ids))
}, numeric(1L))

within_degree <- numeric(length(core_names))
participation <- numeric(length(core_names))
names(within_degree) <- names(participation) <- core_names
for (node in core_names) {
  neighbor_names <- igraph::neighbors(graph_core, node, mode = "all")$name
  neighbor_modules <- membership[neighbor_names]
  node_degree <- length(neighbor_names)
  within_degree[[node]] <- sum(neighbor_modules == membership[[node]])
  module_degree <- table(neighbor_modules)
  participation[[node]] <- 1 - sum(as.numeric(module_degree)^2) / node_degree^2
}
zi <- numeric(length(core_names))
zi_estimable <- logical(length(core_names))
names(zi) <- names(zi_estimable) <- core_names
for (module in unique(membership)) {
  nodes <- names(membership)[membership == module]
  values <- within_degree[nodes]
  estimable <- length(nodes) >= 3L && is.finite(stats::sd(values)) && stats::sd(values) > 0
  if (estimable) {
    zi[nodes] <- (values - mean(values)) / stats::sd(values)
    zi_estimable[nodes] <- TRUE
  } else {
    zi[nodes] <- 0
    zi_estimable[nodes] <- FALSE
  }
}

role <- ifelse(
  zi > 2.5,
  ifelse(
    participation <= 0.30, "Provincial hub",
    ifelse(participation <= 0.75, "Connector hub", "Kinless hub")
  ),
  ifelse(
    participation <= 0.05, "Ultra-peripheral",
    ifelse(
      participation <= 0.62, "Peripheral",
      ifelse(participation <= 0.80, "Connector", "Kinless")
    )
  )
)

node_taxonomy <- filter_audit[match(core_names, filter_audit$MAG), , drop = FALSE]
node_roles <- data.frame(
  MAG = core_names,
  Taxonomy = node_taxonomy$Taxonomy,
  Phylum = taxonomy_rank(node_taxonomy$Taxonomy, 2L),
  Class = taxonomy_rank(node_taxonomy$Taxonomy, 3L),
  Genus = taxonomy_rank(node_taxonomy$Taxonomy, 6L),
  Species = taxonomy_rank(node_taxonomy$Taxonomy, 7L),
  Module = unname(membership[core_names]),
  Degree = as.numeric(degree_value[core_names]),
  Strength = as.numeric(strength_value[core_names]),
  StableDegree70 = as.numeric(stable_degree70[core_names]),
  StableDegree80 = as.numeric(stable_degree80[core_names]),
  BootstrapExpectedDegree = as.numeric(expected_degree[core_names]),
  MeanIncidentStability = as.numeric(mean_incident_stability[core_names]),
  WithinModuleDegree = as.numeric(within_degree[core_names]),
  Zi = as.numeric(zi[core_names]),
  ZiEstimable = as.logical(zi_estimable[core_names]),
  Pi = as.numeric(participation[core_names]),
  Role = role,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
distance_weight <- 1 / pmax(abs(igraph::edge_attr(graph_core, "PartialCorrelation")), 1e-6)
node_roles$Betweenness <- as.numeric(igraph::betweenness(
  graph_core, directed = FALSE, weights = distance_weight, normalized = TRUE
)[node_roles$MAG])
node_roles$HarmonicCentrality <- as.numeric(igraph::harmonic_centrality(
  graph_core, mode = "all", weights = distance_weight, normalized = TRUE
)[node_roles$MAG])
degree_cutoff <- as.numeric(stats::quantile(node_roles$Degree, 0.90, type = 1))
node_roles$TopologyPriority <- node_roles$StableDegree70 >= 2L & (
  node_roles$Zi > 2.5 |
    node_roles$Pi > 0.62 |
    node_roles$Degree >= degree_cutoff
)
node_roles$CandidateEvidence <- vapply(seq_len(nrow(node_roles)), function(i) {
  evidence <- character()
  if (node_roles$Zi[[i]] > 2.5) evidence <- c(evidence, "within-module hub")
  if (node_roles$Pi[[i]] > 0.62) evidence <- c(evidence, "cross-module connector")
  if (node_roles$Degree[[i]] >= degree_cutoff) evidence <- c(evidence, "top-decile degree")
  if (node_roles$StableDegree70[[i]] >= 2L) evidence <- c(evidence, ">=2 consensus edges")
  if (length(evidence) == 0L) "none" else paste(evidence, collapse = "; ")
}, character(1L))
node_roles <- node_roles[order(
  -node_roles$TopologyPriority,
  -node_roles$StableDegree70,
  -node_roles$BootstrapExpectedDegree,
  -node_roles$Degree,
  node_roles$MAG
), ]
rownames(node_roles) <- NULL
topology_candidates <- node_roles[node_roles$TopologyPriority, , drop = FALSE]
topology_candidates$Interpretation <- "Topology-priority candidate; ecological keystone status requires perturbation or longitudinal validation"

component_sizes <- sort(igraph::components(graph_core)$csize, decreasing = TRUE)
observed_modularity <- igraph::modularity(
  graph_core, as.integer(factor(membership, levels = unique(membership))),
  weights = rep(1, igraph::ecount(graph_core))
)
observed_transitivity <- igraph::transitivity(graph_core, type = "global")
observed_assortativity <- igraph::assortativity_degree(graph_core, directed = FALSE)
add_check("Topology", "core-edge-count", igraph::ecount(graph_core) == primary_edges_n, igraph::ecount(graph_core))
add_check("Topology", "core-node-count", igraph::vcount(graph_core) == primary_nodes_n, igraph::vcount(graph_core))
add_check("Topology", "component-accounting", sum(component_sizes) == primary_nodes_n, paste(component_sizes, collapse = "/"))
add_check("Topology", "module-count", length(unique(node_roles$Module)) >= 2L, length(unique(node_roles$Module)))
add_check("Topology", "zi-finite", all(is.finite(node_roles$Zi)), range(node_roles$Zi))
add_check("Topology", "pi-bounded", all(node_roles$Pi >= -1e-12 & node_roles$Pi <= 1 + 1e-12), range(node_roles$Pi))
add_check("Topology", "topology-priority-candidates", nrow(topology_candidates) > 0L, nrow(topology_candidates))

log_msg("Running 1000 degree-preserving topology nulls")
null_replicates <- 1000L
topology_null <- data.frame(
  Replicate = seq_len(null_replicates),
  Modularity = NA_real_,
  Transitivity = NA_real_,
  DegreePreserved = FALSE,
  stringsAsFactors = FALSE
)
original_degree_sorted <- sort(as.numeric(igraph::degree(graph_core)))
for (i in seq_len(null_replicates)) {
  set.seed(primary_seed + 10000L + i)
  null_graph <- igraph::rewire(
    graph_core,
    with = igraph::keeping_degseq(
      niter = 20L * igraph::ecount(graph_core), loops = FALSE
    )
  )
  set.seed(primary_seed + 20000L + i)
  null_community <- igraph::cluster_louvain(
    null_graph, weights = rep(1, igraph::ecount(null_graph))
  )
  topology_null$Modularity[[i]] <- igraph::modularity(
    null_graph, igraph::membership(null_community),
    weights = rep(1, igraph::ecount(null_graph))
  )
  topology_null$Transitivity[[i]] <- igraph::transitivity(null_graph, type = "global")
  topology_null$DegreePreserved[[i]] <- identical(
    sort(as.numeric(igraph::degree(null_graph))), original_degree_sorted
  )
}
topology_null_summary <- data.frame(
  Metric = c("Modularity", "Transitivity"),
  Observed = c(observed_modularity, observed_transitivity),
  NullMean = c(mean(topology_null$Modularity), mean(topology_null$Transitivity)),
  NullSD = c(stats::sd(topology_null$Modularity), stats::sd(topology_null$Transitivity)),
  NullLower95 = c(
    stats::quantile(topology_null$Modularity, 0.025, type = 6),
    stats::quantile(topology_null$Transitivity, 0.025, type = 6)
  ),
  NullUpper95 = c(
    stats::quantile(topology_null$Modularity, 0.975, type = 6),
    stats::quantile(topology_null$Transitivity, 0.975, type = 6)
  ),
  EmpiricalUpperP = c(
    (1 + sum(topology_null$Modularity >= observed_modularity)) / (null_replicates + 1),
    (1 + sum(topology_null$Transitivity >= observed_transitivity)) / (null_replicates + 1)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
add_check("Topology null", "degree-preserved", all(topology_null$DegreePreserved), sum(topology_null$DegreePreserved))
add_check("Topology null", "finite-null-metrics", all(is.finite(topology_null$Modularity)) && all(is.finite(topology_null$Transitivity)), "1000 x 2 metrics")

lcc_size <- function(graph) {
  if (igraph::vcount(graph) == 0L) return(0L)
  max(igraph::components(graph)$csize)
}
trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

log_msg("Comparing adaptive-degree deletion with 1000 random deletion orders")
initial_nodes <- igraph::vcount(graph_core)
removal_steps <- 0:(initial_nodes - 1L)
target_graph <- graph_core
target_curve <- numeric(initial_nodes)
target_removed <- character(initial_nodes)
for (k in removal_steps) {
  target_curve[[k + 1L]] <- lcc_size(target_graph) / initial_nodes
  if (k < initial_nodes - 1L) {
    candidates <- igraph::V(target_graph)$name
    current_degree <- igraph::degree(target_graph)
    tie_expected <- stats::setNames(node_roles$BootstrapExpectedDegree, node_roles$MAG)[candidates]
    remove_order <- order(-current_degree[candidates], -tie_expected, candidates)
    remove_node <- candidates[remove_order[[1L]]]
    target_removed[[k + 1L]] <- remove_node
    target_graph <- igraph::delete_vertices(target_graph, remove_node)
  }
}
random_replicates <- 1000L
random_curve_matrix <- matrix(
  NA_real_, nrow = random_replicates, ncol = initial_nodes
)
random_auc <- numeric(random_replicates)
set.seed(primary_seed + 30000L)
for (i in seq_len(random_replicates)) {
  deletion_order <- sample(igraph::V(graph_core)$name, initial_nodes, replace = FALSE)
  random_graph <- graph_core
  for (k in removal_steps) {
    random_curve_matrix[i, k + 1L] <- lcc_size(random_graph) / initial_nodes
    if (k < initial_nodes - 1L) {
      random_graph <- igraph::delete_vertices(random_graph, deletion_order[[k + 1L]])
    }
  }
  random_auc[[i]] <- trapz(removal_steps / initial_nodes, random_curve_matrix[i, ])
}
target_auc <- trapz(removal_steps / initial_nodes, target_curve)
random_curve_summary <- data.frame(
  RemovedNodes = removal_steps,
  RemovedFraction = removal_steps / initial_nodes,
  RandomMedian = apply(random_curve_matrix, 2L, stats::median),
  RandomLower95 = apply(random_curve_matrix, 2L, stats::quantile, probs = 0.025, type = 6),
  RandomUpper95 = apply(random_curve_matrix, 2L, stats::quantile, probs = 0.975, type = 6),
  TargetedAdaptiveDegree = target_curve,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
robustness_summary <- data.frame(
  Nodes = initial_nodes,
  InitialLargestComponent = lcc_size(graph_core),
  InitialLargestComponentFraction = lcc_size(graph_core) / initial_nodes,
  TargetedAUC = target_auc,
  RandomMedianAUC = stats::median(random_auc),
  RandomLower95AUC = stats::quantile(random_auc, 0.025, type = 6),
  RandomUpper95AUC = stats::quantile(random_auc, 0.975, type = 6),
  RandomOrders = random_replicates,
  EmpiricalLowerP = (1 + sum(random_auc <= target_auc)) / (random_replicates + 1),
  Normalization = "largest connected component / initial non-isolate nodes",
  Interpretation = "topological deletion diagnostic; not ecological resilience",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
targeted_deletion_order <- data.frame(
  Step = seq_len(initial_nodes - 1L),
  RemovedMAG = target_removed[seq_len(initial_nodes - 1L)],
  stringsAsFactors = FALSE
)
add_check("Robustness", "random-orders", length(random_auc) == 1000L, length(random_auc))
add_check("Robustness", "curve-bounds", all(random_curve_matrix >= 0 & random_curve_matrix <= 1) && all(target_curve >= 0 & target_curve <= 1), range(c(random_curve_matrix, target_curve)))
add_check("Robustness", "targeted-more-damaging", target_auc < stats::median(random_auc), paste(target_auc, stats::median(random_auc), sep = "/"))

edge_keys_from_adjacency <- function(adjacency) {
  ij <- which(adjacency & upper.tri(adjacency), arr.ind = TRUE)
  if (nrow(ij) == 0L) return(character())
  paste(
    rownames(adjacency)[ij[, 1L]],
    colnames(adjacency)[ij[, 2L]],
    sep = "||"
  )
}
primary_edge_keys <- edge_keys_from_adjacency(primary$adjacency)

branch_specs <- list(
  list(Name = "Pseudocount 1e-7", Features = primary_features, Pseudocount = 1e-7, Adjust = TRUE, Method = "glasso", Threshold = 0.05),
  list(Name = "Pseudocount 1e-5", Features = primary_features, Pseudocount = 1e-5, Adjust = TRUE, Method = "glasso", Threshold = 0.05),
  list(Name = "Unadjusted CLR", Features = primary_features, Pseudocount = 1e-6, Adjust = FALSE, Method = "glasso", Threshold = 0.05),
  list(Name = "MB neighborhood (OR)", Features = primary_features, Pseudocount = 1e-6, Adjust = TRUE, Method = "mb", Threshold = 0.05),
  list(Name = "Relaxed filter", Features = relaxed_features, Pseudocount = 1e-6, Adjust = TRUE, Method = "glasso", Threshold = 0.05),
  list(Name = "Strict filter", Features = strict_features, Pseudocount = 1e-6, Adjust = TRUE, Method = "glasso", Threshold = 0.05),
  list(Name = "StARS 0.025", Features = primary_features, Pseudocount = 1e-6, Adjust = TRUE, Method = "glasso", Threshold = 0.025),
  list(Name = "StARS 0.10", Features = primary_features, Pseudocount = 1e-6, Adjust = TRUE, Method = "glasso", Threshold = 0.10)
)

primary_sensitivity_row <- data.frame(
  Branch = "Primary",
  Features = length(primary_features),
  Method = "glasso",
  EnvironmentAdjusted = TRUE,
  Pseudocount = 1e-6,
  StARSThreshold = 0.05,
  StARSReplicates = 100L,
  SelectedLambda = primary$lambda,
  SelectedPathIndex = primary$opt_index,
  Edges = primary_edges_n,
  Nonisolates = primary_nodes_n,
  SharedFeaturesWithPrimary = length(primary_features),
  PrimaryEligibleEdges = primary_edges_n,
  EdgeIntersection = primary_edges_n,
  EdgeUnion = primary_edges_n,
  EdgeJaccardWithPrimary = 1,
  PrimaryEdgeRetention = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
sensitivity_rows <- list(primary_sensitivity_row)
log_msg("Running eight preregistered network sensitivity branches")
for (i in seq_along(branch_specs)) {
  spec <- branch_specs[[i]]
  log_msg("Sensitivity branch: ", spec$Name)
  branch_x <- abundance[, spec$Features, drop = FALSE]
  branch <- fit_stars_network(
    branch_x, metadata,
    pseudocount = spec$Pseudocount,
    adjust = spec$Adjust,
    method = spec$Method,
    stars_threshold = spec$Threshold,
    rep_num = 50L,
    seed = primary_seed + 40000L + i
  )
  branch_keys <- edge_keys_from_adjacency(branch$adjacency)
  shared_features <- intersect(primary_features, spec$Features)
  primary_eligible <- primary_edge_keys[vapply(
    strsplit(primary_edge_keys, "\\|\\|"),
    function(z) all(z %in% shared_features), logical(1L)
  )]
  branch_eligible <- branch_keys[vapply(
    strsplit(branch_keys, "\\|\\|"),
    function(z) all(z %in% shared_features), logical(1L)
  )]
  intersection <- intersect(primary_eligible, branch_eligible)
  union <- union(primary_eligible, branch_eligible)
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
    Branch = spec$Name,
    Features = length(spec$Features),
    Method = spec$Method,
    EnvironmentAdjusted = spec$Adjust,
    Pseudocount = spec$Pseudocount,
    StARSThreshold = spec$Threshold,
    StARSReplicates = 50L,
    SelectedLambda = branch$lambda,
    SelectedPathIndex = branch$opt_index,
    Edges = length(branch_keys),
    Nonisolates = sum(rowSums(branch$adjacency) > 0),
    SharedFeaturesWithPrimary = length(shared_features),
    PrimaryEligibleEdges = length(primary_eligible),
    EdgeIntersection = length(intersection),
    EdgeUnion = length(union),
    EdgeJaccardWithPrimary = if (length(union) > 0L) length(intersection) / length(union) else NA_real_,
    PrimaryEdgeRetention = if (length(primary_eligible) > 0L) length(intersection) / length(primary_eligible) else NA_real_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
network_sensitivity <- do.call(rbind, sensitivity_rows)
add_check("Sensitivity", "branch-count", nrow(network_sensitivity) == 9L, nrow(network_sensitivity))
add_check("Sensitivity", "finite-jaccard", all(is.finite(network_sensitivity$EdgeJaccardWithPrimary)), network_sensitivity$EdgeJaccardWithPrimary)
add_check("Sensitivity", "all-branches-nonempty", all(network_sensitivity$Edges > 0L), network_sensitivity$Edges)

data_lineage <- data.frame(
  Stage = c(
    "Public study", "Frozen source", "Inference unit", "Feature universe",
    "Primary transform", "Primary graph", "Uncertainty", "Interpretation"
  ),
  ResourceOrOperation = c(
    "Korchagina et al. 2026; Figshare 30284068 v2",
    "500 sample-relative MAG profiles with checksum-locked metadata",
    "equal sample-relative mean within each hot spring, then closure",
    "780 catalog MAGs to 93 by outcome-free prevalence and abundance filter",
    "reclosure + 1e-6 zero replacement + CLR + region/temperature/pH residuals",
    "huge 1.5.1 graphical lasso with StARS 0.05",
    "1000 region-stratified bootstraps + sensitivity + degree-preserving null",
    "conditional association and topology hypotheses only"
  ),
  Units = c("500 metagenomes", "500 metagenomes", "56 hot springs", "56 hot springs", "56 x 93", "93 nodes", "resampling", "observational"),
  Boundary = c(
    "non-probability sampling designed to capture local diversity",
    "catalog-relative abundance; mean recruitment approximately 13%",
    "one spring is one inferential unit despite 1-33 local samples",
    "filter changes the network estimand",
    "linear adjustment cannot remove unmeasured habitat filtering",
    "sparse Gaussian graphical model assumptions apply",
    "bootstrap is conditional on chosen nodes and primary lambda",
    "no direct interaction, causality, direction, or ecological keystone claim"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

primary_network_summary <- data.frame(
  Springs = nrow(primary_x),
  CandidateMAGs = ncol(primary_x),
  ResidualDesignRank = attr(primary$transformed, "design_rank"),
  ResidualDF = attr(primary$transformed, "residual_df"),
  SelectedLambda = primary$lambda,
  SelectedPathIndex = primary$opt_index,
  SelectedStARSVariability = primary$selected_variability,
  Edges = primary_edges_n,
  PositiveEdges = sum(primary_edge_table$Sign == "Positive"),
  NegativeEdges = sum(primary_edge_table$Sign == "Negative"),
  Nonisolates = primary_nodes_n,
  Isolates = length(primary_features) - primary_nodes_n,
  ComponentsAmongAllCandidates = igraph::components(graph_all)$no,
  ComponentsAmongNonisolates = igraph::components(graph_core)$no,
  LargestComponent = max(component_sizes),
  Modules = length(unique(node_roles$Module)),
  Modularity = observed_modularity,
  Transitivity = observed_transitivity,
  DegreeAssortativity = observed_assortativity,
  ConsensusEdges70 = sum(primary_edge_table$Consensus70),
  HighStabilityEdges80 = sum(primary_edge_table$HighStability80),
  TopologyPriorityCandidates = nrow(topology_candidates),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

interpretation_boundaries <- data.frame(
  TemptingClaim = c(
    "An edge is a direct ecological interaction",
    "A positive edge proves cooperation",
    "A negative edge proves competition",
    "A high-degree or Zi-Pi node is a keystone",
    "Targeted deletion proves ecosystem fragility",
    "Environment adjustment removed habitat filtering",
    "A stable graph is a population-level truth"
  ),
  SupportedClaim = c(
    "The two CLR-residual features are conditionally associated under this sparse model",
    "The signed partial association is positive after measured linear adjustment",
    "The signed partial association is negative after measured linear adjustment",
    "The node is topology-priority and merits perturbation or longitudinal follow-up",
    "The inferred graph fragments faster under an adaptive degree attack than random deletion",
    "Linear effects of measured region, median temperature, and median pH were residualized",
    "Edges show resampling stability conditional on this dataset, feature gate, estimator, and lambda"
  ),
  MissingEvidence = c(
    "absolute abundance, temporal order, perturbation, spatial co-localization or culture",
    "metabolic exchange or reciprocal growth assay",
    "resource competition or inhibition experiment",
    "removal experiment showing disproportionate community change",
    "time-series recovery or perturbation-response measurements",
    "unmeasured chemistry, fine-scale geography, nonlinear effects and measurement error",
    "independent ecosystem replication and profiler/catalog sensitivity"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pal_pub <- c(
  blue = "#0072B2", orange = "#D55E00", green = "#009E73",
  purple = "#CC79A7", gold = "#E69F00", sky = "#56B4E9",
  brown = "#A6761D", grey = "#777777"
)
theme_pub <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.25),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey45"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(color = "grey35")
    )
}
save_pub <- function(plot, stem, width = 210, height = 140, dpi = 350) {
  base <- file.path(figure_dir, stem)
  ggplot2::ggsave(
    paste0(base, ".pdf"), plot,
    width = width, height = height, units = "mm",
    device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".png"), plot,
    width = width, height = height, units = "mm",
    dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".tiff"), plot,
    width = width, height = height, units = "mm",
    dpi = dpi, compression = "lzw", bg = "white"
  )
  invisible(plot)
}

filter_plot_data <- filter_audit
filter_plot_data$Status <- ifelse(filter_plot_data$PrimaryFilter, "Primary set", "Excluded")
filter_plot_data$Status <- factor(filter_plot_data$Status, levels = c("Excluded", "Primary set"))
p_filter_scatter <- ggplot2::ggplot(
  filter_plot_data,
  ggplot2::aes(Prevalence * 100, MeanRelativeAbundance * 100, color = Status)
) +
  ggplot2::geom_hline(yintercept = 0.1, linetype = 2, color = "grey45") +
  ggplot2::geom_vline(xintercept = 70, linetype = 2, color = "grey45") +
  ggplot2::geom_point(size = 1.55, alpha = 0.70) +
  ggplot2::scale_y_log10(
    labels = function(x) formatC(x, format = "fg", digits = 3),
    breaks = c(0.001, 0.01, 0.1, 1, 10)
  ) +
  ggplot2::scale_color_manual(values = c("Excluded" = "grey72", "Primary set" = pal_pub[["orange"]])) +
  ggplot2::annotate(
    "label", x = 98, y = max(filter_plot_data$MeanRelativeAbundance * 100) * 0.75,
    hjust = 1, label = "93 MAGs retained", size = 3.2,
    label.size = 0.2, fill = "white"
  ) +
  ggplot2::labs(
    title = "The feature gate is part of the network estimand",
    subtitle = "Each point is one recovered MAG across 56 equal-weight hot springs",
    x = "Hot-spring prevalence (%)",
    y = "Mean catalog-relative abundance (%)",
    color = NULL
  ) +
  theme_pub(9.5) +
  ggplot2::theme(legend.position = "bottom")

filter_counts <- data.frame(
  Gate = factor(
    c("Raw catalog", "Relaxed", "Primary", "Strict"),
    levels = c("Raw catalog", "Relaxed", "Primary", "Strict")
  ),
  MAGs = c(780L, 183L, 93L, 63L),
  Definition = c(
    "780 catalog MAGs", "Prev >=60%; mean >=0.1%",
    "Prev >=70%; mean >=0.1%", "Prev >=70%; mean >=0.2%"
  ),
  stringsAsFactors = FALSE
)
p_filter_counts <- ggplot2::ggplot(
  filter_counts,
  ggplot2::aes(Gate, MAGs, fill = Gate)
) +
  ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(label = MAGs), vjust = -0.35, fontface = "bold", size = 3.5) +
  ggplot2::scale_fill_manual(values = c(
    "Raw catalog" = "grey65", "Relaxed" = pal_pub[["gold"]],
    "Primary" = pal_pub[["orange"]], "Strict" = pal_pub[["blue"]]
  )) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
  ggplot2::labs(
    title = "Node filters change graph size",
    subtitle = "Filters are label-free and fixed before graph fitting",
    x = NULL, y = "Candidate MAGs"
  ) +
  theme_pub(9.5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

p_filter <- p_filter_scatter + p_filter_counts +
  patchwork::plot_layout(widths = c(1.45, 0.85)) +
  patchwork::plot_annotation(
    title = "From 780 catalog MAGs to a prespecified 93-node candidate universe",
    caption = "Abundance is closed within the recovered-MAG catalog; thresholds do not imply whole-community prevalence.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.caption = ggplot2::element_text(color = "grey35", size = 8.5)
    )
  )
save_pub(p_filter, "29-feature-contract", width = 235, height = 125)

set.seed(primary_seed + 50000L)
layout_xy <- igraph::layout_with_fr(
  graph_core,
  weights = pmax(abs(igraph::edge_attr(graph_core, "PartialCorrelation")), 1e-4),
  niter = 3000L,
  grid = "nogrid"
)
layout_xy[, 1L] <- scales::rescale(layout_xy[, 1L], to = c(-1, 1))
layout_xy[, 2L] <- scales::rescale(layout_xy[, 2L], to = c(-1, 1))
network_layout <- data.frame(
  MAG = igraph::V(graph_core)$name,
  X = layout_xy[, 1L],
  Y = layout_xy[, 2L],
  stringsAsFactors = FALSE
)
network_layout <- merge(network_layout, node_roles, by = "MAG", sort = FALSE)
network_layout <- network_layout[match(igraph::V(graph_core)$name, network_layout$MAG), ]
network_layout$Label <- ""
label_order <- node_roles$MAG[order(
  -node_roles$TopologyPriority,
  -node_roles$StableDegree70,
  -node_roles$BootstrapExpectedDegree,
  -node_roles$Degree,
  node_roles$MAG
)]
label_nodes <- head(label_order, 10L)
label_index <- match(label_nodes, network_layout$MAG)
network_layout$Label[label_index] <- ifelse(
  network_layout$Genus[label_index] == "Unclassified",
  network_layout$MAG[label_index],
  paste0(network_layout$Genus[label_index], "\n", network_layout$MAG[label_index])
)

network_edges <- data.frame(
  Source = core_edge_ends[, 1L],
  Target = core_edge_ends[, 2L],
  PartialCorrelation = igraph::edge_attr(graph_core, "PartialCorrelation"),
  BootstrapFrequency = igraph::edge_attr(graph_core, "BootstrapFrequency"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
network_edges$Sign <- ifelse(network_edges$PartialCorrelation > 0, "Positive", "Negative")
network_edges$AbsolutePartialCorrelation <- abs(network_edges$PartialCorrelation)
network_edges$X <- network_layout$X[match(network_edges$Source, network_layout$MAG)]
network_edges$Y <- network_layout$Y[match(network_edges$Source, network_layout$MAG)]
network_edges$Xend <- network_layout$X[match(network_edges$Target, network_layout$MAG)]
network_edges$Yend <- network_layout$Y[match(network_edges$Target, network_layout$MAG)]

module_levels <- sort(unique(network_layout$Module))
module_colors <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#A6761D", "#999999",
  "#332288", "#88CCEE", "#44AA99", "#DDCC77"
)[seq_along(module_levels)]
names(module_colors) <- module_levels
p_network <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = network_edges,
    ggplot2::aes(
      X, Y, xend = Xend, yend = Yend,
      color = Sign, alpha = BootstrapFrequency,
      linewidth = AbsolutePartialCorrelation
    ),
    lineend = "round"
  ) +
  ggplot2::geom_point(
    data = network_layout,
    ggplot2::aes(
      X, Y, fill = Module,
      size = Degree, shape = TopologyPriority
    ),
    color = "grey15", stroke = 0.35
  ) +
  ggrepel::geom_text_repel(
    data = network_layout[network_layout$Label != "", ],
    ggplot2::aes(X, Y, label = Label),
    seed = primary_seed + 50001L,
    size = 2.55, color = "black",
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    segment.color = "grey55",
    segment.size = 0.3
  ) +
  ggplot2::scale_color_manual(values = c("Positive" = pal_pub[["orange"]], "Negative" = pal_pub[["blue"]])) +
  ggplot2::scale_fill_manual(values = module_colors) +
  ggplot2::scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 23),
    labels = c(`FALSE` = "Other node", `TRUE` = "Topology-priority")
  ) +
  ggplot2::scale_alpha_continuous(
    range = c(0.20, 0.95), limits = c(0, 1),
    breaks = c(0.4, 0.7, 1.0),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::scale_linewidth_continuous(
    range = c(0.25, 1.6),
    breaks = c(0.05, 0.10, 0.15),
    labels = c("0.05", "0.10", "0.15")
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.2, 7.0), labels = function(x) formatC(x, format = "f", digits = 0)
  ) +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::labs(
    title = "Environment-adjusted conditional-association network",
    subtitle = paste0(
      primary_nodes_n, " non-isolates and ", primary_edges_n,
      " edges; isolates omitted only from the drawing"
    ),
    color = "Partial sign", fill = "Module",
    alpha = "Bootstrap\nfrequency", linewidth = "|Partial r|",
    size = "Degree", shape = NULL,
    caption = "Edge opacity reports 1,000 region-stratified bootstraps at the locked primary lambda."
  ) +
  ggplot2::theme_void(base_size = 9.5, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.title.position = "plot",
    plot.caption = ggplot2::element_text(color = "grey35"),
    legend.position = "right",
    plot.margin = ggplot2::margin(10, 18, 8, 8)
  )
save_pub(p_network, "29-conditional-network", width = 225, height = 175)

role_colors <- c(
  "Ultra-peripheral" = "#999999",
  "Peripheral" = "#56B4E9",
  "Connector" = "#E69F00",
  "Kinless" = "#CC79A7",
  "Provincial hub" = "#009E73",
  "Connector hub" = "#D55E00",
  "Kinless hub" = "#332288"
)
z_limit <- max(3.0, max(node_roles$Zi) + 0.35)
p_zipi <- ggplot2::ggplot(
  node_roles,
  ggplot2::aes(Pi, Zi, color = Role, size = StableDegree70)
) +
  ggplot2::geom_vline(xintercept = 0.62, linetype = 2, color = "grey45") +
  ggplot2::geom_hline(yintercept = 2.5, linetype = 2, color = "grey45") +
  ggplot2::geom_point(alpha = 0.88) +
  ggrepel::geom_text_repel(
    data = node_roles[node_roles$TopologyPriority, ],
    ggplot2::aes(label = ifelse(Genus == "Unclassified", MAG, Genus)),
    seed = primary_seed + 50002L,
    size = 2.6, color = "black", min.segment.length = 0,
    max.overlaps = Inf, segment.color = "grey55", segment.size = 0.3
  ) +
  ggplot2::scale_color_manual(values = role_colors, drop = FALSE) +
  ggplot2::scale_size_continuous(range = c(1.8, 6.2), breaks = c(0, 2, 4, 6)) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(min(-0.5, min(node_roles$Zi) - 0.2), z_limit)) +
  ggplot2::labs(
    title = "Zi-Pi is a topology map",
    subtitle = "Classic thresholds: Zi > 2.5; connector Pi > 0.62",
    x = "Participation coefficient (Pi)",
    y = "Within-module degree z-score (Zi)",
    color = "Topological role", size = "Consensus\ndegree"
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "bottom")

candidate_display <- head(topology_candidates, 10L)
candidate_display$Display <- ifelse(
  candidate_display$Genus == "Unclassified",
  candidate_display$MAG,
  paste(candidate_display$Genus, candidate_display$MAG, sep = " · ")
)
candidate_display$Display <- factor(candidate_display$Display, levels = rev(candidate_display$Display))
p_candidates <- ggplot2::ggplot(
  candidate_display,
  ggplot2::aes(BootstrapExpectedDegree, Display)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = BootstrapExpectedDegree, yend = Display),
    color = "grey70", linewidth = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = Role, size = StableDegree70), alpha = 0.9
  ) +
  ggplot2::scale_color_manual(values = role_colors, drop = FALSE) +
  ggplot2::scale_size_continuous(range = c(2.2, 6.0)) +
  ggplot2::labs(
    title = "Candidates pass a stability gate",
    subtitle = "Ranked for follow-up; none is called a keystone",
    x = "Bootstrap-expected degree", y = NULL,
    color = "Role", size = "Consensus\ndegree"
  ) +
  theme_pub(8.5) +
  ggplot2::theme(legend.position = "bottom")

p_zipi_combined <- p_zipi + p_candidates +
  patchwork::plot_layout(widths = c(1.25, 1)) +
  patchwork::plot_annotation(
    title = "Network cartography separates module hubs, connectors and peripheral nodes",
    caption = "A topology-priority label requires >=2 edges with bootstrap frequency >=0.70 plus a preregistered role or degree criterion.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.caption = ggplot2::element_text(color = "grey35", size = 8.5)
    )
  )
save_pub(p_zipi_combined, "29-zippi-roles", width = 235, height = 145)

edge_stability_plot <- primary_edge_table
edge_stability_plot <- edge_stability_plot[order(-edge_stability_plot$SelectionFrequency), ]
edge_stability_plot$EdgeRank <- seq_len(nrow(edge_stability_plot))
p_edge_stability <- ggplot2::ggplot(
  edge_stability_plot,
  ggplot2::aes(EdgeRank, SelectionFrequency, color = Sign)
) +
  ggplot2::geom_hline(yintercept = 0.70, linetype = 2, color = pal_pub[["gold"]]) +
  ggplot2::geom_hline(yintercept = 0.80, linetype = 2, color = pal_pub[["green"]]) +
  ggplot2::geom_segment(
    ggplot2::aes(xend = EdgeRank, y = 0, yend = SelectionFrequency),
    linewidth = 0.35, alpha = 0.45
  ) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::scale_color_manual(values = c("Positive" = pal_pub[["orange"]], "Negative" = pal_pub[["blue"]])) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Bootstrap edge frequency",
    subtitle = paste0(
      sum(primary_edge_table$Consensus70), "/", nrow(primary_edge_table),
      " primary edges reach 70%; ",
      sum(primary_edge_table$HighStability80), " reach 80%"
    ),
    x = "Primary edges ranked by bootstrap frequency",
    y = "Selection frequency", color = "Partial sign"
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "bottom")

sensitivity_plot <- network_sensitivity[network_sensitivity$Branch != "Primary", ]
sensitivity_plot$Branch <- factor(
  sensitivity_plot$Branch,
  levels = sensitivity_plot$Branch[order(sensitivity_plot$EdgeJaccardWithPrimary)]
)
sensitivity_plot$Label <- paste0(sensitivity_plot$Features, " nodes · ", sensitivity_plot$Edges, " edges")
p_sensitivity <- ggplot2::ggplot(
  sensitivity_plot,
  ggplot2::aes(EdgeJaccardWithPrimary, Branch)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = EdgeJaccardWithPrimary, yend = Branch),
    color = "grey72", linewidth = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = PrimaryEdgeRetention), size = 3.2
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label), hjust = -0.10, size = 2.8
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1.25), breaks = seq(0, 1, 0.25),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  ggplot2::scale_color_gradient(
    low = pal_pub[["orange"]], high = pal_pub[["blue"]],
    limits = c(0, 1), labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Specification sensitivity",
    subtitle = "Jaccard is calculated only over nodes shared with the primary branch",
    x = "Edge Jaccard with primary", y = NULL,
    color = "Primary-edge\nretention"
  ) +
  theme_pub(8.7) +
  ggplot2::theme(legend.position = "bottom")

p_stability_combined <- p_edge_stability + p_sensitivity +
  patchwork::plot_layout(widths = c(1.05, 1.25)) +
  patchwork::plot_annotation(
    title = "Edge uncertainty has two layers: resampling and analysis specification",
    caption = "Bootstrap frequencies hold the primary feature set and lambda fixed; branch comparisons deliberately change one analysis choice.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.caption = ggplot2::element_text(color = "grey35", size = 8.5)
    )
  )
save_pub(p_stability_combined, "29-edge-stability-sensitivity", width = 240, height = 145)

p_attack <- ggplot2::ggplot(
  random_curve_summary,
  ggplot2::aes(RemovedFraction)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = RandomLower95, ymax = RandomUpper95),
    fill = "grey75", alpha = 0.65
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = RandomMedian, color = "Random deletion"),
    linewidth = 0.9
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = TargetedAdaptiveDegree, color = "Adaptive degree attack"),
    linewidth = 1.0
  ) +
  ggplot2::scale_color_manual(values = c(
    "Random deletion" = "grey35",
    "Adaptive degree attack" = pal_pub[["orange"]]
  )) +
  ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1), labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Adaptive deletion fragments the graph faster",
    subtitle = paste0(
      "AUC targeted ", sprintf("%.3f", target_auc),
      " vs random median ", sprintf("%.3f", stats::median(random_auc))
    ),
    x = "Nodes removed",
    y = "Largest component / initial nodes",
    color = NULL
  ) +
  theme_pub(9) +
  ggplot2::theme(legend.position = "bottom")

null_long <- rbind(
  data.frame(Replicate = topology_null$Replicate, Metric = "Modularity", Value = topology_null$Modularity),
  data.frame(Replicate = topology_null$Replicate, Metric = "Transitivity", Value = topology_null$Transitivity)
)
null_observed <- data.frame(
  Metric = c("Modularity", "Transitivity"),
  Observed = c(observed_modularity, observed_transitivity),
  Label = paste0("Observed\np = ", formatC(topology_null_summary$EmpiricalUpperP, format = "f", digits = 3)),
  stringsAsFactors = FALSE
)
p_null <- ggplot2::ggplot(null_long, ggplot2::aes(Value)) +
  ggplot2::geom_histogram(bins = 32, fill = pal_pub[["sky"]], color = "white", linewidth = 0.25) +
  ggplot2::geom_vline(
    data = null_observed,
    ggplot2::aes(xintercept = Observed),
    color = pal_pub[["orange"]], linewidth = 0.9
  ) +
  ggplot2::geom_text(
    data = null_observed,
    ggplot2::aes(x = Observed, y = Inf, label = Label),
    color = pal_pub[["orange"]], hjust = 1.08, vjust = 1.15, size = 2.8
  ) +
  ggplot2::facet_wrap(~Metric, scales = "free_x", ncol = 1) +
  ggplot2::scale_x_continuous(labels = function(x) sprintf("%.2f", x)) +
  ggplot2::labs(
    title = "Degree-preserving topology nulls",
    subtitle = "1,000 rewires preserve every node degree",
    x = "Null statistic", y = "Rewired graphs"
  ) +
  theme_pub(8.7)

p_robustness <- p_attack + p_null +
  patchwork::plot_layout(widths = c(1.25, 0.95)) +
  patchwork::plot_annotation(
    title = "Graph robustness is a property of the inferred topology—not proof of ecosystem resilience",
    caption = "Deletion curves exclude the 47 primary isolates; both attacks begin from the same 46-node, three-component subgraph.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.caption = ggplot2::element_text(color = "grey35", size = 8.5)
    )
  )
save_pub(p_robustness, "29-robustness-null", width = 235, height = 145)

random_auc_table <- data.frame(
  Replicate = seq_len(random_replicates),
  RandomDeletionAUC = random_auc,
  stringsAsFactors = FALSE
)
node_role_counts <- as.data.frame(table(node_roles$Role), stringsAsFactors = FALSE)
names(node_role_counts) <- c("Role", "Nodes")
node_role_counts <- node_role_counts[node_role_counts$Nodes > 0, , drop = FALSE]

write_tsv(data_lineage, file.path(output_dir, "data-lineage.tsv"))
write_tsv(primary_network_summary, file.path(output_dir, "primary-network-summary.tsv"))
write_tsv(primary_edge_table, file.path(output_dir, "primary-edges.tsv"))
write_tsv_gz(edge_bootstrap, file.path(output_dir, "edge-bootstrap.tsv.gz"))
write_tsv(node_roles, file.path(output_dir, "node-roles.tsv"))
write_tsv(node_role_counts, file.path(output_dir, "node-role-counts.tsv"))
write_tsv(topology_candidates, file.path(output_dir, "topology-candidates.tsv"))
write_tsv(network_layout, file.path(output_dir, "network-layout.tsv"))
write_tsv(network_edges, file.path(output_dir, "network-edges-layout.tsv"))
write_tsv_gz(topology_null, file.path(output_dir, "degree-preserving-null.tsv.gz"))
write_tsv(topology_null_summary, file.path(output_dir, "degree-preserving-null-summary.tsv"))
write_tsv(random_curve_summary, file.path(output_dir, "deletion-robustness-curves.tsv"))
write_tsv(robustness_summary, file.path(output_dir, "deletion-robustness-summary.tsv"))
write_tsv(targeted_deletion_order, file.path(output_dir, "targeted-deletion-order.tsv"))
write_tsv_gz(random_auc_table, file.path(output_dir, "random-deletion-auc.tsv.gz"))
write_tsv(network_sensitivity, file.path(output_dir, "network-sensitivity.tsv"))
write_tsv(interpretation_boundaries, file.path(output_dir, "interpretation-boundaries.tsv"))

figure_stems <- c(
  "29-feature-contract", "29-conditional-network", "29-zippi-roles",
  "29-edge-stability-sensitivity", "29-robustness-null"
)
figure_audit <- do.call(rbind, lapply(figure_stems, function(stem) {
  do.call(rbind, lapply(c("pdf", "png", "tiff"), function(extension) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    data.frame(
      Figure = stem,
      Format = extension,
      Exists = file.exists(path),
      Bytes = if (file.exists(path)) as.numeric(file.info(path)$size) else NA_real_,
      SHA256 = if (file.exists(path)) sha256_file(path) else NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}))
write_tsv(figure_audit, file.path(output_dir, "figure-audit.tsv"))
add_check("Figures", "figure-files", nrow(figure_audit) == 15L && all(figure_audit$Exists), paste0(sum(figure_audit$Exists), "/15"))
add_check("Figures", "figure-nonempty", all(figure_audit$Bytes > 10000), paste(range(figure_audit$Bytes), collapse = "/"))
anchor_path <- file.path(figure_dir, "29-kurtz-fig2-original.png")
anchor_hash <- if (file.exists(anchor_path)) sha256_file(anchor_path) else "missing"
add_check("Figures", "anchor-present", file.exists(anchor_path) && file.info(anchor_path)$size > 100000, anchor_path)
add_check(
  "Figures", "anchor-sha256",
  identical(anchor_hash, "a67e5d64b526d35d0d4e77645c1c1494fcd83b842177fe01458fae5c78b4a4ee"),
  anchor_hash
)

result_files <- c(
  "data-lineage.tsv", "primary-network-summary.tsv", "primary-edges.tsv",
  "edge-bootstrap.tsv.gz", "node-roles.tsv", "node-role-counts.tsv",
  "topology-candidates.tsv", "network-layout.tsv", "network-edges-layout.tsv",
  "degree-preserving-null.tsv.gz", "degree-preserving-null-summary.tsv",
  "deletion-robustness-curves.tsv", "deletion-robustness-summary.tsv",
  "targeted-deletion-order.tsv", "random-deletion-auc.tsv.gz",
  "network-sensitivity.tsv", "interpretation-boundaries.tsv", "figure-audit.tsv"
)
result_paths <- file.path(output_dir, result_files)
add_check("Outputs", "analysis-result-files", all(file.exists(result_paths)), paste0(sum(file.exists(result_paths)), "/", length(result_paths)))
writeLines(
  paste(vapply(result_paths, sha256_file, character(1L)), result_files),
  file.path(output_dir, "result-checksums.sha256"),
  useBytes = TRUE
)
add_check("Outputs", "result-checksum-lines", length(readLines(file.path(output_dir, "result-checksums.sha256"))) == length(result_files), length(result_files))

chapter_text <- paste(readLines(chapter_path, warn = FALSE), collapse = "\n")
chapter_tokens <- c(
  "draft: false", "eval: true", "freeze: auto", "expected_images: 6",
  "huge 1.5.1", "StARS", "BroadRegion", "graphical lasso",
  "Zi", "Pi", "不等于生态 keystone",
  "a67e5d64b526d35d0d4e77645c1c1494fcd83b842177fe01458fae5c78b4a4ee",
  "## 这一步对应论文里的哪张图", "## 理论：",
  "## 准备工作", "## 可复制代码", "## 审计与升级",
  "## 出版级美化", "## 常见坑", "## 这段 Methods 怎么写",
  "## 换成你自己的数据怎么做", "## 参考"
)
for (token in chapter_tokens) {
  add_check(
    "Chapter", paste0("chapter-", gsub("[^a-z0-9]+", "-", tolower(token))),
    grepl(token, chapter_text, fixed = TRUE), token
  )
}
add_check("Chapter", "no-source-theme-dependency", !grepl("source(\"R/theme_pub.R\")", chapter_text, fixed = TRUE), "inline plotting functions")
add_check("Chapter", "no-meta-copy-self-contained", !grepl("本篇可独立跑通", chapter_text, fixed = TRUE), "reader-facing copy")
add_check("Chapter", "no-interaction-overclaim", !grepl("证明了菌群互作", chapter_text, fixed = TRUE), "association language")
add_check("Chapter", "seed-inline", grepl("set.seed(20260729)", chapter_text, fixed = TRUE), "deterministic code")

write_tsv(checks, file.path(output_dir, "validation-checks.tsv"))
failures <- checks[checks$Status != "PASS", , drop = FALSE]
summary <- list(
  status = if (nrow(failures) == 0L) "passed" else "failed",
  article = 29L,
  seed = primary_seed,
  source_metagenomes = 500L,
  inference_units = nrow(primary_x),
  raw_mag_features = ncol(abundance),
  primary_mag_features = ncol(primary_x),
  residual_design_rank = attr(primary$transformed, "design_rank"),
  residual_df = attr(primary$transformed, "residual_df"),
  selected_lambda = primary$lambda,
  selected_path_index = primary$opt_index,
  primary_edges = primary_edges_n,
  positive_edges = sum(primary_edge_table$Sign == "Positive"),
  negative_edges = sum(primary_edge_table$Sign == "Negative"),
  nonisolates = primary_nodes_n,
  isolates = length(primary_features) - primary_nodes_n,
  nonisolate_components = igraph::components(graph_core)$no,
  largest_component = max(component_sizes),
  modules = length(unique(node_roles$Module)),
  consensus_edges_70 = sum(primary_edge_table$Consensus70),
  high_stability_edges_80 = sum(primary_edge_table$HighStability80),
  bootstrap_replicates = bootstrap_replicates,
  bootstrap_attempts = bootstrap_attempts,
  topology_priority_candidates = nrow(topology_candidates),
  observed_modularity = observed_modularity,
  modularity_null_upper_p = topology_null_summary$EmpiricalUpperP[topology_null_summary$Metric == "Modularity"],
  observed_transitivity = observed_transitivity,
  transitivity_null_upper_p = topology_null_summary$EmpiricalUpperP[topology_null_summary$Metric == "Transitivity"],
  targeted_deletion_auc = target_auc,
  random_deletion_median_auc = stats::median(random_auc),
  sensitivity_branches = nrow(network_sensitivity),
  checksum_entries = checksum_entries,
  checks = nrow(checks),
  passed = sum(checks$Status == "PASS"),
  failed = nrow(failures),
  generated_figure_files = nrow(figure_audit),
  package_versions = list(
    huge = as.character(utils::packageVersion("huge")),
    igraph = as.character(utils::packageVersion("igraph")),
    ggplot2 = as.character(utils::packageVersion("ggplot2")),
    R = paste(R.version$major, R.version$minor, sep = ".")
  ),
  generated_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
jsonlite::write_json(
  summary, file.path(output_dir, "validation-summary.json"),
  auto_unbox = TRUE, pretty = TRUE
)

log_msg(
  "checks=", nrow(checks), "; passed=", sum(checks$Status == "PASS"),
  "; failed=", nrow(failures)
)
if (nrow(failures) > 0L) {
  print(failures, row.names = FALSE)
  stop("Article 29 validation failed.", call. = FALSE)
}
log_msg("Article 29 validation passed.")
cat("Article 29 validation passed:", nrow(checks), "checks.\n")
