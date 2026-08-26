#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(TZ = "UTC")
set.seed(20260722)

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
required <- c("project-root", "input-dir", "output-dir", "figure-dir")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c(
  "Matrix", "vegan", "ape", "ggplot2", "patchwork", "scales",
  "jsonlite", "digest", "magick"
)
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

project_root <- normalizePath(args[["project-root"]], mustWork = TRUE)
input_dir <- normalizePath(args[["input-dir"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
figure_dir <- normalizePath(args[["figure-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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
    add_check("Frozen input", paste0("sha256-", relative), identical(observed, expected), observed)
    expected_files <- c(expected_files, relative)
  }
  payloads <- basename(list.files(directory, full.names = TRUE))
  payloads <- sort(setdiff(payloads, "file-checksums.sha256"))
  add_check(
    "Frozen input", "checksum-manifest-complete",
    identical(payloads, sort(expected_files)),
    paste0("payloads=", length(payloads), "; entries=", length(expected_files))
  )
  length(lines)
}

checksum_entries <- verify_checksum_manifest(input_dir)
notice_path <- file.path(project_root, "data", "small", "22-data-NOTICE.txt")
notice <- readLines(notice_path, warn = FALSE)
add_check("Frozen input", "notice-present", length(notice) > 10L, notice_path)
add_check(
  "Frozen input", "notice-separates-lineages",
  any(grepl("independent datasets", notice, fixed = TRUE)),
  "Human and MAG branches are explicitly independent"
)

resource_manifest <- utils::read.delim(
  file.path(input_dir, "resource-manifest.tsv"), check.names = FALSE, quote = ""
)
filter_contract <- utils::read.delim(
  file.path(input_dir, "filter-contract.tsv"), check.names = FALSE, quote = ""
)
feature_audit <- utils::read.delim(
  file.path(input_dir, "human-feature-audit.tsv"), check.names = FALSE, quote = ""
)
expected_ids <- c("EH7091", "EH7086", "61153444", "61153429", "61153471")
add_check(
  "Lineage", "resource-identifiers",
  identical(as.character(resource_manifest$RepositoryID), expected_ids),
  paste(resource_manifest$RepositoryID, collapse = ",")
)
add_check(
  "Lineage", "resource-manifest-rows",
  nrow(resource_manifest) == 5L,
  nrow(resource_manifest)
)
add_check(
  "Lineage", "filter-contract-spaces",
  identical(filter_contract$FeatureSpace, c("Species", "Gene family", "MAG catalog")),
  paste(filter_contract$FeatureSpace, collapse = ",")
)
add_check(
  "Lineage", "gene-raw-publisher-md5",
  identical(resource_manifest$RawMD5[resource_manifest$RepositoryID == "EH7086"],
            "49caa4e88cbd51e7f5700cbcf4590e55"),
  resource_manifest$RawMD5[resource_manifest$RepositoryID == "EH7086"]
)
add_check(
  "Lineage", "figshare-publisher-md5",
  identical(
    resource_manifest$RawMD5[resource_manifest$RepositoryID %in% c("61153444", "61153429", "61153471")],
    c(
      "0c1762f31a5c4473dec78571a7b74287",
      "402095b04e2c5518cbec462f41528d4f",
      "51288407e9927f23a25b210edcda7b47"
    )
  ),
  "Three Figshare MD5 values"
)

read_wide_profile <- function(path, annotation_columns) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  tab <- utils::read.delim(
    con, check.names = FALSE, quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
  stopifnot(all(annotation_columns %in% names(tab)))
  sample_columns <- setdiff(names(tab), annotation_columns)
  abundance <- t(data.matrix(tab[, sample_columns, drop = FALSE]))
  rownames(abundance) <- sample_columns
  colnames(abundance) <- as.character(tab[[annotation_columns[[1L]]]])
  annotations <- tab[, annotation_columns, drop = FALSE]
  rm(tab)
  gc(verbose = FALSE)
  list(abundance = abundance, annotations = annotations)
}

species_input <- read_wide_profile(
  file.path(input_dir, "species-relative-abundance.tsv.gz"),
  c("Feature", "Species")
)
gene_input <- read_wide_profile(
  file.path(input_dir, "gene-family-prevalence10.tsv.gz"),
  "GeneFamily"
)
mag_input <- read_wide_profile(
  file.path(input_dir, "mag-relative-abundance.tsv.gz"),
  c("MAG", "Taxonomy")
)
human_metadata <- utils::read.delim(
  file.path(input_dir, "human-sample-metadata.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
mag_metadata <- utils::read.delim(
  file.path(input_dir, "hot-spring-sample-metadata.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
mag_recruitment <- utils::read.delim(
  file.path(input_dir, "mag-recruitment.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)

species_percent <- species_input$abundance
gene10_native <- gene_input$abundance
mag_native <- mag_input$abundance

add_check("Input shape", "species-shape", identical(dim(species_percent), c(24L, 298L)), paste(dim(species_percent), collapse = "x"))
add_check("Input shape", "gene-prevalence10-shape", identical(dim(gene10_native), c(24L, 415581L)), paste(dim(gene10_native), collapse = "x"))
add_check("Input shape", "mag-shape", identical(dim(mag_native), c(500L, 780L)), paste(dim(mag_native), collapse = "x"))
add_check("Input shape", "human-profiles", nrow(human_metadata) == 24L, nrow(human_metadata))
add_check("Input shape", "human-subjects", length(unique(human_metadata$subject_id)) == 15L, length(unique(human_metadata$subject_id)))
add_check("Input shape", "mag-samples", nrow(mag_metadata) == 500L, nrow(mag_metadata))
add_check("Input shape", "hot-springs", length(unique(mag_metadata$hotspring)) == 56L, length(unique(mag_metadata$hotspring)))
add_check("Alignment", "human-sample-order", identical(rownames(species_percent), human_metadata$sample_id) && identical(rownames(gene10_native), human_metadata$sample_id), "species/gene/metadata")
add_check("Alignment", "mag-sample-order", identical(rownames(mag_native), mag_metadata$sample) && identical(rownames(mag_native), mag_recruitment$sample_id), "BIOM/metadata/recruitment")
add_check("Value audit", "finite-nonnegative-species", all(is.finite(species_percent)) && min(species_percent) >= 0, range(species_percent))
add_check("Value audit", "finite-nonnegative-gene", all(is.finite(gene10_native)) && min(gene10_native) >= 0, range(gene10_native))
add_check("Value audit", "finite-nonnegative-mag", all(is.finite(mag_native)) && min(mag_native) >= 0, range(mag_native))
add_check("Value audit", "gene-unstratified", !any(grepl("|", colnames(gene10_native), fixed = TRUE)), "No taxon-stratified IDs")
add_check("Value audit", "gene-no-special", !any(colnames(gene10_native) %in% c("UNMAPPED", "UNGROUPED", "UNINTEGRATED")), "No special rows")

species_native_sums <- rowSums(species_percent)
mag_native_sums <- rowSums(mag_native)
add_check(
  "Denominator", "species-native-percent-sums",
  min(species_native_sums) > 94 && max(species_native_sums) < 101,
  sprintf("%.6f..%.6f", min(species_native_sums), max(species_native_sums))
)
add_check(
  "Denominator", "mag-native-catalog-closure",
  max(abs(mag_native_sums - 1)) < 2e-5,
  sprintf("max deviation %.8g", max(abs(mag_native_sums - 1)))
)

close_rows <- function(x) {
  totals <- rowSums(x)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Every sample must have a positive finite denominator.", call. = FALSE)
  }
  x / totals
}

species <- close_rows(species_percent)
gene_prevalence <- colMeans(gene10_native > 0)
gene_matrices <- list(
  `10%` = close_rows(gene10_native[, gene_prevalence >= 0.10, drop = FALSE]),
  `20%` = close_rows(gene10_native[, gene_prevalence >= 0.20, drop = FALSE]),
  `50%` = close_rows(gene10_native[, gene_prevalence >= 0.50, drop = FALSE])
)
gene <- gene_matrices[["20%"]]
mag <- close_rows(mag_native)

add_check("Filter", "gene-features-10", ncol(gene_matrices[["10%"]]) == 415581L, ncol(gene_matrices[["10%"]]))
add_check("Filter", "gene-features-20", ncol(gene) == 178928L, ncol(gene))
add_check("Filter", "gene-features-50", ncol(gene_matrices[["50%"]]) == 1896L, ncol(gene_matrices[["50%"]]))
add_check("Denominator", "species-closed", max(abs(rowSums(species) - 1)) < 1e-12, max(abs(rowSums(species) - 1)))
add_check("Denominator", "gene-closed", max(abs(rowSums(gene) - 1)) < 1e-12, max(abs(rowSums(gene) - 1)))
add_check("Denominator", "mag-reclosed", max(abs(rowSums(mag) - 1)) < 1e-12, max(abs(rowSums(mag) - 1)))

hill_values <- function(x) {
  richness <- rowSums(x > 0)
  entropy_terms <- ifelse(x > 0, x * log(x), 0)
  q1 <- exp(-rowSums(entropy_terms))
  q2 <- 1 / rowSums(x^2)
  data.frame(q0 = richness, q1 = q1, q2 = q2, check.names = FALSE)
}

hill_long <- function(x, feature_space, dataset) {
  h <- hill_values(x)
  n <- nrow(h)
  data.frame(
    Dataset = dataset,
    FeatureSpace = feature_space,
    SampleID = rep(rownames(h), times = 3L),
    Q = rep(c(0L, 1L, 2L), each = n),
    Metric = rep(c("Observed richness", "exp(Shannon)", "Inverse Simpson"), each = n),
    Value = c(h$q0, h$q1, h$q2),
    stringsAsFactors = FALSE
  )
}

alpha <- rbind(
  hill_long(species, "Species", "AsnicarF_2017"),
  hill_long(gene, "Gene families", "AsnicarF_2017"),
  hill_long(mag, "MAG catalog", "Western US hot springs")
)
alpha$SubjectID <- human_metadata$subject_id[
  match(alpha$SampleID, human_metadata$sample_id)
]
alpha$AgeCategory <- human_metadata$age_category[
  match(alpha$SampleID, human_metadata$sample_id)
]
alpha$HotSpring <- mag_metadata$hotspring[
  match(alpha$SampleID, mag_metadata$sample)
]
alpha$Temperature <- mag_metadata$temperature[
  match(alpha$SampleID, mag_metadata$sample)
]
alpha$DetectionBoundary <- ifelse(
  alpha$Q == 0L,
  "Detection/filter dependent",
  "Abundance-weighted effective number"
)

for (space in unique(alpha$FeatureSpace)) {
  wide <- reshape(
    alpha[alpha$FeatureSpace == space, c("SampleID", "Q", "Value")],
    idvar = "SampleID", timevar = "Q", direction = "wide"
  )
  add_check(
    "Alpha diversity", paste0("hill-monotonic-", gsub(" ", "-", tolower(space))),
    all(wide$Value.0 + 1e-10 >= wide$Value.1) && all(wide$Value.1 + 1e-10 >= wide$Value.2),
    paste0("samples=", nrow(wide))
  )
}
add_check("Alpha diversity", "alpha-row-count", nrow(alpha) == 1644L, nrow(alpha))

aitchison_distance <- function(x, pseudocount) {
  logx <- log(x + pseudocount)
  clr <- logx - rowMeans(logx)
  stats::dist(clr, method = "euclidean")
}

distance_row <- function(feature_space, distance_name, variant, features, d) {
  values <- as.vector(d)
  data.frame(
    FeatureSpace = feature_space,
    Distance = distance_name,
    Variant = variant,
    Features = features,
    Samples = attr(d, "Size"),
    PairCount = length(values),
    Minimum = min(values),
    Median = stats::median(values),
    Maximum = max(values),
    stringsAsFactors = FALSE
  )
}

distances <- list(
  species_bray = vegan::vegdist(species, method = "bray"),
  species_jaccard = vegan::vegdist(species > 0, method = "jaccard", binary = TRUE),
  species_aitchison_1e6 = aitchison_distance(species, 1e-6),
  species_aitchison_1e5 = aitchison_distance(species, 1e-5),
  gene_bray = vegan::vegdist(gene, method = "bray"),
  gene_jaccard = vegan::vegdist(gene > 0, method = "jaccard", binary = TRUE),
  gene_aitchison_1e8 = aitchison_distance(gene, 1e-8),
  gene_aitchison_1e7 = aitchison_distance(gene, 1e-7),
  mag_bray = vegan::vegdist(mag, method = "bray"),
  mag_jaccard = vegan::vegdist(mag > 0, method = "jaccard", binary = TRUE)
)

beta_audit <- rbind(
  distance_row("Species", "Bray-Curtis", "Primary closure", ncol(species), distances$species_bray),
  distance_row("Species", "Binary Jaccard", "Presence/absence", ncol(species), distances$species_jaccard),
  distance_row("Species", "Aitchison", "Pseudocount 1e-06", ncol(species), distances$species_aitchison_1e6),
  distance_row("Species", "Aitchison", "Pseudocount 1e-05", ncol(species), distances$species_aitchison_1e5),
  distance_row("Gene families", "Bray-Curtis", "Prevalence >=20%", ncol(gene), distances$gene_bray),
  distance_row("Gene families", "Binary Jaccard", "Prevalence >=20%", ncol(gene), distances$gene_jaccard),
  distance_row("Gene families", "Aitchison", "Pseudocount 1e-08", ncol(gene), distances$gene_aitchison_1e8),
  distance_row("Gene families", "Aitchison", "Pseudocount 1e-07", ncol(gene), distances$gene_aitchison_1e7),
  distance_row("MAG catalog", "Bray-Curtis", "Catalog closure", ncol(mag), distances$mag_bray),
  distance_row("MAG catalog", "Binary Jaccard", "Presence/absence", ncol(mag), distances$mag_jaccard)
)
add_check("Beta diversity", "beta-audit-rows", nrow(beta_audit) == 10L, nrow(beta_audit))
add_check("Beta diversity", "human-pair-counts", all(beta_audit$PairCount[beta_audit$Samples == 24L] == 276L), unique(beta_audit$PairCount[beta_audit$Samples == 24L]))
add_check("Beta diversity", "mag-pair-counts", all(beta_audit$PairCount[beta_audit$Samples == 500L] == 124750L), unique(beta_audit$PairCount[beta_audit$Samples == 500L]))
add_check("Beta diversity", "distances-finite-nonnegative", all(is.finite(beta_audit$Minimum)) && min(beta_audit$Minimum) >= 0, min(beta_audit$Minimum))

pcoa_scores <- function(d, feature_space, dataset) {
  fit <- ape::pcoa(d, correction = "cailliez")
  coordinates <- if (!is.null(fit$vectors.cor)) fit$vectors.cor else fit$vectors
  relative <- if ("Rel_corr_eig" %in% names(fit$values)) {
    fit$values$Rel_corr_eig
  } else {
    fit$values$Relative_eig
  }
  raw_eigen <- fit$values$Eigenvalues
  data <- data.frame(
    Dataset = dataset,
    FeatureSpace = feature_space,
    SampleID = rownames(coordinates),
    Axis1 = coordinates[, 1L],
    Axis2 = coordinates[, 2L],
    Axis1Variance = relative[[1L]],
    Axis2Variance = relative[[2L]],
    NegativeEigenvalues = sum(raw_eigen < -1e-10),
    NegativeEigenvalueMass = sum(abs(raw_eigen[raw_eigen < -1e-10])),
    Correction = "Cailliez",
    CorrectionNote = fit$note,
    stringsAsFactors = FALSE
  )
  data
}

ordination <- rbind(
  pcoa_scores(distances$species_bray, "Species", "AsnicarF_2017"),
  pcoa_scores(distances$gene_bray, "Gene families", "AsnicarF_2017"),
  pcoa_scores(distances$mag_bray, "MAG catalog", "Western US hot springs")
)
ordination$SubjectID <- human_metadata$subject_id[match(ordination$SampleID, human_metadata$sample_id)]
ordination$FamilyRole <- human_metadata$family_role[match(ordination$SampleID, human_metadata$sample_id)]
ordination$AgeCategory <- human_metadata$age_category[match(ordination$SampleID, human_metadata$sample_id)]
ordination$DaysFromFirstCollection <- human_metadata$days_from_first_collection[match(ordination$SampleID, human_metadata$sample_id)]
ordination$HotSpring <- mag_metadata$hotspring[match(ordination$SampleID, mag_metadata$sample)]
ordination$Temperature <- mag_metadata$temperature[match(ordination$SampleID, mag_metadata$sample)]
add_check("Ordination", "ordination-row-count", nrow(ordination) == 548L, nrow(ordination))
add_check("Ordination", "ordination-finite", all(is.finite(ordination$Axis1)) && all(is.finite(ordination$Axis2)), "Axis1/Axis2")
add_check("Ordination", "ordination-correction-recorded", all(ordination$Correction == "Cailliez"), unique(ordination$Correction))

gene_threshold_alpha <- do.call(
  rbind,
  lapply(names(gene_matrices), function(threshold) {
    out <- hill_long(gene_matrices[[threshold]], "Gene families", "AsnicarF_2017")
    out$Threshold <- threshold
    out$Features <- ncol(gene_matrices[[threshold]])
    out
  })
)
gene_bray10 <- vegan::vegdist(gene_matrices[["10%"]], method = "bray")
gene_bray50 <- vegan::vegdist(gene_matrices[["50%"]], method = "bray")
rho_gene_bray10 <- stats::cor(as.vector(gene_bray10), as.vector(distances$gene_bray), method = "spearman")
rho_gene_bray50 <- stats::cor(as.vector(gene_bray50), as.vector(distances$gene_bray), method = "spearman")
rho_gene_aitchison <- stats::cor(as.vector(distances$gene_aitchison_1e8), as.vector(distances$gene_aitchison_1e7), method = "spearman")
rho_species_aitchison <- stats::cor(as.vector(distances$species_aitchison_1e6), as.vector(distances$species_aitchison_1e5), method = "spearman")
mag_richness <- alpha[alpha$FeatureSpace == "MAG catalog" & alpha$Q == 0L, c("SampleID", "Value")]
mag_richness$Recruitment <- mag_recruitment$total_hit_rate[match(mag_richness$SampleID, mag_recruitment$sample_id)]
rho_mag_recruitment <- stats::cor(mag_richness$Value, mag_richness$Recruitment, method = "spearman")

sensitivity_rows <- list()
for (threshold in names(gene_matrices)) {
  alpha_part <- gene_threshold_alpha[gene_threshold_alpha$Threshold == threshold, ]
  for (q in 0:2) {
    values <- alpha_part$Value[alpha_part$Q == q]
    sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
      SensitivityType = "Prevalence filter",
      FeatureSpace = "Gene families",
      Variant = threshold,
      Features = ncol(gene_matrices[[threshold]]),
      Metric = paste0("Median Hill q=", q),
      Value = stats::median(values),
      Reference = "20% primary",
      Interpretation = "Feature-universe sensitivity; not a biological group test",
      stringsAsFactors = FALSE
    )
  }
}
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
  SensitivityType = "Distance geometry", FeatureSpace = "Gene families",
  Variant = "Bray prevalence 10% vs 20%", Features = ncol(gene_matrices[["10%"]]),
  Metric = "Pairwise-distance Spearman rho", Value = rho_gene_bray10,
  Reference = "20% primary", Interpretation = "Pairs are non-independent; no p-value",
  stringsAsFactors = FALSE
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
  SensitivityType = "Distance geometry", FeatureSpace = "Gene families",
  Variant = "Bray prevalence 50% vs 20%", Features = ncol(gene_matrices[["50%"]]),
  Metric = "Pairwise-distance Spearman rho", Value = rho_gene_bray50,
  Reference = "20% primary", Interpretation = "Pairs are non-independent; no p-value",
  stringsAsFactors = FALSE
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
  SensitivityType = "Zero replacement", FeatureSpace = "Gene families",
  Variant = "Aitchison 1e-08 vs 1e-07", Features = ncol(gene),
  Metric = "Pairwise-distance Spearman rho", Value = rho_gene_aitchison,
  Reference = "1e-08 primary", Interpretation = "Pairs are non-independent; no p-value",
  stringsAsFactors = FALSE
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
  SensitivityType = "Zero replacement", FeatureSpace = "Species",
  Variant = "Aitchison 1e-06 vs 1e-05", Features = ncol(species),
  Metric = "Pairwise-distance Spearman rho", Value = rho_species_aitchison,
  Reference = "1e-06 primary", Interpretation = "Pairs are non-independent; no p-value",
  stringsAsFactors = FALSE
)
sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
  SensitivityType = "Catalog recovery", FeatureSpace = "MAG catalog",
  Variant = "Richness vs read recruitment", Features = ncol(mag),
  Metric = "Descriptive Spearman rho", Value = rho_mag_recruitment,
  Reference = "500 samples", Interpretation = "No p-value and no causal interpretation",
  stringsAsFactors = FALSE
)
sensitivity <- do.call(rbind, sensitivity_rows)
add_check("Sensitivity", "sensitivity-row-count", nrow(sensitivity) == 14L, nrow(sensitivity))
add_check("Sensitivity", "distance-rhos-finite", all(is.finite(c(rho_gene_bray10, rho_gene_bray50, rho_gene_aitchison, rho_species_aitchison))), c(rho_gene_bray10, rho_gene_bray50, rho_gene_aitchison, rho_species_aitchison))
add_check("Sensitivity", "recruitment-rho-finite", is.finite(rho_mag_recruitment), rho_mag_recruitment)
add_check(
  "Catalog recovery", "mean-read-recruitment",
  abs(mean(mag_recruitment$total_hit_rate) - 0.12983) < 1e-4,
  sprintf("%.8f", mean(mag_recruitment$total_hit_rate))
)
add_check(
  "Catalog recovery", "read-recruitment-range",
  min(mag_recruitment$total_hit_rate) > 0.008 && max(mag_recruitment$total_hit_rate) < 0.779,
  sprintf("%.6f..%.6f", min(mag_recruitment$total_hit_rate), max(mag_recruitment$total_hit_rate))
)

human_repeats <- do.call(
  rbind,
  lapply(split(human_metadata, human_metadata$subject_id), function(part) {
    data.frame(
      UnitType = "Human subject",
      UnitID = part$subject_id[[1L]],
      Profiles = nrow(part),
      Repeated = if (nrow(part) > 1L) "Yes" else "No",
      SampleIDs = paste(part$sample_id, collapse = ";"),
      Dataset = "AsnicarF_2017",
      Boundary = "Repeated profiles are not independent biological replicates",
      stringsAsFactors = FALSE
    )
  })
)
spring_repeats <- do.call(
  rbind,
  lapply(split(mag_metadata, mag_metadata$hotspring), function(part) {
    data.frame(
      UnitType = "Hot spring",
      UnitID = part$hotspring[[1L]],
      Profiles = nrow(part),
      Repeated = if (nrow(part) > 1L) "Yes" else "No",
      SampleIDs = paste(part$sample, collapse = ";"),
      Dataset = "Western US hot springs",
      Boundary = "Site clustering requires blocked inference in later analyses",
      stringsAsFactors = FALSE
    )
  })
)
repeated_audit <- rbind(human_repeats, spring_repeats)
add_check("Design", "repeated-audit-units", nrow(repeated_audit) == 71L, nrow(repeated_audit))
add_check("Design", "human-profile-ledger", sum(human_repeats$Profiles) == 24L, sum(human_repeats$Profiles))
add_check("Design", "mag-profile-ledger", sum(spring_repeats$Profiles) == 500L, sum(spring_repeats$Profiles))

data_lineage <- data.frame(
  Dataset = c("AsnicarF_2017", "AsnicarF_2017", "Western US hot springs"),
  FeatureSpace = c("Species", "Gene families", "MAG catalog"),
  RepositoryID = c("EH7091", "EH7086", "Figshare 61153444"),
  Samples = c(24L, 24L, 500L),
  IndependentUnits = c(15L, 15L, 56L),
  FrozenFeatures = c(298L, 415581L, 780L),
  PrimaryFeatures = c(298L, 178928L, 780L),
  NativeScale = c("Percent", "Relative fraction", "Catalog-relative fraction"),
  PrimaryDenominator = c("Observed species", "Filtered ordinary gene families", "Recovered MAG catalog"),
  Relationship = c("Sample-matched human branch", "Sample-matched human branch", "Independent environmental branch"),
  stringsAsFactors = FALSE
)

interpretation_boundaries <- data.frame(
  Topic = c(
    "Hill q=0", "Hill q=1", "Hill q=2", "Species versus genes",
    "MAG closure", "Bray-Curtis", "Binary Jaccard", "Aitchison",
    "PCoA", "Inference"
  ),
  AuthorizedInterpretation = c(
    "Observed richness conditional on detection and filtering",
    "Effective number weighted by common features",
    "Effective number dominated by abundant features",
    "Same human samples but different feature universes and measurement models",
    "Composition within 780 recovered MAGs only; recruitment remains separate",
    "Difference in closed abundance profiles",
    "Difference in detected feature membership",
    "Log-ratio geometry after declared filtering and zero replacement",
    "Descriptive coordinates with raw negative eigenvalues and Cailliez correction recorded",
    "No group test, permutation, or pairwise-distance p-value in Article 22"
  ),
  ProhibitedClaim = c(
    "Complete biological richness", "Raw feature count", "Raw feature count",
    "Which omics layer is more diverse", "Percent of all community reads",
    "Taxonomic turnover without denominator context", "Abundance difference",
    "Pseudocount-free truth", "Significant group separation", "Biological significance"
  ),
  stringsAsFactors = FALSE
)

# Publication graphics -------------------------------------------------------
pal_space <- c("Species" = "#0072B2", "Gene families" = "#D55E00", "MAG catalog" = "#009E73")
pal_q <- c("q = 0" = "#0072B2", "q = 1" = "#E69F00", "q = 2" = "#009E73")
pal_age <- c("newborn" = "#56B4E9", "child" = "#E69F00", "adult" = "#CC79A7")
theme_pub <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.25),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey45"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold")
    )
}
save_pub <- function(plot, stem, width_mm, height_mm) {
  base <- file.path(figure_dir, stem)
  ggplot2::ggsave(
    paste0(base, ".pdf"), plot, width = width_mm, height = height_mm,
    units = "mm", device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".png"), plot, width = width_mm, height = height_mm,
    units = "mm", dpi = 350, bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".tiff"), plot, width = width_mm, height = height_mm,
    units = "mm", dpi = 350, compression = "lzw", bg = "white"
  )
}

cards <- data.frame(
  FeatureSpace = c("Species", "Gene families", "MAG catalog"),
  X = 1:3,
  Y = 1,
  Label = c(
    "Human species\n24 profiles | 15 subjects\n298 MetaPhlAn features",
    "Human gene families\n24 profiles | 15 subjects\n178,928 UniRef90 families",
    "Hot-spring MAG catalog\n500 samples | 56 springs\n780 recovered MAGs"
  )
)
p_cards <- ggplot2::ggplot(cards, ggplot2::aes(X, Y, fill = FeatureSpace)) +
  ggplot2::geom_tile(width = 0.90, height = 0.72, color = "white", linewidth = 1, alpha = 0.84) +
  ggplot2::geom_text(ggplot2::aes(label = Label), size = 3.15, lineheight = 1.05) +
  ggplot2::scale_fill_manual(values = pal_space, guide = "none") +
  ggplot2::scale_x_continuous(limits = c(0.45, 3.55), expand = c(0, 0)) +
  ggplot2::scale_y_continuous(limits = c(0.58, 1.42), expand = c(0, 0)) +
  ggplot2::labs(
    title = "A. Input universes",
    subtitle = "Species and genes are sample matched; the MAG catalog is an independent environmental dataset",
    x = NULL,
    y = NULL
  ) +
  theme_pub(9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )

flow <- data.frame(
  FeatureSpace = rep(c("Species", "Gene families", "MAG catalog"), each = 3L),
  Stage = rep(1:3, times = 3L),
  Label = c(
    "MetaPhlAn\npercentage", "298 species\n(no prevalence filter)", "Observed-species\nclosure",
    "HUMAnN\nrelative fraction", ">=5/24 profiles\n(20% prevalence)", "Filtered-gene\nclosure",
    "Genome coverage\nrelative fraction", "780 recovered MAGs\n(no cohort filter)", "Recovered-MAG\ncatalog closure"
  )
)
flow$FeatureSpace <- factor(flow$FeatureSpace, levels = rev(c("Species", "Gene families", "MAG catalog")))
arrows <- expand.grid(
  FeatureSpace = levels(flow$FeatureSpace),
  Stage = c(1.46, 2.46),
  stringsAsFactors = FALSE
)
arrows$FeatureSpace <- factor(arrows$FeatureSpace, levels = levels(flow$FeatureSpace))
p_flow <- ggplot2::ggplot(flow, ggplot2::aes(Stage, FeatureSpace, fill = FeatureSpace)) +
  ggplot2::geom_tile(width = 0.88, height = 0.68, color = "white", linewidth = 0.8, alpha = 0.82) +
  ggplot2::geom_segment(
    data = arrows,
    ggplot2::aes(x = Stage, xend = Stage + 0.08, y = FeatureSpace, yend = FeatureSpace),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(2.2, "mm")),
    linewidth = 0.45,
    color = "grey30"
  ) +
  ggplot2::geom_text(ggplot2::aes(label = Label), size = 3, lineheight = 0.95) +
  ggplot2::scale_fill_manual(values = pal_space, guide = "none") +
  ggplot2::scale_x_continuous(
    breaks = 1:3,
    labels = c("Native scale", "Feature rule", "Diversity denominator"),
    limits = c(0.5, 3.5)
  ) +
  ggplot2::labs(
    title = "B. Analysis denominator",
    subtitle = "Read recruitment to the MAG catalog remains a separate sample-level quantity",
    x = NULL,
    y = NULL
  ) +
  theme_pub(9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
p_resolution <- p_cards / p_flow + patchwork::plot_layout(heights = c(0.43, 0.57)) +
  patchwork::plot_annotation(
    title = "Diversity begins with a feature universe",
    subtitle = "Two matched human profiles plus one independent recovered-genome catalog"
  )

alpha_plot <- alpha
alpha_plot$QLabel <- factor(
  paste0("q = ", alpha_plot$Q),
  levels = c("q = 0", "q = 1", "q = 2")
)
alpha_plot$FeatureSpace <- factor(alpha_plot$FeatureSpace, levels = c("Species", "Gene families", "MAG catalog"))
p_alpha <- ggplot2::ggplot(alpha_plot, ggplot2::aes(QLabel, Value, fill = QLabel, color = QLabel)) +
  ggplot2::geom_violin(scale = "width", trim = TRUE, alpha = 0.28, linewidth = 0.35) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75, linewidth = 0.35) +
  ggplot2::geom_jitter(width = 0.11, height = 0, size = 0.55, alpha = 0.35) +
  ggplot2::facet_wrap(~FeatureSpace, scales = "free_y", nrow = 1L) +
  ggplot2::scale_fill_manual(values = pal_q, guide = "none") +
  ggplot2::scale_color_manual(values = pal_q, guide = "none") +
  ggplot2::scale_y_log10(labels = scales::label_number(big.mark = ",", accuracy = 1)) +
  ggplot2::labs(
    title = "Hill numbers separate richness from dominance",
    subtitle = "Each panel has its own y range; absolute values are conditional on the declared feature universe",
    x = NULL,
    y = "Effective number of features (log scale)",
    caption = "q=0: observed richness | q=1: exp(Shannon) | q=2: inverse Simpson"
  ) +
  theme_pub(10)

human_ord <- ordination[ordination$Dataset == "AsnicarF_2017", ]
human_ord$AgeCategory <- factor(human_ord$AgeCategory, levels = c("newborn", "child", "adult"))
make_human_ordination <- function(space, show_legend = TRUE) {
  dat <- human_ord[human_ord$FeatureSpace == space, ]
  dat <- dat[order(dat$SubjectID, dat$DaysFromFirstCollection), ]
  x_pct <- 100 * unique(dat$Axis1Variance)[[1L]]
  y_pct <- 100 * unique(dat$Axis2Variance)[[1L]]
  ggplot2::ggplot(dat, ggplot2::aes(Axis1, Axis2)) +
    ggplot2::geom_path(ggplot2::aes(group = SubjectID), color = "grey72", linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(color = AgeCategory), size = 2.3, alpha = 0.9) +
    ggplot2::scale_color_manual(values = pal_age, drop = FALSE, guide = if (show_legend) "legend" else "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(space, " (n=24)"),
      subtitle = "Same human profiles",
      x = sprintf("PCoA1 (%.1f%%)", x_pct),
      y = sprintf("PCoA2 (%.1f%%)", y_pct),
      color = "Age category"
    ) +
    theme_pub(9)
}
mag_ord <- ordination[ordination$FeatureSpace == "MAG catalog", ]
mag_x_pct <- 100 * unique(mag_ord$Axis1Variance)[[1L]]
mag_y_pct <- 100 * unique(mag_ord$Axis2Variance)[[1L]]
p_mag_ord <- ggplot2::ggplot(mag_ord, ggplot2::aes(Axis1, Axis2, color = Temperature)) +
  ggplot2::geom_point(size = 1.45, alpha = 0.75) +
  ggplot2::scale_color_viridis_c(option = "C", na.value = "grey75") +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "MAG catalog (n=500)",
    subtitle = "Independent hot-spring samples",
    x = sprintf("PCoA1 (%.1f%%)", mag_x_pct),
    y = sprintf("PCoA2 (%.1f%%)", mag_y_pct),
    color = "Temperature (C)"
  ) +
  theme_pub(9)
p_beta <- make_human_ordination("Species", TRUE) +
  make_human_ordination("Gene families", FALSE) +
  p_mag_ord +
  patchwork::plot_layout(widths = c(1, 1, 1.08)) +
  patchwork::plot_annotation(
    title = "Bray-Curtis PCoA after explicit closure",
    subtitle = "Cailliez-corrected coordinates; raw negative eigenvalues remain in the audit table"
  )

gene_q0 <- gene_threshold_alpha[gene_threshold_alpha$Q == 0L, ]
gene_q0$Threshold <- factor(gene_q0$Threshold, levels = c("10%", "20%", "50%"))
p_filter <- ggplot2::ggplot(gene_q0, ggplot2::aes(Threshold, Value, fill = Threshold)) +
  ggplot2::geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.75) +
  ggplot2::geom_jitter(width = 0.10, size = 1, alpha = 0.55) +
  ggplot2::scale_fill_manual(values = c("10%" = "#56B4E9", "20%" = "#D55E00", "50%" = "#009E73"), guide = "none") +
  ggplot2::scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  ggplot2::labs(
    title = "A. Prevalence filter",
    x = "Minimum sample prevalence",
    y = "Observed gene families (log scale)"
  ) +
  theme_pub(9)

stability <- data.frame(
  Comparison = factor(
    c("Gene Bray: 10% vs 20%", "Gene Bray: 50% vs 20%", "Gene Aitchison: 1e-8 vs 1e-7", "Species Aitchison: 1e-6 vs 1e-5"),
    levels = rev(c("Gene Bray: 10% vs 20%", "Gene Bray: 50% vs 20%", "Gene Aitchison: 1e-8 vs 1e-7", "Species Aitchison: 1e-6 vs 1e-5"))
  ),
  Rho = c(rho_gene_bray10, rho_gene_bray50, rho_gene_aitchison, rho_species_aitchison),
  FeatureSpace = c("Gene families", "Gene families", "Gene families", "Species")
)
p_stability <- ggplot2::ggplot(stability, ggplot2::aes(Rho, Comparison, color = FeatureSpace)) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = Rho, yend = Comparison), color = "grey75", linewidth = 1) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", Rho)), hjust = 1.2, size = 2.8, color = "black") +
  ggplot2::scale_color_manual(values = pal_space, guide = "none") +
  ggplot2::scale_x_continuous(limits = c(0, 1.03), breaks = seq(0, 1, 0.25)) +
  ggplot2::labs(
    title = "B. Distance stability",
    x = "Pairwise-distance Spearman rho",
    y = NULL
  ) +
  theme_pub(9)

p_recruit <- ggplot2::ggplot(mag_richness, ggplot2::aes(100 * Recruitment, Value)) +
  ggplot2::geom_point(size = 1.2, alpha = 0.45, color = "#009E73") +
  ggplot2::geom_hline(
    yintercept = stats::median(mag_richness$Value),
    color = "grey35", linetype = "dashed", linewidth = 0.55
  ) +
  ggplot2::annotate(
    "label", x = Inf, y = Inf,
    label = sprintf("Descriptive rho = %.3f\nNo p-value", rho_mag_recruitment),
    hjust = 1.05, vjust = 1.15, size = 2.8, label.size = 0.25
  ) +
  ggplot2::labs(
    title = "C. Read recruitment",
    x = "Reads recruited to 780 MAGs (%)",
    y = "Observed MAG richness"
  ) +
  theme_pub(9)
p_sensitivity <- p_filter + p_stability + p_recruit +
  patchwork::plot_layout(widths = c(0.82, 1.18, 1.00)) +
  patchwork::plot_annotation(
    title = "Sensitivity checks are part of the diversity result",
    subtitle = "Distance-pair rho is descriptive (no p-value); dashed line marks median MAG richness"
  )

figure_specs <- data.frame(
  Stem = c(
    "22-resolution-boundaries", "22-alpha-hill-numbers",
    "22-beta-ordination", "22-sensitivity-recruitment"
  ),
  WidthMM = c(190, 180, 210, 220),
  HeightMM = c(112, 105, 100, 110),
  stringsAsFactors = FALSE
)
plots <- list(p_resolution, p_alpha, p_beta, p_sensitivity)
for (i in seq_len(nrow(figure_specs))) {
  save_pub(plots[[i]], figure_specs$Stem[[i]], figure_specs$WidthMM[[i]], figure_specs$HeightMM[[i]])
}

visible_labels <- c(
  "Diversity begins with a feature universe",
  "Hill numbers separate richness from dominance",
  "Bray-Curtis PCoA after explicit closure",
  "Filtering, zero replacement, and catalog recovery are analysis inputs",
  cards$Label, flow$Label, stability$Comparison
)
add_check("Graphics", "visible-labels-ascii", all(iconv(visible_labels, to = "ASCII", sub = NA) == visible_labels), "English ASCII labels")
for (i in seq_len(nrow(figure_specs))) {
  stem <- figure_specs$Stem[[i]]
  width_px <- round(figure_specs$WidthMM[[i]] / 25.4 * 350)
  height_px <- round(figure_specs$HeightMM[[i]] / 25.4 * 350)
  for (extension in c("pdf", "png", "tiff")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    add_check("Graphics", paste0(stem, "-", extension, "-exists"), file.exists(path) && file.info(path)$size > 10000, if (file.exists(path)) file.info(path)$size else "missing")
  }
  for (extension in c("png", "tiff")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    info <- magick::image_info(magick::image_read(path))
    add_check(
      "Graphics", paste0(stem, "-", extension, "-350dpi-pixels"),
      abs(info$width[[1L]] - width_px) <= 2L && abs(info$height[[1L]] - height_px) <= 2L,
      paste0(info$width[[1L]], "x", info$height[[1L]], "; expected ", width_px, "x", height_px)
    )
  }
}

# Machine-readable outputs --------------------------------------------------
write_tsv(data_lineage, file.path(output_dir, "data-lineage.tsv"))
write_tsv(alpha, file.path(output_dir, "alpha-diversity.tsv"))
write_tsv(beta_audit, file.path(output_dir, "beta-distance-audit.tsv"))
write_tsv(ordination, file.path(output_dir, "ordination-scores.tsv"))
write_tsv(sensitivity, file.path(output_dir, "sensitivity-audit.tsv"))
write_tsv(repeated_audit, file.path(output_dir, "repeated-measures-audit.tsv"))
write_tsv(interpretation_boundaries, file.path(output_dir, "interpretation-boundaries.tsv"))
write_tsv(checks, file.path(output_dir, "validation-audit.tsv"))

checks_failed <- sum(checks$Status == "FAIL")
checks_passed <- sum(checks$Status == "PASS")
summary <- list(
  status = if (checks_failed == 0L) "passed" else "failed",
  human_profiles = nrow(human_metadata),
  human_subjects = length(unique(human_metadata$subject_id)),
  species_features = ncol(species),
  gene_features_prevalence10 = ncol(gene_matrices[["10%"]]),
  gene_features_primary = ncol(gene),
  gene_features_prevalence50 = ncol(gene_matrices[["50%"]]),
  mag_samples = nrow(mag_metadata),
  hot_springs = length(unique(mag_metadata$hotspring)),
  mag_features = ncol(mag),
  mean_mag_read_recruitment = mean(mag_recruitment$total_hit_rate),
  mag_richness_recruitment_rho = rho_mag_recruitment,
  biological_group_tests = 0L,
  permutations_run = 0L,
  checksum_entries = checksum_entries,
  alpha_rows = nrow(alpha),
  beta_audit_rows = nrow(beta_audit),
  ordination_rows = nrow(ordination),
  sensitivity_rows = nrow(sensitivity),
  checks_passed = checks_passed,
  checks_failed = checks_failed,
  random_seed = 20260722L,
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  package_versions = as.list(vapply(packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1L)))
)
jsonlite::write_json(
  summary,
  file.path(output_dir, "validation-summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 10
)

log_lines <- c(
  "Article 22 alpha/beta diversity validation",
  paste0("Status: ", summary$status),
  paste0("Random seed: ", summary$random_seed),
  paste0("Frozen checksum entries: ", checksum_entries),
  paste0("Human profiles/subjects: ", summary$human_profiles, "/", summary$human_subjects),
  paste0("Species features: ", summary$species_features),
  paste0("Gene features (10%/20%/50%): ", summary$gene_features_prevalence10, "/", summary$gene_features_primary, "/", summary$gene_features_prevalence50),
  paste0("MAG samples/hot springs/features: ", summary$mag_samples, "/", summary$hot_springs, "/", summary$mag_features),
  paste0("Mean MAG read recruitment: ", sprintf("%.8f", summary$mean_mag_read_recruitment)),
  paste0("Descriptive MAG richness-recruitment rho: ", sprintf("%.6f", summary$mag_richness_recruitment_rho)),
  "Biological group tests: 0",
  "Permutations run: 0",
  paste0("Checks passed/failed: ", checks_passed, "/", checks_failed)
)
writeLines(log_lines, file.path(output_dir, "validation.log"), useBytes = TRUE)

cat(paste(log_lines, collapse = "\n"), "\n", sep = "")
if (checks_failed > 0L) {
  failed <- checks[checks$Status == "FAIL", , drop = FALSE]
  print(failed, row.names = FALSE)
  stop("Article 22 validation failed.", call. = FALSE)
}
