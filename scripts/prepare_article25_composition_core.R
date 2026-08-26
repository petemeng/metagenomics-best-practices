#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, digits = 17)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
primary_seed <- 20260725L
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
required <- c("source-rda", "output-dir", "notice")
missing <- setdiff(required, names(args))
if (length(missing) > 0L) {
  stop("Missing arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
}

packages <- c("curatedMetagenomicData", "digest")
unavailable <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(unavailable) > 0L) {
  stop("Missing R packages: ", paste(unavailable, collapse = ", "), call. = FALSE)
}

source_rda <- normalizePath(args[["source-rda"]], mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
notice_path <- normalizePath(args[["notice"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(notice_path), recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}
write_profile <- function(feature_map, matrix, path) {
  con <- gzfile(path, open = "wt", compression = 9L)
  on.exit(close(con), add = TRUE)
  out <- data.frame(feature_map, matrix, check.names = FALSE)
  utils::write.table(
    out, con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}

expected_source_sha256 <- "f42713b0d59d15f65172ead2375c2ea2c1b6309b44cb2b29231edb190e56ab23"
observed_source_sha256 <- sha256_file(source_rda)
if (!identical(observed_source_sha256, expected_source_sha256)) {
  stop("HMP_2012 relative-abundance RDA checksum mismatch.", call. = FALSE)
}
if (!identical(as.character(utils::packageVersion("curatedMetagenomicData")), "3.12.0")) {
  stop("curatedMetagenomicData 3.12.0 is required.", call. = FALSE)
}

source_environment <- new.env(parent = emptyenv())
loaded_names <- load(source_rda, envir = source_environment)
if (!identical(loaded_names, "2021-03-31.HMP_2012.relative_abundance")) {
  stop("Unexpected object name in HMP_2012 RDA.", call. = FALSE)
}
source_profile <- source_environment[[loaded_names]]
stopifnot(
  is.matrix(source_profile),
  identical(dim(source_profile), c(1188L, 748L)),
  !anyDuplicated(rownames(source_profile)),
  !anyDuplicated(colnames(source_profile)),
  all(is.finite(source_profile)),
  min(source_profile) >= 0,
  max(source_profile) <= 100 + 1e-8
)

metadata_all <- curatedMetagenomicData::sampleMetadata
metadata <- metadata_all[metadata_all$study_name == "HMP_2012", , drop = FALSE]
metadata <- metadata[match(colnames(source_profile), metadata$sample_id), , drop = FALSE]
stopifnot(
  nrow(metadata) == 748L,
  !anyNA(metadata$sample_id),
  identical(metadata$sample_id, colnames(source_profile)),
  !anyDuplicated(metadata$sample_id),
  all(metadata$study_condition == "control"),
  all(metadata$disease == "healthy")
)

body_levels_raw <- c("oralcavity", "nasalcavity", "skin", "stool", "vagina")
broad_body_labels <- c(
  oralcavity = "Oral cavity",
  nasalcavity = "Nasal cavity",
  skin = "Skin",
  stool = "Stool",
  vagina = "Vagina"
)
stopifnot(setequal(unique(metadata$body_site), body_levels_raw))

# Match the seven habitat columns displayed in HMP Figure 3. Left and right
# retroauricular crease are combined before selecting one sample per subject.
habitat_labels <- c(
  anterior_nares = "Anterior nares",
  retroauricular_crease = "Retroauricular crease",
  buccal_mucosa = "Buccal mucosa",
  tongue_dorsum = "Tongue dorsum",
  supragingival_plaque = "Supragingival plaque",
  stool = "Stool",
  posterior_fornix = "Posterior fornix"
)
analytic_habitat_raw <- metadata$body_subsite
analytic_habitat_raw[analytic_habitat_raw %in% c(
  "l_retroauricular_crease", "r_retroauricular_crease"
)] <- "retroauricular_crease"
included_habitat <- analytic_habitat_raw %in% names(habitat_labels)
analytic_habitat <- unname(habitat_labels[analytic_habitat_raw])
analytic_habitat[!included_habitat] <- NA_character_

read_depth <- suppressWarnings(as.numeric(metadata$number_reads))
read_depth_sort <- ifelse(is.finite(read_depth), read_depth, -1)
group_key <- ifelse(
  included_habitat,
  paste(metadata$subject_id, analytic_habitat_raw, sep = "::"),
  NA_character_
)
selection_order <- which(included_habitat)
selection_order <- selection_order[order(
  metadata$subject_id[selection_order],
  match(analytic_habitat_raw[selection_order], names(habitat_labels)),
  -read_depth_sort[selection_order],
  metadata$sample_id[selection_order]
)]
candidate_rank <- rep(NA_integer_, nrow(metadata))
candidate_rank[selection_order] <- ave(
  seq_along(selection_order),
  group_key[selection_order],
  FUN = seq_along
)
representative <- !is.na(candidate_rank) & candidate_rank == 1L
stopifnot(
  sum(representative) == length(unique(group_key[included_habitat])),
  all(table(group_key[representative]) == 1L),
  length(unique(metadata$subject_id)) == 103L
)

metadata_frozen <- data.frame(
  SampleID = metadata$sample_id,
  SubjectID = metadata$subject_id,
  BroadBodySite = unname(broad_body_labels[metadata$body_site]),
  BodySiteRaw = metadata$body_site,
  BodySubsite = metadata$body_subsite,
  Habitat = analytic_habitat,
  PrimaryHabitatIncluded = included_habitat,
  Age = metadata$age,
  Gender = metadata$gender,
  Country = metadata$country,
  NumberReads = read_depth,
  NumberBases = suppressWarnings(as.numeric(metadata$number_bases)),
  SequencingPlatform = metadata$sequencing_platform,
  NCBIAccession = metadata$NCBI_accession,
  Representative = representative,
  RepresentativeRule = "highest number_reads per subject x selected HMP habitat; sample_id tie-break",
  SourceColumn = seq_len(nrow(metadata)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

selection_audit <- data.frame(
  SampleID = metadata$sample_id,
  SubjectID = metadata$subject_id,
  BroadBodySite = unname(broad_body_labels[metadata$body_site]),
  BodySubsite = metadata$body_subsite,
  Habitat = analytic_habitat,
  PrimaryHabitatIncluded = included_habitat,
  NumberReads = read_depth,
  CandidateRank = candidate_rank,
  Representative = representative,
  SelectionReason = ifelse(
    representative,
    "maximum finite number_reads; lexical sample_id tie-break",
    ifelse(
      included_habitat,
      "lower number_reads or lexical tie-break within subject x selected habitat",
      "body subsite is outside the seven HMP Figure 3 habitat columns"
    )
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

lineages <- rownames(source_profile)
terminal_token <- sub("^.*\\|", "", lineages)
terminal_prefix <- substr(terminal_token, 1L, 3L)

clean_token <- function(x) {
  x <- sub("^[a-z]__", "", x)
  gsub("_", " ", x, fixed = TRUE)
}
extract_token <- function(lineage, prefix) {
  tokens <- strsplit(lineage, "\\|", fixed = FALSE)[[1L]]
  hit <- tokens[startsWith(tokens, prefix)]
  if (length(hit) == 0L) return("Unclassified")
  clean_token(hit[[1L]])
}
make_display_labels <- function(rank_lineages, rank_tokens) {
  labels <- clean_token(rank_tokens)
  generic <- tolower(labels) %in% c("unclassified", "unknown", "")
  if (any(generic)) {
    parents <- vapply(
      rank_lineages[generic],
      function(x) {
        tokens <- strsplit(x, "\\|", fixed = FALSE)[[1L]]
        if (length(tokens) <= 1L) return("Unclassified")
        clean_token(tokens[[length(tokens) - 1L]])
      },
      character(1L)
    )
    labels[generic] <- paste(parents, "unclassified")
  }
  if (anyDuplicated(labels)) {
    labels <- make.unique(labels, sep = " [lineage ")
    labels <- ifelse(grepl(" \\[lineage ", labels), paste0(labels, "]"), labels)
  }
  labels
}

extract_rank <- function(prefix, rank_name, id_prefix) {
  keep <- terminal_prefix == prefix
  rank_lineages <- lineages[keep]
  rank_tokens <- terminal_token[keep]
  values_percent <- source_profile[keep, , drop = FALSE]
  raw_sum_percent <- colSums(values_percent)
  if (any(!is.finite(raw_sum_percent)) || any(raw_sum_percent <= 0)) {
    stop(rank_name, " contains an invalid sample denominator.", call. = FALSE)
  }
  values_fraction <- sweep(values_percent / 100, 2L, raw_sum_percent / 100, "/")
  labels <- make_display_labels(rank_lineages, rank_tokens)
  feature_map <- data.frame(
    FeatureID = sprintf("%s%04d", id_prefix, seq_along(rank_lineages)),
    Rank = rank_name,
    Label = labels,
    Lineage = rank_lineages,
    Phylum = vapply(rank_lineages, extract_token, character(1L), prefix = "p__"),
    CoreEligible = rank_name == "Species" & !grepl("unclassified", rank_tokens, ignore.case = TRUE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(values_fraction) <- feature_map$FeatureID
  list(
    feature_map = feature_map,
    matrix = values_fraction,
    raw_sum_percent = raw_sum_percent
  )
}

rank_profiles <- list(
  Phylum = extract_rank("p__", "Phylum", "PH"),
  Genus = extract_rank("g__", "Genus", "GE"),
  Species = extract_rank("s__", "Species", "SP")
)
stopifnot(
  all(vapply(rank_profiles, function(x) nrow(x$matrix) > 0L, logical(1L))),
  all(vapply(rank_profiles, function(x) max(abs(colSums(x$matrix) - 1)) < 1e-12, logical(1L)))
)

rank_closure_audit <- do.call(
  rbind,
  lapply(names(rank_profiles), function(rank_name) {
    x <- rank_profiles[[rank_name]]
    data.frame(
      Rank = rank_name,
      Features = nrow(x$matrix),
      Samples = ncol(x$matrix),
      RawSumPercentMin = min(x$raw_sum_percent),
      RawSumPercentMedian = stats::median(x$raw_sum_percent),
      RawSumPercentMax = max(x$raw_sum_percent),
      ReclosedSumMaxError = max(abs(colSums(x$matrix) - 1)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)

profile_paths <- c(
  Phylum = file.path(output_dir, "phylum-relative-abundance.tsv.gz"),
  Genus = file.path(output_dir, "genus-relative-abundance.tsv.gz"),
  Species = file.path(output_dir, "species-relative-abundance.tsv.gz")
)
for (rank_name in names(rank_profiles)) {
  write_profile(
    rank_profiles[[rank_name]]$feature_map,
    rank_profiles[[rank_name]]$matrix,
    profile_paths[[rank_name]]
  )
}

metadata_path <- file.path(output_dir, "sample-metadata.tsv")
selection_path <- file.path(output_dir, "sample-selection-audit.tsv")
closure_path <- file.path(output_dir, "rank-closure-audit.tsv")
contract_path <- file.path(output_dir, "analysis-contract.tsv")
resource_path <- file.path(output_dir, "resource-manifest.tsv")

write_tsv(metadata_frozen, metadata_path)
write_tsv(selection_audit, selection_path)
write_tsv(rank_closure_audit, closure_path)

representative_counts <- table(factor(metadata_frozen$Habitat[representative], levels = unname(habitat_labels)))
analysis_contract <- data.frame(
  Item = c(
    "seed", "source_release", "source_profile_unit", "source_samples",
    "source_subjects", "primary_inference_unit", "representative_rule",
    "habitat_order", "composition_closure", "composition_summary",
    "primary_rank", "primary_detection_threshold", "primary_prevalence_threshold",
    "prevalence_interval", "bootstrap_replicates", "bootstrap_unit",
    "sensitivity_detection_thresholds", "sensitivity_prevalence_thresholds",
    "unclassified_core_rule", "causal_scope"
  ),
  Value = c(
    as.character(primary_seed), "curatedMetagenomicData resource 2021-03-31",
    "within-rank relative abundance fraction", as.character(nrow(metadata_frozen)),
    as.character(length(unique(metadata_frozen$SubjectID))),
    paste0("one subject x selected HMP habitat; n=", sum(representative)),
    "highest number_reads; lexical sample_id tie-break",
    paste(unname(habitat_labels), collapse = " | "),
    "reclose separately within phylum, genus, and species ranks",
    "mean of sample-wise closed compositions; never pool reads",
    "species", "0.0001 fraction (0.01%)", "0.80",
    "Wilson score 95% interval", "1000", "subject within habitat",
    "0 | 0.00001 | 0.0001 | 0.001",
    "0.50 | 0.80 | 0.90 | 1.00",
    "exclude terminal species labels containing unclassified",
    "descriptive habitat-specific carriage and composition; no causal interpretation"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(analysis_contract, contract_path)

resource_manifest <- data.frame(
  Resource = c(
    "HMP_2012 relative abundance", "HMP_2012 sample metadata",
    "HMP Nature Figure 3", "Taxonomic profile processing"
  ),
  Version = c("2021-03-31", "curatedMetagenomicData 3.12.0", "Nature 486 (2012)", "bioBakery 3 reprocessing"),
  Source = c(
    "ExperimentHub curatedMetagenomicData/HMP_2012",
    "curatedMetagenomicData::sampleMetadata",
    "https://doi.org/10.1038/nature11234",
    "MetaPhlAn3 relative_abundance"
  ),
  SHA256 = c(observed_source_sha256, "frozen-subset-checksum-in-file-manifest", "5afd411dee3ad8604263da4a3069f429942b6605674ce7cc8a3d0a161dca66c2", "not-applicable"),
  License = c(
    "Artistic-2.0 package plus original-study terms",
    "Artistic-2.0 package plus original-study terms",
    "CC BY-NC-SA 3.0",
    "input provenance; see curatedMetagenomicData resource metadata"
  ),
  Notes = c(
    "1188 lineage rows x 748 samples; values supplied as percent",
    paste0("748 samples; 103 subjects; representatives across seven HMP Figure 3 habitats: ", paste(names(representative_counts), representative_counts, sep = "=", collapse = ", ")),
    "Original figure is an attribution anchor; generated figures use reprocessed profiles",
    "Rank tables are reclosed separately and expressed as fractions"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(resource_manifest, resource_path)

notice <- c(
  "Article 25 data notice — HMP_2012 composition and operational core microbiome",
  "",
  "The source is the 2021-03-31 HMP_2012 relative_abundance resource distributed by curatedMetagenomicData 3.12.0.",
  paste0("The source RDA SHA256 is ", observed_source_sha256, "."),
  "The profile was generated by uniform bioBakery 3 reprocessing and is not the original 2012 paper's exact taxonomic table.",
  "All 748 source samples and 103 subjects are retained in the frozen ledger.",
  "The primary inferential unit is one subject x selected HMP Figure 3 habitat, represented by the sample with the largest finite number_reads and a lexical sample_id tie-break.",
  "The seven habitats are anterior nares, retroauricular crease, buccal mucosa, tongue dorsum, supragingival plaque, stool, and posterior fornix.",
  "Left and right retroauricular crease are combined before one sample per subject is selected.",
  "This representative rule prevents repeated samples from one subject and habitat from being counted as independent prevalence observations.",
  "Other sparsely sampled body subsites remain in the frozen source ledger but are excluded from the primary habitat analysis.",
  "Phylum, genus, and species tables are extracted from explicit terminal MetaPhlAn ranks and reclosed separately to sum to one per sample.",
  "A zero means not detected by this versioned profiling workflow; it does not prove biological absence.",
  "The operational primary core uses species abundance >=0.01% and prevalence >=80% within each selected habitat.",
  "Wilson intervals, subject bootstrap inclusion frequency, and a detection-by-prevalence threshold grid quantify core-definition uncertainty.",
  "Terminal species labels containing 'unclassified' are excluded from core membership but retained in composition denominators.",
  "Body-site composition and carriage are descriptive; they do not identify causal host-microbe effects.",
  "Original HMP Figure 3 is CC BY-NC-SA 3.0 and is included with article-level attribution.",
  ""
)
writeLines(notice, notice_path, useBytes = TRUE)

payloads <- sort(basename(list.files(output_dir, full.names = TRUE)))
payloads <- setdiff(payloads, "file-checksums.sha256")
checksums <- vapply(file.path(output_dir, payloads), sha256_file, character(1L))
writeLines(
  paste(checksums, payloads),
  file.path(output_dir, "file-checksums.sha256"),
  useBytes = TRUE
)

cat("Article 25 frozen data created.\n")
cat("Source samples:", nrow(metadata_frozen), "\n")
cat("Independent subjects:", length(unique(metadata_frozen$SubjectID)), "\n")
cat("Representative subject-habitat units:", sum(representative), "\n")
cat("Representative counts:", paste(names(representative_counts), representative_counts, sep = "=", collapse = ", "), "\n")
cat("Features:", paste(names(rank_profiles), vapply(rank_profiles, function(x) nrow(x$matrix), integer(1L)), sep = "=", collapse = ", "), "\n")
