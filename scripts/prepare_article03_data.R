args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: prepare_article03_data.R <meta_all.tsv> ",
    "<design_audit.tsv> <summary.json>"
  )
}

input_path <- args[[1]]
audit_path <- args[[2]]
summary_path <- args[[3]]

metadata <- utils::read.delim(
  input_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)

required_columns <- c(
  "Sample_ID", "Study", "Group", "Age", "Gender", "BMI", "Library_Size"
)
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(metadata$Sample_ID)) {
  stop("Sample_ID values must be unique")
}

discovery_studies <- c("FR-CRC", "AT-CRC", "CN-CRC", "US-CRC", "DE-CRC")
validation_studies <- c("IT-CRC", "IT-CRC-2", "JP-CRC")
study_order <- c(discovery_studies, validation_studies)
group_order <- c("CTR", "CRC")

crc <- metadata[
  metadata$Study %in% study_order & metadata$Group %in% group_order,
  required_columns,
  drop = FALSE
]
crc$Role <- ifelse(
  crc$Study %in% discovery_studies,
  "Discovery / meta-analysis",
  "Independent validation"
)

safe_median <- function(x, scale = 1, digits = 1) {
  observed <- x[!is.na(x)]
  if (length(observed) == 0L) {
    return(NA_real_)
  }
  round(stats::median(observed) / scale, digits)
}

summarise_arm <- function(study, group) {
  arm <- crc[crc$Study == study & crc$Group == group, , drop = FALSE]
  n <- nrow(arm)
  gender_observed <- sum(!is.na(arm$Gender))
  female_pct <- if (gender_observed == 0L) {
    NA_real_
  } else {
    round(100 * sum(arm$Gender == "F", na.rm = TRUE) / gender_observed, 1)
  }

  data.frame(
    Study = study,
    Role = unique(arm$Role),
    Group = ifelse(group == "CTR", "Control", "CRC"),
    N = n,
    Age_observed = sum(!is.na(arm$Age)),
    Age_complete_pct = round(100 * mean(!is.na(arm$Age)), 1),
    Age_median_years = safe_median(arm$Age),
    Gender_observed = gender_observed,
    Gender_complete_pct = round(100 * mean(!is.na(arm$Gender)), 1),
    Female_pct_observed = female_pct,
    BMI_observed = sum(!is.na(arm$BMI)),
    BMI_complete_pct = round(100 * mean(!is.na(arm$BMI)), 1),
    BMI_median = safe_median(arm$BMI),
    Library_observed = sum(!is.na(arm$Library_Size)),
    Library_complete_pct = round(100 * mean(!is.na(arm$Library_Size)), 1),
    Library_median_m_reads = safe_median(
      arm$Library_Size,
      scale = 1e6,
      digits = 2
    ),
    stringsAsFactors = FALSE
  )
}

audit <- do.call(
  rbind,
  lapply(
    study_order,
    function(study) {
      do.call(
        rbind,
        lapply(group_order, function(group) summarise_arm(study, group))
      )
    }
  )
)
rownames(audit) <- NULL

arm_minima <- vapply(
  split(crc, crc$Study),
  function(x) min(table(x$Group)),
  numeric(1)
)

summary_values <- list(
  analysis_samples = nrow(crc),
  analysis_studies = length(unique(crc$Study)),
  crc_cases = sum(crc$Group == "CRC"),
  controls = sum(crc$Group == "CTR"),
  smallest_group_n = min(audit$N),
  largest_group_n = max(audit$N),
  median_smallest_arm = stats::median(arm_minima),
  age_missing = sum(is.na(crc$Age)),
  gender_missing = sum(is.na(crc$Gender)),
  bmi_missing = sum(is.na(crc$BMI)),
  library_size_missing = sum(is.na(crc$Library_Size)),
  cohorts_with_complete_age = sum(
    vapply(split(crc$Age, crc$Study), function(x) !anyNA(x), logical(1))
  ),
  cohorts_with_complete_gender = sum(
    vapply(split(crc$Gender, crc$Study), function(x) !anyNA(x), logical(1))
  ),
  cohorts_with_complete_bmi = sum(
    vapply(split(crc$BMI, crc$Study), function(x) !anyNA(x), logical(1))
  ),
  median_library_size = stats::median(crc$Library_Size),
  minimum_library_size = min(crc$Library_Size),
  maximum_library_size = max(crc$Library_Size)
)

stopifnot(
  summary_values$analysis_samples == 768L,
  summary_values$analysis_studies == 8L,
  summary_values$crc_cases == 386L,
  summary_values$controls == 382L,
  summary_values$smallest_group_n == 24L,
  summary_values$largest_group_n == 74L,
  summary_values$median_smallest_arm == 49,
  summary_values$age_missing == 60L,
  summary_values$gender_missing == 60L,
  summary_values$bmi_missing == 71L,
  summary_values$library_size_missing == 0L,
  summary_values$cohorts_with_complete_age == 7L,
  summary_values$cohorts_with_complete_gender == 7L,
  summary_values$cohorts_with_complete_bmi == 3L,
  nrow(audit) == 16L
)

dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  audit,
  audit_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

json_number <- function(value) {
  if (length(value) != 1L || is.na(value)) {
    stop("Summary JSON values must be scalar and non-missing")
  }
  format(value, scientific = FALSE, trim = TRUE)
}
json_lines <- vapply(
  names(summary_values),
  function(name) {
    sprintf('  "%s": %s', name, json_number(summary_values[[name]]))
  },
  character(1)
)
json <- paste0(
  "{\n",
  paste(json_lines, collapse = ",\n"),
  "\n}\n"
)
writeLines(json, summary_path, useBytes = TRUE)
cat(json)
