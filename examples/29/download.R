options(timeout = max(600, getOption("timeout")))
base_url <- "https://raw.githubusercontent.com/petemeng/metagenomics-best-practices/e219a2cd01609cc208a89e291cbd363ed7f2fbd8/"
files <- read.delim(paste0(base_url, "examples/29/files.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(!grepl("(^/|(^|/)\\.\\.(/|$))", files$path)))
for (i in seq_len(nrow(files))) {
  target <- files$path[i]
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(target)) {
    temporary <- tempfile(tmpdir = dirname(target))
    download.file(paste0(base_url, files$source[i]), temporary, mode = "wb", quiet = TRUE)
    if (unname(tools::md5sum(temporary)) != files$md5[i]) {
      unlink(temporary)
      stop("Downloaded file checksum mismatch: ", target)
    }
    if (!file.rename(temporary, target)) stop("Cannot save: ", target)
  }
  if (unname(tools::md5sum(target)) != files$md5[i]) {
    stop("Existing input differs from this version: ", target,
         ". Use a new folder; do not overwrite your own data.")
  }
}
message("Inputs downloaded; saved results have NOT been recomputed.")
