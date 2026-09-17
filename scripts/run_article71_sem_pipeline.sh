#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${project_root}"

python_bin="${PYTHON_BIN:-python3}"
rscript_bin="${RSCRIPT_BIN:-Rscript}"
run_root="${ARTICLE71_RUN_ROOT:-results/article71-rebuild}"
cache_dir="${ARTICLE71_CACHE_DIR:-db/sem/article71}"

export R_LIBS_USER="${project_root}/.r-lib:${HOME}/R/library${R_LIBS_USER:+:${R_LIBS_USER}}"

"${python_bin}" - <<'PY'
import importlib
for package in ("numpy", "pandas"):
    importlib.import_module(package)
PY

"${rscript_bin}" -e '
required <- c("sandwich", "lmtest", "car", "jsonlite", "mediation", "ggplot2", "patchwork")
available <- vapply(
  required,
  function(package) suppressWarnings(requireNamespace(package, quietly = TRUE)),
  logical(1)
)
missing <- required[!available]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "))
}
'

mkdir -p "${run_root}/prepared" "${run_root}/analysis" "${run_root}/bundle" "${run_root}/figures"

"${python_bin}" scripts/download_article71_sem_data.py \
  --cache-dir "${cache_dir}"

"${python_bin}" scripts/prepare_article71_sem.py \
  --cache-dir "${cache_dir}" \
  --output-dir "${run_root}/prepared"

"${rscript_bin}" scripts/run_article71_sem_models.R \
  --input-dir "${run_root}/prepared" \
  --output-dir "${run_root}/analysis"

cp "${run_root}/prepared/"* "${run_root}/bundle/"
cp "${run_root}/analysis/"* "${run_root}/bundle/"

"${rscript_bin}" scripts/plot_article71_sem.R \
  --input-dir "${run_root}/bundle" \
  --figure-dir "${run_root}/figures"

"${python_bin}" scripts/freeze_article71_sem.py \
  --project-root "${project_root}" \
  --work-dir "${run_root}/bundle" \
  --output-dir "${run_root}/frozen"

"${python_bin}" scripts/validate_article71_sem.py \
  --project-root "${project_root}" \
  --frozen-dir "${run_root}/frozen" \
  --qa-dir "${run_root}/qa"

printf 'article71_rebuild\t%s\n' "${project_root}/${run_root}"
