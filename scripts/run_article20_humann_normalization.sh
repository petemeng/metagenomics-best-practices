#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
biobakery_prefix="${BIOBAKERY_ENV_PREFIX:-${HOME}/miniconda3/envs/metagenome-biobakery-2026.07}"
article19_dir="${ARTICLE19_DIR:-${root}/data/small/19-humann3-frozen}"
output_dir="${ARTICLE20_FROZEN_DIR:-${root}/data/small/20-functional-normalization-frozen}"

export PYTHONHASHSEED=0
"${biobakery_prefix}/bin/python" \
  "${root}/scripts/freeze_article20_humann_normalization.py" \
  --environment-prefix "${biobakery_prefix}" \
  --article19-dir "${article19_dir}" \
  --output-dir "${output_dir}"
