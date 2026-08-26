#!/usr/bin/env bash
set -euo pipefail

PIPELINE_RELEASE="5.5.0"
NEXTFLOW_VERSION="26.04.0"

export NXF_APPTAINER_CACHEDIR="/shared/containers/nf-core-mag-${PIPELINE_RELEASE}"
export NXF_VER="${NEXTFLOW_VERSION}"

nextflow run nf-core/mag \
  -r "${PIPELINE_RELEASE}" \
  -profile apptainer \
  -params-file params.publication.yml \
  -c hpc.slurm.config \
  -resume
