#!/usr/bin/env python3
"""Run visible chapter-24 diagnostic code using only its listed public assets.

This validates the explanatory tables, not a substitute for the separate
MaAsLin3/ANCOM-BC2/ALDEx2 full-model validation.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import yaml


def executable_chunks(source):
    chunks = []
    for block in re.findall(r'^```\{r[^}]*\}\n(.*?)^```\s*$', source, re.M | re.S):
        options = yaml.safe_load('\n'.join(re.findall(r'^#\| ?(.*)$', block, re.M))) or {}
        if options.get('eval') is False:
            continue
        if options.get('include') is False or options.get('echo') is False:
            raise ValueError('A required diagnostic definition is hidden')
        chunks.append(re.sub(r'^#\|.*\n?', '', block, flags=re.M))
    return chunks


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, required=True)
    p.add_argument('--work-root', type=Path, required=True)
    p.add_argument('--rscript', default='Rscript')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    root = a.project_root.resolve()
    work = a.work_root.resolve()
    work.mkdir(parents=True, exist_ok=False)
    rows = list(csv.DictReader((root / 'examples/24/files.tsv').open(), delimiter='\t'))
    for row in rows:
        source = root / row['source']
        content = source.read_bytes()
        if hashlib.sha256(content).hexdigest() != row['sha256']:
            raise ValueError('Reader inventory hash mismatch: ' + row['source'])
        destination = work / row['path']
        if not destination.resolve().is_relative_to(work):
            raise ValueError('Unsafe reader path')
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    qmd = (root / 'chapters/24-differential-abundance.qmd').read_text()
    chunks = executable_chunks(qmd)
    assertions = r'''
stopifnot(
  identical(as.integer(detection_table$Detected), c(0L, 16L, 5L, 25L)),
  identical(as.integer(detection_table$NotDetected), c(59L, 35L, 54L, 26L)),
  abs(detection_table$Upper95[1] - 6.060889) < 0.0001,
  identical(dim(species_tss), c(212L, 110L)),
  identical(dim(pathway_tss), c(394L, 110L)),
  max(abs(colSums(species_tss) - 1)) < 1e-10,
  max(abs(colSums(pathway_tss) - 1)) < 1e-10,
  abs(min(retained_mass$SpeciesMass) - 0.5828294) < 1e-6,
  abs(median(retained_mass$SpeciesMass) - 0.9731469) < 1e-6,
  abs(min(retained_mass$PathwayMass) - 0.9961496) < 1e-6,
  identical(as.integer(aldex_q_comparison$Below005), c(4L, 4L)),
  identical(as.integer(summary_table$Value), c(110L, 212L, 394L, 144L, 7L, 0L, 50L, 1L, 4L, 0L)),
  all(!subset(species_results, Label %in% selected_species & Component == "Abundance")$Reportable),
  qr(model.matrix(~ disease + z_age + gender + z_BMI + z_log10_reads, primary_metadata))$rank == 6L,
  !anyDuplicated(primary_raw$subject_id)
)
write.table(detection_table, "detection-evidence.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(mass_summary, "denominator-evidence.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("DIAGNOSTICS_PASS: 15 scientific assertions\n")
'''
    (work / 'reader-diagnostics.R').write_text('\n\n'.join(chunks) + assertions)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with (work / 'run.log').open('w') as log:
        result = subprocess.run([a.rscript, '--vanilla', 'reader-diagnostics.R'], cwd=work,
                                env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
    ok = result.returncode == 0 and 'DIAGNOSTICS_PASS' in (work / 'run.log').read_text()
    report = {'status': 'passed' if ok else 'failed',
              'scope': 'public-input-only visible diagnostic code; separate full-model run required',
              'asset_count': len(rows), 'executed_chunks': len(chunks),
              'scientific_assertions': 15, 'exit_code': result.returncode}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))
    raise SystemExit(not ok)
