#!/usr/bin/env python3
"""Freeze compact, checksum-covered Article 47 phylogenomics evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

from article41_44_utils import freeze_compact_bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root, work, output = (
        args.project_root.resolve(), args.work_dir.resolve(), args.output_dir.resolve()
    )
    sources = [
        root / "scripts/article41_44_utils.py",
        root / "scripts/article42_44_validation_utils.py",
        root / "scripts/prepare_article47_phylogenomics.py",
        root / "scripts/run_article47_phylogenomics.py",
        root / "scripts/summarize_article47_phylogenomics.py",
        root / "scripts/plot_article47_phylogenomics.R",
        root / "scripts/freeze_article47_phylogenomics.py",
        root / "scripts/validate_article47_phylogenomics.py",
        root / "env/phylogeny.yml",
        root / "env/phylogeny-linux-64.lock",
    ]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Article 47 source files: " + ", ".join(missing))

    extra = [
        (work / "genome-ledger.tsv", Path("genome-ledger.tsv")),
        (work / "reference-request-ledger.tsv", Path("reference-request-ledger.tsv")),
        (work / "novelty-pretree-audit.tsv", Path("novelty-pretree-audit.tsv")),
        (work / "alignment-paths.tsv", Path("alignment-paths.tsv")),
        (
            root / "data/small/45-drep-dereplication-frozen/file-checksums.sha256",
            Path("upstream/article45-file-checksums.sha256"),
        ),
        (
            root / "data/small/46-gtdbtk-taxonomy-frozen/file-checksums.sha256",
            Path("upstream/article46-file-checksums.sha256"),
        ),
    ]
    for domain in ("bacteria", "archaea"):
        gtt = work / "gtotree" / domain
        iq = work / "iqtree" / domain
        extra.extend(
            [
                (gtt / "Aligned_SCGs_mod_names.faa", Path(f"raw/{domain}/alignment.faa")),
                (gtt / "Genomes_summary_info.tsv", Path(f"raw/{domain}/gtotree-genomes.tsv")),
                (gtt / "SCG_hit_counts.tsv", Path(f"raw/{domain}/gtotree-scg-hits.tsv")),
                (gtt / "citations.txt", Path(f"raw/{domain}/gtotree-citations.txt")),
                (gtt / "gtotree-runlog.txt", Path(f"raw/{domain}/gtotree-runlog.txt")),
                (iq / f"article47-{domain}.treefile", Path(f"raw/{domain}/iqtree.treefile")),
                (iq / f"article47-{domain}.iqtree", Path(f"raw/{domain}/iqtree-report.txt")),
                (iq / f"article47-{domain}.log", Path(f"raw/{domain}/iqtree.log")),
            ]
        )
    missing_extra = [str(source) for source, _ in extra if not source.is_file()]
    if missing_extra:
        raise FileNotFoundError("Missing Article 47 evidence files: " + ", ".join(missing_extra))

    freeze_compact_bundle(
        root=root,
        work=work,
        output=output,
        article=47,
        slug="novel-mag-phylogenomics",
        source_files=sources,
        extra_files=extra,
    )
    print(f"Frozen Article 47 evidence: {output}")


if __name__ == "__main__":
    main()
