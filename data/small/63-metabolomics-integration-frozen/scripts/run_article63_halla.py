#!/usr/bin/env python3
"""Run restart-safe HAllA discovery analyses for Article 63."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import pandas as pd


EXPECTED_MICROBES = 166
EXPECTED_METABOLITES = 153
EXPECTED_PAIRS = EXPECTED_MICROBES * EXPECTED_METABOLITES


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--halla-cli", type=Path, required=True)
    parser.add_argument("--r-library", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=63001)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def input_shape(path: Path) -> tuple[int, int]:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    if frame.index.duplicated().any() or frame.columns.duplicated().any():
        raise RuntimeError(f"Duplicate feature or sample IDs in {path}")
    if not frame.notna().all().all():
        raise RuntimeError(f"Missing values in {path}")
    return frame.shape


def output_is_complete(path: Path) -> bool:
    association_path = path / "all_associations.txt"
    cluster_path = path / "sig_clusters.txt"
    if not association_path.exists() or not cluster_path.exists():
        return False
    associations = pd.read_csv(association_path, sep="\t")
    required = {"X_features", "Y_features", "association", "p-values", "q-values"}
    return len(associations) == EXPECTED_PAIRS and required.issubset(associations.columns)


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    halla_cli = args.halla_cli.resolve()
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if not halla_cli.exists():
        raise FileNotFoundError(halla_cli)
    output_dir.mkdir(parents=True, exist_ok=True)

    environment = os.environ.copy()
    existing_r_libraries = environment.get("R_LIBS_USER", "")
    environment["R_LIBS_USER"] = ":".join(
        item for item in (str(args.r_library.resolve()), existing_r_libraries) if item
    )
    environment.setdefault("MPLCONFIGDIR", "/tmp/article63-mpl")
    Path(environment["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    version_result = subprocess.run(
        [str(halla_cli), "--version"], env=environment,
        text=True, capture_output=True, check=True
    )
    version_text = (version_result.stdout + version_result.stderr).strip()

    manifest: dict[str, object] = {
        "article": 63,
        "halla_version_output": version_text,
        "seed": args.seed,
        "threads": args.threads,
        "association": "spearman",
        "linkage": "average",
        "fdr_method": "fdr_bh",
        "fdr_alpha": 0.05,
        "fnr_threshold": 0.2,
        "pvalue_mode": (
            "HAllA 0.8.40 default analytic Spearman p-values; "
            "permutation settings are not invoked for metrics returning p-values"
        ),
        "expected_pairs_per_branch": EXPECTED_PAIRS,
        "branches": {},
    }

    for branch in ("adjusted", "raw"):
        x_path = input_dir / "halla" / f"prism-microbiome-{branch}.tsv"
        y_path = input_dir / "halla" / f"prism-metabolome-{branch}.tsv"
        x_shape = input_shape(x_path)
        y_shape = input_shape(y_path)
        if x_shape != (EXPECTED_MICROBES, 155) or y_shape != (EXPECTED_METABOLITES, 155):
            raise RuntimeError(f"Unexpected {branch} HAllA shapes: {x_shape}, {y_shape}")

        branch_dir = output_dir / branch
        if args.force and branch_dir.exists():
            shutil.rmtree(branch_dir)
        command = [
            str(halla_cli),
            "-x", str(x_path),
            "-y", str(y_path),
            "-o", str(branch_dir),
            "-m", "spearman",
            "--linkage_method", "average",
            "--fdr_alpha", "0.05",
            "--fdr_method", "fdr_bh",
            "--fnr_thresh", "0.2",
            "--rank_cluster", "best",
            "--seed", str(args.seed + (0 if branch == "adjusted" else 1)),
            "--num_threads", str(args.threads),
            "--dont_copy",
            "--no_hallagram",
            "--no_progress",
        ]
        log_path = output_dir / f"{branch}-halla.log"
        if not output_is_complete(branch_dir):
            completed = subprocess.run(
                command, env=environment, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False
            )
            log_path.write_text(completed.stdout, encoding="utf-8")
            if completed.returncode != 0:
                raise RuntimeError(
                    f"HAllA {branch} failed with exit code {completed.returncode}; "
                    f"see {log_path}"
                )
        if not output_is_complete(branch_dir):
            raise RuntimeError(f"Incomplete HAllA output in {branch_dir}")

        associations = pd.read_csv(branch_dir / "all_associations.txt", sep="\t")
        clusters = pd.read_csv(branch_dir / "sig_clusters.txt", sep="\t")
        if set(associations["X_features"]) != {
            f"MB{index:03d}" for index in range(1, EXPECTED_MICROBES + 1)
        }:
            raise RuntimeError(f"Unexpected microbiome feature universe in {branch}")
        if set(associations["Y_features"]) != {
            f"MT{index:03d}" for index in range(1, EXPECTED_METABOLITES + 1)
        }:
            raise RuntimeError(f"Unexpected metabolome feature universe in {branch}")
        manifest["branches"][branch] = {
            "microbiome_input": str(x_path),
            "metabolome_input": str(y_path),
            "microbiome_sha256": sha256(x_path),
            "metabolome_sha256": sha256(y_path),
            "associations": len(associations),
            "marginal_q_lt_0_05": int(associations["q-values"].lt(0.05).sum()),
            "significant_hierarchical_blocks": len(clusters),
            "command": command,
        }

    (output_dir / "run-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
