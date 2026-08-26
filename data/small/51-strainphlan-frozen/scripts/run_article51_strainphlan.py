#!/usr/bin/env python3
"""Run a checksum-gated and explicitly seeded StrainPhlAn analysis for Article 51."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path

from article41_44_utils import parse_time, read_tsv, sha256, write_tsv


ASSET_FOLDERS = {
    "consensus marker profile": "consensus_markers",
    "clade marker FASTA": "clade_markers",
    "reference genome": "reference_genomes",
    "official precomputed output": "official_baseline",
    "MetaPhlAn metadata extracted by HTTP range": "database",
}


def run_timed(
    label: str, command: list[str], work: Path, env: dict[str, str]
) -> dict[str, object]:
    stdout_path = work / "logs" / f"{label}.stdout.log"
    stderr_path = work / "logs" / f"{label}.stderr.log"
    time_path = work / "logs" / f"{label}.time.txt"
    timed = ["/usr/bin/time", "-v", "-o", str(time_path), *command]
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        completed = subprocess.run(
            timed, stdout=stdout, stderr=stderr, env=env, check=False
        )
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed; inspect {stderr_path}")
    return {
        "Label": label,
        "ExitStatus": completed.returncode,
        "Command": shlex.join(command),
        "Stdout": str(stdout_path),
        "Stderr": str(stderr_path),
        "TimeLog": str(time_path),
    }


def version(environment: str, command: list[str], env: dict[str, str]) -> str:
    completed = subprocess.run(
        ["conda", "run", "-n", environment, *command],
        text=True,
        capture_output=True,
        env=env,
        check=False,
    )
    lines = [
        line.strip()
        for line in (completed.stdout + "\n" + completed.stderr).splitlines()
        if line.strip()
    ]
    return lines[0] if lines else "VERSION_NOT_REPORTED"


def verify_assets(work: Path) -> list[dict[str, object]]:
    manifest = read_tsv(work / "asset-manifest.tsv")
    if len(manifest) != 37:
        raise ValueError(f"Expected 37 official assets, observed {len(manifest)}")
    audit: list[dict[str, object]] = []
    for row in manifest:
        asset_type = row["AssetType"]
        if asset_type not in ASSET_FOLDERS:
            raise ValueError(f"Unrecognized asset type: {asset_type}")
        path = work / ASSET_FOLDERS[asset_type] / row["Name"]
        exists = path.is_file()
        observed_bytes = path.stat().st_size if exists else -1
        observed_sha = sha256(path) if exists else "MISSING"
        passed = (
            exists
            and observed_bytes == int(row["Bytes"])
            and observed_sha == row["SHA256"]
        )
        audit.append(
            {
                "AssetType": asset_type,
                "Name": row["Name"],
                "ExpectedBytes": int(row["Bytes"]),
                "ObservedBytes": observed_bytes,
                "ExpectedSHA256": row["SHA256"],
                "ObservedSHA256": observed_sha,
                "ChecksumPass": passed,
            }
        )
    write_tsv(work / "asset-check-audit.tsv", audit)
    if not all(row["ChecksumPass"] for row in audit):
        failed = [row["Name"] for row in audit if not row["ChecksumPass"]]
        raise ValueError(f"Refusing to deserialize failed assets: {failed}")
    return audit


def require_one(paths: list[Path], label: str) -> Path:
    paths = [path for path in paths if path.is_file() and path.stat().st_size > 0]
    if len(paths) != 1:
        raise ValueError(f"Expected one non-empty {label}, observed {paths}")
    return paths[0]


def replace_option(command: list[str], option: str, value: str) -> None:
    try:
        index = command.index(option)
    except ValueError as error:
        raise ValueError(f"Missing command option: {option}") from error
    command[index + 1] = value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--environment", default="metagenome-biobakery-2026.07")
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article51-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article51_strainphlan.py first")
    contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))
    verify_assets(work)

    output = work / "output"
    temporary = work / "tmp"
    runtime_tmp = work / "runtime-tmp"
    for path in (output, temporary, runtime_tmp):
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)
    (work / "logs").mkdir(exist_ok=True)
    for path in (work / "logs").glob("strainphlan-*"):
        path.unlink()

    env = os.environ.copy()
    env.update(
        {
            "TMPDIR": str(runtime_tmp),
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
    )
    commands: list[dict[str, object]] = []

    configuration = work / "phylophlan-seeded.cfg"
    config_command = [
        "conda",
        "run",
        "-n",
        args.environment,
        "phylophlan_write_config_file",
        "-o",
        str(configuration),
        "-d",
        "n",
        "--db_dna",
        "makeblastdb",
        "--map_dna",
        "blastn",
        "--msa",
        "mafft",
        "--tree1",
        "raxml",
        "--overwrite",
    ]
    commands.append(
        run_timed("strainphlan-config", config_command, work, env)
    )
    config_text = configuration.read_text(encoding="utf-8")
    if config_text.count("-p 1989") != 1:
        raise ValueError("Could not locate the single RAxML default seed in configuration")
    config_text = config_text.replace("-p 1989", f"-p {int(contract['seed'])}")
    configuration.write_text(config_text, encoding="utf-8")
    if f"params = -p {int(contract['seed'])} -m GTRCAT" not in config_text:
        raise ValueError("Explicit seeded RAxML configuration was not written")

    command = [
        "conda",
        "run",
        "-n",
        args.environment,
        "strainphlan",
        "--sample_list",
        str(work / "sample-list.txt"),
        "--reference_list",
        str(work / "reference-list.txt"),
        "-m",
        str(work / "clade_markers/t__SGB4933_group.fna"),
        "-d",
        str(work / "database/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl"),
        "-c",
        str(contract["target_clade"]),
        "-o",
        str(output),
        "-n",
        str(args.threads),
        "--trim_sequences",
        str(contract["trim_sequences"]),
        "--sample_with_n_markers",
        str(contract["sample_with_n_markers"]),
        "--sample_with_n_markers_perc",
        str(contract["sample_with_n_markers_perc"]),
        "--marker_in_n_samples_perc",
        str(contract["marker_in_n_samples_perc"]),
        "--sample_with_n_markers_after_filt",
        str(contract["sample_with_n_markers_after_filt"]),
        "--sample_with_n_markers_after_filt_perc",
        str(contract["sample_with_n_markers_after_filt_perc"]),
        "--breadth_thres",
        str(contract["breadth_thres"]),
        "--phylophlan_mode",
        str(contract["phylophlan_mode"]),
        "--phylophlan_configuration",
        str(configuration),
        "--tmp",
        str(temporary),
        "--debug",
    ]
    commands.append(run_timed("strainphlan-analysis", command, work, env))

    info_file = output / f"{contract['target_clade']}.info"
    polymorphic = output / f"{contract['target_clade']}.polymorphic"
    tree = require_one(list(output.glob("RAxML_bestTree.*.tre")), "RAxML tree")
    alignment = require_one(list(output.glob("*_concatenated.aln")), "concatenated alignment")
    for path in (info_file, polymorphic):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing required StrainPhlAn output: {path}")
    if not tree.read_text(encoding="utf-8", errors="replace").strip().endswith(";"):
        raise ValueError("RAxML output is not a complete Newick tree")
    if not alignment.read_text(encoding="utf-8", errors="replace").startswith(">"):
        raise ValueError("Concatenated alignment is not FASTA")

    legacy = contract["official_baseline_threshold_sensitivity"]
    legacy_output = work / "output-official-thresholds"
    legacy_temporary = work / "tmp-official-thresholds"
    for path in (legacy_output, legacy_temporary):
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)
    legacy_command = list(command)
    replace_option(legacy_command, "-o", str(legacy_output))
    replace_option(legacy_command, "--tmp", str(legacy_temporary))
    for option, key in (
        ("--sample_with_n_markers", "sample_with_n_markers"),
        ("--sample_with_n_markers_perc", "sample_with_n_markers_perc"),
        ("--marker_in_n_samples_perc", "marker_in_n_samples_perc"),
        ("--sample_with_n_markers_after_filt", "sample_with_n_markers_after_filt"),
        ("--sample_with_n_markers_after_filt_perc", "sample_with_n_markers_after_filt_perc"),
        ("--breadth_thres", "breadth_thres"),
    ):
        replace_option(legacy_command, option, str(legacy[key]))
    commands.append(
        run_timed("strainphlan-official-thresholds", legacy_command, work, env)
    )
    legacy_info = legacy_output / f"{contract['target_clade']}.info"
    legacy_polymorphic = legacy_output / f"{contract['target_clade']}.polymorphic"
    legacy_tree = require_one(
        list(legacy_output.glob("RAxML_bestTree.*.tre")),
        "official-threshold-sensitivity RAxML tree",
    )
    legacy_alignment = require_one(
        list(legacy_output.glob("*_concatenated.aln")),
        "official-threshold-sensitivity concatenated alignment",
    )
    for path in (legacy_info, legacy_polymorphic):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing threshold-sensitivity output: {path}")

    write_tsv(work / "command-log.tsv", commands)
    write_tsv(
        work / "resource-summary.tsv",
        [parse_time(Path(row["TimeLog"])) for row in commands],
    )
    write_tsv(
        work / "determinism-audit.tsv",
        [
            {
                "RandomComponent": "RAxML starting parsimony tree",
                "Control": f"-p {int(contract['seed'])}",
                "Evidence": str(configuration),
                "Pass": f"-p {int(contract['seed'])}" in config_text,
            },
            {
                "RandomComponent": "Python hash iteration",
                "Control": "PYTHONHASHSEED=0",
                "Evidence": "run environment",
                "Pass": env["PYTHONHASHSEED"] == "0",
            },
            {
                "RandomComponent": "MAFFT marker alignment",
                "Control": "--thread 1 in generated configuration",
                "Evidence": str(configuration),
                "Pass": "--thread 1" in config_text,
            },
        ],
    )
    write_tsv(
        work / "output-paths.tsv",
        [
            {
                "TargetClade": contract["target_clade"],
                "Info": str(info_file),
                "Polymorphic": str(polymorphic),
                "Alignment": str(alignment),
                "Tree": str(tree),
                "Configuration": str(configuration),
                "ConfigurationSHA256": sha256(configuration),
                "OfficialThresholdInfo": str(legacy_info),
                "OfficialThresholdPolymorphic": str(legacy_polymorphic),
                "OfficialThresholdAlignment": str(legacy_alignment),
                "OfficialThresholdTree": str(legacy_tree),
            }
        ],
    )
    write_tsv(
        work / "tool-versions.tsv",
        [
            {"Tool": "StrainPhlAn", "Version": version(args.environment, ["strainphlan", "-v"], env), "Role": "marker filtering and strain phylogeny workflow"},
            {"Tool": "MetaPhlAn", "Version": version(args.environment, ["metaphlan", "--version"], env), "Role": "Jan21 marker metadata interpretation"},
            {"Tool": "PhyloPhlAn", "Version": version(args.environment, ["phylophlan", "-v"], env), "Role": "marker alignment and phylogeny controller"},
            {"Tool": "RAxML", "Version": version(args.environment, ["raxmlHPC-PTHREADS-SSE3", "-v"], env), "Role": "seeded maximum-likelihood tree inference"},
            {"Tool": "MAFFT", "Version": version(args.environment, ["mafft", "--version"], env), "Role": "single-thread marker alignment"},
            {"Tool": "BLAST+", "Version": version(args.environment, ["blastn", "-version"], env), "Role": "reference-marker mapping"},
            {"Tool": "Python", "Version": version(args.environment, ["python", "--version"], env), "Role": "workflow runtime"},
        ],
    )
    (work / ".article51-run-complete").write_text("complete\n", encoding="utf-8")
    print(f"Article 51 completed: {tree.name}; {alignment.name}")


if __name__ == "__main__":
    main()
