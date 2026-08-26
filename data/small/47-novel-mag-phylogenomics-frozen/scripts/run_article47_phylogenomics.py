#!/usr/bin/env python3
"""Run GToTree alignments and seeded IQ-TREE inference for Article 47."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path

from article41_44_utils import write_tsv


def run_timed(label: str, command: list[str], work: Path, env: dict[str, str]) -> dict[str, object]:
    stdout_path = work / "logs" / f"{label}.stdout.log"
    stderr_path = work / "logs" / f"{label}.stderr.log"
    time_path = work / "logs" / f"{label}.time.txt"
    timed = ["/usr/bin/time", "-v", "-o", str(time_path), *command]
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(timed, stdout=stdout, stderr=stderr, env=env, check=False)
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
    result = subprocess.run(
        ["conda", "run", "-n", environment, *command],
        text=True, capture_output=True, env=env, check=False,
    )
    lines = [line.strip() for line in (result.stdout + "\n" + result.stderr).splitlines() if line.strip()]
    if command and command[0] == "hmmsearch":
        for line in lines:
            if line.startswith("# HMMER "):
                return line.removeprefix("# ")
    return lines[0] if lines else "VERSION_NOT_REPORTED"


def find_alignment(output: Path) -> Path:
    preferred = output / "Aligned_SCGs_mod_names.faa"
    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred
    candidates = [
        path for path in output.rglob("*")
        if path.is_file() and "aligned" in path.name.lower()
        and "scg" in path.name.lower() and "run_files" not in path.parts
        and path.suffix.lower() in {".faa", ".fa", ".fasta"}
    ]
    if len(candidates) != 1:
        raise ValueError(f"Expected one GToTree concatenated alignment under {output}, observed {candidates}")
    return candidates[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--environment", default="metagenome-phylogeny-2026.07")
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / ".article47-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article47_phylogenomics.py first")
    contract = json.loads((work / "run-contract.json").read_text(encoding="utf-8"))
    env = os.environ.copy()
    env.update({"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": str(args.threads)})
    commands: list[dict[str, object]] = []
    alignments: list[dict[str, object]] = []
    for domain_index, domain in enumerate(contract["domains"]):
        key = domain.lower()
        gtt_output = work / "gtotree" / key
        iq_output = work / "iqtree" / key
        if iq_output.exists():
            shutil.rmtree(iq_output)
        iq_output.mkdir(parents=True)
        gtt_command = [
            "conda", "run", "-n", args.environment, "GToTree",
            "-f", str(work / f"inputs/trees/{key}-query-fastas.txt"),
            "-a", str(work / f"inputs/trees/{key}-reference-accessions.txt"),
            "-H", domain,
            "-m", str(work / f"inputs/trees/{key}-labels.tsv"),
            "-o", str(gtt_output),
            "-F",
            "-N", "-k", "-d", "-P",
            "-n", "4", "-M", "4", "-j", "4",
        ]
        try:
            alignment = find_alignment(gtt_output)
            commands.append({
                "Label": f"gtotree-{key}",
                "ExitStatus": 0,
                "Command": shlex.join(gtt_command),
                "Stdout": str(work / "logs" / f"gtotree-{key}.stdout.log"),
                "Stderr": str(work / "logs" / f"gtotree-{key}.stderr.log"),
                "TimeLog": str(work / "logs" / f"gtotree-{key}.time.txt"),
            })
        except ValueError:
            if gtt_output.exists():
                shutil.rmtree(gtt_output)
            gtt_output.mkdir(parents=True)
            commands.append(run_timed(f"gtotree-{key}", gtt_command, work, env))
            alignment = find_alignment(gtt_output)
        prefix = iq_output / f"article47-{key}"
        iq_command = [
            "conda", "run", "-n", args.environment, "iqtree3",
            "-s", str(alignment), "-st", "AA", "-m", "MFP",
            "-B", "1000", "--alrt", "1000",
            "-T", "AUTO", "--threads-max", str(args.threads),
            "--seed", str(int(contract["seed"]) + domain_index),
            "--prefix", str(prefix), "--safe",
        ]
        commands.append(run_timed(f"iqtree-{key}", iq_command, work, env))
        alignments.append({
            "Domain": domain,
            "Alignment": str(alignment),
            "TreeFile": str(Path(f"{prefix}.treefile")),
            "IQTreeReport": str(Path(f"{prefix}.iqtree")),
            "Seed": int(contract["seed"]) + domain_index,
        })
    write_tsv(work / "command-log.tsv", commands)
    write_tsv(work / "alignment-paths.tsv", alignments)
    tools = [
        {"Tool": "GToTree", "Version": version(args.environment, ["GToTree", "-v"], env), "Role": "SCG extraction and alignment"},
        {"Tool": "IQ-TREE", "Version": version(args.environment, ["iqtree3", "--version"], env), "Role": "maximum-likelihood phylogeny"},
        {"Tool": "HMMER", "Version": version(args.environment, ["hmmsearch", "-h"], env), "Role": "SCG search"},
        {"Tool": "MUSCLE", "Version": version(args.environment, ["muscle", "-version"], env), "Role": "marker alignment"},
        {"Tool": "trimAl", "Version": version(args.environment, ["trimal", "--version"], env), "Role": "alignment trimming"},
        {"Tool": "Prodigal", "Version": version(args.environment, ["prodigal", "-v"], env), "Role": "gene prediction"},
    ]
    write_tsv(work / "tool-versions.tsv", tools)
    (work / ".article47-run-complete").write_text("complete\n", encoding="utf-8")
    print("Article 47 GToTree/IQ-TREE inference completed")


if __name__ == "__main__":
    main()
