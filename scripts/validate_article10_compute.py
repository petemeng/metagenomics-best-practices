#!/usr/bin/env python3
"""Validate Article 10 compute planning, local arrays, and evidence boundaries.

This validator intentionally separates three planes:

1. a native-Linux resource and deterministic local-task smoke test;
2. syntax and contract checks for an HPC SLURM/Apptainer template;
3. control planes that are not attached to this host and remain NOT_RUN.

The real input is the checksum-identified Article 08 MOCK1 evidence.  The
5,000-record tasks exercise staging, accounting, atomic output, and resume
semantics only; they are not an assembler or full-FASTQ performance benchmark.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Sequence


EXPECTED_SHA256 = {
    "08-ena-fastq-sources.tsv": (
        "96638ae0d16ce5953760b596887c9f8f443c1dc33e0f82bef2e6ac040e69d962"
    ),
    "08-read-prefix-metrics.tsv": (
        "6c426e0c89a95b1c1cc7ec631629e3635a7fd83c61c90d4086801b28495b5a61"
    ),
    "10-job-array.tsv": (
        "368baec61944a605305723971f0442f8fc568ce6f7107452d44299e4cb76776f"
    ),
    "10-runtime-contract.tsv": (
        "dad2c5154b3e1f703bab76698adc04607668a4ba0bcc136b53191143cab62e6a"
    ),
    "10-container-smoke.log": (
        "aeb469f3084bb8062d8c9538a8585deb0aa6ec90db1ff3f7563e5ac1662d6750"
    ),
}
EXPECTED_PLATFORM_BYTES = {
    "Illumina": 3_845_199_421,
    "ONT": 3_117_261_341,
    "PacBio": 3_982_506_052,
}
EXPECTED_TOTAL_BYTES = 10_944_966_814
EXPECTED_METRICS_ROWS = 15_000
EXPECTED_ROWS_PER_TASK = 5_000
EXPECTED_APPTAINER_RELEASE = "1.5.2"
EXPECTED_OCI_DIGEST = (
    "sha256:eafc1edb577d2e9b458664a15f23ea1c370214193226069eb22921169fc7e43f"
)
EXPECTED_SIF_SHA256 = (
    "5dd8b9d84e1dfec614fdece44123833d2acc6330a44538bf6b9ab2f0b76902d8"
)
SMOKE_SEED = 20260719


@dataclass(frozen=True)
class CommandResult:
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    return parser.parse_args()


def run_command(
    command: Sequence[str],
    *,
    env: dict[str, str] | None = None,
    timeout: int = 300,
) -> CommandResult:
    merged_env = {**os.environ, "LC_ALL": "C", "LANG": "C"}
    if env:
        merged_env.update(env)
    completed = subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=merged_env,
    )
    return CommandResult(
        command=tuple(command),
        returncode=completed.returncode,
        stdout=completed.stdout.strip(),
        stderr=completed.stderr.strip(),
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path,
    rows: Iterable[dict[str, object]],
    fields: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def redact(text: str, project_root: Path) -> str:
    redacted = text
    replacements = {
        str(project_root.resolve()): "<PROJECT_ROOT>",
        str(Path.home()): "<HOME>",
    }
    for variable in (
        "http_proxy",
        "https_proxy",
        "HTTP_PROXY",
        "HTTPS_PROXY",
        "ALL_PROXY",
        "all_proxy",
    ):
        value = os.environ.get(variable)
        if value:
            replacements[value] = "<PROXY_REDACTED>"
    for source, replacement in sorted(
        replacements.items(), key=lambda item: len(item[0]), reverse=True
    ):
        redacted = redacted.replace(source, replacement)
    redacted = re.sub(r"/tmp/article10-[^/\s]+", "<TASK_SCRATCH>", redacted)
    redacted = re.sub(
        r"(https?://)([^/@\s]+)@",
        r"\1<CREDENTIALS_REDACTED>@",
        redacted,
    )
    return redacted


def add_check(
    rows: list[dict[str, object]],
    category: str,
    check: str,
    expected: object,
    observed: object,
    status: str,
    scope: str,
) -> None:
    rows.append(
        {
            "category": category,
            "check": check,
            "expected": expected,
            "observed": observed,
            "status": status,
            "scope": scope,
        }
    )


def os_release() -> str:
    path = Path("/etc/os-release")
    if not path.exists():
        return platform.system()
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key] = value.strip().strip('"')
    return values.get("PRETTY_NAME", values.get("NAME", platform.system()))


def memory_snapshot() -> dict[str, int]:
    values: dict[str, int] = {}
    path = Path("/proc/meminfo")
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^([A-Za-z_()]+):\s+(\d+)\s+kB$", line)
        if match:
            values[match.group(1)] = int(match.group(2)) * 1024
    return values


def physical_core_count() -> int:
    path = Path("/proc/cpuinfo")
    if not path.is_file():
        return os.cpu_count() or 1
    pairs: set[tuple[str, str]] = set()
    physical_id = ""
    core_id = ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines() + [""]:
        if line.startswith("physical id"):
            physical_id = line.split(":", 1)[1].strip()
        elif line.startswith("core id"):
            core_id = line.split(":", 1)[1].strip()
        elif not line.strip():
            if physical_id or core_id:
                pairs.add((physical_id, core_id))
            physical_id = ""
            core_id = ""
    return len(pairs) if pairs else (os.cpu_count() or 1)


def filesystem_type(path: Path) -> str:
    findmnt = shutil.which("findmnt")
    if not findmnt:
        return "unknown"
    result = run_command([findmnt, "-n", "-o", "FSTYPE", "--target", str(path)])
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def host_inventory(project_root: Path) -> list[dict[str, object]]:
    memory = memory_snapshot()
    stat = os.statvfs(project_root)
    total_bytes = stat.f_blocks * stat.f_frsize
    available_bytes = stat.f_bavail * stat.f_frsize
    rows = [
        {
            "metric": "operating_system",
            "value": os_release(),
            "unit": "text",
            "evidence_scope": "native Linux read-only observation",
        },
        {
            "metric": "kernel_release",
            "value": platform.release(),
            "unit": "text",
            "evidence_scope": "native Linux read-only observation",
        },
        {
            "metric": "architecture",
            "value": platform.machine(),
            "unit": "text",
            "evidence_scope": "native Linux read-only observation",
        },
        {
            "metric": "logical_cpus",
            "value": os.cpu_count() or 1,
            "unit": "count",
            "evidence_scope": "native Linux read-only observation",
        },
        {
            "metric": "physical_cores",
            "value": physical_core_count(),
            "unit": "count",
            "evidence_scope": "native Linux read-only observation",
        },
        {
            "metric": "memory_total",
            "value": memory.get("MemTotal", 0),
            "unit": "bytes",
            "evidence_scope": "native Linux read-only observation",
        },
        {
            "metric": "memory_available_snapshot",
            "value": memory.get("MemAvailable", 0),
            "unit": "bytes",
            "evidence_scope": "dynamic snapshot; not a frozen capacity promise",
        },
        {
            "metric": "swap_total",
            "value": memory.get("SwapTotal", 0),
            "unit": "bytes",
            "evidence_scope": "report only; swap is not a RAM substitute",
        },
        {
            "metric": "project_filesystem_type",
            "value": filesystem_type(project_root),
            "unit": "text",
            "evidence_scope": "mount path redacted",
        },
        {
            "metric": "project_filesystem_total",
            "value": total_bytes,
            "unit": "bytes",
            "evidence_scope": "mount path redacted",
        },
        {
            "metric": "project_filesystem_available_snapshot",
            "value": available_bytes,
            "unit": "bytes",
            "evidence_scope": "dynamic snapshot; mount path redacted",
        },
    ]
    return rows


def parse_elapsed(value: str) -> float:
    parts = value.split(":")
    try:
        if len(parts) == 2:
            return int(parts[0]) * 60 + float(parts[1])
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    except ValueError:
        return -1.0
    return -1.0


def parse_gnu_time(stderr: str) -> dict[str, float]:
    metrics: dict[str, float] = {}
    patterns = {
        "max_rss_kib": r"Maximum resident set size \(kbytes\):\s*(\d+)",
        "filesystem_inputs": r"File system inputs:\s*(\d+)",
        "filesystem_outputs": r"File system outputs:\s*(\d+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, stderr)
        metrics[key] = float(match.group(1)) if match else -1.0
    elapsed_match = re.search(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9:.]+)",
        stderr,
    )
    metrics["elapsed_seconds"] = (
        parse_elapsed(elapsed_match.group(1)) if elapsed_match else -1.0
    )
    return metrics


def input_footprint(
    source_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    grouped: dict[str, int] = defaultdict(int)
    accessions: dict[str, set[str]] = defaultdict(set)
    for row in source_rows:
        platform_key = row["PlatformKey"]
        grouped[platform_key] += int(row["ENABytes"])
        accessions[platform_key].add(row["RunAccession"])
    rows: list[dict[str, object]] = []
    for platform_key in ("Illumina", "ONT", "PacBio"):
        value = grouped[platform_key]
        rows.append(
            {
                "platform": platform_key,
                "run_accessions": ",".join(sorted(accessions[platform_key])),
                "compressed_bytes": value,
                "compressed_gb_decimal": round(value / 1e9, 6),
                "compressed_gib_binary": round(value / (1024**3), 6),
                "transfer_minutes_at_1_gbps": round(value * 8 / 1e9 / 60, 3),
                "evidence_scope": "ENA-reported compressed source bytes",
            }
        )
    return rows


def command_version(command: str) -> str:
    path = shutil.which(command)
    if not path:
        return "NOT_FOUND"
    result = run_command([path, "version"])
    if result.returncode != 0:
        result = run_command([path, "--version"])
    text = (result.stdout or result.stderr).splitlines()
    return text[0].strip() if text else f"exit={result.returncode}"


def save_publication_figure(fig, stem: Path) -> list[str]:
    outputs: list[str] = []
    for suffix in (".pdf", ".png", ".tiff"):
        path = stem.with_suffix(suffix)
        path.parent.mkdir(parents=True, exist_ok=True)
        kwargs: dict[str, object] = {"bbox_inches": "tight"}
        if suffix == ".png":
            kwargs["dpi"] = 350
        elif suffix == ".tiff":
            kwargs["dpi"] = 350
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **kwargs)
        outputs.append(path.name)
    return outputs


def configure_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 13,
            "axes.labelsize": 10,
            "axes.edgecolor": "#314452",
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.color": "#314452",
            "ytick.color": "#314452",
            "text.color": "#24323F",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    return plt


def plot_input_footprint(rows: list[dict[str, object]], figure_dir: Path) -> list[str]:
    plt = configure_matplotlib()
    labels = ["Illumina paired", "ONT R9", "PacBio CCS"]
    values = [float(row["compressed_gb_decimal"]) for row in rows]
    colors = ["#20639B", "#E07A5F", "#2A9D8F"]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    bars = ax.barh(labels, values, color=colors, height=0.58)
    ax.invert_yaxis()
    ax.set_xlabel("Compressed source size (GB, decimal)")
    fig.suptitle(
        "One public benchmark sample: compressed FASTQ footprint",
        x=0.13,
        y=0.975,
        ha="left",
        fontsize=13,
    )
    fig.text(
        0.13,
        0.915,
        f"ENA-reported bytes • {sum(values):.2f} GB total • not a RAM forecast",
        fontsize=9,
        color="#5B6870",
    )
    ax.grid(axis="x", color="#DDE5EA", linewidth=0.7)
    ax.set_axisbelow(True)
    for bar, value in zip(bars, values):
        ax.text(
            value + max(values) * 0.025,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.2f} GB",
            va="center",
            fontsize=9,
        )
    ax.set_xlim(0, max(values) * 1.28)
    fig.tight_layout(rect=(0, 0, 1, 0.84))
    outputs = save_publication_figure(
        fig, figure_dir / "10-input-resource-budget"
    )
    plt.close(fig)
    return outputs


def draw_box(ax, x: float, y: float, width: float, height: float, title: str, body: str, color: str):
    from matplotlib.patches import FancyBboxPatch

    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        linewidth=1.1,
        edgecolor=color,
        facecolor="#FFFFFF",
    )
    ax.add_patch(patch)
    ax.text(x + 0.025, y + height - 0.06, title, weight="bold", fontsize=9.5)
    ax.text(
        x + 0.025,
        y + height - 0.135,
        body,
        fontsize=8.2,
        va="top",
        linespacing=1.35,
        color="#42525D",
    )


def plot_control_loop(figure_dir: Path) -> list[str]:
    plt = configure_matplotlib()
    from matplotlib.patches import FancyArrowPatch

    fig, ax = plt.subplots(figsize=(8.3, 5.1))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Resource requests are a measured control loop", loc="left", pad=14)
    ax.text(
        0,
        0.955,
        "Calibrate within the same tool, database, read set, and assembly design",
        fontsize=9.3,
        color="#5B6870",
    )
    boxes = [
        (0.03, 0.59, "1  Pilot", "One representative task\nFixed inputs + versions", "#20639B"),
        (0.285, 0.59, "2  Measure", "MaxRSS • Elapsed\nRead/write • exit state", "#2A9D8F"),
        (0.54, 0.59, "3  Add headroom", "Declare the margin\nDo not hide OOM retries", "#E9A23B"),
        (0.795, 0.59, "4  Bound array", "CPU • memory • time\nConcurrency cap", "#C44536"),
        (0.41, 0.24, "5  Freeze evidence", "JobID + sacct table\nImage + database digest", "#6C5CE7"),
    ]
    for x, y, title, body, color in boxes:
        draw_box(ax, x, y, 0.18, 0.22, title, body, color)
    arrow_pairs = [
        ((0.21, 0.70), (0.28, 0.70)),
        ((0.465, 0.70), (0.535, 0.70)),
        ((0.72, 0.70), (0.79, 0.70)),
        ((0.885, 0.58), (0.59, 0.45)),
        ((0.41, 0.35), (0.15, 0.58)),
    ]
    for start, end in arrow_pairs:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=12,
                linewidth=1.25,
                color="#6B7A86",
                connectionstyle="arc3,rad=0.0",
            )
        )
    ax.text(
        0.5,
        0.08,
        "A literature number is a prior. Your measured pilot is the request contract.",
        ha="center",
        fontsize=10,
        weight="bold",
        color="#24323F",
    )
    fig.tight_layout()
    outputs = save_publication_figure(
        fig, figure_dir / "10-resource-control-loop"
    )
    plt.close(fig)
    return outputs


def plot_resume_contract(
    completed: int,
    skipped: int,
    figure_dir: Path,
) -> list[str]:
    plt = configure_matplotlib()
    from matplotlib.patches import FancyArrowPatch

    fig, ax = plt.subplots(figsize=(8.1, 4.7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Restart-safe arrays publish only complete outputs", loc="left", pad=12)
    draw_box(
        ax,
        0.04,
        0.51,
        0.24,
        0.27,
        f"First local smoke: {completed} completed",
        "Task-specific scratch\nValidate JSON\nNo sentinel on signal",
        "#20639B",
    )
    draw_box(
        ax,
        0.38,
        0.51,
        0.24,
        0.27,
        "Atomic publish",
        "Move final output\nWrite SHA-256 sentinel\nKeep per-task logs",
        "#2A9D8F",
    )
    draw_box(
        ax,
        0.72,
        0.51,
        0.24,
        0.27,
        f"Second local smoke: {skipped} skipped",
        "Verify output checksum\nSkip valid task\nRecompute missing task only",
        "#6C5CE7",
    )
    for start, end in [((0.28, 0.645), (0.38, 0.645)), ((0.62, 0.645), (0.72, 0.645))]:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.35,
                color="#6B7A86",
            )
        )
    ax.text(
        0.5,
        0.29,
        "SLURM control plane: NOT RUN on this native-Linux host",
        ha="center",
        fontsize=10,
        weight="bold",
        color="#9A6700",
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "#FFF4D6",
            "edgecolor": "#E9A23B",
        },
    )
    ax.text(
        0.5,
        0.13,
        "The same body was syntax-checked and executed with RUN_MODE=native-smoke.",
        ha="center",
        fontsize=9,
        color="#42525D",
    )
    fig.tight_layout()
    outputs = save_publication_figure(
        fig, figure_dir / "10-restart-safe-array"
    )
    plt.close(fig)
    return outputs


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else project_root / "results" / "10-computing-hpc-cloud"
    )
    figure_dir = (
        args.figure_dir.resolve()
        if args.figure_dir
        else project_root / "figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    small = project_root / "data" / "small"
    source_path = small / "08-ena-fastq-sources.tsv"
    metrics_path = small / "08-read-prefix-metrics.tsv"
    task_table_path = small / "10-job-array.tsv"
    runtime_contract_path = small / "10-runtime-contract.tsv"
    smoke_log_path = small / "10-container-smoke.log"
    slurm_script = project_root / "scripts" / "10_slurm_array_smoke.sh"
    task_script = project_root / "scripts" / "article10_task.py"

    checks: list[dict[str, object]] = []
    log_lines = [
        "Article 10 compute validation",
        f"date={date.today().isoformat()}",
        "scope=native-linux-local-smoke",
        "slurm_control_plane=NOT_RUN unless commands are attached",
        "performance_scope=workflow-plumbing-smoke-not-assembler-benchmark",
    ]

    for path in (
        source_path,
        metrics_path,
        task_table_path,
        runtime_contract_path,
        smoke_log_path,
        slurm_script,
        task_script,
    ):
        add_check(
            checks,
            "input",
            f"{path.name}_exists",
            True,
            path.is_file(),
            "PASS" if path.is_file() else "FAIL",
            "repository input",
        )
    if any(row["status"] == "FAIL" for row in checks):
        raise SystemExit("Article 10 required inputs are missing")

    for name, expected_hash in EXPECTED_SHA256.items():
        path = small / name
        observed_hash = sha256(path)
        add_check(
            checks,
            "input",
            f"{name}_sha256",
            expected_hash,
            observed_hash,
            "PASS" if observed_hash == expected_hash else "FAIL",
            "byte-locked input",
        )

    source_rows = read_tsv(source_path)
    task_rows = read_tsv(task_table_path)
    runtime_rows = read_tsv(runtime_contract_path)
    footprint_rows = input_footprint(source_rows)
    platform_bytes = {
        str(row["platform"]): int(row["compressed_bytes"]) for row in footprint_rows
    }
    total_bytes = sum(platform_bytes.values())
    add_check(
        checks,
        "real-data",
        "ena_source_rows",
        4,
        len(source_rows),
        "PASS" if len(source_rows) == 4 else "FAIL",
        "Article 08 checksum-identified sources",
    )
    for key, expected in EXPECTED_PLATFORM_BYTES.items():
        observed = platform_bytes.get(key, 0)
        add_check(
            checks,
            "real-data",
            f"{key.lower()}_compressed_bytes",
            expected,
            observed,
            "PASS" if observed == expected else "FAIL",
            "ENA-reported compressed bytes",
        )
    add_check(
        checks,
        "real-data",
        "total_compressed_bytes",
        EXPECTED_TOTAL_BYTES,
        total_bytes,
        "PASS" if total_bytes == EXPECTED_TOTAL_BYTES else "FAIL",
        "four public FASTQ files",
    )
    add_check(
        checks,
        "real-data",
        "job_array_task_count",
        3,
        len(task_rows),
        "PASS" if len(task_rows) == 3 else "FAIL",
        "three platform metric strata",
    )
    task_ids = [int(row["TaskID"]) for row in task_rows]
    add_check(
        checks,
        "real-data",
        "job_array_ids",
        "1,2,3",
        ",".join(map(str, task_ids)),
        "PASS" if task_ids == [1, 2, 3] else "FAIL",
        "deterministic SLURM_ARRAY_TASK_ID mapping",
    )
    expected_task_rows = sum(int(row["ExpectedRows"]) for row in task_rows)
    add_check(
        checks,
        "real-data",
        "expected_metric_rows",
        EXPECTED_METRICS_ROWS,
        expected_task_rows,
        "PASS" if expected_task_rows == EXPECTED_METRICS_ROWS else "FAIL",
        "three deterministic prefixes",
    )

    write_tsv(
        output_dir / "input-footprint.tsv",
        footprint_rows,
        [
            "platform",
            "run_accessions",
            "compressed_bytes",
            "compressed_gb_decimal",
            "compressed_gib_binary",
            "transfer_minutes_at_1_gbps",
            "evidence_scope",
        ],
    )

    inventory_rows = host_inventory(project_root)
    inventory = {str(row["metric"]): row["value"] for row in inventory_rows}
    write_tsv(
        output_dir / "host-inventory.tsv",
        inventory_rows,
        ["metric", "value", "unit", "evidence_scope"],
    )
    host_checks = [
        (
            "linux_host",
            "Linux",
            platform.system(),
            platform.system() == "Linux",
            "native execution plane",
        ),
        (
            "logical_cpu_floor",
            ">=2",
            inventory["logical_cpus"],
            int(inventory["logical_cpus"]) >= 2,
            "smoke-task capacity only",
        ),
        (
            "memory_floor_bytes",
            ">=1 GiB",
            inventory["memory_total"],
            int(inventory["memory_total"]) >= 1024**3,
            "smoke-task capacity only",
        ),
        (
            "filesystem_available_floor_bytes",
            f">={metrics_path.stat().st_size * 4}",
            inventory["project_filesystem_available_snapshot"],
            int(inventory["project_filesystem_available_snapshot"])
            >= metrics_path.stat().st_size * 4,
            "dynamic local snapshot",
        ),
    ]
    for name, expected, observed, passed, scope in host_checks:
        add_check(
            checks,
            "host",
            name,
            expected,
            observed,
            "PASS" if passed else "FAIL",
            scope,
        )

    runtime_contract = {
        row["ContractKey"]: row["ExpectedValue"] for row in runtime_rows
    }
    runtime_expectations = {
        "recommended_apptainer_release": EXPECTED_APPTAINER_RELEASE,
        "oci_reference": f"docker://alpine@{EXPECTED_OCI_DIGEST}",
        "oci_platform": "linux/amd64",
        "legacy_sif_sha256": EXPECTED_SIF_SHA256,
        "legacy_sif_exec": "EXPECTED_FAIL",
        "local_slurm_control_plane": "NOT_FOUND",
    }
    for key, expected in runtime_expectations.items():
        observed = runtime_contract.get(key, "")
        add_check(
            checks,
            "container-contract",
            key,
            expected,
            observed,
            "PASS" if observed == expected else "FAIL",
            "frozen one-time evidence",
        )
    smoke_text = smoke_log_path.read_text(encoding="utf-8")
    for token in (
        "apptainer_release=1.5.2",
        EXPECTED_OCI_DIGEST,
        "docker_payload_smoke=PASS",
        "version=3.7.2",
        f"sif_sha256={EXPECTED_SIF_SHA256}",
        "sif_exec=EXPECTED_FAIL",
        "sbatch=NOT_FOUND",
        "sacct=NOT_FOUND",
    ):
        add_check(
            checks,
            "container-contract",
            f"log_token_{hashlib.sha256(token.encode()).hexdigest()[:10]}",
            token,
            token if token in smoke_text else "missing",
            "PASS" if token in smoke_text else "FAIL",
            "normalized one-time container/scheduler log",
        )

    apptainer_version = command_version("apptainer")
    singularity_version = command_version("singularity")
    sbatch_path = shutil.which("sbatch")
    sacct_path = shutil.which("sacct")
    shellcheck_path = shutil.which("shellcheck")
    runtime_audit_rows = [
        {
            "component": "recommended_apptainer",
            "expected": EXPECTED_APPTAINER_RELEASE,
            "observed": runtime_contract["recommended_apptainer_release"],
            "status": "PASS",
            "scope": "current official release contract; not locally executed",
        },
        {
            "component": "current_apptainer_cli",
            "expected": "cluster-provided runtime",
            "observed": apptainer_version,
            "status": "PASS" if apptainer_version != "NOT_FOUND" else "NOT_RUN",
            "scope": "native host command availability",
        },
        {
            "component": "current_singularity_cli",
            "expected": "report if present; do not treat as current Apptainer",
            "observed": singularity_version,
            "status": "PASS" if singularity_version != "NOT_FOUND" else "NOT_RUN",
            "scope": "legacy compatibility observation",
        },
        {
            "component": "legacy_sif_exec",
            "expected": "EXPECTED_FAIL",
            "observed": runtime_contract["legacy_sif_exec"],
            "status": "EXPECTED_FAIL",
            "scope": "frozen local installation boundary",
        },
        {
            "component": "slurm_sbatch",
            "expected": "cluster control plane",
            "observed": "FOUND" if sbatch_path else "NOT_FOUND",
            "status": "PASS" if sbatch_path else "NOT_RUN",
            "scope": "native host command availability",
        },
        {
            "component": "slurm_sacct",
            "expected": "cluster accounting plane",
            "observed": "FOUND" if sacct_path else "NOT_FOUND",
            "status": "PASS" if sacct_path else "NOT_RUN",
            "scope": "native host command availability",
        },
    ]
    write_tsv(
        output_dir / "container-runtime-audit.tsv",
        runtime_audit_rows,
        ["component", "expected", "observed", "status", "scope"],
    )
    for row in runtime_audit_rows:
        add_check(
            checks,
            "runtime",
            str(row["component"]),
            row["expected"],
            row["observed"],
            str(row["status"]),
            str(row["scope"]),
        )

    shell_result = run_command(["bash", "-n", str(slurm_script)])
    add_check(
        checks,
        "scheduler-template",
        "bash_syntax",
        "exit 0",
        f"exit {shell_result.returncode}",
        "PASS" if shell_result.returncode == 0 else "FAIL",
        "syntax only; not a SLURM submission",
    )
    script_text = slurm_script.read_text(encoding="utf-8")
    required_tokens = (
        "#SBATCH --array=1-3%3",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=2",
        "#SBATCH --mem=4G",
        "#SBATCH --time=00:10:00",
        "%A_%a",
        "SLURM_TMPDIR",
        "trap cleanup EXIT",
        "trap on_signal USR1 TERM INT",
        "sha256sum",
        "apptainer exec",
        "--cleanenv",
        "--bind",
        "--pwd",
        "RUN_MODE",
        "native-smoke",
        "mv \"${tmp_final}\" \"${final_json}\"",
        "--seed \"${SEED}\"",
    )
    scheduler_rows: list[dict[str, object]] = []
    for token in required_tokens:
        present = token in script_text
        row = {
            "check": token,
            "expected": "present",
            "observed": "present" if present else "missing",
            "status": "PASS" if present else "FAIL",
            "scope": "static SLURM/Apptainer template contract",
        }
        scheduler_rows.append(row)
        add_check(
            checks,
            "scheduler-template",
            f"token_{hashlib.sha256(token.encode()).hexdigest()[:10]}",
            "present",
            row["observed"],
            row["status"],
            row["scope"],
        )
    if shellcheck_path:
        shellcheck_result = run_command([shellcheck_path, str(slurm_script)])
        shellcheck_status = "PASS" if shellcheck_result.returncode == 0 else "FAIL"
        shellcheck_observed = f"exit {shellcheck_result.returncode}"
        log_lines.extend(
            [
                "[shellcheck]",
                redact(shellcheck_result.stdout, project_root),
                redact(shellcheck_result.stderr, project_root),
            ]
        )
    else:
        shellcheck_status = "NOT_RUN"
        shellcheck_observed = "NOT_FOUND"
    scheduler_rows.append(
        {
            "check": "shellcheck",
            "expected": "run when installed",
            "observed": shellcheck_observed,
            "status": shellcheck_status,
            "scope": "optional static linter",
        }
    )
    add_check(
        checks,
        "scheduler-template",
        "shellcheck",
        "run when installed",
        shellcheck_observed,
        shellcheck_status,
        "optional static linter",
    )
    write_tsv(
        output_dir / "scheduler-template-audit.tsv",
        scheduler_rows,
        ["check", "expected", "observed", "status", "scope"],
    )

    time_command = shutil.which("time")
    if time_command != "/usr/bin/time" and Path("/usr/bin/time").is_file():
        time_command = "/usr/bin/time"
    if not time_command:
        raise SystemExit("GNU time is required for the Article 10 telemetry smoke")

    array_dir = output_dir / "array-output"
    if array_dir.exists():
        shutil.rmtree(array_dir)
    array_dir.mkdir(parents=True, exist_ok=True)
    telemetry_rows: list[dict[str, object]] = []
    resume_rows: list[dict[str, object]] = []
    first_completed = 0
    second_skipped = 0
    observed_rows = 0

    common_env = {
        "PROJECT_ROOT": str(project_root),
        "TASK_TABLE": str(task_table_path),
        "METRICS": str(metrics_path),
        "OUT_DIR": str(array_dir),
        "RUN_MODE": "native-smoke",
        "PYTHONHASHSEED": str(SMOKE_SEED),
    }
    for task in task_rows:
        task_id = task["TaskID"]
        task_env = {**common_env, "SLURM_ARRAY_TASK_ID": task_id}
        result = run_command(
            [time_command, "-v", "bash", str(slurm_script)],
            env=task_env,
        )
        log_lines.extend(
            [
                f"[first-run task={task_id}]",
                redact(result.stdout, project_root),
                redact(result.stderr, project_root),
            ]
        )
        completed = result.returncode == 0 and "ACTION=COMPLETED" in result.stdout
        first_completed += int(completed)
        time_metrics = parse_gnu_time(result.stderr)
        output_path = (
            array_dir
            / f"{int(task_id):02d}-{task['PlatformKey'].lower()}.json"
        )
        done_path = Path(f"{output_path}.done")
        payload: dict[str, object] = {}
        if output_path.is_file():
            payload = json.loads(output_path.read_text(encoding="utf-8"))
            observed_rows += int(payload.get("rows", 0))
        output_hash = sha256(output_path) if output_path.is_file() else ""
        sentinel_hash = ""
        if done_path.is_file():
            sentinel_hash = done_path.read_text(encoding="utf-8").split()[0]
        valid_sentinel = bool(output_hash) and output_hash == sentinel_hash
        telemetry_rows.append(
            {
                "task_id": task_id,
                "platform": task["PlatformKey"],
                "run_accession": task["RunAccession"],
                "rows": payload.get("rows", 0),
                "max_rss_kib": int(time_metrics["max_rss_kib"]),
                "elapsed_seconds": time_metrics["elapsed_seconds"],
                "filesystem_inputs": int(time_metrics["filesystem_inputs"]),
                "filesystem_outputs": int(time_metrics["filesystem_outputs"]),
                "output_sha256": output_hash,
                "status": "PASS"
                if completed
                and valid_sentinel
                and int(payload.get("rows", 0)) == EXPECTED_ROWS_PER_TASK
                else "FAIL",
                "scope": "native local workflow smoke; not an assembler benchmark",
            }
        )
        resume_rows.append(
            {
                "task_id": task_id,
                "phase": "first_run",
                "expected_action": "COMPLETED",
                "observed_action": "COMPLETED" if completed else "FAILED",
                "output_checksum_valid": valid_sentinel,
                "status": "PASS" if completed and valid_sentinel else "FAIL",
            }
        )

    for task in task_rows:
        task_id = task["TaskID"]
        task_env = {**common_env, "SLURM_ARRAY_TASK_ID": task_id}
        result = run_command(["bash", str(slurm_script)], env=task_env)
        log_lines.extend(
            [
                f"[second-run task={task_id}]",
                redact(result.stdout, project_root),
                redact(result.stderr, project_root),
            ]
        )
        skipped = result.returncode == 0 and "ACTION=SKIPPED" in result.stdout
        second_skipped += int(skipped)
        resume_rows.append(
            {
                "task_id": task_id,
                "phase": "second_run",
                "expected_action": "SKIPPED",
                "observed_action": "SKIPPED" if skipped else "FAILED",
                "output_checksum_valid": skipped,
                "status": "PASS" if skipped else "FAIL",
            }
        )

    write_tsv(
        output_dir / "probe-telemetry.tsv",
        telemetry_rows,
        [
            "task_id",
            "platform",
            "run_accession",
            "rows",
            "max_rss_kib",
            "elapsed_seconds",
            "filesystem_inputs",
            "filesystem_outputs",
            "output_sha256",
            "status",
            "scope",
        ],
    )
    write_tsv(
        output_dir / "resume-audit.tsv",
        resume_rows,
        [
            "task_id",
            "phase",
            "expected_action",
            "observed_action",
            "output_checksum_valid",
            "status",
        ],
    )
    telemetry_valid = all(
        row["status"] == "PASS"
        and int(row["max_rss_kib"]) > 0
        and float(row["elapsed_seconds"]) >= 0
        for row in telemetry_rows
    )
    resume_valid = all(row["status"] == "PASS" for row in resume_rows)
    smoke_checks = [
        (
            "first_run_completed",
            3,
            first_completed,
            first_completed == 3,
        ),
        (
            "second_run_skipped",
            3,
            second_skipped,
            second_skipped == 3,
        ),
        (
            "observed_metric_rows",
            EXPECTED_METRICS_ROWS,
            observed_rows,
            observed_rows == EXPECTED_METRICS_ROWS,
        ),
        (
            "telemetry_rows_valid",
            True,
            telemetry_valid,
            telemetry_valid,
        ),
        (
            "resume_rows_valid",
            True,
            resume_valid,
            resume_valid,
        ),
    ]
    for name, expected, observed, passed in smoke_checks:
        add_check(
            checks,
            "local-smoke",
            name,
            expected,
            observed,
            "PASS" if passed else "FAIL",
            "native local workflow plumbing only",
        )

    figure_outputs: list[str] = []
    figure_outputs.extend(plot_input_footprint(footprint_rows, figure_dir))
    figure_outputs.extend(plot_control_loop(figure_dir))
    figure_outputs.extend(
        plot_resume_contract(first_completed, second_skipped, figure_dir)
    )
    for figure_name in figure_outputs:
        figure_path = figure_dir / figure_name
        add_check(
            checks,
            "figure",
            f"{figure_name}_exists",
            True,
            figure_path.is_file() and figure_path.stat().st_size > 0,
            "PASS"
            if figure_path.is_file() and figure_path.stat().st_size > 0
            else "FAIL",
            "publication artifact",
        )

    status_counts = Counter(str(row["status"]) for row in checks)
    final_status = "passed" if status_counts.get("FAIL", 0) == 0 else "failed"
    summary = {
        "status": final_status,
        "generated_on": date.today().isoformat(),
        "validation_scope": "native-linux-local-smoke",
        "performance_scope": "workflow-plumbing-smoke-not-assembler-benchmark",
        "source_project": "PRJEB52977",
        "source_sample": "SAMEA14435832",
        "source_fastq_rows": len(source_rows),
        "source_compressed_bytes": total_bytes,
        "source_compressed_gb_decimal": round(total_bytes / 1e9, 6),
        "array_task_count": len(task_rows),
        "first_run_completed": first_completed,
        "second_run_skipped": second_skipped,
        "metric_rows_validated": observed_rows,
        "smoke_seed": SMOKE_SEED,
        "host": {
            "operating_system": inventory["operating_system"],
            "kernel_release": inventory["kernel_release"],
            "architecture": inventory["architecture"],
            "logical_cpus": inventory["logical_cpus"],
            "physical_cores": inventory["physical_cores"],
            "memory_total_bytes": inventory["memory_total"],
            "memory_available_snapshot_bytes": inventory[
                "memory_available_snapshot"
            ],
            "project_filesystem_type": inventory["project_filesystem_type"],
            "project_filesystem_total_bytes": inventory[
                "project_filesystem_total"
            ],
            "project_filesystem_available_snapshot_bytes": inventory[
                "project_filesystem_available_snapshot"
            ],
        },
        "container": {
            "recommended_apptainer_release": EXPECTED_APPTAINER_RELEASE,
            "fixed_oci_reference": f"docker://alpine@{EXPECTED_OCI_DIGEST}",
            "local_apptainer_version": apptainer_version,
            "local_singularity_version": singularity_version,
            "legacy_sif_exec": "EXPECTED_FAIL",
            "legacy_sif_sha256": EXPECTED_SIF_SHA256,
        },
        "scheduler": {
            "slurm_control_plane_validated": bool(sbatch_path and sacct_path),
            "sbatch_available": bool(sbatch_path),
            "sacct_available": bool(sacct_path),
            "template_bash_syntax_validated": shell_result.returncode == 0,
            "template_local_body_validated": first_completed == 3
            and second_skipped == 3,
            "shellcheck_status": shellcheck_status,
        },
        "check_counts": dict(sorted(status_counts.items())),
        "input_sha256": {
            name: sha256(small / name) for name in EXPECTED_SHA256
        },
        "figures": sorted(figure_outputs),
    }
    (output_dir / "compute-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_tsv(
        output_dir / "compute-audit.tsv",
        checks,
        ["category", "check", "expected", "observed", "status", "scope"],
    )
    log_lines.extend(
        [
            "[summary]",
            f"status={final_status}",
            f"PASS={status_counts.get('PASS', 0)}",
            f"NOT_RUN={status_counts.get('NOT_RUN', 0)}",
            f"EXPECTED_FAIL={status_counts.get('EXPECTED_FAIL', 0)}",
            f"FAIL={status_counts.get('FAIL', 0)}",
            f"first_run_completed={first_completed}",
            f"second_run_skipped={second_skipped}",
            f"metric_rows_validated={observed_rows}",
        ]
    )
    (output_dir / "compute-validation.log").write_text(
        "\n".join(line for line in log_lines if line is not None) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if final_status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
