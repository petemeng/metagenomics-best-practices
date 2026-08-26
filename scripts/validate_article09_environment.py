#!/usr/bin/env python3
"""Validate the Article 09 Linux execution plane and create audit figures.

The current machine may be native Linux or the Linux guest inside WSL 2.
Windows-side WSL installation commands are intentionally never inferred from
native-Linux evidence.  A successful native-Linux run therefore reports the
Windows control-plane check as NOT_RUN while still validating the portable
Linux, conda/mamba, package, and real-data contracts.
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
import tempfile
from collections import Counter
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Sequence


EXPECTED_ENVIRONMENT = "metagenome-platform-smoke-2026.07"
EXPECTED_ENVIRONMENT_SHA256 = (
    "8f758b6ffdcf1561ece7d187ff34bc3f5a174fd8c6da66a101b206fcc869d20c"
)
EXPECTED_MINIFORGE_RELEASE = "26.3.2-2"
EXPECTED_MINIFORGE_SHA256 = (
    "42260ffe3830fb953d5eee1bbb32229ff06aa7c3833c1ed7a9a0420a95685d94"
)
EXPECTED_PACKAGES = {
    "python": "3.12.13",
    "seqkit": "2.10.0",
    "pigz": "2.8",
    "matplotlib-base": "3.10.5",
}
EXPECTED_INPUT_SHA256 = {
    "08-ena-fastq-sources.tsv": (
        "96638ae0d16ce5953760b596887c9f8f443c1dc33e0f82bef2e6ac040e69d962"
    ),
    "08-prefix-source-summary.json": (
        "f7bc59e303f26f723aacb2ff8a1d065411c1dcface225ce6d16969baefa66d2d"
    ),
    "08-read-prefix-metrics.tsv": (
        "6c426e0c89a95b1c1cc7ec631629e3635a7fd83c61c90d4086801b28495b5a61"
    ),
}
EXPECTED_PROJECT = "PRJEB52977"
EXPECTED_SAMPLE = "SAMEA14435832"
EXPECTED_SOURCE_ROWS = 4
EXPECTED_METRICS_ROWS = 15_000
EXPECTED_PLATFORMS = ("Illumina", "ONT", "PacBio")
EXPECTED_ROWS_PER_PLATFORM = 5_000


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
    parser.add_argument("--environment", default=EXPECTED_ENVIRONMENT)
    parser.add_argument("--data-smoke-only", action="store_true")
    parser.add_argument("--data-smoke-output", type=Path)
    return parser.parse_args()


def run_command(command: Sequence[str], timeout: int = 180) -> CommandResult:
    completed = subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
        env={**os.environ, "LC_ALL": "C", "LANG": "C"},
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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_version(text: str) -> str:
    match = re.search(r"\d+(?:\.\d+)+(?:[-+._A-Za-z0-9]*)?", text)
    return match.group(0) if match else text.strip()


def detect_host_scope() -> tuple[str, str]:
    proc_text = ""
    proc_version = Path("/proc/version")
    if proc_version.exists():
        proc_text = proc_version.read_text(encoding="utf-8", errors="replace")
    is_wsl = bool(
        re.search(r"microsoft|wsl", proc_text, flags=re.IGNORECASE)
        or os.environ.get("WSL_INTEROP")
        or os.environ.get("WSL_DISTRO_NAME")
    )
    if is_wsl:
        return "wsl-linux-execution-plane", "WSL detected from the Linux guest"
    return "native-linux", "Native Linux; Windows control plane is outside scope"


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


def redact(text: str, project_root: Path) -> str:
    """Remove local paths, proxy values, credentials, and token-like URL parts."""
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
    redacted = re.sub(
        r"(https?://)([^/@\s]+)@",
        r"\1<CREDENTIALS_REDACTED>@",
        redacted,
    )
    redacted = re.sub(
        r"(https?://[^/\s]+)/(?:t|token)/[^/\s]+/",
        r"\1/<TOKEN_REDACTED>/",
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


def data_smoke(project_root: Path) -> dict[str, object]:
    small = project_root / "data" / "small"
    source_path = small / "08-ena-fastq-sources.tsv"
    summary_path = small / "08-prefix-source-summary.json"
    metrics_path = small / "08-read-prefix-metrics.tsv"

    for path in (source_path, summary_path, metrics_path):
        if not path.exists():
            raise FileNotFoundError(path)

    hashes = {
        path.name: sha256(path)
        for path in (source_path, summary_path, metrics_path)
    }
    sources = read_tsv(source_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))

    platform_counts: Counter[str] = Counter()
    sample_counts: Counter[str] = Counter()
    metrics_rows = 0
    with metrics_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metrics_rows += 1
            platform_counts[row["PlatformKey"]] += 1
            sample_counts[row["SampleAccession"]] += 1

    checks: list[dict[str, object]] = []

    def record(check: str, expected: object, observed: object) -> None:
        checks.append(
            {
                "metric": check,
                "expected": expected,
                "observed": observed,
                "status": "PASS" if observed == expected else "FAIL",
            }
        )

    for filename, expected_hash in EXPECTED_INPUT_SHA256.items():
        record(f"sha256:{filename}", expected_hash, hashes[filename])
    record("project_accession", EXPECTED_PROJECT, summary.get("source_project"))
    record("sample_accession", EXPECTED_SAMPLE, summary.get("sample_accession"))
    record("source_rows", EXPECTED_SOURCE_ROWS, len(sources))
    record(
        "source_sample_identity",
        [EXPECTED_SAMPLE],
        sorted({row["SampleAccession"] for row in sources}),
    )
    record("metrics_rows", EXPECTED_METRICS_ROWS, metrics_rows)
    record(
        "platforms",
        list(EXPECTED_PLATFORMS),
        sorted(platform_counts),
    )
    for platform_name in EXPECTED_PLATFORMS:
        record(
            f"rows:{platform_name}",
            EXPECTED_ROWS_PER_PLATFORM,
            platform_counts.get(platform_name, 0),
        )
    record(
        "metrics_sample_identity",
        {EXPECTED_SAMPLE: EXPECTED_METRICS_ROWS},
        dict(sample_counts),
    )
    record("raw_fastq_stored", False, summary.get("raw_fastq_stored"))
    record(
        "illumina_prefix_mates_synchronized",
        True,
        summary.get("illumina_prefix_mates_synchronized"),
    )
    included = [
        row["PlatformKey"]
        for row in sources
        if row["IncludeInMetrics"].upper() == "TRUE"
    ]
    record("metric_source_platforms", list(EXPECTED_PLATFORMS), included)

    failures = [row for row in checks if row["status"] == "FAIL"]
    return {
        "status": "passed" if not failures else "failed",
        "python_version": platform.python_version(),
        "project_accession": summary.get("source_project"),
        "sample_accession": summary.get("sample_accession"),
        "source_rows": len(sources),
        "metrics_rows": metrics_rows,
        "platform_counts": dict(platform_counts),
        "raw_fastq_stored": summary.get("raw_fastq_stored"),
        "checks": checks,
    }


def save_figure(fig: object, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        figure_dir / f"{stem}.pdf",
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        figure_dir / f"{stem}.png",
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        figure_dir / f"{stem}.tiff",
        dpi=350,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def draw_layer_map(
    figure_dir: Path,
    host_scope: str,
    environment_passed: bool,
    data_passed: bool,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titleweight": "bold",
        }
    )
    colors = {
        "navy": "#173F5F",
        "blue": "#20639B",
        "teal": "#2A9D8F",
        "amber": "#E9C46A",
        "ink": "#24323F",
        "pale": "#F4F7FA",
        "white": "#FFFFFF",
    }
    layers = [
        (
            "Windows control plane",
            "Install / update / verify WSL 2",
            "VERIFY ON WINDOWS",
            colors["amber"],
        ),
        (
            "Linux execution plane",
            "Bash, Linux filesystem, x86_64",
            "WSL GUEST" if host_scope.startswith("wsl") else "NATIVE-LINUX EVIDENCE",
            colors["blue"],
        ),
        (
            "Isolated tool environment",
            "Miniforge policy + pinned conda environment",
            "PASS" if environment_passed else "FAIL",
            colors["teal"] if environment_passed else "#C44536",
        ),
        (
            "Real-data smoke",
            "PRJEB52977 / SAMEA14435832 frozen evidence",
            "PASS" if data_passed else "FAIL",
            colors["teal"] if data_passed else "#C44536",
        ),
    ]

    fig, ax = plt.subplots(figsize=(10.2, 7.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title(
        "A reproducible Windows workflow has two evidence planes",
        loc="left",
        fontsize=16,
        color=colors["navy"],
        pad=18,
    )
    ax.text(
        0,
        9.45,
        "Official Windows instructions and locally executed Linux evidence are reported separately.",
        color=colors["ink"],
        fontsize=10.5,
    )

    y_positions = [7.55, 5.6, 3.65, 1.7]
    for index, ((title, subtitle, badge, color), y) in enumerate(
        zip(layers, y_positions)
    ):
        box = FancyBboxPatch(
            (0.35, y - 0.65),
            9.25,
            1.28,
            boxstyle="round,pad=0.025,rounding_size=0.12",
            linewidth=1.25,
            edgecolor=color,
            facecolor=colors["pale"],
        )
        ax.add_patch(box)
        ax.add_patch(
            FancyBboxPatch(
                (0.58, y - 0.35),
                2.25,
                0.7,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                linewidth=0,
                facecolor=color,
            )
        )
        ax.text(
            1.705,
            y,
            badge,
            ha="center",
            va="center",
            color=colors["white"] if index != 0 else colors["ink"],
            fontsize=8.5,
            weight="bold",
        )
        ax.text(
            3.15,
            y + 0.18,
            title,
            ha="left",
            va="center",
            color=colors["navy"],
            fontsize=12,
            weight="bold",
        )
        ax.text(
            3.15,
            y - 0.22,
            subtitle,
            ha="left",
            va="center",
            color=colors["ink"],
            fontsize=9.5,
        )
        if index < len(layers) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (5.0, y - 0.7),
                    (5.0, y_positions[index + 1] + 0.7),
                    arrowstyle="-|>",
                    mutation_scale=14,
                    linewidth=1.1,
                    color=colors["navy"],
                )
            )

    ax.text(
        0.35,
        0.55,
        "Publication rule: do not describe native-Linux evidence as a Windows or WSL execution.",
        color=colors["ink"],
        fontsize=9.5,
        style="italic",
    )
    save_figure(fig, figure_dir, "09-wsl2-layer-map")
    plt.close(fig)


def draw_validation_figure(
    figure_dir: Path,
    environment_rows: list[dict[str, object]],
    data_rows: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titleweight": "bold",
        }
    )
    palette = {
        "PASS": "#2A9D8F",
        "FAIL": "#C44536",
        "NOT_RUN": "#E9C46A",
        "ink": "#24323F",
        "navy": "#173F5F",
        "grid": "#DDE5EC",
    }
    categories = [
        "Windows control",
        "Linux platform",
        "Manager policy",
        "Pinned packages",
        "Real-data smoke",
    ]
    row_groups: dict[str, list[str]] = {category: [] for category in categories}
    for row in environment_rows:
        category = str(row["category"])
        target = {
            "Windows control": "Windows control",
            "Linux platform": "Linux platform",
            "Manager policy": "Manager policy",
            "Pinned environment": "Pinned packages",
        }.get(category)
        if target:
            row_groups[target].append(str(row["status"]))
    row_groups["Real-data smoke"].extend(
        str(row["status"]) for row in data_rows
    )

    counts = {
        category: Counter(row_groups[category])
        for category in categories
    }
    fig, ax = plt.subplots(figsize=(10.2, 6.2))
    left = [0] * len(categories)
    for status_name in ("PASS", "NOT_RUN", "FAIL"):
        values = [counts[category].get(status_name, 0) for category in categories]
        ax.barh(
            categories,
            values,
            left=left,
            color=palette[status_name],
            height=0.62,
            label=status_name.replace("_", " "),
        )
        left = [a + b for a, b in zip(left, values)]

    for index, category in enumerate(categories):
        total = sum(counts[category].values())
        ax.text(
            total + 0.25,
            index,
            f"{total} check" if total == 1 else f"{total} checks",
            va="center",
            color=palette["ink"],
            fontsize=9,
        )

    fig.suptitle(
        "Environment acceptance audit",
        x=0.105,
        y=0.98,
        ha="left",
        fontsize=16,
        color=palette["navy"],
        weight="bold",
    )
    ax.set_xlabel("Number of checks")
    ax.set_ylabel("")
    ax.grid(axis="x", color=palette["grid"], linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        ncol=3,
        loc="upper right",
        bbox_to_anchor=(0.95, 0.982),
    )
    ax.invert_yaxis()
    ax.text(
        0,
        -0.18,
        "NOT RUN is an explicit evidence boundary, not a failed Linux validation.",
        transform=ax.transAxes,
        color=palette["ink"],
        fontsize=9.3,
        style="italic",
    )
    fig.subplots_adjust(top=0.82, bottom=0.18, left=0.24, right=0.93)
    save_figure(fig, figure_dir, "09-environment-validation")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()

    if args.data_smoke_only:
        if args.data_smoke_output is None:
            raise SystemExit("--data-smoke-output is required with --data-smoke-only")
        result = data_smoke(root)
        args.data_smoke_output.parent.mkdir(parents=True, exist_ok=True)
        args.data_smoke_output.write_text(
            json.dumps(result, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        return 0 if result["status"] == "passed" else 1

    if args.output_dir is None or args.figure_dir is None:
        raise SystemExit("--output-dir and --figure-dir are required")

    output_dir = args.output_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    matplotlib_cache = Path(tempfile.mkdtemp(prefix="article09-matplotlib-"))
    os.environ["MPLCONFIGDIR"] = str(matplotlib_cache)

    environment_rows: list[dict[str, object]] = []
    host_scope, host_scope_note = detect_host_scope()
    is_wsl = host_scope.startswith("wsl")
    today = date.today().isoformat()

    add_check(
        environment_rows,
        "Windows control",
        "WSL 2 installation and VERSION 2",
        "Capture wsl.exe evidence on Windows",
        "Linux guest detected; Windows capture still required"
        if is_wsl
        else "Not a Windows/WSL host",
        "NOT_RUN",
        "Windows host only",
    )
    add_check(
        environment_rows,
        "Linux platform",
        "validation scope",
        "Explicit host classification",
        host_scope,
        "PASS",
        host_scope_note,
    )
    add_check(
        environment_rows,
        "Linux platform",
        "kernel",
        "Linux",
        platform.system(),
        "PASS" if platform.system() == "Linux" else "FAIL",
        "Executed locally",
    )
    machine = platform.machine()
    add_check(
        environment_rows,
        "Linux platform",
        "architecture",
        "x86_64",
        machine,
        "PASS" if machine == "x86_64" else "FAIL",
        "Pinned Linux x86_64 environment",
    )

    findmnt = run_command(
        ("findmnt", "-T", str(root), "-no", "FSTYPE,TARGET")
    )
    filesystem = (
        findmnt.stdout.split(maxsplit=1)[0]
        if findmnt.returncode == 0 and findmnt.stdout
        else "unknown"
    )
    under_windows_mount = bool(
        re.match(r"^/mnt/[A-Za-z](?:/|$)", str(root))
    )
    storage_ok = findmnt.returncode == 0 and (
        not is_wsl or not under_windows_mount
    )
    add_check(
        environment_rows,
        "Linux platform",
        "project storage",
        "Linux filesystem; WSL projects not under /mnt/<drive>",
        f"{filesystem}; "
        + (
            "WSL Windows-mounted path"
            if under_windows_mount
            else "not a WSL Windows-mounted path"
        ),
        "PASS" if storage_ok else "FAIL",
        "Executed locally; mount target path redacted",
    )

    conda_path = shutil.which("conda")
    mamba_path = shutil.which("mamba")
    conda_version = (
        run_command((conda_path, "--version")) if conda_path else None
    )
    mamba_version = (
        run_command((mamba_path, "--version")) if mamba_path else None
    )
    add_check(
        environment_rows,
        "Manager policy",
        "conda entrypoint",
        "Available",
        normalize_version(conda_version.stdout)
        if conda_version and conda_version.returncode == 0
        else "Unavailable",
        "PASS"
        if conda_version and conda_version.returncode == 0
        else "FAIL",
        "Host environment manager",
    )
    add_check(
        environment_rows,
        "Manager policy",
        "mamba entrypoint",
        "Available",
        normalize_version(mamba_version.stdout)
        if mamba_version and mamba_version.returncode == 0
        else "Unavailable",
        "PASS"
        if mamba_version and mamba_version.returncode == 0
        else "FAIL",
        "Host environment manager",
    )

    config: dict[str, object] = {}
    if conda_path:
        config_result = run_command(
            (
                conda_path,
                "config",
                "--show",
                "channels",
                "channel_priority",
                "solver",
                "--json",
            )
        )
        if config_result.returncode == 0:
            config = json.loads(config_result.stdout)
    channels = config.get("channels", [])
    priority = config.get("channel_priority")
    solver = config.get("solver")
    add_check(
        environment_rows,
        "Manager policy",
        "host channels",
        ["conda-forge"],
        channels,
        "PASS" if channels == ["conda-forge"] else "FAIL",
        "Only safe config fields queried; full .condarc omitted",
    )
    add_check(
        environment_rows,
        "Manager policy",
        "channel priority",
        "strict",
        priority,
        "PASS" if priority == "strict" else "FAIL",
        "Host manager policy",
    )
    add_check(
        environment_rows,
        "Manager policy",
        "solver",
        "libmamba",
        solver,
        "PASS" if solver == "libmamba" else "FAIL",
        "Host manager policy",
    )

    environment_path = root / "env" / "platform-smoke.yml"
    environment_hash = sha256(environment_path) if environment_path.exists() else ""
    add_check(
        environment_rows,
        "Pinned environment",
        "environment YAML SHA-256",
        EXPECTED_ENVIRONMENT_SHA256,
        environment_hash or "missing",
        "PASS"
        if environment_hash == EXPECTED_ENVIRONMENT_SHA256
        else "FAIL",
        "Repository contract",
    )
    notice_text = (
        (root / "data" / "small" / "09-data-NOTICE.txt").read_text(
            encoding="utf-8"
        )
        if (root / "data" / "small" / "09-data-NOTICE.txt").exists()
        else ""
    )
    installer_pin_ok = (
        EXPECTED_MINIFORGE_RELEASE in notice_text
        and EXPECTED_MINIFORGE_SHA256 in notice_text
    )
    add_check(
        environment_rows,
        "Pinned environment",
        "Miniforge release record",
        f"{EXPECTED_MINIFORGE_RELEASE}; {EXPECTED_MINIFORGE_SHA256}",
        "Release and checksum recorded"
        if installer_pin_ok
        else "Pin missing",
        "PASS" if installer_pin_ok else "FAIL",
        "Documentation pin only; installer not downloaded on this host",
    )

    package_rows: list[dict[str, object]] = []
    packages: dict[str, dict[str, object]] = {}
    explicit_text = ""
    env_available = False
    if conda_path:
        list_result = run_command(
            (conda_path, "list", "-n", args.environment, "--json")
        )
        if list_result.returncode == 0:
            env_available = True
            for item in json.loads(list_result.stdout):
                packages[str(item["name"])] = item
        explicit_result = run_command(
            (conda_path, "list", "--explicit", "-n", args.environment)
        )
        if explicit_result.returncode == 0:
            explicit_text = redact(explicit_result.stdout, root)

    add_check(
        environment_rows,
        "Pinned environment",
        "named environment",
        args.environment,
        args.environment if env_available else "Unavailable",
        "PASS" if env_available else "FAIL",
        "Executed locally",
    )

    for package_name, expected_version in EXPECTED_PACKAGES.items():
        item = packages.get(package_name, {})
        installed_version = str(item.get("version", "missing"))
        package_status = (
            "PASS" if installed_version == expected_version else "FAIL"
        )
        package_rows.append(
            {
                "package": package_name,
                "requested_version": expected_version,
                "installed_version": installed_version,
                "build": item.get("build_string", item.get("build", "")),
                "channel": item.get("channel", ""),
                "status": package_status,
            }
        )
        add_check(
            environment_rows,
            "Pinned environment",
            f"{package_name} package",
            expected_version,
            installed_version,
            package_status,
            "Exact version from conda package metadata",
        )

    cli_versions: dict[str, str] = {}
    cli_commands = {
        "python": ("python", "--version"),
        "seqkit": ("seqkit", "version"),
        "pigz": ("pigz", "--version"),
        "matplotlib-base": (
            "python",
            "-c",
            "import matplotlib; print(matplotlib.__version__)",
        ),
    }
    for name, command in cli_commands.items():
        if not conda_path or not env_available:
            cli_versions[name] = "missing"
            continue
        result = run_command(
            (conda_path, "run", "-n", args.environment, *command)
        )
        cli_versions[name] = (
            normalize_version(f"{result.stdout}\n{result.stderr}")
            if result.returncode == 0
            else "command failed"
        )
        add_check(
            environment_rows,
            "Pinned environment",
            f"{name} command",
            EXPECTED_PACKAGES[name],
            cli_versions[name],
            "PASS"
            if result.returncode == 0
            and cli_versions[name] == EXPECTED_PACKAGES[name]
            else "FAIL",
            "Executed through conda run",
        )

    explicit_path = output_dir / "explicit-linux-64.txt"
    if explicit_text:
        explicit_path.write_text(explicit_text + "\n", encoding="utf-8")
    add_check(
        environment_rows,
        "Pinned environment",
        "explicit package lock",
        "Generated from the validated environment",
        "Generated" if explicit_path.exists() else "Missing",
        "PASS" if explicit_path.exists() else "FAIL",
        "Linux x86_64 platform lock; URLs contain no credentials",
    )

    data_result: dict[str, object]
    data_command_status = 1
    with tempfile.TemporaryDirectory(prefix="article09-smoke-") as temp_dir:
        temp_json = Path(temp_dir) / "data-smoke.json"
        if conda_path and env_available:
            data_command = run_command(
                (
                    conda_path,
                    "run",
                    "-n",
                    args.environment,
                    "python",
                    str(root / "scripts" / "validate_article09_environment.py"),
                    "--project-root",
                    str(root),
                    "--data-smoke-only",
                    "--data-smoke-output",
                    str(temp_json),
                ),
                timeout=300,
            )
            data_command_status = data_command.returncode
        if temp_json.exists():
            data_result = json.loads(temp_json.read_text(encoding="utf-8"))
        else:
            data_result = {
                "status": "failed",
                "python_version": "",
                "checks": [
                    {
                        "metric": "conda_data_smoke",
                        "expected": "completed",
                        "observed": "not completed",
                        "status": "FAIL",
                    }
                ],
            }

    data_rows = list(data_result["checks"])
    data_rows.append(
        {
            "metric": "smoke_python_version",
            "expected": EXPECTED_PACKAGES["python"],
            "observed": data_result.get("python_version", ""),
            "status": "PASS"
            if data_result.get("python_version") == EXPECTED_PACKAGES["python"]
            and data_command_status == 0
            else "FAIL",
        }
    )

    write_tsv(
        output_dir / "environment-audit.tsv",
        environment_rows,
        ("category", "check", "expected", "observed", "status", "scope"),
    )
    write_tsv(
        output_dir / "package-audit.tsv",
        package_rows,
        (
            "package",
            "requested_version",
            "installed_version",
            "build",
            "channel",
            "status",
        ),
    )
    write_tsv(
        output_dir / "data-smoke.tsv",
        data_rows,
        ("metric", "expected", "observed", "status"),
    )

    environment_failures = [
        row for row in environment_rows if row["status"] == "FAIL"
    ]
    data_failures = [row for row in data_rows if row["status"] == "FAIL"]
    overall_status = (
        "passed"
        if not environment_failures and not data_failures
        else "failed"
    )
    status_counts = Counter(str(row["status"]) for row in environment_rows)
    data_status_counts = Counter(str(row["status"]) for row in data_rows)
    summary = {
        "status": overall_status,
        "audit_date": today,
        "validation_scope": host_scope,
        "windows_control_plane_validated": False,
        "linux_execution_plane_validated": not any(
            row["status"] == "FAIL"
            for row in environment_rows
            if row["category"] == "Linux platform"
        ),
        "environment_name": args.environment,
        "environment_yaml_sha256": environment_hash,
        "miniforge_release_recommended": EXPECTED_MINIFORGE_RELEASE,
        "miniforge_installer_downloaded_on_host": False,
        "host_os": os_release(),
        "kernel_release": platform.release(),
        "architecture": machine,
        "filesystem_type": filesystem,
        "conda_version": (
            normalize_version(conda_version.stdout)
            if conda_version and conda_version.returncode == 0
            else ""
        ),
        "mamba_version": (
            normalize_version(mamba_version.stdout)
            if mamba_version and mamba_version.returncode == 0
            else ""
        ),
        "channel_policy": {
            "channels": channels,
            "channel_priority": priority,
            "solver": solver,
        },
        "cli_versions": cli_versions,
        "package_count": len(packages),
        "environment_pass_checks": status_counts.get("PASS", 0),
        "environment_not_run_checks": status_counts.get("NOT_RUN", 0),
        "environment_fail_checks": status_counts.get("FAIL", 0),
        "data_pass_checks": data_status_counts.get("PASS", 0),
        "data_fail_checks": data_status_counts.get("FAIL", 0),
        "environment_checks": dict(status_counts),
        "data_checks": dict(data_status_counts),
        "project_accession": data_result.get("project_accession"),
        "sample_accession": data_result.get("sample_accession"),
        "metrics_rows": data_result.get("metrics_rows"),
        "platform_counts": data_result.get("platform_counts"),
        "raw_fastq_stored": data_result.get("raw_fastq_stored"),
        "redaction": {
            "project_paths": True,
            "home_paths": True,
            "proxy_values": True,
            "full_condarc_exported": False,
        },
    }
    (output_dir / "environment-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    log_lines = [
        "Article 09 environment validation",
        f"date: {today}",
        f"status: {overall_status}",
        f"validation_scope: {host_scope}",
        "windows_control_plane_validated: false",
        f"host_os: {os_release()}",
        f"kernel: Linux {platform.release()} {machine}",
        f"filesystem_type: {filesystem}",
        f"conda_version: {summary['conda_version']}",
        f"mamba_version: {summary['mamba_version']}",
        f"channels: {', '.join(str(value) for value in channels)}",
        f"channel_priority: {priority}",
        f"solver: {solver}",
        f"environment: {args.environment}",
        f"environment_yaml_sha256: {environment_hash}",
        "package_versions: "
        + ", ".join(
            f"{row['package']}={row['installed_version']}"
            for row in package_rows
        ),
        f"real_data_project: {data_result.get('project_accession')}",
        f"real_data_sample: {data_result.get('sample_accession')}",
        f"metrics_rows: {data_result.get('metrics_rows')}",
        f"raw_fastq_stored: {data_result.get('raw_fastq_stored')}",
        (
            "environment_checks: "
            f"PASS={status_counts.get('PASS', 0)}, "
            f"NOT_RUN={status_counts.get('NOT_RUN', 0)}, "
            f"FAIL={status_counts.get('FAIL', 0)}"
        ),
        (
            "data_checks: "
            f"PASS={data_status_counts.get('PASS', 0)}, "
            f"FAIL={data_status_counts.get('FAIL', 0)}"
        ),
        "sensitive_config: proxy values, user paths, and full .condarc omitted",
    ]
    (output_dir / "environment-validation.log").write_text(
        redact("\n".join(log_lines) + "\n", root),
        encoding="utf-8",
    )

    draw_layer_map(
        figure_dir,
        host_scope,
        not environment_failures,
        not data_failures,
    )
    draw_validation_figure(figure_dir, environment_rows, data_rows)
    shutil.rmtree(matplotlib_cache, ignore_errors=True)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if overall_status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
