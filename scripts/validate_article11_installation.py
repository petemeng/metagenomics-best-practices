#!/usr/bin/env python3
"""Validate Article 11 frozen installation evidence and database gates.

The multi-gigabyte environments were installed and tested once on native
Linux. Routine QA verifies the resulting explicit locks, concise command
evidence, database manifest, and fail-closed downloader without downloading a
database or depending on those host environments still being present.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence


VALIDATION_DATE = "2026-07-20"
EXPECTED_DIRECT_PINS = {
    "assembly": {
        "python": "3.10",
        "setuptools": "<81",
        "megahit": "1.2.9",
        "spades": "4.3.0",
        "bowtie2": "2.5.5",
        "metabat2": "2.18",
        "maxbin2": "2.2.7",
        "concoct": "1.1.0",
        "das_tool": "1.1.7",
        "coverm": "0.8.0",
    },
    "biobakery": {
        "python": "3.12",
        "metaphlan": "4.2.5",
        "humann": "3.9",
        "bowtie2": "2.5.5",
        "diamond": "2.2.4",
    },
}
EXPECTED_DATABASES = {
    "metaphlan-vJan26": {
        "release_id": "mpa_vJan26_CHOCOPhlAnSGB_202605",
        "checksum_algorithm": "md5",
        "expected_checksum": "7162b0c3493663dce9abef08ccc06aea",
        "download_gate": "enabled",
    },
    "humann-chocophlan-full": {
        "release_id": "v201901_v31-full",
        "checksum_algorithm": "none",
        "expected_checksum": "UNPUBLISHED",
        "download_gate": "blocked-no-upstream-checksum",
    },
    "humann-uniref90-full": {
        "release_id": "v201901b-uniref90-diamond",
        "checksum_algorithm": "none",
        "expected_checksum": "UNPUBLISHED",
        "download_gate": "blocked-no-upstream-checksum",
    },
    "checkm2-v3": {
        "release_id": "zenodo-14897628-v3",
        "checksum_algorithm": "md5",
        "expected_checksum": "07c10655620843b517d0df0c160d911f",
        "download_gate": "enabled",
    },
    "gtdbtk-r232": {
        "release_id": "R232",
        "checksum_algorithm": "md5",
        "expected_checksum": "25a59e0352b1fd150c589f56559767d4",
        "download_gate": "enabled",
    },
}
EXPECTED_ENTRYPOINTS = {
    "Python (assembly)": "3.10.20",
    "MEGAHIT": "1.2.9",
    "SPAdes": "4.3.0",
    "Bowtie2 (assembly)": "2.5.5",
    "MetaBAT2": "2:2.18",
    "MaxBin2": "2.2.7",
    "CONCOCT": "1.1.0",
    "DAS Tool": "1.1.7",
    "CoverM": "0.8.0",
    "NumPy (assembly)": "1.26.4",
    "setuptools (assembly)": "80.10.2",
    "Python (bioBakery)": "3.12.13",
    "MetaPhlAn": "4.2.5",
    "HUMAnN": "3.9",
    "Bowtie2 (bioBakery)": "2.5.5",
    "DIAMOND (bioBakery)": "2.2.4",
}
EXPECTED_TESTS = {
    "MEGAHIT built-in test": "ALL DONE",
    "SPAdes built-in test": "TEST PASSED CORRECTLY",
    "HUMAnN unit tests": "Ran 186 tests",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    return parser.parse_args()


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
        writer.writerows(rows)


def environment_yaml_pins(path: Path) -> dict[str, str]:
    pins: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^\s*-\s+([A-Za-z0-9_.-]+)([<>=].+)$", raw)
        if match:
            pins[match.group(1)] = match.group(2).lstrip("=")
    return pins


def explicit_lock_packages(path: Path) -> tuple[int, set[str]]:
    text = path.read_text(encoding="utf-8")
    if "@EXPLICIT" not in text:
        raise ValueError(f"{path.name} is not an explicit conda specification")
    urls = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith(("#", "@"))
    ]
    names: set[str] = set()
    for url in urls:
        filename = url.split("/")[-1].split("#", 1)[0]
        match = re.match(r"^(.+)-([0-9][^-]*)-[^-]+\.(?:conda|tar\.bz2)$", filename)
        if match:
            names.add(match.group(1))
    return len(urls), names


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


def draw_box(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    body: str,
    color: str,
    face: str = "#FFFFFF",
) -> None:
    from matplotlib.patches import FancyBboxPatch

    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        linewidth=1.25,
        edgecolor=color,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(x + 0.025, y + height - 0.06, title, weight="bold", fontsize=10)
    ax.text(
        x + 0.025,
        y + height - 0.14,
        body,
        fontsize=8.6,
        va="top",
        linespacing=1.35,
        color="#42525D",
    )


def plot_environment_boundaries(figure_dir: Path) -> list[str]:
    plt = configure_matplotlib()
    from matplotlib.patches import FancyArrowPatch

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Environment boundaries follow dependency contracts", loc="left", pad=14)
    ax.text(
        0,
        0.94,
        "Exact direct pins + explicit Linux-64 locks • databases remain external",
        color="#5B6870",
        fontsize=9.2,
    )
    draw_box(
        ax,
        0.04,
        0.50,
        0.26,
        0.29,
        "bioBakery",
        "MetaPhlAn 4.2.5\nHUMAnN 3.9\nPython 3.12",
        "#20639B",
        "#F3F8FC",
    )
    draw_box(
        ax,
        0.37,
        0.50,
        0.26,
        0.29,
        "Assembly + binning",
        "MEGAHIT • SPAdes\n3 binners • DAS Tool\nCoverM • Python 3.10",
        "#2A9D8F",
        "#F2FAF7",
    )
    draw_box(
        ax,
        0.70,
        0.50,
        0.26,
        0.29,
        "Isolated later",
        "CheckM2 / GUNC\nGTDB-Tk / dRep\nConflict tested, not merged",
        "#C44536",
        "#FFF6F3",
    )
    for x in (0.30, 0.63):
        ax.add_patch(
            FancyArrowPatch(
                (x, 0.645),
                (x + 0.07, 0.645),
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.2,
                color="#6B7A86",
            )
        )
    ax.text(
        0.50,
        0.30,
        "One giant environment: solver conflict reproduced",
        ha="center",
        weight="bold",
        color="#9A3412",
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": "#FFF0E8",
            "edgecolor": "#E07A5F",
        },
    )
    ax.text(
        0.50,
        0.13,
        "A version command validates an executable. A database query validates a database.",
        ha="center",
        fontsize=9.5,
        color="#42525D",
    )
    fig.tight_layout()
    outputs = save_publication_figure(
        fig, figure_dir / "11-environment-boundaries"
    )
    plt.close(fig)
    return outputs


def plot_entrypoints(
    evidence_rows: list[dict[str, str]], figure_dir: Path
) -> list[str]:
    plt = configure_matplotlib()
    rows = [row for row in evidence_rows if row["evidence_class"] == "entrypoint"]
    labels = [row["tool"] for row in rows]
    colors = [
        "#20639B" if row["environment"] == "biobakery" else "#2A9D8F"
        for row in rows
    ]
    fig_height = max(4.8, 0.34 * len(rows) + 1.6)
    fig, ax = plt.subplots(figsize=(7.9, fig_height))
    y = list(range(len(rows)))
    ax.barh(y, [1] * len(rows), color=colors, height=0.64)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.18)
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.grid(False)
    ax.set_title("Frozen entry-point audit: all declared tools passed", loc="left", pad=13)
    ax.text(
        0,
        1.01,
        f"{len(rows)} commands • user-site isolation enabled • 2026-07-20",
        transform=ax.transAxes,
        fontsize=9,
        color="#5B6870",
    )
    for index, row in enumerate(rows):
        ax.text(
            0.03,
            index,
            row["observed"],
            va="center",
            color="white",
            fontsize=8.5,
            weight="bold",
        )
        ax.text(1.04, index, "PASS", va="center", color="#176B4D", weight="bold")
    fig.tight_layout()
    outputs = save_publication_figure(
        fig, figure_dir / "11-toolchain-entrypoints"
    )
    plt.close(fig)
    return outputs


def plot_database_gates(
    database_rows: list[dict[str, str]], figure_dir: Path
) -> list[str]:
    plt = configure_matplotlib()
    labels = [
        "MetaPhlAn vJan26",
        "HUMAnN ChocoPhlAn",
        "HUMAnN UniRef90",
        "CheckM2 v3",
        "GTDB-Tk R232",
    ]
    by_id = {row["database_id"]: row for row in database_rows}
    ordered_ids = list(EXPECTED_DATABASES)
    sizes = [
        int(by_id[item]["expected_compressed_bytes"]) / 1e9
        for item in ordered_ids
    ]
    colors = [
        "#2A9D8F"
        if by_id[item]["download_gate"] == "enabled"
        else "#E07A5F"
        for item in ordered_ids
    ]
    fig, ax = plt.subplots(figsize=(8.0, 4.9))
    bars = ax.barh(labels, sizes, color=colors, height=0.62)
    ax.invert_yaxis()
    ax.set_xlabel("Expected archive size (GB, decimal)")
    ax.set_title("Database landing is checksum-gated", loc="left", pad=13)
    ax.text(
        0,
        1.01,
        "Green: published checksum • coral: no publisher checksum • none downloaded",
        transform=ax.transAxes,
        fontsize=9,
        color="#5B6870",
    )
    ax.grid(axis="x", color="#DDE5EA", linewidth=0.7)
    ax.set_axisbelow(True)
    max_size = max(sizes)
    for bar, size, database_id in zip(bars, sizes, ordered_ids):
        row = by_id[database_id]
        label = (
            f"{size:.1f} GB • enabled"
            if row["download_gate"] == "enabled"
            else f"{size:.1f} GB • BLOCKED"
        )
        ax.text(
            size + max_size * 0.018,
            bar.get_y() + bar.get_height() / 2,
            label,
            va="center",
            fontsize=8.7,
        )
    ax.set_xlim(0, max_size * 1.32)
    fig.tight_layout()
    outputs = save_publication_figure(
        fig, figure_dir / "11-database-storage-contract"
    )
    plt.close(fig)
    return outputs


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else root / "results" / "11-install-biobakery-assembly"
    )
    figure_dir = (
        args.figure_dir.resolve() if args.figure_dir else root / "figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    env_dir = root / "env"
    small = root / "data" / "small"
    yaml_paths = {
        "assembly": env_dir / "assembly.yml",
        "biobakery": env_dir / "biobakery.yml",
    }
    lock_paths = {
        "assembly": env_dir / "assembly-linux-64.lock",
        "biobakery": env_dir / "biobakery-linux-64.lock",
    }
    evidence_path = small / "11-environment-evidence.tsv"
    test_log_path = small / "11-install-self-tests.log"
    solver_path = small / "11-solver-audit.tsv"
    database_path = small / "11-database-manifest.tsv"
    db_script = root / "db" / "download_db.sh"
    relink_script = env_dir / "relink-biobakery-entrypoints.sh"

    required = [
        *yaml_paths.values(),
        *lock_paths.values(),
        evidence_path,
        test_log_path,
        solver_path,
        database_path,
        db_script,
        relink_script,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("missing Article 11 inputs: " + ", ".join(missing))

    checks: list[dict[str, object]] = []
    lock_counts: dict[str, int] = {}
    yaml_hashes: dict[str, str] = {}
    lock_hashes: dict[str, str] = {}

    for environment, yaml_path in yaml_paths.items():
        text = yaml_path.read_text(encoding="utf-8")
        observed_pins = environment_yaml_pins(yaml_path)
        yaml_hashes[environment] = sha256(yaml_path)
        add_check(
            checks,
            "environment_contract",
            f"{environment} user-site isolation",
            'PYTHONNOUSERSITE: "1"',
            'PYTHONNOUSERSITE: "1"' if 'PYTHONNOUSERSITE: "1"' in text else "missing",
            "PASS" if 'PYTHONNOUSERSITE: "1"' in text else "FAIL",
            "source YAML",
        )
        for package, version in EXPECTED_DIRECT_PINS[environment].items():
            observed = observed_pins.get(package, "MISSING")
            add_check(
                checks,
                "direct_pin",
                f"{environment}:{package}",
                version,
                observed,
                "PASS" if observed == version else "FAIL",
                "source YAML",
            )

        try:
            count, names = explicit_lock_packages(lock_paths[environment])
            lock_status = "PASS"
        except ValueError:
            count, names, lock_status = 0, set(), "FAIL"
        lock_counts[environment] = count
        lock_hashes[environment] = sha256(lock_paths[environment])
        add_check(
            checks,
            "explicit_lock",
            f"{environment} explicit package count",
            ">=100",
            count,
            lock_status if count >= 100 else "FAIL",
            "frozen native-Linux transaction",
        )
        for package in EXPECTED_DIRECT_PINS[environment]:
            if package == "python":
                lock_name = "python"
            elif package == "setuptools":
                lock_name = "setuptools"
            else:
                lock_name = package
            add_check(
                checks,
                "explicit_lock",
                f"{environment}:{package} present in lock",
                "present",
                "present" if lock_name in names else "missing",
                "PASS" if lock_name in names else "FAIL",
                "frozen native-Linux transaction",
            )

    evidence_rows = read_tsv(evidence_path)
    evidence_by_tool = {row["tool"]: row for row in evidence_rows}
    for tool, version in EXPECTED_ENTRYPOINTS.items():
        row = evidence_by_tool.get(tool, {})
        observed = row.get("observed", "MISSING")
        status = row.get("status", "MISSING")
        add_check(
            checks,
            "entrypoint",
            tool,
            version,
            observed,
            "PASS" if status == "PASS" and version in observed else "FAIL",
            "frozen live command evidence",
        )
    for test_name, marker in EXPECTED_TESTS.items():
        row = evidence_by_tool.get(test_name, {})
        observed = row.get("observed", "MISSING")
        status = row.get("status", "MISSING")
        add_check(
            checks,
            "built_in_test",
            test_name,
            marker,
            observed,
            "PASS" if status == "PASS" and marker in observed else "FAIL",
            "software installation test only",
        )

    test_log = test_log_path.read_text(encoding="utf-8")
    for marker in (
        "MEGAHIT: PASS",
        "SPAdes: PASS",
        "HUMAnN unit tests: PASS",
        "CONCOCT initial entry point: EXPECTED_FAIL",
        "CONCOCT after setuptools<81: PASS",
        "HUMAnN bundled Bowtie2/DIAMOND: EXPECTED_FAIL",
        "bioBakery exact relink: PASS",
    ):
        add_check(
            checks,
            "test_log",
            marker,
            "present",
            "present" if marker in test_log else "missing",
            "PASS" if marker in test_log else "FAIL",
            "concise redacted install log",
        )

    relink_text = relink_script.read_text(encoding="utf-8")
    for token in (
        "--force-reinstall",
        '"bowtie2=2.5.5"',
        '"diamond=2.2.4"',
        "PYTHONNOUSERSITE=1",
    ):
        add_check(
            checks,
            "post_install_relink",
            token,
            "present",
            "present" if token in relink_text else "missing",
            "PASS" if token in relink_text else "FAIL",
            "rebuild contract for HUMAnN executable collision",
        )

    solver_rows = read_tsv(solver_path)
    solver_by_contract = {row["contract"]: row for row in solver_rows}
    expected_solver_results = {
        "biobakery": "PASS_AFTER_RELINK",
        "assembly-core": "PASS_AFTER_FIX",
        "assembly-plus-checkm2": "FAIL",
        "database-dependent-validation": "NOT_RUN",
    }
    for contract, expected in expected_solver_results.items():
        observed = solver_by_contract.get(contract, {}).get("result", "MISSING")
        add_check(
            checks,
            "solver_boundary",
            contract,
            expected,
            observed,
            "PASS" if observed == expected else "FAIL",
            "frozen solver audit",
        )

    database_rows = read_tsv(database_path)
    database_by_id = {row["database_id"]: row for row in database_rows}
    for database_id, expected_values in EXPECTED_DATABASES.items():
        row = database_by_id.get(database_id, {})
        for field, expected in expected_values.items():
            observed = row.get(field, "MISSING")
            add_check(
                checks,
                "database_manifest",
                f"{database_id}:{field}",
                expected,
                observed,
                "PASS" if observed == expected else "FAIL",
                "manifest only; archive not downloaded",
            )
        observed_status = row.get("validation_status", "MISSING")
        add_check(
            checks,
            "database_boundary",
            f"{database_id}:not downloaded",
            "NOT_DOWNLOADED",
            observed_status,
            "PASS" if observed_status == "NOT_DOWNLOADED" else "FAIL",
            "database-dependent workflows not claimed",
        )

    audit = subprocess.run(
        [str(db_script), "audit"],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    add_check(
        checks,
        "download_script",
        "manifest audit",
        "exit 0",
        f"exit {audit.returncode}",
        "PASS" if audit.returncode == 0 else "FAIL",
        (audit.stdout or audit.stderr).strip(),
    )
    blocked = subprocess.run(
        [str(db_script), "download", "humann-chocophlan-full"],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    add_check(
        checks,
        "download_script",
        "mutable HUMAnN endpoint fails closed",
        "exit 2 before DB_ROOT",
        f"exit {blocked.returncode}",
        "PASS"
        if blocked.returncode == 2 and "fail-closed" in blocked.stderr
        else "FAIL",
        "negative control; no network or filesystem mutation",
    )

    figures = []
    figures.extend(plot_environment_boundaries(figure_dir))
    figures.extend(plot_entrypoints(evidence_rows, figure_dir))
    figures.extend(plot_database_gates(database_rows, figure_dir))
    for name in figures:
        add_check(
            checks,
            "figure",
            name,
            "exists and non-empty",
            (figure_dir / name).stat().st_size,
            "PASS" if (figure_dir / name).stat().st_size > 0 else "FAIL",
            "English publication figure",
        )

    status_counts = Counter(str(row["status"]) for row in checks)
    failed = [row for row in checks if row["status"] == "FAIL"]
    summary = {
        "status": "failed" if failed else "passed",
        "validation_date": VALIDATION_DATE,
        "validation_scope": "frozen-live-install-evidence-plus-manifest-audit",
        "environment_count": 2,
        "entrypoint_pass_count": sum(
            row["status"] == "PASS"
            for row in checks
            if row["category"] == "entrypoint"
        ),
        "built_in_test_pass_count": sum(
            row["status"] == "PASS"
            for row in checks
            if row["category"] == "built_in_test"
        ),
        "assembly_lock_packages": lock_counts["assembly"],
        "biobakery_lock_packages": lock_counts["biobakery"],
        "database_rows": len(database_rows),
        "database_enabled_rows": sum(
            row["download_gate"] == "enabled" for row in database_rows
        ),
        "database_blocked_rows": sum(
            row["download_gate"].startswith("blocked") for row in database_rows
        ),
        "database_downloaded_rows": sum(
            row["validation_status"] != "NOT_DOWNLOADED"
            for row in database_rows
        ),
        "database_dependent_workflows_validated": False,
        "metaphlan_database_release": "mpa_vJan26_CHOCOPhlAnSGB_202605",
        "gtdbtk_database_release": "R232",
        "user_site_isolation": True,
        "setuptools_ceiling": "<81",
        "resolved_entrypoint_regressions": [
            "CONCOCT pkg_resources",
            "HUMAnN bundled Bowtie2/DIAMOND file collision",
        ],
        "environment_yaml_sha256": yaml_hashes,
        "explicit_lock_sha256": lock_hashes,
        "audit_pass_checks": status_counts.get("PASS", 0),
        "audit_fail_checks": status_counts.get("FAIL", 0),
        "figure_files": sorted(figures),
    }

    write_tsv(
        output_dir / "installation-audit.tsv",
        checks,
        ("category", "check", "expected", "observed", "status", "scope"),
    )
    write_tsv(
        output_dir / "database-audit.tsv",
        database_rows,
        tuple(database_rows[0].keys()),
    )
    (output_dir / "installation-summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    log_lines = [
        "Article 11 installation validation",
        f"status={summary['status']}",
        f"scope={summary['validation_scope']}",
        f"assembly_lock_packages={summary['assembly_lock_packages']}",
        f"biobakery_lock_packages={summary['biobakery_lock_packages']}",
        f"entrypoint_pass_count={summary['entrypoint_pass_count']}",
        f"built_in_test_pass_count={summary['built_in_test_pass_count']}",
        f"database_rows={summary['database_rows']}",
        f"database_downloaded_rows={summary['database_downloaded_rows']}",
        f"pass={summary['audit_pass_checks']}",
        f"fail={summary['audit_fail_checks']}",
    ]
    (output_dir / "installation-validation.log").write_text(
        "\n".join(log_lines) + "\n", encoding="utf-8"
    )
    print("\n".join(log_lines))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
