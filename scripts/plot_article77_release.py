#!/usr/bin/env python3
"""Create deterministic English-only publication figures for Article 77."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


SEED = 20_260_777
COLORS = {
    "blue": "#0072B2",
    "sky": "#56B4E9",
    "green": "#009E73",
    "orange": "#E69F00",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "gray": "#8A8A8A",
    "light": "#EEF3F7",
    "dark": "#263238",
}
STEMS = (
    "77-repository-routing",
    "77-accession-graph",
    "77-release-layers",
    "77-artifact-readiness",
    "77-identifier-state",
    "77-container-digest",
    "77-database-manifest",
    "77-release-gates",
    "77-availability-statement",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def setup() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "svg.fonttype": "none",
        }
    )


def save(fig: plt.Figure, stem: str, output: Path) -> None:
    fig.savefig(output / f"{stem}.png", dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(output / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def read(input_dir: Path, filename: str) -> pd.DataFrame:
    return pd.read_csv(input_dir / filename, sep="\t")


def repository_routing(routes: pd.DataFrame, output: Path) -> None:
    assets = [
        "Raw reads", "Study description", "Source material", "Primary metagenome assembly",
        "MAG", "UViG / vOTU", "Gene catalog", "Functional annotation",
        "Figure source tables", "Workflow source", "Container image", "Database manifest",
    ]
    destinations = ["SRA / ENA", "INSDC assembly", "DOI repository", "GitHub", "OCI registry"]
    matrix = np.zeros((len(assets), len(destinations)), dtype=int)
    assignments = {
        "Raw reads": [0], "Study description": [0], "Source material": [0],
        "Primary metagenome assembly": [1], "MAG": [1], "UViG / vOTU": [1, 2],
        "Gene catalog": [2], "Functional annotation": [2], "Figure source tables": [2],
        "Workflow source": [2, 3], "Container image": [4], "Database manifest": [2, 3],
    }
    for i, asset in enumerate(assets):
        for j in assignments[asset]:
            matrix[i, j] = 2 if len(assignments[asset]) == 1 else 1
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    ax.imshow(matrix, cmap=ListedColormap(["#F2F2F2", "#F0E442", "#009E73"]), vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(range(len(destinations)), destinations)
    ax.set_yticks(range(len(assets)), assets)
    ax.tick_params(axis="x", rotation=22)
    for i in range(len(assets)):
        for j in range(len(destinations)):
            label = "Primary" if matrix[i, j] == 2 else ("Linked" if matrix[i, j] == 1 else "")
            ax.text(j, i, label, ha="center", va="center", fontsize=7.5, color=COLORS["dark"])
    ax.set_title("Route each release object to the repository that preserves its semantics")
    ax.set_xlabel("Preservation destination")
    ax.set_ylabel("Release object")
    ax.set_xticks(np.arange(-0.5, len(destinations), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(assets), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    save(fig, "77-repository-routing", output)


def box(ax: plt.Axes, xy: tuple[float, float], text: str, color: str, width: float = 1.7, height: float = 0.65) -> None:
    x, y = xy
    patch = FancyBboxPatch((x - width / 2, y - height / 2), width, height,
                           boxstyle="round,pad=0.04,rounding_size=0.08",
                           linewidth=1.4, edgecolor=color, facecolor=color + "22")
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=9, color=COLORS["dark"], weight="bold")


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = "#66757F") -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13,
                                linewidth=1.4, color=color, shrinkA=5, shrinkB=5))


def accession_graph(output: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.8))
    positions = {
        "Study\nPRJEB52977": (0.8, 2.6),
        "Sample\nSAMEA14435832": (3.0, 3.35),
        "Sample\nSAMEA14435833": (3.0, 1.85),
        "Experiment\nresolve from ENA": (5.2, 3.35),
        "Experiment\nresolve from ENA ": (5.2, 1.85),
        "Run\nERR9765746": (7.4, 3.35),
        "Run\nERR9765747": (7.4, 1.85),
        "Analysis / derived data\nlocal or future accession": (9.8, 2.6),
    }
    for label, pos in positions.items():
        state_color = COLORS["green"] if any(token in label for token in ("PRJ", "SAMEA", "ERR")) else COLORS["orange"]
        box(ax, pos, label.rstrip(), state_color, width=1.85 if "Analysis" not in label else 2.25)
    for start, end in [
        ((1.65, 2.7), (2.1, 3.2)), ((1.65, 2.5), (2.1, 2.0)),
        ((3.9, 3.35), (4.3, 3.35)), ((3.9, 1.85), (4.3, 1.85)),
        ((6.1, 3.35), (6.5, 3.35)), ((6.1, 1.85), (6.5, 1.85)),
        ((8.3, 3.25), (8.65, 2.85)), ((8.3, 1.95), (8.65, 2.35)),
    ]:
        arrow(ax, start, end)
    ax.text(4.6, 0.65, "Never infer an accession\nfrom a naming pattern", ha="center", color=COLORS["red"], weight="bold")
    ax.text(9.2, 0.65, "Submission IDs are not\npublication IDs", ha="center", color=COLORS["red"], weight="bold")
    ax.set_xlim(-0.3, 11.1)
    ax.set_ylim(0.2, 4.2)
    ax.axis("off")
    ax.set_title("Persistent identifiers follow the archive object graph", pad=10)
    save(fig, "77-accession-graph", output)


def release_layers(output: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.8))
    layers = [
        ("Primary archive", "Study · Sample · Experiment · Run", COLORS["blue"]),
        ("Derived scientific data", "Assembly · MAG · UViG · catalog · profiles", COLORS["green"]),
        ("Reproducible computation", "Code release · parameters · container digest · databases", COLORS["purple"]),
        ("Paper contract", "Methods · Data Availability · Code Availability", COLORS["orange"]),
    ]
    y_positions = [3.7, 2.7, 1.7, 0.7]
    for (title, content, color), y in zip(layers, y_positions, strict=True):
        patch = FancyBboxPatch((0.7, y - 0.35), 8.6, 0.7, boxstyle="round,pad=0.04",
                               edgecolor=color, facecolor=color + "20", linewidth=1.6)
        ax.add_patch(patch)
        ax.text(1.0, y, title, va="center", weight="bold", color=color, fontsize=10)
        ax.text(4.05, y, content, va="center", color=COLORS["dark"], fontsize=9.2)
    for y1, y2 in zip(y_positions[:-1], y_positions[1:], strict=True):
        arrow(ax, (5.0, y1 - 0.37), (5.0, y2 + 0.37))
    ax.text(9.65, 2.2, "Links must\nremain resolvable", ha="center", va="center", color=COLORS["red"], weight="bold")
    ax.set_xlim(0.4, 10.4)
    ax.set_ylim(0.15, 4.25)
    ax.axis("off")
    ax.set_title("A paper is reproducible only when all four release layers connect")
    save(fig, "77-release-layers", output)


def artifact_readiness(readiness: pd.DataFrame, output: Path) -> None:
    order = readiness.sort_values(["StatusCode", "ReleaseObject"], ascending=[True, True]).reset_index(drop=True)
    colors = {
        "Missing": COLORS["red"], "Blocked": COLORS["red"], "Local only": COLORS["orange"],
        "Draft ready": COLORS["orange"], "Local package ready": COLORS["blue"],
        "Existing third-party": COLORS["green"],
    }
    fig, ax = plt.subplots(figsize=(9.4, 5.4))
    y = np.arange(len(order))
    ax.hlines(y, 0, order["StatusCode"], color="#D7DEE3", linewidth=3)
    ax.scatter(order["StatusCode"], y, s=95, c=[colors[value] for value in order["Status"]], edgecolor="white", linewidth=0.8, zorder=3)
    for i, row in order.iterrows():
        ax.text(row["StatusCode"] + 0.08, i, row["Status"], va="center", fontsize=7.5, color=colors[row["Status"]])
    ax.set_yticks(y, order["ReleaseObject"])
    ax.set_xticks([0, 1, 2, 3], ["Blocked / missing", "Local / draft", "Package ready", "Existing archive"])
    ax.set_xlim(-0.05, 3.55)
    ax.set_xlabel("Evidence state (not a quality score)")
    ax.set_title("The worked example is locally auditable but not externally releasable")
    ax.grid(axis="x", color="#E7ECEF")
    save(fig, "77-artifact-readiness", output)


def identifier_state(identifiers: pd.DataFrame, output: Path) -> None:
    counts = identifiers["State"].value_counts().reindex(
        ["EXISTING_THIRD_PARTY", "LOCAL_ONLY", "RESOLVE_BEFORE_RELEASE", "BLOCKED"], fill_value=0
    )
    labels = ["Existing third-party", "Local provenance", "Resolve before release", "Pending / blocked"]
    colors = [COLORS["green"], COLORS["blue"], COLORS["orange"], COLORS["red"]]
    fig, ax = plt.subplots(figsize=(8.7, 4.6))
    bars = ax.bar(labels, counts.values, color=colors, width=0.68)
    ax.bar_label(bars, padding=4, fontsize=11, weight="bold")
    ax.set_ylabel("Identifier records")
    ax.set_title("Existing identifiers are separated from identifiers that do not yet exist")
    ax.tick_params(axis="x", rotation=12)
    ax.grid(axis="y", color="#E7ECEF")
    ax.text(0.5, -0.26, "PENDING_NOT_INVENTED is a blocking token, never a citation", transform=ax.transAxes,
            ha="center", color=COLORS["red"], weight="bold")
    save(fig, "77-identifier-state", output)


def container_digest(output: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.6))
    nodes = [
        ((1.0, 2.5), "Build recipe\n+ source commit", COLORS["blue"], 2.0),
        ((3.6, 2.5), "Mutable tag\nregistry/image:version", COLORS["orange"], 2.0),
        ((6.3, 2.5), "Index digest\nsha256:...", COLORS["green"], 1.8),
        ((8.8, 3.25), "linux/amd64\nmanifest digest", COLORS["purple"], 1.8),
        ((8.8, 1.75), "linux/arm64\nmanifest digest", COLORS["purple"], 1.8),
    ]
    for pos, label, color, width in nodes:
        box(ax, pos, label, color, width=width, height=0.8)
    for start, end in [
        ((2.0, 2.5), (2.55, 2.5)), ((4.6, 2.5), (5.35, 2.5)),
        ((7.2, 2.65), (7.9, 3.05)), ((7.2, 2.35), (7.9, 1.95)),
    ]:
        arrow(ax, start, end)
    ax.text(3.6, 0.65, "Tag can move", ha="center", color=COLORS["red"], weight="bold")
    ax.text(7.55, 0.65, "Record index + tested-platform digest", ha="center", color=COLORS["green"], weight="bold")
    ax.set_xlim(-0.2, 10.0)
    ax.set_ylim(0.2, 4.1)
    ax.axis("off")
    ax.set_title("Pin container content, not only a mutable image tag")
    save(fig, "77-container-digest", output)


def database_manifest(databases: pd.DataFrame, output: Path) -> None:
    fields = ["Release", "Artifact", "Checksum", "SourceURL"]
    matrix = databases[fields].notna().astype(int).values
    labels = [name.replace(" reference database", "").replace(" reference package", "") for name in databases["Database"]]
    fig, ax = plt.subplots(figsize=(8.8, 4.2))
    ax.imshow(matrix, cmap=ListedColormap([COLORS["red"], COLORS["green"]]), vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(fields)), ["Release", "Artifact", "Checksum", "Source URL"])
    ax.set_yticks(range(len(labels)), labels)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, "Recorded" if matrix[i, j] else "Missing", ha="center", va="center",
                    color="white", weight="bold", fontsize=8.5)
    ax.set_title("A public database manifest excludes machine-local paths")
    ax.set_xlabel("Required provenance field")
    ax.set_ylabel("Reference database")
    save(fig, "77-database-manifest", output)


def release_gates(gates: pd.DataFrame, output: Path) -> None:
    status_order = ["Pass", "Review", "Not applicable", "Blocked"]
    colors = {"Pass": COLORS["green"], "Review": COLORS["orange"], "Not applicable": COLORS["gray"], "Blocked": COLORS["red"]}
    counts = gates["Status"].value_counts().reindex(status_order, fill_value=0)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), gridspec_kw={"width_ratios": [0.9, 1.55]})
    bars = axes[0].bar(status_order, counts.values, color=[colors[x] for x in status_order])
    axes[0].bar_label(bars, padding=3, weight="bold")
    axes[0].set_ylabel("Release gates")
    axes[0].tick_params(axis="x", rotation=20)
    axes[0].grid(axis="y", color="#E7ECEF")
    order = gates.sort_values(["StatusCode", "GateOrder"]).reset_index(drop=True)
    y = np.arange(len(order))
    axes[1].scatter(order["StatusCode"], y, c=[colors[x] for x in order["Status"]], s=45)
    axes[1].set_yticks(y, order["Gate"], fontsize=7)
    axes[1].set_xticks([0, 1, 2, 3], ["Blocked", "Review", "N/A", "Pass"])
    axes[1].set_xlim(-0.2, 3.2)
    axes[1].grid(axis="x", color="#E7ECEF")
    fig.suptitle("External release remains blocked until every mandatory gate passes", weight="bold", fontsize=13)
    fig.tight_layout()
    save(fig, "77-release-gates", output)


def availability_statement(data: pd.DataFrame, code: pd.DataFrame, output: Path) -> None:
    status_code = {"Missing": 0, "Blocked": 0, "Pending": 1, "Local only": 1, "Local package ready": 2, "Complete": 3}
    data_status = data.set_index("Component")["Status"].to_dict()
    code_status = code.set_index("Component")["Status"].to_dict()
    columns = [
        "Reused data", "Generated sequences", "Processed tables", "Restrictions",
        "Repository URL", "Version ID", "License", "Container digest", "Database manifest",
    ]
    statuses = [
        [data_status["Reused study"], data_status["MAGs"], data_status["Figure source data"], data_status["Restrictions"], "Not applicable", data_status["Gene catalog"], "Not applicable", "Not applicable", "Not applicable"],
        ["Not applicable", "Not applicable", "Not applicable", "Not applicable", code_status["Repository URL"], code_status["Version DOI"], code_status["License"], code_status["Container digest"], code_status["Database manifest"]],
    ]
    matrix = np.full((2, len(columns)), np.nan)
    for i in range(2):
        for j, status in enumerate(statuses[i]):
            if status != "Not applicable":
                matrix[i, j] = status_code[status]
    cmap = ListedColormap([COLORS["red"], COLORS["orange"], COLORS["blue"], COLORS["green"]])
    fig, ax = plt.subplots(figsize=(11.2, 4.3))
    masked = np.ma.masked_invalid(matrix)
    ax.imshow(masked, cmap=cmap, vmin=0, vmax=3, aspect="auto")
    ax.set_xticks(range(len(columns)), columns, rotation=55, ha="right", fontsize=7)
    ax.set_yticks([0, 1], ["Data Availability", "Code Availability"])
    for i in range(2):
        for j in range(len(columns)):
            if statuses[i][j] != "Not applicable":
                ax.text(j, i, statuses[i][j].replace(" ", "\n"), ha="center", va="center", fontsize=5.8,
                        color="white" if matrix[i, j] in (0, 3) else COLORS["dark"], weight="bold")
    ax.set_title("Availability statements must cover reused data, generated data, code and immutable execution records")
    ax.set_xlabel("Statement component")
    save(fig, "77-availability-statement", output)


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output = args.figures_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    np.random.seed(SEED)
    setup()
    routes = read(input_dir, "repository-routing-matrix.tsv")
    readiness = read(input_dir, "release-readiness.tsv")
    identifiers = read(input_dir, "identifier-registry.tsv")
    databases = read(input_dir, "database-manifest-public.tsv")
    gates = read(input_dir, "release-gate-ledger.tsv")
    data = read(input_dir, "data-availability-components.tsv")
    code = read(input_dir, "code-availability-components.tsv")

    repository_routing(routes, output)
    accession_graph(output)
    release_layers(output)
    artifact_readiness(readiness, output)
    identifier_state(identifiers, output)
    container_digest(output)
    database_manifest(databases, output)
    release_gates(gates, output)
    availability_statement(data, code, output)

    anchor_source = input_dir / "tenhoopen-figure2-original.jpg"
    anchor_target = output / "77-tenhoopen-figure2-original.jpg"
    shutil.copy2(anchor_source, anchor_target)
    files = [output / f"{stem}.{extension}" for stem in STEMS for extension in ("png", "svg")]
    files.append(anchor_target)
    manifest = {
        "article": 77,
        "plot_seed": SEED,
        "figures": [
            {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in files
        ],
    }
    (input_dir / "figure-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"figures\t{output}\t{len(files)} files")


if __name__ == "__main__":
    main()
