#!/usr/bin/env python3
"""Prepare Article 74's release, parameter, database, and runtime audit packet."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import re
import shutil
import tarfile
from pathlib import Path

import pandas as pd


ARTICLE = 74
ANALYSIS_SEED = 74_001
PLOT_SEED = 20_260_774
SNAPSHOT_DATE = "2026-08-24"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--runtime-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n", float_format="%.10g")


def extract_member(archive: Path, member: str, destination: Path) -> None:
    with tarfile.open(archive, "r:gz") as handle:
        source = handle.extractfile(member)
        if source is None:
            raise FileNotFoundError(f"Missing {member} in {archive}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(source.read())


def read_fastq_records(path: Path) -> int:
    lines = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for lines, _ in enumerate(handle, start=1):
            pass
    if lines % 4:
        raise ValueError(f"Incomplete FASTQ record in {path}")
    return lines // 4


def parse_duration(value: object) -> float:
    text = str(value).strip()
    total = 0.0
    for number, unit in re.findall(r"([0-9.]+)\s*(ms|us|ns|d|h|m|s)", text):
        factor = {
            "d": 86400.0,
            "h": 3600.0,
            "m": 60.0,
            "s": 1.0,
            "ms": 1e-3,
            "us": 1e-6,
            "ns": 1e-9,
        }[unit]
        total += float(number) * factor
    return total


def parse_memory_mib(value: object) -> float:
    text = str(value).strip()
    match = re.fullmatch(r"([0-9.]+)\s*([KMGT]?B)", text, flags=re.I)
    if not match:
        return 0.0
    number, unit = match.groups()
    factor = {"B": 1 / 2**20, "KB": 1 / 2**10, "MB": 1, "GB": 2**10, "TB": 2**20}[unit.upper()]
    return float(number) * factor


def select_traces(runtime: Path) -> tuple[Path, Path]:
    traces = []
    for path in sorted((runtime / "stub-output/pipeline_info").glob("execution_trace_*.txt")):
        try:
            frame = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        traces.append((path, frame))
    first = [pair for pair in traces if len(pair[1]) >= 5 and not pair[1]["status"].eq("CACHED").any()]
    resumed = [pair for pair in traces if pair[1]["status"].eq("CACHED").any()]
    if not first or not resumed:
        raise ValueError("Could not identify successful first and resumed traces")
    return first[-1][0], resumed[-1][0]


def stable_release_records(evidence: Path) -> dict[str, object]:
    mag_release = json.loads((evidence / "mag-release.json").read_text())
    mag_tag = json.loads((evidence / "mag-tag.json").read_text())
    func_release = json.loads((evidence / "funcscan-release.json").read_text())
    func_tag = json.loads((evidence / "funcscan-tag.json").read_text())
    nf_release = json.loads((evidence / "nextflow-release.json").read_text())
    return {
        "nf-core/mag": {
            "release": mag_release["tag_name"],
            "name": mag_release["name"],
            "published_at": mag_release["published_at"],
            "commit": mag_tag["object"]["sha"],
            "source_url": "https://github.com/nf-core/mag/releases/tag/5.5.0",
        },
        "nf-core/funcscan": {
            "release": func_release["tag_name"],
            "name": func_release["name"],
            "published_at": func_release["published_at"],
            "commit": func_tag["object"]["sha"],
            "source_url": "https://github.com/nf-core/funcscan/releases/tag/4.0.0",
        },
        "Nextflow": {
            "release": nf_release["tag_name"].removeprefix("v"),
            "name": nf_release["name"],
            "published_at": nf_release["published_at"],
            "source_url": "https://github.com/nextflow-io/nextflow/releases/tag/v26.04.0",
            "assets": [
                {
                    "name": item["name"],
                    "bytes": item["size"],
                    "digest": item.get("digest", ""),
                    "url": item["browser_download_url"],
                }
                for item in nf_release["assets"]
                if item["name"] in {"nextflow", "nextflow-26.04.0-dist"}
            ],
        },
    }


def prepare_sources(evidence: Path, output: Path) -> dict[str, object]:
    manifest = json.loads((evidence / "download-manifest.json").read_text())
    if manifest["article"] != ARTICLE or manifest["resource_count"] != 13:
        raise ValueError("Unexpected Article 74 download manifest")
    for name, record in manifest["resources"].items():
        path = evidence / name
        if path.stat().st_size != record["bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"Source checksum mismatch: {name}")

    source = output / "source"
    source.mkdir(parents=True, exist_ok=True)
    mag_archive = evidence / "mag-5.5.0.tar.gz"
    func_archive = evidence / "funcscan-4.0.0.tar.gz"
    selected = {
        "mag-nextflow.config": (mag_archive, "mag-5.5.0/nextflow.config"),
        "mag-input-schema.json": (mag_archive, "mag-5.5.0/assets/schema_input.json"),
        "mag-test-minimal.config": (mag_archive, "mag-5.5.0/conf/test_minimal.config"),
        "mag-resource-guidance.md": (mag_archive, "mag-5.5.0/docs/usage/resource_guidance.md"),
        "mag-usage.md": (mag_archive, "mag-5.5.0/docs/usage.md"),
        "mag-output.md": (mag_archive, "mag-5.5.0/docs/output.md"),
        "mag-changelog.md": (mag_archive, "mag-5.5.0/CHANGELOG.md"),
        "mag-metromap-original.png": (mag_archive, "mag-5.5.0/docs/images/mag_metromap_light.png"),
        "funcscan-nextflow.config": (func_archive, "funcscan-4.0.0/nextflow.config"),
        "funcscan-schema.json": (func_archive, "funcscan-4.0.0/nextflow_schema.json"),
        "funcscan-readme.md": (func_archive, "funcscan-4.0.0/README.md"),
        "funcscan-usage.md": (func_archive, "funcscan-4.0.0/docs/usage.md"),
        "funcscan-output.md": (func_archive, "funcscan-4.0.0/docs/output.md"),
        "funcscan-metromap-original.png": (func_archive, "funcscan-4.0.0/docs/images/funcscan_metro_workflow.png"),
    }
    for name, (archive, member) in selected.items():
        extract_member(archive, member, source / name)

    for name in (
        "checkm2-record.json",
        "gunc-progenomes2.1.dmnd.gz.md5",
        "gunc-progenomes2.1.dmnd.md5",
        "gunc-database-v1.1.0.py",
    ):
        shutil.copy2(evidence / name, source / name)

    stable = stable_release_records(evidence)
    (source / "release-records.stable.json").write_text(
        json.dumps(stable, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    source_records = {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(source.iterdir())
        if path.is_file()
    }
    normalized_manifest = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "download_manifest_sha256": sha256(evidence / "download-manifest.json"),
        "immutable_downloads": manifest["resources"],
        "selected_source_files": source_records,
    }
    (output / "source-manifest.json").write_text(
        json.dumps(normalized_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    shutil.copy2(source / "mag-metromap-original.png", output / "mag-metromap-original.png")
    shutil.copy2(source / "funcscan-metromap-original.png", output / "funcscan-metromap-original.png")
    return stable


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    evidence = args.evidence_dir.resolve()
    runtime = args.runtime_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    releases = prepare_sources(evidence, output)
    release_rows = []
    for component, record in releases.items():
        archive_name = {
            "nf-core/mag": "mag-5.5.0.tar.gz",
            "nf-core/funcscan": "funcscan-4.0.0.tar.gz",
            "Nextflow": "nextflow-26.04.0-dist",
        }[component]
        archive = evidence / archive_name
        release_rows.append(
            {
                "Component": component,
                "Release": record["release"],
                "PublishedUTC": record["published_at"],
                "Commit": record.get("commit", "release asset"),
                "Artifact": archive_name,
                "Bytes": archive.stat().st_size,
                "SHA256": sha256(archive),
                "SourceURL": record["source_url"],
            }
        )
    release_lock = pd.DataFrame(release_rows)
    write_tsv(release_lock, output / "release-lock.tsv")

    checkm2 = json.loads((evidence / "checkm2-record.json").read_text())
    databases = pd.DataFrame(
        [
            {
                "Database": "GTDB-Tk reference package",
                "Release": "R11-RS232",
                "Artifact": "gtdbtk_r232_data.tar.gz",
                "Bytes": 60806405195,
                "ChecksumType": "MD5",
                "Checksum": "25a59e0352b1fd150c589f56559767d4",
                "SourceURL": "https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz",
                "LocalPath": "/shared/db/gtdb/release232/gtdbtk_r232_data.tar.gz",
            },
            {
                "Database": "CheckM2 reference database",
                "Release": "Zenodo 14897628; dataset version 3",
                "Artifact": checkm2["files"][0]["key"],
                "Bytes": checkm2["files"][0]["size"],
                "ChecksumType": "MD5",
                "Checksum": checkm2["files"][0]["checksum"].split(":", 1)[1],
                "SourceURL": checkm2["files"][0]["links"]["self"],
                "LocalPath": "/shared/db/checkm2/zenodo-14897628/checkm2_db_v14897628.dmnd",
            },
            {
                "Database": "GUNC ProGenomes",
                "Release": "ProGenomes 2.1; GUNC 1.1.0 download endpoint",
                "Artifact": "gunc_db_progenomes2.1.dmnd.gz",
                "Bytes": 7185710760,
                "ChecksumType": "MD5",
                "Checksum": "bc93a855e0760aad5c4e5f2d0e26da46",
                "SourceURL": "https://black.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz",
                "LocalPath": "/shared/db/gunc/progenomes2.1/gunc_db_progenomes2.1.dmnd",
            },
        ]
    )
    write_tsv(databases, output / "database-lock.tsv")

    clean_paths = [
        root / "data/raw/article13/ERR9765746_clean_R1.fastq.gz",
        root / "data/raw/article13/ERR9765746_clean_R2.fastq.gz",
    ]
    input_rows = []
    for mate, path in zip(("R1", "R2"), clean_paths, strict=True):
        input_rows.append(
            {
                "Project": "PRJEB52977",
                "BioSample": "SAMEA14435832",
                "Run": "ERR9765746",
                "Sample": "ERR9765746",
                "Group": "MOCK1",
                "Mate": mate,
                "Records": read_fastq_records(path),
                "RelativePath": path.relative_to(root).as_posix(),
                "Bytes": path.stat().st_size,
                "SHA256": sha256(path),
            }
        )
    input_manifest = pd.DataFrame(input_rows)
    if input_manifest["Records"].nunique() != 1 or input_manifest["Records"].iloc[0] != 99991:
        raise ValueError("Unexpected clean paired-read count")
    write_tsv(input_manifest, output / "input-manifest.tsv")
    shutil.copy2(runtime / "samplesheet.csv", output / "samplesheet.csv")
    shutil.copy2(runtime / "params.stub.yml", output / "params.stub.yml")
    shutil.copy2(runtime / "params.publication.yml", output / "params.publication.yml")
    shutil.copy2(runtime / "hpc.slurm.config", output / "hpc.slurm.config")
    shutil.copy2(runtime / "production-command.sh", output / "production-command.sh")

    defaults = pd.DataFrame(
        [
            ("coassemble_group", "false", "false", "One assembly per biological sample; group remains metadata"),
            ("binning_map_mode", "group", "own", "Keep the abundance/coverage unit aligned to each sample"),
            ("skip_spades", "false", "true", "Avoid silently comparing two assemblers inside the primary run"),
            ("skip_megahit", "false", "false", "MEGAHIT is the prespecified short-read assembler"),
            ("skip_metabat2", "false", "false", "MetaBAT2 is the prespecified primary binner"),
            ("skip_maxbin2", "false", "true", "Alternative binners require a separate benchmark plan"),
            ("skip_concoct", "false", "true", "Alternative binners require a separate benchmark plan"),
            ("run_busco", "true", "false", "Use the project MAG acceptance pair below"),
            ("run_checkm2", "false", "true", "Estimate completeness and contamination"),
            ("run_gunc", "false", "true", "Audit lineage inconsistency/chimerism independently"),
            ("gtdb_db", "GTDB R232 URL", "local checksum-locked R232", "Prevent implicit re-download or release drift"),
            ("metabat_rng_seed", "1", "1", "Preserve deterministic MetaBAT2 execution"),
            ("semibin_rng_seed", "1", "1", "Preserve deterministic setting even while SemiBin is disabled"),
        ],
        columns=["Parameter", "PipelineDefault5.5.0", "PublicationOverride", "Reason"],
    )
    write_tsv(defaults, output / "pipeline-defaults-audit.tsv")

    precedence = pd.DataFrame(
        [
            (1, "Pipeline script defaults", "nextflow.config", "Lowest"),
            (2, "Configuration files and profiles", "institution + -profile apptainer + -c hpc.slurm.config", "Infrastructure"),
            (3, "Parameter file", "-params-file params.publication.yml", "Scientific contract"),
            (4, "Command-line pipeline parameters", "--parameter value", "Highest; reserve for emergency overrides"),
        ],
        columns=["Rank", "Layer", "Example", "Role"],
    )
    write_tsv(precedence, output / "parameter-precedence.tsv")

    units = pd.DataFrame(
        [
            ("Separate assembly + own mapping", False, "own", "sample", "sample", "Independent samples; primary analysis"),
            ("Separate assembly + group mapping", False, "group", "sample", "all runs in group -> each assembly", "Longitudinal/group-aware coverage; justify cross-mapping"),
            ("Group co-assembly + group mapping", True, "group", "group", "all runs in group -> pooled assembly", "Prespecified shared-community recovery"),
            ("Group co-assembly + all mapping", True, "all", "group", "all runs -> every assembly", "Special comparative design; high compute and leakage risk"),
        ],
        columns=["Strategy", "coassemble_group", "binning_map_mode", "AssemblyUnit", "CoverageUnit", "UseCase"],
    )
    write_tsv(units, output / "execution-unit-matrix.tsv")

    profiles = pd.DataFrame(
        [
            ("apptainer", "Software isolation", "One engine only; shared immutable image cache", "Configuration parsed locally; engine not executed"),
            ("hpc.slurm.config", "Scheduler/resources", "SLURM executor, queue, ceilings, mounts", "Configuration contract only; no sbatch on QA host"),
            ("params.publication.yml", "Scientific choices", "Assembler, binner, QC, databases, seeds", "Version controlled and peer reviewed"),
            ("test_minimal", "Control-plane QA", "Input schema and a minimal DAG", "Executed with -stub-run; not a biological benchmark"),
        ],
        columns=["Layer", "Purpose", "LockedContent", "LocalEvidence"],
    )
    write_tsv(profiles, output / "profile-contract.tsv")

    profile_parse = pd.DataFrame(
        [
            ("manifest.version", "5.5.0"),
            ("manifest.nextflowVersion", "!>=26.04.0"),
            ("plugins", "nf-schema@2.7.2"),
            ("process.executor", "slurm"),
            ("process.queue", "compute"),
            ("process.resourceLimits.cpus", "32"),
            ("process.resourceLimits.memory", "220 GB"),
            ("apptainer.enabled", "true"),
            ("apptainer.cacheDir", "/shared/containers/nf-core-mag-5.5.0"),
            ("docker.enabled", "false"),
            ("singularity.enabled", "false"),
            ("conda.enabled", "false"),
        ],
        columns=["EffectiveConfigKey", "ParsedValue"],
    )
    write_tsv(profile_parse, output / "profile-parse-evidence.tsv")

    runtime_environment = pd.DataFrame(
        [
            ("Nextflow", "26.04.0 build 12031", "Official distribution asset; SHA256 in release-lock.tsv", "executed"),
            ("Java", "OpenJDK 17.0.18", "Temporary conda-forge runtime used only for control-plane QA", "executed"),
            ("nf-core/mag", "5.5.0; local Nextflow revision 71af2049ef", "Tag commit 56abab5b023ce953c9c43fe21090d156ad0e18af", "executed in stub mode"),
            ("nf-schema", "2.7.2", "Pipeline plugin lock", "executed"),
            ("Apptainer", "profile enabled", "No Apptainer executable on local QA host", "config parse only"),
            ("SLURM", "executor=slurm", "No sbatch executable on local QA host", "config parse only"),
        ],
        columns=["Component", "VersionOrSetting", "Identity", "LocalStatus"],
    )
    write_tsv(runtime_environment, output / "runtime-environment.tsv")

    first_path, resume_path = select_traces(runtime)
    runtime_rows = []
    for run_label, path in (("first-success", first_path), ("resume", resume_path)):
        frame = pd.read_csv(path, sep="\t")
        for row in frame.itertuples(index=False):
            runtime_rows.append(
                {
                    "Run": run_label,
                    "TaskID": int(row.task_id),
                    "Hash": row.hash,
                    "Process": row.name,
                    "Status": row.status,
                    "Exit": int(row.exit),
                    "DurationSeconds": parse_duration(row.duration),
                    "RealtimeSeconds": parse_duration(row.realtime),
                    "PeakRSSMiB": parse_memory_mib(row.peak_rss),
                }
            )
    runtime_trace = pd.DataFrame(runtime_rows)
    write_tsv(runtime_trace, output / "stub-runtime-trace.tsv")

    first = runtime_trace[runtime_trace["Run"].eq("first-success")]
    resumed = runtime_trace[runtime_trace["Run"].eq("resume")]
    runtime_summary = pd.DataFrame(
        [
            ("offline-schema", "FAILED_AS_EXPECTED", 0, 0, "Remote CheckM URL cannot be validated under NXF_OFFLINE; failure retained as a boundary test"),
            ("first-success", "COMPLETED", len(first), int(first["Status"].eq("CACHED").sum()), "Online URL validation; official test_minimal profile; -stub-run"),
            ("resume", "COMPLETED", len(resumed), int(resumed["Status"].eq("CACHED").sum()), "Same work directory and parameters; MultiQC rebuilt after aggregate inputs changed"),
        ],
        columns=["Attempt", "Outcome", "Tasks", "CachedTasks", "Interpretation"],
    )
    write_tsv(runtime_summary, output / "runtime-summary.tsv")

    stub_scope = pd.DataFrame(
        [
            ("Input samplesheet conforms to the locked schema", True, "nf-schema accepted sample/run/group and paired FASTQ paths"),
            ("Pinned engine can parse and launch nf-core/mag 5.5.0", True, "Nextflow 26.04.0; release archive checksum recorded"),
            ("Minimal DAG schedules and completes", True, "Five stub tasks completed in the first successful run"),
            ("Task cache can be reused with -resume", True, "Four upstream tasks cached; MultiQC rebuilt"),
            ("Apptainer image execution works on the QA host", False, "Profile parsed only; Apptainer is not installed on this host"),
            ("SLURM submission works on the QA host", False, "Configuration only; sbatch is absent"),
            ("Assembly, bins, CheckM2, GUNC, or GTDB results are biologically valid", False, "Stub outputs are zero-byte/sentinel files and cannot support scientific claims"),
        ],
        columns=["Claim", "SupportedByLocalRun", "EvidenceBoundary"],
    )
    write_tsv(stub_scope, output / "stub-scope.tsv")

    hardware = pd.DataFrame(
        [
            ("Local schema/DAG stub", "2-4", "15 GB ceiling", "5 GB", "<5 min after plugin cache", "Executed here; not a biological run"),
            ("Short-read per-sample MAG run", "16-32", ">=72 GB; 220 GB queue ceiling", ">=500 GB working space", "12-72 h per sample, cohort dependent", "Planning range; benchmark on a pilot sample"),
            ("GTDB-Tk classification", "10", "140 GB first request", ">=150 GB including expanded R232", "hours per MAG set", "Official nf-core/mag 5.5 resource request"),
            ("MEGAHIT", "8", "40 GB first request", "input-dependent scratch", "up to 16 h first request", "Official nf-core/mag 5.5 resource request"),
            ("CheckM2 or GUNC", "6", "36 GB first request", "1.74 GB / 7.19 GB compressed DB plus outputs", "up to 8 h first request", "Official nf-core/mag 5.5 resource request"),
        ],
        columns=["Mode", "Cores", "RAM", "Disk", "ExpectedTime", "Basis"],
    )
    write_tsv(hardware, output / "hardware-envelope.tsv")

    funcscan = pd.DataFrame(
        [
            ("ARG", "run_arg_screening", "ABRicate; AMRFinderPlus; fARGene; RGI; DeepARG", "hAMRonization + argNorm", "Lock each tool database and preserve ontology version"),
            ("AMP", "run_amp_screening", "ampir; Macrel; AMPlify; optional HMMER", "AMPcombi", "Consensus is not experimental activity validation"),
            ("BGC", "run_bgc_screening", "antiSMASH; BiG-SLiCE; DeepBGC; GECCO; optional HMMER", "comBGC", "Contig-edge truncation and minimum length affect recovery"),
            ("CAZyme/CGC", "run_cazyme_screening", "dbCAN", "CGC and substrate tables when GFF is supplied", "Release usage.md says run_cazyme_annotation, but the executable schema/config use run_cazyme_screening"),
        ],
        columns=["Branch", "ActivationFlag", "PrimaryTools", "AggregateOutput", "PublicationAudit"],
    )
    write_tsv(funcscan, output / "funcscan-branch-contract.tsv")

    bundle = pd.DataFrame(
        [
            ("Pipeline identity", "release tag + full commit + source archive SHA256", True),
            ("Engine identity", "Nextflow release asset SHA256 + Java major", True),
            ("Scientific parameters", "params.publication.yml", True),
            ("Infrastructure", "Apptainer profile + hpc.slurm.config", True),
            ("Input identity", "samplesheet + per-FASTQ SHA256", True),
            ("Database identity", "release, URL, bytes, checksum, local mount", True),
            ("Runtime lineage", "trace/report/timeline/DAG + command + work cache policy", True),
            ("Container identity", "SHA256 of every pulled SIF", False),
            ("Scientific acceptance", "non-stub output QC and predefined inclusion table", False),
        ],
        columns=["Layer", "RequiredArtifact", "PresentInLocalPacket"],
    )
    write_tsv(bundle, output / "provenance-bundle.tsv")

    methods_contract = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "pipeline": {"name": "nf-core/mag", "release": "5.5.0", "commit": releases["nf-core/mag"]["commit"]},
        "engine": {"name": "Nextflow", "release": "26.04.0", "java": "17.0.18"},
        "execution": {"production_profile": "apptainer + hpc.slurm.config", "local_qa": "test_minimal + -stub-run"},
        "unit": {"coassemble_group": False, "binning_map_mode": "own"},
        "determinism": {"metabat_rng_seed": 1, "semibin_rng_seed": 1},
        "database_releases": {"GTDB": "R11-RS232", "CheckM2": "Zenodo 14897628 v3", "GUNC": "ProGenomes 2.1"},
        "scientific_boundary": "Stub execution validates input, DAG, and cache behavior only; publication claims require non-stub HPC outputs and predefined biological acceptance.",
    }
    (output / "methods-contract.json").write_text(
        json.dumps(methods_contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    metrics = {
        "article": ARTICLE,
        "snapshot_date": SNAPSHOT_DATE,
        "analysis_seed": ANALYSIS_SEED,
        "plot_seed": PLOT_SEED,
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "input_pairs": int(input_manifest["Records"].iloc[0]),
        "first_success_tasks": int(len(first)),
        "resume_tasks": int(len(resumed)),
        "resume_cached_tasks": int(resumed["Status"].eq("CACHED").sum()),
        "resume_rebuilt_tasks": int(resumed["Status"].ne("CACHED").sum()),
        "database_lock_count": len(databases),
        "source_snapshot_count": len(json.loads((output / "source-manifest.json").read_text())["selected_source_files"]),
    }
    (output / "analysis-metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"prepared\t{output}\t{len(first)} first tasks\t{len(resumed)} resumed tasks")


if __name__ == "__main__":
    main()
