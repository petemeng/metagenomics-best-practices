#!/usr/bin/env python3
"""Run the Article 34 Prodigal, MMseqs2, CD-HIT, and truth-audit workflow."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterator, TextIO


SEED = 20260734
PRIMARY_IDENTITY = 0.95
PRIMARY_COVERAGE = 0.95
BRANCHES = {
    "megahit-m1": ("MEGAHIT", "Individual", "MOCK1"),
    "megahit-m2": ("MEGAHIT", "Individual", "MOCK2"),
    "megahit-co": ("MEGAHIT", "Co-assembly", "MOCK1+MOCK2"),
    "metaspades-m1": ("metaSPAdes", "Individual", "MOCK1"),
    "metaspades-m2": ("metaSPAdes", "Individual", "MOCK2"),
    "metaspades-co": ("metaSPAdes", "Co-assembly", "MOCK1+MOCK2"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            elif header is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_record(handle: TextIO, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def deterministic_gzip_text(path: Path) -> TextIO:
    raw = path.open("wb")
    gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(gz, encoding="utf-8", newline="")


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = deterministic_gzip_text(path) if path.suffix == ".gz" else path.open("w", encoding="utf-8", newline="")
    with handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def pool_fasta(inputs: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for source in inputs:
            for header, sequence in fasta_records(source):
                identifier = header.split()[0]
                if identifier in seen:
                    raise ValueError(f"Duplicate identifier {identifier} while creating {output}")
                seen.add(identifier)
                write_record(handle, header, sequence)
    if not seen:
        raise ValueError(f"Empty FASTA pool: {output}")


def extract_fasta(source: Path, identifiers: set[str], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    found: set[str] = set()
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in fasta_records(source):
            identifier = header.split()[0]
            if identifier in identifiers:
                write_record(handle, header, sequence)
                found.add(identifier)
    missing = identifiers - found
    if missing:
        raise ValueError(f"Failed to extract {len(missing)} representatives into {output}")


def split_truth_fasta(source: Path, output_dir: Path) -> list[tuple[str, Path]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, list[tuple[str, str]]] = {}
    for header, sequence in fasta_records(source):
        identifier = header.split()[0]
        fields = identifier.split("__")
        if len(fields) < 3 or fields[0] != "truth":
            raise ValueError(f"Unexpected normalized truth identifier: {identifier}")
        grouped.setdefault(fields[1], []).append((header, sequence))
    assets: list[tuple[str, Path]] = []
    for accession in sorted(grouped):
        output = output_dir / f"{accession}.fna"
        with output.open("w", encoding="utf-8", newline="\n") as handle:
            for header, sequence in grouped[accession]:
                write_record(handle, header, sequence)
        assets.append((accession, output))
    if len(assets) != 87:
        raise ValueError(f"Expected 87 truth genome FASTAs, found {len(assets)}")
    return assets


def combine_gff(inputs: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##gff-version  3\n")
        for source in inputs:
            for raw in source.read_text(encoding="utf-8").splitlines():
                if raw and not raw.startswith("#"):
                    handle.write(raw + "\n")


def fasta_ids(path: Path) -> set[str]:
    return {header.split()[0] for header, _ in fasta_records(path)}


class Runner:
    def __init__(self, work: Path, env_prefix: Path, threads: int, force: bool):
        self.work = work
        self.env_prefix = env_prefix
        self.threads = threads
        self.force = force
        self.logs = work / "logs"
        self.markers = work / "markers"
        self.logs.mkdir(parents=True, exist_ok=True)
        self.markers.mkdir(parents=True, exist_ok=True)
        self.records: list[dict] = []
        self.environment = os.environ.copy()
        self.environment.update(
            {
                "LC_ALL": "C",
                "LANG": "C",
                "OMP_NUM_THREADS": str(threads),
                "PYTHONHASHSEED": "0",
                "PYTHONNOUSERSITE": "1",
            }
        )

    def binary(self, name: str) -> Path:
        path = self.env_prefix / "bin" / name
        if not path.is_file():
            raise FileNotFoundError(path)
        return path

    def run(self, label: str, command: list[str], outputs: list[Path], inputs: list[Path] | None = None) -> None:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "-", label)
        marker = self.markers / f"{safe}.json"
        payload = {
            "label": label,
            "command": command,
            "inputs": {str(path): sha256(path) for path in (inputs or []) if path.is_file()},
        }
        cached = False
        if not self.force and marker.is_file() and all(path.is_file() and path.stat().st_size > 0 for path in outputs):
            previous = json.loads(marker.read_text(encoding="utf-8"))
            cached = previous.get("contract") == payload
        if not cached:
            for output in outputs:
                output.parent.mkdir(parents=True, exist_ok=True)
            stdout_path = self.logs / f"{safe}.stdout.log"
            stderr_path = self.logs / f"{safe}.stderr.log"
            time_path = self.logs / f"{safe}.time.txt"
            timed_command = ["/usr/bin/time", "-v", "-o", str(time_path), *command]
            with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
                subprocess.run(timed_command, check=True, env=self.environment, stdout=stdout, stderr=stderr)
            missing = [str(path) for path in outputs if not path.is_file() or path.stat().st_size == 0]
            if missing:
                raise RuntimeError(f"Step {label} did not create nonempty outputs: {missing}")
            marker.write_text(
                json.dumps(
                    {
                        "contract": payload,
                        "outputs": {str(path): {"bytes": path.stat().st_size, "sha256": sha256(path)} for path in outputs},
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
        self.records.append(
            {
                "Step": label,
                "Cached": "yes" if cached else "no",
                "Command": shlex.join(command),
                "Outputs": ";".join(str(path.relative_to(self.work)) for path in outputs),
            }
        )


def parse_gene_metadata(fna: Path, faa: Path, branch: str) -> list[dict]:
    aa_lengths = {header.split()[0]: len(sequence.rstrip("*")) for header, sequence in fasta_records(faa)}
    rows: list[dict] = []
    for header, sequence in fasta_records(fna):
        parts = header.split(" # ")
        identifier = parts[0].split()[0]
        if len(parts) < 5:
            raise ValueError(f"Unexpected Prodigal header: {header}")
        partial_match = re.search(r"(?:^|;)partial=([01]{2})(?:;|$)", parts[4])
        if partial_match is None:
            raise ValueError(f"No Prodigal partial code in header: {header}")
        partial = partial_match.group(1)
        completeness = "Complete" if partial == "00" else "Incomplete" if partial == "11" else "Partial"
        contig = identifier.rsplit("_", 1)[0]
        if branch == "truth":
            assembler, strategy, mock = "Exact reference", "Truth", "MOCK2"
            genome = contig.split("__")[1]
        else:
            assembler, strategy, mock = BRANCHES[branch]
            genome = ""
        rows.append(
            {
                "GeneID": identifier,
                "ContigID": contig,
                "Branch": branch,
                "Assembler": assembler,
                "OriginStrategy": strategy,
                "Mock": mock,
                "GenBankAssembly": genome,
                "Start": int(parts[1]),
                "End": int(parts[2]),
                "Strand": "+" if parts[3] == "1" else "-",
                "PartialCode": partial,
                "Completeness": completeness,
                "NtLength": len(sequence),
                "AaLength": aa_lengths[identifier],
            }
        )
    if set(aa_lengths) != {row["GeneID"] for row in rows}:
        raise ValueError(f"Protein/nucleotide identifiers differ for {branch}")
    return rows


def mmseq_cluster(
    runner: Runner,
    label: str,
    input_faa: Path,
    output_dir: Path,
    identity: float,
    coverage: float,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    db = output_dir / "sequenceDB"
    clu = output_dir / "clusterDB"
    rep_db = output_dir / "representativeDB"
    tmp = output_dir / "tmp"
    tmp.mkdir(exist_ok=True)
    membership = output_dir / "membership.tsv"
    reps_faa = output_dir / "representatives.faa"
    mmseqs = str(runner.binary("mmseqs"))
    runner.run(
        f"{label}-createdb",
        [mmseqs, "createdb", str(input_faa), str(db)],
        [db, Path(str(db) + ".index")],
        [input_faa],
    )
    cluster_label = f"{label}-cluster"
    cluster_marker = runner.markers / f"{re.sub(r'[^A-Za-z0-9_.-]+', '-', cluster_label)}.json"
    if runner.force or not cluster_marker.is_file():
        for path in output_dir.glob(f"{clu.name}*"):
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir()
    runner.run(
        cluster_label,
        [
            mmseqs,
            "cluster",
            str(db),
            str(clu),
            str(tmp),
            "--min-seq-id",
            f"{identity:.2f}",
            "-c",
            f"{coverage:.2f}",
            "--cov-mode",
            "1",
            "--cluster-mode",
            "2",
            "--alignment-mode",
            "3",
            "--threads",
            str(runner.threads),
        ],
        [Path(str(clu) + ".index"), Path(str(clu) + ".dbtype")],
        [input_faa],
    )
    runner.run(
        f"{label}-createtsv",
        [mmseqs, "createtsv", str(db), str(db), str(clu), str(membership)],
        [membership],
    )
    runner.run(
        f"{label}-createsubdb",
        [mmseqs, "createsubdb", str(clu), str(db), str(rep_db)],
        [Path(str(rep_db) + ".index"), Path(str(rep_db) + ".dbtype")],
    )
    runner.run(
        f"{label}-createsubdb-headers",
        [mmseqs, "createsubdb", str(clu), str(db) + "_h", str(rep_db) + "_h"],
        [Path(str(rep_db) + "_h.index"), Path(str(rep_db) + "_h.dbtype")],
    )
    runner.run(
        f"{label}-convert2fasta",
        [mmseqs, "convert2fasta", str(rep_db), str(reps_faa)],
        [reps_faa],
    )
    return {"db": db, "clu": clu, "rep_db": rep_db, "membership": membership, "faa": reps_faa}


def direct_catalog(source_faa: Path, source_fna: Path, output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    faa = output_dir / "representatives.faa"
    fna = output_dir / "representatives.fna"
    membership = output_dir / "membership.tsv"
    shutil.copyfile(source_faa, faa)
    shutil.copyfile(source_fna, fna)
    with membership.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier in sorted(fasta_ids(faa)):
            handle.write(f"{identifier}\t{identifier}\n")
    return {"faa": faa, "fna": fna, "membership": membership}


def parse_cdhit_membership(path: Path, output: Path) -> None:
    clusters: list[list[tuple[str, bool]]] = []
    current: list[tuple[str, bool]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">Cluster"):
            if current:
                clusters.append(current)
            current = []
            continue
        match = re.search(r">([^\.]+)\.\.\.", raw)
        if match:
            current.append((match.group(1), raw.rstrip().endswith("*")))
    if current:
        clusters.append(current)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for cluster in clusters:
            representatives = [identifier for identifier, is_rep in cluster if is_rep]
            if len(representatives) != 1:
                raise ValueError(f"CD-HIT cluster has {len(representatives)} representatives")
            representative = representatives[0]
            for member, _ in cluster:
                handle.write(f"{representative}\t{member}\n")


def make_pools(work: Path, assembler: str, individual_reps: dict[str, Path] | None = None, variant: str = "primary") -> dict[str, Path]:
    key = assembler.lower().replace("spades", "spades")
    prefix = "megahit" if assembler == "MEGAHIT" else "metaspades"
    pool_dir = work / "pools" / key
    pool_dir.mkdir(parents=True, exist_ok=True)
    individual_faa = pool_dir / "individual.raw.faa"
    individual_fna = pool_dir / "individual.raw.fna"
    co_faa = pool_dir / "co.raw.faa"
    co_fna = pool_dir / "co.raw.fna"
    pool_fasta([work / f"genes/{prefix}-m1.faa", work / f"genes/{prefix}-m2.faa"], individual_faa)
    pool_fasta([work / f"genes/{prefix}-m1.fna", work / f"genes/{prefix}-m2.fna"], individual_fna)
    pool_fasta([work / f"genes/{prefix}-co.faa"], co_faa)
    pool_fasta([work / f"genes/{prefix}-co.fna"], co_fna)
    result = {"individual_faa": individual_faa, "individual_fna": individual_fna, "co_faa": co_faa, "co_fna": co_fna}
    if individual_reps is not None:
        mix_faa = pool_dir / f"mix-stage2-{variant}.faa"
        mix_fna = pool_dir / f"mix-stage2-{variant}.fna"
        pool_fasta([individual_reps["faa"], co_faa], mix_faa)
        pool_fasta([individual_reps["fna"], co_fna], mix_fna)
        result.update({"mix_faa": mix_faa, "mix_fna": mix_fna})
    return result


def add_contract(
    rows: list[dict],
    work: Path,
    catalog_id: str,
    assembler: str,
    strategy: str,
    method: str,
    identity: float,
    coverage: float,
    cascaded: bool,
    assets: dict[str, Path],
    input_faa: Path,
    stage1_membership: Path | None = None,
) -> None:
    rows.append(
        {
            "CatalogID": catalog_id,
            "Assembler": assembler,
            "Strategy": strategy,
            "Method": method,
            "MinimumAminoAcidIdentity": f"{identity:.2f}",
            "MinimumMemberCoverage": f"{coverage:.2f}",
            "CoverageMode": "target/member",
            "ClusterMode": "greedy incremental by length" if method != "No clustering" else "not applicable",
            "Cascaded": "yes" if cascaded else "no",
            "InputFAA": str(input_faa.relative_to(work)),
            "RepresentativesFAA": str(assets["faa"].relative_to(work)),
            "RepresentativesFNA": str(assets["fna"].relative_to(work)),
            "MembershipTSV": str(assets["membership"].relative_to(work)),
            "Stage1MembershipTSV": "" if stage1_membership is None else str(stage1_membership.relative_to(work)),
        }
    )


def create_catalogs(runner: Runner) -> list[dict]:
    work = runner.work
    contracts: list[dict] = []
    for assembler in ("MEGAHIT", "metaSPAdes"):
        pools = make_pools(work, assembler)
        stem = assembler.lower()
        individual_dir = work / f"catalogs/{stem}/individual-primary"
        individual = mmseq_cluster(runner, f"{stem}-individual-primary", pools["individual_faa"], individual_dir, 0.95, 0.95)
        individual["fna"] = individual_dir / "representatives.fna"
        extract_fasta(pools["individual_fna"], fasta_ids(individual["faa"]), individual["fna"])
        add_contract(
            contracts,
            work,
            f"{stem}-individual-primary",
            assembler,
            "Individual",
            "MMseqs2 9.d36de",
            0.95,
            0.95,
            False,
            individual,
            pools["individual_faa"],
        )

        co_dir = work / f"catalogs/{stem}/co-primary"
        co = direct_catalog(pools["co_faa"], pools["co_fna"], co_dir)
        add_contract(
            contracts,
            work,
            f"{stem}-co-primary",
            assembler,
            "Co-assembly",
            "No clustering",
            1.0,
            1.0,
            False,
            co,
            pools["co_faa"],
        )

        pools = make_pools(work, assembler, individual, "primary")
        mix_dir = work / f"catalogs/{stem}/mix-primary"
        mix = mmseq_cluster(runner, f"{stem}-mix-primary", pools["mix_faa"], mix_dir, 0.95, 0.95)
        mix["fna"] = mix_dir / "representatives.fna"
        extract_fasta(pools["mix_fna"], fasta_ids(mix["faa"]), mix["fna"])
        add_contract(
            contracts,
            work,
            f"{stem}-mix-primary",
            assembler,
            "Mix",
            "MMseqs2 9.d36de",
            0.95,
            0.95,
            True,
            mix,
            pools["mix_faa"],
            individual["membership"],
        )

    for variant, identity, coverage in (("id90-cov80", 0.90, 0.80), ("id99-cov95", 0.99, 0.95)):
        pools = make_pools(work, "MEGAHIT")
        ind_dir = work / f"catalogs/megahit/individual-{variant}"
        individual = mmseq_cluster(runner, f"megahit-individual-{variant}", pools["individual_faa"], ind_dir, identity, coverage)
        individual["fna"] = ind_dir / "representatives.fna"
        extract_fasta(pools["individual_fna"], fasta_ids(individual["faa"]), individual["fna"])
        pools = make_pools(work, "MEGAHIT", individual, variant)
        mix_dir = work / f"catalogs/megahit/mix-{variant}"
        mix = mmseq_cluster(runner, f"megahit-mix-{variant}", pools["mix_faa"], mix_dir, identity, coverage)
        mix["fna"] = mix_dir / "representatives.fna"
        extract_fasta(pools["mix_fna"], fasta_ids(mix["faa"]), mix["fna"])
        add_contract(
            contracts,
            work,
            f"megahit-mix-{variant}",
            "MEGAHIT",
            "Mix",
            "MMseqs2 9.d36de",
            identity,
            coverage,
            True,
            mix,
            pools["mix_faa"],
            individual["membership"],
        )

    primary = next(row for row in contracts if row["CatalogID"] == "megahit-mix-primary")
    cdhit_dir = work / "catalogs/megahit/mix-cdhit"
    cdhit_dir.mkdir(parents=True, exist_ok=True)
    input_faa = work / primary["InputFAA"]
    input_fna = work / "pools/megahit/mix-stage2-primary.fna"
    reps_faa = cdhit_dir / "representatives.faa"
    clstr = Path(str(reps_faa) + ".clstr")
    runner.run(
        "megahit-mix-cdhit",
        [
            str(runner.binary("cd-hit")),
            "-i",
            str(input_faa),
            "-o",
            str(reps_faa),
            "-c",
            "0.95",
            "-G",
            "0",
            "-aS",
            "0.95",
            "-g",
            "1",
            "-n",
            "5",
            "-d",
            "0",
            "-T",
            str(runner.threads),
            "-M",
            "0",
        ],
        [reps_faa, clstr],
        [input_faa],
    )
    membership = cdhit_dir / "membership.tsv"
    parse_cdhit_membership(clstr, membership)
    reps_fna = cdhit_dir / "representatives.fna"
    extract_fasta(input_fna, fasta_ids(reps_faa), reps_fna)
    add_contract(
        contracts,
        work,
        "megahit-mix-cdhit",
        "MEGAHIT",
        "Mix",
        "CD-HIT 4.8.1 local identity",
        0.95,
        0.95,
        True,
        {"faa": reps_faa, "fna": reps_fna, "membership": membership},
        input_faa,
        work / primary["Stage1MembershipTSV"],
    )
    return contracts


def search_truth(runner: Runner, contracts: list[dict], truth: dict[str, Path]) -> None:
    mmseqs = str(runner.binary("mmseqs"))
    truth_db = truth["db"]
    jobs: list[tuple[str, Path, Path, Path, Path, Path, str]] = []
    for contract in contracts:
        catalog_id = contract["CatalogID"]
        audit_dir = runner.work / "truth-audit" / catalog_id
        audit_dir.mkdir(parents=True, exist_ok=True)
        faa = runner.work / contract["RepresentativesFAA"]
        catalog_db = audit_dir / "catalogDB"
        runner.run(
            f"truth-{catalog_id}-createdb",
            [mmseqs, "createdb", str(faa), str(catalog_db)],
            [catalog_db, Path(str(catalog_db) + ".index")],
            [faa],
        )
        for direction, coverage_mode in (
            ("catalog-to-truth-recovery", "1"),
            ("catalog-to-truth-support", "2"),
        ):
            jobs.append(
                (
                    catalog_id,
                    audit_dir,
                    catalog_db,
                    audit_dir / f"{direction}DB",
                    audit_dir / f"tmp-{direction}",
                    audit_dir / f"{direction}.tsv",
                    coverage_mode,
                )
            )

    workers = min(2, len(jobs))
    default_job_threads = max(1, runner.threads // workers)

    def run_search(job: tuple[str, Path, Path, Path, Path, Path, str]) -> None:
        catalog_id, audit_dir, query_db, result_db, tmp, table, coverage_mode = job
        direction = table.stem
        search_label = f"truth-{catalog_id}-{direction}-search"
        search_marker = runner.markers / f"{re.sub(r'[^A-Za-z0-9_.-]+', '-', search_label)}.json"
        job_threads = default_job_threads
        if search_marker.is_file():
            previous = json.loads(search_marker.read_text(encoding="utf-8"))
            previous_command = previous.get("contract", {}).get("command", [])
            if "--threads" in previous_command:
                job_threads = int(previous_command[previous_command.index("--threads") + 1])
        if runner.force or not search_marker.is_file():
            for path in audit_dir.glob(f"{result_db.name}*"):
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
            if tmp.exists():
                shutil.rmtree(tmp)
            tmp.mkdir()
        runner.run(
            search_label,
            [
                mmseqs,
                "search",
                str(query_db),
                str(truth_db),
                str(result_db),
                str(tmp),
                "--min-seq-id",
                "0.95",
                "-c",
                "0.95",
                "--cov-mode",
                coverage_mode,
                "--alignment-mode",
                "3",
                "--max-seqs",
                "100",
                "--threads",
                str(job_threads),
            ],
            [Path(str(result_db) + ".index"), Path(str(result_db) + ".dbtype")],
        )
        runner.run(
            f"truth-{catalog_id}-{direction}-convertalis",
            [
                mmseqs,
                "convertalis",
                str(query_db),
                str(truth_db),
                str(result_db),
                str(table),
                "--format-output",
                "query,target,pident,alnlen,qcov,tcov,evalue,bits",
            ],
            [table],
        )

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        list(executor.map(run_search, jobs))


def tool_version(command: list[str]) -> str:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    return (result.stdout + "\n" + result.stderr).strip().replace("\n", " | ")[:500]


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    env_prefix = args.env_prefix.resolve()
    work = args.work_dir.resolve()
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if not (work / ".article34-inputs-complete").is_file():
        raise FileNotFoundError("Run prepare_article34_gene_catalog_inputs.py first")
    runner = Runner(work, env_prefix, args.threads, args.force)

    genes = work / "genes"
    genes.mkdir(parents=True, exist_ok=True)
    all_metadata: list[dict] = []
    prodigal = str(runner.binary("prodigal"))

    def call_genes(branch: str) -> list[dict]:
        source = work / ("truth/mock2-exact-genomes.fna" if branch == "truth" else f"assemblies/{branch}.fna")
        faa = genes / f"{branch}.faa"
        fna = genes / f"{branch}.fna"
        gff = genes / f"{branch}.gff"
        runner.run(
            f"prodigal-{branch}",
            [prodigal, "-i", str(source), "-a", str(faa), "-d", str(fna), "-o", str(gff), "-f", "gff", "-p", "meta", "-q"],
            [faa, fna, gff],
            [source],
        )
        return parse_gene_metadata(fna, faa, branch)

    branches = [*BRANCHES]
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(len(branches), args.threads)) as executor:
        for rows in executor.map(call_genes, branches):
            all_metadata.extend(rows)

    truth_inputs = split_truth_fasta(work / "truth/mock2-exact-genomes.fna", work / "truth/genomes")
    truth_parts = work / "genes/truth-parts"
    truth_parts.mkdir(parents=True, exist_ok=True)

    def call_truth(accession_and_source: tuple[str, Path]) -> tuple[Path, Path, Path]:
        accession, source = accession_and_source
        faa = truth_parts / f"{accession}.faa"
        fna = truth_parts / f"{accession}.fna"
        gff = truth_parts / f"{accession}.gff"
        runner.run(
            f"prodigal-truth-{accession}",
            [prodigal, "-i", str(source), "-a", str(faa), "-d", str(fna), "-o", str(gff), "-f", "gff", "-p", "meta", "-q"],
            [faa, fna, gff],
            [source],
        )
        return faa, fna, gff

    with concurrent.futures.ThreadPoolExecutor(max_workers=min(len(truth_inputs), args.threads)) as executor:
        truth_outputs = list(executor.map(call_truth, truth_inputs))
    pool_fasta([assets[0] for assets in truth_outputs], genes / "truth.faa")
    pool_fasta([assets[1] for assets in truth_outputs], genes / "truth.fna")
    combine_gff([assets[2] for assets in truth_outputs], genes / "truth.gff")
    all_metadata.extend(parse_gene_metadata(genes / "truth.fna", genes / "truth.faa", "truth"))

    assembly_metadata = [row for row in all_metadata if row["Branch"] != "truth"]
    truth_metadata = [row for row in all_metadata if row["Branch"] == "truth"]
    write_tsv(work / "summary/gene-metadata.tsv.gz", assembly_metadata)
    write_tsv(work / "summary/truth-gene-metadata.tsv.gz", truth_metadata)

    contracts = create_catalogs(runner)
    write_tsv(work / "summary/catalog-contracts.tsv", contracts)

    truth_cluster = mmseq_cluster(
        runner,
        "truth-nr-primary",
        genes / "truth.faa",
        work / "truth/nonredundant-primary",
        PRIMARY_IDENTITY,
        PRIMARY_COVERAGE,
    )
    search_truth(runner, contracts, {"db": truth_cluster["db"], "rep_db": truth_cluster["rep_db"]})

    versions = [
        {"Tool": "Python", "Version": tool_version([str(runner.binary("python")), "--version"])},
        {"Tool": "Prodigal", "Version": tool_version([prodigal, "-v"])},
        {"Tool": "MMseqs2", "Version": tool_version([str(runner.binary("mmseqs")), "version"])},
        {"Tool": "CD-HIT", "Version": tool_version([str(runner.binary("cd-hit")), "-h"])},
        {"Tool": "SeqKit", "Version": tool_version([str(runner.binary("seqkit")), "version"])},
    ]
    write_tsv(work / "summary/tool-versions.tsv", versions)
    write_tsv(work / "summary/command-log.tsv", sorted(runner.records, key=lambda row: row["Step"]))
    contract = {
        "seed": SEED,
        "threads": args.threads,
        "environment_prefix": str(env_prefix),
        "environment_yaml_sha256": sha256(root / "env/gene-catalog.yml"),
        "environment_lock_sha256": sha256(root / "env/gene-catalog-linux-64.lock"),
        "prodigal_mode": "meta",
        "primary_clustering": {
            "software": "MMseqs2 9.d36de",
            "minimum_amino_acid_identity": 0.95,
            "minimum_member_coverage": 0.95,
            "coverage_mode": 1,
            "cluster_mode": 2,
            "alignment_mode": 3,
        },
        "truth_audit": {
            "identity": 0.95,
            "target_coverage": 0.95,
            "alignment_mode": 3,
            "max_seqs": 100,
            "denominator": "Prodigal-meta ORF clusters predicted from 87 exact MOCK2 reference genomes",
        },
        "catalogs": len(contracts),
    }
    (work / "summary/run-contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    (work / ".article34-run-complete").write_text(json.dumps(contract, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(contract, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as error:
        print(f"Command failed with exit code {error.returncode}: {error.cmd}", file=sys.stderr)
        raise
