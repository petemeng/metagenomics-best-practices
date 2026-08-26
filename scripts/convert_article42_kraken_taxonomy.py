#!/usr/bin/env python3
"""Convert Kraken2 contig calls to the seven-rank TaxVAMB input contract."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from article41_44_utils import write_tsv


RANKS = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kraken-output", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--names", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def load_nodes(path: Path) -> tuple[dict[int, int], dict[int, str]]:
    parent, rank = {}, {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = [field.strip() for field in raw.split("|")]
            taxid = int(fields[0])
            parent[taxid] = int(fields[1])
            rank[taxid] = fields[2]
    return parent, rank


def load_names(path: Path) -> dict[int, str]:
    names = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = [field.strip() for field in raw.split("|")]
            if len(fields) >= 4 and fields[3] == "scientific name":
                names[int(fields[0])] = fields[1].replace(";", ",").replace("\t", " ")
    return names


def lineage(taxid: int, parents: dict[int, int], ranks: dict[int, str], names: dict[int, str]) -> tuple[str, str]:
    if taxid <= 0 or taxid not in parents:
        return "", "unclassified"
    selected: dict[str, str] = {}
    seen = set()
    current = taxid
    while current not in seen and current in parents:
        seen.add(current)
        rank = ranks.get(current, "")
        # NCBI Taxonomy renamed the old "superkingdom" rank to "domain".
        # TaxVAMB's canonical seven-rank contract still uses that first slot.
        if rank == "domain":
            rank = "superkingdom"
        if rank in RANKS and current in names:
            selected[rank] = names[current]
        next_taxid = parents[current]
        if next_taxid == current:
            break
        current = next_taxid
    values = []
    deepest = "unclassified"
    for rank in RANKS:
        if rank not in selected:
            break
        values.append(selected[rank])
        deepest = "domain" if rank == "superkingdom" else rank
    if not values:
        return "", "unclassified"
    return ";".join(values), deepest


def main() -> int:
    args = parse_args()
    parents, ranks = load_nodes(args.nodes)
    names = load_names(args.names)
    rows = []
    status = Counter()
    with args.kraken_output.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Malformed Kraken2 row: {raw[:200]}")
            call, contig, taxid_text = fields[:3]
            taxid = int(taxid_text) if call == "C" else 0
            prediction, deepest = lineage(taxid, parents, ranks, names)
            rows.append({"contigs": contig, "predictions": prediction})
            status[(call, deepest)] += 1
    if not rows:
        raise ValueError("No Kraken2 contig calls")
    write_tsv(args.output, rows)
    write_tsv(
        args.summary,
        [
            {"KrakenStatus": call, "DeepestRank": deepest, "Contigs": count}
            for (call, deepest), count in sorted(status.items())
        ],
    )
    print(f"PASS TaxVAMB taxonomy: {len(rows)} contigs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
