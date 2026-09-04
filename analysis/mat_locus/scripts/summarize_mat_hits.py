#!/usr/bin/env python3
"""Summarize MAT-locus BLAST hits across genome scans."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


BLAST_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovs",
    "stitle",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organized-samples", required=True, type=Path)
    parser.add_argument("--query-map", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_query_map(path: Path) -> dict[str, dict[str, str]]:
    queries = {}
    for row in read_tsv(path):
        queries[row["query_id"]] = row
    return queries


def read_blast(path: Path, search_type: str) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=BLAST_FIELDS)
        for row in reader:
            row["search_type"] = search_type
            rows.append(row)
    return rows


def best_hits(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    best = {}
    for row in rows:
        query = row["qseqid"]
        bitscore = float(row["bitscore"])
        if query not in best or bitscore > float(best[query]["bitscore"]):
            best[query] = row
    return best


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    samples = read_tsv(args.organized_samples)
    query_map = load_query_map(args.query_map)

    output_fields = [
        "genome_id",
        "tree_label",
        "species_name",
        "search_type",
        "query_id",
        "gene_symbol",
        "query_mating_type",
        "subject_id",
        "percent_identity",
        "alignment_length",
        "query_coverage",
        "evalue",
        "bitscore",
        "subject_start",
        "subject_end",
        "subject_title",
    ]

    output_rows = []
    for sample in samples:
        genome_id = sample["genome_id"]
        rows = []
        rows.extend(read_blast(args.results_dir / "blastp" / f"{genome_id}.blastp.tsv", "blastp"))
        rows.extend(read_blast(args.results_dir / "tblastn" / f"{genome_id}.tblastn.tsv", "tblastn"))
        for query_id, hit in sorted(best_hits(rows).items()):
            query = query_map.get(query_id, {})
            output_rows.append(
                {
                    "genome_id": genome_id,
                    "tree_label": sample.get("tree_label", ""),
                    "species_name": sample.get("species_name", ""),
                    "search_type": hit["search_type"],
                    "query_id": query_id,
                    "gene_symbol": query.get("gene_symbol", ""),
                    "query_mating_type": query.get("mating_type", ""),
                    "subject_id": hit["sseqid"],
                    "percent_identity": hit["pident"],
                    "alignment_length": hit["length"],
                    "query_coverage": hit["qcovs"],
                    "evalue": hit["evalue"],
                    "bitscore": hit["bitscore"],
                    "subject_start": hit["sstart"],
                    "subject_end": hit["send"],
                    "subject_title": hit["stitle"],
                }
            )

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {args.output}")
    print(f"Rows: {len(output_rows)}")


if __name__ == "__main__":
    main()
