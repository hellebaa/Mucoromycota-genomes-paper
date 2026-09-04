#!/usr/bin/env python3
"""Summarize MAT-locus and sex-pathway blastp hits from protein databases."""

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
    "qlen",
    "slen",
    "stitle",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_query_metadata(mat_catalog: Path, sex_catalog: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    for row in read_tsv(mat_catalog):
        metadata[row["query_id"]] = {
            "evidence_track": row.get("evidence_track", ""),
            "locus_role": row.get("locus_role", ""),
            "gene_symbol": row.get("gene_symbol", ""),
            "mating_type_or_allele": row.get("mating_type_or_allele", ""),
            "source_species": row.get("source_species", ""),
            "priority_tier": row.get("priority_tier", ""),
            "use_for_mat_calling": row.get("use_for_mat_calling", ""),
            "notes": row.get("notes", ""),
        }
    for row in read_tsv(sex_catalog):
        metadata[row["query_id"]] = {
            "evidence_track": row.get("evidence_track", ""),
            "locus_role": row.get("pathway", ""),
            "gene_symbol": row.get("gene_symbol", ""),
            "mating_type_or_allele": "",
            "source_species": row.get("source_species", ""),
            "priority_tier": row.get("priority_tier", ""),
            "use_for_mat_calling": row.get("use_for_mat_calling", ""),
            "notes": row.get("notes", ""),
        }
    return metadata


def parse_blast_file(path: Path, genome_id: str) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=BLAST_FIELDS)
        for row in reader:
            row["genome_id"] = genome_id
            rows.append(row)
    return rows


def best_rows(rows: list[dict[str, str]], key_fields: tuple[str, ...]) -> list[dict[str, str]]:
    best: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(row[field] for field in key_fields)
        bitscore = float(row["bitscore"])
        if key not in best or bitscore > float(best[key]["bitscore"]):
            best[key] = row
    return [best[key] for key in sorted(best)]


def add_query_metadata(row: dict[str, str], query_metadata: dict[str, dict[str, str]]) -> dict[str, str]:
    annotated = dict(row)
    query = query_metadata.get(row["qseqid"], {})
    for field in [
        "evidence_track",
        "locus_role",
        "gene_symbol",
        "mating_type_or_allele",
        "source_species",
        "priority_tier",
        "use_for_mat_calling",
        "notes",
    ]:
        annotated[field] = query.get(field, "")
    return annotated


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blastdb-manifest", required=True, type=Path)
    parser.add_argument("--blastp-dir", required=True, type=Path)
    parser.add_argument("--mat-query-catalog", required=True, type=Path)
    parser.add_argument("--sex-pathway-catalog", required=True, type=Path)
    parser.add_argument("--out-all", required=True, type=Path)
    parser.add_argument("--out-best-per-query", required=True, type=Path)
    parser.add_argument("--out-best-per-gene", required=True, type=Path)
    args = parser.parse_args()

    query_metadata = load_query_metadata(args.mat_query_catalog, args.sex_pathway_catalog)
    manifest_rows = read_tsv(args.blastdb_manifest)
    all_rows: list[dict[str, str]] = []

    for manifest in manifest_rows:
        genome_id = manifest["genome_id"]
        blast_path = args.blastp_dir / f"{genome_id}.mat_locus_plus_sex_pathway.blastp.tsv"
        all_rows.extend(parse_blast_file(blast_path, genome_id))

    annotated_rows = [add_query_metadata(row, query_metadata) for row in all_rows]
    best_per_query = best_rows(annotated_rows, ("genome_id", "qseqid"))
    best_per_gene = best_rows(annotated_rows, ("genome_id", "gene_symbol", "evidence_track"))

    output_fields = [
        "genome_id",
        "evidence_track",
        "locus_role",
        "gene_symbol",
        "mating_type_or_allele",
        "priority_tier",
        "use_for_mat_calling",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "qcovs",
        "qlen",
        "slen",
        "evalue",
        "bitscore",
        "qstart",
        "qend",
        "sstart",
        "send",
        "source_species",
        "stitle",
        "notes",
    ]

    write_rows(args.out_all, annotated_rows, output_fields)
    write_rows(args.out_best_per_query, best_per_query, output_fields)
    write_rows(args.out_best_per_gene, best_per_gene, output_fields)

    print(f"Wrote all hits: {args.out_all} ({len(annotated_rows)} rows)")
    print(f"Wrote best per query: {args.out_best_per_query} ({len(best_per_query)} rows)")
    print(f"Wrote best per gene: {args.out_best_per_gene} ({len(best_per_gene)} rows)")


if __name__ == "__main__":
    main()
