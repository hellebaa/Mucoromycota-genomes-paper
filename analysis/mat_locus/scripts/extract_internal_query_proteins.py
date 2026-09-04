#!/usr/bin/env python3
"""Extract internal query proteins listed in the MAT query catalog."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())
    if header:
        records.append((header, "".join(seq)))
    return records


def find_manifest_row(manifest_rows: list[dict[str, str]], genome_id: str) -> dict[str, str] | None:
    for row in manifest_rows:
        if row["genome_id"] == genome_id:
            return row
    return None


def find_record(records: list[tuple[str, str]], protein_id: str) -> tuple[str, str] | None:
    for header, seq in records:
        if protein_id == header.split()[0]:
            return header, seq
        fields = header.replace("|", " ").replace(";", " ").split()
        if protein_id in fields or protein_id in header:
            return header, seq
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-catalog", required=True, type=Path)
    parser.add_argument("--genome-manifest", required=True, type=Path)
    parser.add_argument("--out-fasta", required=True, type=Path)
    parser.add_argument("--missing-out", required=True, type=Path)
    args = parser.parse_args()

    catalog_rows = read_tsv(args.query_catalog)
    manifest_rows = read_tsv(args.genome_manifest)
    args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.missing_out.parent.mkdir(parents=True, exist_ok=True)

    protein_cache: dict[Path, list[tuple[str, str]]] = {}
    written = 0
    missing_rows = []

    with args.out_fasta.open("w") as out:
        for query in catalog_rows:
            if query["sequence_source"] != "internal_protein_id":
                continue

            genome_id = query["source_genome_id"]
            protein_id = query["protein_id"]
            manifest_row = find_manifest_row(manifest_rows, genome_id)
            if manifest_row is None or not manifest_row.get("proteins_faa"):
                missing_rows.append({**query, "missing_reason": "missing proteome in manifest"})
                continue

            protein_path = Path(manifest_row["proteins_faa"])
            if protein_path not in protein_cache:
                protein_cache[protein_path] = read_fasta(protein_path)

            record = find_record(protein_cache[protein_path], protein_id)
            if record is None:
                missing_rows.append({**query, "missing_reason": "protein_id not found in proteome"})
                continue

            original_header, seq = record
            out.write(f">{query['query_id']} {query['gene_symbol']} {query['source_species']} original={original_header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")
            written += 1

    with args.missing_out.open("w", newline="") as handle:
        fieldnames = list(catalog_rows[0]) + ["missing_reason"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(missing_rows)

    print(f"Wrote: {args.out_fasta}")
    print(f"Internal query proteins written: {written}")
    print(f"Missing/internal unresolved: {len(missing_rows)}")


if __name__ == "__main__":
    main()
