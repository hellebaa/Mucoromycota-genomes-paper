#!/usr/bin/env python3
"""Scan FASTA records for exact motif occurrences."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def fasta_records(path: Path):
    record_id = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, "".join(chunks).upper()
                record_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if record_id is not None:
        yield record_id, "".join(chunks).upper()


def find_all(seq: str, motif: str) -> list[int]:
    positions = []
    start = 0
    while True:
        hit = seq.find(motif, start)
        if hit == -1:
            return positions
        positions.append(hit)
        start = hit + 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--motif", default="GCATAGGT")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    motif = args.motif.upper()
    rows = []
    for record_id, seq in fasta_records(args.fasta):
        positions = find_all(seq, motif)
        for pos in positions:
            rows.append(
                {
                    "query_id": record_id,
                    "motif": motif,
                    "motif_start_0based": pos,
                    "motif_end_0based": pos + len(motif),
                    "query_length_bp": len(seq),
                    "relative_position": round(pos / len(seq), 6) if seq else "",
                }
            )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        fieldnames = ["query_id", "motif", "motif_start_0based", "motif_end_0based", "query_length_bp", "relative_position"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote motif hits: {args.out} ({len(rows)} hits)")


if __name__ == "__main__":
    main()
