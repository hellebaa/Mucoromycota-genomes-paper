#!/usr/bin/env python3
"""Build combined MAT-locus and sex-pathway protein query FASTA files."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header:
        records.append((header, "".join(seq)))
    return records


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mucor-mat-fasta",
        required=True,
        type=Path,
        help="Historical Mucor circinelloides MAT protein FASTA.",
    )
    parser.add_argument(
        "--sex-pathway-fasta",
        required=True,
        type=Path,
        help="Protein FASTA containing TSP1/TSP2 or other sex-pathway queries.",
    )
    parser.add_argument(
        "--internal-mat-fasta",
        type=Path,
        help="Optional internal MAT candidate protein FASTA extracted from the genome mirror.",
    )
    parser.add_argument(
        "--out-mat-fasta",
        required=True,
        type=Path,
        help="Output MAT-locus query protein FASTA.",
    )
    parser.add_argument(
        "--out-all-fasta",
        required=True,
        type=Path,
        help="Output combined MAT-locus plus sex-pathway protein FASTA.",
    )
    args = parser.parse_args()

    mat_records = read_fasta(args.mucor_mat_fasta)
    if args.internal_mat_fasta and args.internal_mat_fasta.exists():
        mat_records.extend(read_fasta(args.internal_mat_fasta))

    sex_pathway_records = read_fasta(args.sex_pathway_fasta)
    all_records = mat_records + sex_pathway_records

    write_fasta(mat_records, args.out_mat_fasta)
    write_fasta(all_records, args.out_all_fasta)

    print(f"Wrote MAT-locus queries: {args.out_mat_fasta} ({len(mat_records)} records)")
    print(f"Wrote all protein queries: {args.out_all_fasta} ({len(all_records)} records)")


if __name__ == "__main__":
    main()
