#!/usr/bin/env python3
"""Map query proteins to local proteomes and OrthoFinder groups.

This is a lightweight pre-BLAST helper for cases where close homologs are
expected, such as the Mucor circinelloides MAT query proteins against local
Mucor circinelloides proteomes. Use BLAST/DIAMOND for the full cross-species
scan.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument(
        "--proteome",
        action="append",
        required=True,
        help="Proteome as label:path. Can be supplied more than once.",
    )
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--hogs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--kmer-size", default=6, type=int)
    parser.add_argument("--top-candidates", default=20, type=int)
    return parser.parse_args()


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    seq_parts: list[str] = []
    with path.open(errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def normalize_gene_id(header_or_token: str) -> str:
    gene = header_or_token.split()[0]
    for suffix in ("-T1", "-T2", "-T3", "-T4", "-T5"):
        gene = gene.replace(suffix, "")
    return gene


def load_group_lookup(path: Path, group_column: str) -> dict[str, str]:
    lookup = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            group = row[group_column]
            for column, value in row.items():
                if column in {group_column, "OG", "Gene Tree Parent Clade"} or not value:
                    continue
                for token in value.split(", "):
                    first = token.split()[0]
                    lookup[first] = group
                    lookup[normalize_gene_id(first)] = group
    return lookup


def shared_kmers(query: str, subject: str, kmer_size: int) -> int:
    kmers = {
        query[i : i + kmer_size]
        for i in range(len(query) - kmer_size + 1)
        if len(set(query[i : i + kmer_size])) >= 3
    }
    return sum(1 for kmer in kmers if kmer in subject)


def smith_waterman_score(query: str, subject: str) -> int:
    previous = [0] * (len(subject) + 1)
    best = 0
    for query_aa in query:
        current = [0]
        for j, subject_aa in enumerate(subject, 1):
            score = max(
                0,
                previous[j - 1] + (2 if query_aa == subject_aa else -1),
                previous[j] - 2,
                current[j - 1] - 2,
            )
            current.append(score)
            best = max(best, score)
        previous = current
    return best


def parse_proteome_arg(value: str) -> tuple[str, Path]:
    if ":" not in value:
        raise ValueError(f"Proteome argument must be label:path, got {value}")
    label, path = value.split(":", 1)
    return label, Path(path)


def main() -> None:
    args = parse_args()
    og_lookup = load_group_lookup(args.orthogroups, "Orthogroup")
    hog_lookup = load_group_lookup(args.hogs, "HOG")
    queries = parse_fasta(args.query_fasta)
    proteomes = [(label, parse_fasta(path)) for label, path in map(parse_proteome_arg, args.proteome)]

    fields = [
        "query_id",
        "query_description",
        "query_length",
        "proteome",
        "best_gene_id",
        "best_header",
        "best_length",
        "shared_kmers",
        "smith_waterman_score",
        "orthogroup",
        "hierarchical_orthogroup",
    ]
    rows = []
    for query_header, query_seq in queries:
        query_id, _, query_description = query_header.partition(" ")
        for label, records in proteomes:
            candidates = []
            for subject_header, subject_seq in records:
                kmer_count = shared_kmers(query_seq, subject_seq, args.kmer_size)
                if kmer_count:
                    candidates.append((kmer_count, subject_header, subject_seq))
            candidates.sort(reverse=True, key=lambda item: item[0])
            best = None
            for kmer_count, subject_header, subject_seq in candidates[: args.top_candidates]:
                score = smith_waterman_score(query_seq, subject_seq)
                if best is None or score > best[0]:
                    best = (score, kmer_count, subject_header, subject_seq)

            if best is None:
                rows.append(
                    {
                        "query_id": query_id,
                        "query_description": query_description,
                        "query_length": len(query_seq),
                        "proteome": label,
                        "best_gene_id": "",
                        "best_header": "",
                        "best_length": "",
                        "shared_kmers": 0,
                        "smith_waterman_score": 0,
                        "orthogroup": "",
                        "hierarchical_orthogroup": "",
                    }
                )
                continue

            score, kmer_count, subject_header, subject_seq = best
            gene_id = normalize_gene_id(subject_header)
            rows.append(
                {
                    "query_id": query_id,
                    "query_description": query_description,
                    "query_length": len(query_seq),
                    "proteome": label,
                    "best_gene_id": gene_id,
                    "best_header": subject_header,
                    "best_length": len(subject_seq),
                    "shared_kmers": kmer_count,
                    "smith_waterman_score": score,
                    "orthogroup": og_lookup.get(gene_id, ""),
                    "hierarchical_orthogroup": hog_lookup.get(gene_id, ""),
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.output}")
    print(f"Rows: {len(rows)}")


if __name__ == "__main__":
    main()
