#!/usr/bin/env python3
"""Summarize contig lengths and gene-bearing status for genome assemblies."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable


GENE_FEATURE_PRIORITY = ("gene", "mRNA", "transcript", "CDS")


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_selected(row: dict[str, str], include: str) -> bool:
    if include == "all":
        return True
    value = row.get("include_core_analyses", "").strip().lower()
    if include == "core":
        return value == "yes"
    value = row.get("include_final_ultrametric_tree", "").strip().lower()
    if include == "final_tree":
        return value == "yes"
    raise ValueError(f"Unsupported include mode: {include}")


def parse_fasta_lengths(path: Path) -> list[dict[str, int | str]]:
    records: list[dict[str, int | str]] = []
    current_id: str | None = None
    current_len = 0

    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append({"contig_id": current_id, "length_bp": current_len})
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())

    if current_id is not None:
        records.append({"contig_id": current_id, "length_bp": current_len})

    return records


def parse_gff_gene_stats(path: Path) -> dict[str, dict[str, int]]:
    by_feature: dict[str, dict[str, dict[str, int]]] = {
        feature: defaultdict(lambda: {"feature_count": 0, "feature_bp": 0})
        for feature in GENE_FEATURE_PRIORITY
    }

    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            seqid, _, feature_type, start, end = parts[:5]
            if feature_type not in by_feature:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            span = abs(end_i - start_i) + 1
            by_feature[feature_type][seqid]["feature_count"] += 1
            by_feature[feature_type][seqid]["feature_bp"] += span

    for feature in GENE_FEATURE_PRIORITY:
        if by_feature[feature]:
            return dict(by_feature[feature])
    return {}


def nxx(lengths: list[int], fraction: float) -> tuple[int, int]:
    total = sum(lengths)
    if total == 0:
        return 0, 0
    threshold = total * fraction
    running = 0
    for index, length in enumerate(sorted(lengths, reverse=True), start=1):
        running += length
        if running >= threshold:
            return length, index
    return 0, 0


def infer_family_order(row: dict[str, str]) -> tuple[str, str]:
    notes = row.get("notes", "")
    family = row.get("family", "")
    order = row.get("order", "")
    if family or order:
        return family, order
    # Keep this intentionally light; the manifest may not carry taxonomy.
    match = re.search(r"family=([^;]+);order=([^;]+)", notes)
    if match:
        return match.group(1), match.group(2)
    return "", ""


def summarize_genome(row: dict[str, str]) -> tuple[list[dict[str, str | int | float]], dict[str, str | int | float]]:
    assembly = Path(row["assembly_fasta"])
    gff = Path(row["gff"])
    contigs = parse_fasta_lengths(assembly)
    gene_stats = parse_gff_gene_stats(gff)
    family, order = infer_family_order(row)

    contig_rows: list[dict[str, str | int | float]] = []
    for rank, contig in enumerate(sorted(contigs, key=lambda item: int(item["length_bp"]), reverse=True), start=1):
        contig_id = str(contig["contig_id"])
        length_bp = int(contig["length_bp"])
        stats = gene_stats.get(contig_id, {"feature_count": 0, "feature_bp": 0})
        gene_count = int(stats["feature_count"])
        gene_bp = int(stats["feature_bp"])
        gene_bearing = "yes" if gene_count > 0 else "no"
        contig_rows.append(
            {
                "genome_id": row["genome_id"],
                "species_name": row.get("species_name", ""),
                "source": row.get("source", ""),
                "family": family,
                "order": order,
                "contig_rank": rank,
                "contig_id": contig_id,
                "length_bp": length_bp,
                "length_mb": round(length_bp / 1_000_000, 6),
                "gene_feature_count": gene_count,
                "gene_feature_bp": gene_bp,
                "gene_bearing_contig": gene_bearing,
                "putative_chromosome_current": gene_bearing,
            }
        )

    lengths = [int(item["length_bp"]) for item in contigs]
    n50, l50 = nxx(lengths, 0.5)
    n90, l90 = nxx(lengths, 0.9)
    gene_bearing_contigs = sum(1 for item in contig_rows if item["gene_bearing_contig"] == "yes")
    assembly_bp = sum(lengths)
    gene_bearing_bp = sum(int(item["length_bp"]) for item in contig_rows if item["gene_bearing_contig"] == "yes")
    genome_row = {
        "genome_id": row["genome_id"],
        "species_name": row.get("species_name", ""),
        "source": row.get("source", ""),
        "family": family,
        "order": order,
        "assembly_bp": assembly_bp,
        "assembly_mb": round(assembly_bp / 1_000_000, 6),
        "n_contigs": len(contigs),
        "max_contig_bp": max(lengths) if lengths else 0,
        "n50_bp": n50,
        "l50": l50,
        "n90_bp": n90,
        "l90": l90,
        "gene_bearing_contigs": gene_bearing_contigs,
        "gene_bearing_bp": gene_bearing_bp,
        "gene_bearing_assembly_percent": round(100 * gene_bearing_bp / assembly_bp, 4) if assembly_bp else 0,
    }
    return contig_rows, genome_row


def write_tsv(path: Path, rows: Iterable[dict[str, str | int | float]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("/cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv"))
    parser.add_argument("--out-contigs", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/assembly_contig_gene_summary.tsv"))
    parser.add_argument("--out-genomes", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/assembly_genome_summary.tsv"))
    parser.add_argument("--include", choices=["core", "final_tree", "all"], default="core")
    args = parser.parse_args()

    manifest = [row for row in read_manifest(args.manifest) if is_selected(row, args.include)]
    contig_rows: list[dict[str, str | int | float]] = []
    genome_rows: list[dict[str, str | int | float]] = []
    skipped: list[str] = []

    for row in manifest:
        assembly = row.get("assembly_fasta", "")
        gff = row.get("gff", "")
        if not assembly or not gff or not Path(assembly).exists() or not Path(gff).exists():
            skipped.append(row.get("genome_id", ""))
            continue
        genome_contigs, genome_summary = summarize_genome(row)
        contig_rows.extend(genome_contigs)
        genome_rows.append(genome_summary)

    if not contig_rows:
        raise SystemExit("No contig rows produced. Check manifest paths and include mode.")

    write_tsv(args.out_contigs, contig_rows, list(contig_rows[0].keys()))
    write_tsv(args.out_genomes, genome_rows, list(genome_rows[0].keys()))

    print(f"Wrote contig table: {args.out_contigs} ({len(contig_rows)} rows)")
    print(f"Wrote genome table: {args.out_genomes} ({len(genome_rows)} rows)")
    if skipped:
        print(f"Skipped genomes with missing assembly/GFF: {len(skipped)}")
        print(";".join(skipped))


if __name__ == "__main__":
    main()
