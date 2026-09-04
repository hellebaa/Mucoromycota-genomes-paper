#!/usr/bin/env python3
"""Scan candidate centromere BLAST-hit intervals for a diagnostic motif."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def read_manifest(path: Path) -> dict[str, Path]:
    assemblies: dict[str, Path] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            assembly = row.get("assembly_fasta", "")
            if assembly:
                assemblies[row["genome_id"]] = Path(assembly)
    return assemblies


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current_id: str | None = None
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                sequences[current_id] = []
            elif current_id is not None:
                sequences[current_id].append(line)
    return {key: "".join(value).upper() for key, value in sequences.items()}


def find_all(seq: str, motif: str) -> list[int]:
    positions = []
    start = 0
    while True:
        hit = seq.find(motif, start)
        if hit == -1:
            return positions
        positions.append(hit)
        start = hit + 1


def parse_bed_name(name: str) -> tuple[str, str, str]:
    parts = name.split("|", 2)
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    return "", "", name


def read_bed(path: Path) -> list[dict[str, str | int | float]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            genome_id, query_class, qseqid = parse_bed_name(parts[3])
            rows.append(
                {
                    "contig_id": parts[0],
                    "start_0based": int(parts[1]),
                    "end_0based": int(parts[2]),
                    "name": parts[3],
                    "score": parts[4] if len(parts) > 4 else "",
                    "strand": parts[5] if len(parts) > 5 else "",
                    "genome_id": genome_id,
                    "query_class": query_class,
                    "qseqid": qseqid,
                }
            )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str | int | float]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("/cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv"))
    parser.add_argument("--centromere-bed", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits.bed"))
    parser.add_argument("--motif", default="GCATAGGT")
    parser.add_argument("--out-hits", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits_motif_scan.tsv"))
    parser.add_argument("--out-retained-bed", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits_GCATAGGT_retained.bed"))
    parser.add_argument("--out-contig-summary", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_motif_retained_contig_summary.tsv"))
    parser.add_argument("--out-genome-summary", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_motif_retained_genome_summary.tsv"))
    args = parser.parse_args()

    motif = args.motif.upper()
    assemblies = read_manifest(args.manifest)
    bed_rows = read_bed(args.centromere_bed)

    fasta_cache: dict[str, dict[str, str]] = {}
    hit_rows = []
    retained_bed_rows = []
    skipped = []

    for row in bed_rows:
        genome_id = str(row["genome_id"])
        contig_id = str(row["contig_id"])
        if not genome_id or genome_id not in assemblies:
            skipped.append(str(row["name"]))
            continue
        if genome_id not in fasta_cache:
            fasta_cache[genome_id] = read_fasta(assemblies[genome_id])
        contigs = fasta_cache[genome_id]
        if contig_id not in contigs:
            skipped.append(str(row["name"]))
            continue

        contig_seq = contigs[contig_id]
        start = max(0, int(row["start_0based"]))
        end = min(len(contig_seq), int(row["end_0based"]))
        interval_seq = contig_seq[start:end]
        motif_positions = find_all(interval_seq, motif)
        motif_genomic_positions = [start + pos for pos in motif_positions]
        motif_count = len(motif_positions)
        retained = motif_count > 0
        hit_row = {
            "genome_id": genome_id,
            "contig_id": contig_id,
            "start_0based": start,
            "end_0based": end,
            "length_bp": end - start,
            "query_class": row["query_class"],
            "qseqid": row["qseqid"],
            "score": row["score"],
            "strand": row["strand"],
            "motif": motif,
            "motif_count": motif_count,
            "motif_retained": "yes" if retained else "no",
            "motif_positions_within_hit_0based": ",".join(str(pos) for pos in motif_positions),
            "motif_positions_genomic_0based": ",".join(str(pos) for pos in motif_genomic_positions),
        }
        hit_rows.append(hit_row)
        if retained:
            retained_bed_rows.append(row)

    write_tsv(args.out_hits, hit_rows, list(hit_rows[0].keys()) if hit_rows else [
        "genome_id",
        "contig_id",
        "start_0based",
        "end_0based",
        "length_bp",
        "query_class",
        "qseqid",
        "score",
        "strand",
        "motif",
        "motif_count",
        "motif_retained",
        "motif_positions_within_hit_0based",
        "motif_positions_genomic_0based",
    ])

    args.out_retained_bed.parent.mkdir(parents=True, exist_ok=True)
    with args.out_retained_bed.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for row in retained_bed_rows:
            writer.writerow([row["contig_id"], row["start_0based"], row["end_0based"], row["name"], row["score"], row["strand"]])

    by_contig: dict[tuple[str, str], list[dict[str, str | int | float]]] = defaultdict(list)
    by_genome: dict[str, list[dict[str, str | int | float]]] = defaultdict(list)
    for row in hit_rows:
        by_genome[str(row["genome_id"])].append(row)
        if row["motif_retained"] == "yes":
            by_contig[(str(row["genome_id"]), str(row["contig_id"]))].append(row)

    contig_summary = []
    for (genome_id, contig_id), rows in sorted(by_contig.items()):
        starts = [int(row["start_0based"]) for row in rows]
        ends = [int(row["end_0based"]) for row in rows]
        contig_summary.append(
            {
                "genome_id": genome_id,
                "contig_id": contig_id,
                "motif_retained_hit_count": len(rows),
                "total_motif_instances": sum(int(row["motif_count"]) for row in rows),
                "query_classes": ",".join(sorted({str(row["query_class"]) for row in rows})),
                "candidate_region_start": min(starts),
                "candidate_region_end": max(ends),
                "candidate_region_span_bp": max(ends) - min(starts),
                "max_score": max(float(row["score"]) for row in rows if str(row["score"]) != ""),
            }
        )

    genome_summary = []
    for genome_id, rows in sorted(by_genome.items()):
        retained_rows = [row for row in rows if row["motif_retained"] == "yes"]
        genome_summary.append(
            {
                "genome_id": genome_id,
                "blast_hit_count": len(rows),
                "motif_retained_hit_count": len(retained_rows),
                "motif_retained_contigs": len({str(row["contig_id"]) for row in retained_rows}),
                "total_motif_instances": sum(int(row["motif_count"]) for row in retained_rows),
            }
        )

    write_tsv(args.out_contig_summary, contig_summary, list(contig_summary[0].keys()) if contig_summary else [
        "genome_id",
        "contig_id",
        "motif_retained_hit_count",
        "total_motif_instances",
        "query_classes",
        "candidate_region_start",
        "candidate_region_end",
        "candidate_region_span_bp",
        "max_score",
    ])
    write_tsv(args.out_genome_summary, genome_summary, list(genome_summary[0].keys()) if genome_summary else [
        "genome_id",
        "blast_hit_count",
        "motif_retained_hit_count",
        "motif_retained_contigs",
        "total_motif_instances",
    ])

    print(f"Wrote motif scan for BLAST hits: {args.out_hits} ({len(hit_rows)} hits)")
    print(f"Wrote motif-retained BED: {args.out_retained_bed} ({len(retained_bed_rows)} hits)")
    print(f"Wrote motif-retained contig summary: {args.out_contig_summary} ({len(contig_summary)} contigs)")
    print(f"Wrote motif-retained genome summary: {args.out_genome_summary} ({len(genome_summary)} genomes)")
    if skipped:
        print(f"Skipped intervals with missing genome/contig: {len(skipped)}")


if __name__ == "__main__":
    main()
