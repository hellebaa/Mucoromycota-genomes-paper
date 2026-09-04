#!/usr/bin/env python3
"""Summarize proximity between CEN-core BLAST hits and GremLINE BLAST hits."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_name(name: str) -> tuple[str, str, str]:
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
            genome_id, evidence_class, query_id = parse_name(parts[3])
            rows.append(
                {
                    "genome_id": genome_id,
                    "contig_id": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "name": parts[3],
                    "score": float(parts[4]) if len(parts) > 4 and parts[4] else 0.0,
                    "strand": parts[5] if len(parts) > 5 else "",
                    "evidence_class": evidence_class,
                    "query_id": query_id,
                }
            )
    return rows


def interval_distance(a: dict[str, str | int | float], b: dict[str, str | int | float]) -> int:
    a_start = int(a["start"])
    a_end = int(a["end"])
    b_start = int(b["start"])
    b_end = int(b["end"])
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0


def midpoint(row: dict[str, str | int | float]) -> float:
    return (int(row["start"]) + int(row["end"])) / 2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cen-bed", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits.bed"))
    parser.add_argument("--gremline-bed", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/gremline_blast_hits.bed"))
    parser.add_argument("--max-distance-bp", type=int, default=100_000)
    parser.add_argument("--out-cen-summary", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/cen_gremline_proximity_cen_hit_summary.tsv"))
    parser.add_argument("--out-contig-summary", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/cen_gremline_proximity_contig_summary.tsv"))
    args = parser.parse_args()

    cen_hits = read_bed(args.cen_bed)
    grem_hits = read_bed(args.gremline_bed)

    grem_by_contig: dict[tuple[str, str], list[dict[str, str | int | float]]] = defaultdict(list)
    for hit in grem_hits:
        grem_by_contig[(str(hit["genome_id"]), str(hit["contig_id"]))].append(hit)

    cen_summary = []
    for cen in cen_hits:
        key = (str(cen["genome_id"]), str(cen["contig_id"]))
        candidates = grem_by_contig.get(key, [])
        nearby = []
        nearest = None
        for grem in candidates:
            distance = interval_distance(cen, grem)
            if nearest is None or distance < int(nearest["distance_bp"]):
                nearest = {
                    "distance_bp": distance,
                    "gremline_start": grem["start"],
                    "gremline_end": grem["end"],
                    "gremline_name": grem["name"],
                    "gremline_strand": grem["strand"],
                    "midpoint_delta_bp": int(midpoint(grem) - midpoint(cen)),
                }
            if distance <= args.max_distance_bp:
                nearby.append(grem)
        cen_summary.append(
            {
                "genome_id": cen["genome_id"],
                "contig_id": cen["contig_id"],
                "cen_start": cen["start"],
                "cen_end": cen["end"],
                "cen_name": cen["name"],
                "cen_score": cen["score"],
                "gremline_hits_same_contig": len(candidates),
                "gremline_hits_within_threshold": len(nearby),
                "has_nearby_gremline": "yes" if nearby else "no",
                "nearest_gremline_distance_bp": nearest["distance_bp"] if nearest else "",
                "nearest_gremline_midpoint_delta_bp": nearest["midpoint_delta_bp"] if nearest else "",
                "nearest_gremline_name": nearest["gremline_name"] if nearest else "",
                "nearest_gremline_strand": nearest["gremline_strand"] if nearest else "",
            }
        )

    grouped: dict[tuple[str, str], list[dict[str, str | int | float]]] = defaultdict(list)
    for row in cen_summary:
        grouped[(str(row["genome_id"]), str(row["contig_id"]))].append(row)

    contig_summary = []
    for (genome_id, contig_id), rows in sorted(grouped.items()):
        nearest_distances = [
            int(row["nearest_gremline_distance_bp"])
            for row in rows
            if str(row["nearest_gremline_distance_bp"]) != ""
        ]
        contig_summary.append(
            {
                "genome_id": genome_id,
                "contig_id": contig_id,
                "cen_hit_count": len(rows),
                "cen_hits_with_nearby_gremline": sum(1 for row in rows if row["has_nearby_gremline"] == "yes"),
                "has_cen_and_gremline_context": "yes" if any(row["has_nearby_gremline"] == "yes" for row in rows) else "no",
                "gremline_hits_same_contig_max": max(int(row["gremline_hits_same_contig"]) for row in rows),
                "nearest_gremline_distance_bp_min": min(nearest_distances) if nearest_distances else "",
            }
        )

    args.out_cen_summary.parent.mkdir(parents=True, exist_ok=True)
    with args.out_cen_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cen_summary[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(cen_summary)

    with args.out_contig_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(contig_summary[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(contig_summary)

    print(f"Wrote CEN-hit proximity summary: {args.out_cen_summary} ({len(cen_summary)} CEN hits)")
    print(f"Wrote contig proximity summary: {args.out_contig_summary} ({len(contig_summary)} contigs)")


if __name__ == "__main__":
    main()
