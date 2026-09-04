#!/usr/bin/env python3
"""Estimate rough chromosome-scale contig counts from assembly summaries."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/assembly_contig_gene_summary.tsv"))
    parser.add_argument("--telomeres", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/telomere_contig_summary_multi_motif.tsv"))
    parser.add_argument("--out", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/rough_karyotype_summary.tsv"))
    parser.add_argument("--min-chromosome-mb", type=float, default=0.5)
    parser.add_argument("--drop-ratio", type=float, default=0.25)
    args = parser.parse_args()

    contigs = read_tsv(args.contigs)
    telomeres = read_tsv(args.telomeres) if args.telomeres.exists() else []
    tel_by_key = {(row["genome_id"], row["contig_id"]): row for row in telomeres}

    by_genome: dict[str, list[dict[str, str]]] = {}
    for row in contigs:
        by_genome.setdefault(row["genome_id"], []).append(row)

    out_rows = []
    for genome_id, rows in sorted(by_genome.items()):
        rows = sorted(rows, key=lambda row: int(row["contig_rank"]))
        chrom_like = [
            row for row in rows
            if row.get("gene_bearing_contig") == "yes" and to_float(row.get("length_mb", "")) >= args.min_chromosome_mb
        ]

        sharp_drop_after_rank = ""
        for current, next_row in zip(rows, rows[1:]):
            current_len = to_float(current.get("length_mb", ""))
            next_len = to_float(next_row.get("length_mb", ""))
            if current_len > 0 and next_len / current_len <= args.drop_ratio:
                sharp_drop_after_rank = current["contig_rank"]
                break

        tel_supported = 0
        two_telomeres = 0
        for row in chrom_like:
            tel = tel_by_key.get((genome_id, row["contig_id"]), {})
            has_left = tel.get("left_telomere") == "yes"
            has_right = tel.get("right_telomere") == "yes"
            if has_left or has_right:
                tel_supported += 1
            if has_left and has_right:
                two_telomeres += 1

        out_rows.append(
            {
                "genome_id": genome_id,
                "species_name": rows[0].get("species_name", ""),
                "n_contigs_total": len(rows),
                "assembly_mb": round(sum(to_float(row.get("length_mb", "")) for row in rows), 3),
                "chromosome_scale_contigs_min_0p5mb_gene_bearing": len(chrom_like),
                "chromosome_scale_bp": int(sum(to_float(row.get("length_bp", "")) for row in chrom_like)),
                "chromosome_scale_assembly_percent": round(
                    100 * sum(to_float(row.get("length_bp", "")) for row in chrom_like) /
                    max(1, sum(to_float(row.get("length_bp", "")) for row in rows)),
                    3,
                ),
                "chromosome_scale_contigs_with_telomere": tel_supported,
                "chromosome_scale_contigs_with_two_telomeres": two_telomeres,
                "sharp_drop_after_contig_rank": sharp_drop_after_rank,
                "rough_karyotype_note": "draft estimate; confirm with telomeres, Hi-C, centromeres, and manual inspection",
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"Wrote rough karyotype summary: {args.out} ({len(out_rows)} genomes)")


if __name__ == "__main__":
    main()
