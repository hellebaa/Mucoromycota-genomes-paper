#!/usr/bin/env python3
"""Scan assembly contig ends for canonical TTAGGG telomere-repeat arrays."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def read_motifs(path: Path | None, motif: str, reverse_motif: str) -> list[dict[str, str]]:
    if path is None:
        return [
            {
                "motif": motif.upper(),
                "reverse_motif": reverse_motif.upper() if reverse_motif else reverse_complement(motif),
                "label": motif.upper(),
                "notes": "",
            }
        ]

    motifs = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            motif_value = row["motif"].strip().upper()
            motifs.append(
                {
                    "motif": motif_value,
                    "reverse_motif": reverse_complement(motif_value),
                    "label": row.get("label", motif_value).strip() or motif_value,
                    "notes": row.get("notes", "").strip(),
                }
            )
    return motifs


def selected(row: dict[str, str], include: str) -> bool:
    if include == "all":
        return True
    column = "include_core_analyses" if include == "core" else "include_final_ultrametric_tree"
    return row.get(column, "").strip().lower() == "yes"


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


def scan_window(seq: str, motif: str, min_copies: int) -> list[tuple[int, int, str, int]]:
    pattern = re.compile(f"(?:{motif})" + "{" + f"{min_copies}," + "}")
    hits = []
    for match in pattern.finditer(seq):
        start, end = match.span()
        hits.append((start, end, motif, (end - start) // len(motif)))
    return hits


def best_terminal_hit(
    seq: str,
    side: str,
    motif: str,
    reverse_motif: str,
    window_bp: int,
    min_copies: int,
    max_distance_bp: int,
) -> dict[str, int | str] | None:
    contig_len = len(seq)
    if side == "left":
        offset = 0
        window = seq[:window_bp]
    else:
        offset = max(0, contig_len - window_bp)
        window = seq[offset:]

    hits = []
    for current_motif in (motif, reverse_motif):
        for start, end, motif_name, copies in scan_window(window, current_motif, min_copies):
            genomic_start = offset + start
            genomic_end = offset + end
            distance = genomic_start if side == "left" else contig_len - genomic_end
            if distance <= max_distance_bp:
                hits.append(
                    {
                        "start": genomic_start,
                        "end": genomic_end,
                        "motif": motif_name,
                        "copies": copies,
                        "distance_to_contig_end_bp": distance,
                    }
                )

    if not hits:
        return None
    return sorted(hits, key=lambda hit: (hit["distance_to_contig_end_bp"], -hit["copies"]))[0]


def write_tsv(path: Path, rows: list[dict[str, int | str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("/cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv"))
    parser.add_argument("--out-bed", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/telomere_hits.bed"))
    parser.add_argument("--out-summary", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/telomere_contig_summary.tsv"))
    parser.add_argument("--include", choices=["core", "final_tree", "all"], default="core")
    parser.add_argument("--motif", default="TTAGGG")
    parser.add_argument("--reverse-motif", default="CCCTAA")
    parser.add_argument("--motifs", type=Path, help="Optional TSV with columns motif, label, notes.")
    parser.add_argument("--window-bp", type=int, default=10_000)
    parser.add_argument("--max-distance-bp", type=int, default=1_000)
    parser.add_argument("--min-copies", type=int, default=3)
    args = parser.parse_args()
    motifs = read_motifs(args.motifs, args.motif, args.reverse_motif)

    bed_rows: list[dict[str, int | str]] = []
    summary_rows: list[dict[str, int | str]] = []

    for row in read_manifest(args.manifest):
        if not selected(row, args.include):
            continue
        assembly = row.get("assembly_fasta", "")
        if not assembly or not Path(assembly).exists():
            continue
        for contig_id, seq in fasta_records(Path(assembly)):
            side_hits = {"left": [], "right": []}
            for motif_info in motifs:
                for side in ("left", "right"):
                    hit = best_terminal_hit(
                        seq,
                        side,
                        motif_info["motif"],
                        motif_info["reverse_motif"],
                        args.window_bp,
                        args.min_copies,
                        args.max_distance_bp,
                    )
                    if hit is not None:
                        hit["label"] = motif_info["label"]
                        side_hits[side].append(hit)

            best_hits = {}
            for side, hits in side_hits.items():
                best_hits[side] = sorted(
                    hits,
                    key=lambda hit: (int(hit["distance_to_contig_end_bp"]), -int(hit["copies"]), -len(str(hit["motif"]))),
                )[0] if hits else None

            for side, hit in best_hits.items():
                if hit is None:
                    continue
                bed_rows.append(
                    {
                        "chrom": contig_id,
                        "start": hit["start"],
                        "end": hit["end"],
                        "name": f"{row['genome_id']}|{side}_telomere|{hit['label']}|{hit['motif']}x{hit['copies']}",
                        "score": hit["copies"],
                        "strand": "+" if hit["motif"] in {item["motif"] for item in motifs} else "-",
                    }
                )
            summary_rows.append(
                {
                    "genome_id": row["genome_id"],
                    "species_name": row.get("species_name", ""),
                    "contig_id": contig_id,
                    "contig_length_bp": len(seq),
                    "left_telomere": "yes" if best_hits["left"] else "no",
                    "left_motif_label": best_hits["left"]["label"] if best_hits["left"] else "",
                    "left_motif": best_hits["left"]["motif"] if best_hits["left"] else "",
                    "left_copies": best_hits["left"]["copies"] if best_hits["left"] else 0,
                    "left_distance_to_end_bp": best_hits["left"]["distance_to_contig_end_bp"] if best_hits["left"] else "",
                    "right_telomere": "yes" if best_hits["right"] else "no",
                    "right_motif_label": best_hits["right"]["label"] if best_hits["right"] else "",
                    "right_motif": best_hits["right"]["motif"] if best_hits["right"] else "",
                    "right_copies": best_hits["right"]["copies"] if best_hits["right"] else 0,
                    "right_distance_to_end_bp": best_hits["right"]["distance_to_contig_end_bp"] if best_hits["right"] else "",
                    "both_ends_telomeric": "yes" if best_hits["left"] and best_hits["right"] else "no",
                    "left_motif_candidates": ";".join(f"{hit['label']}:{hit['motif']}x{hit['copies']}" for hit in side_hits["left"]),
                    "right_motif_candidates": ";".join(f"{hit['label']}:{hit['motif']}x{hit['copies']}" for hit in side_hits["right"]),
                }
            )

    args.out_bed.parent.mkdir(parents=True, exist_ok=True)
    with args.out_bed.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["chrom", "start", "end", "name", "score", "strand"], delimiter="\t")
        writer.writerows(bed_rows)
    if summary_rows:
        write_tsv(args.out_summary, summary_rows, list(summary_rows[0].keys()))
    print(f"Wrote telomere BED: {args.out_bed} ({len(bed_rows)} hits)")
    print(f"Wrote telomere summary: {args.out_summary} ({len(summary_rows)} contigs)")


if __name__ == "__main__":
    main()
