#!/usr/bin/env python3
"""Summarize candidate centromere/GremLINE BLAST hits by genome and contig."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
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
]


def query_class(qseqid: str) -> str:
    if qseqid.startswith("CEN_"):
        return "CEN_motif"
    if "GremLINE" in qseqid:
        return "GremLINE"
    return "other_centromere_query"


def genome_from_file(path: Path, suffix: str) -> str:
    name = path.name
    if suffix and name.endswith(suffix):
        return name[: -len(suffix)]
    blastout_match = re.search(r"blastout_CEN_([^_]+)(?:_table)?$", name)
    if blastout_match:
        return blastout_match.group(1)
    for token in (".blast.tsv", ".blast.out", ".tsv", ".out", "_table"):
        if name.endswith(token):
            return name[: -len(token)]
    match = re.search(r"(gz|gg)?[A-Za-z]+[A-Za-z0-9]+", name)
    return match.group(0) if match else path.stem


def parse_blast(path: Path, genome_id: str, min_pident: float, min_length: int, max_evalue: float):
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            row = dict(zip(BLAST_FIELDS, parts[:12]))
            try:
                row["pident"] = float(row["pident"])
                row["length"] = int(row["length"])
                row["sstart"] = int(row["sstart"])
                row["send"] = int(row["send"])
                row["evalue"] = float(row["evalue"])
                row["bitscore"] = float(row["bitscore"])
            except ValueError:
                continue
            if row["pident"] < min_pident or row["length"] < min_length or row["evalue"] > max_evalue:
                continue
            start = min(row["sstart"], row["send"]) - 1
            end = max(row["sstart"], row["send"])
            row.update(
                {
                    "genome_id": genome_id,
                    "query_class": query_class(str(row["qseqid"])),
                    "bed_start": start,
                    "bed_end": end,
                    "strand": "+" if row["sstart"] <= row["send"] else "-",
                }
            )
            yield row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast-dir", type=Path, required=True)
    parser.add_argument("--pattern", default="*")
    parser.add_argument("--file-suffix", default="")
    parser.add_argument("--out-bed", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits.bed"))
    parser.add_argument("--out-summary", type=Path, default=Path("/cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_contig_summary.tsv"))
    parser.add_argument("--min-pident", type=float, default=80.0)
    parser.add_argument("--min-length", type=int, default=80)
    parser.add_argument("--max-evalue", type=float, default=1e-10)
    args = parser.parse_args()

    hits = []
    for path in sorted(args.blast_dir.glob(args.pattern)):
        if not path.is_file():
            continue
        genome_id = genome_from_file(path, args.file_suffix)
        hits.extend(parse_blast(path, genome_id, args.min_pident, args.min_length, args.max_evalue))

    args.out_bed.parent.mkdir(parents=True, exist_ok=True)
    with args.out_bed.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for hit in hits:
            writer.writerow(
                [
                    hit["sseqid"],
                    hit["bed_start"],
                    hit["bed_end"],
                    f"{hit['genome_id']}|{hit['query_class']}|{hit['qseqid']}",
                    hit["bitscore"],
                    hit["strand"],
                ]
            )

    grouped: dict[tuple[str, str], list[dict[str, str | int | float]]] = defaultdict(list)
    for hit in hits:
        grouped[(str(hit["genome_id"]), str(hit["sseqid"]))].append(hit)

    summary_rows = []
    for (genome_id, contig_id), contig_hits in sorted(grouped.items()):
        classes = sorted({str(hit["query_class"]) for hit in contig_hits})
        starts = [int(hit["bed_start"]) for hit in contig_hits]
        ends = [int(hit["bed_end"]) for hit in contig_hits]
        summary_rows.append(
            {
                "genome_id": genome_id,
                "contig_id": contig_id,
                "n_centromere_query_hits": len(contig_hits),
                "query_classes": ",".join(classes),
                "has_cen_motif_hit": "yes" if "CEN_motif" in classes else "no",
                "has_gremline_hit": "yes" if "GremLINE" in classes else "no",
                "candidate_region_start": min(starts),
                "candidate_region_end": max(ends),
                "candidate_region_span_bp": max(ends) - min(starts),
                "max_bitscore": max(float(hit["bitscore"]) for hit in contig_hits),
                "max_pident": max(float(hit["pident"]) for hit in contig_hits),
            }
        )

    with args.out_summary.open("w", newline="") as handle:
        if summary_rows:
            writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)

    print(f"Wrote centromere BLAST BED: {args.out_bed} ({len(hits)} hits)")
    print(f"Wrote centromere contig summary: {args.out_summary} ({len(summary_rows)} contigs)")


if __name__ == "__main__":
    main()
