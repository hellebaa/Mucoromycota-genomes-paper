#!/usr/bin/env python3
"""Extract FASTA sequences from BED intervals."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_assembly(manifest: Path, genome_id: str) -> Path:
    rows = read_manifest(manifest)
    for row in rows:
        if row.get("genome_id") == genome_id:
            assembly = row.get("assembly_fasta", "")
            if assembly and Path(assembly).exists():
                return Path(assembly)
            raise SystemExit(f"Assembly FASTA missing for {genome_id}: {assembly}")
    raise SystemExit(f"Genome ID not found in manifest: {genome_id}")


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


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[index : index + width] for index in range(0, len(seq), width))


def read_bed(path: Path) -> list[dict[str, str | int]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rows.append(
                {
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "name": parts[3] if len(parts) >= 4 else f"{parts[0]}:{parts[1]}-{parts[2]}",
                }
            )
    return rows


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{wrap_sequence(sequence)}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--assembly", type=Path)
    parser.add_argument("--manifest", type=Path, default=Path("/cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv"))
    parser.add_argument("--genome-id", default="MucCirc3")
    parser.add_argument("--out-all", type=Path, required=True)
    parser.add_argument("--out-cen", type=Path, required=True)
    parser.add_argument("--out-gremline", type=Path, required=True)
    parser.add_argument("--out-index", type=Path, required=True)
    args = parser.parse_args()

    assembly = args.assembly or find_assembly(args.manifest, args.genome_id)
    sequences = read_fasta(assembly)
    bed_rows = read_bed(args.bed)

    all_records = []
    cen_records = []
    gremline_records = []
    index_rows = []
    for row in bed_rows:
        chrom = str(row["chrom"])
        start = int(row["start"])
        end = int(row["end"])
        name = str(row["name"])
        if chrom not in sequences:
            raise SystemExit(f"BED chromosome {chrom} not found in {assembly}")
        sequence = sequences[chrom][start:end]
        header = f"{name}::{chrom}:{start}-{end}"
        all_records.append((header, sequence))
        if name.startswith("CEN_"):
            cen_records.append((header, sequence))
            feature_type = "CEN_core"
        elif "GremLINE" in name:
            gremline_records.append((header, sequence))
            feature_type = "GremLINE"
        else:
            feature_type = "other"
        index_rows.append(
            {
                "query_id": header,
                "feature_name": name,
                "feature_type": feature_type,
                "source_genome": args.genome_id,
                "source_contig": chrom,
                "start_0based": start,
                "end_0based": end,
                "length_bp": len(sequence),
                "assembly_fasta": str(assembly),
            }
        )

    write_fasta(args.out_all, all_records)
    write_fasta(args.out_cen, cen_records)
    write_fasta(args.out_gremline, gremline_records)
    args.out_index.parent.mkdir(parents=True, exist_ok=True)
    with args.out_index.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(index_rows)

    print(f"Wrote all interval FASTA: {args.out_all} ({len(all_records)} records)")
    print(f"Wrote CEN core FASTA: {args.out_cen} ({len(cen_records)} records)")
    print(f"Wrote GremLINE FASTA: {args.out_gremline} ({len(gremline_records)} records)")
    print(f"Wrote query index: {args.out_index}")


if __name__ == "__main__":
    main()
