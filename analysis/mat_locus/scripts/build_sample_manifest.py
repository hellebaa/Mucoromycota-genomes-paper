#!/usr/bin/env python3
"""Build a MAT-locus sample manifest from metadata and a data/genomes mirror."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ASSEMBLY_SUFFIXES = (
    ".contigs.fasta",
    ".decon.masked.sorted.fasta",
    ".fasta",
    "_genomic.fna",
    ".fna",
    ".fa",
)
GFF_SUFFIXES = (".gff3", ".gff")
PROTEIN_SUFFIXES = (".proteins.faa", ".faa", ".proteins.fa", ".pep")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--genomes-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--include-excluded",
        action="store_true",
        help="Include rows excluded from final core analyses.",
    )
    return parser.parse_args()


def is_candidate(path: Path, suffixes: tuple[str, ...]) -> bool:
    name = path.name
    if name.startswith("."):
        return False
    if "cds" in name.lower() or "transcript" in name.lower():
        return False
    return any(name.endswith(suffix) for suffix in suffixes)


def first_existing(paths: list[Path]) -> str:
    existing = [p for p in paths if p.exists()]
    return str(existing[0]) if existing else ""


def find_file(genome_dir: Path, genome_id: str, suffixes: tuple[str, ...]) -> str:
    preferred = [genome_dir / f"{genome_id}{suffix}" for suffix in suffixes]
    hit = first_existing(preferred)
    if hit:
        return hit
    candidates = sorted(
        p for p in genome_dir.iterdir() if p.is_file() and is_candidate(p, suffixes)
    )
    return str(candidates[0]) if candidates else ""


def find_proteins(genome_dir: Path, genome_id: str) -> str:
    preferred = [genome_dir / f"{genome_id}{suffix}" for suffix in PROTEIN_SUFFIXES]
    hit = first_existing(preferred)
    if hit:
        return hit
    candidates = sorted(
        p
        for p in genome_dir.iterdir()
        if p.is_file() and any(p.name.endswith(suffix) for suffix in PROTEIN_SUFFIXES)
    )
    return str(candidates[0]) if candidates else ""


def genome_dir_candidates(genomes_dir: Path, row: dict[str, str]) -> list[Path]:
    genome_id = row.get("genome_id", "")
    tree_label = row.get("tree_label", "")
    candidates = [genomes_dir / genome_id, genomes_dir / tree_label]
    if genome_id == "EntParv1":
        candidates.append(genomes_dir / "EntParv")
    return list(dict.fromkeys(candidates))


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.metadata.open(newline="") as handle:
        reader = csv.DictReader(handle)
        metadata_rows = list(reader)

    output_fields = [
        "genome_id",
        "tree_label",
        "species_name",
        "family",
        "order",
        "source",
        "metadata_status",
        "genome_dir",
        "assembly_fasta",
        "annotation_gff",
        "protein_fasta",
        "input_status",
        "notes",
    ]

    rows = []
    for row in metadata_rows:
        if not args.include_excluded:
            if row.get("include_core_analyses") != "yes":
                continue
            if row.get("include_final_ultrametric_tree") != "yes":
                continue

        genome_id = row["genome_id"]
        genome_dir = first_existing(genome_dir_candidates(args.genomes_dir, row))
        notes = []
        assembly = annotation = proteins = ""
        if genome_dir:
            gdir = Path(genome_dir)
            assembly = find_file(gdir, genome_id, ASSEMBLY_SUFFIXES)
            annotation = find_file(gdir, genome_id, GFF_SUFFIXES)
            proteins = find_proteins(gdir, genome_id)
        else:
            notes.append("genome directory not found")

        if not assembly:
            notes.append("assembly FASTA not found")
        if not annotation:
            notes.append("annotation GFF not found")
        if not proteins:
            notes.append("protein FASTA not found")

        rows.append(
            {
                "genome_id": genome_id,
                "tree_label": row.get("tree_label", ""),
                "species_name": row.get("species_name", ""),
                "family": row.get("family", ""),
                "order": row.get("order", ""),
                "source": row.get("source", ""),
                "metadata_status": row.get("status_source_metadata", ""),
                "genome_dir": genome_dir,
                "assembly_fasta": assembly,
                "annotation_gff": annotation,
                "protein_fasta": proteins,
                "input_status": "ready" if assembly and annotation and proteins else "incomplete",
                "notes": "; ".join(notes),
            }
        )

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    ready = sum(1 for row in rows if row["input_status"] == "ready")
    print(f"Wrote {args.output}")
    print(f"Rows: {len(rows)}; ready: {ready}; incomplete: {len(rows) - ready}")


if __name__ == "__main__":
    main()
