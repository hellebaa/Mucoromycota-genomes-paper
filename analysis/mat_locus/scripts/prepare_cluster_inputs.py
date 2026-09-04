#!/usr/bin/env python3
"""Organize MAT-locus inputs into standard per-analysis paths."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--mode", choices=("symlink", "copy"), default="symlink")
    return parser.parse_args()


def link_or_copy(source: Path, dest: Path, mode: str) -> None:
    if not source:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if mode == "copy":
        shutil.copy2(source, dest)
    else:
        os.symlink(source, dest)


def main() -> None:
    args = parse_args()
    inputs = args.work_dir / "inputs"
    assemblies = inputs / "assemblies"
    annotations = inputs / "annotations"
    proteins = inputs / "proteins"
    queries = inputs / "queries"
    for directory in (assemblies, annotations, proteins, queries):
        directory.mkdir(parents=True, exist_ok=True)

    link_or_copy(args.query_fasta, queries / args.query_fasta.name, args.mode)

    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    output_fields = [
        "genome_id",
        "tree_label",
        "species_name",
        "assembly_fasta",
        "annotation_gff",
        "protein_fasta",
        "ready_for_scan",
        "notes",
    ]
    organized_rows = []

    for row in rows:
        genome_id = row["genome_id"]
        notes = []
        assembly_dest = assemblies / f"{genome_id}.fasta"
        gff_dest = annotations / f"{genome_id}.gff"
        protein_dest = proteins / f"{genome_id}.faa"

        if row.get("assembly_fasta"):
            link_or_copy(Path(row["assembly_fasta"]), assembly_dest, args.mode)
        else:
            notes.append("missing assembly")
        if row.get("annotation_gff"):
            link_or_copy(Path(row["annotation_gff"]), gff_dest, args.mode)
        else:
            notes.append("missing annotation")
        if row.get("protein_fasta"):
            link_or_copy(Path(row["protein_fasta"]), protein_dest, args.mode)
        else:
            notes.append("missing proteins")

        organized_rows.append(
            {
                "genome_id": genome_id,
                "tree_label": row.get("tree_label", ""),
                "species_name": row.get("species_name", ""),
                "assembly_fasta": str(assembly_dest) if assembly_dest.exists() else "",
                "annotation_gff": str(gff_dest) if gff_dest.exists() else "",
                "protein_fasta": str(protein_dest) if protein_dest.exists() else "",
                "ready_for_scan": "yes" if not notes else "no",
                "notes": "; ".join(notes),
            }
        )

    output = inputs / "organized_samples.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(organized_rows)

    ready = sum(1 for row in organized_rows if row["ready_for_scan"] == "yes")
    print(f"Wrote {output}")
    print(f"Rows: {len(organized_rows)}; ready: {ready}; incomplete: {len(organized_rows) - ready}")


if __name__ == "__main__":
    main()
