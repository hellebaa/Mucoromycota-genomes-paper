#!/usr/bin/env python3
"""Audit the NeLS data/genomes folder against curated genome metadata."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


ALIASED_REMOTE_FOLDERS = {
    "EntParv1": "EntParv",
}


def classify_file(filename: str) -> str:
    name = filename.lower()
    if name.endswith((".gff", ".gff3")):
        return "gff"
    if name.endswith((".faa", ".faa.gz")) or "protein.faa" in name or "proteins.faa" in name:
        return "proteins"
    if "cds" in name and name.endswith((".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz")):
        return "cds"
    if name.endswith((".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz")) and (
        "genomic" in name
        or "dna" in name
        or "contig" in name
        or "fasta" in name
        or "masked" in name
    ):
        return "assembly"
    if name.endswith("_idtbl.tsv"):
        return "idtbl"
    if name.endswith("_genepostbl.tsv"):
        return "genepostbl"
    return "other"


def read_manifest(path: Path) -> dict[str, list[tuple[str, int]]]:
    files_by_folder: dict[str, list[tuple[str, int]]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            relpath = row["relative_path"]
            if not relpath:
                continue
            if "/" not in relpath and row["type"] == "file":
                continue
            folder = relpath.split("/", 1)[0]
            files_by_folder.setdefault(folder, [])
            if row["type"] == "file" and "/" in relpath:
                filename = relpath.split("/", 1)[1]
                files_by_folder[folder].append((filename, int(row["size_bytes"] or 0)))
    return files_by_folder


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--unmatched-out", required=True, type=Path)
    args = parser.parse_args()

    files_by_folder = read_manifest(args.manifest)
    metadata_rows = read_metadata(args.metadata)
    metadata_ids = [row["genome_id"] for row in metadata_rows]
    core_ids = {
        row["genome_id"]
        for row in metadata_rows
        if row.get("include_core_analyses", "").lower() == "yes"
    }
    final_tree_ids = {
        row["genome_id"]
        for row in metadata_rows
        if row.get("include_final_ultrametric_tree", "").lower() == "yes"
    }

    audit_rows = []
    for row in metadata_rows:
        genome_id = row["genome_id"]
        remote_folder = genome_id if genome_id in files_by_folder else ALIASED_REMOTE_FOLDERS.get(genome_id, genome_id)
        present = remote_folder in files_by_folder
        files = files_by_folder.get(remote_folder, [])

        files_by_type = {
            "assembly": [],
            "gff": [],
            "proteins": [],
            "cds": [],
            "idtbl": [],
            "genepostbl": [],
            "other": [],
        }
        zero_byte_files = []
        for filename, size in files:
            files_by_type[classify_file(filename)].append(filename)
            if size == 0:
                zero_byte_files.append(filename)

        missing = []
        required_classes = [
            ("assembly", "assembly_fasta"),
            ("gff", "gff_or_gff3"),
            ("proteins", "proteins_faa"),
            ("cds", "cds_or_transcripts"),
            ("idtbl", "idtbl"),
            ("genepostbl", "genepostbl"),
        ]
        for key, label in required_classes:
            if not files_by_type[key]:
                missing.append(label)

        status = "ok"
        if not present:
            status = "missing_folder"
        elif missing:
            status = "missing_files"
        if genome_id != remote_folder:
            status = f"{status};alias_folder"
        if zero_byte_files:
            status = f"{status};zero_byte_files"

        audit_rows.append(
            {
                "genome_id": genome_id,
                "species_name": row.get("species_name", ""),
                "source": row.get("source", ""),
                "nels_folder": remote_folder if present else "",
                "in_final_ultrametric_tree": "yes" if genome_id in final_tree_ids else "no",
                "in_core_analyses": "yes" if genome_id in core_ids else "no",
                "status": status,
                "missing_file_classes": "; ".join(missing),
                "zero_byte_files": "; ".join(zero_byte_files),
                "assembly_files": "; ".join(files_by_type["assembly"]),
                "gff_files": "; ".join(files_by_type["gff"]),
                "protein_files": "; ".join(files_by_type["proteins"]),
                "cds_files": "; ".join(files_by_type["cds"]),
                "idtbl_files": "; ".join(files_by_type["idtbl"]),
                "genepostbl_files": "; ".join(files_by_type["genepostbl"]),
                "other_files": "; ".join(files_by_type["other"]),
            }
        )

    consumed_remote_folders = set(metadata_ids) | set(ALIASED_REMOTE_FOLDERS.values())
    unmatched_remote_folders = sorted(set(files_by_folder) - consumed_remote_folders)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)

    with args.unmatched_out.open("w") as handle:
        for folder in unmatched_remote_folders:
            handle.write(f"{folder}\n")

    print(f"metadata_rows\t{len(metadata_rows)}")
    print(f"remote_top_level_folders\t{len(files_by_folder)}")
    print(f"final_ultrametric_expected\t{len(final_tree_ids)}")
    print(f"core_expected\t{len(core_ids)}")
    print(f"status_counts\t{dict(Counter(row['status'] for row in audit_rows))}")
    print("core_or_final_with_issues")
    for row in audit_rows:
        if (
            row["in_core_analyses"] == "yes" or row["in_final_ultrametric_tree"] == "yes"
        ) and row["status"] != "ok":
            print(
                "\t".join(
                    [
                        row["genome_id"],
                        row["nels_folder"],
                        row["status"],
                        row["missing_file_classes"],
                    ]
                )
            )
    print(f"unmatched_remote_folders\t{len(unmatched_remote_folders)}")


if __name__ == "__main__":
    main()
