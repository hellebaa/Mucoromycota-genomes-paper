#!/usr/bin/env python3
"""Build a manifest of canonical per-genome files in a NeLS/Saga data mirror."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ALIASED_FOLDERS = {
    "EntParv1": "EntParv",
}

IGNORED_MISSING_GENOMES = {
    "ggAmyRoux1",
    "gzDisGlob1",
    "gzMucFlav1",
    "gzMucRace1",
}


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def rel_or_empty(path: Path | None) -> str:
    return str(path) if path else ""


def choose_file(candidates: list[Path], prefer: list[str] | None = None) -> Path | None:
    if not candidates:
        return None
    prefer = prefer or []
    candidates = sorted(candidates, key=lambda path: path.name)
    for token in prefer:
        for candidate in candidates:
            if token.lower() in candidate.name.lower():
                return candidate
    return candidates[0]


def top_level_files(folder: Path) -> list[Path]:
    if not folder.exists():
        return []
    return [path for path in folder.iterdir() if path.is_file()]


def repeat_file(folder: Path, repeat_subdir: str, suffixes: tuple[str, ...]) -> Path | None:
    repeat_dir = folder / "repeats" / repeat_subdir
    if not repeat_dir.exists():
        return None
    candidates = [
        path
        for path in repeat_dir.iterdir()
        if path.is_file() and path.name.lower().endswith(suffixes)
    ]
    return choose_file(candidates)


def build_row(data_root: Path, metadata_row: dict[str, str]) -> dict[str, str]:
    genome_id = metadata_row["genome_id"]
    folder_name = genome_id
    folder = data_root / "genomes" / folder_name
    alias_note = ""
    if not folder.exists() and genome_id in ALIASED_FOLDERS:
        folder_name = ALIASED_FOLDERS[genome_id]
        folder = data_root / "genomes" / folder_name
        alias_note = f"folder alias: {genome_id}->{folder_name}"

    files = top_level_files(folder)
    lower = {path: path.name.lower() for path in files}

    gff = choose_file([path for path, name in lower.items() if name.endswith((".gff", ".gff3"))])
    proteins = choose_file(
        [
            path
            for path, name in lower.items()
            if name.endswith((".faa", ".faa.gz"))
            or "protein.faa" in name
            or "proteins.faa" in name
            or name.endswith((".proteins.fa", ".proteins.fasta"))
        ],
        prefer=[genome_id, folder_name, "protein"],
    )
    cds = choose_file(
        [
            path
            for path, name in lower.items()
            if "cds" in name and name.endswith((".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz"))
        ],
        prefer=[genome_id, folder_name, "cds-transcripts", "cds_from_genomic"],
    )
    assembly = choose_file(
        [
            path
            for path, name in lower.items()
            if name.endswith((".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz"))
            and "cds" not in name
            and "protein" not in name
            and "proteins" not in name
        ],
        prefer=[genome_id, folder_name, "contigs", "genomic", "dna"],
    )
    idtbl = choose_file([path for path, name in lower.items() if name.endswith("_idtbl.tsv")])
    genepostbl = choose_file([path for path, name in lower.items() if name.endswith("_genepostbl.tsv")])

    denovo_library = repeat_file(folder, "denovo", (".repeatmodeler-library.fasta",))
    denovo_out = repeat_file(folder, "denovo", (".repeatmasker-denovo.out",))
    denovo_gff = repeat_file(folder, "denovo", (".repeatmasker-denovo.gff3", ".repeatmasker-denovo.out.gff"))
    denovo_tbl = repeat_file(folder, "denovo", (".repeatmasker-denovo.tbl",))

    shared_out = repeat_file(folder, "shared_library", (".sharedlib.repeatmasker.out",))
    shared_gff = repeat_file(
        folder,
        "shared_library",
        (".sharedlib.repeatmasker.out.gff", ".sharedlib.repeatmasker.gff3"),
    )
    shared_tbl = repeat_file(folder, "shared_library", (".sharedlib.repeatmasker.tbl",))
    shared_masked = repeat_file(
        folder,
        "shared_library",
        (".sharedlib.repeatmasker.masked.fasta", ".sharedlib.repeatmasker.masked", ".masked"),
    )

    required = {
        "assembly_fasta": assembly,
        "gff": gff,
        "proteins_faa": proteins,
        "cds_fasta": cds,
        "idtbl": idtbl,
        "genepostbl": genepostbl,
    }
    missing = [name for name, value in required.items() if value is None]
    status = "ok"
    if not folder.exists():
        status = "missing_folder"
    elif missing:
        status = "missing_core_files"
    if genome_id in IGNORED_MISSING_GENOMES and status == "missing_folder":
        status = "ignored_excluded_missing_folder"
    if alias_note and status != "ok":
        status = f"{status};alias_folder"

    shared_te_library = data_root / "repeats" / "shared_library" / "library" / "mucoro_TE_shared.fa"
    shared_te_note = "" if shared_te_library.exists() else "shared TE library not found at expected path"

    notes = "; ".join(
        value
        for value in [
            alias_note,
            "excluded messy/subgenome case; ignored for canonical core mirror"
            if genome_id in IGNORED_MISSING_GENOMES
            else "",
            shared_te_note,
            metadata_row.get("notes", ""),
        ]
        if value
    )

    return {
        "genome_id": genome_id,
        "species_name": metadata_row.get("species_name", ""),
        "source": metadata_row.get("source", ""),
        "include_final_ultrametric_tree": metadata_row.get("include_final_ultrametric_tree", ""),
        "include_core_analyses": metadata_row.get("include_core_analyses", ""),
        "data_root": str(data_root),
        "genome_folder": str(folder) if folder.exists() else "",
        "status": status,
        "missing_core_files": ";".join(missing),
        "assembly_fasta": rel_or_empty(assembly),
        "gff": rel_or_empty(gff),
        "proteins_faa": rel_or_empty(proteins),
        "cds_fasta": rel_or_empty(cds),
        "idtbl": rel_or_empty(idtbl),
        "genepostbl": rel_or_empty(genepostbl),
        "denovo_repeatmodeler_library": rel_or_empty(denovo_library),
        "denovo_repeatmasker_out": rel_or_empty(denovo_out),
        "denovo_repeatmasker_gff": rel_or_empty(denovo_gff),
        "denovo_repeatmasker_tbl": rel_or_empty(denovo_tbl),
        "shared_repeatmasker_out": rel_or_empty(shared_out),
        "shared_repeatmasker_gff": rel_or_empty(shared_gff),
        "shared_repeatmasker_tbl": rel_or_empty(shared_tbl),
        "shared_repeatmasker_masked_fasta": rel_or_empty(shared_masked),
        "shared_te_library_fasta": str(shared_te_library) if shared_te_library.exists() else "",
        "notes": notes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        default="/cluster/work/users/hellb/Mucor/data",
        type=Path,
        help="Local mirror of Projects/NMBU_threeD_Sandve_2024/Mucor/data.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="Path to genome_metadata_master.csv. Defaults to <data-root>/metadata/genome_metadata_master.csv.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Output manifest TSV. Defaults to <data-root>/metadata/genome_file_manifest.tsv.",
    )
    args = parser.parse_args()

    metadata = args.metadata or args.data_root / "metadata" / "genome_metadata_master.csv"
    out = args.out or args.data_root / "metadata" / "genome_file_manifest.tsv"

    metadata_rows = read_metadata(metadata)
    manifest_rows = [build_row(args.data_root, row) for row in metadata_rows]

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    issues = [
        row
        for row in manifest_rows
        if row["status"] != "ok" and not row["status"].startswith("ignored_")
    ]
    core_or_final_issues = [
        row
        for row in issues
        if row["include_core_analyses"].lower() == "yes"
        or row["include_final_ultrametric_tree"].lower() == "yes"
    ]
    print(f"Wrote: {out}")
    print(f"Genomes: {len(manifest_rows)}")
    print(f"Issues: {len(issues)}")
    print(f"Core/final issues: {len(core_or_final_issues)}")
    for row in issues:
        print(
            "\t".join(
                [
                    row["genome_id"],
                    row["status"],
                    row["missing_core_files"],
                    row["genome_folder"],
                ]
            )
        )


if __name__ == "__main__":
    main()
