#!/usr/bin/env python3
"""Map MAT BLAST hits back to GFF gene neighborhoods.

The protein homology screen is intentionally broad. This script asks the more
useful MAT-locus question: does the best sexP/sexM-like HMG hit map to a gene
that is locally near the best tptA and rnhA hits?
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def parse_attributes(value: str) -> dict[str, str]:
    attrs = {}
    for part in value.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, val = part.split("=", 1)
        elif " " in part:
            key, val = part.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = val.strip().strip('"')
    return attrs


def split_values(value: str) -> list[str]:
    parts = []
    for item in re.split(r"[,| ]+", value or ""):
        item = item.strip()
        if item:
            parts.append(item)
    return parts


def normalized_ids(raw_id: str) -> set[str]:
    if not raw_id:
        return set()
    ids = {raw_id, raw_id.split()[0]}
    for item in list(ids):
        ids.add(re.sub(r"-T\d+$", "", item))
        ids.add(re.sub(r"\.t\d+$", "", item))
        ids.add(re.sub(r"\.p\d+$", "", item))
        ids.add(item.replace("transcript:", ""))
        ids.add(item.replace("protein:", ""))
    return {item for item in ids if item}


def subject_candidates(row: dict[str, str]) -> set[str]:
    candidates = set()
    for field in ("sseqid", "stitle"):
        value = row.get(field, "")
        candidates.update(normalized_ids(value))
        for token in split_values(value.replace(";", " ").replace("=", " ")):
            candidates.update(normalized_ids(token))
    return candidates


def parse_gff(path: Path) -> tuple[list[dict[str, object]], dict[str, int]]:
    genes: list[dict[str, object]] = []
    gene_index_by_id: dict[str, int] = {}
    child_to_gene: dict[str, str] = {}
    pending_aliases: dict[str, list[str]] = {}
    pending_text: dict[str, list[str]] = {}

    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attrs_raw = fields
            attrs = parse_attributes(attrs_raw)
            feature_id = attrs.get("ID", "")
            parents = [item for item in attrs.get("Parent", "").split(",") if item]

            if feature_type == "gene":
                aliases = []
                for key in (
                    "ID",
                    "Name",
                    "Alias",
                    "locus_tag",
                    "protein_id",
                    "orig_protein_id",
                    "transcriptId",
                    "proteinId",
                ):
                    aliases.extend(split_values(attrs.get(key, "")))
                gene_id = feature_id or attrs.get("Name", "") or attrs.get("locus_tag", "")
                gene = {
                    "seqid": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene_id": gene_id,
                    "name": attrs.get("Name", ""),
                    "locus_tag": attrs.get("locus_tag", ""),
                    "aliases": aliases + pending_aliases.pop(gene_id, []),
                    "annotation_text_parts": [
                        value
                        for value in [
                            attrs.get("product", ""),
                            attrs.get("product_name", ""),
                            attrs.get("Note", ""),
                            attrs.get("description", ""),
                            attrs_raw,
                            *pending_text.pop(gene_id, []),
                        ]
                        if value
                    ],
                }
                genes.append(gene)
                gene_index_by_id[gene_id] = len(genes) - 1
                continue

            parent_gene_id = ""
            if feature_type in {"mRNA", "transcript"} and parents:
                parent_gene_id = parents[0]
                if feature_id:
                    child_to_gene[feature_id] = parent_gene_id
            elif feature_type == "CDS" and parents:
                parent_gene_id = child_to_gene.get(parents[0], parents[0])

            if not parent_gene_id:
                continue

            aliases = []
            for key in (
                "ID",
                "Name",
                "Alias",
                "locus_tag",
                "protein_id",
                "orig_protein_id",
                "transcriptId",
                "proteinId",
            ):
                aliases.extend(split_values(attrs.get(key, "")))
            aliases.append(feature_id)

            if parent_gene_id in gene_index_by_id:
                gene = genes[gene_index_by_id[parent_gene_id]]
                gene["aliases"] = list(gene["aliases"]) + aliases
                text = " | ".join(
                    value
                    for value in [
                        attrs.get("product", ""),
                        attrs.get("product_name", ""),
                        attrs.get("Note", ""),
                        attrs.get("description", ""),
                        attrs_raw,
                    ]
                    if value
                )
                if text:
                    gene["annotation_text_parts"] = list(gene["annotation_text_parts"]) + [text]
            else:
                pending_aliases.setdefault(parent_gene_id, []).extend(aliases)
                pending_text.setdefault(parent_gene_id, []).append(attrs_raw)

    genes.sort(key=lambda row: (str(row["seqid"]), int(row["start"]), int(row["end"])))
    id_to_index: dict[str, int] = {}
    for index, gene in enumerate(genes):
        gene["order_index"] = index
        deduped_text = []
        for text in gene.get("annotation_text_parts", []):
            if text and text not in deduped_text:
                deduped_text.append(str(text))
        gene["annotation_text"] = " | ".join(deduped_text)
        all_aliases = [str(gene.get("gene_id", "")), str(gene.get("name", "")), str(gene.get("locus_tag", ""))]
        all_aliases.extend(str(alias) for alias in gene.get("aliases", []))
        for alias in all_aliases:
            for norm in normalized_ids(alias):
                id_to_index.setdefault(norm, index)
    return genes, id_to_index


def load_gff_paths(annotation_samples: Path) -> dict[str, Path]:
    paths = {}
    for row in read_tsv(annotation_samples):
        gff = Path(row.get("gff_path", ""))
        if gff.exists():
            paths[row["sample_id"]] = gff
    return paths


def best_hits_by_genome_gene(best_per_gene: Path) -> dict[str, dict[str, dict[str, str]]]:
    hits: dict[str, dict[str, dict[str, str]]] = {}
    for row in read_tsv(best_per_gene):
        gene = row["gene_symbol"]
        genome = row["genome_id"]
        if gene not in {"tptA", "rnhA", "sexP", "sexM"}:
            continue
        hits.setdefault(genome, {})
        current = hits[genome].get(gene)
        if current is None or to_float(row["bitscore"]) > to_float(current["bitscore"]):
            hits[genome][gene] = row
    return hits


def map_hit(row: dict[str, str] | None, id_to_index: dict[str, int]) -> int | None:
    if row is None:
        return None
    for candidate in subject_candidates(row):
        if candidate in id_to_index:
            return id_to_index[candidate]
    return None


def choose_hmg(hit_indices: dict[str, int | None], hits: dict[str, dict[str, str]]) -> str:
    choices = []
    for gene in ("sexP", "sexM"):
        if hit_indices.get(gene) is not None and gene in hits:
            choices.append((to_float(hits[gene]["bitscore"]), gene))
    if not choices:
        return ""
    return sorted(choices, reverse=True)[0][1]


def local_status(
    hmg_gene: str,
    mapped: dict[str, int | None],
    genes: list[dict[str, object]],
    max_neighbor_distance: int,
) -> tuple[str, str, str]:
    if not hmg_gene:
        return "hmg_not_mapped_to_gff", "", ""
    hmg_index = mapped[hmg_gene]
    if hmg_index is None:
        return "hmg_not_mapped_to_gff", "", ""

    hmg = genes[hmg_index]
    present_local = []
    distances = []
    for flank in ("tptA", "rnhA"):
        flank_index = mapped.get(flank)
        if flank_index is None:
            continue
        flank_gene = genes[flank_index]
        if flank_index == hmg_index:
            distances.append(f"{flank}:same_gene_as_hmg")
            continue
        if flank_gene["seqid"] != hmg["seqid"]:
            distances.append(f"{flank}:different_scaffold")
            continue
        distance = int(flank_gene["order_index"]) - int(hmg["order_index"])
        distances.append(f"{flank}:{distance}")
        if abs(distance) <= max_neighbor_distance:
            present_local.append(flank)

    if {"tptA", "rnhA"}.issubset(set(present_local)):
        return "complete_local_context", ",".join(present_local), ";".join(distances)
    if present_local:
        return "partial_local_context", ",".join(present_local), ";".join(distances)
    return "hmg_mapped_no_local_flanks", "", ";".join(distances)


def hit_text(row: dict[str, str] | None) -> str:
    if row is None:
        return ""
    return (
        f"{row['sseqid']}|q={row['qseqid']}|pident={row['pident']}|"
        f"qcov={row['qcovs']}|bitscore={row['bitscore']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--best-per-gene", required=True, type=Path)
    parser.add_argument("--annotation-samples", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--neighbor-genes", default=10, type=int)
    args = parser.parse_args()

    final_genomes = [
        row
        for row in read_csv(args.metadata)
        if row.get("include_final_ultrametric_tree") == "yes" or row.get("include_core_analyses") == "yes"
    ]
    gff_paths = load_gff_paths(args.annotation_samples)
    hits_by_genome = best_hits_by_genome_gene(args.best_per_gene)

    summary_rows = []
    neighbor_rows = []
    for meta in sorted(final_genomes, key=lambda row: row.get("tree_label", "") or row["genome_id"]):
        genome = meta["genome_id"]
        gff = gff_paths.get(genome)
        hits = hits_by_genome.get(genome, {})
        mapped = {"tptA": None, "rnhA": None, "sexP": None, "sexM": None}
        genes: list[dict[str, object]] = []

        if gff is not None:
            genes, id_to_index = parse_gff(gff)
            for gene in mapped:
                mapped[gene] = map_hit(hits.get(gene), id_to_index)

        hmg_gene = choose_hmg(mapped, hits)
        status, local_flanks, distances = local_status(hmg_gene, mapped, genes, args.neighbor_genes)

        def mapped_gene_id(gene_symbol: str) -> str:
            index = mapped.get(gene_symbol)
            if index is None or not genes:
                return ""
            return str(genes[index]["gene_id"])

        summary_rows.append(
            {
                "genome_id": genome,
                "tree_label": meta.get("tree_label", ""),
                "species_name": meta.get("species_name", ""),
                "order": meta.get("order", ""),
                "family": meta.get("family", ""),
                "gff_path": str(gff or ""),
                "context_status": "gff_missing" if gff is None else status,
                "selected_hmg_query_type": hmg_gene,
                "local_flanks": local_flanks,
                "flank_gene_distance_from_hmg": distances,
                "hmg_gene_id": mapped_gene_id(hmg_gene) if hmg_gene else "",
                "tptA_gene_id": mapped_gene_id("tptA"),
                "rnhA_gene_id": mapped_gene_id("rnhA"),
                "sexP_gene_id": mapped_gene_id("sexP"),
                "sexM_gene_id": mapped_gene_id("sexM"),
                "tptA_hit": hit_text(hits.get("tptA")),
                "rnhA_hit": hit_text(hits.get("rnhA")),
                "sexP_hit": hit_text(hits.get("sexP")),
                "sexM_hit": hit_text(hits.get("sexM")),
            }
        )

        hmg_index = mapped.get(hmg_gene) if hmg_gene else None
        if hmg_index is not None and genes:
            hmg = genes[hmg_index]
            scaffold_indices = [i for i, gene in enumerate(genes) if gene["seqid"] == hmg["seqid"]]
            scaffold_position = scaffold_indices.index(hmg_index)
            left = max(0, scaffold_position - args.neighbor_genes)
            right = min(len(scaffold_indices), scaffold_position + args.neighbor_genes + 1)
            for neighbor_index in scaffold_indices[left:right]:
                neighbor = genes[neighbor_index]
                neighbor_rows.append(
                    {
                        "genome_id": genome,
                        "tree_label": meta.get("tree_label", ""),
                        "species_name": meta.get("species_name", ""),
                        "selected_hmg_query_type": hmg_gene,
                        "hmg_gene_id": hmg["gene_id"],
                        "relative_gene_index": int(neighbor["order_index"]) - int(hmg["order_index"]),
                        "seqid": neighbor["seqid"],
                        "start": neighbor["start"],
                        "end": neighbor["end"],
                        "strand": neighbor["strand"],
                        "gene_id": neighbor["gene_id"],
                        "name": neighbor["name"],
                        "locus_tag": neighbor["locus_tag"],
                        "annotation_text": neighbor["annotation_text"],
                    }
                )

    summary_path = args.output_prefix.with_suffix(".summary.tsv")
    neighbors_path = args.output_prefix.with_suffix(".neighbors.tsv")
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_fields = [
        "genome_id",
        "tree_label",
        "species_name",
        "order",
        "family",
        "gff_path",
        "context_status",
        "selected_hmg_query_type",
        "local_flanks",
        "flank_gene_distance_from_hmg",
        "hmg_gene_id",
        "tptA_gene_id",
        "rnhA_gene_id",
        "sexP_gene_id",
        "sexM_gene_id",
        "tptA_hit",
        "rnhA_hit",
        "sexP_hit",
        "sexM_hit",
    ]
    neighbor_fields = [
        "genome_id",
        "tree_label",
        "species_name",
        "selected_hmg_query_type",
        "hmg_gene_id",
        "relative_gene_index",
        "seqid",
        "start",
        "end",
        "strand",
        "gene_id",
        "name",
        "locus_tag",
        "annotation_text",
    ]

    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    with neighbors_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=neighbor_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(neighbor_rows)

    print(f"Wrote {summary_path} ({len(summary_rows)} rows)")
    print(f"Wrote {neighbors_path} ({len(neighbor_rows)} rows)")


if __name__ == "__main__":
    main()
