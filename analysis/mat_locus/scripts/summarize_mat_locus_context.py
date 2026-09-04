#!/usr/bin/env python3
"""Summarize MAT candidate locus context from BLAST hits and GFF annotations."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organized-samples", required=True, type=Path)
    parser.add_argument("--blast-summary", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--min-bitscore", default=50.0, type=float)
    parser.add_argument("--min-query-coverage", default=25.0, type=float)
    parser.add_argument(
        "--neighbor-genes",
        default=10,
        type=int,
        help="Number of neighboring genes to report on each side of each HMG candidate.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def normalized_ids(raw_id: str) -> set[str]:
    ids = {raw_id}
    ids.add(raw_id.split()[0])
    ids.add(re.sub(r"-T\d+$", "", raw_id.split()[0]))
    ids.add(re.sub(r"\.p\d+$", "", raw_id.split()[0]))
    ids.add(re.sub(r"\.t\d+$", "", raw_id.split()[0]))
    return {x for x in ids if x}


def canonical_signal(gene_symbol: str) -> str:
    if gene_symbol in {"sexP", "sexM"}:
        return "HMG"
    return gene_symbol


def parse_gff(path: Path) -> tuple[list[dict[str, str]], dict[str, int]]:
    genes = []
    id_to_index = {}
    transcript_to_gene = {}

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
            parent = attrs.get("Parent", "")

            if feature_type in {"mRNA", "transcript"} and feature_id and parent:
                transcript_to_gene[feature_id] = parent.split(",")[0]

            if feature_type != "gene":
                continue

            gene_id = feature_id or attrs.get("Name", "") or attrs.get("locus_tag", "")
            product = attrs.get("product", "") or attrs.get("Note", "") or attrs.get("description", "")
            gene = {
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "gene_id": gene_id,
                "name": attrs.get("Name", ""),
                "locus_tag": attrs.get("locus_tag", ""),
                "product": product,
                "attributes": attrs_raw,
            }
            genes.append(gene)

    genes.sort(key=lambda row: (row["seqid"], row["start"], row["end"]))
    for index, gene in enumerate(genes):
        for key in ("gene_id", "name", "locus_tag"):
            for ident in normalized_ids(gene.get(key, "")):
                id_to_index.setdefault(ident, index)

    for transcript, gene_id in transcript_to_gene.items():
        if gene_id in id_to_index:
            id_to_index[transcript] = id_to_index[gene_id]
            id_to_index[re.sub(r"-T\d+$", "", transcript)] = id_to_index[gene_id]

    return genes, id_to_index


def index_hits(rows: list[dict[str, str]], min_bitscore: float, min_qcov: float) -> dict[str, list[dict[str, str]]]:
    by_genome = defaultdict(list)
    for row in rows:
        try:
            bitscore = float(row.get("bitscore", "0") or 0)
            qcov = float(row.get("query_coverage", "0") or 0)
        except ValueError:
            continue
        if bitscore < min_bitscore or qcov < min_qcov:
            continue
        row = dict(row)
        row["signal"] = canonical_signal(row.get("gene_symbol", ""))
        by_genome[row["genome_id"]].append(row)
    return by_genome


def best_by_signal(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    best = {}
    for row in rows:
        signal = row["signal"]
        bitscore = float(row["bitscore"])
        current = best.get(signal)
        if current is None or bitscore > float(current["bitscore"]):
            best[signal] = row
    return best


def find_gene_index(hit: dict[str, str], id_to_index: dict[str, int]) -> int | None:
    subject = hit.get("subject_id", "")
    title = hit.get("subject_title", "")
    candidates = set()
    for token in [subject, title, *title.split()]:
        candidates.update(normalized_ids(token))
    for candidate in candidates:
        if candidate in id_to_index:
            return id_to_index[candidate]
    return None


def locus_status(signals: set[str], same_scaffold: bool) -> str:
    if {"HMG", "tptA", "rnhA"}.issubset(signals) and same_scaffold:
        return "complete_candidate_context"
    if "HMG" in signals and ({"tptA", "rnhA"} & signals) and same_scaffold:
        return "partial_candidate_context"
    if "HMG" in signals:
        return "hmg_only_or_no_gff_context"
    return "no_hmg_candidate"


def main() -> None:
    args = parse_args()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    samples = read_tsv(args.organized_samples)
    hits_by_genome = index_hits(
        read_tsv(args.blast_summary),
        min_bitscore=args.min_bitscore,
        min_qcov=args.min_query_coverage,
    )

    summary_rows = []
    neighbor_rows = []

    for sample in samples:
        genome_id = sample["genome_id"]
        gff_path = Path(sample.get("annotation_gff", ""))
        genes = []
        id_to_index = {}
        if gff_path.exists():
            genes, id_to_index = parse_gff(gff_path)

        best = best_by_signal(hits_by_genome.get(genome_id, []))
        mapped = {}
        for signal, hit in best.items():
            index = find_gene_index(hit, id_to_index) if id_to_index else None
            if index is not None:
                mapped[signal] = (hit, genes[index], index)

        scaffolds = {gene["seqid"] for hit, gene, index in mapped.values()}
        same_scaffold = len(scaffolds) == 1 and len(mapped) >= 2
        signals = set(mapped)

        hmg_records = [value for signal, value in mapped.items() if signal == "HMG"]
        for hit, gene, index in hmg_records:
            left = max(0, index - args.neighbor_genes)
            right = min(len(genes), index + args.neighbor_genes + 1)
            for neighbor_index in range(left, right):
                neighbor = genes[neighbor_index]
                neighbor_rows.append(
                    {
                        "genome_id": genome_id,
                        "tree_label": sample.get("tree_label", ""),
                        "species_name": sample.get("species_name", ""),
                        "hmg_subject_id": hit.get("subject_id", ""),
                        "hmg_gene_id": gene["gene_id"],
                        "relative_gene_index": neighbor_index - index,
                        "seqid": neighbor["seqid"],
                        "start": neighbor["start"],
                        "end": neighbor["end"],
                        "strand": neighbor["strand"],
                        "gene_id": neighbor["gene_id"],
                        "name": neighbor["name"],
                        "locus_tag": neighbor["locus_tag"],
                        "product": neighbor["product"],
                    }
                )

        summary_rows.append(
            {
                "genome_id": genome_id,
                "tree_label": sample.get("tree_label", ""),
                "species_name": sample.get("species_name", ""),
                "annotation_gff": str(gff_path) if gff_path.exists() else "",
                "locus_status": locus_status(signals, same_scaffold),
                "same_scaffold_candidate": "yes" if same_scaffold else "no",
                "signals_mapped_to_gff": ",".join(sorted(signals)),
                "hmg_gene_id": mapped.get("HMG", ("", {"gene_id": ""}, ""))[1]["gene_id"] if "HMG" in mapped else "",
                "tptA_gene_id": mapped.get("tptA", ("", {"gene_id": ""}, ""))[1]["gene_id"] if "tptA" in mapped else "",
                "rnhA_gene_id": mapped.get("rnhA", ("", {"gene_id": ""}, ""))[1]["gene_id"] if "rnhA" in mapped else "",
                "hmg_subject_id": best.get("HMG", {}).get("subject_id", ""),
                "tptA_subject_id": best.get("tptA", {}).get("subject_id", ""),
                "rnhA_subject_id": best.get("rnhA", {}).get("subject_id", ""),
                "notes": "" if gff_path.exists() else "annotation GFF missing",
            }
        )

    summary_path = args.output_prefix.with_suffix(".summary.tsv")
    neighbors_path = args.output_prefix.with_suffix(".neighbors.tsv")

    summary_fields = [
        "genome_id",
        "tree_label",
        "species_name",
        "annotation_gff",
        "locus_status",
        "same_scaffold_candidate",
        "signals_mapped_to_gff",
        "hmg_gene_id",
        "tptA_gene_id",
        "rnhA_gene_id",
        "hmg_subject_id",
        "tptA_subject_id",
        "rnhA_subject_id",
        "notes",
    ]
    neighbor_fields = [
        "genome_id",
        "tree_label",
        "species_name",
        "hmg_subject_id",
        "hmg_gene_id",
        "relative_gene_index",
        "seqid",
        "start",
        "end",
        "strand",
        "gene_id",
        "name",
        "locus_tag",
        "product",
    ]

    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    with neighbors_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=neighbor_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(neighbor_rows)

    print(f"Wrote {summary_path}")
    print(f"Wrote {neighbors_path}")


if __name__ == "__main__":
    main()
