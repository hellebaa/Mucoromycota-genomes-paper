#!/usr/bin/env python3
"""Scan local GFF/GFF3 annotations for putative MAT-locus genes.

This annotation-first screen complements sequence homology searches. It looks
for annotated HMG/MAT-like genes and reports neighboring genes, especially
nearby `tptA` and `rnhA` candidates.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


SIGNAL_PATTERNS = {
    "HMG": re.compile(
        r"\b(high mobility group|high-mobility group|HMG[- ]?box|MAT1-2-1|sex plus|sex minus|sexP|sexM|Casanova)\b",
        re.IGNORECASE,
    ),
    "tptA": re.compile(
        r"\b(triose phosphate|triose-phosphate|phosphate transporter|phosphate translocator|TPTA|tptA)\b",
        re.IGNORECASE,
    ),
    "rnhA": re.compile(
        r"\b(RNA helicase|RNH|rnhA|DEAD.?box|DEAH.?box)\b",
        re.IGNORECASE,
    ),
}

FEATURE_TYPES = {"gene", "mRNA", "transcript", "CDS"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-genomes-dir", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument(
        "--neighbor-genes",
        type=int,
        default=10,
        help="Number of neighboring genes to report on each side of each HMG/MAT candidate.",
    )
    return parser.parse_args()


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


def clean_join(values: list[str]) -> str:
    seen = []
    for value in values:
        value = value.strip()
        if value and value not in seen:
            seen.append(value)
    return " | ".join(seen)


def useful_annotation_text(attrs: dict[str, str]) -> str:
    keys = [
        "ID",
        "Name",
        "Alias",
        "locus_tag",
        "protein_id",
        "orig_protein_id",
        "transcriptId",
        "proteinId",
        "product",
        "product_name",
        "Note",
        "note",
        "description",
        "Dbxref",
        "Ontology_term",
        "go_function",
        "go_process",
    ]
    return clean_join([attrs.get(key, "") for key in keys])


def split_parents(parent_value: str) -> list[str]:
    return [part.strip() for part in parent_value.split(",") if part.strip()]


def match_signals(text: str) -> list[str]:
    text_lower = text.lower()
    signals = []
    for signal, pattern in SIGNAL_PATTERNS.items():
        if signal == "HMG" and ("hmg-coa" in text_lower or "hmg coa" in text_lower):
            continue
        if pattern.search(text):
            signals.append(signal)
    return signals


def parse_gff(path: Path) -> list[dict[str, object]]:
    gene_rows: list[dict[str, object]] = []
    gene_index_by_id: dict[str, int] = {}
    child_to_gene: dict[str, str] = {}
    pending_child_text: dict[str, list[str]] = defaultdict(list)

    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attrs_raw = fields
            if feature_type not in FEATURE_TYPES:
                continue
            attrs = parse_attributes(attrs_raw)
            feature_id = attrs.get("ID", "")
            parents = split_parents(attrs.get("Parent", ""))
            text = useful_annotation_text(attrs)

            if feature_type == "gene":
                gene_id = feature_id or attrs.get("Name", "") or attrs.get("locus_tag", "")
                gene = {
                    "seqid": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene_id": gene_id,
                    "name": attrs.get("Name", ""),
                    "locus_tag": attrs.get("locus_tag", ""),
                    "protein_ids": [],
                    "annotation_texts": [text],
                }
                gene_index_by_id[gene_id] = len(gene_rows)
                gene_rows.append(gene)
                if gene_id in pending_child_text:
                    gene["annotation_texts"].extend(pending_child_text.pop(gene_id))
                continue

            parent_gene_id = ""
            if feature_type in {"mRNA", "transcript"} and parents:
                parent_gene_id = parents[0]
                if feature_id:
                    child_to_gene[feature_id] = parent_gene_id
            elif feature_type == "CDS" and parents:
                parent = parents[0]
                parent_gene_id = child_to_gene.get(parent, parent)

            if not parent_gene_id:
                continue
            protein_id = attrs.get("protein_id", "") or attrs.get("orig_protein_id", "") or attrs.get("proteinId", "")
            if parent_gene_id in gene_index_by_id:
                gene = gene_rows[gene_index_by_id[parent_gene_id]]
                gene["annotation_texts"].append(text)
                if protein_id:
                    gene["protein_ids"].append(protein_id)
            else:
                pending_child_text[parent_gene_id].append(text)

    gene_rows.sort(key=lambda row: (str(row["seqid"]), int(row["start"]), int(row["end"])))
    for order_index, gene in enumerate(gene_rows):
        text = clean_join(gene["annotation_texts"])  # type: ignore[arg-type]
        gene["annotation_text"] = text
        gene["signals"] = match_signals(text)
        gene["protein_ids"] = clean_join(gene["protein_ids"])  # type: ignore[arg-type]
        gene["order_index"] = order_index
    return gene_rows


def sample_id_from_gff(path: Path, raw_genomes_dir: Path) -> str:
    rel = path.relative_to(raw_genomes_dir)
    folder = rel.parts[0]
    stem = path.name
    for suffix in (".gff3", ".gff", ".gtf"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    if stem and stem != folder:
        return stem
    return folder


def locus_status(window_signals: set[str]) -> str:
    if {"HMG", "tptA", "rnhA"}.issubset(window_signals):
        return "annotation_complete_candidate_context"
    if "HMG" in window_signals and ({"tptA", "rnhA"} & window_signals):
        return "annotation_partial_candidate_context"
    if "HMG" in window_signals:
        return "annotation_hmg_only"
    return "no_hmg_annotation"


def main() -> None:
    args = parse_args()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    gff_files = sorted(
        path
        for path in args.raw_genomes_dir.rglob("*")
        if path.is_file() and path.name.endswith((".gff", ".gff3", ".gtf"))
    )

    candidate_rows = []
    neighbor_rows = []
    sample_rows = []

    for gff in gff_files:
        sample_id = sample_id_from_gff(gff, args.raw_genomes_dir)
        genes = parse_gff(gff)
        signal_genes = [gene for gene in genes if gene["signals"]]
        hmg_genes = [gene for gene in signal_genes if "HMG" in gene["signals"]]
        sample_best_status = "no_hmg_annotation"

        for gene in signal_genes:
            candidate_rows.append(
                {
                    "sample_id": sample_id,
                    "gff_path": str(gff),
                    "seqid": gene["seqid"],
                    "start": gene["start"],
                    "end": gene["end"],
                    "strand": gene["strand"],
                    "gene_id": gene["gene_id"],
                    "name": gene["name"],
                    "locus_tag": gene["locus_tag"],
                    "protein_ids": gene["protein_ids"],
                    "signals": ",".join(gene["signals"]),
                    "annotation_text": gene["annotation_text"],
                }
            )

        indices_by_scaffold = defaultdict(list)
        for index, gene in enumerate(genes):
            indices_by_scaffold[gene["seqid"]].append(index)

        for hmg in hmg_genes:
            index = int(hmg["order_index"])
            scaffold_indices = indices_by_scaffold[hmg["seqid"]]
            scaffold_position = scaffold_indices.index(index)
            left_pos = max(0, scaffold_position - args.neighbor_genes)
            right_pos = min(len(scaffold_indices), scaffold_position + args.neighbor_genes + 1)
            window_indices = scaffold_indices[left_pos:right_pos]
            window = [genes[window_index] for window_index in window_indices]
            window_signals = {signal for gene in window for signal in gene["signals"]}
            status = locus_status(window_signals)
            if status == "annotation_complete_candidate_context":
                sample_best_status = status
            elif status == "annotation_partial_candidate_context" and sample_best_status != "annotation_complete_candidate_context":
                sample_best_status = status
            elif sample_best_status == "no_hmg_annotation":
                sample_best_status = status

            for neighbor in window:
                neighbor_rows.append(
                    {
                        "sample_id": sample_id,
                        "hmg_gene_id": hmg["gene_id"],
                        "hmg_seqid": hmg["seqid"],
                        "hmg_start": hmg["start"],
                        "relative_gene_index": scaffold_indices.index(int(neighbor["order_index"])) - scaffold_position,
                        "seqid": neighbor["seqid"],
                        "start": neighbor["start"],
                        "end": neighbor["end"],
                        "strand": neighbor["strand"],
                        "gene_id": neighbor["gene_id"],
                        "name": neighbor["name"],
                        "locus_tag": neighbor["locus_tag"],
                        "protein_ids": neighbor["protein_ids"],
                        "signals": ",".join(neighbor["signals"]),
                        "annotation_text": neighbor["annotation_text"],
                        "window_status": status,
                    }
                )

        sample_rows.append(
            {
                "sample_id": sample_id,
                "gff_path": str(gff),
                "gene_count": len(genes),
                "annotated_signal_gene_count": len(signal_genes),
                "hmg_candidate_count": len(hmg_genes),
                "sample_best_status": sample_best_status,
            }
        )

    outputs = {
        args.output_prefix.with_suffix(".samples.tsv"): (
            sample_rows,
            [
                "sample_id",
                "gff_path",
                "gene_count",
                "annotated_signal_gene_count",
                "hmg_candidate_count",
                "sample_best_status",
            ],
        ),
        args.output_prefix.with_suffix(".candidates.tsv"): (
            candidate_rows,
            [
                "sample_id",
                "gff_path",
                "seqid",
                "start",
                "end",
                "strand",
                "gene_id",
                "name",
                "locus_tag",
                "protein_ids",
                "signals",
                "annotation_text",
            ],
        ),
        args.output_prefix.with_suffix(".neighbors.tsv"): (
            neighbor_rows,
            [
                "sample_id",
                "hmg_gene_id",
                "hmg_seqid",
                "hmg_start",
                "relative_gene_index",
                "seqid",
                "start",
                "end",
                "strand",
                "gene_id",
                "name",
                "locus_tag",
                "protein_ids",
                "signals",
                "annotation_text",
                "window_status",
            ],
        ),
    }

    for path, (rows, fields) in outputs.items():
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {path} ({len(rows)} rows)")

    print(f"GFF files scanned: {len(gff_files)}")


if __name__ == "__main__":
    main()
