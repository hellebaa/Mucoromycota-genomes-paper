#!/usr/bin/env python3
"""Create a first-pass MAT/sex evidence review table.

This intentionally produces review categories rather than final biological
calls. HMG-only BLAST hits are noisy, so MAT status still requires manual locus
context inspection.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_delimited(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def load_best_by_gene(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows = read_delimited(path, "\t")
    best: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["genome_id"], row["gene_symbol"])
        bitscore = to_float(row["bitscore"])
        if key not in best or bitscore > to_float(best[key]["bitscore"]):
            best[key] = row
    return best


def load_best_exploratory_hmg(path: Path) -> dict[str, dict[str, str]]:
    rows = read_delimited(path, "\t")
    best: dict[str, dict[str, str]] = {}
    for row in rows:
        if row["locus_role"] != "exploratory_hmg":
            continue
        genome_id = row["genome_id"]
        bitscore = to_float(row["bitscore"])
        if genome_id not in best or bitscore > to_float(best[genome_id]["bitscore"]):
            best[genome_id] = row
    return best


def classify_hit(row: dict[str, str] | None, *, kind: str) -> str:
    if row is None:
        return "absent"

    qcov = to_float(row.get("qcovs", ""))
    pident = to_float(row.get("pident", ""))
    bitscore = to_float(row.get("bitscore", ""))

    if kind == "hmg":
        if qcov >= 60 and bitscore >= 100 and pident >= 25:
            return "strong"
        if qcov >= 30 and bitscore >= 50 and pident >= 25:
            return "moderate"
        return "weak"

    if kind == "flank":
        if qcov >= 60 and bitscore >= 150 and pident >= 30:
            return "strong"
        if qcov >= 30 and bitscore >= 70 and pident >= 25:
            return "moderate"
        return "weak"

    if kind == "sex_pathway":
        if qcov >= 70 and bitscore >= 150 and pident >= 30:
            return "strong"
        if qcov >= 40 and bitscore >= 80 and pident >= 25:
            return "moderate"
        return "weak"

    return "present"


def hit_summary(row: dict[str, str] | None) -> str:
    if row is None:
        return ""
    return (
        f"{row['sseqid']} "
        f"q={row['qseqid']} "
        f"pident={row['pident']} "
        f"qcov={row['qcovs']} "
        f"bitscore={row['bitscore']}"
    )


def infer_review_status(statuses: dict[str, str], annotation_status: str) -> tuple[str, str, str]:
    hmg = [statuses["sexP"], statuses["sexM"]]
    flanks = [statuses["tptA"], statuses["rnhA"]]
    strong_hmg = [name for name in ("sexP", "sexM") if statuses[name] == "strong"]
    moderate_hmg = [name for name in ("sexP", "sexM") if statuses[name] in {"strong", "moderate"}]
    strong_flank_count = sum(status == "strong" for status in flanks)
    moderate_flank_count = sum(status in {"strong", "moderate"} for status in flanks)

    if strong_hmg and strong_flank_count == 2:
        allele = "/".join(strong_hmg)
        confidence = "high" if annotation_status == "annotation_complete_candidate_context" else "medium"
        return "candidate_mat_locus", confidence, f"{allele} with strong tptA and rnhA homology"

    if annotation_status == "annotation_complete_candidate_context" and moderate_hmg:
        allele = "/".join(moderate_hmg)
        return "candidate_mat_locus", "medium", f"{allele} with complete candidate context in prior GFF scan"

    if moderate_hmg and moderate_flank_count >= 1:
        allele = "/".join(moderate_hmg)
        return "homology_hits_context_needed", "low", f"{allele} plus flanking-gene homologs; inspect genomic neighborhood before MAT call"

    if any(status in {"strong", "moderate"} for status in hmg):
        return "hmg_hit_context_needed", "low", "HMG-like hit found; broad HMG homology is not enough for a MAT call"

    if annotation_status in {"annotation_complete_candidate_context", "annotation_partial_candidate_context"}:
        return "annotation_context_candidate", "low", "Prior GFF annotation scan suggested MAT-like context"

    return "no_clear_mat_candidate", "low", "No clear MAT-locus candidate from this first-pass protein BLAST summary"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--best-per-gene", required=True, type=Path)
    parser.add_argument("--annotation-samples", required=True, type=Path)
    parser.add_argument("--context-summary", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    metadata_rows = read_delimited(args.metadata, ",")
    final_rows = [
        row
        for row in metadata_rows
        if row.get("include_final_ultrametric_tree") == "yes"
        or row.get("include_core_analyses") == "yes"
    ]

    best_by_gene = load_best_by_gene(args.best_per_gene)
    exploratory_hmg = load_best_exploratory_hmg(args.best_per_gene)
    annotation_rows = {
        row["sample_id"]: row
        for row in read_delimited(args.annotation_samples, "\t")
    }
    context_rows = {}
    if args.context_summary and args.context_summary.exists():
        context_rows = {
            row["genome_id"]: row
            for row in read_delimited(args.context_summary, "\t")
        }

    output_fields = [
        "genome_id",
        "tree_label",
        "species_name",
        "strain_id",
        "family",
        "order",
        "metadata_sexual_reproduction",
        "metadata_mat_locus",
        "review_mat_status",
        "review_confidence",
        "review_notes",
        "annotation_scan_status",
        "annotation_hmg_candidate_count",
        "blast_gff_context_status",
        "blast_gff_selected_hmg",
        "blast_gff_local_flanks",
        "blast_gff_flank_distance_from_hmg",
        "blast_gff_hmg_gene_id",
        "blast_gff_tptA_gene_id",
        "blast_gff_rnhA_gene_id",
        "tptA_hit_class",
        "rnhA_hit_class",
        "sexP_hit_class",
        "sexM_hit_class",
        "exploratory_hmg_hit_class",
        "TSP1_hit_class",
        "TSP2_hit_class",
        "tptA_best_hit",
        "rnhA_best_hit",
        "sexP_best_hit",
        "sexM_best_hit",
        "exploratory_hmg_best_hit",
        "TSP1_best_hit",
        "TSP2_best_hit",
    ]

    output_rows = []
    for row in sorted(final_rows, key=lambda r: r["tree_label"] or r["genome_id"]):
        genome_id = row["genome_id"]
        annotation = annotation_rows.get(genome_id, {})
        context = context_rows.get(genome_id, {})
        hits = {
            "tptA": best_by_gene.get((genome_id, "tptA")),
            "rnhA": best_by_gene.get((genome_id, "rnhA")),
            "sexP": best_by_gene.get((genome_id, "sexP")),
            "sexM": best_by_gene.get((genome_id, "sexM")),
            "exploratory_hmg": exploratory_hmg.get(genome_id),
            "tdh/TSP1": best_by_gene.get((genome_id, "tdh/TSP1")),
            "TSP2": best_by_gene.get((genome_id, "TSP2")),
        }
        statuses = {
            "tptA": classify_hit(hits["tptA"], kind="flank"),
            "rnhA": classify_hit(hits["rnhA"], kind="flank"),
            "sexP": classify_hit(hits["sexP"], kind="hmg"),
            "sexM": classify_hit(hits["sexM"], kind="hmg"),
            "exploratory_hmg": classify_hit(hits["exploratory_hmg"], kind="hmg"),
            "tdh/TSP1": classify_hit(hits["tdh/TSP1"], kind="sex_pathway"),
            "TSP2": classify_hit(hits["TSP2"], kind="sex_pathway"),
        }
        review_status, confidence, notes = infer_review_status(
            statuses,
            annotation.get("sample_best_status", ""),
        )
        context_status = context.get("context_status", "")
        if context_status == "complete_local_context":
            review_status = "candidate_mat_locus"
            confidence = "high" if annotation.get("sample_best_status") == "annotation_complete_candidate_context" else "medium"
            notes = f"{context.get('blast_gff_selected_hmg', context.get('selected_hmg_query_type', 'HMG'))} with local tptA and rnhA context"
        elif context_status == "partial_local_context" and review_status != "candidate_mat_locus":
            review_status = "possible_mat_locus_partial_context"
            confidence = "low"
            notes = "HMG candidate has one local MAT-flanking gene; inspect locus manually"
        output_rows.append(
            {
                "genome_id": genome_id,
                "tree_label": row.get("tree_label", ""),
                "species_name": row.get("species_name", ""),
                "strain_id": row.get("strain_id", ""),
                "family": row.get("family", ""),
                "order": row.get("order", ""),
                "metadata_sexual_reproduction": row.get("sexual_reproduction", ""),
                "metadata_mat_locus": row.get("mat_locus", ""),
                "review_mat_status": review_status,
                "review_confidence": confidence,
                "review_notes": notes,
                "annotation_scan_status": annotation.get("sample_best_status", ""),
                "annotation_hmg_candidate_count": annotation.get("hmg_candidate_count", ""),
                "blast_gff_context_status": context.get("context_status", ""),
                "blast_gff_selected_hmg": context.get("selected_hmg_query_type", ""),
                "blast_gff_local_flanks": context.get("local_flanks", ""),
                "blast_gff_flank_distance_from_hmg": context.get("flank_gene_distance_from_hmg", ""),
                "blast_gff_hmg_gene_id": context.get("hmg_gene_id", ""),
                "blast_gff_tptA_gene_id": context.get("tptA_gene_id", ""),
                "blast_gff_rnhA_gene_id": context.get("rnhA_gene_id", ""),
                "tptA_hit_class": statuses["tptA"],
                "rnhA_hit_class": statuses["rnhA"],
                "sexP_hit_class": statuses["sexP"],
                "sexM_hit_class": statuses["sexM"],
                "exploratory_hmg_hit_class": statuses["exploratory_hmg"],
                "TSP1_hit_class": statuses["tdh/TSP1"],
                "TSP2_hit_class": statuses["TSP2"],
                "tptA_best_hit": hit_summary(hits["tptA"]),
                "rnhA_best_hit": hit_summary(hits["rnhA"]),
                "sexP_best_hit": hit_summary(hits["sexP"]),
                "sexM_best_hit": hit_summary(hits["sexM"]),
                "exploratory_hmg_best_hit": hit_summary(hits["exploratory_hmg"]),
                "TSP1_best_hit": hit_summary(hits["tdh/TSP1"]),
                "TSP2_best_hit": hit_summary(hits["TSP2"]),
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=output_fields)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {args.out}")
    print(f"Rows: {len(output_rows)}")


if __name__ == "__main__":
    main()
