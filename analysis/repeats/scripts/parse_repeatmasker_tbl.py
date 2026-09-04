#!/usr/bin/env python3
"""Parse RepeatMasker .tbl summaries into one repeat summary table."""

from __future__ import annotations

import argparse
import glob
import os
import re
from pathlib import Path

import pandas as pd


SUMMARY_FIELDS = {
    "Total_length_bp": 0,
    "GC_percent": 0.0,
    "Bases_masked_bp": 0,
    "Bases_masked_percent": 0.0,
    "Retroelements_bp": 0,
    "Retroelements_percent": 0.0,
    "DNA_transposons_bp": 0,
    "DNA_transposons_percent": 0.0,
    "Rolling_circles_bp": 0,
    "Rolling_circles_percent": 0.0,
    "Unclassified_bp": 0,
    "Unclassified_percent": 0.0,
}


def parse_int(value: str) -> int:
    return int(value.replace(",", ""))


def parse_tbl(filepath: str) -> dict[str, object]:
    species = os.path.basename(filepath).split(".")[0]
    data: dict[str, object] = {"Species": species, **SUMMARY_FIELDS}

    text = Path(filepath).read_text()

    match = re.search(r"total length:\s+([0-9,]+)", text)
    if match:
        data["Total_length_bp"] = parse_int(match.group(1))

    match = re.search(r"GC level:\s+([0-9.]+)", text)
    if match:
        data["GC_percent"] = float(match.group(1))

    match = re.search(r"bases masked:\s+([0-9,]+)\s+bp\s*\\(?\\s*([0-9.]+)", text)
    if match:
        data["Bases_masked_bp"] = parse_int(match.group(1))
        data["Bases_masked_percent"] = float(match.group(2))

    repeat_patterns = {
        "Retroelements": ("Retroelements_bp", "Retroelements_percent"),
        "DNA transposons": ("DNA_transposons_bp", "DNA_transposons_percent"),
        "Rolling-circles": ("Rolling_circles_bp", "Rolling_circles_percent"),
        "Unclassified:": ("Unclassified_bp", "Unclassified_percent"),
    }
    for label, (bp_field, percent_field) in repeat_patterns.items():
        match = re.search(rf"{re.escape(label)}\s+[0-9,]+\s+([0-9,]+)\s+bp\s+([0-9.]+)", text)
        if match:
            data[bp_field] = parse_int(match.group(1))
            data[percent_field] = float(match.group(2))

    other_bp = (
        int(data["Bases_masked_bp"])
        - int(data["Retroelements_bp"])
        - int(data["DNA_transposons_bp"])
        - int(data["Rolling_circles_bp"])
        - int(data["Unclassified_bp"])
    )
    other_percent = (
        float(data["Bases_masked_percent"])
        - float(data["Retroelements_percent"])
        - float(data["DNA_transposons_percent"])
        - float(data["Rolling_circles_percent"])
        - float(data["Unclassified_percent"])
    )

    data["Other_bp"] = other_bp
    data["Other_percent"] = round(other_percent, 2)
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        default="/cluster/projects/nn8013k/results/helle/repeats/repeat_tbl",
        help="Directory containing RepeatMasker .tbl files.",
    )
    parser.add_argument(
        "--output",
        default="repeat_summary_table.tsv",
        help="Output TSV path.",
    )
    args = parser.parse_args()

    rows = [parse_tbl(filepath) for filepath in glob.glob(os.path.join(args.input_dir, "*.tbl"))]
    if not rows:
        raise SystemExit(f"No .tbl files found in {args.input_dir}")

    df = pd.DataFrame(rows).sort_values("Species")
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
