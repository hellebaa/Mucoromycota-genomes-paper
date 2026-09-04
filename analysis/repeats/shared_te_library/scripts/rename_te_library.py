#!/usr/bin/env python3
"""Rename filtered TE consensus sequences with standardized family IDs."""

from collections import defaultdict

from Bio import SeqIO


input_fa = "mucoro_TE_clusters.filtered.fa"
output_fa = "mucoro_TE_shared.fa"

counts = defaultdict(int)

with open(output_fa, "w") as out:
    for rec in SeqIO.parse(input_fa, "fasta"):
        header = rec.description

        if "#" in header:
            teclass = header.split("#")[1].split()[0]
        else:
            teclass = "Unknown"

        teclass = teclass.replace("/", "_")
        counts[teclass] += 1

        rec.id = f"MUC_{teclass}_{counts[teclass]:04d}"
        rec.description = ""
        SeqIO.write(rec, out, "fasta")
