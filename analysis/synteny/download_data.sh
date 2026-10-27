#!/usr/bin/env bash
# Download data files required for the synteny analysis.
# Run from any directory; paths are resolved relative to the repo root.
#
# 1. Phylogeny, metadata and orthogroups, which define the genome set.
# 2. R/00_genome_table.R writes genomes.tsv from the phylogeny and lists the
#    gene tables and pairwise orthologue tables needed for the analysed pairs.
# 3. Those files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

scripts/download_files.sh \
    data/orthofinder/results/Species_Tree/SpeciesTree_rooted_node_labels.txt \
    data/orthofinder/results/Orthogroups/Orthogroups.tsv \
    data/species.tsv \
    data/phylogeny/sco_cds_withOutgroups_b500.treefile \
    data/phylogeny/sco_cds_b500.date.nwk \
    data/metadata/README.md \
    data/metadata/genome_metadata_master.csv \
    data/metadata/assembly_qc_summary.csv \
    data/metadata/genome_id_aliases.csv

FILES=()
while IFS= read -r f; do FILES+=("${f}"); done < <(Rscript analysis/synteny/R/00_genome_table.R)
scripts/download_files.sh "${FILES[@]}"
