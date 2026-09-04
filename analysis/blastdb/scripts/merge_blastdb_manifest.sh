#!/usr/bin/env bash
# Merge per-array-task BLAST DB manifest rows after build_blastdb_array.slurm.

set -euo pipefail

BLASTDB_ROOT="${BLASTDB_ROOT:-/cluster/work/users/hellb/Mucor/results/blastdb}"
OUT="${BLASTDB_ROOT}/blastdb_manifest.tsv"

mkdir -p "${BLASTDB_ROOT}"

{
    echo -e "genome_id\tprotein_blastdb_prefix\tassembly_blastdb_prefix"
    find "${BLASTDB_ROOT}/logs" -name "blastdb_manifest.*.tsv" -type f -print \
        | sort -V \
        | xargs awk -F '\t' 'NF == 3 && $1 !~ /^Building/ && $2 != "" && $3 != ""'
} > "${OUT}"

echo "Wrote ${OUT}"
