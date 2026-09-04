#!/bin/bash
#SBATCH --job-name=prep_cen_queries
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/prep_cen_queries_%j.out
#SBATCH --error=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/prep_cen_queries_%j.err

set -euo pipefail

PROJECT_ROOT=/cluster/work/users/hellb/Mucor
SCRIPT_ROOT="${PROJECT_ROOT}/scripts/assembly_structure"
RESULTS_DIR="${PROJECT_ROOT}/results/assembly_structure"
QUERY_DIR="${RESULTS_DIR}/centromere_queries"

mkdir -p "${RESULTS_DIR}/logs" "${QUERY_DIR}"

module --quiet purge
module load Python/3.13.5-GCCcore-14.3.0

python "${SCRIPT_ROOT}/scripts/extract_bed_sequences.py" \
  --bed "${SCRIPT_ROOT}/config/MucCirc3_CEN_repeats.bed" \
  --manifest "${PROJECT_ROOT}/data/metadata/genome_file_manifest.tsv" \
  --genome-id MucCirc3 \
  --out-all "${QUERY_DIR}/MucCirc3_CEN_and_GremLINE_queries.fasta" \
  --out-cen "${QUERY_DIR}/MucCirc3_CEN_core_queries.fasta" \
  --out-gremline "${QUERY_DIR}/MucCirc3_GremLINE_queries.fasta" \
  --out-index "${QUERY_DIR}/MucCirc3_CEN_and_GremLINE_query_index.tsv"

python "${SCRIPT_ROOT}/scripts/scan_fasta_motif.py" \
  --fasta "${QUERY_DIR}/MucCirc3_CEN_core_queries.fasta" \
  --motif GCATAGGT \
  --out "${QUERY_DIR}/MucCirc3_CEN_core_GCATAGGT_motif_hits.tsv"
