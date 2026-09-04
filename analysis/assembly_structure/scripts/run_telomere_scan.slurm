#!/bin/bash
#SBATCH --job-name=telomeres
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/telomeres_%j.out
#SBATCH --error=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/telomeres_%j.err

set -euo pipefail

PROJECT_ROOT=/cluster/work/users/hellb/Mucor
SCRIPT_ROOT="${PROJECT_ROOT}/scripts/assembly_structure/scripts"
RESULTS_DIR="${PROJECT_ROOT}/results/assembly_structure"

mkdir -p "${RESULTS_DIR}/logs"

module --quiet purge
module load Python/3.13.5-GCCcore-14.3.0

python "${SCRIPT_ROOT}/scan_telomere_motifs.py" \
  --manifest "${PROJECT_ROOT}/data/metadata/genome_file_manifest.tsv" \
  --include core \
  --motif TTAGGG \
  --reverse-motif CCCTAA \
  --window-bp 10000 \
  --max-distance-bp 1000 \
  --min-copies 3 \
  --out-bed "${RESULTS_DIR}/telomere_hits.bed" \
  --out-summary "${RESULTS_DIR}/telomere_contig_summary.tsv"
