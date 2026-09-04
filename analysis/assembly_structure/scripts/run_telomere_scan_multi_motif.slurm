#!/bin/bash
#SBATCH --job-name=telomeres_multi
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/telomeres_multi_%j.out
#SBATCH --error=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/telomeres_multi_%j.err

set -euo pipefail

PROJECT_ROOT=/cluster/work/users/hellb/Mucor
SCRIPT_ROOT="${PROJECT_ROOT}/scripts/assembly_structure"
RESULTS_DIR="${PROJECT_ROOT}/results/assembly_structure"

mkdir -p "${RESULTS_DIR}/logs"

module --quiet purge
module load Python/3.13.5-GCCcore-14.3.0

python "${SCRIPT_ROOT}/scripts/scan_telomere_motifs.py" \
  --manifest "${PROJECT_ROOT}/data/metadata/genome_file_manifest.tsv" \
  --include core \
  --motifs "${SCRIPT_ROOT}/config/telomere_motifs.tsv" \
  --window-bp 25000 \
  --max-distance-bp 5000 \
  --min-copies 3 \
  --out-bed "${RESULTS_DIR}/telomere_hits_multi_motif.bed" \
  --out-summary "${RESULTS_DIR}/telomere_contig_summary_multi_motif.tsv"
