#!/bin/bash
#SBATCH --job-name=asm_struct
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/asm_struct_%j.out
#SBATCH --error=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/asm_struct_%j.err

set -euo pipefail

PROJECT_ROOT=/cluster/work/users/hellb/Mucor
DATA_ROOT="${PROJECT_ROOT}/data"
SCRIPT_ROOT="${PROJECT_ROOT}/scripts/assembly_structure/scripts"
RESULTS_DIR="${PROJECT_ROOT}/results/assembly_structure"

mkdir -p "${RESULTS_DIR}/logs"

module --quiet purge
module load Python/3.13.5-GCCcore-14.3.0

python "${SCRIPT_ROOT}/summarize_assembly_contigs.py" \
  --manifest "${DATA_ROOT}/metadata/genome_file_manifest.tsv" \
  --include core \
  --out-contigs "${RESULTS_DIR}/assembly_contig_gene_summary.tsv" \
  --out-genomes "${RESULTS_DIR}/assembly_genome_summary.tsv"
