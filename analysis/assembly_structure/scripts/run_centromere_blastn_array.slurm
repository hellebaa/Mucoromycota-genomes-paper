#!/bin/bash
#SBATCH --job-name=cen_blastn
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=1-61
#SBATCH --output=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/cen_blastn_%A_%a.out
#SBATCH --error=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/cen_blastn_%A_%a.err

set -euo pipefail

PROJECT_ROOT=/cluster/work/users/hellb/Mucor
MANIFEST="${PROJECT_ROOT}/data/metadata/genome_file_manifest.tsv"
QUERY_FASTA="${PROJECT_ROOT}/results/assembly_structure/centromere_queries/MucCirc3_CEN_core_queries.fasta"
RESULTS_DIR="${PROJECT_ROOT}/results/assembly_structure"
BLAST_DIR="${RESULTS_DIR}/centromere_blastn_cen_core"

mkdir -p "${RESULTS_DIR}/logs" "${BLAST_DIR}"

module --quiet purge
module load BLAST+/2.17.0-gompi-2025b

mapfile -t GENOME_ROWS < <(awk -F'\t' 'NR > 1 && $5 == "yes" && $10 != "" {print $1 "\t" $10}' "${MANIFEST}")

ROW="${GENOME_ROWS[$SLURM_ARRAY_TASK_ID-1]}"
GENOME_ID="$(printf '%s\n' "${ROW}" | cut -f1)"
ASSEMBLY_FASTA="$(printf '%s\n' "${ROW}" | cut -f2)"
OUT="${BLAST_DIR}/${GENOME_ID}.MucCirc3_CEN_core.blastn.tsv"

echo "Genome: ${GENOME_ID}"
echo "Assembly: ${ASSEMBLY_FASTA}"
echo "Output: ${OUT}"

blastn \
  -query "${QUERY_FASTA}" \
  -subject "${ASSEMBLY_FASTA}" \
  -out "${OUT}" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -perc_identity 75 \
  -evalue 1e-10 \
  -num_threads "${SLURM_CPUS_PER_TASK}"
