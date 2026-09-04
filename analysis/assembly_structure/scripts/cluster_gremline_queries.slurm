#!/bin/bash
#SBATCH --job-name=cluster_gremline
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/cluster_gremline_%j.out
#SBATCH --error=/cluster/work/users/hellb/Mucor/results/assembly_structure/logs/cluster_gremline_%j.err

set -euo pipefail

PROJECT_ROOT=/cluster/work/users/hellb/Mucor
QUERY_DIR="${PROJECT_ROOT}/results/assembly_structure/centromere_queries"

mkdir -p "${PROJECT_ROOT}/results/assembly_structure/logs"

module --quiet purge
module load CD-HIT/4.8.1-GCC-14.2.0

cd-hit-est \
  -i "${QUERY_DIR}/MucCirc3_GremLINE_queries.fasta" \
  -o "${QUERY_DIR}/MucCirc3_GremLINE_queries.cdhit90_rep.fasta" \
  -c 0.90 \
  -aS 0.80 \
  -T "${SLURM_CPUS_PER_TASK}" \
  -M 0

{
  printf "input_gremline_queries\t"
  grep -c '^>' "${QUERY_DIR}/MucCirc3_GremLINE_queries.fasta"
  printf "cdhit90_representatives\t"
  grep -c '^>' "${QUERY_DIR}/MucCirc3_GremLINE_queries.cdhit90_rep.fasta"
} > "${QUERY_DIR}/MucCirc3_GremLINE_query_counts.tsv"
