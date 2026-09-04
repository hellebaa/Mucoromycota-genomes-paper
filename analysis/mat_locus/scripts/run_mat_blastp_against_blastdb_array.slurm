#!/bin/bash
#SBATCH --job-name=mat_blastp
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=1-61
#SBATCH --output=/cluster/work/users/hellb/Mucor/results/mat_locus/logs/mat_blastp_%A_%a.out
#SBATCH --error=/cluster/work/users/hellb/Mucor/results/mat_locus/logs/mat_blastp_%A_%a.err

set -euo pipefail

module load BLAST+/2.17.0-gompi-2025b

BLASTDB_MANIFEST="${BLASTDB_MANIFEST:-/cluster/work/users/hellb/Mucor/results/blastdb/blastdb_manifest.tsv}"
QUERY_FASTA="${QUERY_FASTA:-/cluster/work/users/hellb/Mucor/results/mat_locus/queries/mat_locus_plus_sex_pathway_queries.faa}"
RESULTS_DIR="${RESULTS_DIR:-/cluster/work/users/hellb/Mucor/results/mat_locus}"
EVALUE="${EVALUE:-1e-5}"
MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-25}"

mkdir -p "${RESULTS_DIR}/blastp" "${RESULTS_DIR}/logs"

row=$(awk -F '\t' -v n="${SLURM_ARRAY_TASK_ID}" 'NR == n + 1 {print; exit}' "${BLASTDB_MANIFEST}")
if [[ -z "${row}" ]]; then
    echo "No BLAST DB manifest row ${SLURM_ARRAY_TASK_ID}; exiting cleanly."
    exit 0
fi

genome_id=$(printf "%s\n" "${row}" | cut -f1)
protein_db=$(printf "%s\n" "${row}" | cut -f2)

if [[ -z "${genome_id}" || -z "${protein_db}" ]]; then
    echo "Malformed BLAST DB manifest row ${SLURM_ARRAY_TASK_ID}: ${row}" >&2
    exit 1
fi

if [[ ! -s "${QUERY_FASTA}" ]]; then
    echo "Missing query FASTA: ${QUERY_FASTA}" >&2
    exit 1
fi

if [[ ! -s "${protein_db}.pin" && ! -s "${protein_db}.pdb" ]]; then
    echo "Missing protein BLAST DB prefix for ${genome_id}: ${protein_db}" >&2
    exit 1
fi

out="${RESULTS_DIR}/blastp/${genome_id}.mat_locus_plus_sex_pathway.blastp.tsv"

blastp \
    -query "${QUERY_FASTA}" \
    -db "${protein_db}" \
    -evalue "${EVALUE}" \
    -max_target_seqs "${MAX_TARGET_SEQS}" \
    -num_threads "${SLURM_CPUS_PER_TASK}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen stitle" \
    -out "${out}"

echo -e "${genome_id}\t${out}"
