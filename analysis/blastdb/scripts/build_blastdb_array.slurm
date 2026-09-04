#!/bin/bash
#SBATCH --job-name=blastdb
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=1-75
#SBATCH --output=/cluster/work/users/hellb/Mucor/results/blastdb/logs/blastdb_%A_%a.out
#SBATCH --error=/cluster/work/users/hellb/Mucor/results/blastdb/logs/blastdb_%A_%a.err

set -euo pipefail

module load BLAST+/2.17.0-gompi-2025b

MANIFEST="${MANIFEST:-/cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv}"
BLASTDB_ROOT="${BLASTDB_ROOT:-/cluster/work/users/hellb/Mucor/results/blastdb}"
SCRIPT_DIR="${SCRIPT_DIR:-/cluster/work/users/hellb/Mucor/scripts/blastdb/scripts}"
APPLY="${APPLY:-0}"

mkdir -p "${BLASTDB_ROOT}/logs"

if [[ "${APPLY}" != "1" ]]; then
    echo "Dry run only. Submit with --export=ALL,APPLY=1 to build BLAST databases."
    awk -F '\t' -v n="${SLURM_ARRAY_TASK_ID}" '
        NR == 1 {
            for (i = 1; i <= NF; i++) idx[$i] = i
            next
        }
        NR == n + 1 {
            print "row=" n
            print "genome_id=" $idx["genome_id"]
            print "status=" $idx["status"]
            print "include_core_analyses=" $idx["include_core_analyses"]
            print "include_final_ultrametric_tree=" $idx["include_final_ultrametric_tree"]
            print "proteins_faa=" $idx["proteins_faa"]
            print "assembly_fasta=" $idx["assembly_fasta"]
        }
    ' "${MANIFEST}"
    exit 0
fi

"${SCRIPT_DIR}/build_blastdb_one_genome.sh" \
    "${MANIFEST}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "${BLASTDB_ROOT}" > "${BLASTDB_ROOT}/logs/blastdb_manifest.${SLURM_ARRAY_TASK_ID}.tsv"
