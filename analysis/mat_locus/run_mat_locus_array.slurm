#!/usr/bin/env bash
#SBATCH --job-name=mat_locus
#SBATCH --account=CHANGE_ME
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-80
#SBATCH --output=logs/mat_locus_%A_%a.out
#SBATCH --error=logs/mat_locus_%A_%a.err

set -euo pipefail

# Edit these paths for the cluster.
PROJECT_DIR="${PROJECT_DIR:-/path/to/Mucoromycota-genomes-paper}"
GENOMES_DIR="${GENOMES_DIR:-/path/to/data/genomes}"
METADATA_CSV="${METADATA_CSV:-${PROJECT_DIR}/data/metadata/genome_metadata_master.csv}"
QUERY_FASTA="${QUERY_FASTA:-/path/to/MAT_MucCirc._prot.fasta}"
MAT_WORK="${MAT_WORK:-/path/to/work/Mucoromycota/MAT_locus}"

mkdir -p "${MAT_WORK}/logs" "${MAT_WORK}/results/summary"
cd "${MAT_WORK}"

# Adjust module names to the cluster environment.
module purge
module load BLAST+/2.15.0 || module load blast
module load Python/3.11 || module load python

MANIFEST="${MAT_WORK}/sample_manifest.tsv"
ORGANIZED="${MAT_WORK}/inputs/organized_samples.tsv"
QUERY_ORGANIZED="${MAT_WORK}/inputs/queries/$(basename "${QUERY_FASTA}")"

if [[ ! -s "${MANIFEST}" ]]; then
    python "${PROJECT_DIR}/analysis/mat_locus/scripts/build_sample_manifest.py" \
        --metadata "${METADATA_CSV}" \
        --genomes-dir "${GENOMES_DIR}" \
        --output "${MANIFEST}"
fi

if [[ ! -s "${ORGANIZED}" ]]; then
    python "${PROJECT_DIR}/analysis/mat_locus/scripts/prepare_cluster_inputs.py" \
        --manifest "${MANIFEST}" \
        --query-fasta "${QUERY_FASTA}" \
        --work-dir "${MAT_WORK}" \
        --mode symlink
fi

"${PROJECT_DIR}/analysis/mat_locus/scripts/scan_one_sample.sh" \
    "${ORGANIZED}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "${QUERY_ORGANIZED}" \
    "${MAT_WORK}/results"

if [[ "${SLURM_ARRAY_TASK_ID}" == "1" ]]; then
    echo "Array jobs started. Run the summarizer after all array jobs finish:"
    echo "python ${PROJECT_DIR}/analysis/mat_locus/scripts/summarize_mat_hits.py \\"
    echo "  --organized-samples ${ORGANIZED} \\"
    echo "  --query-map ${PROJECT_DIR}/analysis/mat_locus/config/mat_query_accessions.tsv \\"
    echo "  --results-dir ${MAT_WORK}/results \\"
    echo "  --output ${MAT_WORK}/results/summary/mat_locus_hit_summary.tsv"
fi
