#!/bin/bash
#SBATCH --job-name=nels2saga
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=nels2saga.%j.out

set -euo pipefail

NELS_HOST="${NELS_HOST:-nels.storage}"
NELS_DATA_DIR="${NELS_DATA_DIR:-Projects/NMBU_threeD_Sandve_2024/Mucor/data/}"
SAGA_DATA_DIR="${SAGA_DATA_DIR:-/cluster/work/users/hellb/Mucor/data/}"
APPLY="${APPLY:-0}"

mkdir -p "${SAGA_DATA_DIR}"

RSYNC_ARGS=(
    -avh
    --progress
    --delete
)

if [[ "${APPLY}" != "1" ]]; then
    RSYNC_ARGS+=(--dry-run)
    echo "Dry run only. Submit with --export=ALL,APPLY=1 to transfer files."
fi

rsync "${RSYNC_ARGS[@]}" \
    "${NELS_HOST}:${NELS_DATA_DIR}" \
    "${SAGA_DATA_DIR}"

cat <<EOF

Next step on Saga:
  python analysis/data_curation/scripts/build_genome_file_manifest.py \\
    --data-root ${SAGA_DATA_DIR%/} \\
    --metadata ${SAGA_DATA_DIR%/}/metadata/genome_metadata_master.csv \\
    --out ${SAGA_DATA_DIR%/}/metadata/genome_file_manifest.tsv
EOF
