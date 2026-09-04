#!/usr/bin/env bash
# Mirror the canonical NeLS data directory to Saga working storage.

set -euo pipefail

NELS_HOST="${NELS_HOST:-nels.storage}"
NELS_DATA_DIR="${NELS_DATA_DIR:-Projects/NMBU_threeD_Sandve_2024/Mucor/data/}"
SAGA_DATA_DIR="${SAGA_DATA_DIR:-/cluster/work/users/hellb/Mucor/data/}"
DRY_RUN=1

usage() {
    cat <<EOF
Usage: $(basename "$0") [--apply]

Mirrors:
  ${NELS_HOST}:${NELS_DATA_DIR}
to:
  ${SAGA_DATA_DIR}

Options:
  --apply   Actually copy/update files. Without this flag, rsync runs as dry run.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply)
            DRY_RUN=0
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

mkdir -p "${SAGA_DATA_DIR}"

RSYNC_ARGS=(
    -avh
    --progress
    --delete
)

if [[ "${DRY_RUN}" -eq 1 ]]; then
    RSYNC_ARGS+=(--dry-run)
    echo "Dry run only. Re-run with --apply to transfer files."
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
