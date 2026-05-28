#!/usr/bin/env bash
# Download specific files from NeLS storage.
#
# Usage:
#   nels_download.sh <local/path/relative/to/repo/root> [<more/files> ...]
#
# Each argument is a path relative to the repo root (e.g. data/orthofinder/results/foo.txt).
# The file is fetched from the corresponding path on NeLS:
#   nels.storage:Projects/NMBU_threeD_Sandve_2024/Mucor/data/<rest of path>
#
# Requires ~/.ssh/config to have a Host entry named "nels.storage" (see NELS_SETUP.md).

set -euo pipefail

NELS_HOST="nels.storage"
NELS_REMOTE_BASE="Projects/NMBU_threeD_Sandve_2024/Mucor"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ $# -eq 0 ]]; then
    echo "Usage: $(basename "$0") <local/path> [<local/path> ...]" >&2
    exit 1
fi

for local_rel in "$@"; do
    local_abs="${REPO_ROOT}/${local_rel}"
    remote_path="${NELS_REMOTE_BASE}/${local_rel}"

    mkdir -p "$(dirname "${local_abs}")"

    echo "==> ${local_rel}"
    rsync -avh --progress "${NELS_HOST}:${remote_path}" "${local_abs}"
done
