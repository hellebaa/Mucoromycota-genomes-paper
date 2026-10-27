#!/usr/bin/env bash
# Download specific files from NeLS storage.
#
# Usage:
#   download_files.sh <local/path/relative/to/repo/root> [<more/files> ...]
#
# Each argument is a path relative to the repo root (e.g. data/orthofinder/results/foo.txt).
# The file is fetched from the corresponding path on NeLS:
#   nels.storage:Projects/NMBU_threeD_Sandve_2024/Mucor/data/<rest of path>
#
# All files are transferred in a single rsync call (one ssh connection).
# -k follows symlinked directories on NeLS (e.g. data/orthofinder/results), so
# they arrive as real directories.
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

echo "==> Downloading $# files from ${NELS_HOST}"
printf '%s\n' "$@" |
    rsync -avhk --files-from=- "${NELS_HOST}:${NELS_REMOTE_BASE}/" "${REPO_ROOT}/"
