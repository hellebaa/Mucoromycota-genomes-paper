#!/usr/bin/env bash
# Download data files required for the synteny analysis.
# Run from any directory; paths are resolved relative to the repo root.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

"${REPO_ROOT}/scripts/download_files.sh" \
    data/orthofinder/results/Species_Tree/SpeciesTree_rooted_node_labels.txt \
    data/orthofinder/results/Orthogroups/Orthogroups.tsv
