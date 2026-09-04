#!/usr/bin/env bash
# Download data files required for the overview figure.
# Run from any directory; paths are resolved relative to the repo root.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

"${REPO_ROOT}/scripts/download_files.sh" \
    data/orthofinder/results/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv \
    data/phylogeny/sco_cds_b500.date.nwk
