#!/usr/bin/env bash
# Make a read-only manifest of the canonical NeLS genomes folder.

set -euo pipefail

NELS_HOST="${NELS_HOST:-nels.storage}"
NELS_GENOMES_DIR="${NELS_GENOMES_DIR:-Projects/NMBU_threeD_Sandve_2024/Mucor/data/genomes}"
OUT="${1:-outputs/nels_genomes_manifest.tsv}"

mkdir -p "$(dirname "${OUT}")"

ssh "${NELS_HOST}" "python3 -c '
import os
base = \"${NELS_GENOMES_DIR}\"
print(\"relative_path\\ttype\\tsize_bytes\")
for root, dirs, files in os.walk(base):
    for name in sorted(dirs):
        path = os.path.join(root, name)
        print(f\"{os.path.relpath(path, base)}\\tdir\\t\")
    for name in sorted(files):
        path = os.path.join(root, name)
        print(f\"{os.path.relpath(path, base)}\\tfile\\t{os.path.getsize(path)}\")
'" > "${OUT}"

echo "Wrote ${OUT}"
