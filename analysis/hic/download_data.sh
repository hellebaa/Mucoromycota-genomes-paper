#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
case "${1:-}" in
 "") MANIFEST="$ROOT/analysis/hic/config/reference_v1_files.tsv" ;;
 --inputs) MANIFEST="$ROOT/analysis/hic/config/reference_v1_inputs.tsv" ;;
 *) echo "Usage: $0 [--inputs]" >&2; exit 1 ;;
esac
bash "$ROOT/scripts/download_files.sh" --manifest "$MANIFEST"
