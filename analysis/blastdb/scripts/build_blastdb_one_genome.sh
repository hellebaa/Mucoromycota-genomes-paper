#!/usr/bin/env bash
# Build protein and assembly nucleotide BLAST databases for one manifest row.

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $(basename "$0") <genome_file_manifest.tsv> <row_number> <blastdb_root>" >&2
    echo "row_number is 1-based, excluding the header." >&2
    exit 1
fi

MANIFEST="$1"
ROW_NUMBER="$2"
BLASTDB_ROOT="$3"

row=$(awk -F '\t' -v n="${ROW_NUMBER}" 'NR == n + 1 {print; exit}' "${MANIFEST}")
if [[ -z "${row}" ]]; then
    echo "No manifest row ${ROW_NUMBER} in ${MANIFEST}" >&2
    exit 1
fi

header=$(head -n 1 "${MANIFEST}")

value_for() {
    local column="$1"
    awk -F '\t' -v col="${column}" -v header="${header}" -v row="${row}" '
        BEGIN {
            n = split(header, h, "\t")
            split(row, r, "\t")
            for (i = 1; i <= n; i++) {
                if (h[i] == col) {
                    print r[i]
                    exit
                }
            }
        }
    '
}

GENOME_ID=$(value_for genome_id)
STATUS=$(value_for status)
IN_CORE=$(value_for include_core_analyses)
IN_TREE=$(value_for include_final_ultrametric_tree)
PROTEINS=$(value_for proteins_faa)
ASSEMBLY=$(value_for assembly_fasta)

if [[ "${STATUS}" != "ok" ]]; then
    echo "Skipping ${GENOME_ID}: status=${STATUS}" >&2
    exit 0
fi

if [[ "${INCLUDE_ALL:-0}" != "1" && "${IN_CORE}" != "yes" && "${IN_TREE}" != "yes" ]]; then
    echo "Skipping ${GENOME_ID}: not core/final; set INCLUDE_ALL=1 to include." >&2
    exit 0
fi

if [[ -z "${PROTEINS}" || ! -s "${PROTEINS}" ]]; then
    echo "Missing proteins_faa for ${GENOME_ID}: ${PROTEINS}" >&2
    exit 1
fi

if [[ -z "${ASSEMBLY}" || ! -s "${ASSEMBLY}" ]]; then
    echo "Missing assembly_fasta for ${GENOME_ID}: ${ASSEMBLY}" >&2
    exit 1
fi

PROT_DIR="${BLASTDB_ROOT}/proteins/${GENOME_ID}"
ASM_DIR="${BLASTDB_ROOT}/assemblies/${GENOME_ID}"
mkdir -p "${PROT_DIR}" "${ASM_DIR}"

PROT_PREFIX="${PROT_DIR}/${GENOME_ID}.proteins"
ASM_PREFIX="${ASM_DIR}/${GENOME_ID}.assembly"

makeblastdb \
    -in "${PROTEINS}" \
    -dbtype prot \
    -out "${PROT_PREFIX}" \
    -title "${GENOME_ID} proteins" >&2

makeblastdb \
    -in "${ASSEMBLY}" \
    -dbtype nucl \
    -out "${ASM_PREFIX}" \
    -title "${GENOME_ID} assembly" >&2

echo -e "${GENOME_ID}\t${PROT_PREFIX}\t${ASM_PREFIX}"
