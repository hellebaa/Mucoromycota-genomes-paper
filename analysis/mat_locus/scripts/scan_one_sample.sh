#!/usr/bin/env bash
# Run MAT-locus homology searches for one sample from organized_samples.tsv.

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $(basename "$0") <organized_samples.tsv> <1-based-row-index> <query.faa> <results_dir>" >&2
    exit 1
fi

SAMPLES_TSV="$1"
ROW_INDEX="$2"
QUERY_FASTA="$3"
RESULTS_DIR="$4"

ROW="$(awk -v n="${ROW_INDEX}" 'BEGIN{FS=OFS="\t"} NR==n+1 {print; found=1} END{if(!found) exit 2}' "${SAMPLES_TSV}" || true)"
if [[ -z "${ROW}" ]]; then
    echo "No sample for row index ${ROW_INDEX}; exiting cleanly."
    exit 0
fi

GENOME_ID="$(printf '%s\n' "${ROW}" | cut -f1)"
ASSEMBLY_FASTA="$(printf '%s\n' "${ROW}" | cut -f4)"
PROTEIN_FASTA="$(printf '%s\n' "${ROW}" | cut -f6)"
READY="$(printf '%s\n' "${ROW}" | cut -f7)"

mkdir -p "${RESULTS_DIR}/blastp" "${RESULTS_DIR}/tblastn" "${RESULTS_DIR}/blastdb" "${RESULTS_DIR}/logs"

if [[ "${READY}" != "yes" ]]; then
    echo "${GENOME_ID}: inputs incomplete; skipping."
    exit 0
fi

PROT_DB="${RESULTS_DIR}/blastdb/${GENOME_ID}.proteins"
NUCL_DB="${RESULTS_DIR}/blastdb/${GENOME_ID}.assembly"

makeblastdb -in "${PROTEIN_FASTA}" -dbtype prot -out "${PROT_DB}" \
    > "${RESULTS_DIR}/logs/${GENOME_ID}.makeblastdb.prot.log" 2>&1

blastp \
    -query "${QUERY_FASTA}" \
    -db "${PROT_DB}" \
    -evalue 1e-5 \
    -max_target_seqs 20 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle" \
    -out "${RESULTS_DIR}/blastp/${GENOME_ID}.blastp.tsv"

makeblastdb -in "${ASSEMBLY_FASTA}" -dbtype nucl -out "${NUCL_DB}" \
    > "${RESULTS_DIR}/logs/${GENOME_ID}.makeblastdb.nucl.log" 2>&1

tblastn \
    -query "${QUERY_FASTA}" \
    -db "${NUCL_DB}" \
    -evalue 1e-5 \
    -max_target_seqs 20 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle" \
    -out "${RESULTS_DIR}/tblastn/${GENOME_ID}.tblastn.tsv"

echo "${GENOME_ID}: done"
