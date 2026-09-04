#!/bin/bash
#SBATCH --job-name=RM
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-63
#SBATCH --output=RM_%A_%a.out

set -euo pipefail

module load RepeatMasker/4.1.4-foss-2022a

cd /cluster/projects/nn8013k/results/helle/repeats

mkdir -p RM_output

GENOMES=(genomes/*.f*a)
GENOME=${GENOMES[$SLURM_ARRAY_TASK_ID-1]}

SCRATCH=/cluster/projects/nn8013k/results/helle/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCRATCH"

export REPEATMASKER_CACHE_DIR=$SCRATCH/RM_cache
mkdir -p "$REPEATMASKER_CACHE_DIR"

RepeatMasker \
  -engine ncbi \
  -no_is \
  -lib consensus/mucoro_TE_shared.fa \
  -gff \
  -xsmall \
  -pa 8 \
  -dir RM_output \
  "$GENOME"

rm -rf "$SCRATCH"
