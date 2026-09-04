#!/bin/bash
#SBATCH --job-name=cdhit
#SBATCH --account=nn8013k
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G

set -euo pipefail

module load CD-HIT/4.8.1-GCC-14.2.0

cd /cluster/projects/nn8013k/results/helle/repeats/consensus/

cd-hit-est \
  -i all_mucoro_repeat_families.fa \
  -o all_mucoro_repeat_familie_clusters.fa \
  -c 0.85 \
  -aS 0.8 \
  -T 12 \
  -M 0
