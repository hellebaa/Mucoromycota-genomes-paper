#!/bin/bash
#SBATCH --job-name=TErename
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=step4_rename.%j.out

set -euo pipefail

cd /cluster/projects/nn8013k/results/helle/repeats/consensus/

module load Python/3.13.5-GCCcore-14.3.0
module load Biopython/1.86-gfbf-2025b

# Remove simple repeats, low-complexity sequence, satellites, and RNA repeats.
grep -v -E "Simple_repeat|Low_complexity|Satellite|rRNA|tRNA" \
  all_mucoro_repeat_familie_clusters.fa.classified > mucoro_TE_clusters.filtered.fa

python /cluster/projects/nn8013k/results/helle/repeats/rename_te_library.py
