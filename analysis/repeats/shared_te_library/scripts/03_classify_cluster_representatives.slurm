#!/bin/bash
#SBATCH --job-name=TEclass
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=step3_classify.%j.out

set -euo pipefail

cd /cluster/projects/nn8013k/results/helle/repeats/consensus/

module load RepeatModeler/2.0.4-foss-2022a

RepeatClassifier \
  -consensi all_mucoro_repeat_familie_clusters.fa \
  -threads 16
