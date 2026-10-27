#!/usr/bin/env bash
# Run the full synteny analysis: download inputs, compute indices, statistics,
# figures and the HTML report. Run from any directory.
#
#   bash analysis/synteny/run_all.sh            # everything
#   SKIP_DOWNLOAD=1 bash analysis/synteny/run_all.sh

set -euo pipefail
# UTF-8 locale; under the C locale knitr escapes non-ASCII characters in figure
# captions (e.g. a less-than-or-equal sign is printed as <U+2264>).
export LC_ALL=C.UTF-8

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"
export REPO_ROOT

if [[ -z "${SKIP_DOWNLOAD:-}" ]]; then
    bash analysis/synteny/download_data.sh   # also writes genomes.tsv
else
    Rscript analysis/synteny/R/00_genome_table.R > /dev/null
fi

Rscript analysis/synteny/R/test_functions.R
Rscript analysis/synteny/R/01_pairs_and_qc.R
Rscript analysis/synteny/R/02_pair_indices.R
Rscript analysis/synteny/R/03_statistics.R
Rscript analysis/synteny/R/04_figures.R
Rscript -e 'rmarkdown::render("analysis/synteny/synteny_report.Rmd", knit_root_dir = Sys.getenv("REPO_ROOT"), quiet = TRUE)'
