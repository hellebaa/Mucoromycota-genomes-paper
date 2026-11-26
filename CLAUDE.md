# Project context for Claude

Comparative genomics paper on Mucoromycota fungi. Analyses are in Python/R/bash.

## Repository conventions

- One directory per analysis under `analysis/`. Each has a `download_data.sh`
  that lists the input files it needs from remote storage.
- Data files live under `data/<analysis>/` and are never committed (`.gitignore`).
- `scripts/download_files.sh` is the single entry point for fetching remote
  files. It currently downloads from NeLS storage; it will be replaced with a
  Zenodo downloader at publication. Analysis scripts must not reference NeLS directly — always go through
  `scripts/download_files.sh`.
- `scripts/nels.R` (`nels_ls`, `nels_read`, `nels_file`, `nels_readRDS`) reads
  NeLS files from R for interactive exploration. Do not use it in pipeline
  scripts; they fetch inputs through `download_data.sh`.

## Adding files to an analysis download list

Edit the relevant `analysis/<name>/download_data.sh` and add paths as
additional arguments to the `scripts/download_files.sh` call. Paths are
relative to the repo root and must start with `data/`.

## Data storage

Remote storage: NeLS, host alias `nels.storage` (see `NELS_SETUP.md`).
Remote base path: `~/Projects/NMBU_threeD_Sandve_2024/Mucor/data/`
Local mirror: `data/`
