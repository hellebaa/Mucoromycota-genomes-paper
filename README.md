# Mucoromycota Genomes Paper

Scripts and analyses for our comparative genomics paper on Mucoromycota fungi genomes.

## Repository structure

```
analysis/          One subdirectory per analysis
  synteny/         Synteny analysis
    download_data.sh   Download input files from storage
scripts/
  download_files.sh  Shared file downloader (currently NeLS, later Zenodo)
  nels.R             R helpers for reading NeLS files interactively
data/              Downloaded data files — not tracked by git
```

## Data access

Data files are not stored in the git repo but rather stored remotely and downloaded on demand. 
During development we will use NeLS. (See [NELS_SETUP.md](NELS_SETUP.md) for first-time SSH setup.)
Before publication we will move the files to Zenodo.

Each analysis directory contains a `download_data.sh` script that fetches the
files it needs:

```bash
bash analysis/synteny/download_data.sh
```

### Reading NeLS files from R

For interactive exploration, [`scripts/nels.R`](scripts/nels.R) reads files
from NeLS directly over SSH, using the same `nels.storage` host. Paths are the
same as in the repo:

```r
source("scripts/nels.R")
nels_ls("data")                                  # list a folder
sp <- nels_read("data/species.tsv")              # stream a table (data.table::fread)
f  <- nels_file("data/phylogeny/sco_cds_b500.date.nwk")  # download once to data/, local path
```

Pipelines should still fetch their inputs with `download_data.sh`, so that
they keep working when the data moves to Zenodo.

## Adding a new analysis

1. Create `analysis/<name>/`
2. Add a `download_data.sh` that calls `scripts/download_files.sh` with the
   required `data/` paths
3. Put analysis scripts in `analysis/<name>/`

## Shared Hi-C reference annotations

The [three-species Hi-C reference release](analysis/hic/README.md) provides candidate-centromere and telomere BED files, evidence tables and checksum-verified NeLS downloads.
