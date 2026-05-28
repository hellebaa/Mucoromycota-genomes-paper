# Mucoromycota Genomes Paper

Scripts and analyses for a comparative genomics paper on Mucoromycota fungi genomes.

## Repository structure

```
analysis/          One subdirectory per analysis
  synteny/         Synteny analysis
    download_data.sh   Download input files from storage
scripts/
  download_files.sh  Shared file downloader (currently NeLS, later Zenodo)
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

## Adding a new analysis

1. Create `analysis/<name>/`
2. Add a `download_data.sh` that calls `scripts/download_files.sh` with the
   required `data/` paths
3. Put analysis scripts in `analysis/<name>/`
