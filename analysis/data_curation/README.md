# Data Curation And Saga Mirror

This module contains small utilities for keeping the project data layout
auditable across NeLS, Saga, and local development.

## Canonical Storage

NeLS is the canonical shared data store during analysis development:

```text
Projects/NMBU_threeD_Sandve_2024/Mucor/data/
```

Saga analyses should run from a disposable working mirror, not directly from
NeLS:

```text
/cluster/work/users/hellb/Mucor/data/
```

## Step 1: Mirror NeLS To Saga

Precondition: the metadata files should also be on NeLS:

```text
Projects/NMBU_threeD_Sandve_2024/Mucor/data/metadata/
  genome_metadata_master.csv
  genome_id_aliases.csv
  assembly_qc_summary.csv
  tree_set_audit.csv
```

If needed, upload the current local metadata staging folder to NeLS:

```bash
rsync -avh --progress \
  /Users/hellebaalsrud/Documents/Projects/Mucormycota_NeLS_staging/data/metadata/ \
  nels.storage:Projects/NMBU_threeD_Sandve_2024/Mucor/data/metadata/
```

Submit on Saga:

```bash
sbatch analysis/data_curation/scripts/sync_nels_to_saga.slurm
```

By default this runs a dry run. If the file list looks right, submit the actual
transfer:

```bash
sbatch --export=ALL,APPLY=1 analysis/data_curation/scripts/sync_nels_to_saga.slurm
```

The SLURM script uses a plain `rsync` command. If NeLS drops the SSH connection,
resubmit the same command; `rsync` will skip files that are already complete.

## Step 2: Build Genome File Manifest

After the mirror is complete, build the file manifest on Saga:

```bash
python analysis/data_curation/scripts/build_genome_file_manifest.py \
  --data-root /cluster/work/users/hellb/Mucor/data \
  --metadata /cluster/work/users/hellb/Mucor/data/metadata/genome_metadata_master.csv \
  --out /cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv
```

This produces one row per curated genome ID with concrete paths to assembly,
GFF, protein, CDS, repeat, ID-table, and gene-position files. Downstream SLURM
arrays should read this manifest instead of globbing directories independently.
