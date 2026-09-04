# BLAST Databases

This module builds reusable BLAST databases from the Saga mirror of the NeLS
data directory. These databases are shared infrastructure for MAT-locus
annotation, kinetochore gene mining, and later targeted gene searches.

## Inputs

Build the genome file manifest first:

```text
/cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv
```

The manifest is made by:

```bash
python /cluster/work/users/hellb/Mucor/scripts/build_genome_file_manifest.py
```

## Outputs

BLAST databases are written outside `data/`, because they are derived temporary
analysis infrastructure:

```text
/cluster/work/users/hellb/Mucor/results/blastdb/
  proteins/<genome_id>/
  assemblies/<genome_id>/
  logs/
```

## Transfer Scripts To Saga

From the laptop/local machine:

```bash
rsync -ravz \
  /Users/hellebaalsrud/.codex/worktrees/a8cb/Mucoromycota-genomes-paper/analysis/blastdb/ \
  hellb@login-1.saga.sigma2.no:/cluster/work/users/hellb/Mucor/scripts/blastdb/
```

## Run On Saga

Submit a dry run first:

```bash
mkdir -p /cluster/work/users/hellb/Mucor/results/blastdb/logs
sbatch /cluster/work/users/hellb/Mucor/scripts/blastdb/scripts/build_blastdb_array.slurm
```

If the selected genomes look right, submit the real build:

```bash
sbatch --export=ALL,APPLY=1 /cluster/work/users/hellb/Mucor/scripts/blastdb/scripts/build_blastdb_array.slurm
```

The SLURM array is currently set to `1-75`, matching the 75-row genome
manifest.

By default, the array builds databases only for genomes marked as included in
core analyses or the final ultrametric tree. To include all genomes with usable
files, submit with:

```bash
sbatch --export=ALL,APPLY=1,INCLUDE_ALL=1 /cluster/work/users/hellb/Mucor/scripts/blastdb/scripts/build_blastdb_array.slurm
```

## Databases Created

For each selected genome:

- `makeblastdb -dbtype prot` from `proteins_faa`
- `makeblastdb -dbtype nucl` from `assembly_fasta`

The databases are built without `-parse_seqids` to avoid failures from mixed
public/private FASTA header formats.

Saga module:

```bash
module load BLAST+/2.17.0-gompi-2025b
```

The file:

```text
/cluster/work/users/hellb/Mucor/results/blastdb/blastdb_manifest.tsv
```

records the expected database prefixes for downstream SLURM jobs.

After the array finishes, merge per-task manifest rows with:

```bash
bash /cluster/work/users/hellb/Mucor/scripts/blastdb/scripts/merge_blastdb_manifest.sh
```
