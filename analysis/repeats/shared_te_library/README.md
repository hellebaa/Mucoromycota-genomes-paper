# Shared TE Library

This analysis documents how the comparative Mucoromycota TE library was built
and then used for uniform RepeatMasker annotation across genomes.

The historical final shared library on Saga is:

```text
/cluster/projects/nn8013k/results/helle/repeat/consensus/mucoro_TE_shared.fa
```

Large input/output files live on NeLS and Saga, not in Git.

## Canonical NeLS Inputs

Species-specific de novo repeat annotations are stored per genome:

```text
Projects/NMBU_threeD_Sandve_2024/Mucor/data/genomes/<genome_id>/repeats/denovo/
  <genome_id>.repeatmodeler-library.fasta
  <genome_id>.repeatmasker-denovo.gff3
  <genome_id>.repeatmasker-denovo.out
  <genome_id>.repeatmasker-denovo.tbl
```

The `<genome_id>.repeatmodeler-library.fasta` files correspond to
RepeatModeler `consensi.fa.classified`/`consensi.fa` outputs from the
species-specific de novo runs.

## Workflow Outline

### Step 1: Collect Species-Specific TE Consensuses

For each genome, de novo repeat libraries were generated with RepeatModeler
v2.0.4. The resulting species-specific consensus FASTA files were extracted
and concatenated manually into one input FASTA:

```bash
cat species*/consensi.fa > all_mucoro_repeat_families.fa
```

No historical SLURM/script file has been identified for this step. The remaining
workflow steps were run as scripts on Saga.

### Step 2: Cluster Consensus Sequences

All consensus sequences were clustered with CD-HIT-EST v4.8.1 at 85%
nucleotide identity and 80% alignment coverage of the shorter sequence:

```text
-c 0.85 -aS 0.8
```

Cluster representatives were retained as provisional shared TE family
sequences.

Historical Saga script:

```bash
#!/bin/bash
#SBATCH --job-name=cdhit
#SBATCH --account=nn8013k
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G

module load CD-HIT/4.8.1-GCC-14.2.0

cd /cluster/projects/nn8013k/results/helle/repeats/consensus/

cd-hit-est \
  -i all_mucoro_repeat_families.fa \
  -o all_mucoro_repeat_familie_clusters.fa \
  -c 0.85 \
  -aS 0.8 \
  -T 12 \
  -M 0
```

### Step 3: Classify Cluster Representatives

Cluster representative sequences were classified with RepeatClassifier from
RepeatModeler v2.0.4 using default parameters.

Historical Saga script:

```bash
#!/bin/bash
#SBATCH --job-name=TEclass
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=step3_classify.%j.out

cd /cluster/projects/nn8013k/results/helle/repeats/consensus/

module load RepeatModeler/2.0.4-foss-2022a

RepeatClassifier \
  -consensi all_mucoro_repeat_familie_clusters.fa \
  -threads 16
```

This produced:

```text
mucoro_TE_clusters.fa.classified
```

### Step 4: Filter and Rename Shared Library

Sequences annotated as `Simple_repeat`, `Low_complexity`, `Satellite`, `rRNA`,
or `tRNA` were removed. Remaining sequences were renamed with standardized TE
family identifiers while preserving classification information, for example:

```text
MUC_Gypsy_0001
```

Filtering and renaming were run on Saga with:

```bash
#!/bin/bash
#SBATCH --job-name=TErename
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=step4_rename.%j.out

cd /cluster/projects/nn8013k/results/helle/repeats/consensus/

module load Python/3.13.5-GCCcore-14.3.0
module load Biopython/1.86-gfbf-2025b

grep -v -E "Simple_repeat|Low_complexity|Satellite|rRNA|tRNA" \
  all_mucoro_repeat_familie_clusters.fa.classified > mucoro_TE_clusters.filtered.fa

python /cluster/projects/nn8013k/results/helle/repeats/rename_te_library.py
```

The Python script read `mucoro_TE_clusters.filtered.fa` and wrote the final
shared TE library:

```text
mucoro_TE_shared.fa
```

Note: Step 3 records `mucoro_TE_clusters.fa.classified` as the produced file,
while the Step 4 script filtered
`all_mucoro_repeat_familie_clusters.fa.classified`. Confirm whether this was a
rename/copy before final methods freeze.

### Step 5: Re-Annotate Genomes with Shared Library

All genomes were re-annotated with RepeatMasker using the curated shared TE
library under identical parameters across species.

Historical Saga script:

```bash
#!/bin/bash
#SBATCH --job-name=RM
#SBATCH --account=nn8013k
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-63
#SBATCH --output=RM_%A_%a.out

module load RepeatMasker/4.1.4-foss-2022a

cd /cluster/projects/nn8013k/results/helle/repeats

mkdir -p RM_output

GENOMES=(genomes/*.f*a)
GENOME=${GENOMES[$SLURM_ARRAY_TASK_ID-1]}

SCRATCH=/cluster/projects/nn8013k/results/helle/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $SCRATCH

export REPEATMASKER_CACHE_DIR=$SCRATCH/RM_cache
mkdir -p $REPEATMASKER_CACHE_DIR

RepeatMasker \
  -engine ncbi \
  -no_is \
  -lib consensus/mucoro_TE_shared.fa \
  -gff \
  -xsmall \
  -pa 8 \
  -dir RM_output \
  $GENOME

rm -rf $SCRATCH
```

This generated shared-library RepeatMasker outputs for 63 genomes in:

```text
/cluster/projects/nn8013k/results/helle/repeats/RM_output/
```

The outputs are being organized on NeLS per genome:

```text
Projects/NMBU_threeD_Sandve_2024/Mucor/data/genomes/<genome_id>/repeats/shared_library/
  <genome_id>.sharedlib.repeatmasker.out
  <genome_id>.sharedlib.repeatmasker.out.gff
  <genome_id>.sharedlib.repeatmasker.tbl
  <genome_id>.sharedlib.repeatmasker.masked.fasta
```

## Notes To Complete

- Confirm Step 3/Step 4 classified FASTA filename rename/copy.
- Add the final NeLS location of `mucoro_TE_shared.fa`.
- Add a manifest of genomes included in the library build.
- Add a manifest of genomes re-annotated with the shared library.
