# MAT Locus Annotation

This analysis annotates candidate MAT locus genes and broader sex-pathway genes
across the Mucoromycota genome set.

The analysis uses two evidence tracks:

1. MAT-locus candidates
   - Mucorales: canonical `tptA` - HMG/`sexP`/`sexM` - `rnhA` context.
   - Mortierellales/Umbelopsidales: exploratory HMG/homeodomain candidates;
     absence of Mucorales-like hits should not be interpreted as MAT absence.
2. Sex-pathway candidates
   - Trisporoid-pathway genes such as `tdh`/`TSP1`, treated as evidence for
     sexual signalling machinery, not as MAT-locus calls.

The query catalogs are:

```text
analysis/mat_locus/config/mat_locus_query_catalog.tsv
analysis/mat_locus/config/sex_pathway_query_accessions.tsv
```

The legacy Mucor circinelloides query accession table is retained at:

```text
analysis/mat_locus/config/mat_query_accessions.tsv
```

The corresponding Mucor circinelloides MAT protein query FASTA is bundled at:

```text
analysis/mat_locus/queries/mucor_mat_locus_queries.faa
```

The local historical working folder is:

`/Users/hellebaalsrud/Documents/Projects/Mucormycota/MATlocus`

The pipeline assumes the cluster has a mirror of the project data folder with
genome inputs under:

`data/genomes/<genome_id-or-public-id>/`

No assemblies, GFFs, proteomes, or BLAST databases should be committed to Git.

## Recommended Saga Layout

On Saga, use the NeLS mirror and shared BLAST databases:

```text
/cluster/work/users/hellb/Mucor/data/
/cluster/work/users/hellb/Mucor/results/blastdb/
/cluster/work/users/hellb/Mucor/results/mat_locus/
/cluster/work/users/hellb/Mucor/scripts/
```

The genome-file manifest should already exist:

```text
/cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv
```

The BLAST database manifest should already exist:

```text
/cluster/work/users/hellb/Mucor/results/blastdb/blastdb_manifest.tsv
```

Transfer this module to Saga with:

```bash
rsync -ravz \
  /Users/hellebaalsrud/.codex/worktrees/a8cb/Mucoromycota-genomes-paper/analysis/mat_locus/ \
  hellb@login-1.saga.sigma2.no:/cluster/work/users/hellb/Mucor/scripts/mat_locus/
```

## Query Preparation

Extract internal query proteins from the Saga mirror:

```bash
python /cluster/work/users/hellb/Mucor/scripts/mat_locus/scripts/extract_internal_query_proteins.py \
  --query-catalog /cluster/work/users/hellb/Mucor/scripts/mat_locus/config/mat_locus_query_catalog.tsv \
  --genome-manifest /cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv \
  --out-fasta /cluster/work/users/hellb/Mucor/results/mat_locus/queries/internal_mat_locus_queries.faa \
  --missing-out /cluster/work/users/hellb/Mucor/results/mat_locus/queries/internal_mat_locus_queries.missing.tsv
```

The core Mucor circinelloides MAT-locus query FASTA contains:

```text
ACH92694.1  Mucor circinelloides tptA plus
ACH92691.1  Mucor circinelloides tptA minus
ACH92695.1  Mucor circinelloides sexP
ACH92692.1  Mucor circinelloides sexM
ACH92696.1  Mucor circinelloides rnhA plus
ACH92693.1  Mucor circinelloides rnhA minus
```

The sex-pathway query FASTA contains:

```text
CAA98021.1  Mucor mucedo tdh/TSP1 protein accession
CAP72566.1  Mucor mucedo TSP2 protein accession
```

The final query files should be placed under:

```text
/cluster/work/users/hellb/Mucor/results/mat_locus/queries/
  mat_locus_queries.faa
  mat_locus_plus_sex_pathway_queries.faa
  sex_pathway_queries.faa
```

Build the query FASTA files with:

```bash
python /cluster/work/users/hellb/Mucor/scripts/mat_locus/scripts/build_query_fastas.py \
  --mucor-mat-fasta /cluster/work/users/hellb/Mucor/scripts/mat_locus/queries/mucor_mat_locus_queries.faa \
  --sex-pathway-fasta /cluster/work/users/hellb/Mucor/scripts/mat_locus/queries/sex_pathway_queries.faa \
  --internal-mat-fasta /cluster/work/users/hellb/Mucor/results/mat_locus/queries/internal_mat_locus_queries.faa \
  --out-mat-fasta /cluster/work/users/hellb/Mucor/results/mat_locus/queries/mat_locus_queries.faa \
  --out-all-fasta /cluster/work/users/hellb/Mucor/results/mat_locus/queries/mat_locus_plus_sex_pathway_queries.faa
```

Expected counts:

```text
mat_locus_queries.faa                    14 records
mat_locus_plus_sex_pathway_queries.faa   16 records
```

## Protein Homology Screen

Run the MAT-locus plus sex-pathway query proteins against the reusable protein
BLAST databases:

```bash
mkdir -p /cluster/work/users/hellb/Mucor/results/mat_locus/logs
sbatch /cluster/work/users/hellb/Mucor/scripts/mat_locus/scripts/run_mat_blastp_against_blastdb_array.slurm
```

The array is currently set to `1-61`, matching the 61 core/final protein
databases in:

```text
/cluster/work/users/hellb/Mucor/results/blastdb/blastdb_manifest.tsv
```

Per-genome BLAST results are written to:

```text
/cluster/work/users/hellb/Mucor/results/mat_locus/blastp/<genome_id>.mat_locus_plus_sex_pathway.blastp.tsv
```

After the array finishes, summarize all hits:

```bash
python /cluster/work/users/hellb/Mucor/scripts/mat_locus/scripts/summarize_mat_blastp_hits.py \
  --blastdb-manifest /cluster/work/users/hellb/Mucor/results/blastdb/blastdb_manifest.tsv \
  --blastp-dir /cluster/work/users/hellb/Mucor/results/mat_locus/blastp \
  --mat-query-catalog /cluster/work/users/hellb/Mucor/scripts/mat_locus/config/mat_locus_query_catalog.tsv \
  --sex-pathway-catalog /cluster/work/users/hellb/Mucor/scripts/mat_locus/config/sex_pathway_query_accessions.tsv \
  --out-all /cluster/work/users/hellb/Mucor/results/mat_locus/summary/mat_locus_plus_sex_pathway.blastp.all_hits.tsv \
  --out-best-per-query /cluster/work/users/hellb/Mucor/results/mat_locus/summary/mat_locus_plus_sex_pathway.blastp.best_per_query.tsv \
  --out-best-per-gene /cluster/work/users/hellb/Mucor/results/mat_locus/summary/mat_locus_plus_sex_pathway.blastp.best_per_gene.tsv
```

## Legacy Run

1. Build a sample manifest from the metadata table and available genome folders:

```bash
python analysis/mat_locus/scripts/build_sample_manifest.py \
  --metadata data/metadata/genome_metadata_master.csv \
  --genomes-dir "${GENOMES_DIR}" \
  --output "${MAT_WORK}/sample_manifest.tsv"
```

2. Organize inputs as symlinks:

```bash
python analysis/mat_locus/scripts/prepare_cluster_inputs.py \
  --manifest "${MAT_WORK}/sample_manifest.tsv" \
  --query-fasta "${QUERY_FASTA}" \
  --work-dir "${MAT_WORK}" \
  --mode symlink
```

3. Submit the SLURM array:

```bash
sbatch analysis/mat_locus/run_mat_locus_array.slurm
```

Adjust paths and module names in the SLURM header before submitting.

If the metadata file is not under `${PROJECT_DIR}/data/metadata/`, set
`METADATA_CSV=/path/to/genome_metadata_master.csv` when submitting.

4. Summarize best homology hits:

```bash
python analysis/mat_locus/scripts/summarize_mat_hits.py \
  --organized-samples "${MAT_WORK}/inputs/organized_samples.tsv" \
  --query-map analysis/mat_locus/config/mat_query_accessions.tsv \
  --results-dir "${MAT_WORK}/results" \
  --output "${MAT_WORK}/results/summary/mat_locus_hit_summary.tsv"
```

5. Use the GFF annotations to inspect locus context:

```bash
python analysis/mat_locus/scripts/summarize_mat_locus_context.py \
  --organized-samples "${MAT_WORK}/inputs/organized_samples.tsv" \
  --blast-summary "${MAT_WORK}/results/summary/mat_locus_hit_summary.tsv" \
  --output-prefix "${MAT_WORK}/results/summary/mat_locus_context"
```

This step asks whether the best HMG candidate maps back to an annotated gene,
then reports whether the best `tptA`, HMG, and `rnhA` candidates map to the
same scaffold. It also writes the neighboring genes around each mapped HMG
candidate so the local MAT-like structure can be checked manually.

## Notes

- The first pass is a homology screen, not a final MAT locus call.
- Best hits should be checked manually for locus structure: HMG gene plus
  neighboring `tptA` and `rnhA`, strand/orientation, and whether the candidate
  is truncated or displaced.
- Broad HMG orthogroups are not enough to call a MAT locus. The GFF context
  check is the important filter after homology search.
- The existing MAT workbook contains useful manual notes and can be merged into
  the supplementary sex/ecology curation table later.
