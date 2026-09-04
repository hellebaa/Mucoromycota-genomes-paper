# Assembly structure supporting figure

This analysis summarizes chromosome-level assembly structure for the
Mucoromycota genome set.

The first pass produces one row per contig/scaffold, with contig length and
whether gene annotation is present on that contig. This supports a draft
supplementary figure with one mini bar plot per genome: one bar per contig,
colored by gene-bearing status.

Later layers can be added to the same contig-level table, including telomere
hits, centromere candidates, repeat density, GC, and gene-density windows.

Current optional evidence layers:

- Telomeres: direct scans for terminal arrays of the canonical fungal/animal
  telomere motif `TTAGGG` and its reverse complement `CCCTAA`.
- Candidate centromeres: BLAST hits to Mucor centromere-repeat/GremLINE query
  sequences. The current reference coordinates come from the MucCirc3 BED file
  made from published Mucor centromere motifs and GremLINE annotations.

## Inputs

Run on the Saga mirror of the NeLS data directory:

```bash
/cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv
```

The manifest should contain at least:

- `genome_id`
- `species_name`
- `include_core_analyses`
- `assembly_fasta`
- `gff`

## Outputs

Default output directory:

```bash
/cluster/work/users/hellb/Mucor/results/assembly_structure
```

Main files:

- `assembly_contig_gene_summary.tsv`: one row per contig
- `assembly_genome_summary.tsv`: one row per genome
- `telomere_hits.bed`: terminal telomere motif arrays
- `telomere_contig_summary.tsv`: one row per contig with left/right telomere calls
- `centromere_blast_hits.bed`: filtered centromere/GremLINE BLAST hits
- `centromere_contig_summary.tsv`: one row per contig with candidate-centromere evidence
- `assembly_contig_bars.pdf`: draft supporting figure
- `assembly_contig_bars.png`: PNG copy of the draft supporting figure

## Saga run

```bash
sbatch /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/summarize_assembly_structure.slurm
```

After the TSVs are generated, they can be downloaded locally and plotted, or the
R plotting script can be run on Saga if the required R packages are available.

## Telomere scan

```bash
sbatch /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/run_telomere_scan.slurm
```

The default scan checks the first and last 10 kb of each contig and records
terminal motif arrays with at least three consecutive motif copies.

For a broader exploratory telomere pass, run the multi-motif scan:

```bash
sbatch /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/run_telomere_scan_multi_motif.slurm
```

This uses `config/telomere_motifs.tsv`, checks the first and last 25 kb of each
contig, and allows hits up to 5 kb from the contig end. The output files are:

- `telomere_hits_multi_motif.bed`
- `telomere_contig_summary_multi_motif.tsv`

Treat missing telomere hits as "not detected with the tested motif/settings",
not as evidence that telomeres are absent.

## Candidate centromere BLAST summary

First prepare query sequences from the MucCirc3 BED coordinates:

```bash
sbatch /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/prepare_centromere_queries.slurm
```

This extracts:

- all BED intervals
- the 9 `CEN_*` core intervals
- all GremLINE intervals

It also scans only the extracted `CEN_*` sequences for the candidate signature
motif `GCATAGGT`.

Then BLAST the 9 CEN-core queries against the core assemblies:

```bash
sbatch /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/run_centromere_blastn_array.slurm
```

The array size is set to 61 for the current core dataset. If the number of core
genomes changes, check it with:

```bash
awk -F'\t' 'NR > 1 && $5 == "yes" && $10 != "" {n++} END {print n}' \
  /cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv
```

Use this after BLASTing Mucor centromere/GremLINE query sequences against the
assemblies to summarize hits:

```bash
python /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/summarize_centromere_blast_hits.py \
  --blast-dir /cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blastn_cen_core \
  --pattern '*.MucCirc3_CEN_core.blastn.tsv' \
  --file-suffix '.MucCirc3_CEN_core.blastn.tsv' \
  --out-bed /cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits.bed \
  --out-summary /cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_contig_summary.tsv
```

This is an evidence table, not a final centromere caller. Final centromere
coordinates should combine motif hits, GremLINE/repeat evidence, gene-poor
regions, and, where available, Hi-C/3D genome information.

To scan the BLAST-hit intervals for the conserved `GCATAGGT` signature motif:

```bash
python /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/scan_centromere_blast_hit_motifs.py \
  --manifest /cluster/work/users/hellb/Mucor/data/metadata/genome_file_manifest.tsv \
  --centromere-bed /cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits.bed \
  --motif GCATAGGT \
  --out-hits /cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits_motif_scan.tsv \
  --out-retained-bed /cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits_GCATAGGT_retained.bed \
  --out-contig-summary /cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_motif_retained_contig_summary.tsv \
  --out-genome-summary /cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_motif_retained_genome_summary.tsv
```

These outputs distinguish all CEN-core BLAST hits from the subset that retains
the short motif inside the aligned candidate interval.

## Rough karyotype summary

After generating assembly and telomere summaries:

```bash
python /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/summarize_rough_karyotype.py \
  --contigs /cluster/work/users/hellb/Mucor/results/assembly_structure/assembly_contig_gene_summary.tsv \
  --telomeres /cluster/work/users/hellb/Mucor/results/assembly_structure/telomere_contig_summary_multi_motif.tsv \
  --out /cluster/work/users/hellb/Mucor/results/assembly_structure/rough_karyotype_summary.tsv
```

The default rough estimate counts gene-bearing contigs at least 0.5 Mb long and
reports telomere support and the first sharp length drop between ranked contigs.
This is a draft assembly-interpretation aid, not a final karyotype call.

## GremLINE context

GremLINE evidence is treated as a pericentromeric/context track, not as a
standalone centromere call. First cluster the extracted MucCirc3 GremLINE
intervals to reduce redundancy:

```bash
sbatch /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/cluster_gremline_queries.slurm
```

Then BLAST the clustered representative GremLINE sequences against the core
assemblies:

```bash
sbatch /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/run_gremline_blastn_array.slurm
```

Summarize GremLINE BLAST hits:

```bash
python /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/summarize_centromere_blast_hits.py \
  --blast-dir /cluster/work/users/hellb/Mucor/results/assembly_structure/gremline_blastn_cdhit90 \
  --pattern '*.MucCirc3_GremLINE_cdhit90.blastn.tsv' \
  --file-suffix '.MucCirc3_GremLINE_cdhit90.blastn.tsv' \
  --out-bed /cluster/work/users/hellb/Mucor/results/assembly_structure/gremline_blast_hits.bed \
  --out-summary /cluster/work/users/hellb/Mucor/results/assembly_structure/gremline_contig_summary.tsv \
  --min-pident 80 \
  --min-length 500 \
  --max-evalue 1e-20
```

Finally, summarize whether CEN-core BLAST hits have nearby GremLINE evidence on
the same contig:

```bash
python /cluster/work/users/hellb/Mucor/scripts/assembly_structure/scripts/summarize_cen_gremline_proximity.py \
  --cen-bed /cluster/work/users/hellb/Mucor/results/assembly_structure/centromere_blast_hits.bed \
  --gremline-bed /cluster/work/users/hellb/Mucor/results/assembly_structure/gremline_blast_hits.bed \
  --max-distance-bp 100000 \
  --out-cen-summary /cluster/work/users/hellb/Mucor/results/assembly_structure/cen_gremline_proximity_cen_hit_summary.tsv \
  --out-contig-summary /cluster/work/users/hellb/Mucor/results/assembly_structure/cen_gremline_proximity_contig_summary.tsv
```
