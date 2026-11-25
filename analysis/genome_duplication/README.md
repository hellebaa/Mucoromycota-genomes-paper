# Genome duplication analysis

Segmental and whole-genome duplications in 61 Mucoromycota genomes
(Mortierellales, Umbelopsidales, Mucorales), from doubletrouble duplicate pairs
and OrthoFinder gene trees.

## Files

| File | Purpose |
|---|---|
| `wgd_paper.Rmd` (+ `.html`) | Analyses for the paper: duplicate classes, differences across genomes, segmental duplicates by lineage, ancestral segmental duplications (order-free neighbourhood synteny), Ks of segmental and ancestral-segmental duplicates |
| `wgd_mucoromycota.Rmd` (+ `.html`) | Extended analyses and checks behind the paper report: assembly QC and near-identical pairs, gene-tree duplication rates, macrosynteny, copy-number asymmetry, expanded families, robustness of the synteny test |
| `PLAN.md` | Analysis plan |

## Running

Knit either Rmd from RStudio, or from the repo root:

```bash
Rscript -e 'rmarkdown::render("analysis/genome_duplication/wgd_paper.Rmd")'
```

Input files are fetched from NeLS on the first run through
`scripts/download_files.sh`. Three slow steps are cached in
`data/genome_duplication/` and reused afterwards:

- `contig_lengths.tsv`: assembly sequence lengths, computed on NeLS
- `ortholog_ks.tsv`: Ks of one-to-one orthologs (MSA2dist)
- `neighbourhood_support.tsv.gz`: neighbourhood synteny support of every paralog pair

R packages: data.table, ape, nlme, ggplot2, ggtree, aplot, ggrepel, mclust,
patchwork, Biostrings and doubletrouble (for MSA2dist).
