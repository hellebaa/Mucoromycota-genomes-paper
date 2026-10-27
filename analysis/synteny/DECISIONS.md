# Decision log: synteny analysis

Decisions made while setting up the analysis of within-chromosome gene-order
scrambling (mesosynteny) across Mortierellales, Umbelopsidales and Mucorales.
Newest entries are added at the end. Literature behind the choices is in
[LITERATURE.md](LITERATURE.md).

## 2026-10-27

### D1. Synteny data source: OrthoFinder orthologues, not GENESPACE

- The project prefers the GENESPACE synteny data on NeLS. No Mucoromycota
  GENESPACE run exists under `Projects/NMBU_threeD_Sandve_2024/Mucor/data`
  (searched with `find` for `*genespace*`, `*synHits*`, `*synten*`,
  `*riparian*`, `*.anchors`). A GENESPACE run elsewhere in the NeLS project
  (`genomes/orthology/genespace`) is for salmonid fish and unrelated. The user
  confirmed the Mucor run may have been removed and asked to only use
  `Mucor/data`.
- Used instead: OrthoFinder 2.5.4 results already in the project
  (`data/orthofinder/results` -> `Results_May21_4`, 69 proteomes), pairwise
  orthologue tables (`Orthologues/Orthologues_<sp>/<sp>__v__<sp2>.tsv`), and
  the per-genome gene position tables (`data/genomes/<g>/<g>_genePosTbl.tsv`,
  `<g>_IDtbl.tsv`).
- Rationale: GENESPACE would add synteny-constrained orthology and MCScanX
  collinear blocks. Re-running it means re-running OrthoFinder/DIAMOND and
  MCScanX on ~60 genomes. Its block-based output is also designed around
  collinearity, which is exactly what breaks down under mesosynteny; the
  indices used here (D6) only need one-to-one orthologue positions, which is
  how odp (Schultz et al. 2023), macrosyntR and GenomicBreaks' orthologue
  indices are applied.
- If a GENESPACE run becomes available, its `syntenicHits` could replace the
  anchors in `build_anchors()`; the indices would not change.

### D2. Orthologue anchors: OrthoFinder pairwise one-to-one orthologues

- Anchors are rows of the pairwise orthologue tables with exactly one gene in
  each genome. These come from OrthoFinder gene trees, so one-to-one
  orthologues inside duplicated families are kept. That matters for
  Mucorales, where whole-genome duplication left many two-copy families that
  an orthogroup-level single-copy filter would drop.
- OrthoFinder separates genes with `", "` but public FASTA headers keep their
  descriptions, which can contain `", "`. Gene counts per cell therefore only
  count elements whose first word is a known ID of that genome.
- Orthologue IDs are mapped to coordinates via gene, protein or transcript
  ID (`*_IDtbl.tsv`). All included genomes map >= 98.8% of their OrthoFinder
  genes (results/tables/genome_qc.tsv).

### D3. Genome set

- All 61 Mucoromycota genomes in `data/genomes` with gene tables
  (Mortierellales 13, Umbelopsidales 6, Mucorales 42). The Dikarya,
  *Conidiobolus* and Glomeromycotina outgroups are not used.
- Excluded: *Mucor lusitanicus* `MucCirc3`. OrthoFinder used JGI protein IDs
  (`jgi|Muccir1_3| <proteinId> |...`) that do not occur in
  `MucCirc3_IDtbl.tsv` (which uses GFF `gene_N`/`mRNA_N` IDs). Mapping would
  need the GFF3 `proteinId` attribute and a special-case parser for one of
  42 Mucorales genomes. Recorded in `genomes.tsv` (`include = no`).
- Assembly contiguity is high for all genomes: 78-100% of genes sit on
  sequences carrying >= 100 genes (lowest: `gzMucFlav4` 78%, `PhaArti` 87%).
  No contiguity filter. All within-chromosome indices are computed per
  homologous sequence pair, so fragmentation mostly reduces the number of
  anchors per pair rather than creating apparent scrambling.
- Public genomes (annotation source "Other") and EBP-Nor genomes (one
  annotation pipeline) are mixed in the main analysis to get enough
  phylogenetically independent pairs in Mortierellales and Umbelopsidales.
  Sensitivity analysis: EBP-Nor pairs only.

### D4. Genome pairs

- All pairs within each order (Mortierellales 78, Umbelopsidales 15,
  Mucorales 820 after excluding MucCirc3) plus all EBP-Nor pairs across
  orders (for completeness; not used in the tests). 1,158 pairs in total.
- One orthologue table per pair (first genome in `genomes.tsv` order vs the
  second); indices are symmetrised (D6).

### D5. Divergence measure

- Primary: patristic distance on the maximum-likelihood tree of concatenated
  single-copy orthologue CDS (`data/phylogeny/sco_cds_withOutgroups_b500.treefile`,
  substitutions/site). This is the "genetic distance" the hypothesis is
  framed in.
- Sensitivity: protein distance on the OrthoFinder species tree; time as
  2 x the relative node age of the MRCA on the dated tree
  (`sco_cds_b500.date.nwk`, root age scaled to 1).
- The preliminary figure's distances (e.g. *M. fragilis* vs *M. plumbeus*
  0.050) do not match the current tree: in the CDS tree that pair is 0.318
  subs/site and *M. zonata* vs *L. hyalina* is 0.040. The preliminary
  figure used an older tree with different taxa (e.g. *A. rouxii*,
  *D. globulifera* inside *Mucor*). Matched-divergence comparisons are
  redone with the current tree; the closest *Mucor* pairs at the same
  distance as *M. zonata* vs *L. hyalina* are strains of *M. circinelloides*
  (0.026) and *M. fragilis* (0.071).

### D6. Indices of gene order (what "scrambling" is measured with)

All indices are re-implementations for orthologue anchors of published
indices (GenomicBreaks; Plessy et al. 2024; Rocha 2003; odp/macrosyntR).
GenomicBreaks was not installed: it is a GitHub-only package built for
whole-genome alignment objects, and the definitions are short enough to
re-implement transparently. `R/test_functions.R` checks each index on
synthetic genomes (collinear, one inversion, one translocation, fully
scrambled). Gene order is the rank of an anchor among anchors on the same
sequence.

| Scale | Index | Definition |
|---|---|---|
| macro | neighbour on same chromosome | fraction of adjacent anchor pairs (target) whose orthologues are on the same query sequence (1 - inter-chromosomal breaks) |
| macro | synteny index | anchors on the best-matching query sequence, summed over target sequences / all anchors |
| macro | macrosynteny fraction | anchors on sequence pairs with Bonferroni-significant one-sided Fisher test (odp, macrosyntR) |
| meso | strand randomisation index | per target sequence and its best query sequence, `|n_same - n_opposite| / n` of relative gene orientation, weighted mean |
| meso | correlation index | `|Pearson r|` of anchor positions on target sequence vs best query sequence, weighted mean |
| meso | Kendall's tau index | `|Kendall tau|` of anchor ranks, as above |
| micro | GOC | fraction of adjacencies whose orthologues are within 4 anchors on the query (Rocha 2003; GenomicBreaks default vicinity) |
| micro | conserved adjacency | as GOC with vicinity 1 |
| micro | genes per collinear block | mean length of runs of anchors adjacent in both genomes with consistent orientation |

- Indices are computed in both directions and averaged.
- Deviation from GenomicBreaks: correlation and strand indices are computed
  per pair of homologous sequences (target sequence and its best-matching
  query sequence) and use absolute values, because sequence orientation in
  an assembly is arbitrary.
- Random expectation: anchor positions are shuffled within each query
  sequence (chromosome content kept) and relative orientation is
  randomised; 20 permutations. Within-chromosome indices are normalised as
  `(observed - random) / (1 - random)`: 1 = conserved order, 0 = as random
  as chromosome content allows.

### D7. Primary index chosen after a pilot (disclosed)

- The initial plan was GOC as primary index, because it is the standard
  orthologue-based gene-order index and links to breakpoint-rate studies.
- A pilot on six pairs (before any rate comparison was run) showed that GOC
  measures local neighbourhoods, which stay intact in all clades, while the
  dot-plot pattern in the hypothesis (dots spread within chromosome boxes)
  is a loss of order and orientation at the scale of whole chromosomes:

  | pair | d | GOC intra (norm.) | strand index | correlation index |
  |---|---|---|---|---|
  | M. zonata vs L. hyalina | 0.040 | 0.99 | 0.96 | 1.00 |
  | M. circinelloides strains | 0.026 | 0.98 | 0.29 | 0.87 |
  | M. fragilis strains | 0.071 | 0.96 | 0.18 | 0.69 |
  | U. isabellina strains | 0.124 | 1.00 | 0.89 | 1.00 |
  | M. fragilis vs M. plumbeus | 0.318 | 0.87 | 0.12 | 0.24 |

- Primary index: the strand randomisation index. Plessy et al. (2024)
  designed it for genomes "where gene order changed a lot but genes stayed
  on the same chromosome", i.e. mesosynteny, and Hane et al. (2011) define
  mesosynteny by randomised order *and orientation*. It also has a simple
  expectation under random inversions: the fraction of genes with
  conserved relative orientation decays as `exp(-2 mu)`, `mu` being the
  number of inversions covering a gene. Secondary meso-scale indices:
  Kendall's tau and correlation index. Micro (GOC) and macro (same
  chromosome) indices are reported as well, so the conclusion can be checked
  at every scale.
- Risk: choosing the index after seeing a few values can favour the
  hypothesis. Mitigation: all indices are reported with the same tests, the
  pilot table above is kept, and the six pilot pairs are only a small part
  of the 913 within-order pairs.

### D8. Statistical model and test

- Decay model per index: `I(d) = c * exp(-lambda * d)`, with the ceiling
  `c` (annotation/assembly noise) estimated as the median index of
  near-identical pairs (`d <= 0.01` subs/site; 14 pairs, mostly *M.
  plumbeus* strains). Per-pair rate `lambda_ij = -log(I_ij / c) / d_ij`;
  normalised values below 0.01 are set to 0.01 and counted as saturated
  (rate is then a lower bound).
- Pairs used for rates: within-order, `0.01 < d <= 0.65` subs/site. 0.65
  is just above the deepest split in Mortierellales (0.63) and
  Umbelopsidales (0.59), so all clades are compared over the same range;
  Mucorales pairs beyond 0.65 are shown in figures but not tested.
- Non-independence: pairs share branches. Tests use random maximal sets of
  phylogenetically independent pairs (paths in the ML tree share no
  branch; Maddison 2000), drawn by random-order greedy selection within each
  clade; 1,000 draws. In each draw, per-pair rates are compared with
  one-sided Wilcoxon rank-sum tests: Mucorales > Mortierellales,
  Umbelopsidales > Mortierellales, Mucorales > Umbelopsidales. Reported:
  median p over draws and the fraction of draws with p < 0.05, and the
  per-clade median rate (median and 95% range over draws).
- Descriptive curves: non-linear least squares fit of the decay model with
  clade-specific `lambda` on all pairs in range (standard errors ignore
  non-independence and are not used for inference).
- The primary test is pre-specified as: strand index, CDS distance, all
  genomes, Mucorales > Mortierellales and Umbelopsidales > Mortierellales.

### D9. Output locations

- The task allowed creating files only in `analysis/synteny/` (plus
  downloads into `data/`). All outputs therefore go to
  `analysis/synteny/results/` (`tables/`, `figures/`) rather than
  `data/synteny/`; per-pair caches go to `results/cache/`, which is
  git-ignored.

### D10. Download script

- `download_data.sh` builds its file list from `genomes.tsv`: gene tables
  for included genomes and one orthologue table per analysed pair (1,330
  files, ~0.7 GB). It still calls only `scripts/download_files.sh`.
- `RSYNC_RSH` is set to reuse one SSH connection (ControlMaster) for all
  files; this cut the download from roughly one SSH login per file to one
  in total. The control socket is in a short `mktemp -d /tmp/cm.XXXXXX`
  directory because Unix socket paths are limited to ~104 characters. It
  does not reference NeLS and is ignored by a non-SSH downloader.

### D11. Interim look at partial results (disclosed)

- While indices were being computed, clade medians for ~260 finished pairs
  were inspected (after D7 fixed the primary index, before any test was
  run). They showed that (i) relative gene orientation randomises in both
  Mortierellales and Mucorales by d ~ 0.4 subs/site, faster in Mucorales,
  whereas (ii) global order (Kendall's tau) stays at ~0.5-0.6 in deep
  Mortierellales pairs but drops to ~0.1 in Mucorales, and (iii) the strand
  index does not decay at a constant exponential rate in Mortierellales.
- The primary index was not changed. Kendall's tau is reported as the main
  secondary index because it corresponds most directly to the dot-plot
  pattern (positions, not orientation). The decay model form is revisited
  in D12 on the full data, using the shape of the decay only (not the clade
  difference).

### D12. Decay model: stretched exponential with a shared shape

- D8 assumed plain exponential decay. The interim look (D11) showed the
  strand index decaying slowly at first and then quickly in Mortierellales
  (0.96 at d = 0.04, ~0.7 at d ~ 0.15, ~0.05 at d ~ 0.45), which a constant
  rate cannot fit; per-pair rates would then depend on d, and clades with
  different divergence ranges among their pairs would not be comparable.
- Model changed to `I(d) = c * exp(-(lambda * d)^beta)`: clade-specific rate
  `lambda`, one shape `beta` shared by all clades, estimated per index and
  divergence measure by non-linear least squares on all pairs in range.
  `beta = 1` is the exponential model; `beta > 1` describes decay that
  accelerates, as expected when an index first tolerates rearrangements
  (e.g. small inversions that leave chromosome-scale order intact) and then
  breaks down. A shared shape cannot by itself create a difference between
  clades. Per-pair rates: `lambda_ij = (-log(I_ij / c))^(1 / beta) / d_ij`.
- The fitted curves are drawn over the data in the figures so the fit can be
  judged per clade.
- Implementation note: `directional_indices()` was rewritten in base R for
  speed (about 2x faster per pair single-threaded) after the first full run
  was started;
  `R/test_functions.R` checks that it reproduces the dplyr reference
  implementation exactly on synthetic genomes and on *M. zonata* vs
  *L. hyalina*. All pairs were recomputed with the new code.

### D13. Independent-pair selection: shortest paths first

- D8 drew independent pairs by greedy selection in uniformly random order.
  The first run showed that this gives few pairs per draw (median 4
  Mortierellales, 1-2 Umbelopsidales, 9 Mucorales), because a deep pair
  picked early blocks many branches. With 2 Umbelopsidales pairs a one-sided
  Wilcoxon test against 4 Mortierellales pairs cannot reach p < 0.05
  (minimum p = 1/15), so the Umbelopsidales comparisons were untestable.
- Changed to greedy selection of the shortest paths (fewest branches) first,
  ties broken at random, which approximates the maximal number of
  independent pairs that Maddison (2000) recommends. Pairs per draw (all
  genomes): 6 Mortierellales, 3 Umbelopsidales, 11-12 Mucorales (the maximum
  possible is 6 and 3 for the first two). Shortest-first favours closely
  related pairs; divergence is accounted for by the per-pair rate.
- Effect on the primary test (strand index, CDS distance, all genomes):
  Mucorales > Mortierellales median p 0.025 -> 0.0004; Umbelopsidales >
  Mortierellales median p 0.88 -> 0.74 (no evidence either way of the
  change). The direction of every comparison was the same under both
  selection rules. This change was made after seeing the first results and
  is disclosed here for that reason.
- Also fixed: the 95% range of per-clade median rates was computed after
  the median had overwritten the column (`summarise()` evaluates
  sequentially), which collapsed the interval to a point.

### D14. Model-free matched-divergence test and descriptive smoothers

- After the first figures, the shared-shape decay model visibly misfits
  parts of the data (Mucorales strand index levels off at ~0.05-0.15 instead
  of 0, because genes move in blocks that keep their orientation; the
  Mortierellales tau index stays ~1 among close pairs). To avoid relying on
  the model, a model-free test was added: raw normalised index values of
  phylogenetically independent pairs with 0.01 < d <= 0.15 (all three clades
  have pairs there; nothing is saturated), one-sided Wilcoxon tests, 1,000
  draws (shortest paths first, as in D13).
- The plain exponential model specified first (D8) is reported as a
  sensitivity analysis next to the stretched model.
- Main-figure curves are descriptive LOESS smoothers (degree 1, span 1)
  instead of the model fits; the model fits are in Supplementary Figure S5.
- Result: same conclusions as the rate-based tests (strand index, tau,
  correlation: Mucorales vs Mortierellales p = 0.004; Umbelopsidales vs
  Mortierellales p = 0.76; the same-chromosome index shows no
  Mucorales-Mortierellales difference, p = 0.60).

### D15. Genome table generated from the phylogeny; one rsync call per download

- `genomes.tsv` was first assembled by a one-off script that was not kept in
  the repo. It is now written by `R/00_genome_table.R` from the dated tree
  (`sco_cds_b500.date.nwk`, whose 61 tips are exactly the Mucoromycota
  genomes):
  - order = all tips below a node given as the MRCA of two tips
    (Mortierellales: gzPodHuml1 + ggMorZona1; Umbelopsidales: gzUmbIsab1 +
    gzUmbRama1; Mucorales: SynRace + gzMucPlum1);
  - CENP-A absent below MRCA(gzUmbIsab1, gzMucPlum1), CENP-C absent below
    the Mucorales node (the losses as described for this project; not read
    from a gene presence/absence table);
  - excluded genomes listed by name with the reason (MucCirc3, D3);
  - species, family and source joined from `genome_metadata_master.csv`;
    OrthoFinder species name matched in the `Orthogroups.tsv` header.
  The script checks that every tip gets an order, that the tree-derived
  order agrees with the metadata, and that every genome has an OrthoFinder
  name. The regenerated file is byte-identical to the earlier one, so pair
  directions and all results are unchanged.
- The pair rule (within-order pairs plus EBP-Nor pairs across orders) now
  lives in one function, `analysis_pairs()`, used both for the download list
  and by `01_pairs_and_qc.R` (previously duplicated in awk and R).
- `scripts/download_files.sh` now passes all paths to a single
  `rsync --files-from` call instead of one rsync per file. `-k` makes
  symlinked directories on NeLS (`data/orthofinder/results`) arrive as real
  directories. This replaces the SSH connection-reuse workaround in D10.
  Checking all 1,287 files takes ~4 s when they are already present.
