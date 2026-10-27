# Literature search log: synteny analysis

Search performed 2026-10-27 (Claude Code, web search + full-text reading where
open access). Goal: find published concepts and methods for quantifying
within-chromosome gene-order scrambling ("mesosynteny") and for relating it to
centromere/kinetochore biology, before choosing a method.

Legend for **Use**: *adopted* = method or index used in this analysis;
*context* = cited for background/interpretation; *considered* = evaluated but
not used (reason given).

## Search queries

| # | Query (web search) | Main hits |
|---|---|---|
| 1 | (direct) Hane et al. 2011 Genome Biology, full text via PMC3219968 | Hane 2011 |
| 2 | Mucor circinelloides lacks centromeric histone CENP-A mosaic point regional centromeres | Navarro-Mendoza 2019 |
| 3 | mesosynteny quantification method Dothideomycetes intrachromosomal inversions gene order | Hane 2011 and follow-ups |
| 4 | loss of CenH3 holocentricity genome rearrangement gene order synteny | Drinnenberg 2014; Luzula 2024 |
| 5 | Fischer Rocha Dujon 2006 highly variable rates of genome rearrangements hemiascomycetous yeast lineages | Fischer 2006 |
| 6 | Coghlan Wolfe fourfold faster rate of genome rearrangement nematodes Drosophila | Coghlan & Wolfe 2002 |
| 7 | odp oxford dot plot software Schultz ctenophore fusion with mixing ALG irreversible | Schultz 2023 (via 8) |
| 8 | Mucorales centromere evolution CENP-C loss kinetochore Mucoromycota Umbelopsis Mortierella (extended) | Navarro-Mendoza 2019, 2023 |
| 9 | Schultz 2023 Nature ancient gene linkages ... Oxford dot plot odp mixing | Schultz 2023 (methods read in full) |
| 10 | Merian elements Lepidoptera chromosome evolution holocentric intrachromosomal rearrangement Wright 2024 | Wright 2024 |
| 11 | Nigon elements nematode chromosome evolution gene order shuffled within chromosomes holocentric Oscheius tipulae | Gonzalez de la Rosa 2021 |
| 12 | Bhutkar 2008 chromosomal rearrangement 12 Drosophila genomes Muller elements synteny blocks | Bhutkar 2008 |
| 13 | synteny index SI gene order conservation measure Shifman Ninyo Gophna Snir Phylo SI | Shifman 2014 |
| 14 | Huynen Bork 1998 measuring genome evolution gene order conservation versus sequence identity | Huynen & Bork 1998 |
| 15 | GENESPACE Lovell 2022 eLife synteny orthology | Lovell 2022 |
| 16 | mesosynteny chromosome-level fungal genomes 2022 2023 2024 ... Pezizomycotina comparative (extended) | Hane 2011; Shen 2020; Cutaneotrichosporon 2023 |
| 17 | holocentric chromosomes higher rate intrachromosomal rearrangements synteny decay Rhynchospora Carex Luzula | Carex 2024; Rhynchospora; Luzula |
| 18 | centromere loss inversion rate fungi intercentromeric recombination Cryptococcus Candida Sanyal | Coelho 2024; Guin 2020 |
| 19 | Navarro-Mendoza kinetochore survey Mucoromycotina 55 species CENP-A CENP-C | Navarro-Mendoza 2019 |
| 20 | Shen 2020 Science Advances contrasting modes genome evolution Ascomycota gene order | Shen 2020 (full text not accessible) |
| 21 | Plessy 2024 Genome Research extreme genome scrambling Oikopleura ... GenomicBreaks | Plessy 2024 (preprint read in full) |
| 22 | d'Alencon 2010 PNAS synteny holocentric Lepidoptera local rearrangements | d'Alencon 2010 |
| 23 | syntenet R Bioconductor package synteny network MCScanX | Almeida-Silva 2023 |
| 24 | Drinnenberg 2014 eLife recurrent loss of CenH3 holocentricity insects | Drinnenberg 2014 |
| 25 | GenomicBreaks reference pages (GOC, synteny_index, correlation_index, strand_randomisation_index, tau_index) and source code | GenomicBreaks docs |
| 26 | macrosyntR CRAN oxford grid orthologs Fisher test | macrosyntR |
| 27 | Ranz 2022 Proc Biol Sci randomization limits chromosomal gene organization Lepidoptera Diptera | Ranz 2022 |
| 28 | Mucorales chromosome-level genome assembly synteny ... Rhizopus Mucor lusitanicus | R. stolonifer 2025; Lichtheimia 2014 |
| 29 | Mortierella chromosome-level genome synteny Mortierellomycotina centromere CENP-A | nothing Mortierella-specific; Malassezia 2020 |
| 30 | van Hooff Tromer Snel Kops 2017 EMBO Reports kinetochore network CENP-A CENP-C | van Hooff 2017 |

No published study was found that relates kinetochore/CENP gene loss in fungi
to the rate of within-chromosome gene-order scrambling. Searches 8, 19, 28 and
29 found no Mortierellales or Umbelopsidales centromere studies and no
comparative synteny study across Mucoromycota orders.

## Sources

### Mesosynteny and gene-order scrambling: concepts and methods

- **Hane JK, Rouxel T, Howlett BJ, Kema GHJ, Goodwin SB, Oliver RP (2011).**
  A novel mode of chromosomal evolution peculiar to filamentous Ascomycete
  fungi. *Genome Biology* 12:R45.
  [link](https://link.springer.com/article/10.1186/gb-2011-12-5-r45),
  [PMC3219968](https://pmc.ncbi.nlm.nih.gov/articles/PMC3219968).
  Defines mesosynteny as conservation of gene content within homologous
  chromosomes with randomised gene order and orientation. Method: PROMER
  (MUMmer) dot plots on scaffolds >= 500 kb; binomial test of conserved
  sequence per scaffold pair; sliding-window detection of collinear regions;
  genome pairs with mean collinear-region length >= 20 kb called
  macrosyntenic, otherwise mesosyntenic. Proposes frequent intrachromosomal
  inversions with very rare translocations; no centromere discussion; no
  model relating decay to divergence.
  **Use:** *adopted* (concept, the macro- vs meso-scale distinction, and
  orientation randomisation as a hallmark); the alignment-based thresholds
  were not reused because we work with orthologue anchors and need a
  continuous index to compare rates.

- **Plessy C, Mansfield MJ, Bliznina A, et al. (2024).** Extreme genome
  scrambling in marine planktonic *Oikopleura dioica* cryptic species.
  *Genome Research* 34:426-440. doi:10.1101/gr.278295.123;
  [preprint](https://www.biorxiv.org/content/10.1101/2023.05.09.539028.full).
  Scrambling quantified with whole-genome alignments (LAST) and the R package
  GenomicBreaks: strand randomisation index, collinear vs isolated
  alignments, breakpoint accumulation rate (breakpoints/Mb/Myr) compared with
  *Ciona* and *Drosophila*. Chromosome arms remain the unit of macrosynteny.
  **Use:** *adopted* (strand randomisation index as primary index, synteny
  and correlation indices, collinear-block definition, rate-per-divergence
  framing).

- **GenomicBreaks R package** (Plessy lab, OIST).
  [GitHub](https://github.com/oist/GenomicBreaks),
  [reference](https://oist.github.io/GenomicBreaks/reference/index.html).
  Documentation and source of `synteny_index`, `correlation_index`
  (weighted Pearson of positions on the best-matching chromosome),
  `strand_randomisation_index` (`|P - M| / (P + M)`, weighted mean),
  `GOC` (Rocha 2003; vicinity 4; recommended for orthologue positions) and
  `tau_index` (Yanai 2005 / Ranz 2022). **Use:** *adopted*, re-implemented
  for orthologue anchors in `R/functions.R` (the package was not installed;
  see DECISIONS.md D6).

- **Rocha EPC (2003).** DNA repeats lead to the accelerated loss of gene order
  in bacteria. *Trends in Genetics* 19:600-603. doi:10.1016/j.tig.2003.09.011.
  Origin of the gene order conservation (GOC) index as implemented in
  GenomicBreaks. **Use:** *adopted* (GOC; micro-scale index).

- **Schultz DT, Haddock SHD, Bredeson JV, Green RE, Simakov O, Rokhsar DS
  (2023).** Ancient gene linkages support ctenophores as sister to other
  animals. *Nature* 618:110-117. doi:10.1038/s41586-023-05936-6;
  [PMC10232365](https://pmc.ncbi.nlm.nih.gov/articles/PMC10232365).
  Oxford dot plots of reciprocal-best-hit or OrthoFinder orthologues;
  Bonferroni-corrected one-sided Fisher's exact tests for chromosome-pair
  associations; "mixing" of fused linkage groups scored by transitions
  between groups against random orderings. Software: odp
  ([GitHub](https://github.com/conchoecia/odp)). **Use:** *adopted* (Oxford
  dot plots; Fisher test for the macrosynteny fraction); *considered* the
  mixing score (designed for fusions of ancestral linkage groups, not for
  order within one homologous chromosome pair).

- **El Hilali S, Copley R. macrosyntR** (CRAN, 2023). Draws ordered Oxford
  grids from orthologue tables; Fisher tests of chromosome associations;
  community-based ordering. [CRAN](https://cran.r-project.org/web/packages/macrosyntR/refman/macrosyntR.html).
  **Use:** *adopted* (same Fisher-test logic; ordering of dot-plot axes by
  main partner chromosome); package itself not used (our dot plots also
  colour gene orientation).

- **Ranz JM, Gonzalez PM, Su RN, Bedford SJ, Abreu-Goodger C, Markow T
  (2022).** Multiscale analysis of the randomization limits of the
  chromosomal gene organization between Lepidoptera and Diptera.
  *Proc R Soc B* 289:20212183. [PMC8767184](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8767184/).
  Holocentric (monarch) vs monocentric (*Drosophila*) genomes; multiscale
  indices including the tau index; chromosome repatterning slower in
  Lepidoptera than Diptera but same order of magnitude. **Use:** *context*
  (multiscale framing; tau index mentioned); full text blocked by a captcha.

- **Huynen MA, Bork P (1998).** Measuring genome evolution. *PNAS*
  95:5849-5856. [PubMed](https://pubmed.ncbi.nlm.nih.gov/9600883/).
  Gene order changes faster than gene content and protein sequence; genome
  distances related to divergence. **Use:** *context* (index vs sequence
  divergence framework).

- **Shifman A, Ninyo N, Gophna U, Snir S (2014).** Phylo SI: a new
  genome-wide approach for prokaryotic phylogeny. *Nucleic Acids Research*.
  [PMC3936750](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3936750/).
  Synteny index (shared genes in k-neighbourhoods). **Use:** *considered*;
  equivalent in spirit to GOC, which we use.

### Rearrangement rates relative to divergence

- **Fischer G, Rocha EPC, Brunet F, Vergassola M, Dujon B (2006).** Highly
  variable rates of genome rearrangements between hemiascomycetous yeast
  lineages. *PLoS Genetics* 2:e32.
  [PMC1391921](https://pmc.ncbi.nlm.nih.gov/articles/PMC1391921/).
  Rearrangement rates estimated per branch, correlated within but highly
  variable between lineages. **Use:** *context* (precedent for comparing
  lineage-specific rearrangement rates, normalised by sequence divergence).

- **Coghlan A, Wolfe KH (2002).** Fourfold faster rate of genome
  rearrangement in nematodes than in *Drosophila*. *Genome Research*
  12:857-867. doi:10.1101/gr.172702.
  0.4-1.0 breaks/Mb/Myr in *Caenorhabditis* (holocentric) vs *Drosophila*.
  **Use:** *context* (rate comparison between lineages; holocentric lineage
  with fast intrachromosomal rearrangement).

- **Bhutkar A, Schaeffer SW, Russo SM, et al. (2008).** Chromosomal
  rearrangement inferred from comparisons of 12 *Drosophila* genomes.
  *Genetics* 179:1657-1680.
  ~95% of genes stay on the same Muller element; paracentric inversions
  shuffle order within elements; scrambling roughly linear with distance.
  **Use:** *context* (meso-like pattern in a monocentric, CENP-A-bearing
  lineage: within-chromosome scrambling is not by itself diagnostic of
  centromere change).

- **Shen X-X, Steenwyk JL, LaBella AL, et al. (2020).** Genome-scale
  phylogeny and contrasting modes of genome evolution in the fungal phylum
  Ascomycota. *Science Advances* 6:eabd0079. Full text not accessible
  (HTTP 403); search snippets report higher gene-order divergence in
  Pezizomycotina than Saccharomycotina. **Use:** *context* only.

- **Cutaneotrichosporon chromosome-level assemblies (2023).** Chromosome-level
  genome assemblies of *Cutaneotrichosporon* spp. reveal imbalanced evolution
  between nucleotide sequences and chromosome synteny. *BMC Genomics*.
  [PMC10568926](https://pmc.ncbi.nlm.nih.gov/articles/PMC10568926/).
  Synteny disruption (number of locally collinear blocks) much higher than
  expected from ANI, compared descriptively with *Saccharomyces*; authors
  speculate holocentric or very short centromeres. **Use:** *context*
  (closest fungal precedent for "more rearrangement than expected for
  sequence divergence" linked to centromere type; descriptive only).

### Centromeres and kinetochores

- **Navarro-Mendoza MI, Perez-Arques C, Panchal S, et al. (2019).** Early
  diverging fungus *Mucor circinelloides* lacks centromeric histone CENP-A and
  displays a mosaic of point and regional centromeres. *Current Biology*
  29:3791-3802. [PubMed 31679929](https://pubmed.ncbi.nlm.nih.gov/31679929/),
  [PMC6925572](https://pmc.ncbi.nlm.nih.gov/articles/PMC6925572).
  CENP-A and CENP-C absent; kinetochore (Mis12, Dsn1) localises to a short
  AT-rich region with a conserved motif flanked by Grem-LINE1 retroelements;
  survey of 55 Mucoromycotina proteomes; CENP-A and CENP-C lost in the
  ancestor of Mucorales. **Use:** *context* (the centromere type in the
  scrambled clade).

- **Heterochromatin and RNAi act independently to ensure genome stability in
  Mucorales human fungal pathogens (2023).** *PNAS* 120:e2220475120 (authors
  not verified in this search).
  [link](https://www.pnas.org/doi/10.1073/pnas.2220475120). **Use:** *context*
  (genome stability mechanisms in Mucorales).

- **van Hooff JJE, Tromer E, van Wijk LM, Snel B, Kops GJPL (2017).**
  Evolutionary dynamics of the kinetochore network in eukaryotes as revealed
  by comparative genomics. *EMBO Reports* 18:1559-1571.
  [PMC5579357](https://pmc.ncbi.nlm.nih.gov/articles/PMC5579357).
  Reported the absence of CENP-A and CENP-C in Mucoromycotina species.
  **Use:** *context*.

- **Drinnenberg IA, deYoung D, Henikoff S, Malik HS (2014).** Recurrent loss
  of CenH3 is associated with independent transitions to holocentricity in
  insects. *eLife* 3:e03676. [PMC4359364](https://pmc.ncbi.nlm.nih.gov/articles/PMC4359364/).
  **Use:** *context* (CenH3 loss and centromere architecture).

- **Coelho MA, et al. (2024).** Comparative genomics of *Cryptococcus* and
  *Kwoniella* (*PLoS Biology*; title paraphrased): karyotype change by
  fusions with pericentric inversions, centromere loss, and
  intercentromeric recombination mediated by shared centromeric transposons.
  [PMC11185503](https://pmc.ncbi.nlm.nih.gov/articles/PMC11185503).
  **Use:** *context* (centromeres as sites and constraints of fungal
  rearrangements).

- **Guin K, et al. (2020).** Spatial inter-centromeric interactions
  facilitated the emergence of evolutionary new centromeres. *eLife*
  9:e58556. **Use:** *context*.

- **Sankaranarayanan SR, et al. (2020).** Loss of centromere function drives
  karyotype evolution in closely related *Malassezia* species. *eLife*
  9:e53944. **Use:** *context*.

### Holocentric lineages (CENP-A-independent or diffuse centromeres)

- **d'Alencon E, Sezutsu H, Legeai F, et al. (2010).** Extensive synteny
  conservation of holocentric chromosomes in Lepidoptera despite high rates
  of local genome rearrangements. *PNAS* 107:7680-7685.
  Conserved chromosome-level synteny but only 1.3-2.0 genes per conserved
  block; ~2 breaks/Mb/Myr, faster than *Drosophila*, possibly linked to
  holocentricity. **Use:** *context* (closest animal analogue of
  mesosynteny; "genes per collinear block" index adopted as a secondary
  index).

- **Wright CJ, Stevens L, Mackintosh A, Lawniczak M, Blaxter M (2024).**
  Comparative genomics reveals the dynamics of chromosome evolution in
  Lepidoptera. *Nature Ecology & Evolution*.
  [PMC11009112](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11009112/).
  32 Merian elements largely intact over ~250 Myr despite holocentricity.
  **Use:** *context*.

- **Gonzalez de la Rosa PM, Thomson M, Trivedi U, et al. (2021).** A
  telomere-to-telomere assembly of *Oscheius tipulae* and the evolution of
  rhabditid nematode chromosomes. *G3* 11:jkaa020.
  [PMC8022731](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8022731/).
  Seven Nigon elements; karyotype stable but gene order within elements not
  conserved (holocentric nematodes). **Use:** *context*.

- Holocentric plants: *Carex* (Genomic hotspots of chromosome rearrangements
  explain conserved synteny despite high rates of chromosome evolution in a
  holocentric lineage, 2024; [PMC11628656](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11628656/));
  *Luzula sylvatica* (2024; [PMC11538461](https://pmc.ncbi.nlm.nih.gov/articles/PMC11538461)).
  Mainly fission/fusion; **Use:** *context*.

### Synteny software considered

- **Lovell JT, Sreedasyam A, Schranz ME, et al. (2022).** GENESPACE tracks
  regions of interest and gene copy number variation across multiple
  genomes. *eLife* 11:e78526. **Use:** *considered* (preferred by the
  project, but no Mucoromycota GENESPACE run exists under
  `Mucor/data` on NeLS; see DECISIONS.md D1).

- **Almeida-Silva F, Zhao T, Ullrich KK, Schranz ME, Van de Peer Y (2023).**
  syntenet: an R/Bioconductor package for the inference and analysis of
  synteny networks. *Bioinformatics* 39:btac806. MCScanX port in Rcpp.
  **Use:** *considered* (collinear-block calling is not needed for the
  indices; blocks are defined directly from anchor adjacency).

### Background not found in this search (cited from prior knowledge)

- Maddison WP (2000). Testing character correlation using pairwise
  comparisons on a phylogeny. *J Theor Biol* 202:195-204. Basis for using
  phylogenetically independent (non-overlapping) pairs. **Use:** *adopted*.
- Emms DM, Kelly S (2019). OrthoFinder: phylogenetic orthology inference for
  comparative genomics. *Genome Biology* 20:238. Source of the pairwise
  orthologues. **Use:** *adopted*.
