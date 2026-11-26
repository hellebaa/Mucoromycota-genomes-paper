# Plan: signatures of genome duplication in Mucoromycota

## Question

Corrochano et al. 2016 (Curr Biol 26:1577) proposed a whole-genome duplication
(WGD) early in the Mucoromycotina. Their phylome had four Mucorales genomes
(*Lichtheimia*, *Phycomyces*, *Mucor*, *Rhizopus*) with *Mortierella alpina* as
the closest outgroup, so the proposed WGD lies after the split from
Mortierellales: somewhere between the Umbelopsidales + Mucorales ancestor
(node N9 in our OrthoFinder tree) and the Mucorales ancestor (N12).

Their evidence is a burst of gene-tree duplications on that branch (0.70-0.96
duplications per gene) and larger gene families than in Dikarya. They found few
duplicated collinear regions (4-13 per genome), no 1:2 conserved synteny with an
outgroup and no Ks peak, and did not compare a WGD with elevated small-scale
duplication. We test the hypothesis with far denser sampling (42 Mucorales,
6 Umbelopsidales, 13 Mortierellales) and Hi-C assemblies in each order.

Inputs: doubletrouble duplicate pairs (SD segmental, TD tandem, PD proximal,
DD dispersed) with Ka/Ks per genome; OrthoFinder orthogroups, gene trees and
duplication events; dated species tree; genome metadata.

## 1. Genome set and QC

- Harmonise genome IDs across doubletrouble, OrthoFinder, the tree and metadata.
- Contig count, N50 and gene-bearing contigs for every genome, computed the same
  way for all (the QC table covers only part of our own assemblies).
  Segmental-duplicate detection depends on contiguity.
- Flag likely artifacts: gzMucFlav4 (2,762 SD pairs at median Ks ~0.01:
  duplicated haplotypes?), strain pairs of one species, ggAmyRoux1 (86.8% BUSCO
  duplicated; excluded from the tree).

## 2. Duplicate landscape per species

- Pairs by type (SD/TD/PD/DD), scaled by gene count, on the dated tree.
- Do genomes with more genes have more duplicates, and of which type?

## 3. Ks age distributions (recent events only)

- Pairs with 0.05 < Ks < 2 (lower bound removes allelic/assembly duplicates).
  An event at N9/N12 is beyond Ks saturation, so Ks cannot date it.
- Peaks per species and type (mixture models on log Ks).
- Ortholog Ks between sister species, to place peaks before or after speciation
  (e.g. shared or separate events in *Absidia* and *Rhizopus*; *R. delemar*
  WGD, Ma et al. 2009).

## 4. Testing a WGD on the Mucorales / Mucoromycotina stem

Mortierellales (and Umbelopsidales, if the event is on N12) are the
unduplicated outgroup.

a. Locate the burst: duplications per branch from OrthoFinder, scaled by branch
   length on the dated tree. Does N9 or N12 stand out among branches of
   similar length? (First pass, >=50% support: N8 86, N9 773, N12 1,264;
   Mortierellales stem N10 1,588.)
b. Collinearity of those duplicates: are duplications mapped to N9/N12 enriched
   for SD pairs in Mucorales genomes, relative to younger duplicates?
c. 1:2 conserved synteny (key test): does a region in ggMorZona1 / gzUmbVina1
   map to two regions in a Mucorales genome, with complementary gene loss?
   Start at the macro level (shared chromosomes), since local gene order erodes
   fast in Mucorales (synteny analysis).
d. Concordant retention: do the same orthogroups keep two copies across
   Mucorales lineages more often than independent small-scale duplication
   would give?
e. Later, if a-d are promising: model-based test (Whale) of WGD vs no WGD on
   the branch, accounting for gene-tree uncertainty.

## 5. Copy number with microsynteny

For gene families expanded in a clade relative to its sister, are the extra
copies in collinear blocks (WGD-like) or dispersed/tandem (small-scale)?

## 6. Order-free synteny (added)

Collinearity (doubletrouble SD: >= 5 anchors in order) is lost when gene order
is scrambled within chromosomes, as in Mucorales. Neighbourhood synteny: a
paralog pair is syntenic if paralog pairs of at least k other gene families
link the +-W gene neighbourhoods of its two copies, ignoring order (W = 10,
k = 2; W = 20/30 as checks), calibrated with shuffled gene order. Applied to
the duplicates of each gene-tree branch, within the same genomes across
branches of different age, and per duplication across Mucorales genomes.

## Separate question

Mortierellales retain 435-841 SD pairs per genome, mostly Ks > 2, and the
Mortierellales stem (N10) has the largest duplication burst. Origin unknown.
