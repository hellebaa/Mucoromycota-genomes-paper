# Shared functions for the synteny / mesosynteny analysis.
#
# All scripts are run from the repository root (see run_all.sh). Paths are
# taken from syn_config().
#
# Terminology
#   anchor   one-to-one orthologue pair (OrthoFinder pairwise orthologues)
#   target   genome "a" of a pair (rows are ordered along its sequences)
#   query    genome "b" of a pair
#   rank     index of an anchor among the anchors on the same sequence
#
# Gene-order indices are re-implementations, for orthologue anchors, of
# published indices (see analysis/synteny/LITERATURE.md):
#   GOC                Rocha 2003 Trends Genet; GenomicBreaks::GOC (Plessy et al. 2024)
#   synteny index      GenomicBreaks::synteny_index (Plessy et al. 2024)
#   correlation index  GenomicBreaks::correlation_index (Plessy et al. 2024)
#   strand index       GenomicBreaks::strand_randomisation_index (Plessy et al. 2024)
#   macrosynteny       Fisher test of chromosome pairs, as in odp (Schultz et al.
#                      2023) and macrosyntR (El Hilali & Copley 2023)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
})

syn_config <- function() {
  repo <- normalizePath(Sys.getenv("REPO_ROOT", "."), mustWork = TRUE)
  if (!dir.exists(file.path(repo, "analysis", "synteny"))) {
    stop("Run from the repository root or set REPO_ROOT")
  }
  out <- file.path(repo, "analysis", "synteny", "results")
  list(
    repo       = repo,
    data       = file.path(repo, "data"),
    genomes    = file.path(repo, "analysis", "synteny", "genomes.tsv"),
    out        = out,
    tables     = file.path(out, "tables"),
    figures    = file.path(out, "figures"),
    cache      = file.path(out, "cache"),
    tree_ml    = file.path(repo, "data", "phylogeny", "sco_cds_withOutgroups_b500.treefile"),
    tree_dated = file.path(repo, "data", "phylogeny", "sco_cds_b500.date.nwk"),
    tree_prot  = file.path(repo, "data", "orthofinder", "results", "Species_Tree",
                           "SpeciesTree_rooted_node_labels.txt"),
    orthologues = file.path(repo, "data", "orthofinder", "results", "Orthologues"),
    vicinity   = 4L,     # GOC neighbourhood (GenomicBreaks default)
    min_pair_n = 5L,     # min anchors on a chromosome pair for correlation/strand indices
    n_perm     = 20L,    # permutations for the random-order expectation
    d_max      = 0.65,   # max subs/site for between-clade rate comparisons
    d_min_rate = 0.01,   # pairs closer than this only inform the annotation-noise ceiling
    n_draws    = 1000L,  # random sets of phylogenetically independent pairs
    seed       = 20261027L
  )
}

clade_levels <- c("Mortierellales", "Umbelopsidales", "Mucorales")
clade_colours <- c(Mortierellales = "#3B6FB6", Umbelopsidales = "#C27A0E", Mucorales = "#A23B3B")

read_genome_table <- function(cfg = syn_config(), all = FALSE) {
  read_tsv(cfg$genomes, show_col_types = FALSE) %>%
    filter(all | include == "yes") %>%
    mutate(order = factor(order, levels = clade_levels))
}

# Genome pairs analysed: all pairs within an order, plus all pairs of EBP-Nor
# genomes across orders. Genome a precedes genome b in genomes.tsv, which also
# fixes which OrthoFinder orthologue table (a__v__b) is used.
analysis_pairs <- function(genomes) {
  idx <- combn(nrow(genomes), 2)
  g <- function(col, i) genomes[[col]][idx[i, ]]
  tibble(a = g("genome", 1), b = g("genome", 2),
         of_a = g("of_name", 1), of_b = g("of_name", 2),
         order_a = g("order", 1), order_b = g("order", 2),
         source_a = g("source", 1), source_b = g("source", 2)) %>%
    filter(order_a == order_b | (source_a == "EBPNOR" & source_b == "EBPNOR"))
}

# ---------------------------------------------------------------------------
# Gene positions and ID aliases
# ---------------------------------------------------------------------------

read_gene_positions <- function(genome, cfg = syn_config()) {
  f <- file.path(cfg$data, "genomes", genome, paste0(genome, "_genePosTbl.tsv"))
  read_tsv(f, col_types = "ciicc") %>%
    transmute(gene = geneID, seq = seqname, start, end, strand,
              pos = (start + end) / 2) %>%
    distinct(gene, .keep_all = TRUE)
}

# Lookup table from any identifier OrthoFinder may have used (gene, protein or
# transcript ID) to the gene ID used in the position table.
read_gene_aliases <- function(genome, genes, cfg = syn_config()) {
  f <- file.path(cfg$data, "genomes", genome, paste0(genome, "_IDtbl.tsv"))
  ids <- read_tsv(f, col_types = "ccc")
  bind_rows(
    tibble(id = genes$gene, gene = genes$gene),
    tibble(id = ids$proteinID, gene = ids$geneID),
    tibble(id = ids$rnaID, gene = ids$geneID)
  ) %>%
    filter(!is.na(id), gene %in% genes$gene) %>%
    distinct(id, .keep_all = TRUE)
}

load_genome <- function(genome, cfg = syn_config()) {
  genes <- read_gene_positions(genome, cfg)
  list(genes = genes, alias = read_gene_aliases(genome, genes, cfg))
}

# ---------------------------------------------------------------------------
# One-to-one orthologue anchors
# ---------------------------------------------------------------------------

orthologue_file <- function(of_a, of_b, cfg = syn_config()) {
  file.path(cfg$orthologues, paste0("Orthologues_", of_a),
            paste0(of_a, "__v__", of_b, ".tsv"))
}

# OrthoFinder separates genes with ", ", but public FASTA headers can carry
# descriptions that themselves contain ", ". Count only the elements whose
# first word is a known identifier of that genome.
parse_gene_cell <- function(x, known_ids) {
  parts <- str_split(x, ", ")
  first <- map(parts, ~ word(.x, 1))
  hits <- map(first, ~ .x[.x %in% known_ids])
  tibble(n = lengths(hits), id = map_chr(hits, ~ if (length(.x) == 1) .x else NA_character_))
}

# Returns one row per one-to-one orthologue with positions in both genomes.
build_anchors <- function(ga, gb, of_a, of_b, cfg = syn_config()) {
  raw <- read_tsv(orthologue_file(of_a, of_b, cfg), col_types = "ccc",
                  col_names = c("og", "a", "b"), skip = 1)
  pa <- parse_gene_cell(raw$a, ga$alias$id)
  pb <- parse_gene_cell(raw$b, gb$alias$id)
  one2one <- tibble(og = raw$og, id_a = pa$id, id_b = pb$id) %>%
    filter(pa$n == 1, pb$n == 1)
  one2one %>%
    inner_join(ga$alias, by = c(id_a = "id")) %>% rename(gene_a = gene) %>%
    inner_join(gb$alias, by = c(id_b = "id")) %>% rename(gene_b = gene) %>%
    distinct(gene_a, .keep_all = TRUE) %>% distinct(gene_b, .keep_all = TRUE) %>%
    inner_join(ga$genes %>% select(gene, seq_a = seq, pos_a = pos, strand_a = strand),
               by = c(gene_a = "gene")) %>%
    inner_join(gb$genes %>% select(gene, seq_b = seq, pos_b = pos, strand_b = strand),
               by = c(gene_b = "gene")) %>%
    add_ranks() %>%
    mutate(rel = if_else(strand_a == strand_b, 1L, -1L)) %>%
    select(og, gene_a, gene_b, seq_a, pos_a, rank_a, seq_b, pos_b, rank_b, rel)
}

add_ranks <- function(an) {
  an %>%
    group_by(seq_a) %>% mutate(rank_a = rank(pos_a, ties.method = "first")) %>%
    group_by(seq_b) %>% mutate(rank_b = rank(pos_b, ties.method = "first")) %>%
    ungroup()
}

swap_anchors <- function(an) {
  an %>% rename(gene_a = gene_b, gene_b = gene_a, seq_a = seq_b, seq_b = seq_a,
                pos_a = pos_b, pos_b = pos_a, rank_a = rank_b, rank_b = rank_a)
}

# ---------------------------------------------------------------------------
# Gene-order indices (directional: target a, query b)
# ---------------------------------------------------------------------------

# Written in base R on plain vectors because it runs ~40 times per pair
# (2 directions x (1 + n_perm shuffles)); R/test_functions.R checks it against
# a dplyr reference implementation.
#
# Adjacency: each anchor and the next anchor along the same target sequence
# (anchors at the end of a target sequence have no follower and are dropped,
# as in GenomicBreaks::GOC).
# Collinear block: maximal run of anchors that follow each other in both
# genomes with consistent orientation ("collinear alignments" of Plessy et al.
# 2024 / GenomicBreaks::flagColinearAlignments, applied to genes).
directional_indices <- function(an, v = 4L, min_pair_n = 5L) {
  o <- order(an$seq_a, an$rank_a)
  sa <- an$seq_a[o]; sb <- an$seq_b[o]; rel <- an$rel[o]
  ra <- an$rank_a[o]; rb <- an$rank_b[o]; pa <- an$pos_a[o]; pb <- an$pos_b[o]
  n <- length(sa)

  follows <- sa[-1] == sa[-n]                 # anchor i+1 follows anchor i on the target
  same_b <- sb[-1] == sb[-n]
  step <- rb[-1] - rb[-n]
  intra <- same_b[follows]
  near <- (same_b & abs(step) <= v)[follows]
  adj1 <- (same_b & abs(step) == 1L)[follows]

  linked <- c(FALSE, follows & same_b & step == rel[-1] & rel[-n] == rel[-1])
  block_size <- tabulate(cumsum(!linked))

  pair_n <- table(sa, sb)
  best_n <- apply(pair_n, 1, max)
  best_b <- setNames(colnames(pair_n)[max.col(pair_n, ties.method = "first")], rownames(pair_n))
  on_best <- sb == best_b[sa] & best_n[sa] >= min_pair_n
  bp <- vapply(split(which(on_best), sa[on_best]), function(i) {
    c(length(i), abs(cor(pa[i], pb[i])), abs(cor(ra[i], rb[i], method = "kendall")),
      abs(sum(rel[i])) / length(i))
  }, numeric(4))

  tibble(
    n_adj         = length(near),
    goc           = mean(near),                    # GOC, vicinity v (Rocha 2003)
    goc1          = mean(adj1),                    # strict conserved adjacency
    inter_break   = mean(!intra),                  # neighbour moved to another chromosome
    goc_intra     = mean(near[intra]),             # GOC among intra-chromosomal adjacencies
    goc1_intra    = mean(adj1[intra]),
    synteny_index = sum(best_n) / n,               # GenomicBreaks::synteny_index
    corr_index    = weighted.mean(bp[2, ], bp[1, ]),  # GenomicBreaks::correlation_index (|r|)
    tau_index     = weighted.mean(bp[3, ], bp[1, ]),  # |Kendall tau|, rank-based variant
    strand_index  = weighted.mean(bp[4, ], bp[1, ]),  # GenomicBreaks::strand_randomisation_index
    block_mean    = mean(block_size),              # genes per collinear block (d'Alencon 2010)
    frac_block3   = sum(block_size[block_size >= 3]) / n
  )
}

# Random gene order within chromosomes: shuffle anchor positions (and ranks,
# jointly) within each query sequence, keeping chromosome content, and
# randomise relative orientation.
shuffle_within_query <- function(an) {
  idx <- seq_len(nrow(an))
  for (i in split(idx, an$seq_b)) idx[i] <- i[sample.int(length(i))]
  an$pos_b <- an$pos_b[idx]
  an$rank_b <- an$rank_b[idx]
  an$rel <- sample(c(-1L, 1L), nrow(an), replace = TRUE)
  an
}

order_indices <- c("goc", "goc1", "goc_intra", "goc1_intra", "corr_index", "tau_index",
                   "strand_index", "block_mean", "frac_block3")

# Symmetric indices for a pair: mean of the two directions. Random-order
# expectations from n_perm shuffles, and normalised values
# (obs - random) / (1 - random), where 1 = fully conserved order and
# 0 = order as random as chromosome content allows.
pair_indices <- function(an, cfg = syn_config()) {
  both <- function(x) {
    bind_rows(directional_indices(x, cfg$vicinity, cfg$min_pair_n),
              directional_indices(swap_anchors(x), cfg$vicinity, cfg$min_pair_n)) %>%
      summarise(across(everything(), mean))
  }
  obs <- both(an)
  null <- map_dfr(seq_len(cfg$n_perm), function(i) both(shuffle_within_query(an))) %>%
    summarise(across(all_of(order_indices), mean))
  norm <- map2_dfc(obs[order_indices], null[order_indices],
                   ~ if (.y < 1) (.x - .y) / (1 - .y) else NA_real_)
  bind_cols(
    tibble(n_anchors = nrow(an),
           n_seq_a = n_distinct(an$seq_a), n_seq_b = n_distinct(an$seq_b)),
    obs,
    rename_with(null, ~ paste0(.x, "_rand")),
    rename_with(norm, ~ paste0(.x, "_norm")),
    macrosynteny(an)
  )
}

# Fraction of anchors on significantly associated chromosome pairs
# (one-sided Fisher / hypergeometric test, Bonferroni over all pairs),
# as used for Oxford dot plots by odp and macrosyntR.
macrosynteny_pairs <- function(an, alpha = 0.05) {
  N <- nrow(an)
  na <- count(an, seq_a, name = "n_a")
  nb <- count(an, seq_b, name = "n_b")
  m <- nrow(na) * nrow(nb)
  count(an, seq_a, seq_b) %>%
    left_join(na, by = "seq_a") %>% left_join(nb, by = "seq_b") %>%
    mutate(p = phyper(n - 1, n_b, N - n_b, n_a, lower.tail = FALSE),
           p_adj = pmin(1, p * m),
           sig = p_adj < alpha)
}

macrosynteny <- function(an) {
  mp <- macrosynteny_pairs(an)
  tibble(macro_frac = sum(mp$n[mp$sig]) / sum(mp$n),
         macro_pairs = sum(mp$sig))
}

# ---------------------------------------------------------------------------
# Phylogeny helpers
# ---------------------------------------------------------------------------

pair_distance_table <- function(genomes, cfg = syn_config()) {
  ml <- read.tree(cfg$tree_ml)
  dated <- read.tree(cfg$tree_dated)
  prot <- read.tree(cfg$tree_prot)
  prot$tip.label <- sub("(\\.proteins|_prot)$", "", prot$tip.label)
  d_ml <- cophenetic(ml)
  d_dt <- cophenetic(dated)
  d_pr <- cophenetic(prot)
  root_age <- max(node.depth.edgelength(dated))
  lookup <- function(d, a, b) {
    ok <- a %in% rownames(d) & b %in% rownames(d)
    out <- rep(NA_real_, length(a))
    out[ok] <- d[cbind(a[ok], b[ok])]
    out
  }
  function(a, b) {
    tibble(d_sub = lookup(d_ml, a, b),
           d_prot = lookup(d_pr, a, b),
           # node age of the MRCA as a fraction of the root age of the dated tree
           age_rel = lookup(d_dt, a, b) / 2 / root_age)
  }
}

# Edge indices on the path between two tips of a tree (unrooted is fine).
tip_paths <- function(tree, pairs) {
  g <- igraph::graph_from_edgelist(tree$edge, directed = FALSE)
  tip <- setNames(seq_along(tree$tip.label), tree$tip.label)
  map2(pairs$a, pairs$b, function(a, b) {
    as.integer(igraph::shortest_paths(g, tip[[a]], tip[[b]], output = "epath")$epath[[1]])
  })
}

# One maximal set of phylogenetically independent pairs (paths share no
# branch; Maddison 2000), by greedy selection of the shortest paths first
# (number of branches), ties broken at random. Shortest-first approximates the
# largest possible number of pairs; the random tie-breaks give different sets
# across draws.
sample_independent_pairs <- function(paths, n_edges) {
  used <- logical(n_edges)
  keep <- logical(length(paths))
  for (i in order(lengths(paths), runif(length(paths)))) {
    e <- paths[[i]]
    if (!any(used[e])) {
      keep[i] <- TRUE
      used[e] <- TRUE
    }
  }
  which(keep)
}
