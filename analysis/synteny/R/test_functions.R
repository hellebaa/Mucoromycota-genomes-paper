# Sanity checks of the gene-order indices on synthetic genomes with known
# structure. Run: Rscript analysis/synteny/R/test_functions.R

source("analysis/synteny/R/functions.R")
set.seed(1)

make_anchors <- function(seq_b, pos_b, rel, seq_a = rep("A1", length(pos_b))) {
  tibble(og = paste0("og", seq_along(pos_b)), gene_a = paste0("a", seq_along(pos_b)),
         gene_b = paste0("b", seq_along(pos_b)),
         seq_a = seq_a, pos_a = seq_along(pos_b), seq_b = seq_b, pos_b = pos_b, rel = rel) %>%
    add_ranks()
}
check <- function(label, cond) {
  cat(sprintf("%-60s %s\n", label, if (isTRUE(cond)) "ok" else "FAILED"))
  if (!isTRUE(cond)) quit(status = 1)
}

n <- 200
colinear <- make_anchors(rep("B1", n), 1:n, rep(1L, n))
ix <- directional_indices(colinear)
check("colinear: GOC = 1", ix$goc == 1)
check("colinear: correlation, strand, synteny = 1",
      ix$corr_index == 1 && ix$strand_index == 1 && ix$synteny_index == 1)
check("colinear: one block of n genes", ix$block_mean == n)

# One large inversion (genes 51-150 reversed): order still conserved locally
inv_pos <- c(1:50, 150:51, 151:200)
inv_rel <- c(rep(1L, 50), rep(-1L, 100), rep(1L, 50))
ix <- directional_indices(make_anchors(rep("B1", n), inv_pos, inv_rel))
check("inversion: GOC stays high (> 0.98)", ix$goc > 0.98)
check("inversion: strict adjacency breaks only at 2 ends", isTRUE(all.equal(ix$n_adj * (1 - ix$goc1), 2)))
check("inversion: strand index = |150 - 50| / 200 = 0", ix$strand_index == 0)
check("inversion: three collinear blocks", ix$block_mean == n / 3)

# Translocation: second half on another query chromosome
tr <- make_anchors(rep(c("B1", "B2"), each = n / 2), c(1:(n / 2), 1:(n / 2)), rep(1L, n))
ix <- directional_indices(tr)
check("translocation: one inter-chromosomal break", ix$inter_break == 1 / (n - 1))
check("translocation: synteny index 0.5", ix$synteny_index == 0.5)
check("translocation: intra-chromosomal GOC = 1", ix$goc_intra == 1)

# Fully scrambled order within one chromosome (mesosynteny)
meso <- make_anchors(rep("B1", n), sample(n), sample(c(-1L, 1L), n, TRUE))
ix <- directional_indices(meso)
check("mesosynteny: synteny index = 1", ix$synteny_index == 1)
check("mesosynteny: GOC near random (< 0.15)", ix$goc < 0.15)
check("mesosynteny: correlation near 0 (< 0.2)", ix$corr_index < 0.2)

# Normalised indices: colinear = 1, scrambled ~ 0
cfg <- modifyList(syn_config(), list(n_perm = 10L))
pc <- pair_indices(colinear, cfg)
pm <- pair_indices(meso, cfg)
check("normalised GOC: colinear = 1", pc$goc_intra_norm == 1)
check("normalised GOC: scrambled ~ 0 (|x| < 0.1)", abs(pm$goc_intra_norm) < 0.1)
meso2 <- make_anchors(rep(c("B1", "B2"), each = n / 2), c(sample(n / 2), sample(n / 2)),
                      sample(c(-1L, 1L), n, TRUE), seq_a = rep(c("A1", "A2"), each = n / 2))
check("macrosynteny: scrambled chromosomes still significant", macrosynteny(meso2)$macro_frac == 1)
# Fast base-R directional_indices() against a straightforward dplyr reference
reference_indices <- function(an, v = 4L, min_pair_n = 5L) {
  adj <- an %>%
    arrange(seq_a, rank_a) %>%
    group_by(seq_a) %>%
    mutate(next_seq_b = lead(seq_b), next_rank_b = lead(rank_b)) %>%
    ungroup() %>%
    filter(!is.na(next_seq_b)) %>%
    mutate(same_chr = next_seq_b == seq_b, step = next_rank_b - rank_b)
  near <- adj$same_chr & abs(adj$step) <= v
  adj1 <- adj$same_chr & abs(adj$step) == 1L
  tab <- count(an, seq_a, seq_b)
  best <- tab %>% group_by(seq_a) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup()
  bp <- an %>%
    semi_join(filter(best, n >= min_pair_n), by = c("seq_a", "seq_b")) %>%
    group_by(seq_a, seq_b) %>%
    summarise(n = n(), r = abs(cor(pos_a, pos_b)),
              tau = abs(cor(rank_a, rank_b, method = "kendall")),
              s = abs(sum(rel)) / n(), .groups = "drop")
  blocks <- an %>%
    arrange(seq_a, rank_a) %>%
    group_by(seq_a) %>%
    mutate(linked = coalesce(seq_b == lag(seq_b) & rank_b - lag(rank_b) == rel &
                               lag(rel) == rel, FALSE),
           block = cumsum(!linked)) %>%
    count(seq_a, block, name = "size")
  tibble(n_adj = length(near), goc = mean(near), goc1 = mean(adj1),
         inter_break = mean(!adj$same_chr), goc_intra = mean(near[adj$same_chr]),
         goc1_intra = mean(adj1[adj$same_chr]), synteny_index = sum(best$n) / sum(tab$n),
         corr_index = weighted.mean(bp$r, bp$n), tau_index = weighted.mean(bp$tau, bp$n),
         strand_index = weighted.mean(bp$s, bp$n), block_mean = mean(blocks$size),
         frac_block3 = sum(blocks$size[blocks$size >= 3]) / sum(blocks$size))
}
same_as_reference <- function(an) {
  isTRUE(all.equal(as.data.frame(directional_indices(an)), as.data.frame(reference_indices(an))))
}
check("fast = reference: inversion", same_as_reference(make_anchors(rep("B1", n), inv_pos, inv_rel)))
check("fast = reference: two scrambled chromosomes", same_as_reference(meso2))
rnd <- make_anchors(sample(c("B1", "B2", "B3"), 600, TRUE), runif(600) * 1e6,
                    sample(c(-1L, 1L), 600, TRUE), seq_a = rep(c("A1", "A2", "A3"), each = 200))
check("fast = reference: random three-chromosome genome", same_as_reference(rnd))
cfg <- syn_config()
if (file.exists(orthologue_file("ggMorZona1.proteins", "gzLinHyal1.proteins", cfg))) {
  real <- build_anchors(load_genome("ggMorZona1", cfg), load_genome("gzLinHyal1", cfg),
                        "ggMorZona1.proteins", "gzLinHyal1.proteins", cfg)
  check("fast = reference: M. zonata vs L. hyalina", same_as_reference(real))
  check("fast = reference: same pair, query as target", same_as_reference(swap_anchors(real)))
}
cat("All tests passed\n")
