# Rates of gene-order decay per clade and tests of the hypothesis that
# within-chromosome gene order decays faster in clades that lost CENP-A
# (Umbelopsidales, Mucorales) and CENP-A + CENP-C (Mucorales) than in
# Mortierellales.
#
# Model: an index normalised to 1 = conserved and 0 = random order decays with
# divergence d as a stretched exponential
#   I(d) = c * exp(-(lambda * d)^beta),
# where c (<= 1) is the ceiling set by annotation/assembly noise, estimated
# from near-identical genome pairs (d <= d_min_rate), lambda is a
# clade-specific rate and beta a shape shared by all clades (beta = 1 is plain
# exponential decay; see DECISIONS.md D12). beta is estimated per index and
# divergence measure from all pairs in range. Per pair,
#   lambda_ij = (-log(I_ij / c))^(1 / beta) / d_ij.
# Under random inversions, the strand index decays as exp(-2 * mu) where mu is
# the expected number of inversions covering a gene.
#
# Inference: pairs share branches of the phylogeny, so tests use random
# maximal sets of phylogenetically independent pairs (paths share no branch;
# Maddison 2000). In each draw, per-pair rates are compared between clades
# with one-sided Wilcoxon rank-sum tests. Results are summarised over draws.
#
# Output
#   results/tables/noise_ceiling.tsv
#   results/tables/decay_fits.tsv
#   results/tables/rates_by_clade.tsv
#   results/tables/rate_tests.tsv
#   results/tables/matched_values_by_clade.tsv
#   results/tables/matched_tests.tsv
#   results/cache/independent_draws.rds

source("analysis/synteny/R/functions.R")
cfg <- syn_config()
set.seed(cfg$seed)

pi_tab <- read_tsv(file.path(cfg$tables, "pair_indices.tsv"), show_col_types = FALSE)

# Indices tested, by scale. Values are normalised to 1 = conserved order and
# 0 = random order within chromosomes (or, for the macro index, 1 = no
# neighbour moved to another chromosome).
indices <- tribble(
  ~index,              ~scale,  ~label,
  "strand_index_norm", "meso",  "Strand randomisation index",
  "tau_index_norm",    "meso",  "Kendall's tau index",
  "corr_index_norm",   "meso",  "Correlation index",
  "goc_intra_norm",    "micro", "GOC (vicinity 4), intra-chromosomal",
  "goc1_intra_norm",   "micro", "Conserved adjacencies, intra-chromosomal",
  "same_chr",          "macro", "Neighbours on the same chromosome"
)
primary_index <- "strand_index_norm"

divergences <- tribble(
  ~divergence, ~div_label,
  "d_sub",     "CDS substitutions/site (ML tree)",
  "d_prot",    "Protein distance (OrthoFinder species tree)",
  "d_time",    "Time (2 x relative node age, dated tree)"
)

within <- pi_tab %>%
  filter(within) %>%
  mutate(clade = factor(clade, levels = clade_levels),
         same_chr = 1 - inter_break,
         d_time = 2 * age_rel)

long <- within %>%
  select(a, b, clade, both_ebpnor, d_sub, d_prot, d_time, all_of(indices$index)) %>%
  pivot_longer(all_of(indices$index), names_to = "index", values_to = "value")

# --- Noise ceiling c from near-identical pairs --------------------------------
ceiling_tab <- long %>%
  filter(d_sub <= cfg$d_min_rate) %>%
  group_by(index) %>%
  summarise(n_pairs = n(), c = min(1, median(value)), .groups = "drop")
write_tsv(ceiling_tab, file.path(cfg$tables, "noise_ceiling.tsv"))

# Pairs used for rate comparisons: within-order, d_min_rate < d_sub <= d_max
# (d_max ~ the deepest split within Mortierellales and Umbelopsidales, so all
# clades contribute pairs over the same divergence range).
rate_pairs <- within %>%
  filter(d_sub > cfg$d_min_rate, d_sub <= cfg$d_max) %>%
  select(a, b, clade, both_ebpnor, d_sub)

floor_value <- 0.01  # normalised values below this are treated as saturated

in_range <- long %>%
  semi_join(rate_pairs, by = c("a", "b")) %>%
  left_join(select(ceiling_tab, index, c), by = "index") %>%
  pivot_longer(all_of(divergences$divergence), names_to = "divergence", values_to = "d")

# --- Decay fits on all pairs in range (non-independent; for curves and beta) --
fit_decay <- function(df) {
  df <- filter(df, !is.na(value), d > 0)
  k <- nlevels(df$clade)
  scale <- median(df$d)  # fit on d / median(d) for numerical stability
  df$x <- df$d / scale
  f <- try(nls(value ~ c * exp(-(lambda[as.integer(clade)] * x)^beta), data = df,
               start = list(lambda = rep(1, k), beta = 1), algorithm = "port",
               lower = c(rep(1e-6, k), 0.2), upper = c(rep(1e3, k), 5)), silent = TRUE)
  if (inherits(f, "try-error")) return(tibble())
  co <- summary(f)$coefficients
  tibble(clade = levels(df$clade),
         lambda = co[seq_len(k), 1] / scale,
         lambda_se_naive = co[seq_len(k), 2] / scale,
         beta = co[k + 1, 1],
         beta_se_naive = co[k + 1, 2],
         half_life = log(2)^(1 / beta) / lambda,
         n_pairs = as.integer(table(df$clade)))
}
decay_fits <- in_range %>%
  group_by(index, divergence) %>%
  group_modify(~ fit_decay(.x)) %>%
  ungroup() %>%
  left_join(select(ceiling_tab, index, c), by = "index")
write_tsv(decay_fits, file.path(cfg$tables, "decay_fits.tsv"))

# Per-pair rates under the stretched-exponential model (D12) and, as a
# sensitivity, the plain exponential model (beta = 1) specified first (D8).
pair_rates <- in_range %>%
  left_join(distinct(decay_fits, index, divergence, beta), by = c("index", "divergence")) %>%
  crossing(model = c("stretched", "exponential")) %>%
  mutate(beta = if_else(model == "exponential", 1, beta),
         saturated = value / c < floor_value,
         rel = pmin(1, pmax(value / c, floor_value)),
         rate = (-log(rel))^(1 / beta) / d)

# --- Phylogenetically independent pairs ---------------------------------------
tree <- read.tree(cfg$tree_ml)
draw_sets <- function(pairs_df, n_draws) {
  paths <- tip_paths(tree, pairs_df)
  map(seq_len(n_draws), function(i) {
    unlist(map(split(seq_len(nrow(pairs_df)), pairs_df$clade), function(rows) {
      rows[sample_independent_pairs(paths[rows], nrow(tree$edge))]
    }), use.names = FALSE)
  })
}

# Matched-divergence window for the model-free comparison: closely related
# pairs, where all three clades have pairs and no index is saturated.
d_matched <- 0.15
genome_sets <- list(
  all = rate_pairs,
  ebpnor = filter(rate_pairs, both_ebpnor),
  matched = filter(rate_pairs, d_sub <= d_matched)
)
draws <- map(genome_sets, draw_sets, n_draws = cfg$n_draws)
saveRDS(list(pairs = genome_sets, draws = draws), file.path(cfg$cache, "independent_draws.rds"))

comparisons <- tribble(
  ~comparison,                          ~faster,          ~slower,
  "Mucorales > Mortierellales",         "Mucorales",      "Mortierellales",
  "Umbelopsidales > Mortierellales",    "Umbelopsidales", "Mortierellales",
  "Mucorales > Umbelopsidales",         "Mucorales",      "Umbelopsidales"
)

one_sided_p <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, alternative = "greater", exact = FALSE)$p.value)
}

# For each draw of independent pairs: per-clade medians of `stat` and
# one-sided Wilcoxon tests that `score` (larger = more scrambled) is larger
# in the faster clade.
run_draws <- function(df, set_name, keys_by) {
  pairs_df <- genome_sets[[set_name]]
  keys <- pairs_df %>% transmute(a, b, row = row_number())
  df <- semi_join(df, pairs_df, by = c("a", "b"))
  imap_dfr(draws[[set_name]], function(rows, draw) {
    df %>%
      inner_join(filter(keys, row %in% rows), by = c("a", "b")) %>%
      group_by(across(all_of(keys_by))) %>%
      group_modify(function(x, key) {
        med <- x %>% group_by(clade) %>%
          summarise(n = n(), median_stat = median(stat), n_saturated = sum(saturated),
                    .groups = "drop")
        tests <- comparisons %>%
          mutate(p = map2_dbl(faster, slower, ~ one_sided_p(x$score[x$clade == .x],
                                                           x$score[x$clade == .y])))
        tibble(med = list(med), tests = list(tests))
      }) %>%
      ungroup() %>%
      mutate(draw = draw, genome_set = set_name)
  })
}

summarise_draws <- function(res, keys_by) {
  meds <- res %>%
    select(genome_set, all_of(keys_by), draw, med) %>%
    unnest(med) %>%
    group_by(genome_set, across(all_of(keys_by)), clade) %>%
    summarise(pairs_per_draw = median(n),
              lo = quantile(median_stat, 0.025),
              hi = quantile(median_stat, 0.975),
              median = median(median_stat),
              saturated_per_draw = median(n_saturated),
              .groups = "drop")
  tests <- res %>%
    select(genome_set, all_of(keys_by), draw, tests) %>%
    unnest(tests) %>%
    group_by(genome_set, across(all_of(keys_by)), comparison) %>%
    summarise(median_p = median(p, na.rm = TRUE),
              frac_draws_p05 = mean(p < 0.05, na.rm = TRUE),
              n_draws = sum(!is.na(p)),
              .groups = "drop") %>%
    left_join(indices, by = "index")
  list(meds = meds, tests = tests)
}

# Rate comparisons (all genomes, EBP-Nor genomes)
rate_keys <- c("index", "divergence", "model")
rate_res <- map_dfr(c("all", "ebpnor"), ~ run_draws(
  mutate(pair_rates, stat = rate, score = rate), .x, rate_keys))
rate_sum <- summarise_draws(rate_res, rate_keys)

rates_by_clade <- rate_sum$meds %>%
  rename(median_rate = median, rate_lo = lo, rate_hi = hi)
write_tsv(rates_by_clade, file.path(cfg$tables, "rates_by_clade.tsv"))

rate_tests <- rate_sum$tests %>%
  mutate(primary = index == primary_index & divergence == "d_sub" &
           genome_set == "all" & model == "stretched")
write_tsv(rate_tests, file.path(cfg$tables, "rate_tests.tsv"))

# Model-free comparison at matched divergence: index values (lower = more
# scrambled) of independent pairs with d_min_rate < d_sub <= d_matched.
matched_values <- long %>%
  semi_join(genome_sets$matched, by = c("a", "b")) %>%
  mutate(stat = value, score = -value, saturated = FALSE)
matched_sum <- summarise_draws(run_draws(matched_values, "matched", "index"), "index")
write_tsv(matched_sum$meds %>% rename(median_value = median, value_lo = lo, value_hi = hi),
          file.path(cfg$tables, "matched_values_by_clade.tsv"))
write_tsv(matched_sum$tests, file.path(cfg$tables, "matched_tests.tsv"))

message("Primary test (", primary_index, ", d_sub, all genomes, stretched model):")
print(rate_tests %>% filter(primary) %>% select(comparison, median_p, frac_draws_p05))
message("Matched divergence (", cfg$d_min_rate, " < d <= ", d_matched, "):")
print(matched_sum$tests %>% filter(index == primary_index) %>%
        select(comparison, median_p, frac_draws_p05))
