# Exploration: why does the strand index drop while Kendall's tau barely moves
# in closely related Mucor pairs? Dot plots and runs of inverted genes for
# selected pairs. Not part of run_all.sh; run after 02_pair_indices.R.
#
# Per chromosome pair (each target sequence and its main partner, as for the
# indices), the raw strand index is 1 - 2 * (fraction of genes in the minority
# orientation): every inverted gene counts the same, whatever the size of the
# inversion it sits in. Kendall's tau counts gene pairs whose order is
# reversed: a clean inversion of k genes on a chromosome pair with n genes
# reverses k(k - 1)/2 of n(n - 1)/2 gene pairs, so tau falls by about
# 2 * (k / n)^2. Small inversions therefore lower the strand index but hardly
# change tau.
#
# An inverted run is a maximal run of consecutive genes (along genome a) in the
# minority orientation of their chromosome pair. span_b is the number of genes
# the run spans in genome b; span_b == n_genes means the run is also
# contiguous there, as expected for a simple inversion.
#
# Output
#   results/figures/explore_strand_tau_dotplots.{pdf,png}
#   results/figures/explore_strand_tau_runs.{pdf,png}
#   results/tables/explore_inverted_runs.tsv

source("analysis/synteny/R/functions.R")
source("analysis/synteny/R/plot_functions.R")
cfg <- syn_config()
genomes <- read_genome_table(cfg)
pi_tab <- read_tsv(file.path(cfg$tables, "pair_indices.tsv"), show_col_types = FALSE)

pairs <- tribble(
  ~a,           ~b,           ~title,
  "gzMucPlum3", "gzMucPlum6", "Mucor, near-identical",
  "gzMucPlum3", "gzMucPlum1", "Mucor, control",  # fragmented gzMucPlum1 on the y axis
  "gzMucPlum3", "gzMucPlum5", "Mucor, near-identical",
  "gzMucPlum2", "gzMucPlum3", "Mucor, circled",
  "gzMucPlum2", "gzMucPlum5", "Mucor, circled",
  "gzMucFrag2", "gzMucFrag3", "Mucor",
  "MorAlpi3",   "gzMorAlpi1", "Mortierellales",
  "ggMorZona1", "gzLinHyal1", "Mortierellales"
) %>%
  mutate(pr = map2(a, b, pair_row),
         d = map_dbl(pr, "d_sub"),
         strand = map_dbl(pr, "strand_index_norm"),
         tau = map_dbl(pr, "tau_index_norm"),
         tau_raw = map_dbl(pr, "tau_index"),
         label = sprintf("%s - %s (d = %.3f)", a, b, d)) %>%
  select(-pr)

# --- Runs of inverted genes ---------------------------------------------------
inverted_runs <- function(an) {
  best <- an %>% count(seq_a, seq_b) %>%
    group_by(seq_a) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
    filter(n >= cfg$min_pair_n) %>%
    select(seq_a, seq_b, n_chr = n)
  an %>%
    inner_join(best, by = c("seq_a", "seq_b")) %>%
    arrange(seq_a, rank_a) %>%
    group_by(seq_a) %>%
    mutate(inverted = rel != sign(sum(rel)), run = consecutive_id(inverted)) %>%
    group_by(seq_a, seq_b, n_chr, run, inverted) %>%
    summarise(n_genes = n(), from_a = min(pos_a), to_a = max(pos_a),
              span_b = max(rank_b) - min(rank_b) + 1L, .groups = "drop")
}

runs <- pairs %>%
  mutate(runs = map2(a, b, ~ inverted_runs(pair_anchors(.x, .y)))) %>%
  select(a, b, label, runs) %>%
  unnest(runs)
write_tsv(runs %>% filter(inverted) %>% select(-inverted, -label),
          file.path(cfg$tables, "explore_inverted_runs.tsv"))

# Raw strand index (direction a -> b) from the inverted genes, and the tau
# expected if every inverted run were a clean inversion and nothing else moved
per_chr <- runs %>%
  group_by(label, seq_a, n_chr) %>%
  summarise(inv = sum(n_genes[inverted]),
            reversed_pairs = sum(choose(n_genes[inverted], 2)), .groups = "drop") %>%
  mutate(tau_inv = 1 - 2 * reversed_pairs / choose(n_chr, 2))
run_summary <- runs %>%
  filter(inverted) %>%
  group_by(label) %>%
  summarise(inverted_runs = n(),
            single_gene_runs = sum(n_genes == 1),
            largest_run = max(n_genes),
            runs_contiguous_in_b = mean(span_b == n_genes)) %>%
  right_join(per_chr %>% group_by(label) %>%
               summarise(genes = sum(n_chr), inverted_genes = sum(inv),
                         strand_raw = 1 - 2 * sum(inv) / sum(n_chr),
                         tau_if_inversions = weighted.mean(tau_inv, n_chr)),
             by = "label") %>%
  left_join(select(pairs, label, strand, tau, tau_raw), by = "label") %>%
  mutate(across(where(is.double), ~ round(.x, 3))) %>%
  arrange(match(label, pairs$label))
options(width = 200)
print(run_summary)

# --- Figure 1: whole-genome dot plots -----------------------------------------
dots <- pmap(pairs, function(a, b, title, ...) dotplot(a, b, title = title))
save_fig(wrap_plots(dots, ncol = 4, guides = "collect") & theme(legend.position = "bottom"),
         "explore_strand_tau_dotplots", width = 180, height = 110)

# --- Figure 2: inverted genes by run size, and one chromosome pair per pair --
size_levels <- c("1", "2-4", "5-19", "20-99", "100+")
size_colours <- setNames(c("#BFE3D3", "#7CC4A6", "#3FA07D", "#1B7A5C", "#0B4A37"), size_levels)
by_size <- runs %>%
  filter(inverted) %>%
  mutate(size = cut(n_genes, c(0, 1, 4, 19, 99, Inf), labels = size_levels)) %>%
  group_by(label, size) %>%
  summarise(genes = sum(n_genes), .groups = "drop") %>%
  left_join(per_chr %>% group_by(label) %>% summarise(total = sum(n_chr)), by = "label") %>%
  mutate(frac = genes / total, label = factor(label, levels = rev(pairs$label)))
bar_text <- run_summary %>%
  transmute(label = factor(label, levels = levels(by_size$label)),
            frac = inverted_genes / genes,
            text = sprintf("%d genes in %d runs (largest %d)", inverted_genes, inverted_runs,
                           largest_run))
bars <- ggplot(by_size, aes(frac, label, fill = size)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(data = bar_text, aes(frac, label, label = text), inherit.aes = FALSE,
            hjust = -0.05, size = 2.1, colour = "grey30") +
  scale_fill_manual(values = size_colours, name = "Genes per inverted run", drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.5))) +
  labs(x = "Fraction of genes in the minority orientation (raw strand index = 1 - 2 x fraction)",
       y = NULL) +
  theme_syn() + theme(legend.position = "bottom")

# Rank-rank dot plot of the chromosome pair with most inverted genes. As in the
# whole-genome dot plots, the partner is flipped when gene order runs
# backwards overall (negative rank correlation), and orientation is relative
# to that overall direction.
zoom_panel <- function(a, b) {
  an <- pair_anchors(a, b)
  chr <- inverted_runs(an) %>% filter(inverted) %>%
    count(seq_a, seq_b, wt = n_genes) %>% slice_max(n, n = 1, with_ties = FALSE)
  d <- an %>%
    semi_join(chr, by = c("seq_a", "seq_b")) %>%
    mutate(dir = sign(cor(pos_a, pos_b, method = "spearman")),
           x = rank(pos_a), y = rank(dir * pos_b),
           orientation = factor(if_else(rel == dir, "same", "inverted"),
                                levels = c("same", "inverted"))) %>%
    arrange(orientation)
  ggplot(d, aes(x, y, colour = orientation)) +
    geom_point(size = 0.35, stroke = 0) +
    scale_colour_manual(values = orientation_colours, guide = "none") +
    coord_equal() +
    labs(x = sprintf("%s, %s (gene rank)", a, chr$seq_a),
         y = sprintf("%s, %s (gene rank)", b, chr$seq_b),
         title = sprintf("%s - %s\n%d of %d genes inverted", a, b,
                         sum(d$orientation == "inverted"), nrow(d))) +
    theme_syn(7) +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.3, colour = "grey30"),
          axis.line = element_blank())
}
zooms <- wrap_plots(zoom_panel("gzMucPlum3", "gzMucPlum5"), zoom_panel("gzMucPlum2", "gzMucPlum3"),
                    zoom_panel("gzMucFrag2", "gzMucFrag3"), nrow = 1)
save_fig(bars / zooms + plot_layout(heights = c(1, 1.05)),
         "explore_strand_tau_runs", width = 180, height = 140)
