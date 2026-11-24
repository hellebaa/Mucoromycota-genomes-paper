# Publication figures for the synteny / mesosynteny analysis.
#
# Output (PDF + PNG) in results/figures/
#   fig_mesosynteny_main      dot plots at matched divergence + decay + rates
#   figS_dotplots_shallow     more dot plots of closely related pairs
#   figS_indices_vs_divergence all indices vs divergence, full range
#   figS_divergence_measures  primary index vs the three divergence measures
#   figS_rates_forest         rate estimates, all indices and sensitivity settings
#   figS_decay_model_fits     stretched-exponential fits used for the rates
#   figS_mesosynteny_map      same-chromosome index vs strand index per pair
#   figS_heatmap_strand       strand index for all pairs, genomes in tree order
#   figS_decay_mucor_split    Figure 1B with Mucorales split into Mucor and the
#                             other genera (D16)

source("analysis/synteny/R/functions.R")
source("analysis/synteny/R/plot_functions.R")
cfg <- syn_config()
set.seed(cfg$seed)  # jittered points
dir.create(cfg$figures, recursive = TRUE, showWarnings = FALSE)

genomes <- read_genome_table(cfg)
pi_tab <- read_tsv(file.path(cfg$tables, "pair_indices.tsv"), show_col_types = FALSE) %>%
  mutate(same_chr = 1 - inter_break, d_time = 2 * age_rel)
fits <- read_tsv(file.path(cfg$tables, "decay_fits.tsv"), show_col_types = FALSE)
rates <- read_tsv(file.path(cfg$tables, "rates_by_clade.tsv"), show_col_types = FALSE)
tests <- read_tsv(file.path(cfg$tables, "rate_tests.tsv"), show_col_types = FALSE)

index_labels <- c(
  strand_index_norm = "Strand index",
  tau_index_norm    = "Kendall's tau index",
  corr_index_norm   = "Correlation index",
  goc_intra_norm    = "GOC (vicinity 4)",
  goc1_intra_norm   = "Conserved adjacency",
  same_chr          = "Neighbour on same chromosome",
  synteny_index     = "Synteny index",
  macro_frac        = "Genes on syntenic chromosome pairs",
  block_mean        = "Genes per collinear block"
)
div_labels <- c(d_sub = "Divergence (CDS subs/site)",
                d_prot = "Divergence (protein, OrthoFinder tree)",
                d_time = "Divergence (time; 1 = 2 x root age)")

scale_clade <- list(
  scale_colour_manual(values = clade_colours, name = NULL, drop = FALSE),
  scale_fill_manual(values = clade_colours, name = NULL, drop = FALSE),
  scale_shape_manual(values = c(Mortierellales = 16, Umbelopsidales = 17, Mucorales = 15),
                     name = NULL, drop = FALSE)
)

# Mucorales split into Mucor and the other genera (D16)
mucor_group_labels <- expression("Mortierellales", "Umbelopsidales", italic("Mucor"),
                                 "Other Mucorales")
scale_mucor <- list(
  scale_colour_manual(values = mucor_group_colours, breaks = mucor_group_levels,
                      labels = mucor_group_labels, name = NULL, drop = FALSE),
  scale_shape_manual(values = c(Mortierellales = 16, Umbelopsidales = 17, Mucor = 15,
                                `Other Mucorales` = 18),
                     breaks = mucor_group_levels, labels = mucor_group_labels,
                     name = NULL, drop = FALSE)
)

# Matched-divergence examples, EBP-Nor genomes where possible
main_pairs <- tribble(
  ~clade,           ~level,    ~a,           ~b,
  "Mortierellales", "shallow", "ggMorZona1", "gzLinHyal1",
  "Umbelopsidales", "shallow", "gzUmbIsab1", "UmbIsab2",
  "Mucorales",      "shallow", "gzMucFrag1", "gzMucFrag2",
  "Mortierellales", "deep",    "gzLinHyal1", "gzMorAlpi1",
  "Umbelopsidales", "deep",    "gzUmbRama1", "gzUmbVina1",
  "Mucorales",      "deep",    "gzMucPlum1", "ggMucLanc1"
)
main_dots <- pmap(main_pairs, function(clade, level, a, b) dotplot(a, b, title = clade))

# ---------------------------------------------------------------------------
# Decay of indices with divergence
# ---------------------------------------------------------------------------

within <- pi_tab %>% filter(within) %>% mutate(clade = factor(clade, levels = clade_levels))

# curve = "loess": descriptive smoother per clade; "model": stretched-exponential
# fit used for the rate estimates (03_statistics.R)
decay_panel <- function(index, divergence = "d_sub", xmax = cfg$d_max, curve = "loess",
                        df = within, scale = scale_clade) {
  dd <- df %>% filter(.data[[divergence]] <= xmax) %>%
    transmute(clade, d = .data[[divergence]], value = .data[[index]])
  p <- ggplot(dd, aes(d, value, colour = clade, shape = clade)) +
    geom_hline(yintercept = 0, linewidth = 0.2, linetype = 2, colour = "grey60") +
    geom_point(size = 1, alpha = 0.7, stroke = 0)
  ft <- fits %>% filter(.data$index == !!index, .data$divergence == !!divergence)
  if (curve == "model" && nrow(ft) > 0) {
    curves <- ft %>%
      mutate(clade = factor(clade, levels = clade_levels)) %>%
      crossing(d = seq(0, xmax, length.out = 200)) %>%
      mutate(value = c * exp(-(lambda * d)^beta))
    p <- p + geom_line(data = curves, linewidth = 0.6, show.legend = FALSE)
  } else if (curve == "loess") {
    p <- p + geom_smooth(aes(group = clade), method = "loess", formula = y ~ x, span = 1,
                         method.args = list(degree = 1), se = FALSE, linewidth = 0.6,
                         show.legend = FALSE)
  }
  p + scale +
    scale_y_continuous(limits = c(min(0, min(dd$value, na.rm = TRUE)), 1.02),
                       breaks = seq(0, 1, 0.25)) +
    labs(x = div_labels[[divergence]], y = index_labels[[index]]) +
    theme_syn()
}

# Model-free comparison at matched divergence: index values of all pairs in
# `df` with d_min_rate < d <= 0.15 (points), the median over independent-pair
# draws (bar; `med`) and one-sided p from 03_statistics.R (`pv`, printed at
# x position `x`).
matched_vals <- read_tsv(file.path(cfg$tables, "matched_values_by_clade.tsv"), show_col_types = FALSE)
matched_tests <- read_tsv(file.path(cfg$tables, "matched_tests.tsv"), show_col_types = FALSE)
d_matched <- 0.15
matched_idx <- c(strand_index_norm = "Strand", tau_index_norm = "Tau", same_chr = "Same chr.")
matched_panel <- function(df = within, med = matched_vals,
                          pv = matched_tests %>%
                            filter(comparison == "Mucorales > Mortierellales") %>% mutate(x = 2),
                          scale = scale_clade,
                          xlabs = c(Mortierellales = "Mor", Umbelopsidales = "Umb",
                                    Mucorales = "Muc")) {
  as_index <- function(x) factor(matched_idx[x], levels = matched_idx)
  pts <- df %>%
    filter(d_sub > cfg$d_min_rate, d_sub <= d_matched) %>%
    select(clade, all_of(names(matched_idx))) %>%
    pivot_longer(all_of(names(matched_idx)), names_to = "index") %>%
    mutate(index = as_index(index))
  med <- med %>% filter(index %in% names(matched_idx)) %>%
    mutate(index = as_index(index), clade = factor(clade, levels = levels(df$clade)))
  pv <- pv %>% filter(index %in% names(matched_idx)) %>%
    mutate(index = as_index(index),
           label = if_else(median_p < 0.001, "p < 0.001", sprintf("p = %.3f", median_p)))
  ggplot(pts, aes(clade, value, colour = clade, shape = clade)) +
    geom_jitter(width = 0.2, height = 0, size = 1, alpha = 0.7, stroke = 0) +
    geom_errorbar(data = med, aes(y = median_value, ymin = median_value, ymax = median_value),
                  width = 0.6, linewidth = 0.7, show.legend = FALSE) +
    geom_text(data = pv, aes(x = x, y = 0.02, label = label), inherit.aes = FALSE,
              size = 2.2, colour = "grey30") +
    facet_wrap(~index, nrow = 1) +
    scale +
    scale_x_discrete(labels = xlabs) +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = sprintf("Index value, %.2f < d <= %.2f", cfg$d_min_rate, d_matched)) +
    theme_syn() +
    theme(panel.spacing = unit(1.5, "mm"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
}

# Mesosynteny map: chromosome-scale order vs chromosome membership
meso_map <- within %>% filter(d_sub <= cfg$d_max, d_sub > cfg$d_min_rate) %>%
  ggplot(aes(same_chr, strand_index_norm, colour = clade, shape = clade)) +
  geom_point(size = 1.1, alpha = 0.75, stroke = 0) +
  annotate("text", x = 0.985, y = 0.03, label = "mesosynteny", hjust = 1, vjust = 0,
           size = 2.4, fontface = "italic", colour = "grey35") +
  annotate("text", x = 0.985, y = 1.0, label = "macrosynteny", hjust = 1, vjust = 1,
           size = 2.4, fontface = "italic", colour = "grey35") +
  scale_clade +
  labs(x = index_labels[["same_chr"]], y = index_labels[["strand_index_norm"]]) +
  theme_syn()

# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

dots <- wrap_plots(main_dots, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
bottom <- (decay_panel("strand_index_norm") + decay_panel("tau_index_norm") +
             decay_panel("same_chr") + matched_panel()) +
  plot_layout(ncol = 4, guides = "collect", widths = c(1, 1, 1, 1.5)) &
  theme(legend.position = "bottom")
main <- wrap_elements(dots) / bottom + plot_layout(heights = c(2.1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B")))
save_fig(main, "fig_mesosynteny_main", width = 180, height = 230)

# ---------------------------------------------------------------------------
# Supplementary figures
# ---------------------------------------------------------------------------

shallow_pairs <- tribble(
  ~a,           ~b,
  "MorAlpi3",   "gzMorAlpi1",
  "LinElon",    "gzLinHyal1",
  "MorAlpi2",   "gzMorAlpi1",
  "gzMucCirc1", "gzMucCirc2",
  "gzMucFrag2", "gzMucFrag3",
  "gzMucPlum2", "gzMucPlum7",
  "AbsRepe",    "gzAbsCyli1",
  "gzMucFlav2", "gzMucFlav3",
  "CunBert",    "gzCunBlak1",
  "LicRamo",    "gzLicCory1",
  "UmbVina3",   "gzUmbRama1",
  "gzUmbVina1", "gzUmbVina2"
) %>%
  left_join(select(pi_tab, a, b, d_sub), by = c("a", "b")) %>%
  mutate(d_sub = coalesce(d_sub, map2_dbl(a, b, ~ pair_row(.x, .y)$d_sub))) %>%
  mutate(clade = as.character(genomes$order[match(a, genomes$genome)])) %>%
  arrange(factor(clade, levels = clade_levels), d_sub)
s_dots <- pmap(shallow_pairs, function(a, b, d_sub, clade) dotplot(a, b, title = clade))
save_fig(wrap_plots(s_dots, ncol = 4, guides = "collect") & theme(legend.position = "bottom"),
         "figS_dotplots_shallow", width = 180, height = 170)

all_idx <- c("strand_index_norm", "tau_index_norm", "corr_index_norm", "goc_intra_norm",
             "goc1_intra_norm", "same_chr", "synteny_index", "macro_frac", "block_mean")
s_idx <- within %>%
  select(clade, d_sub, all_of(all_idx)) %>%
  pivot_longer(all_of(all_idx), names_to = "index") %>%
  mutate(index = factor(index_labels[index], levels = index_labels[all_idx])) %>%
  ggplot(aes(d_sub, value, colour = clade, shape = clade)) +
  geom_vline(xintercept = cfg$d_max, linewidth = 0.2, linetype = 2, colour = "grey60") +
  geom_point(size = 0.8, alpha = 0.6, stroke = 0) +
  facet_wrap(~index, scales = "free_y", ncol = 3) +
  scale_clade +
  labs(x = div_labels[["d_sub"]], y = "Index value") +
  theme_syn() + theme(legend.position = "bottom")
save_fig(s_idx, "figS_indices_vs_divergence", width = 180, height = 170)

s_div <- map(names(div_labels), ~ decay_panel("strand_index_norm", .x,
                                             xmax = max(within[[.x]][within$d_sub <= cfg$d_max]))) %>%
  wrap_plots(nrow = 1, guides = "collect") & theme(legend.position = "bottom")
save_fig(s_div, "figS_divergence_measures", width = 180, height = 70)

settings <- tribble(
  ~genome_set, ~divergence, ~model,        ~setting,
  "all",       "d_sub",     "stretched",   "All genomes, CDS distance\n(primary)",
  "all",       "d_prot",    "stretched",   "All genomes, protein distance",
  "all",       "d_time",    "stretched",   "All genomes, time",
  "ebpnor",    "d_sub",     "stretched",   "EBP-Nor genomes, CDS distance",
  "all",       "d_sub",     "exponential", "All genomes, CDS distance,\nexponential model"
)
rate_floor <- 1e-3  # rates of 0 (index at the noise ceiling) drawn at this value
forest <- rates %>%
  inner_join(settings, by = c("genome_set", "divergence", "model")) %>%
  mutate(setting = factor(setting, levels = settings$setting),
         index = factor(index_labels[index], levels = rev(index_labels)),
         clade = factor(clade, levels = clade_levels),
         across(c(median_rate, rate_lo, rate_hi), ~ pmax(.x, rate_floor))) %>%
  ggplot(aes(median_rate, index, colour = clade, shape = clade)) +
  geom_linerange(aes(xmin = rate_lo, xmax = rate_hi), position = position_dodge(0.6),
                 linewidth = 0.4) +
  geom_point(position = position_dodge(0.6), size = 1.2) +
  facet_wrap(~setting, ncol = 3, scales = "free_x") +
  scale_x_log10() +
  scale_clade +
  labs(x = "Decay rate lambda (median over independent-pair draws, 95% range)", y = NULL) +
  theme_syn() + theme(legend.position = "bottom")
save_fig(forest, "figS_rates_forest", width = 180, height = 150)

fit_idx <- c("strand_index_norm", "tau_index_norm", "corr_index_norm", "goc_intra_norm", "same_chr")
s_fits <- map(fit_idx, ~ decay_panel(.x, curve = "model")) %>%
  wrap_plots(ncol = 3, guides = "collect") & theme(legend.position = "bottom")
save_fig(s_fits, "figS_decay_model_fits", width = 180, height = 120)

save_fig(meso_map + theme(legend.position = "bottom"), "figS_mesosynteny_map",
         width = 90, height = 90)

# Heatmap of the strand index, genomes in phylogenetic order
tree <- read.tree(cfg$tree_dated)
tree <- keep.tip(tree, intersect(tree$tip.label, genomes$genome)) %>% ladderize(right = FALSE)
is_tip <- tree$edge[, 2] <= length(tree$tip.label)
tip_order <- tree$tip.label[tree$edge[is_tip, 2]]
hm <- bind_rows(pi_tab %>% select(a, b, value = strand_index_norm),
                pi_tab %>% select(a = b, b = a, value = strand_index_norm)) %>%
  mutate(a = factor(a, levels = tip_order), b = factor(b, levels = rev(tip_order)))
clade_bar <- genomes %>% mutate(genome = factor(genome, levels = tip_order))
s_hm <- ggplot(hm, aes(a, b, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "#F2F0EB", high = "#1F3A5F", limits = c(0, 1),
                      name = "Strand index", na.value = "white") +
  labs(x = NULL, y = NULL) +
  theme_syn(6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5),
        axis.text.y = element_text(size = 4.5), axis.line = element_blank(),
        axis.ticks = element_blank(), aspect.ratio = 1)
bar <- ggplot(clade_bar, aes(genome, 1, fill = order)) + geom_tile() + scale_clade +
  theme_void() + theme(legend.position = "none")
save_fig(bar / s_hm + plot_layout(heights = c(0.03, 1)), "figS_heatmap_strand",
         width = 180, height = 190)

# Figure 1B with Mucorales split into Mucor and the other genera (D16); p in
# the matched panel: each group vs Mortierellales
split_within <- split_mucor(within, read.tree(cfg$tree_ml))
split_vals <- read_tsv(file.path(cfg$tables, "mucor_split_by_group.tsv"), show_col_types = FALSE) %>%
  filter(genome_set == "mucor_matched") %>%
  rename(median_value = median)
split_pv <- read_tsv(file.path(cfg$tables, "mucor_split_tests.tsv"), show_col_types = FALSE) %>%
  filter(genome_set == "mucor_matched", str_ends(comparison, " > Mortierellales")) %>%
  mutate(x = match(sub(" > Mortierellales$", "", comparison), mucor_group_levels))
split_decay <- map(names(matched_idx), ~ decay_panel(.x, df = split_within, scale = scale_mucor)) %>%
  wrap_plots(nrow = 1)
split_matched <- matched_panel(split_within, split_vals, split_pv, scale_mucor,
                               xlabs = expression(Mor, Umb, italic(Mucor), "Other Muc."))
save_fig(split_decay / split_matched + plot_layout(guides = "collect") &
           theme(legend.position = "bottom"),
         "figS_decay_mucor_split", width = 180, height = 130)

message("Figures written to ", cfg$figures)
