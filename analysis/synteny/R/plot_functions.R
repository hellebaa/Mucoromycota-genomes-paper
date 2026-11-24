# Plotting helpers shared by 04_figures.R and the exploration scripts: theme,
# figure export and Oxford dot plots.
#
# Functions use `cfg`, `genomes` (read_genome_table()) and `pi_tab`
# (results/tables/pair_indices.tsv) from the calling script.

suppressPackageStartupMessages(library(patchwork))

theme_syn <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
          axis.text = element_text(colour = "grey20"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0),
          legend.key.size = unit(3, "mm"),
          plot.title = element_text(size = base_size, face = "plain"),
          plot.tag = element_text(face = "bold", size = base_size + 2))
}

save_fig <- function(p, name, width, height) {
  ggsave(file.path(cfg$figures, paste0(name, ".pdf")), p, width = width, height = height,
         units = "mm", device = pdf, useDingbats = FALSE)
  ggsave(file.path(cfg$figures, paste0(name, ".png")), p, width = width, height = height,
         units = "mm", dpi = 300, bg = "white")
}

short_name <- function(g) {
  sp <- genomes$species[match(g, genomes$genome)]
  paste0(str_sub(sp, 1, 1), ". ", word(sp, 2))
}
italic_label <- function(g) sprintf("italic('%s')~'(%s)'", short_name(g), g)

# ---------------------------------------------------------------------------
# Oxford dot plots
# ---------------------------------------------------------------------------

# Keep chromosomes (contigs) holding >= 1.5% of anchors, order query contigs by
# their main target partner, and flip query contigs so that their main
# diagonal runs upwards. Relative orientation is recomputed after flipping.
prep_dotplot <- function(an, min_frac = 0.015) {
  min_n <- max(30, min_frac * nrow(an))
  a_keep <- count(an, seq_a) %>% filter(n >= min_n) %>% arrange(desc(n))
  b_keep <- count(an, seq_b) %>% filter(n >= min_n)
  an <- filter(an, seq_a %in% a_keep$seq_a, seq_b %in% b_keep$seq_b)
  b_info <- an %>%
    group_by(seq_b, seq_a) %>%
    summarise(n = n(), mid_a = median(pos_a),
              r = if (n() > 2) cor(pos_a, pos_b, method = "spearman") else 1,
              .groups = "drop") %>%
    group_by(seq_b) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(ai = match(seq_a, a_keep$seq_a)) %>%
    arrange(ai, mid_a) %>%
    mutate(flip = r < 0)
  len_a <- an %>% group_by(seq_a) %>% summarise(len = max(pos_a)) %>%
    arrange(match(seq_a, a_keep$seq_a)) %>%
    mutate(off = lag(cumsum(len), default = 0), idx = row_number())
  len_b <- an %>% group_by(seq_b) %>% summarise(len = max(pos_b)) %>%
    arrange(match(seq_b, b_info$seq_b)) %>%
    mutate(off = lag(cumsum(len), default = 0), idx = row_number())
  pts <- an %>%
    left_join(select(b_info, seq_b, flip), by = "seq_b") %>%
    left_join(select(len_a, seq_a, off_a = off), by = "seq_a") %>%
    left_join(select(len_b, seq_b, off_b = off, len_b = len), by = "seq_b") %>%
    mutate(pos_b = if_else(flip, len_b - pos_b, pos_b),
           rel = if_else(flip, -rel, rel),
           x = (off_a + pos_a) / 1e6, y = (off_b + pos_b) / 1e6,
           orientation = factor(if_else(rel == 1L, "same", "inverted"),
                                levels = c("same", "inverted")))
  list(points = pts, len_a = len_a, len_b = len_b)
}

genome_cache <- new.env()
get_genome <- function(g) {
  if (is.null(genome_cache[[g]])) genome_cache[[g]] <- load_genome(g, cfg)
  genome_cache[[g]]
}

pair_anchors <- function(a, b) {
  ia <- match(a, genomes$genome); ib <- match(b, genomes$genome)
  if (ia > ib) return(swap_anchors(pair_anchors(b, a)))
  build_anchors(get_genome(a), get_genome(b), genomes$of_name[ia], genomes$of_name[ib], cfg)
}

pair_row <- function(a, b) {
  pi_tab %>% filter((.data$a == !!a & .data$b == !!b) | (.data$a == !!b & .data$b == !!a))
}

orientation_colours <- c(same = "#262626", inverted = "#1B9E77")

dotplot <- function(a, b, title = NULL) {
  dp <- prep_dotplot(pair_anchors(a, b))
  pr <- pair_row(a, b)
  clade <- as.character(genomes$order[match(a, genomes$genome)])
  col <- clade_colours[[clade]]
  sub <- sprintf("d = %.3f | strand %.2f | tau %.2f", pr$d_sub, pr$strand_index_norm,
                 pr$tau_index_norm)
  brk <- function(len) (len$off + len$len / 2) / 1e6
  # label only sequences wide enough for a label (>= 3% of the axis)
  lab <- function(len) if_else(len$len / sum(len$len) >= 0.03, as.character(len$idx), "")
  ggplot(dp$points, aes(x, y, colour = orientation)) +
    geom_vline(xintercept = dp$len_a$off[-1] / 1e6, linewidth = 0.2, colour = "grey85") +
    geom_hline(yintercept = dp$len_b$off[-1] / 1e6, linewidth = 0.2, colour = "grey85") +
    geom_point(size = 0.25, alpha = 0.8, stroke = 0) +
    scale_colour_manual(values = orientation_colours,
                        name = "Gene orientation relative to the other genome", drop = FALSE) +
    scale_x_continuous(breaks = brk(dp$len_a), labels = lab(dp$len_a), expand = c(0, 0)) +
    scale_y_continuous(breaks = brk(dp$len_b), labels = lab(dp$len_b), expand = c(0, 0)) +
    labs(x = parse(text = italic_label(a)), y = parse(text = italic_label(b)),
         title = title, subtitle = sub) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_syn() +
    theme(aspect.ratio = 1,
          axis.text = element_text(size = 5),
          axis.ticks = element_blank(),
          panel.border = element_rect(fill = NA, linewidth = 0.3, colour = "grey30"),
          axis.line = element_blank(),
          plot.subtitle = element_text(size = 6, colour = "grey30"),
          plot.title = element_text(face = "bold", colour = col))
}
