#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) {
    return(default)
  }
  args[[hit + 1]]
}

results_dir <- get_arg("--results-dir", "/Users/hellebaalsrud/Documents/Projects/Mucormycota_NeLS_staging/results/assembly_structure")
metadata_file <- get_arg("--metadata", "/Users/hellebaalsrud/Documents/Projects/Mucormycota_NeLS_staging/data/metadata/genome_metadata_master.csv")
out_prefix <- get_arg("--out-prefix", "/Users/hellebaalsrud/.codex/worktrees/a8cb/Mucoromycota-genomes-paper/outputs/01a03cff-390b-7840-a99b-2d8185721340/assembly_structure_tracks_draft")
ncol_facets <- as.integer(get_arg("--ncol", "5"))
order_filter <- get_arg("--order", "")
genome_prefix_regex <- get_arg("--genome-prefix-regex", "")

read_optional_tsv <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

contigs <- read_tsv(file.path(results_dir, "assembly_contig_gene_summary.tsv"), show_col_types = FALSE, progress = FALSE)
telomeres <- read_optional_tsv(file.path(results_dir, "telomere_contig_summary.tsv"))
cen <- read_optional_tsv(file.path(results_dir, "centromere_contig_summary.tsv"))
grem <- read_optional_tsv(file.path(results_dir, "gremline_contig_summary.tsv"))
metadata <- read_csv(metadata_file, show_col_types = FALSE, progress = FALSE)

metadata_order <- metadata %>%
  filter(include_core_analyses == "yes") %>%
  mutate(order_group = factor(order, levels = c("Mucorales", "Umbelopsidales", "Mortierellales"))) %>%
  arrange(order_group, family, genome_id) %>%
  select(genome_id, order, family)

contigs_plot <- contigs %>%
  left_join(metadata_order, by = "genome_id", suffix = c("", "_metadata")) %>%
  mutate(
    order = coalesce(order_metadata, order),
    family = coalesce(family_metadata, family),
    genome_id = factor(genome_id, levels = metadata_order$genome_id),
    contig_rank = as.integer(contig_rank),
    y = -contig_rank,
    gene_bearing_contig = factor(gene_bearing_contig, levels = c("yes", "no"))
  )

if (nzchar(order_filter)) {
  keep_genomes <- metadata_order %>%
    filter(order == order_filter) %>%
    pull(genome_id)
  contigs_plot <- contigs_plot %>% filter(as.character(genome_id) %in% keep_genomes)
  metadata_order <- metadata_order %>% filter(genome_id %in% keep_genomes)
}

if (nzchar(genome_prefix_regex)) {
  keep_genomes <- metadata_order %>%
    filter(grepl(genome_prefix_regex, genome_id)) %>%
    pull(genome_id)
  contigs_plot <- contigs_plot %>% filter(as.character(genome_id) %in% keep_genomes)
  metadata_order <- metadata_order %>% filter(genome_id %in% keep_genomes)
}

telomere_marks <- telomeres %>%
  select(genome_id, contig_id, left_telomere, right_telomere) %>%
  pivot_longer(
    cols = c(left_telomere, right_telomere),
    names_to = "side",
    values_to = "telomere"
  ) %>%
  filter(telomere == "yes") %>%
  left_join(contigs_plot %>% select(genome_id, contig_id, length_mb, y), by = c("genome_id", "contig_id")) %>%
  mutate(
    x = if_else(side == "left_telomere", 0, length_mb),
    side = recode(side, left_telomere = "Left telomere", right_telomere = "Right telomere"),
    genome_id = factor(genome_id, levels = metadata_order$genome_id)
  )

cen_marks <- cen %>%
  mutate(
    cen_mid_mb = ((candidate_region_start + candidate_region_end) / 2) / 1e6,
    cen_span_mb = pmax((candidate_region_end - candidate_region_start) / 1e6, 0.03)
  ) %>%
  left_join(contigs_plot %>% select(genome_id, contig_id, y), by = c("genome_id", "contig_id")) %>%
  filter(!is.na(y)) %>%
  mutate(genome_id = factor(genome_id, levels = metadata_order$genome_id))

grem_marks <- grem %>%
  mutate(
    grem_mid_mb = ((candidate_region_start + candidate_region_end) / 2) / 1e6,
    grem_span_mb = pmax((candidate_region_end - candidate_region_start) / 1e6, 0.03)
  ) %>%
  left_join(contigs_plot %>% select(genome_id, contig_id, y), by = c("genome_id", "contig_id")) %>%
  filter(!is.na(y)) %>%
  mutate(genome_id = factor(genome_id, levels = metadata_order$genome_id))

label_data <- contigs_plot %>%
  group_by(genome_id) %>%
  summarize(
    order = first(order),
    family = first(family),
    n_contigs = n(),
    gene_bearing = sum(gene_bearing_contig == "yes"),
    assembly_mb = sum(length_mb),
    .groups = "drop"
  ) %>%
  mutate(facet_label = paste0(genome_id, "\n", round(assembly_mb, 1), " Mb; ", n_contigs, " contigs"))

label_lookup <- setNames(as.character(label_data$facet_label), as.character(label_data$genome_id))
contigs_plot$facet_id <- factor(label_lookup[as.character(contigs_plot$genome_id)], levels = label_data$facet_label)
telomere_marks$facet_id <- factor(label_lookup[as.character(telomere_marks$genome_id)], levels = label_data$facet_label)
cen_marks$facet_id <- factor(label_lookup[as.character(cen_marks$genome_id)], levels = label_data$facet_label)
grem_marks$facet_id <- factor(label_lookup[as.character(grem_marks$genome_id)], levels = label_data$facet_label)

telomere_marks <- telomere_marks %>% filter(!is.na(facet_id), !is.na(y), !is.na(x))
cen_marks <- cen_marks %>% filter(!is.na(facet_id), !is.na(y), !is.na(cen_mid_mb))
grem_marks <- grem_marks %>% filter(!is.na(facet_id), !is.na(y), !is.na(grem_mid_mb))

plot_obj <- ggplot(contigs_plot, aes(y = y)) +
  geom_segment(
    aes(x = 0, xend = length_mb, yend = y, color = gene_bearing_contig),
    linewidth = 2.2,
    lineend = "butt"
  ) +
  geom_segment(
    data = grem_marks,
    aes(x = pmax(grem_mid_mb - grem_span_mb / 2, 0), xend = grem_mid_mb + grem_span_mb / 2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 1.1,
    color = "#B6782C",
    alpha = 0.55
  ) +
  geom_segment(
    data = cen_marks,
    aes(x = pmax(cen_mid_mb - cen_span_mb / 2, 0), xend = cen_mid_mb + cen_span_mb / 2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 1.3,
    color = "#7A4EA3"
  ) +
  geom_point(
    data = telomere_marks,
    aes(x = x, y = y, shape = side),
    inherit.aes = FALSE,
    size = 1.2,
    color = "#16817A",
    fill = "#16817A"
  ) +
  facet_wrap(~ facet_id, scales = "free_y", ncol = ncol_facets) +
  scale_color_manual(
    values = c(yes = "#5F5752", no = "#D0D5DD"),
    labels = c(yes = "Gene-bearing contig", no = "No gene annotation")
  ) +
  scale_shape_manual(values = c("Left telomere" = 17, "Right telomere" = 25)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.04))) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(
    x = "Position along contig (Mb)",
    color = NULL,
    shape = NULL,
    caption = "Purple: CEN-core BLAST candidate interval. Ochre: GremLINE-like interval. Teal triangles: terminal TTAGGG/CCCTAA telomere motif arrays."
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F2F4F7", color = "#D0D5DD"),
    strip.text = element_text(size = 5.8, lineheight = 0.95),
    axis.text.x = element_text(size = 5.5),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.caption = element_text(size = 7, hjust = 0),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(paste0(out_prefix, ".pdf"), plot_obj, width = 14, height = 18, units = "in", limitsize = FALSE)
ggsave(paste0(out_prefix, ".png"), plot_obj, width = 14, height = 18, units = "in", dpi = 300, limitsize = FALSE)

message("Wrote: ", paste0(out_prefix, ".pdf"))
message("Wrote: ", paste0(out_prefix, ".png"))
