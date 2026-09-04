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
out_prefix <- get_arg("--out-prefix", "/Users/hellebaalsrud/.codex/worktrees/a8cb/Mucoromycota-genomes-paper/outputs/01a03cff-390b-7840-a99b-2d8185721340/assembly_structure_vertical_bars_draft")
ncol_facets <- as.integer(get_arg("--ncol", "4"))
order_filter <- get_arg("--order", "")
genome_prefix_regex <- get_arg("--genome-prefix-regex", "")
show_candidate_markers <- get_arg("--show-markers", "yes")

read_optional_tsv <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

contigs <- read_tsv(file.path(results_dir, "assembly_contig_gene_summary.tsv"), show_col_types = FALSE, progress = FALSE)
telomere_multi_file <- file.path(results_dir, "telomere_contig_summary_multi_motif.tsv")
telomere_first_pass_file <- file.path(results_dir, "telomere_contig_summary.tsv")
telomeres <- if (file.exists(telomere_multi_file)) {
  read_optional_tsv(telomere_multi_file)
} else {
  read_optional_tsv(telomere_first_pass_file)
}
cen <- read_optional_tsv(file.path(results_dir, "centromere_contig_summary.tsv"))
grem <- read_optional_tsv(file.path(results_dir, "gremline_contig_summary.tsv"))
metadata <- read_csv(metadata_file, show_col_types = FALSE, progress = FALSE)

metadata_order <- metadata %>%
  filter(include_core_analyses == "yes") %>%
  mutate(order_group = factor(order, levels = c("Mucorales", "Umbelopsidales", "Mortierellales"))) %>%
  arrange(order_group, family, genome_id) %>%
  select(genome_id, order, family)

if (nzchar(order_filter)) {
  metadata_order <- metadata_order %>% filter(order == order_filter)
}
if (nzchar(genome_prefix_regex)) {
  metadata_order <- metadata_order %>% filter(grepl(genome_prefix_regex, genome_id))
}

contigs_plot <- contigs %>%
  semi_join(metadata_order, by = "genome_id") %>%
  left_join(metadata_order, by = "genome_id", suffix = c("", "_metadata")) %>%
  mutate(
    genome_id = factor(genome_id, levels = metadata_order$genome_id),
    contig_rank = as.integer(contig_rank),
    gene_bearing_contig = factor(gene_bearing_contig, levels = c("yes", "no"))
  )

label_data <- contigs_plot %>%
  group_by(genome_id) %>%
  summarize(
    n_contigs = n(),
    assembly_mb = sum(length_mb),
    .groups = "drop"
  ) %>%
  mutate(facet_label = paste0(genome_id, "\n", round(assembly_mb, 1), " Mb; ", n_contigs, " contigs"))

label_lookup <- setNames(as.character(label_data$facet_label), as.character(label_data$genome_id))
contigs_plot <- contigs_plot %>%
  mutate(facet_id = factor(label_lookup[as.character(genome_id)], levels = label_data$facet_label))

telomere_marks <- telomeres %>%
  semi_join(contigs_plot %>% select(genome_id, contig_id), by = c("genome_id", "contig_id")) %>%
  mutate(
    n_telomere_ends = (left_telomere == "yes") + (right_telomere == "yes"),
    telomere_label = case_when(
      n_telomere_ends == 2 ~ "Both telomeres",
      n_telomere_ends == 1 ~ "One telomere",
      TRUE ~ "No telomere"
    )
  ) %>%
  filter(n_telomere_ends > 0) %>%
  left_join(contigs_plot %>% select(genome_id, contig_id, contig_rank, length_mb, facet_id), by = c("genome_id", "contig_id")) %>%
  filter(!is.na(facet_id))

candidate_marks <- bind_rows(
  if (nrow(cen) > 0) cen %>% mutate(marker_type = "CEN-core BLAST") else tibble(),
  if (nrow(grem) > 0) grem %>% mutate(marker_type = "GremLINE-like") else tibble()
) %>%
  semi_join(contigs_plot %>% select(genome_id, contig_id), by = c("genome_id", "contig_id")) %>%
  left_join(contigs_plot %>% select(genome_id, contig_id, contig_rank, length_mb, facet_id), by = c("genome_id", "contig_id")) %>%
  mutate(
    marker_type = factor(marker_type, levels = c("CEN-core BLAST", "GremLINE-like")),
    marker_x = contig_rank + case_when(
      marker_type == "CEN-core BLAST" ~ -0.18,
      marker_type == "GremLINE-like" ~ 0.18,
      TRUE ~ 0
    ),
    marker_y = ((candidate_region_start + candidate_region_end) / 2) / 1e6,
    marker_y = pmin(marker_y, length_mb)
  ) %>%
  filter(!is.na(facet_id), !is.na(marker_y))

plot_obj <- ggplot(contigs_plot, aes(x = contig_rank, y = length_mb)) +
  geom_col(aes(fill = gene_bearing_contig), width = 0.82, color = NA)

if (show_candidate_markers == "yes" && nrow(candidate_marks) > 0) {
  plot_obj <- plot_obj +
    geom_point(
      data = candidate_marks,
      aes(x = marker_x, y = marker_y, color = marker_type),
      inherit.aes = FALSE,
      size = 0.75,
      alpha = 0.85
    )
}

if (nrow(telomere_marks) > 0) {
  plot_obj <- plot_obj +
    geom_point(
      data = telomere_marks,
      aes(x = contig_rank, y = length_mb, shape = telomere_label),
      inherit.aes = FALSE,
      size = 1.0,
      color = "#16817A",
      fill = "#16817A"
    )
}

plot_obj <- plot_obj +
  facet_wrap(~ facet_id, scales = "free_x", ncol = ncol_facets) +
  scale_fill_manual(
    values = c(yes = "#B7B2AD", no = "#E7EAEE"),
    labels = c(yes = "Gene-bearing contig", no = "No gene annotation")
  ) +
  scale_color_manual(values = c("CEN-core BLAST" = "#7A4EA3", "GremLINE-like" = "#B6782C")) +
  scale_shape_manual(values = c("One telomere" = 24, "Both telomeres" = 25)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.04))) +
  labs(
    x = "Contigs ranked by length",
    y = "Contig length (Mb)",
    fill = NULL,
    color = NULL,
    shape = NULL,
    caption = "Points show candidate evidence layers: purple, CEN-core BLAST intervals; ochre, GremLINE-like intervals; teal triangles, terminal TTAGGG/CCCTAA telomere motif arrays."
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F2F4F7", color = "#D0D5DD"),
    strip.text = element_text(size = 5.8, lineheight = 0.95),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 5.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.caption = element_text(size = 7, hjust = 0),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(paste0(out_prefix, ".pdf"), plot_obj, width = 14, height = 18, units = "in", limitsize = FALSE)
ggsave(paste0(out_prefix, ".png"), plot_obj, width = 14, height = 18, units = "in", dpi = 300, limitsize = FALSE)

message("Wrote: ", paste0(out_prefix, ".pdf"))
message("Wrote: ", paste0(out_prefix, ".png"))
