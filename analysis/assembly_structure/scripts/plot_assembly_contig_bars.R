#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) {
    return(default)
  }
  args[[hit + 1]]
}

contig_file <- get_arg("--contigs", "/cluster/work/users/hellb/Mucor/results/assembly_structure/assembly_contig_gene_summary.tsv")
out_pdf <- get_arg("--out-pdf", "/cluster/work/users/hellb/Mucor/results/assembly_structure/assembly_contig_bars.pdf")
out_png <- get_arg("--out-png", "/cluster/work/users/hellb/Mucor/results/assembly_structure/assembly_contig_bars.png")
tree_file <- get_arg("--tree", "")
ncol_facets <- as.integer(get_arg("--ncol", "6"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.")
}

contigs <- read.delim(contig_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("genome_id", "contig_rank", "length_mb", "gene_bearing_contig") %in% names(contigs))) {
  stop("Input table must include genome_id, contig_rank, length_mb, and gene_bearing_contig.")
}

genome_order <- unique(contigs$genome_id)
if (nzchar(tree_file) && requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::read.tree(tree_file)
  tree_order <- tree$tip.label[tree$tip.label %in% genome_order]
  genome_order <- c(tree_order, setdiff(genome_order, tree_order))
}

contigs$genome_id <- factor(contigs$genome_id, levels = genome_order)
contigs$gene_bearing_contig <- factor(contigs$gene_bearing_contig, levels = c("yes", "no"))

palette <- c(yes = "#5F5752", no = "#D0D5DD")

p <- ggplot2::ggplot(
  contigs,
  ggplot2::aes(x = contig_rank, y = length_mb, fill = gene_bearing_contig)
) +
  ggplot2::geom_col(width = 0.82, color = NA) +
  ggplot2::facet_wrap(~ genome_id, scales = "free_x", ncol = ncol_facets) +
  ggplot2::scale_fill_manual(
    values = palette,
    breaks = c("yes", "no"),
    labels = c("Gene-bearing contig", "No gene annotation")
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.04))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::labs(
    x = "Contigs ranked by length",
    y = "Contig length (Mb)",
    fill = NULL
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F2F4F7", color = "#D0D5DD"),
    strip.text = ggplot2::element_text(size = 6.5),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(out_pdf, p, width = 14, height = 10, units = "in", limitsize = FALSE)
ggplot2::ggsave(out_png, p, width = 14, height = 10, units = "in", dpi = 300, limitsize = FALSE)

message("Wrote: ", out_pdf)
message("Wrote: ", out_png)
