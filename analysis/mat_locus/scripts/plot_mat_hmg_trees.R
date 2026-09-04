#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: plot_mat_hmg_trees.R <tree.nwk> <N0.tsv> <output-dir> <tree-tag>\n",
    "Example: Rscript plot_mat_hmg_trees.R data/.../OG0000175_tree.txt data/.../N0.tsv outputs/mat_trees resolved"
  )
}

tree_file <- args[[1]]
hog_file <- args[[2]]
output_dir <- args[[3]]
tree_tag <- args[[4]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sexm_hog <- "N0.HOG0000705"
sexp_hog <- "N0.HOG0000708"

read_hog_members <- function(path, hog_id) {
  rows <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  row <- rows[rows$HOG == hog_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Expected exactly one row for ", hog_id, ", found ", nrow(row))
  }

  meta_cols <- c("HOG", "OG", "Gene Tree Parent Clade")
  species_cols <- setdiff(names(row), meta_cols)
  members <- character()
  for (column in species_cols) {
    value <- row[[column]]
    if (is.na(value) || !nzchar(value)) {
      next
    }
    genes <- strsplit(value, ", ", fixed = TRUE)[[1]]
    genes <- genes[nzchar(genes)]
    members <- c(members, paste0(gsub("[.]", "_", column), "_", genes))
  }
  members
}

compact_one_label <- function(label) {
  species <- sub("_(proteins|prot|pep)_.*$", "", label)
  gene_patterns <- c(
    "[A-Z]{2,}[A-Z0-9]*EN_[0-9]+",
    "XP_[0-9]+[.][0-9]+",
    "KAI[0-9]+[.][0-9]+",
    "KAG[0-9]+[.][0-9]+",
    "KAL[0-9]+[.][0-9]+",
    "ORZ[0-9]+[.][0-9]+",
    "ANIA_[0-9]+",
    "NCU[0-9]+",
    "SP[A-Z0-9]+[.][0-9]+",
    "UMAG_[0-9]+",
    "CNG[0-9]+",
    "CNA[0-9]+",
    "jgi[|]Muccir1_3[|][0-9]+"
  )
  gene <- NA_character_
  for (pattern in gene_patterns) {
    hit <- regmatches(label, regexpr(pattern, label, perl = TRUE))
    if (length(hit) == 1 && nzchar(hit)) {
      gene <- hit
      break
    }
  }
  if (is.na(gene)) {
    gene <- sub("^.*_(proteins|prot|pep)_", "", label)
    gene <- substr(gene, 1, 24)
  }
  paste(species, gene, sep = "|")
}

compact_label <- function(labels) {
  vapply(labels, compact_one_label, character(1))
}

member_matches <- function(tip_labels, members) {
  normalize_match <- function(value) {
    gsub("[^A-Za-z0-9]+", "", value)
  }
  fallback_member <- function(member) {
    parts <- strsplit(member, " ")[[1]]
    if (length(parts) >= 2 && grepl("_jgi[|]", parts[[1]])) {
      return(paste(parts[[1]], parts[[2]]))
    }
    parts[[1]]
  }
  normalized_tips <- normalize_match(tip_labels)
  matched <- rep(FALSE, length(tip_labels))
  for (member in members) {
    normalized_member <- normalize_match(member)
    member_hits <- grepl(normalized_member, normalized_tips, fixed = TRUE)
    if (!any(member_hits)) {
      fallback <- normalize_match(fallback_member(member))
      member_hits <- grepl(fallback, normalized_tips, fixed = TRUE)
    }
    matched <- matched | member_hits
  }
  tip_labels[matched]
}

plot_tree <- function(tree, tip_table, title, path_base, show_labels = TRUE, width = 8, height = 10) {
  point_cols <- c(
    "SexM-like HOG" = "#2c7fb8",
    "SexP-like HOG" = "#7b3294",
    "Other OG0000175 HMG" = "grey78"
  )
  candidate_cols <- c("MAT HMG query match" = "#d7301f", "Other" = "grey30")

  p <- ggtree(tree, size = 0.35) %<+% tip_table +
    geom_tippoint(aes(color = hog_group, shape = candidate), size = 2.1, alpha = 0.95) +
    scale_color_manual(values = point_cols, drop = FALSE) +
    scale_shape_manual(values = c("MAT HMG query match" = 17, "Other" = 16), drop = FALSE) +
    labs(title = title, color = "HOG", shape = "Candidate") +
    theme_tree2() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  if (show_labels) {
    p <- p + geom_tiplab(aes(label = plot_label), size = 2.4, align = FALSE, offset = 0.01)
  }

  ggsave(paste0(path_base, ".pdf"), p, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(path_base, ".png"), p, width = width, height = height, dpi = 300, limitsize = FALSE)
}

tree <- read.tree(tree_file)
sexm_members <- read_hog_members(hog_file, sexm_hog)
sexp_members <- read_hog_members(hog_file, sexp_hog)

sexm_tips <- member_matches(tree$tip.label, sexm_members)
sexp_tips <- member_matches(tree$tip.label, sexp_members)

tip_table <- data.frame(
  label = tree$tip.label,
  plot_label = compact_label(tree$tip.label),
  hog_group = "Other OG0000175 HMG",
  candidate = "Other",
  stringsAsFactors = FALSE
)
tip_table$hog_group[tip_table$label %in% sexm_tips] <- "SexM-like HOG"
tip_table$hog_group[tip_table$label %in% sexp_tips] <- "SexP-like HOG"
tip_table$candidate[
  grepl("GZMUCI1EN_006548|GZMUCI2EN_007105", tip_table$label)
] <- "MAT HMG query match"

message("Tree tips: ", length(tree$tip.label))
message(sexm_hog, " matched tips: ", length(sexm_tips))
message(sexp_hog, " matched tips: ", length(sexp_tips))

full_base <- file.path(output_dir, paste0("OG0000175_HMG_full_", tree_tag))
plot_tree(
  tree,
  tip_table,
  "OG0000175 HMG-box family, final OrthoFinder run",
  full_base,
  show_labels = FALSE,
  width = 8,
  height = 8
)

for (hog_id in c(sexm_hog, sexp_hog)) {
  keep <- if (hog_id == sexm_hog) sexm_tips else sexp_tips
  if (length(keep) < 2) {
    warning("Skipping ", hog_id, ": fewer than two matched tips")
    next
  }
  pruned <- keep.tip(tree, keep)
  pruned_table <- tip_table[tip_table$label %in% pruned$tip.label, , drop = FALSE]
  group_name <- if (hog_id == sexm_hog) "SexM-like HOG" else "SexP-like HOG"
  base <- file.path(output_dir, paste0(hog_id, "_", gsub("[^A-Za-z0-9]+", "_", group_name), "_", tree_tag))
  plot_tree(
    pruned,
    pruned_table,
    paste0(hog_id, " / ", group_name, " within OG0000175"),
    base,
    show_labels = TRUE,
    width = 10,
    height = max(4, length(keep) * 0.28)
  )
}

membership <- tip_table[tip_table$hog_group != "Other OG0000175 HMG" | tip_table$candidate != "Other", ]
write.table(
  membership,
  file = file.path(output_dir, paste0("OG0000175_HOG_membership_", tree_tag, ".tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
