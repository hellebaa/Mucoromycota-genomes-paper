# Genome QC (contiguity, ID mapping) and the table of genome pairs with
# sequence and time divergence.
#
# Output
#   results/tables/genome_qc.tsv
#   results/tables/pairs.tsv

source("analysis/synteny/R/functions.R")
cfg <- syn_config()
dir.create(cfg$tables, recursive = TRUE, showWarnings = FALSE)

genomes <- read_genome_table(cfg)

# --- Genome QC --------------------------------------------------------------
# Gene-based contiguity: sequences carrying genes, gene N50 (number of genes
# on the sequence that brings the cumulative count to half of all genes) and
# fraction of genes on sequences with >= 100 genes.
of_ids <- read_tsv(file.path(cfg$data, "orthofinder", "results", "Orthogroups", "Orthogroups.tsv"),
                   col_types = cols(.default = "c"), n_max = Inf)

genome_qc <- genomes %>%
  mutate(qc = map2(genome, of_name, function(g, of) {
    gd <- load_genome(g, cfg)
    n <- sort(table(gd$genes$seq), decreasing = TRUE)
    of_genes <- of_ids[[of]] %>% discard(is.na) %>% str_split(", ") %>% unlist() %>% word(1)
    tibble(genes = sum(n),
           seqs = length(n),
           gene_n50 = as.integer(n[which(cumsum(n) >= sum(n) / 2)[1]]),
           frac_genes_seq100 = sum(n[n >= 100]) / sum(n),
           orthofinder_genes = length(of_genes),
           frac_mapped = mean(of_genes %in% gd$alias$id))
  })) %>%
  unnest(qc)

write_tsv(genome_qc, file.path(cfg$tables, "genome_qc.tsv"))
message("Genomes: ", nrow(genome_qc),
        "; min fraction of OrthoFinder genes mapped to coordinates: ",
        round(min(genome_qc$frac_mapped), 3))

# --- Pairs ------------------------------------------------------------------
# analysis_pairs(): all within-order pairs, plus all pairs of EBP-Nor genomes
# across orders (the same rule decides what download_data.sh fetches).
dist_fun <- pair_distance_table(genomes, cfg)
pairs <- analysis_pairs(genomes) %>%
  mutate(within = order_a == order_b,
         clade = if_else(within, as.character(order_a), "between"),
         both_ebpnor = source_a == "EBPNOR" & source_b == "EBPNOR",
         file = orthologue_file(of_a, of_b, cfg),
         file_exists = file.exists(file)) %>%
  bind_cols(dist_fun(.$a, .$b))

stopifnot(all(pairs$file_exists), !anyNA(pairs$d_sub))
write_tsv(pairs %>% select(-file), file.path(cfg$tables, "pairs.tsv"))
message("Pairs: ", nrow(pairs), " (", sum(pairs$within), " within-order)")
