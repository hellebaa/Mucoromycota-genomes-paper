# Genome table derived from the dated phylogeny, and the list of data files the
# analysis needs.
#
# Orders and the CENP-A / CENP-C losses are defined by tree nodes, each given as
# the MRCA of two tips. Every genome below an order's node belongs to that order;
# every genome below a loss node lacks the protein.
#
# Output
#   analysis/synteny/genomes.tsv
#   stdout: data files to download (gene tables per genome, orthologue tables
#           per analysed pair), one path per line, relative to the repo root

source("analysis/synteny/R/functions.R")
cfg <- syn_config()

order_nodes <- list(
  Mortierellales = c("gzPodHuml1", "ggMorZona1"),
  Umbelopsidales = c("gzUmbIsab1", "gzUmbRama1"),
  Mucorales      = c("SynRace", "gzMucPlum1")
)
cenpA_loss_node <- c("gzUmbIsab1", "gzMucPlum1")  # Umbelopsidales + Mucorales
cenpC_loss_node <- order_nodes$Mucorales
exclude <- c(MucCirc3 = paste("OrthoFinder used JGI protein IDs that are not in",
                              "MucCirc3_IDtbl.tsv; cannot map orthologues to coordinates"))

tree <- read.tree(cfg$tree_dated)
tips_below <- function(node) extract.clade(tree, getMRCA(tree, node))$tip.label
of_names <- scan(file.path(cfg$data, "orthofinder", "results", "Orthogroups", "Orthogroups.tsv"),
                 what = "", sep = "\t", nlines = 1, quiet = TRUE)[-1]
meta <- read_csv(file.path(cfg$data, "metadata", "genome_metadata_master.csv"),
                 show_col_types = FALSE)

genomes <- imap_dfr(order_nodes, ~ tibble(genome = tips_below(.x), order = .y)) %>%
  left_join(meta %>% select(genome = tree_label, species = species_name, family, source,
                            meta_order = order), by = "genome") %>%
  mutate(of_name = of_names[match(genome, sub("(\\.proteins|_prot)$", "", of_names))],
         cenpA = if_else(genome %in% tips_below(cenpA_loss_node), "absent", "present"),
         cenpC = if_else(genome %in% tips_below(cenpC_loss_node), "absent", "present"),
         include = if_else(genome %in% names(exclude), "no", "yes"),
         note = coalesce(exclude[genome], "")) %>%
  arrange(genome)

stopifnot(setequal(genomes$genome, tree$tip.label),
          !anyNA(genomes$of_name), !anyNA(genomes$species),
          genomes$order == genomes$meta_order)

genomes <- select(genomes, genome, of_name, species, family, order, source, cenpA, cenpC,
                  include, note)
write_tsv(genomes, cfg$genomes)
message("Wrote ", cfg$genomes, ": ", nrow(genomes), " genomes (",
        sum(genomes$include == "yes"), " included)")

included <- filter(genomes, include == "yes")
pairs <- analysis_pairs(included)
files <- c(
  file.path("data", "genomes", rep(included$genome, each = 2),
            paste0(rep(included$genome, each = 2), c("_genePosTbl.tsv", "_IDtbl.tsv"))),
  sub(paste0(cfg$repo, "/"), "", orthologue_file(pairs$of_a, pairs$of_b, cfg), fixed = TRUE)
)
writeLines(files)
