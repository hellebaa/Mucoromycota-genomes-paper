# Gene-order indices for every genome pair.
#
# Output
#   results/tables/pair_indices.tsv
#   results/cache/pair_indices/<a>__<b>.rds  (per-pair cache; delete to recompute)

source("analysis/synteny/R/functions.R")
cfg <- syn_config()
cache_dir <- file.path(cfg$cache, "pair_indices")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

genomes <- read_genome_table(cfg)
pairs <- read_tsv(file.path(cfg$tables, "pairs.tsv"), show_col_types = FALSE)

genome_data <- set_names(genomes$genome) %>% map(load_genome, cfg = cfg)

run_pair <- function(k) {
  p <- pairs[k, ]
  f <- file.path(cache_dir, paste0(p$a, "__", p$b, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  set.seed(cfg$seed + k)
  an <- build_anchors(genome_data[[p$a]], genome_data[[p$b]], p$of_a, p$of_b, cfg)
  res <- bind_cols(p %>% select(a, b), pair_indices(an, cfg))
  saveRDS(res, f)
  res
}

n_cores <- max(1L, parallel::detectCores() - 1L)
message("Computing indices for ", nrow(pairs), " pairs on ", n_cores, " cores")
res <- parallel::mclapply(seq_len(nrow(pairs)), run_pair, mc.cores = n_cores)
failed <- map_lgl(res, inherits, "try-error")
if (any(failed)) stop("Failed pairs: ", paste(which(failed), collapse = ", "),
                      "\n", res[failed][[1]])

pair_indices <- pairs %>% left_join(bind_rows(res), by = c("a", "b"))
write_tsv(pair_indices, file.path(cfg$tables, "pair_indices.tsv"))
message("Done: ", nrow(pair_indices), " pairs; median anchors per pair ",
        median(pair_indices$n_anchors))
