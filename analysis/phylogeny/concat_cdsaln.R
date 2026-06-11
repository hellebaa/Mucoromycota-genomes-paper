.libPaths(c("~/R/library-4.3/","~/R/library-4.4/","~/R/x86_64-pc-linux-gnu-library/4.4/"))
library(tidyverse)
library(Biostrings)
cds_files <- list.files("/mnt/ScratchProjects/ebp-nor/Fungi/phylogeny/cds_aln_withOutgroups/",full.names = T)
all_cds <- lapply(cds_files, function(cds_file){
  cds <- readBStringSet(cds_file)
  species = gsub(".+_([^_]+)$","\\1",names(cds)) 
  sortSpecies <- sort(species)
  index <- data.frame("species"=sortSpecies) %>% 
      left_join(data.frame("order" = 1:length(cds),"species" = species))
  cds <- cds[index$order]
})
all_cds <- all_cds[!sapply(all_cds,is.null)]
all_cds <- do.call(xscat,all_cds)
names(all_cds) <- sortSpecies
writeXStringSet(all_cds,"/mnt/ScratchProjects/ebp-nor/Fungi/phylogeny/sco_cds_withOutgroups.fasta")