.libPaths(c("~/R/library-4.3/","~/R/library-4.4/","~/R/x86_64-pc-linux-gnu-library/4.4/"))
library(readr)
library(magrittr)
library(tidyverse)
library(dplyr)
library(ape)
library(Biostrings)
options(stringsAsFactors = FALSE)

cds_fasta <- list.files("data/genomes/",pattern = "*cds*",recursive = T,full.names = T)
geneID_files <- list.files("data/genomes/",pattern = "_IDtbl.tsv",recursive = T,full.names = T)

geneIDs <- lapply(geneID_files, function(gf){
  read_tsv(gf) %>%
    mutate(spc = gsub("_IDtbl.tsv","",basename(gf)))
}) %>% bind_rows()

all_cds <- lapply(cds_fasta,function(sp){
  cds <- readBStringSet(sp)
  spc = gsub(".fna","",basename(sp))
  ID = gsub("[|] ","|",names(cds))
  ID = gsub("^([^ ]+).*","\\1",ID)
  r <- data.frame("spc"=spc,"ID"=ID,"Sequence"=cds)
  rownames(r) = NULL
  r
}) %>% bind_rows()


r <- left_join(SCO,all_cds)
b <- filter(r,!is.na(Sequence))
a <- filter(r,is.na(Sequence)) %>% select(-Sequence)
a <- left_join(a, select(geneIDs,geneID,proteinID,spc), by = c("ID"="proteinID","spc"))
a <- select(a,-ID) %>%
  rename("geneID"="ID")
test_a <- left_join(a,all_cds)
all_cds <- rbind(b,test_a)

split_orthogroup <- split(all_cds,all_cds$OG)
message(length(split_orthogroup))
cds_outFolder_macse <- "/mnt/ScratchProjects/ebp-nor/Fungi/phylogeny/cds_sequences_withOutgroups"
dir.create(cds_outFolder_macse, recursive = T)
for (i in 1:length(split_orthogroup)){
  seq <- BStringSet(split_orthogroup[[i]]$Sequence)
  names(seq) <- paste0(split_orthogroup[[i]]$ID,"_",split_orthogroup[[i]]$spc)
  if (length(seq)>1){
    og <- split_orthogroup[[i]]$HOG[1]
      writeXStringSet(seq,filepath = paste0(cds_outFolder_macse,"/",og,".fasta"))    
  }
}
