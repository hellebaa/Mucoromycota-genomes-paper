.libPaths(c("~/R/library-4.3/","~/R/library-4.4/","~/R/x86_64-pc-linux-gnu-library/4.4/"))
library(tidyverse)
library(Biostrings)
N0 <- read_tsv("data/orthofinder/results/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
#outgroups <- c("UstMayd","CryNeof","SchPomb","AspNid","NeuCras","CanAlbi","SacCere","ConCoro")
tbl <- dplyr::select(N0,-`Gene Tree Parent Clade`) %>%
  gather(key="spc",value="ID",-HOG,-OG) %>%
  drop_na() %>% 
  mutate(ID=strsplit(ID,split = ", ")) %>% 
  unnest(cols = "ID") %>%
  mutate(ID = gsub("[|] ","|",ID)) %>%
  mutate(ID = gsub("^([^ ]+).*","\\1",ID)) %>%
  mutate(spc = gsub("_prot","",spc)) %>%
  mutate(spc = gsub(".pep","",spc)) %>%
  mutate(spc = gsub(".proteins","",spc))
  #filter(spc %in% outgroups == F) 
species <- unique(tbl$spc)
counts_HOG_spc <- dplyr::count(tbl,HOG,spc) %>%
  filter(n==1)
counts_HOG <- table(counts_HOG_spc$HOG)
select_HOG <-  names(counts_HOG[counts_HOG==length(species)])
SCO <- filter(tbl,HOG %in% select_HOG)

