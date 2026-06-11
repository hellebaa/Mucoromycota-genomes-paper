.libPaths(c("~/R/library-4.3/","~/R/library-4.4/","~/R/x86_64-pc-linux-gnu-library/4.4/"))
library(magrittr)
library(tidyverse)

gffSplitAttributes <- function(tbl){
  tbl %>% 
    # Create a column with unique values outside the attribute column
    # or else the spread will not work.
    mutate(gffSplitAttributes_tmpidx = 1:n()) %>% 
    mutate(attribute = strsplit(attribute,split = ";(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)",perl = T)) %>% 
    unnest(attribute) %>%
    mutate(key = str_extract(attribute,"^[^=]+")) %>%
    mutate(attribute = gsub("^[^=]+=","",attribute)) %>%
    spread(key = key,value = attribute) %>% 
    select(-gffSplitAttributes_tmpidx) 
}

convertGFF <- function(gff){
  gffTbl <-
    read_tsv(gff,col_names = c("seqname","feature","start","end","strand","attribute"),comment="#", col_types="c-cii-c-c")
  
  geneTbl <- 
    gffTbl %>% 
    filter(feature == "gene") %>% 
    gffSplitAttributes() %>%
    select(seqname,start,end,strand,ID) %>% 
    dplyr::rename(geneID = ID)
  
  # get proteinID and parentID(to link with transcriptID or geneID) from CDS features
  CDSs <-
    gffTbl %>% 
    filter(feature == "CDS") %>% 
    select(attribute) %>% 
    # add temporary index
    mutate(idx = 1:n()) %>% 
    # split the attributes
    mutate(attribute = strsplit(attribute,split = ";")) %>% 
    unnest(cols = "attribute") %>% 
    separate(attribute,c("key","value"), sep="=") 
  if ("protein_id" %in% CDSs$key){
    CDSs <- filter(CDSs, key %in% c("protein_id", "Parent")) %>% 
      spread(key = key,value = value) %>%
      dplyr::select(-idx) %>%
      dplyr::rename(proteinID=protein_id) %>% 
      # Ther may be some CDSs with no proteinID. Don't include them.
      filter(!is.na(proteinID)) %>% 
      distinct()
  } else{
    CDSs <- filter(CDSs, key == "Parent") %>% 
      spread(key = key,value = value) %>%
      dplyr::select(-idx) %>%
      distinct() %>%
      mutate(proteinID=Parent)
  }
  
  features <- unique(gffTbl$feature)
  if ("mRNA" %in% features){
    mRNAs <-
      gffTbl %>% 
      filter(feature == "mRNA") %>% 
      select(attribute) %>%
      gffSplitAttributes() %>%
      select(ID,Parent) %>%
      dplyr::rename(rnaID = ID, geneID=Parent)
    IDtbl <- inner_join(CDSs %>% dplyr::rename(rnaID=Parent),mRNAs) 
  } else{
    IDtbl <- CDSs %>% dplyr::rename(geneID=Parent) %>%
      mutate(rnaID = proteinID)
  }
  return( 
    list(
      genePosTbl=geneTbl,
      IDtbl=IDtbl %>%
        select(proteinID,rnaID,geneID) 
    )
  )  
}

gffFiles <- list.files("data/genomes",pattern = "*.gff*",full.names = T,recursive = T)
for (f in gffFiles){
  r <- convertGFF(f)
  spc <- basename(dirname(f))
  write_tsv(r$IDtbl,file.path(dirname(f),paste0(spc,"_IDtbl.tsv")))
  write_tsv(r$genePosTbl,file.path(dirname(f),paste0(spc,"_genePosTbl.tsv")))
}
