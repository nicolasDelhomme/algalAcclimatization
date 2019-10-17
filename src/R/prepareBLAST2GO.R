#' ---
#' title: "Extract the sequences for B2G"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' 
#' Load libraries
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(parallel))

#' Read the trinity transcript sequences
mRNA <- readDNAStringSet(here("data/trinity/Trinity.fasta"))

#' Modify the IDs (so as to be able to merge with the existing Transdecoder diamond results)
metadata <- readRDS(here("data/analysis/metadata.rds"))

stopifnot(all(sub(" .*","",names(mRNA)) == metadata$tID))

names(mRNA) <- sub(" .*","",names(mRNA))
names(mRNA)[!is.na(metadata$pID)] <- metadata$pID[!is.na(metadata$pID)]

#' # Export
#' Prepare the blast input
m_chk <- breakInChunks(length(mRNA),400)

dir.create(here("data/blastn/tmp"),recursive=TRUE,showWarnings=FALSE)
dev.null <- mclapply(1:length(m_chk),function(i){
  writeXStringSet(mRNA[start(m_chk)[i]:end(m_chk)[i]],
                  file=here(paste0("data/blastn/tmp/Trinity.fasta.",i)))  
},mc.cores=16L)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

