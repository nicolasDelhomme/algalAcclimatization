#' ---
#' title: "Assembly metadata creation"
#' author: "Nicolas Delhomme && Amit Bajhaiya"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(tidyverse))

#' # Data
#' ## Transdecoder
transIDs <- read_delim(here("data/TransDecoder/Trinity.fasta.transdecoder.pep.ID.txt"),
                       delim=" ",col_names=c("pID","tdID","ORF",
                                             "type","length","score","pos")) %>% 
    mutate(pID=sub(">","",pID),
           type=sub("type:","",type)) %>% 
    mutate(tID=sub("\\.p[0-9]+","",pID))

allIDs <- gsub("^>| .*","",scan(here("data/trinity/TrinityIDs.txt"),
               sep="\n",what="character"))

message(sprintf("Out of %s transcripts, %s are putatively protein-coding",
                length(allIDs),
                nrow(transIDs)))

barplot(table(transIDs$type))

#' Initiating the metadata
metadata <- tibble(tID=allIDs) %>% left_join(transIDs) %>% select(-ORF)

#' ## Diamond
#' Diamond was used to align the reads to UniRef90
#' 
#' The associated taxonomy IDs and Names are retrieved
diamond <- read_table2("data/diamond/uniref90.dmnd_Trinity.blt",
                       col_names = c("TRINITY_ID", "reference_ID", 
                                     "identity", "alignment_length", "mismatch",
                                     "gap_open", "start_trinity", "end_trinity", 
                                     "start_ref", "end_ref",
                                     "evalue", "bitscore", "trinity_length", 
                                     "ref_length", "taxonomy")) %>% 
  mutate(trinity_cov=alignment_length/trinity_length,
         ref_cov=alignment_length/ref_length) %>% 
  group_by(TRINITY_ID) %>% 
  arrange(desc(identity),desc(trinity_cov),desc(ref_cov)) %>% 
  slice(1) %>% ungroup()

message(sprintf("Out of %s protein coding transcripts, %s align to UniRef90",
                sum(!is.na(metadata$pID)),
                sum(!is.na(metadata$reference_ID))))

f<-function(lines,pos){
  return(lines[lines$V1 %in% diamond$reference_ID,])
}

tax_map <- read_delim_chunked(
  here("uniref/annotation/uniref90_id-table.txt"),
  chunk_size=1e6,delim=" ",callback=DataFrameCallback$new(f))

colnames(tax_map) <- c("reference_ID","TaxName","TaxID")

mar <- par("mar")
par(mar=c(12.1,4.1,0.1,0.1))
barplot(sort(table(tax_map$V2),decreasing=TRUE)[1:20],las=2,cex.names=.8)
par(mar=mar)

metadata %<>% left_join(left_join(diamond,tax_map,by="reference_ID"),
                        by=c("pID"="TRINITY_ID"))

message(sprintf("The %s UniRef90 alignment identify %s unique taxa",
                sum(!is.na(metadata$reference_ID)),
                length(unique(na.omit(metadata$TaxID)))
          ))

#' Read the taxonomy sql to report the Division


#' Add expression and GC content

#' Add time point specificity

#' # Export
write_rds(metadata,path=here("data/analysis/metadata.rds"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
