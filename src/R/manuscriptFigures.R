#' ---
#' title: "Manuscript figures"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
#' 
suppressPackageStartupMessages({
  library(here)
  library(KEGGREST)
  library(magrittr)
  library(S4Vectors)
  library(tidyverse)
})

#' * Palette
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * DE genes
degList <- lapply(sort(list.files(here("data/analysis/DE/Response"),pattern="*_genes.csv",full.names=TRUE)),
                  read_csv,show_col_types = FALSE,skip=1,
                  col_names=c("ID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

names(degList) <- c("acute","early","late")

#' * Annotation
annot <- read_tsv(here("data/analysis/cold-response/algae_GO_annotation.tsv"),
                  show_col_types=FALSE)

#' # Figures
#' 
#' ## Figure 3
#'
#' ### Lipids

lipidTerms <- tibble(Term=scan(here("doc/Algae-selected-lipid-GO-terms.tsv"),what="character",sep="\n"))

#' Data manipulation
#' 
#' 1. we find the genes that are annotated with the terms
lipidTerms$Position <- lapply(lipidTerms$Term,function(trm,ant){ant$TxID[grep(trm,ant$Term)]},annot)

#' 2. we count them and remove terms that are not found (most likely Aman used the whole transcriptome and not the algae subset to look up the terms)
lipidTerms %<>% mutate(Ntx=elementNROWS(lipidTerms$Position)) %>% filter(Ntx>0)

#' 3. we find how many are in the DEG lists
lipidTerms %<>% bind_cols(sapply(sapply(degList,"[[","ID"),function(d,p){sapply(lapply(p,"%in%",d),sum)},lipidTerms$Position))

#' 4. We removed the absent ones
lipidTerms %<>% filter(rowSums(lipidTerms[,names(degList)])>0)

#' 5. We calculate a frequency
freq <- lipidTerms[,names(degList)] / lipidTerms$Ntx

names(freq) <- paste0("percent_",names(freq))

lipidTerms %<>% bind_cols(freq)

#' 6. We plot the frequency (0 to 0.5 to 1 - blue to white to red) as a heatmap, removed the
#' dendrograms and clustering. Add the numbers as cell labels.
gplots::heatmap.2(as.matrix(lipidTerms[,c("percent_acute","percent_early","percent_late")]),
                  trace = "none",dendrogram="none",
                  Colv=FALSE,Rowv=FALSE,key=FALSE,
                  labCol = names(degList),margins=c(5,20),
                  labRow = lipidTerms$Term,
                  cexCol = 1,cexRow=1,notecol="black",
                  cellnote=as.matrix(lipidTerms[,c("acute","early","late")]),
                  col=hpal)

#' ### Carbohydrates
carbohydrateTerms <- tibble(Term=scan(here("doc/Algae-selected-carbohydrate-GO-terms.tsv"),what="character",sep="\n"))

#' Data manipulation
#' 
#' 1. we find the genes that are annotated with the terms
carbohydrateTerms$Position <- lapply(carbohydrateTerms$Term,function(trm,ant){ant$TxID[grep(trm,ant$Term)]},annot)

#' 2. we count them and remove terms that are not found (most likely Aman used the whole transcriptome and not the algae subset to look up the terms)
carbohydrateTerms %<>% mutate(Ntx=elementNROWS(carbohydrateTerms$Position)) %>% filter(Ntx>0)

#' 3. we find how many are in the DEG lists
carbohydrateTerms %<>% bind_cols(sapply(sapply(degList,"[[","ID"),function(d,p){sapply(lapply(p,"%in%",d),sum)},carbohydrateTerms$Position))

#' 4. We removed the absent ones
carbohydrateTerms %<>% filter(rowSums(carbohydrateTerms[,names(degList)])>0)

#' 5. We calculate a frequency
freq <- carbohydrateTerms[,names(degList)] / carbohydrateTerms$Ntx

names(freq) <- paste0("percent_",names(freq))

carbohydrateTerms %<>% bind_cols(freq)

#' 6. We plot the frequency (0 to 0.5 to 1 - blue to white to red) as a heatmap, removed the
#' dendrograms and clustering. Add the numbers as cell labels.
gplots::heatmap.2(as.matrix(carbohydrateTerms[,c("percent_acute","percent_early","percent_late")]),
                  trace = "none",dendrogram="none",
                  Colv=FALSE,Rowv=FALSE,key=FALSE,
                  labCol = names(degList),margins=c(5,20),
                  labRow = carbohydrateTerms$Term,
                  cexCol = 1,cexRow=1,notecol="black",
                  cellnote=as.matrix(carbohydrateTerms[,c("acute","early","late")]),
                  col=hpal)

#' ## Figure 4
#' 
#' ### Build the annotation
#' 
#' First we read the annotation and retrieve the transcript ID and enzyme code
eCannot <- read_tsv(here("data/analysis/cold-response/algae_annotation.tsv"),
                  show_col_types=FALSE) %>% 
  select(`Sequence Name`,`Enzyme Code`) %>%
  dplyr::rename(TxID=`Sequence Name`,EC=`Enzyme Code`) %>% 
  filter(!is.na(EC))

#' We convert to a longer format 
eCannot %<>% separate_longer_delim(EC,"|")

#' We use the REST API to fetch the corresponding pathways
uEC <- eCannot %>% select(EC) %>% mutate(EC=sub("EC","enzyme",EC)) %>% distinct() %>% unlist(use.names=FALSE)

#' We can only get 10 entries at a time, so we iterate. It is very slow, several hours.
pthws <- do.call("c",
                    lapply(breakInChunks(totalsize = length(uEC),chunksize = 10),
                           function(rng,ecs){
                             message(sprintf("Running range starting at %s",start(rng)))
                             lapply(keggGet(uEC[rng]),"[[","PATHWAY")
                           },uEC))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
