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
  library(pathview)
  library(S4Vectors)
  library(stringr)
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

#' * Expression
load(here("data/analysis/salmon/dds-sample-swap-corrected.rda"))
load(here("data/analysis/DE/vst-aware.rda"))
dds$Response <- factor(levels=c("control","acute","early","late"),
                       sub(".*hrs","late",
                           gsub("12hrs|^4hrs","early",
                                sub("60min","acute",
                                    sub("std","control",dds$Time)))))

scaledAvgExp <- t(scale(t(sapply(split.data.frame(t(vst),dds$Response),colMeans))))
scaledAvgExp[is.nan(scaledAvgExp)] <- 0
scaledAvgExp %<>% as.data.frame() %>% rownames_to_column("TxID") %>% as_tibble()

#' # Figures
#' 
#' ## Figure 3
#'
#' ### Lipids

lipidTerms <- tibble(Term=scan(here("doc/Algae-selected-lipid-GO-terms.tsv"),what="character",sep="\n"))

#' Data manipulation
#' 
#' 1. we find the genes that are annotated with the terms
lipidTerms$Tx <- lapply(lipidTerms$Term,function(trm,ant){ant$TxID[grep(trm,ant$Term)]},annot)

#' 2. we count them and remove terms that are not found (most likely Aman used the whole transcriptome and not the algae subset to look up the terms)
lipidTerms %<>% mutate(Ntx=elementNROWS(lipidTerms$Tx)) %>% filter(Ntx>0)

#' 3. we find how many are in the DEG lists
lipidTerms %<>% bind_cols(sapply(sapply(degList,"[[","ID"),function(d,p){sapply(lapply(p,"%in%",d),sum)},lipidTerms$Tx))

#' 4. We removed the absent ones
lipidTerms %<>% filter(rowSums(lipidTerms[,names(degList)])>0)

#' 5. We calculate a frequency
freq <- lipidTerms[,names(degList)] / lipidTerms$Ntx

names(freq) <- paste0("percent_",names(freq))

lipidTerms %<>% bind_cols(freq)

#' 6. We plot the frequency (0 to 0.5 to 1 - blue to white to red) as a heatmap, removed the
#' dendrograms and clustering. Add the numbers as cell labels.
gplots::heatmap.2(rbind(c(percent_acute=0,percent_early=0.5,percent_late=1),
                        as.matrix(lipidTerms[,c("percent_acute","percent_early","percent_late")])),
                  trace = "none",dendrogram="none",
                  Colv=FALSE,Rowv=FALSE,key=FALSE,
                  labCol = names(degList),margins=c(5,20),
                  labRow = c("colour key",lipidTerms$Term),
                  breaks = seq(0,1,0.01),
                  cexCol = 1,cexRow=1,notecol="black",
                  cellnote=rbind(c("0%","50%","100%"),
                                 as.matrix(lipidTerms[,c("acute","early","late")])),
                  col=hpal,rowsep=1,sepcolor="white")

write_tsv(lipidTerms %>% mutate(Tx=sapply(lipidTerms$Tx,paste,collapse="|")),
          here("data/analysis/lipid/lipid_go_term_summary.tsv"))

#' ### Carbohydrates
carbohydrateTerms <- tibble(Term=scan(here("doc/Algae-selected-carbohydrate-GO-terms.tsv"),what="character",sep="\n"))

#' Data manipulation
#' 
#' 1. we find the genes that are annotated with the terms
carbohydrateTerms$Tx <- lapply(carbohydrateTerms$Term,function(trm,ant){ant$TxID[grep(trm,ant$Term)]},annot)

#' 2. we count them and remove terms that are not found (most likely Aman used the whole transcriptome and not the algae subset to look up the terms)
carbohydrateTerms %<>% mutate(Ntx=elementNROWS(carbohydrateTerms$Tx)) %>% filter(Ntx>0)

#' 3. we find how many are in the DEG lists
carbohydrateTerms %<>% bind_cols(sapply(sapply(degList,"[[","ID"),function(d,p){sapply(lapply(p,"%in%",d),sum)},carbohydrateTerms$Tx))

#' 4. We removed the absent ones
carbohydrateTerms %<>% filter(rowSums(carbohydrateTerms[,names(degList)])>0)

#' 5. We calculate a frequency
freq <- carbohydrateTerms[,names(degList)] / carbohydrateTerms$Ntx

names(freq) <- paste0("percent_",names(freq))

carbohydrateTerms %<>% bind_cols(freq)

#' 6. We plot the frequency (0 to 0.5 to 1 - blue to white to red) as a heatmap, removed the
#' dendrograms and clustering. Add the numbers as cell labels.
gplots::heatmap.2(rbind(c(percent_acute=0,percent_early=0.5,percent_late=1),
                        as.matrix(carbohydrateTerms[,c("percent_acute","percent_early","percent_late")])),
                  trace = "none",dendrogram="none",
                  Colv=FALSE,Rowv=FALSE,key=FALSE,
                  labCol = names(degList),margins=c(5,20),
                  labRow = c("colour key",carbohydrateTerms$Term),
                  breaks = seq(0,1,0.01),
                  cexCol = 1,cexRow=1,notecol="black",
                  cellnote=rbind(c("0%","50%","100%"),
                                 as.matrix(carbohydrateTerms[,c("acute","early","late")])),
                  col=hpal,rowsep=1,sepcolor="white")

write_tsv(carbohydrateTerms %>% mutate(Tx=sapply(carbohydrateTerms$Tx,paste,collapse="|")),
          here("data/analysis/carbs/carbohydrates_go_term_summary.tsv"))

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

#' We remove entries that have only a 3rd level, no 4th level EC code
uEC <- uEC[elementNROWS(sapply(uEC,gregexpr,pattern="\\."))==3]

#' We can only get 10 entries at a time, so we iterate. It is very slow, several hours.
pthws <- do.call("c",
                    lapply(breakInChunks(totalsize = length(uEC),chunksize = 10),
                           function(rng,ecs){
                             message(sprintf("Running range starting at %s",rng[1]))
                             lst <- lapply(keggGet(ecs[rng]),"[[","PATHWAY")
                             names(lst) <- ecs[rng]
                             return(lst)
                           },uEC))
stopifnot(all(names(pthws) == uEC))

#' Report as a dictionary of pathways
upthw <- unlist(pthws)
names(upthw) <- sub(".*\\.","",names(upthw))
upthw <- upthw[!duplicated(upthw)]
pthw <- tibble(ID=names(upthw),
               Description=upthw)

#' Export the dictionary
write_tsv(pthw,here("data/analysis/annotation/kegg-pathways-dictionary.tsv"))

#' Combine the dictionary with the transcript and remove NA PID rows (EC code that have been made obsolete)
eCannot %<>% left_join(tibble(EC=sub("enzyme","EC",rep(names(pthws),elementNROWS(pthws))),
       PID=unlist(sapply(pthws,names),use.names=FALSE),
       PDESC=unlist(pthws,use.names=FALSE)
       ),by="EC",relationship = "many-to-many") %>% filter(!is.na(PID))

#' Export the kegg annotation
write_tsv(eCannot,here("data/analysis/annotation/kegg-transcript-pathway-mapping.tsv"))

#' Export the transcript annotation
kAnnot <- t(sapply(split.data.frame(as.data.frame(eCannot)[,-1],eCannot$TxID),
       function(df){apply(df,2,paste,collapse="|")})) %>% 
  as.data.frame() %>% rownames_to_column("TxID") %>% as_tibble()

write_tsv(kAnnot,here("data/analysis/annotation/algae_KEGG_annotation.tsv"))

#' ### Mock figure for fatty acids
#' 
#' 1. we find the terms
fAcids <- eCannot[grepl("fatty",eCannot$PDESC,ignore.case=TRUE),]

#' 2. we get the transcript lists and their length
fAcids %<>% group_by(PDESC) %>% summarise(TX=list(TxID)) %>% mutate(Ntx=elementNROWS(TX))

#' 3. we find how many are in the DEG lists
fAcids %<>% bind_cols(sapply(sapply(degList,"[[","ID"),function(d,p){sapply(lapply(p,"%in%",d),sum)},fAcids$TX))

#' 4. We calculate a frequency
freq <- fAcids[,names(degList)] / fAcids$Ntx

names(freq) <- paste0("percent_",names(freq))

fAcids %<>% bind_cols(freq)

#' 5. We plot the frequency (0 to 0.5 to 1 - blue to white to red) as a heatmap, removed the
#' dendrograms and clustering. Add the numbers as cell labels.
gplots::heatmap.2(rbind(c(percent_acute=0,percent_early=0.5,percent_late=1),
  as.matrix(fAcids[,c("percent_acute","percent_early","percent_late")])),
                  trace = "none",dendrogram="none",
                  Colv=FALSE,Rowv=FALSE,key=FALSE,
                  labCol = names(degList),margins=c(5,20),
                  labRow = c("colour key",fAcids$PDESC),
                  breaks = seq(0,1,0.01),
                  cexCol = 1,cexRow=1,notecol="black",
                  cellnote=rbind(c("0%","50%","100%"),
                                 as.matrix(fAcids[,c("acute","early","late")])),
                  col=hpal,rowsep=1,sepcolor="white")

#' 6. Write the results
dir.create(here("data/analysis/kegg"),showWarnings=FALSE)
write_tsv(fAcids %>% mutate(TX=sapply(fAcids$TX,paste,collapse="|")),
          here("data/analysis/kegg/fatty_acids_pathways_summary.tsv"))

#' ## Figure 6 (and 7)
#' 
#' Plotting the Fatty acid biosynthesis pathway: ec00061 - we actually use the one for KEGG orthologs (ko)
pspecies <- "ko"
pnumber <- "00061"

#' 1. we get the pathway information
pwy <- keggGet(paste0(pspecies,pnumber))

#' 2. Get the enzymes
#'
#' From the pathway
#' 
ko <- pwy[[1]]$ORTHOLOGY
ko <- ko[grepl("EC",ko)]
tb <- tibble(KO=names(ko),EC=ko) %>%
  mutate(EC=sub("]","",sub(".*EC:","",EC))) %>% 
  separate_longer_delim(EC," ")

#' From the annotation
ez <- eCannot %>% mutate(EC=str_sub(EC,4)) %>% 
  filter(EC %in% (tb %>% select(EC) %>% distinct() %>% unlist(use.names=FALSE)) & PID==paste0("ec",pnumber))

#' 3. Collate the expression values
ez %<>% left_join(scaledAvgExp,by="TxID")

#' 4. Summarise
ez %<>% left_join(tb,by="EC",relationship = "many-to-many") %>% 
  group_by(KO) %>% summarise(across(where(is_double),median))
gdat <- ez %>% column_to_rownames("KO") %>% as.matrix()

#' 5. plot
pathview(
  gene.data=gdat,
  species = pspecies,
  pathway.id = pnumber,
  out.suffix = "fatty-acid-biosynthesis",
  kegg.dir=here("data/analysis/kegg"))

#' 6. mv the result file in the result folder
outfile=list.files(here(),pattern=paste0(pspecies,pnumber,".*"))
file.rename(outfile,here("data/analysis/kegg",outfile))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
