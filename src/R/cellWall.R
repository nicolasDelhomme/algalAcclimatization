#' ---
#' title: "Algae cell wall putative gene of interest"
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
suppressPackageStartupMessages({
  library(tidyverse)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(magrittr)
  library(matrixStats)
})

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Annotation
annotation <- read_tsv(here("data/b2g/blast2go_20190117_export.txt.gz"),show_col_types=FALSE) %>% 
  mutate(`Sequence Name`=sub("\\.p\\d+$","",`Sequence Name`))

#' * Expression
expression <- read_tsv(here("data/analysis/salmon/variance-stabilised_model-aware_gene-expression_data.tsv"))

#' * Function
"line_plot" <- function(samples=samples,vst=vst,genes=genes){

  meds <- sapply(split.data.frame(scale(t(vst[genes,])),
                                  factor(samples$Time,levels=c("std","60min",paste0(c(4,12,24,72,120),"hrs")))),
                 colMedians)
  
  meanSd <- tibble(Means=colMeans(meds),Sd=colSds(meds),Timepoint=colnames(meds))
  
  p <- ggplot(meds %>% as_tibble() %>% mutate(geneID=genes) %>% 
                pivot_longer(1:ncol(meds),names_to="Timepoint") %>% 
                mutate(Timepoint=factor(Timepoint,levels=c("std","60min",paste0(c(4,12,24,72,120),"hrs")))),
              aes(x=Timepoint,y=value,col=Timepoint))+
    geom_violin() +
    geom_line(data=meanSd,aes(x=Timepoint,y=Means,group=1),col="darkblue") + 
    geom_line(data=meanSd,aes(x=Timepoint,y=Means+Sd,group=1),col="darkblue",linetype="dashed") + 
    geom_line(data=meanSd,aes(x=Timepoint,y=Means-Sd,group=1),col="darkblue",linetype="dashed") + 
    scale_y_continuous(name="VST expression") +
    ggtitle(label=paste("Z-score values per time point"))
  
  p
  suppressMessages(suppressWarnings(plot(p)))
  return(NULL)
}

#' * Output
dir.create(here("data/analysis/cell-wall"),showWarnings=FALSE)

#' Reverse sample swap
expression %<>% rename(B2_2_24hrs_D="B2_2_12hrs_D",B2_2_12hrs_D="B2_2_24hrs_D")

#' # Gene of interest
#' ## Lookup
#' We use a combination of the Sequence Description and Gene Ontology (GO) term to fetch candidates
goi <- unlist(annotation[grepl("cell wall",annotation$`Sequence Description`) | 
  grepl("cell wall",annotation$`Annotation GO Term`),"Sequence Name"],use.names=FALSE)

#' replace (\\.p\\d+$) .p followed by any number of digits. The $ the sign that indicates the end of the text
goi <- sub("\\.p\\d+$","",goi)

#' sanity
stopifnot(length(goi)==length(unique(goi)))

message(sprintf("There are %s candidates",length(goi)))

#' * samples information
samples <-read_tsv(here("doc/Samples-swap-corrected.tsv"),show_col_types = FALSE)

#' ## Expression
vst <- expression %>% filter(rowname %in% goi) %>% column_to_rownames() %>% as.matrix()
sel <- rowSums(vst) > 0
vst <- vst[sel,]

samples <- samples[match(colnames(vst),samples$SampleID),]

#' ### Heatmap
#' Clustered by rows and columns
hm <- heatmap.2(t(scale(t(vst))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = sub("B2_2_","",colnames(vst)),
                col=hpal)

#' Clustered by rows only
vst <- vst[,c(grep("std",colnames(vst)),
grep("60min",colnames(vst)),
grep("_4h",colnames(vst)),
grep("12h",colnames(vst)),
grep("24h",colnames(vst)),
grep("72h",colnames(vst)),
grep("120h",colnames(vst)))]

hm <- heatmap.2(t(scale(t(vst))),
                Colv=FALSE, dendrogram="row",
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = sub("B2_2_","",colnames(vst)),
                col=hpal)

#' Cut the tree
pos <- cutree(as.hclust(hm$rowDendrogram),k=5)

#' First cluster
sel <- pos == 1 

hm <- heatmap.2(t(scale(t(vst[sel,]))),
                Colv=FALSE, dendrogram="row",
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = sub("B2_2_","",colnames(vst)),
                col=hpal)

#' transcripts in that cluster
write_tsv(annotation %>% filter(`Sequence Name` %in% rownames(vst[sel,])),file=here("data/analysis/cell-wall/cluster-12-120h-annotation.txt"))
line_plot(samples,vst,rownames(vst[sel,]))
ggsave(file=here("data/analysis/cell-wall/cluster-12-120h-annotation.pdf"),width=12,height=8)

#' Second cluster (24, 72, 124)
sel <- pos == 2

hm <- heatmap.2(t(scale(t(vst[sel,]))),
                Colv=FALSE, dendrogram="row",
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = sub("B2_2_","",colnames(vst)),
                col=hpal)

write_tsv(annotation %>% filter(`Sequence Name` %in% rownames(vst[sel,])),file=here("data/analysis/cell-wall/cluster-24-120h-annotation.txt"))
line_plot(samples,vst,rownames(vst[sel,]))
ggsave(file=here("data/analysis/cell-wall/cluster-24-120h-annotation.pdf"),width=12,height=8)

#' fourth cluster (4h)
sel <- pos == 4 

hm <- heatmap.2(t(scale(t(vst[sel,]))),
                Colv=FALSE, dendrogram="row",
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")}, trace = "none",
                labCol = sub("B2_2_","",colnames(vst)),
                col=hpal,margins=c(5,10),cexRow=0.7)

write_tsv(annotation %>% filter(`Sequence Name` %in% rownames(vst[sel,])),file=here("data/analysis/cell-wall/cluster-4h-annotation.txt"))
line_plot(samples,vst,rownames(vst[sel,]))
ggsave(file=here("data/analysis/cell-wall/cluster-4h-annotation.pdf"),width=12,height=8)

#' third cluster (std)
sel <- pos == 3
hm <- heatmap.2(t(scale(t(vst[sel,]))),
                Colv=FALSE, dendrogram="row",
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")}, trace = "none",
                labCol = sub("B2_2_","",colnames(vst)),
                col=hpal,margins=c(5,10),cexRow=0.7)

write_tsv(annotation %>% filter(`Sequence Name` %in% rownames(vst[sel,])),file=here("data/analysis/cell-wall/cluster-std-annotation.txt"))
line_plot(samples,vst,rownames(vst[sel,]))
ggsave(file=here("data/analysis/cell-wall/cluster-std-annotation.pdf"),width=12,height=8)

#' fifth cluster - aspecific, uninteresting
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```



