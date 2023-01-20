#' ---
#' title: "Algae sucrose synthase putative gene of interest"
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

#' * GO dictionary
go_annot <- annotation %>% select(`Annotation GO ID`,`Annotation GO Term`,`Annotation GO Category`) %>% 
  rename(ID=`Annotation GO ID`,Term=`Annotation GO Term`,Category=`Annotation GO Category`) %>% 
  filter(!is.na(ID)) 

go_annot <- tibble(ID=unlist(str_split(go_annot$ID,"\\|"),use.names=FALSE),
                   Term=unlist(str_split(go_annot$Term,"\\|"),use.names=FALSE),
                   Category=unlist(str_split(go_annot$Category,"\\|"),use.names=FALSE)) %>% distinct()
#' * Expression
expression <- read_tsv(here("data/analysis/salmon/variance-stabilised_model-aware_gene-expression_data.tsv"),
                       show_col_types=FALSE)

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
dir.create(here("data/analysis/sucrose-synthase"),showWarnings=FALSE)

#' Reverse sample swap
expression %<>% rename(B2_2_24hrs_D="B2_2_12hrs_D",B2_2_12hrs_D="B2_2_24hrs_D")

#' # Gene of interest
#' ## Lookup
#' We use a combination of the Sequence Description and Gene Ontology (GO) term to fetch candidates
goi <- unlist(annotation[grepl("sucrose synthase",annotation$`Sequence Description`) | 
  grepl("sucrose synthase",annotation$`Annotation GO Term`),"Sequence Name"],use.names=FALSE)

annotation[annotation$`Sequence Name` %in% goi,"Annotation GO Term"]
annotation[annotation$`Sequence Name` %in% goi,"Sequence Description"]

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

#' * Export
write_tsv(vst %>% as_tibble() %>% rownames_to_column("ID"),
          here("data/analysis/sucrose-synthase/goi-vst-expression.tsv"))

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

#' Analysis that could be done

#' #' Cut the tree
#' pos <- cutree(as.hclust(hm$rowDendrogram),k=5)
#' 
#' #' First cluster
#' sel <- pos == 1 
#' 
#' hm <- heatmap.2(t(scale(t(vst[sel,]))),
#'                 Colv=FALSE, dendrogram="row",
#'                 distfun=pearson.dist,
#'                 hclustfun=function(X){hclust(X,method="ward.D2")},
#'                 labRow = NA,trace = "none",
#'                 labCol = sub("B2_2_","",colnames(vst)),
#'                 col=hpal)
#' 
#' #' transcripts in that cluster
#' write_tsv(annotation %>% filter(`Sequence Name` %in% rownames(vst[sel,])),file=here("data/analysis/cell-wall/cluster-12-120h-annotation.txt"))
#' line_plot(samples,vst,rownames(vst[sel,]))
#' ggsave(file=here("data/analysis/cell-wall/cluster-12-120h-annotation.pdf"),width=12,height=8)
#' 
#' #' Second cluster (24, 72, 124)
#' sel <- pos == 2
#' 
#' hm <- heatmap.2(t(scale(t(vst[sel,]))),
#'                 Colv=FALSE, dendrogram="row",
#'                 distfun=pearson.dist,
#'                 hclustfun=function(X){hclust(X,method="ward.D2")},
#'                 labRow = NA,trace = "none",
#'                 labCol = sub("B2_2_","",colnames(vst)),
#'                 col=hpal)
#' 
#' write_tsv(annotation %>% filter(`Sequence Name` %in% rownames(vst[sel,])),file=here("data/analysis/cell-wall/cluster-24-120h-annotation.txt"))
#' line_plot(samples,vst,rownames(vst[sel,]))
#' ggsave(file=here("data/analysis/cell-wall/cluster-24-120h-annotation.pdf"),width=12,height=8)
#' 
#' #' fourth cluster (4h)
#' sel <- pos == 4 
#' 
#' hm <- heatmap.2(t(scale(t(vst[sel,]))),
#'                 Colv=FALSE, dendrogram="row",
#'                 distfun=pearson.dist,
#'                 hclustfun=function(X){hclust(X,method="ward.D2")}, trace = "none",
#'                 labCol = sub("B2_2_","",colnames(vst)),
#'                 col=hpal,margins=c(5,10),cexRow=0.7)
#' 
#' write_tsv(annotation %>% filter(`Sequence Name` %in% rownames(vst[sel,])),file=here("data/analysis/cell-wall/cluster-4h-annotation.txt"))
#' line_plot(samples,vst,rownames(vst[sel,]))
#' ggsave(file=here("data/analysis/cell-wall/cluster-4h-annotation.pdf"),width=12,height=8)
#' 
#' #' third cluster (std)
#' sel <- pos == 3
#' hm <- heatmap.2(t(scale(t(vst[sel,]))),
#'                 Colv=FALSE, dendrogram="row",
#'                 distfun=pearson.dist,
#'                 hclustfun=function(X){hclust(X,method="ward.D2")}, trace = "none",
#'                 labCol = sub("B2_2_","",colnames(vst)),
#'                 col=hpal,margins=c(5,10),cexRow=0.7)
#' 
#' write_tsv(annotation %>% filter(`Sequence Name` %in% rownames(vst[sel,])),file=here("data/analysis/cell-wall/cluster-std-annotation.txt"))
#' line_plot(samples,vst,rownames(vst[sel,]))
#' ggsave(file=here("data/analysis/cell-wall/cluster-std-annotation.pdf"),width=12,height=8)
#' 
#' #' fifth cluster - aspecific, uninteresting
#' #'
#' #' # GO terms
#' #' Extract the go IDs from the annotation, associate them with the genes in the different clusters,
#' #' before summarising their occurences
#' goFreq <- 
#'   tibble(ID=names(pos),Cluster=factor(pos)) %>% left_join(
#'     annotation %>% 
#'       select(`Annotation GO ID`,`Sequence Name`) %>% 
#'       rename(ID=`Sequence Name`,
#'              GO_ID=`Annotation GO ID`) %>% 
#'       separate_rows(GO_ID,sep="\\|"),by="ID") %>% 
#'   group_by(Cluster,GO_ID) %>% dplyr::count()
#' 
#' #' convert to a matrix
#' mat <- goFreq %>% pivot_wider(names_from=Cluster,values_from=n) %>% column_to_rownames("GO_ID") %>% as.matrix()
#' 
#' #' Impute NAs to be 0
#' mat[is.na(mat)] <- 0
#' 
#' #' Plot the matrix for rows that have more than 5 occurrences
#' pdf(here("data/analysis/cell-wall/GO-terms-occurence-heatmap.pdf"),width=8,height=12)
#' heatmap.2(mat[rowSums(mat)>5,],trace="none",
#'           Colv=FALSE, dendrogram="row",cexRow=0.3,
#'           col=RColorBrewer::brewer.pal(9,"Reds"))
#' dev.off()
#' 
#' #' Extract the annotation for a GO_ID
#' #' 
#' #' _e.g._ GO:0005618
#' go_annot %>% filter(ID=="GO:0005618") %>% select(Term,Category)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```



