#' ---
#' title: "Algae cold response putative gene of interest"
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
  library(emoji)
  library(tidyverse)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(magrittr)
  library(matrixStats)
  library(pvclust)
})

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

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

#' # Collect 
#' 
#' ## Annotation
suppressWarnings(annotation <- read_tsv(here("data/b2g/blast2go_20190117_export.txt.gz"),show_col_types=FALSE) %>% 
  mutate(`Sequence Name`=sub("\\.p\\d+$","",`Sequence Name`)))

#' ## Terms of interest
go_annot <- annotation %>% select(`Annotation GO ID`,`Annotation GO Term`,`Annotation GO Category`) %>% 
  rename(ID=`Annotation GO ID`,Term=`Annotation GO Term`,Category=`Annotation GO Category`) %>% 
  filter(!is.na(ID))

go_annot <- tibble(ID=unlist(str_split(go_annot$ID,"\\|"),use.names=FALSE),
                   Term=unlist(str_split(go_annot$Term,"\\|"),use.names=FALSE),
                   Category=unlist(str_split(go_annot$Category,"\\|"),use.names=FALSE)) %>% distinct()

toi <- annotation[grepl("cold",annotation$`Annotation GO Term`),c("Annotation GO ID","Annotation GO Term")]

toi %<>% separate_longer_delim(everything(),delim="|") %>% distinct() %>% filter(grepl("cold",`Annotation GO Term`))

#' ## Gene of interest
goi <- unlist(annotation[Reduce("|",lapply(toi$`Annotation GO ID`,grepl,annotation$`Annotation GO ID`)),"Sequence Name"])

#' sanity
stopifnot(length(goi)==length(unique(goi)))
message(sprintf("There are %s candidates in the assembly",length(goi)))

#' ## Expression
#' ```{r IMPORTANT NOTE, echo=FALSE, eval=FALSE}
#' message("In the dds object with the corrected sample swap, the file name has no been changed! Only the metadata has been corrected (_i.e._ B2_2_12hrs_D has a Time 24hrs)")
#' ```
suppressPackageStartupMessages(load(here("data/analysis/DE/vst-aware.rda"),verbose=FALSE))
suppressPackageStartupMessages(load(here("data/analysis/salmon/dds-sample-swap-corrected.rda"),verbose=FALSE))
colnames(vst) <- dds$SampleID

#' # Visualise
outdir=here("data/analysis/cold-response")
dir.create(outdir,showWarnings=FALSE)

#' ## Expression
sel <- goi %in% rownames(vst)
message(sprintf("There are %s candidates in the algae",sum(sel)))
goi <- goi[sel]
expression <- vst[goi,]
sel <- rowSums(expression) > 0
message(sprintf("There are %s expressed genes",sum(sel)))
expression <- expression[sel,]

#' ### Heatmap
#' Clustered by rows and columns
hm <- heatmap.2(t(scale(t(expression))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = dds$Time,
                col=hpal)

hm.pvclust <- pvclust(data = t(scale(t(expression))),
                      method.hclust = "ward.D2", 
                      nboot = 100, parallel = TRUE)

plot(hm.pvclust, labels = dds$Time)
pvrect(hm.pvclust)

#' `r emoji("point_right")` **We observe a clear response between 1 and 12h, while all later time points are closest to standard, clearly highlighting an acclimation process.**
#'
#' Clustered by rows only
expression <- expression[,c(grep("std",colnames(expression)),
                            grep("60min",colnames(expression)),
                            grep("_4h",colnames(expression)),
                            grep("12h",colnames(expression)),
                            grep("24h",colnames(expression)),
                            grep("72h",colnames(expression)),
                            grep("120h",colnames(expression)))]

hm <- heatmap.2(t(scale(t(expression))),
                Colv=FALSE, dendrogram="row",
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = dds$Time[match(colnames(expression),dds$SampleID)],
                col=hpal)

#' `r emoji("point_right")` **We observe three clusters, one showing an early response, one showing a late response and one showing more diffuse differences. The latter could be Sam's non-responsive cold-response genes.**
#'
#' #### Subsets
#' 
#' To look at the three cluster of interests, we cut the tree
k <- 3
pos <- cutree(as.hclust(hm$rowDendrogram),k=k)

#' * Process the clusters
#' Note that a rerun can change the following order of the clusters:
#' 1. unresponsive genes
#' 2. acclimation
#' 3. acute response

dev.null <- sapply(1:k,function(i,h,d,a,e){
  sel <- pos == i
  heatmap.2(t(scale(t(e[sel,]))),
            Colv=FALSE, dendrogram="row",
            distfun=pearson.dist,
            hclustfun=function(X){hclust(X,method="ward.D2")},
            labRow = NA,trace = "none",
            labCol = d$Time[match(colnames(e),d$SampleID)],
            col=hpal)
  
  write_tsv(a %>% filter(`Sequence Name` %in% rownames(e[sel,])),
            file=here(outdir,paste0("cluster-",i,"-annotation.txt")))
  
  line_plot(colData(d)[match(colnames(e),d$SampleID),],
            e,rownames(e[sel,]))
  ggsave(file=here(outdir,paste0("cluster-",i,".pdf")),width=12,height=8)
  return(NULL)
},hm,dds,annotation,expression)

#' # GO terms
#' Extract the go IDs from the annotation, associate them with the genes in the different clusters,
#' before summarising their occurences
goFreq <-
  tibble(ID=names(pos),Cluster=factor(pos)) %>% left_join(
    annotation %>%
      select(`Annotation GO ID`,`Sequence Name`) %>%
      rename(`Annotation GO ID`="GOID",`Sequence Name`="ID") %>% 
      separate_longer_delim(GOID,delim="|") %>% 
      filter(GOID %in% toi$`Annotation GO ID`),by="ID") %>%
  group_by(Cluster,GOID) %>% dplyr::count()

#' convert to a matrix
goFreq %<>% pivot_wider(names_from=Cluster,values_from=n) %>% column_to_rownames("GOID")
goFreq %>% rownames_to_column("GOID") %>% left_join(toi %>% rename(`Annotation GO ID`="GOID"),by="GOID")

#' plot the matrix
pdf(here("data/analysis/cell-wall/GO-terms-occurence-heatmap.pdf"),width=8,height=12)
heatmap.2(as.matrix(goFreq),trace="none",
          Colv=FALSE, dendrogram="row",cexRow=0.8,
          col=RColorBrewer::brewer.pal(9,"Reds"))
dev.off()

#' ```{r extract,echo=FALSE,eval=FALSE}
#' # Extract the annotation for a GO_ID _e.g._ GO:0005618
#' go_annot %>% filter(ID=="GO:0005618") %>% select(Term,Category)
#' ```
#'
#' # Export
write_tsv(expression %>% as_tibble() %>% rownames_to_column("ID"),
          here(outdir,"goi-vst-expression.tsv"))

algae_annot <- annotation %>% 
  select(`Sequence Name`,`Annotation GO ID`,`Annotation GO Term`,`Annotation GO Category`) %>%
  filter(`Sequence Name` %in% rownames(vst)) %>% 
  dplyr::rename(TxID=`Sequence Name`,GOID=`Annotation GO ID`,
         Term=`Annotation GO Term`,Category=`Annotation GO Category`) %>% 
  filter(!is.na(GOID))

write_tsv(tibble(ID=unlist(str_split(algae_annot$GOID,"\\|"),use.names=FALSE),
                   Term=unlist(str_split(algae_annot$Term,"\\|"),use.names=FALSE),
                   Category=unlist(str_split(algae_annot$Category,"\\|"),use.names=FALSE)) %>% distinct(),
          here(outdir,"algae_GO_ID-to-Term.tsv"))

write_tsv(annotation %>% filter(`Sequence Name` %in% rownames(vst)),
          here(outdir,"algae_annotation.tsv"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
