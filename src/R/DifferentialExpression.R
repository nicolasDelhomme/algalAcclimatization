#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme & Amit Bajhaiya"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Working directory

#' * Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(VennDiagram))

#' * Helper files
suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
"line_plot" <- function(dds,vst,gene_id){
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    return(
        ggplot(bind_cols(as.data.frame(colData(dds)),
                         melt(vst[sel,])),
               aes(x=Time,y=value,group=Time)) +
            geom_point() + geom_smooth() +
            scale_y_continuous(name="VST expression") + 
            ggtitle(label=paste("Expression for: ",gene_id))
    )
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds)){
    
    if(length(contrast)==1){
        res <- results(dds,name=contrast)
    } else {
        res <- results(dds,contrast=contrast)
    }
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj)
    
    if(verbose){
        message(sprintf("There are %s genes that are DE",sum(sel)))
    }
            
    if(export){
        if(!dir.exists(default_dir)){
            dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
        }
        write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
        write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
    }
    if(plot){
        heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                  distfun = pearson.dist,
                  hclustfun = function(X){hclust(X,method="ward.D2")},
                  trace="none",col=hpal,labRow = FALSE,
                  labCol=labels[sample_sel]
        )
    }
    return(rownames(res[sel,]))
}

#' # Data
#' * DESeq object
load("../analysis/salmon/dds.rda")

#' * Metadata
metadata <- read_rds(path="/mnt/picea/projects/algae/cfunk/algal-acclimatization/analysis/metadata.rds")

#' * Swap the samples
#' B2_2_12hrs_D and B2_2_24hrs_D need swapping
d12 <- which(dds$SampleID == "B2_2_12hrs_D")
d24 <- which(dds$SampleID == "B2_2_24hrs_D")
dds$Time[d12] <- "24hrs"
dds$Time[d24] <- "12hrs"

dds$Time <- factor(dds$Time,levels=c("std","60min","4hrs","12hrs","24hrs","72hrs","120hrs"))

#' ## Select
#' * only the complete proteins
tx.sel <- metadata$type=="complete" & ! is.na(metadata$type)

dds <- dds[tx.sel,]

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE,fitType='local')
vst <- assay(vsd)
vst <- vst - min(vst)
write_tsv(as.data.frame(vst) %>% rownames_to_column(),
          path=here("data/analysis/salmon/variance-stabilised_model-aware_gene-expression_data.tsv"))

#' ## Gene of interest
#' * 
#' This is not relevant for now and I selected the gene at random, but it looks cool :-D
# line_plot(dds,vst,"TRINITY_DN78_c3_g1_i1")

#' ## Differential Expression
dds <- DESeq(dds,fitType="local")

#' * Results
#' This shows all possible contrasts
resultsNames(dds)

h1.vs.std <- extract_results(dds,vst,"Time_60min_vs_std",
                                 default_prefix="60min_vs_std_",
                                 labels=dds$Time,
                                 sample_sel=dds$Time %in% c("std","60min"),
                             export=FALSE,plot=FALSE)

h4.vs.h1 <- extract_results(dds,vst,c("Time","4hrs","60min"),
                             default_prefix="4hrs_vs_60min_",
                             labels=dds$Time,
                             sample_sel=dds$Time %in% c("4hrs","60min"),
                            export=FALSE,plot=FALSE)

h12.vs.h4 <- extract_results(dds,vst,c("Time","12hrs","4hrs"),
                             default_prefix="12hrs_vs_4hrs_",
                             labels=dds$Time,
                             sample_sel=dds$Time %in% c("4hrs","12hrs"),
                             export=FALSE,plot=FALSE)

h24.vs.h12 <- extract_results(dds,vst,c("Time","24hrs","12hrs"),
                             default_prefix="24hrs_vs_12hrs_",
                             labels=dds$Time,
                             sample_sel=dds$Time %in% c("24hrs","12hrs"),
                             export=FALSE,plot=FALSE)

h72.vs.h24 <- extract_results(dds,vst,c("Time","72hrs","24hrs"),
                             default_prefix="72hrs_vs_24hrs_",
                             labels=dds$Time,
                             sample_sel=dds$Time %in% c("72hrs","24hrs"),
                             export=FALSE,plot=FALSE)

h120.vs.h72 <- extract_results(dds,vst,c("Time","120hrs","72hrs"),
                             default_prefix="120hrs_vs_72hrs_",
                             labels=dds$Time,
                             sample_sel=dds$Time %in% c("120hrs","72hrs"),
                             export=FALSE,plot=FALSE)


#' ### Barplot
barplot(sapply(list(H1vsStd=h1.vs.std,
             H4vsH1=h4.vs.h1,
             H12vsH4=h12.vs.h4,
             H24vsH12=h24.vs.h12,
             H72vs24=h72.vs.h24,
             H120vsH72=h120.vs.h72),length),las=2)

#' There are about 14 k genes affected at any time point            
length(Reduce(union,list(H1=h1.vs.std,
            H4=h4.vs.h1,
            H12=h12.vs.h4,
            H24=h24.vs.h12,
            H72=h72.vs.h24,
            H120=h120.vs.h72)))

#' ### Venn Diagram
#' Create a VennDiagram comparing the different Time Point
#' 
#' Many genes are common between at least 2 time points - about 30 k
grid.newpage()
grid.draw(venn.diagram(list(H1=h1.vs.std,
                            H4=h4.vs.h1,
                            H12=h12.vs.h4,
                            H24=h24.vs.h12,
                            H72=h72.vs.h24),
                       NULL,
                       fill=pal[1:5]))


#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


