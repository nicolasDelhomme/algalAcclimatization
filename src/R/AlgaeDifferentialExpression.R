#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' * Libraries
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/plotMA.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
    message(paste("Plotting",gene_id))
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                          data.frame(value=vst[sel,])),
                aes(x=Time,y=value,col=Time,group=Time)) +
        geom_point() + geom_smooth() +
        scale_y_continuous(name="VST expression") + 
        ggtitle(label=paste("Expression for: ",gene_id))
    
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE){
    
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
    
    # a look at independent filtering
    if(plot){
        plot(metadata(res)$filterNumRej,
             type="b", ylab="number of rejections",
             xlab="quantiles of filter")
        lines(metadata(res)$lo.fit, col="red")
        abline(v=metadata(res)$filterTheta)
    }
    
    if(verbose){
        message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                        round(metadata(res)$filterThreshold,digits=5),
                        names(metadata(res)$filterThreshold)))
    
        max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
        message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                        round(max.theta*100,digits=5),
                        round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
    }
    
    if(plot){
        qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
        dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                          qtl.exp=qtl.exp,
                          number.degs=sapply(lapply(qtl.exp,function(qe){
                              res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                                  ! is.na(res$padj) & res$baseMean >= qe
                          }),sum))
        
        plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
            geom_line() + geom_point() +
            scale_x_continuous("quantiles of expression") + 
            scale_y_continuous("base mean expression") +
            geom_hline(yintercept=expression_cutoff,
                       linetype="dotted",col="red"))
        
        if(debug){
            p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                geom_line() + geom_point() +
                scale_x_continuous("quantiles of expression") + 
                scale_y_log10("base mean expression") + 
                geom_hline(yintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
                     geom_line() + geom_point() +
                     geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes"))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Cumulative number of DE genes"))
            
            plot(ggplot(data.frame(x=dat$thetas[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("base mean of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
        }
        p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                          y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
            geom_line() + geom_point() +
            scale_x_log10("base mean of expression") + 
            scale_y_continuous("Number of DE genes per interval") + 
            geom_vline(xintercept=expression_cutoff,
                       linetype="dotted",col="red")
        suppressMessages(suppressWarnings(plot(p)))
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
        res$baseMean >= expression_cutoff
    
    if(verbose){
        message(sprintf("There are %s genes that are DE with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s",
                        sum(sel),
                        lfc,padj,expression_cutoff))
    }
    
    val <- rowSums(vst[sel,sample_sel])==0
    if (sum(val) >0){
        warning(sprintf("There are %s DE genes that have no vst expression in the selected samples",sum(val)))
        sel[sel][val] <- FALSE
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
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel & res$log2FoldChange > 0,]),
                dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' * Data
load(here("data/analysis/salmon/dds-sample-swap-corrected.rda"))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis/DE"),showWarnings=FALSE)
save(vst,file=here("data/analysis/DE/vst-aware.rda"))

#' ## Gene of interests
#' We do not have any goi at that time
# goi <- read_lines(here("doc/goi.txt"))
# stopifnot(all(goi %in% rownames(vst)))
# dev.null <- lapply(goi,line_plot,dds=dds,vst=vst)

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

#' ## Results
# fail
T1vsCtl <- extract_results(dds=dds,vst=vst,
                           contrast="Time_60min_vs_std",
                           default_prefix="DE-Time-1h-vs-std",
                           labels=paste(dds$Time,dds$Replicate,sep="-"),
                           sample_sel=dds$Time %in% c("60min","std"))
T4vs1 <- extract_results(dds=dds,vst=vst,
                         contrast=list("Time_4hrs_vs_std","Time_60min_vs_std"),
                         default_prefix="DE-Time-4h-vs-1h",
                         labels=paste(dds$Time,dds$Replicate,sep="-"),
                         sample_sel=dds$Time %in% c("4hrs","60min"))
T12vs4 <- extract_results(dds=dds,vst=vst,
                          contrast=list("Time_12hrs_vs_std","Time_4hrs_vs_std"),
                          default_prefix="DE-Time-12h-vs-4h",
                          labels=paste(dds$Time,dds$Replicate,sep="-"),
                          sample_sel=dds$Time %in% c("12hrs","4hrs"))
T24vs12 <- extract_results(dds=dds,vst=vst,contrast=list("Time_24hrs_vs_std","Time_12hrs_vs_std"),
                           default_prefix="DE-Time-24h-vs-12h",
                           labels=paste(dds$Time,dds$Replicate,sep="-"),
                           sample_sel=dds$Time %in% c("24hrs","12hrs"))
T72vs24 <- extract_results(dds=dds,vst=vst,contrast=list("Time_72hrs_vs_std","Time_24hrs_vs_std"),
                           default_prefix="DE-Time-72h-vs-24h",
                           labels=paste(dds$Time,dds$Replicate,sep="-"),
                           sample_sel=dds$Time %in% c("72hrs","24hrs"))
T120vs72 <- extract_results(dds=dds,vst=vst,contrast=list("Time_120hrs_vs_std","Time_72hrs_vs_std"),
                            default_prefix="DE-Time-120h-vs-72h",
                            labels=paste(dds$Time,dds$Replicate,sep="-"),
                            sample_sel=dds$Time %in% c("120hrs","72hrs"))

#' ### Venn Diagram
res.list <- list(T1vsCtl=T1vsCtl,
                 T4vs1=T4vs1,
                 T12vs4=T12vs4,
                 T24vs12=T24vs12,
                 T72vs24=T72vs24,
                 T120vs72=T120vs72)

#' #### Fast response
grid.newpage()
grid.draw(venn.diagram(lapply(res.list[1:2],"[[","all"),
                       NULL,
                       fill=pal[1:2]))

grid.newpage()
grid.draw(venn.diagram(lapply(res.list[1:2],"[[","up"),
                       NULL,
                       fill=pal[1:2]))

grid.newpage()
grid.draw(venn.diagram(lapply(res.list[1:2],"[[","dn"),
                       NULL,
                       fill=pal[1:2]))

#' #### Switch
grid.newpage()
grid.draw(venn.diagram(lapply(res.list[2:4],"[[","all"),
                       NULL,
                       fill=pal[1:3]))

grid.newpage()
grid.draw(venn.diagram(lapply(res.list[2:4],"[[","up"),
                       NULL,
                       fill=pal[1:3]))

grid.newpage()
grid.draw(venn.diagram(lapply(res.list[2:4],"[[","dn"),
                       NULL,
                       fill=pal[1:3]))

#' #### Acclimatization
grid.newpage()
grid.draw(venn.diagram(lapply(res.list[4:6],"[[","all"),
                       NULL,
                       fill=pal[1:3]))

grid.newpage()
grid.draw(venn.diagram(lapply(res.list[4:6],"[[","up"),
                       NULL,
                       fill=pal[1:3]))

grid.newpage()
grid.draw(venn.diagram(lapply(res.list[4:6],"[[","dn"),
                       NULL,
                       fill=pal[1:3]))

#' ### Gene Ontology enrichment
background <- rownames(vst)[featureSelect(vst,dds$Time,exp=0.1)]

enr.list <- lapply(res.list,function(r){
    lapply(r,gopher,background=background,task="go",url="algae")
})

dev.null <- lapply(names(enr.list),function(n){
    r <- enr.list[[n]]
    if(!is.null(r$all$go)){
        write_delim(r$all$go,path=file.path(file.path(here("data/analysis/DE",
                                                           paste0(n,"-all-DE-genes_GO-enrichment.txt")))))
        write_delim(r$all$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                            paste0(n,"-all-DE-genes_GO-enrichment_for-REVIGO.txt")))))
    }
    if(!is.null(r$up$go)){
        write_csv(r$up$go,path=file.path(file.path(here("data/analysis/DE",
                                                        paste0(n,"-up-DE-genes_GO-enrichment.txt")))))
        write_delim(r$up$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                           paste0(n,"-up-DE-genes_GO-enrichment_for-REVIGO.txt")))))    
    }
    if(!is.null(r$dn$go)){
        write_csv(r$dn$go,path=file.path(file.path(here("data/analysis/DE",
                                                        paste0(n,"-down-DE-genes_GO-enrichment.txt")))))
        write_delim(r$dn$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                           paste0(n,"-down-DE-genes_GO-enrichment_for-REVIGO.txt")))))    
    }
})

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


