#' ---
#' title: "Differential Expression"
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
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(parallel)
    library(RColorBrewer)
    library(tidyverse)
    library(UpSetR)
    library(VennDiagram)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * parallel default
threads=1

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
                              debug=FALSE,filter=c("median",NULL),...){
    
    # get the filter
    if(!is.null(match.arg(filter))){
        filter <- rowMedians(counts(dds,normalized=TRUE))
        message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
    }
    
    # validation
    if(length(contrast)==1){
        res <- results(dds,name=contrast,filter = filter)
    } else {
        res <- results(dds,contrast=contrast,filter = filter)
    }
    
    stopifnot(length(sample_sel)==ncol(vst))
    
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
        if(debug){
            plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("base mean expression") +
                     geom_hline(yintercept=expression_cutoff,
                                linetype="dotted",col="red"))
        
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
            
            p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                geom_line() + geom_point() +
                scale_x_log10("base mean of expression") + 
                scale_y_continuous("Number of DE genes per interval") + 
                geom_vline(xintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
        }
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
        res$baseMean >= expression_cutoff
    
    if(verbose){
      message(sprintf(paste(
        ifelse(sum(sel)==1,
               "There is %s gene that is DE",
               "There are %s genes that are DE"),
        "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s"),
        sum(sel),padj,
        lfc,expression_cutoff))
    }
    
    # proceed only if there are DE genes
    if(sum(sel) > 0){
        val <- rowSums(vst[sel,sample_sel,drop=FALSE])==0
        if (sum(val) >0){
          warning(sprintf(paste(
            ifelse(sum(val)==1,
                   "There is %s DE gene that has",
                   "There are %s DE genes that have"),
            "no vst expression in the selected samples"),sum(val)))
          sel[sel][val] <- FALSE
        } 

        if(export){
            if(!dir.exists(default_dir)){
                dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
            }
            write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
            write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
        }
        if(plot & sum(sel)>1){
            heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                      distfun = pearson.dist,
                      hclustfun = function(X){hclust(X,method="ward.D2")},
                      trace="none",col=hpal,labRow = FALSE,
                      labCol=labels[sample_sel],...
            )
        }
    }
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel & res$log2FoldChange > 0,]),
                dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' 3. extract and plot the enrichment results
extractEnrichmentResults <- function(enrichment,task="go",
                                     diff.exp=c("all","up","dn"),
                                     go.namespace=c("BP","CC","MF"),
                                     genes=NULL,export=TRUE,plot=TRUE,
                                     default_dir=here("data/analysis/DE"),
                                     default_prefix="DE",
                                     url="athaliana"){
    # process args
    diff.exp <- match.arg(diff.exp)
    de <- ifelse(diff.exp=="all","none",
                 ifelse(diff.exp=="dn","down",diff.exp))

    # sanity
    if( is.null(enrichment[[task]]) | length(enrichment[[task]]) == 0){
        message(paste("No enrichment for",task))
    } else {

        # write out
        if(export){
            write_tsv(enrichment[[task]],
                      file=here(default_dir,
                                paste0(default_prefix,"-genes_GO-enrichment.tsv")))
            if(!is.null(genes)){
                write_tsv(
                    enrichedTermToGenes(genes=genes,terms=enrichment[[task]]$id,url=url,mc.cores=16L),
                    file=here(default_dir,
                              paste0(default_prefix,"-enriched-term-to-genes.tsv"))
                )
            }
        }
        
        if(plot){
            sapply(go.namespace,function(ns){
                titles <- c(BP="Biological Process",
                            CC="Cellular Component",
                            MF="Molecular Function")
                suppressWarnings(tryCatch({plotEnrichedTreemap(enrichment,enrichment=task,
                                                               namespace=ns,
                                                               de=de,title=paste(default_prefix,titles[ns]))},
                                          error = function(e) {
                                              message(paste("Treemap plot failed for",ns, 
                                                            "because of:",e))
                                              return(NULL)
                                          }))
            })
        }
    }
}

#' 4. get the background genes from the dds
getBgGenes <- function(dds,threads=1L){
  ifilt <- mclapply(resultsNames(dds),
                    function(nam,dds){
                      r<-results(dds,name=nam)
                      return(rownames(r)[!is.na(r$padj)])
                    },dds,mc.cores=threads)
  names(ifilt) <- resultsNames(dds)
  nobg <- Reduce(union,ifilt)
  upset(fromList(lapply(ifilt,function(i){which(nobg %in% i)})))
  return(nobg)
}

#' * Data
load(here("data/analysis/salmon/dds-sample-swap-corrected.rda"))

#' ## Update
#' * Creating the response variable
dds$Response <- factor(levels=c("control","acute","early","late"),
                       sub(".*hrs","late",
                           gsub("12hrs|^4hrs","early",
                                sub("60min","acute",
                                    sub("std","control",dds$Time)))))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis/DE"),showWarnings=FALSE)
save(vst,file=here("data/analysis/DE/vst-aware.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("data/analysis/DE/vst-aware.tsv"))

#' ## Gene of interests
#goi <- read_lines(here("doc/goi.txt"))
#stopifnot(all(goi %in% rownames(vst)))
#dev.null <- lapply(goi,line_plot,dds=dds,vst=vst)

#' ## Differential Expression
#' * Updating the model
design(dds) <- ~Response
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

#' ## Results
dir.create(here("data/analysis/DE/Response"))

#' ### Acute _vs._ control
acute_vs_control <- extract_results(dds,vst,"Response_acute_vs_control",
                                    default_dir=here("data/analysis/DE/Response"),
                                    default_prefix="acute_vs_control_",
                                    sample_sel=dds$Response %in% c("acute","control"),
                                    labels=dds$Response)


#' ### Early _vs._ control
early_vs_control <- extract_results(dds,vst,"Response_early_vs_control",
                                    default_dir=here("data/analysis/DE/Response"),
                                    default_prefix="early_vs_control_",
                                    sample_sel=dds$Response %in% c("early","control"),
                                    labels=dds$Response)

#' ### Late _vs._ control
late_vs_control <- extract_results(dds,vst,"Response_late_vs_control",
                                    default_dir=here("data/analysis/DE/Response"),
                                    default_prefix="late_vs_control_",
                                    sample_sel=dds$Response %in% c("late","control"),
                                    labels=dds$Response)


#' ### Venn Diagram
res.list <- list(acute=acute_vs_control,
                 early=early_vs_control,
                 late=late_vs_control)

#' #### All DE genes
grid.newpage()
grid.draw(venn.diagram(lapply(res.list,"[[","all"),
                      NULL,
                      fill=pal[1:3]))

#' #### DE genes (up in mutant)
grid.newpage()
grid.draw(venn.diagram(lapply(res.list,"[[","up"),
                      NULL,
                      fill=pal[1:3]))

#' #### DE genes (up in control)
grid.newpage()
grid.draw(venn.diagram(lapply(res.list,"[[","dn"),
                      NULL,
                      fill=pal[1:3]))

#' ### Gene Ontology enrichment
background <- getBgGenes(dds,threads=length(resultsNames(dds)))

enr.list <- lapply(res.list,function(r){
    lapply(r,gopher,background=background,task="go",url="algae")
})

dev.null <- lapply(names(enr.list),function(n){
    lapply(names(enr.list[[n]]),function(de){
        extractEnrichmentResults(enr.list[[n]][[de]],
                                 diff.exp=de,
                                 genes=res.list[[n]][[de]],
                                 default_prefix=paste(n,de,sep="-"),
                                 url="algae")
    })
})

#' ### Visualisation
#' #### Export
dev.null <- lapply(list.files(here("data/analysis/DE"),pattern="*-genes_GO-enrichment.tsv",full.names=TRUE),
                   function(fil){
                     write_tsv(read_tsv(fil,col_select=c("id","padj"),
                                        show_col_types=FALSE),
                               here(sub("\\.tsv","_for-REVIGO.tsv",fil)),
                               col_names=FALSE)
                   })

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


