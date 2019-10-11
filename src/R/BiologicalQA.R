#' ---
#' title: "Biological QA"
#' author: "Nicolas Delhomme && Amit Bajhaiya"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Working directory
setwd("/mnt/picea/projects/algae/cfunk/algal-acclimatization/Salmon")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/algae/cfunk/algal-acclimatization/Salmon")
#' ```

#' * Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))

#' * Helper functions
source("~/Git/UPSCb/src/R/featureSelection.R")

#' * Graphics
pal <- brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Metadata
#' Sample information
samples <- read_delim(file = "~/Git/UPSCb/projects/algal-acclimatization/doc/Samples.tsv", delim="\t") %>% 
  mutate(Time=relevel(factor(Time),"std")) %>% 
  mutate(Replicate=factor(Replicate))

#' # Raw data
filelist <- list.files(".", 
                    recursive = TRUE, 
                    pattern = "quant.sf",
                    full.names = TRUE)

#' Select the samples containing fungi
names(filelist) <- basename(dirname(filelist))

stopifnot(all(names(filelist) %in% samples$SampleID))

filelist <- filelist[match(samples$SampleID,names(filelist))]

#' Read the transcript expression 
tx <- suppressMessages(tximport(files = filelist, 
                                type = "salmon",
                                txOut=TRUE))

counts <- round(tx$counts)

#' ## Raw Data QC analysis
#' Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s transcripts are probable chimeraes",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' * Sequencing depth
#' The sequencing depth is  very similar between samples,
#' around 15M reads
ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
         bind_cols(samples),
       aes(x,y,fill=Time)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is shifted towards
#' low expression, which may be due to both technical (chimearaes)
#' and biological (_e.g._ bacteria) aspects
#' 
ggplot(melt(log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by condition. The gene coverage 
#' across samples is extremely similar
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Time=samples$Time[match(ind,samples$SampleID)])

#' * Color by Time
#' We observe no bias for any given time point and all
#' samples behave similarly
ggplot(dat,aes(x=values,group=ind,col=Time)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(file.path("..","analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file="../analysis/salmon/raw-unormalised-gene-expression_data.csv")

#' ## Data normalisation 
#' ### Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~ Time)

save(dds,file="../analysis/salmon/dds.rda")

#' Check the size factors (i.e. the sequencing library size effect)
#' 
#' The sequencing depth is not variable -/+ 25%
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
#' The transformation is done blind (ignoring the experimental design), 
#' as we are performing the quality assessment.
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
#' The normalised values are on an approx. log2 scale
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked but not optimally, as we 
#' have a lot of very lowly expressed (presumably artefactual)
#' transcripts. On the other hand, transcripts that are more 
#' highly expressed show a relatively constant variance.
#' There is also evidence that some transcripts display a large 
#' variance. If that variance is linked to the study design, we 
#' will have sufficient power to detect differential expression.
#' To assess that we can perform a Principal Component Analysis (PCA).
meanSdPlot(vst[rowSums(counts)>0,])

#' ## QC on the normalised data
#' ### PCA
#' The PCA is conducted by sample, hence the matrix transposition below.
pc <- prcomp(t(vst))
  
percent <- round(summary(pc)$importance[2,]*100)
  
#' ### 3 first dimensions
#' This looks interesting as the sample separate clearly Time
#' There may be a sample swap between 12 and 24 hours
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Time)],pch=19)
legend("topleft",
       fill=pal[1:nlevels(samples$Time)],
       legend=levels(samples$Time))
par(mar=mar)

#' ### 2D
#' In all likelihood, B2_2_12hrs_D and B2_2_24hrs_D have been sample-swapped
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samples)

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")


ggplotly(p) %>% layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
                       yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' ### Heatmap
#' Filter for noise
#' A cutoff at a VST value of 1 leaves about 15000 genes, which is adequate for the QA
sels <- rangeFeatureSelect(counts=vst,
                   conditions=samples$Time,
                   nrep=3)

#' Let's zoom in a bit
pander(sapply(sels,sum))

#' * Heatmap of "all" genes
#' We know that a lot of transcripts have low expression and may be artefactual.
#' 
#' As such, we select an higher cutoff (8)
ctf=8
 
#' Taking into account all the genes (above that noise threshold), the data
#' clearly cluster according to the time, in two main clusters, an early (std until 12h) and a 
#' late one (24 to 120h). The late cluster is relatively homogeneous with almost no change 
#' between 72 and 120h. Here, the sample swap is very visible.
#' 
heatmap.2(t(scale(t(vst[sels[[ctf]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = samples$Time,
          col=hpal)

#' ## Conclusion
#' The quality of the data is good. The PCA shows that the samples cluster bytime, 
#' however, there has been in all likelihood a sample swap. 
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' 
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
