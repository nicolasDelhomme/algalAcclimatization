#' ---
#' title: "Algae time series Biological QA"
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
  library(parallel)
  library(pander)
  library(plotly)
  library(RColorBrewer)
  library(scatterplot3d)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
pal <- brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Metadata
#' Sample information
samples <- read_delim(file = here("doc/Samples.tsv"), 
                      delim="\t",
                      col_types=cols(
                        SampleID=col_factor(),
                        Time=col_factor(),
                        Replicate=col_factor())
                      ) %>% 
  mutate(Time=relevel(Time,"std"))

#' # Raw data
filelist <- list.files(here("data/salmon"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' Sort the raw data and samples 
samples <- samples[match(basename(dirname(filelist)),samples$SampleID),]
stopifnot(all(samples$SampleID == basename(dirname(filelist))))

#' name the file list vector
names(filelist) <- samples$SampleID

#' Read the expression at the transcript level (we have no gene information)
counts <- suppressMessages(round(tximport(files = filelist, 
                                  type = "salmon",txOut = TRUE)$counts))

#' Read the algae IDs
IDs <- scan(here("data/analysis/annotation/algae-IDs.txt"),
            what = "character")

#' Subset the data
counts <- counts[IDs,]

#' ## Quality Control
#' * Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' * Let us take a look at the sequencing depth, colouring by Time
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=Time)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),
        axis.title.x=element_blank())

#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by Time 
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Time=samples$Time[match(ind,samples$SampleID)])

ggplot(dat,aes(x=values,group=ind,col=Time)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~ Time)

save(dds,file=here("data/analysis/salmon/dds.rda"))

#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked adequately
#' 
meanSdPlot(vst[rowSums(vst)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model
nvar=1

#' An the number of possible combinations
nlevel=nlevels(dds$Time)

#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)
  
#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)

#' The PCA shows that a large fraction of the variance is 
#' explained by both variables.
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(dds$Time)],pch=19)
legend("topleft",
       fill=pal[1:nlevels(dds$Time)],
       legend=levels(dds$Time))
par(mar=mar)

#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samples)

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Replicate,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' ### Heatmap
#' 
#' Filter for noise
#' 
sels <- rangeFeatureSelect(counts=vst,
                           conditions=dds$Time,
                           nrep=3)
vst.cutoff <- 3

#' * Heatmap of "all" genes
#' 
hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = dds$Time,
          col=hpal)

plot(as.hclust(hm$colDendrogram),xlab="",sub="")

#' ## Conclusion
#' The data looks good. Clearly, there has been a sample swap 
#' between B2_2_24hrs_D and B2_2_12hrs_D.

D12p <- which(dds$SampleID == "B2_2_12hrs_D")
D24p <- which(dds$SampleID == "B2_2_24hrs_D")
colData(dds)[D12p,"Time"] <- "24hrs"
colData(dds)[D24p,"Time"] <- "12hrs"
samples[D12p,"Time"] <- "24hrs"
samples[D24p,"Time"] <- "12hrs"

save(dds,file=here("data/analysis/salmon/dds-sample-swap-corrected.rda"))
write_tsv(samples,path = here("doc/Samples-swap-corrected.tsv"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
