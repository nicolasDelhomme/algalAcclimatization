#' ---
#' title: "Algae time series Biological QA"
#' subtitle: "projectID: u2019005"
#' author: "Nicolas Delhomme & Aman Zare"
#' date: "`r Sys.Date()`"
#' output: 
#'  html_document:
#'   fig_width: 9
#'   fig_height: 6
#'   toc: TRUE
#'   number_sections: TRUE
#'   toc_depth: 3
#'   toc_float:
#'    collapsed: TRUE
#'    smooth_scroll: TRUE
#'   code_folding: hide
#'   theme: "flatly" 
#'   highlight: pygments
#'   includes:
#'    before_body: header.html
#'    after_body: footer.html
#'   css: style.css    
#' ---
#' 
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(emoji)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(parallel)
  library(pander)
  library(pheatmap)
  library(plotly)
  library(pvclust)
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
samples <- read_delim(file = here("doc/Samples-swap-corrected.tsv"), 
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
levels =  c("std", "60min", "4hrs", "12hrs", "24hrs", "72hrs", "120hrs")
ggplot(dat,aes(x,y,fill=factor(Time, levels))) +
  geom_col() + 
  scale_y_continuous(name="reads") + 
  scale_fill_manual(values = pal) + 
  theme(axis.text.x=element_text(angle=90,size=10),
        axis.title.x=element_blank()) +
  labs(fill = "Time")

#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by Time 
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Time=samples$Time[match(ind,samples$SampleID)])
dat$Time <- factor(dat$Time, levels)

ggplot(dat,aes(x=values,group=ind,col=Time)) +
  geom_density(na.rm = TRUE) +
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  scale_color_manual(values = pal[1:7])

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
  design = ~ Time) %>%
  suppressMessages()

save(dds,file=here("data/analysis/salmon/dds-sample-swap-corrected.rda"))

#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
# pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked adequately
#'
#' * before: 
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])

#' * after:  
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
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)
  


#' `r emoji("star")` The PCA shows that a large fraction of the variance is explained by first 3 components. 

pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    samples)

#PC1 vs PC2
p1 <- pc.dat %>% 
  ggplot(mapping = aes(x=PC1,y=PC2,text=gsub(replacement = "", x = SampleID, pattern = "_B2_2"))) + 
  geom_point(size = 2) +
  aes(color=Time) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

p1 %<>% ggplotly() %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep=""))) %>% suppressWarnings()

#PC1 vs PC3
p2 <- pc.dat %>% 
  ggplot(mapping = aes(x=PC1,y=PC3,text=gsub(replacement = "", x = SampleID, pattern = "_B2_2"))) + 
  geom_point(size = 2) +
  aes(color=Time) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

p2 %<>% ggplotly() %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[3],"%)",sep=""))) %>% suppressWarnings()


#' ```{r subplot tr, out.width = '100%'}
#' subplot(style(p1, showlegend = F), p2,
#'         titleX = TRUE, titleY = TRUE, nrows = 1, margin = c(0.05, 0.05, 0, 0))
#' ```


#' ### Heatmap
#' 
#' Filter for noise
#' 
conds <- factor(dds$Time)
sels <- rangeFeatureSelect(counts=vst,
                           conditions=dds$Time,
                           nrep=3)
vst.cutoff <- 2

#' * Heatmap of "all" genes
#'
nn <- nrow(vst[sels[[vst.cutoff+1]],])
tn <- nrow(vst)
pn <- round(nn * 100/ tn, digits=1)



#'  `r emoji("warning")` **`r pn`%** (`r nn`) of total `r tn` genes are plotted below:

mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))

#annotation column
colnames(mat) <- gsub("B2_2_","",colnames(mat))
col_anno <- samples %>% select(Time) %>% as.data.frame()
col_anno$Time <- factor(col_anno$Time, levels)
rownames(col_anno) <- colnames(mat)

#annotation colors
col_anno_col = pal[1:7]
names(col_anno_col) <- levels
col_anno_col=list(Time = col_anno_col)

hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_cols = "correlation",
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               annotation_col = col_anno,
               show_rownames = F,
               labels_col = conds,
               angle_col = 90,
               annotation_colors = col_anno_col,
               legend = F)


plot(as.hclust(hm$tree_col),xlab="",sub="")

hm.pvclust <- suppressMessages(
  pvclust(data = t(scale(t(vst[sels[[vst.cutoff]],]))),
          method.hclust = "ward.D2", 
          nboot = 100, parallel = TRUE, quiet = T)
)

plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)

#' <details><summary>bootstrapping results as a table</summary>
#' ```{r bootstrapping results as a table transcript}
#' print(hm.pvclust, digits=3)
#' ```
#' </details>
#' 
#' <hr />
#' &nbsp;
#' 
#' # Summary
#' `r emoji("star")` The data looks very good after fixing the earlier samples swap. Both PCA and heatmap
#' show separation between different times except for 72 and 120 hours which are mixed together. 
#' In PCA plots, the first component explains the majority of variation (27%) and separates samples 
#' at 1,4 and 12 hours from the rest (std, 24,72 and 120 hrs). The second component explains
#' further variation (16%) separating early stages (std, 1 and 4 hours) from later stages. Together the PCA 
#' and Heatmap plots suggest that 4 and 12hrs conditions have the strongest effect. 
#' This is less in 60 min but still there is a clear effect. While 72hr and 120hr clustered separately
#' from other times, they do not show much difference between themselves. 
#' 
#' # Session Info
#' <details><summary>Session Info</summary>
#' ```{r session info}
#' sessionInfo()
#' ```
#' </details>
#'   
#' &nbsp;
#' 
