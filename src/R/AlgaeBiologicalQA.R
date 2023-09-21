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
##' <hr />
#' &nbsp;
#' 
#' # Setup
#' This section and the next are relevant for reproducibility purposes. For results, please skip to section 3 (Quality Control)
#' 
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(emoji)
  library(gplots)
  library(gtools)
  library(here)
  library(hyperSpec)
  library(limma)
  library(magrittr)
  library(parallel)
  library(patchwork)
  library(PCAtools)
  library(pheatmap)
  library(plotly)
  library(pvclust)
  library(RColorBrewer)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)
pal <- brewer.pal(9,"Blues")

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

#' * add filelist to samples as a new column
samples %<>% mutate(Filenames = filelist)

#' Read the expression at the transcript level (we have no gene information)
txi <- suppressMessages(tximport(files = samples$Filenames,
                                 type = "salmon",
                                 txOut=TRUE))

#' Read the algae IDs
IDs <- scan(here("data/analysis/annotation/algae-IDs.txt"),
            what = "character")

#' Subset the data
txi[!names(txi) == "countsFromAbundance"] <- lapply(txi[!names(txi) == "countsFromAbundance"],function(l){l[IDs,]})
stopifnot(all(sapply(lapply(lapply(txi[!names(txi) == "countsFromAbundance"],rownames),"==",IDs),all)))

counts <- txi$counts
colnames(counts) <- samples$SampleID

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
  theme(axis.text.x=element_text(angle=90,size=10),
        axis.title.x=element_blank()) +
  labs(fill = "Time")

#' `r emoji("point_right")` **We observe almost no difference in the raw sequencing depth**
#'
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

#' `r emoji("point_right")` **The cumulative gene coverage is as expected**
#'
#' The same is done for the individual samples colored by Time 
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Time=samples$Time[match(ind,samples$SampleID)])
dat$Time <- factor(dat$Time, levels)

ggplot(dat,aes(x=values,group=ind,col=Time)) +
  geom_density(na.rm = TRUE) +
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' `r emoji("point_right")` **All samples have the same sequencing depth characteristics and there is no deviation when we look at the Time variable**
#' 
#' ## Export
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data.csv"))

#' 
#' <hr />
#' &nbsp;
#' 
#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples,
  design = ~ Time) %>%
  suppressMessages()

save(dds,file=here("data/analysis/salmon/dds-sample-swap-corrected.rda"))

colnames(dds) <- paste(samples$Time,samples$Replicate,sep="_")

#' Check the size factors (_i.e._ the sequencing library size effect)
dds <- estimateSizeFactors(dds) %>% 
  suppressMessages()

boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)

#' and without outliers:
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)

#' `r emoji("point_right")` **There is almost no differences in the libraries' size factors. They are all within +/- 30% of the average library size.**
#'
#' Assess whether there might be a difference in library size linked to a
#' given metadata
boxplot(split(t(normalizationFactors(dds)),dds$Time),las=2,
        main="Sequencing libraries size factor by Time",
        outline=FALSE)

#' `r emoji("point_right")` **The scaling factor distribution is not independent from the Time variable.**

plot(colMeans(normalizationFactors(dds)),
     log10(colSums(counts(dds))),ylab="log10 raw depth",
     xlab="average scaling factor",
     col=rainbow(n=nlevels(dds$Time))[as.integer(dds$Time)],pch=19)
legend("bottomright",fill=rainbow(n=nlevels(dds$Time)),
       legend=levels(dds$Time),cex=0.6)

#' `r emoji("point_right")` **The scaling factor appear linearly proportional to the sequencing depth for all samples but Time 4hrs. This might indicate that the number of genes expressed at 4hrs differs from the other samples.**
#' 
#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE) %>% suppressMessages()
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' let's look at standard deviations before and after VST normalization. 
#' This plot is to see whether there is a dependency of SD on the mean. 
#'
#' Before: 
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])

#' After:  
meanSdPlot(vst[rowSums(vst)>0,])

#' After VST normalization, the red line is almost horizontal which indicates
#' no dependency of variance on mean (homoscedastic).
#' 
#' `r emoji("point_right")` **We can conclude that the variance stabilization worked adequately**
#' 
#' <hr />
#' &nbsp;
#' 
#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Using PCAtools
p <- pca(vst,colData(dds))

#' ### Scree plot
#' 
#' We define the number of variable of the model
nvar <- 1

#' An the number of possible combinations
nlevel <- nlevels(dds$Time)

#' We devise the optimal number of components using two methods
horn <- suppressWarnings(parallelPCA(vst))
elbow <- findElbowPoint(p$variance)

#' We plot the percentage explained by different components and try to empirically assess whether
#' the observed number of components would be in agreement with our model's assumptions.
#' 
#' * the red line represents number of variables in the model  
#' * the orange line represents number of variable combinations.
#' * the black dotted, annotate lines represent the optimal number of components 
#' reported by the horn and elbow methods.
#' 
ggplot(tibble(x=1:length(percent),y=cumsum(percent),p=percent),aes(x=x,y=y)) +
  geom_line() + geom_col(aes(x,p)) + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component",breaks=1:length(percent),minor_breaks=NULL) + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=c(horn$n,elbow),colour="black",linetype="dotted",linewidth=0.5) +
  geom_label(aes(x = horn$n + 1, y = cumsum(percent)[horn$n],label = 'Horn', vjust = 1)) +
  geom_label(aes(x = elbow + 1, y = cumsum(percent)[elbow],label = 'Elbow', vjust = 1))

#' `r emoji("point_right")` **The first component explains 27% of the data variance. Both metrics, Horn and Elbow suggest that five or six components are those that are informative. Indeed the slope of the curve is fairly linear past PC7 and that would indicate that the remaining PCs only capture sample specific noise. While this is only empirical, the scree plot support having only few variables of importance in the dataset.**
#'

#' ### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    as.data.frame(colData(dds)))

#' * PC1 vs PC2
p1 <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p1) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,
       colby = "Time",
       colLegendTitle = "Time",
       lab=samples$Replicate,
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` **The first dimension separates the control and later samples, from the very early (1h centered) and early ones (4-12h).**
#'
#' * PC1 vs PC3
p2 <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=Time,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p2) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[3],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,x = 'PC1', y = 'PC3',
       colby = 'Time',
       colLegendTitle = 'Time',
       lab=samples$Replicate,
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' ```{r subplot tr, out.width = '100%'}
#' subplot(style(p1, showlegend = F), p2,
#'         titleX = TRUE, titleY = TRUE, nrows = 1, margin = c(0.05, 0.05, 0, 0))
#' ```
#' `r emoji("point_right")` **The third dimension separates the very late samples (72-120h) from the 24h ones**
#'
#' ### Pairs plot
#' This allows for looking at more dimensions, five by default
#' 
suppressMessages(pairsplot(p,colby='Time'))

#' `r emoji("point_right")` **All the first five PCs can be explained by an effect of the Time variable**
#'
#' ### Loadings
#' Loadings, _i.e._ the individual effect of every gene in a component can be studied. Here the most important ones are visualized throughout the different PCs
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1 to PC5',
             caption = 'Top 1% variables',
             drawConnectors = TRUE)

#' `r emoji("point_right")` **It could be interesting to lookup the annotations of some of these genes.**
#'
#' ### Correlation
#' This is a plot showing the correlation between the PC and the model variables. Note that while this might be relevant 
#' for a linear variable, it is less so for categorical variables. Sorting categorical variables in a linear order according to the PCs above might help.
#' 
p$metadata$Minutes <- as.integer(sub("std",0,sub("60min",1,sub("hrs","",samples$Time))))
p$metadata$Response <- factor(levels=c("control","acute","early","late"),
                              sub(".*hrs","late",
                                  gsub("12hrs|^4hrs","early",
                                       sub("60min","acute",
                                           sub("std","control",samples$Time)))))

suppressWarnings(eigencorplot(p,metavars=c('Minutes','Response')))

#' `r emoji("point_right")` **The first three components are associated with either variables.**
#'
#' ### Samples Distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- dds$SampleID
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=pal,labels_row=samples$Time,
         labels_col=samples$Time)

#' `r emoji("point_right")` **4h and 12h are the most distant time points, with 4h being the furthest away. 72 and 120 cannot be separated. 1h and 24h are close to the control, in that order.**
#' 
#' ## Sequencing depth
#' The figures show the number of genes expressed per condition at different expression cutoffs. The scale on the lower plot is the same as on the upper.
#' The first plot is a heatmap showing the number of genes above a given cutoff. The second plot shows it as a ratio of the number of genes expressed for (a)
#' given variable(s) divided by the average number of genes expressed at that cutoff across all variable(s). The latter plot is of course biased at higher cutoff 
#' as the number of genes becomes smaller and smaller.
#' The point of these two plots is to assert whether the number of genes expressed varies between conditions, as this would break some assumptions for normalisation and differential expression.
conds <- dds$Time
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

#' `r emoji("point_right")` **As suggested by the initial QC, 4hrs is somewhat of an outlier in the numbers of genes expressed.**
#"
#' Plotting the number of genes that are expressed (at any level)
do.call(rbind,split(t(nrow(vst) - colSums(vst==0)),samples$Time)) %>% as.data.frame() %>% 
  rownames_to_column("Time") %>% pivot_longer(starts_with("V")) %>% 
  ggplot(aes(x=Time, y=value,fill=Time)) + geom_dotplot(binaxis = "y", stackdir = "center") +
  scale_y_continuous("# of expressed genes")

#' `r emoji("point_right")` **There is a clear effect on both 4h and 12h. Whether that explains part of PC1 is to be kept in mind.**
#"
#' ### Heatmap
#' 
#' Here we want to visualise all the informative genes as a heatmap. We first filter the genes to remove those below the selected noise/signal cutoff. 
#' The method employed here is naive, and relies on observing a sharp decrease in the number of genes within the first few low level of expression. 
#' Using an independent filtering method, such as implemented in DESeq2 would be more accurate, but for the purpose of QA validation, a naive approach is sufficient.
#' Note that a sweet spot for computation is around 20 thousand genes, as building the hierarchical clustering for the heatmap scales non-linearly.
#'
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3) %>% 
  suppressWarnings()

#' `r emoji("point_right")` **Here a cutoff of 1 is selected**
vst.cutoff <- 1

nn <- nrow(vst[sels[[vst.cutoff+1]],])
tn <- nrow(vst)
pn <- round(nn * 100/ tn, digits=1)

#'  `r emoji("warning")` **`r pn`%** (`r nn`) of total `r tn` genes are plotted below:

mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))

#' annotation column
colnames(mat) <- gsub("B2_2_","",colnames(mat))
col_anno <- samples %>% select(Time) %>% as.data.frame()
col_anno$Time <- factor(col_anno$Time, levels)
rownames(col_anno) <- colnames(mat)

#' annotation colors
col_anno_col = brewer.pal(nlevels(conds),"Dark2")
names(col_anno_col) <- levels
col_anno_col=list(Time = col_anno_col)

hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_cols = "correlation",
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               annotation_col = col_anno,
               show_rownames = FALSE,
               labels_col = conds,
               angle_col = 90,
               annotation_colors = col_anno_col,
               legend = FALSE)

#' `r emoji("point_right")` **The heatmap shows a similar pattern to the PCA. Even if 4h could be biased as it has 10% less expressed genes, the changes observed are 4h are too drastic to be due to that alone.**
#'
#' ## Clustering of samples
#' ```{r echo=FALSE,eval=FALSE}
#' # Developer: This wouldonly works with the gplots heatmap.2, not the pheatmap
#' plot(as.hclust(hm$colDendrogram),xlab="",sub="")
#' ```
#'
#' Below we assess the previous dendrogram's reproducibility and plot the clusters with au and bp where:
#' 
#' * __au (Approximately Unbiased): computed by multiscale bootstrap resampling__ `r emoji("point_left")` a better approximation
#' * __bp (Bootstrap Probability): computed by normal bootstrap resampling__
#' 
#' `r emoji("warning")`Clusters with AU larger than 95% are highlighted by rectangles, which are strongly supported by data
#' 
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                      method.hclust = "ward.D2", 
                      nboot = 100, parallel = TRUE)

plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)

#' `r emoji("point_right")` **The clustering reveals three high level groups: std+60min, 4hrs+12hrs, 24hrs+72hrs_120hrs and six lower level groups, one per time point at the exception of 72hrs and 120hrs that are mixed.**

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
#' `r emoji("star")` The data looks very good after fixing the earlier samples swap. 
#' 
#' `r emoji("star")` Both PCA and heatmap show separation between different times 
#' except for 72 and 120 hours which are mixed together. 
#' 
#' `r emoji("star")` The 4h samples and to a lesser degree the 12h samples have up 
#' to 10% less genes expressed than the other samples. This could introduce a bias
#' in the expression quantification, albeit the changes observed in the heatmap suggests
#' that this effect has a technical origin rather than a biological origin: due to the higher 
#' expression levels at 4h, the sequencing depth would be shallower and lowly expressed genes 
#' would be drop-outs (while the gene is expressed in the sample, it is not being recorded)
#' 
#' `r emoji("star")` In PCA plots, the first component explains the majority of variation (27%) and separates samples 
#' at 1,4 and 12 hours from the rest (std, 24,72 and 120 hrs). The second component explains
#' further variation (16%) separating early stages (std, 1 and 4 hours) from later stages. 
#' 
#' `r emoji("star")` Together the PCA and heatmap plots suggest that 4 and 12hrs conditions have the strongest effect. 
#' This is less pronounced at 60 min but still there is nonetheless a clear effect. While 72hr and 120hr clustered separately
#' from other time points, they do not show much difference between themselves. 
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
