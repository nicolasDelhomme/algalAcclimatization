library(Mfuzz)
library(vsn)

setwd("/mnt/picea/projects/algae/cfunk/algal-acclimatization/analysis/salmon")

source("~/Git/UPSCb/src/R/gopher.R")

# load the data
load("dds.rda")

# the metadata
metadata <- readRDS("../metadata.rds")

# select for onlye the complete protein coding transcripts
dds <- dds[metadata$type == "complete" & ! is.na(metadata$type),]

# transform the expression values to stabilise the variance
vst <- varianceStabilizingTransformation(dds,blind=FALSE,fitType="local")
vst <- assay(vst)
vst <- vst - min(vst)

meanSdPlot(vst)

# combine the biological replicates to keep only their mean values
res <- sapply(split.data.frame(t(vst),dds$Time),colMeans)

# sort by time
res <- res[,c(1,6,5,3,4,7,2)]

# do the clustering
eset <- ExpressionSet(res)

# filter for low variance
f.eset <- filter.std(eset,0.1)

# standardise the data (https://en.wikipedia.org/wiki/Standard_score)
s.eset <- standardise(f.eset)

# estimate the fuzziness parameter
m1 <- mestimate(s.eset)

# perform the clustering
cl <- mfuzz(s.eset,c=48,m=m1)

save(cl,s.eset,file="../mfuzz.rda")

goi <- names(cl$cluster)[cl$cluster == 5]

enr <- gopher(goi,background=rownames(s.eset),task="go")


