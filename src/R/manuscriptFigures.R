#' ---
#' title: "Manuscript figures"
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
#' 
suppressPackageStartupMessages({
  library(here)
  library(KEGGREST)
  library(magrittr)
  library(pathview)
  library(S4Vectors)
  library(stringr)
  library(tidyverse)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Palette
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * DE genes
degList <- lapply(sort(list.files(here("data/analysis/DE/Response"),pattern="*_genes.csv",full.names=TRUE)),
                  read_csv,show_col_types = FALSE,skip=1,
                  col_names=c("ID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

names(degList) <- c("acute","early","late")

degListPerTimePoint <- lapply(sort(list.files(here("data/analysis/DE"),pattern="*vs_std_genes.csv",full.names=TRUE)),
									read_csv,show_col_types = FALSE,skip=1,
									col_names=c("ID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))

names(degListPerTimePoint) <- c("120hrs_vs_ctrl","12hrs_vs_ctrl","24hrs_vs_ctrl","4hrs_vs_ctrl","1hr_vs_ctrl","72hrs_vs_ctrl")

#' * Annotation
annot <- read_tsv(here("data/analysis/cold-response/algae_GO_annotation.tsv"),
                  show_col_types=FALSE)

#' * Expression
load(here("data/analysis/salmon/dds-sample-swap-corrected.rda"))
load(here("data/analysis/DE/vst-aware.rda"))
dds$Response <- factor(levels=c("control","acute","early","late"),
                       sub(".*hrs","late",
                           gsub("12hrs|^4hrs","early",
                                sub("60min","acute",
                                    sub("std","control",dds$Time)))))

dds$RenameTime <- factor(levels=c("control","1hr","4hrs","12hrs","24hrs","72hrs","120hrs"),
                         sub("60min","1hr",
                             sub("std","control",dds$Time)))

scaledAvgExp <- t(scale(t(sapply(split.data.frame(t(vst),dds$Response),colMeans))))
scaledAvgExp[is.nan(scaledAvgExp)] <- 0
scaledAvgExp %<>% as.data.frame() %>% rownames_to_column("TxID") %>% as_tibble()

#' # Figures
#' 
#' ## Figure 2
#' 
conds <- dds$RenameTime
sels <- rangeFeatureSelect(counts=vst,
													 conditions=conds,
													 nrep=3,
													 plot=FALSE)
vst.cutoff <- 1
mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))
colnames(mat) <- sub("60min","1hr",sub("std","control",paste0(dds$Time,"_",dds$Replicate)))
#' Rearrange the column
mat <- mat[,c("control_A","control_D","control_B","control_C",
							"1hr_D","1hr_A","1hr_B","1hr_C",
							"4hrs_D","4hrs_B","4hrs_A","4hrs_C",
							"12hrs_D","12hrs_B","12hrs_A","12hrs_C",
							"24hrs_D","24hrs_B","24hrs_A","24hrs_C",
							"72hrs_D","72hrs_C","72hrs_B","72hrs_A",
							"120hrs_C","120hrs_D","120hrs_A","120hrs_B")]

#' annotation column
col_anno <- str_sub(colnames(mat),end = -3) %>% as.data.frame()
colnames(col_anno) <- "Time"
col_anno$Time <- factor(col_anno$Time, levels=c("control","1hr","4hrs","12hrs","24hrs","72hrs","120hrs"))
rownames(col_anno) <- colnames(mat)

#' annotation colors
col_anno_col = brewer.pal(nlevels(conds),"Dark2")
names(col_anno_col) <- c("control","1hr","4hrs","12hrs","24hrs","72hrs","120hrs")
col_anno_col=list(Time = col_anno_col)

#' ### Heatmap
hm2 <- pheatmap(mat,
								color = hpal,
								border_color = NA,
								cluster_cols = FALSE,
								clustering_distance_rows = "correlation",
								clustering_method = "ward.D2",
								annotation_col = col_anno,
								show_rownames = FALSE,
								labels_col = col_anno$Time,
								angle_col = 90,
								annotation_colors = col_anno_col,
								legend = FALSE)

#' ## Figure 3
#'
#' ### Lipids

lipidTerms <- tibble(Term=scan(here("doc/Algae-selected-lipid-GO-terms.tsv"),what="character",sep="\n"))

#' Data manipulation
#' 
#' 1. we find the genes that are annotated with the terms
lipidTerms$Tx <- lapply(lipidTerms$Term,function(trm,ant){ant$TxID[grep(trm,ant$Term)]},annot)

#' 2. we count them and remove terms that are not found (most likely Aman used the whole transcriptome and not the algae subset to look up the terms)
lipidTerms %<>% mutate(Ntx=elementNROWS(lipidTerms$Tx)) %>% filter(Ntx>0)

#' 3. we find how many are in the DEG lists
lipidTerms %<>% bind_cols(sapply(sapply(degList,"[[","ID"),function(d,p){sapply(lapply(p,"%in%",d),sum)},lipidTerms$Tx))

#' 4. We removed the absent ones
lipidTerms %<>% filter(rowSums(lipidTerms[,names(degList)])>0)

#' 5. We calculate a frequency
freq <- lipidTerms[,names(degList)] / lipidTerms$Ntx

names(freq) <- paste0("percent_",names(freq))

lipidTerms %<>% bind_cols(freq)

#' 6. We plot the frequency (0 to 0.5 to 1 - blue to white to red) as a heatmap, removed the
#' dendrograms and clustering. Add the numbers as cell labels.
gplots::heatmap.2(rbind(c(percent_acute=0,percent_early=0.5,percent_late=1),
                        as.matrix(lipidTerms[,c("percent_acute","percent_early","percent_late")])),
                  trace = "none",dendrogram="none",
                  Colv=FALSE,Rowv=FALSE,key=FALSE,
                  labCol = names(degList),margins=c(5,20),
                  labRow = c("colour key",lipidTerms$Term),
                  breaks = seq(0,1,0.01),
                  cexCol = 1,cexRow=1,notecol="black",
                  cellnote=rbind(c("0%","50%","100%"),
                                 as.matrix(lipidTerms[,c("acute","early","late")])),
                  col=hpal,rowsep=1,sepcolor="white")

write_tsv(lipidTerms %>% mutate(Tx=sapply(lipidTerms$Tx,paste,collapse="|")),
          here("data/analysis/lipid/lipid_go_term_summary.tsv"))

#' ### Carbohydrates
carbohydrateTerms <- tibble(Term=scan(here("doc/Algae-selected-carbohydrate-GO-terms.tsv"),what="character",sep="\n"))

#' Data manipulation
#' 
#' 1. we find the genes that are annotated with the terms
carbohydrateTerms$Tx <- lapply(carbohydrateTerms$Term,function(trm,ant){ant$TxID[grep(trm,ant$Term)]},annot)

#' 2. we count them and remove terms that are not found (most likely Aman used the whole transcriptome and not the algae subset to look up the terms)
carbohydrateTerms %<>% mutate(Ntx=elementNROWS(carbohydrateTerms$Tx)) %>% filter(Ntx>0)

#' 3. we find how many are in the DEG lists
carbohydrateTerms %<>% bind_cols(sapply(sapply(degList,"[[","ID"),function(d,p){sapply(lapply(p,"%in%",d),sum)},carbohydrateTerms$Tx))

#' 4. We removed the absent ones
carbohydrateTerms %<>% filter(rowSums(carbohydrateTerms[,names(degList)])>0)

#' 5. We calculate a frequency
freq <- carbohydrateTerms[,names(degList)] / carbohydrateTerms$Ntx

names(freq) <- paste0("percent_",names(freq))

carbohydrateTerms %<>% bind_cols(freq)

#' 6. We plot the frequency (0 to 0.5 to 1 - blue to white to red) as a heatmap, removed the
#' dendrograms and clustering. Add the numbers as cell labels.
gplots::heatmap.2(rbind(c(percent_acute=0,percent_early=0.5,percent_late=1),
                        as.matrix(carbohydrateTerms[,c("percent_acute","percent_early","percent_late")])),
                  trace = "none",dendrogram="none",
                  Colv=FALSE,Rowv=FALSE,key=FALSE,
                  labCol = names(degList),margins=c(5,20),
                  labRow = c("colour key",carbohydrateTerms$Term),
                  breaks = seq(0,1,0.01),
                  cexCol = 1,cexRow=1,notecol="black",
                  cellnote=rbind(c("0%","50%","100%"),
                                 as.matrix(carbohydrateTerms[,c("acute","early","late")])),
                  col=hpal,rowsep=1,sepcolor="white")

write_tsv(carbohydrateTerms %>% mutate(Tx=sapply(carbohydrateTerms$Tx,paste,collapse="|")),
          here("data/analysis/carbs/carbohydrates_go_term_summary.tsv"))

#' ## Figure 4
#' 
#' ### Build the annotation
#' 
#' First we read the annotation and retrieve the transcript ID and enzyme code
eCannot <- read_tsv(here("data/analysis/cold-response/algae_annotation.tsv"),
                  show_col_types=FALSE) %>% 
  select(`Sequence Name`,`Enzyme Code`) %>%
  dplyr::rename(TxID=`Sequence Name`,EC=`Enzyme Code`) %>% 
  filter(!is.na(EC))

#' We convert to a longer format 
eCannot %<>% separate_longer_delim(EC,"|")

#' We use the REST API to fetch the corresponding pathways
uEC <- eCannot %>% select(EC) %>% mutate(EC=sub("EC","enzyme",EC)) %>% distinct() %>% unlist(use.names=FALSE)

#' We remove entries that have only a 3rd level, no 4th level EC code
uEC <- uEC[elementNROWS(sapply(uEC,gregexpr,pattern="\\."))==3]

#' We can only get 10 entries at a time, so we iterate. It is very slow, several hours.
pthws <- do.call("c",
                    lapply(breakInChunks(totalsize = length(uEC),chunksize = 10),
                           function(rng,ecs){
                             message(sprintf("Running range starting at %s",rng[1]))
                             lst <- lapply(keggGet(ecs[rng]),"[[","PATHWAY")
                             names(lst) <- ecs[rng]
                             return(lst)
                           },uEC))
stopifnot(all(names(pthws) == uEC))

#' Report as a dictionary of pathways
upthw <- unlist(pthws)
names(upthw) <- sub(".*\\.","",names(upthw))
upthw <- upthw[!duplicated(upthw)]
pthw <- tibble(ID=names(upthw),
               Description=upthw)

#' Export the dictionary
write_tsv(pthw,here("data/analysis/annotation/kegg-pathways-dictionary.tsv"))

#' Combine the dictionary with the transcript and remove NA PID rows (EC code that have been made obsolete)
eCannot %<>% left_join(tibble(EC=sub("enzyme","EC",rep(names(pthws),elementNROWS(pthws))),
       PID=unlist(sapply(pthws,names),use.names=FALSE),
       PDESC=unlist(pthws,use.names=FALSE)
       ),by="EC",relationship = "many-to-many") %>% filter(!is.na(PID))

#' Export the kegg annotation
write_tsv(eCannot,here("data/analysis/annotation/kegg-transcript-pathway-mapping.tsv"))

#' Export the transcript annotation
kAnnot <- t(sapply(split.data.frame(as.data.frame(eCannot)[,-1],eCannot$TxID),
       function(df){apply(df,2,paste,collapse="|")})) %>% 
  as.data.frame() %>% rownames_to_column("TxID") %>% as_tibble()

write_tsv(kAnnot,here("data/analysis/annotation/algae_KEGG_annotation.tsv"))

#' ### Mock figure for fatty acids
#' 
#' 1. we find the terms
fAcids <- eCannot[grepl("fatty",eCannot$PDESC,ignore.case=TRUE),]

#' 2. we get the transcript lists and their length
fAcids %<>% group_by(PDESC) %>% summarise(TX=list(TxID)) %>% mutate(Ntx=elementNROWS(TX))

#' 3. we find how many are in the DEG lists
fAcids %<>% bind_cols(sapply(sapply(degList,"[[","ID"),function(d,p){sapply(lapply(p,"%in%",d),sum)},fAcids$TX))

#' 4. We calculate a frequency
freq <- fAcids[,names(degList)] / fAcids$Ntx

names(freq) <- paste0("percent_",names(freq))

fAcids %<>% bind_cols(freq)

#' 5. We plot the frequency (0 to 0.5 to 1 - blue to white to red) as a heatmap, removed the
#' dendrograms and clustering. Add the numbers as cell labels.
gplots::heatmap.2(rbind(c(percent_acute=0,percent_early=0.5,percent_late=1),
  as.matrix(fAcids[,c("percent_acute","percent_early","percent_late")])),
                  trace = "none",dendrogram="none",
                  Colv=FALSE,Rowv=FALSE,key=FALSE,
                  labCol = names(degList),margins=c(5,20),
                  labRow = c("colour key",fAcids$PDESC),
                  breaks = seq(0,1,0.01),
                  cexCol = 1,cexRow=1,notecol="black",
                  cellnote=rbind(c("0%","50%","100%"),
                                 as.matrix(fAcids[,c("acute","early","late")])),
                  col=hpal,rowsep=1,sepcolor="white")

#' 6. Write the results
dir.create(here("data/analysis/kegg"),showWarnings=FALSE)
write_tsv(fAcids %>% mutate(TX=sapply(fAcids$TX,paste,collapse="|")),
          here("data/analysis/kegg/fatty_acids_pathways_summary.tsv"))

#' ## Figure 5
#' 
#' 1. Getting GO related to cold
toicold <- annot[grepl("cold",annot$Term),c("GOID","Term")]
toicold %<>% separate_longer_delim(everything(),delim="|") %>% distinct() %>% filter(grepl("cold",Term))

message(sprintf("There are %s GO terms related to cold in the assembly",nrow(toicold)))
knitr::kable(toicold)

#' 2. Getting genes with GO related to cold
goicold <- unlist(as.vector(annot[Reduce("|",lapply(toicold$GOID,grepl,annot$GOID)),"TxID"]))
message(sprintf("There are %s candidates in the assembly",length(goicold)))

#' 3. Counting cold-related genes in DEG
countcold <- as.data.frame(matrix(data = 0, nrow = length(degList), ncol = 2, dimnames = list(names(degList),c("Up","Down"))))
countcold$Up <- as.numeric(lapply(names(degList), function(x) sum((degList[[x]][["ID"]] %in% goicold)*(degList[[x]][["log2FoldChange"]] > 0))))
countcold$Down <- as.numeric(lapply(names(degList), function(x) sum((degList[[x]][["ID"]] %in% goicold)*(degList[[x]][["log2FoldChange"]] < 0))))
countcold$NonSig <- length(goicold) - countcold$Up - countcold$Down
countcold$time <- rownames(countcold)

countcoldtimepoint <- as.data.frame(matrix(data = 0, nrow = length(degListPerTimePoint), ncol = 2, dimnames = list(names(degListPerTimePoint),c("Up","Down"))))
countcoldtimepoint$Up <- as.numeric(lapply(names(degListPerTimePoint), function(x) sum((degListPerTimePoint[[x]][["ID"]] %in% goicold)*(degListPerTimePoint[[x]][["log2FoldChange"]] > 0))))
countcoldtimepoint$Down <- as.numeric(lapply(names(degListPerTimePoint), function(x) sum((degListPerTimePoint[[x]][["ID"]] %in% goicold)*(degListPerTimePoint[[x]][["log2FoldChange"]] < 0))))
countcoldtimepoint$NonSig <- length(goicold) - countcoldtimepoint$Up - countcoldtimepoint$Down
countcoldtimepoint$time <- rownames(countcoldtimepoint)

#' 4. Plotting
dat <- melt(countcold)
dat$DE <- factor(dat$variable, levels = c("Up","NonSig","Down"))
ggplot(dat,aes(y=value,x=time,fill=DE)) +
	geom_col(position = "dodge") +
	xlab("Response") +
	ylab("Number of cold-related differentially expressed genes") +
	ggtitle("Cold related genes: per response")

dat <- melt(countcoldtimepoint)
dat$DE <- factor(dat$variable, levels = c("Up","NonSig","Down"))
dat$time <- factor(dat$time,levels = c("1hr_vs_ctrl","4hrs_vs_ctrl","12hrs_vs_ctrl","24hrs_vs_ctrl","72hrs_vs_ctrl","120hrs_vs_ctrl"))
ggplot(dat,aes(y=value,x=time,fill=DE)) +
	geom_col(position = "dodge") +
	xlab("Time") +
	ylab("Number of cold-related differentially expressed genes") +
	ggtitle("Cold related genes: per timepoint")

#' ## Figure 6 (and 7)
#' 
#' Plotting the Fatty acid biosynthesis pathway: ec00061 - we actually use the one for KEGG orthologs (ko)
pspecies <- "ko"
pnumber <- "00061"

#' 1. we get the pathway information
pwy <- keggGet(paste0(pspecies,pnumber))

#' 2. Get the enzymes
#'
#' From the pathway
#' 
ko <- pwy[[1]]$ORTHOLOGY
ko <- ko[grepl("EC",ko)]
tb <- tibble(KO=names(ko),EC=ko) %>%
  mutate(EC=sub("]","",sub(".*EC:","",EC))) %>% 
  separate_longer_delim(EC," ")

#' From the annotation
ez <- eCannot %>% mutate(EC=str_sub(EC,4)) %>% 
  filter(EC %in% (tb %>% select(EC) %>% distinct() %>% unlist(use.names=FALSE)) & PID==paste0("ec",pnumber))

#' 3. Collate the expression values
ez %<>% left_join(scaledAvgExp,by="TxID")

#' 4. Summarise
ez %<>% left_join(tb,by="EC",relationship = "many-to-many") %>% 
  group_by(KO) %>% summarise(across(where(is_double),median))
gdat <- ez %>% column_to_rownames("KO") %>% as.matrix()

#' 5. plot
pathview(
  gene.data=gdat,
  species = pspecies,
  pathway.id = pnumber,
  out.suffix = "fatty-acid-biosynthesis",
  kegg.dir=here("data/analysis/kegg"),
  low=c("blue","blue"))

#' 6. mv the result file in the result folder
outfile=list.files(here(),pattern=paste0(pspecies,pnumber,".*"))
file.rename(outfile,here("data/analysis/kegg",outfile))

#' ## Figure 8
#' 
#' 1. Getting GO related to cell wall
toicellwall <- annot[grepl("cell wall",annot$Term),c("GOID","Term")]
toicellwall %<>% separate_longer_delim(everything(),delim="|") %>% distinct() %>% filter(grepl("cell wall",Term))

message(sprintf("There are %s GO terms related to cold in the assembly",nrow(toicellwall)))
knitr::kable(toicellwall)

#' 2. Getting genes with GO related to cold
goicellwall <- unlist(as.vector(annot[Reduce("|",lapply(toicellwall$GOID,grepl,annot$GOID)),"TxID"]))
message(sprintf("There are %s candidates in the assembly",length(goicellwall)))

#' 3. Counting cold-related genes in DEG
countcellwall <- as.data.frame(matrix(data = 0, nrow = length(degList), ncol = 2, dimnames = list(names(degList),c("Up","Down"))))
countcellwall$Up <- as.numeric(lapply(names(degList), function(x) sum((degList[[x]][["ID"]] %in% goicellwall)*(degList[[x]][["log2FoldChange"]] > 0))))
countcellwall$Down <- as.numeric(lapply(names(degList), function(x) sum((degList[[x]][["ID"]] %in% goicellwall)*(degList[[x]][["log2FoldChange"]] < 0))))
countcellwall$NonSig <- length(goicellwall) - countcellwall$Up - countcellwall$Down
countcellwall$time <- rownames(countcellwall)

countcellwalltimepoint <- as.data.frame(matrix(data = 0, nrow = length(degListPerTimePoint), ncol = 2, dimnames = list(names(degListPerTimePoint),c("Up","Down"))))
countcellwalltimepoint$Up <- as.numeric(lapply(names(degListPerTimePoint), function(x) sum((degListPerTimePoint[[x]][["ID"]] %in% goicellwall)*(degListPerTimePoint[[x]][["log2FoldChange"]] > 0))))
countcellwalltimepoint$Down <- as.numeric(lapply(names(degListPerTimePoint), function(x) sum((degListPerTimePoint[[x]][["ID"]] %in% goicellwall)*(degListPerTimePoint[[x]][["log2FoldChange"]] < 0))))
countcellwalltimepoint$NonSig <- length(goicellwall) - countcellwalltimepoint$Up - countcellwalltimepoint$Down
countcellwalltimepoint$time <- rownames(countcellwalltimepoint)

countcellwall$sample <- factor(rownames(countcellwall),levels = c("60min_vs_std","4hrs_vs_std","12hrsvs_std","24hrs_vs_std","72hrs_vs_std","120hrs_vs_std"))
#' 4. Plotting the expression
dat <- melt(countcellwall)
dat$DE <- factor(dat$variable, levels = c("Up","NonSig","Down"))
ggplot(dat,aes(y=value,x=time,fill=DE)) +
	geom_col(position = "dodge") +
	xlab("Response") +
	ylab("Number of cell wall-related differentially expressed genes") +
	ggtitle("Cellwall related genes: per response")

dat <- melt(countcellwalltimepoint)
dat$DE <- factor(dat$variable, levels = c("Up","NonSig","Down"))
dat$time <- factor(dat$time,levels = c("1hr_vs_ctrl","4hrs_vs_ctrl","12hrs_vs_ctrl","24hrs_vs_ctrl","72hrs_vs_ctrl","120hrs_vs_ctrl"))
ggplot(dat,aes(y=value,x=time,fill=DE)) +
	geom_col(position = "dodge") +
	xlab("Time") +
	ylab("Number of cell wall-related differentially expressed genes") +
	ggtitle("Cellwall related genes: per timepoint")

#' 5. Reading cell wall thickness
cellwallthickness <- reshape(read.table(here("doc/Cell-wall-measurement.tsv"),header = TRUE), idvar = "time_in_hours", timevar = "treatment", direction = "wide")
cellwallthickness$minus <- cellwallthickness$width_in_um.cold - cellwallthickness$width_in_um.control
cellwallthickness$divide <- cellwallthickness$width_in_um.cold / cellwallthickness$width_in_um.control
cellwallthickness$time_in_hours <- factor(cellwallthickness$time_in_hours, levels = c(1,4,12,24,48,72,120,240))

#' 6. Plotting expression with cellwall thickness
dat <- countcellwalltimepoint %>% 
	select(Up, Down, time) %>%
	mutate(Down = (-1)*Down) %>%
	melt()
dat$time <- factor(dat$time,levels = c("1hr_vs_ctrl","4hrs_vs_ctrl","12hrs_vs_ctrl","24hrs_vs_ctrl","72hrs_vs_ctrl","120hrs_vs_ctrl"))
dat$DE <- factor(dat$variable, levels = c("Down","Up"))

dat2 <- cellwallthickness[(cellwallthickness$time_in_hours != 240)&(cellwallthickness$time_in_hours != 48),1:3]
dat2$time <- factor(c("1hr_vs_ctrl","4hrs_vs_ctrl","12hrs_vs_ctrl","24hrs_vs_ctrl","72hrs_vs_ctrl","120hrs_vs_ctrl"),levels = c("1hr_vs_ctrl","4hrs_vs_ctrl","12hrs_vs_ctrl","24hrs_vs_ctrl","72hrs_vs_ctrl","120hrs_vs_ctrl"))

ggplot() + geom_col(data=dat,mapping = aes(y=value,x=time,fill=DE)) +
	geom_line(data=dat2,group=1,mapping=aes(x=time,y=width_in_um.control*200,colour="control"),size = 1) +
	geom_line(data=dat2,group=1,mapping=aes(x=time,y=width_in_um.cold*200,colour="cold"),size = 1) +
	xlab("Time") +
	scale_y_continuous(name = "Number of cell wall-related differentially expressed genes",
										 sec.axis = sec_axis(trans = ~ . / 200,
										 										name = "Mean cell wall thickness (um)")) +
	scale_color_manual(name = "Cell wall thickness", values = c("control" = "brown", "cold" = "darkgreen"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
