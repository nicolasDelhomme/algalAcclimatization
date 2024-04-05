#' ---
#' title: "Fig 5 anf Fig 8"
#' author: "Fai"
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
  library(tidyverse)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(magrittr)
  library(matrixStats)
	library(ggVennDiagram)
	library(UpSetR)
})

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Annotation
annotation <- read_tsv(here("data/b2g/blast2go_20190117_export.txt.gz"),show_col_types=FALSE) %>% 
  mutate(`Sequence Name`=sub("\\.p\\d+$","",`Sequence Name`))

#' * GO dictionary
go_annot <- annotation %>% select(`Annotation GO ID`,`Annotation GO Term`,`Annotation GO Category`) %>% 
  rename(ID=`Annotation GO ID`,Term=`Annotation GO Term`,Category=`Annotation GO Category`) %>% 
  filter(!is.na(ID)) 

go_annot <- tibble(ID=unlist(str_split(go_annot$ID,"\\|"),use.names=FALSE),
                   Term=unlist(str_split(go_annot$Term,"\\|"),use.names=FALSE),
                   Category=unlist(str_split(go_annot$Category,"\\|"),use.names=FALSE)) %>% distinct()

#' # Gene of interest Cell wall
#' ## Lookup
#' We use a combination of the Sequence Description and Gene Ontology (GO) term to fetch candidates
goicellwall <- unlist(annotation[grepl("cell wall",annotation$`Sequence Description`) | 
  grepl("cell wall",annotation$`Annotation GO Term`),"Sequence Name"],use.names=FALSE)

#' replace (\\.p\\d+$) .p followed by any number of digits. The $ the sign that indicates the end of the text
goicellwall <- sub("\\.p\\d+$","",goi)

#' sanity
stopifnot(length(goicellwall)==length(unique(goicellwall)))
#write.table(goicellwall,here("doc/GOI_cellwall.txt"),row.names = FALSE, col.names = FALSE, quote = FALSE)

#' # Gene of interest Cold
#' ## Lookup
toi <- annotation[grepl("cold",annotation$`Annotation GO Term`),c("Annotation GO ID","Annotation GO Term")]

toi %<>% separate_longer_delim(everything(),delim="|") %>% distinct() %>% filter(grepl("cold",`Annotation GO Term`))

#' ## Gene of interest
goi <- unlist(annotation[Reduce("|",lapply(toi$`Annotation GO ID`,grepl,annotation$`Annotation GO ID`)),"Sequence Name"])

#' sanity
stopifnot(length(goi)==length(unique(goi)))


#' * Get the DE result
load(here("data/analysis/salmon/dds-sample-swap-corrected.rda"))

#' ## Select
#' * only the complete proteins
metadata <- read_rds("/mnt/picea/projects/algae/cfunk/algal-acclimatization/analysis/metadata.rds")
tx.sel <- metadata$type=="complete" & ! is.na(metadata$type)
sum(tx.sel) #34650
nrow(dds) #49477
nrow(metadata) #285908
tx.sel2 <- rownames(dds) %in% metadata$tID[tx.sel]
dds <- dds[tx.sel2,]

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE,fitType='local')
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Differential Expression
dds <- DESeq(dds,fitType="local")

#' * Results
#' This shows all possible contrasts
resultsNames(dds)

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
															padj=0.01,lfc=0.5,
															plot=TRUE,verbose=TRUE,
															export=TRUE,default_dir=here("data/analysis/DE"),
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

#This is what available from previous DE
#h1.vs.std
#h4.vs.h1
#h12.vs.h4
#h24.vs.h12
#h72.vs.h24
#h120.vs.h72

#Olivia asked for everything vs std
h4.vs.std <- extract_results(dds,vst,"Time_4hrs_vs_std",
														 default_prefix="4hrs_vs_std_",
														 labels=dds$Time,
														 sample_sel=dds$Time %in% c("std","4hrs"),
														 export=FALSE,plot=FALSE)

h12.vs.std <- extract_results(dds,vst,"Time_12hrs_vs_std",
														 default_prefix="12hrsvs_std_",
														 labels=dds$Time,
														 sample_sel=dds$Time %in% c("std","12hrs"),
														 export=FALSE,plot=FALSE)

h24.vs.std <- extract_results(dds,vst,"Time_24hrs_vs_std",
														 default_prefix="24hrs_vs_std_",
														 labels=dds$Time,
														 sample_sel=dds$Time %in% c("std","24hrs"),
														 export=FALSE,plot=FALSE)

h72.vs.std <- extract_results(dds,vst,"Time_72hrs_vs_std",
														 default_prefix="72hrs_vs_std_",
														 labels=dds$Time,
														 sample_sel=dds$Time %in% c("std","72hrs"),
														 export=FALSE,plot=FALSE)

h120.vs.std <- extract_results(dds,vst,"Time_120hrs_vs_std",
															default_prefix="120hrs_vs_std_",
															labels=dds$Time,
															sample_sel=dds$Time %in% c("std","120hrs"),
															export=FALSE,plot=FALSE)

h1.vs.std <- extract_results(dds,vst,"Time_60min_vs_std",
														 default_prefix="60min_vs_std_",
														 labels=dds$Time,
														 sample_sel=dds$Time %in% c("std","60min"),
														 export=FALSE,plot=FALSE)

#DE files
files <- list.files(here("data/analysis/DE"), pattern = "_std_genes.csv$", full.names = TRUE)
degs <- lapply(files,fread)
names(degs) <- sub("_genes","",sub(".csv","",basename(files)))

countcellwall <- as.data.frame(matrix(data = 0, nrow = length(degs), ncol = 2))
colnames(countcellwall) <- c("Up","Down")
rownames(countcellwall) <- c("60min_vs_std","4hrs_vs_std","12hrsvs_std","24hrs_vs_std","72hrs_vs_std","120hrs_vs_std")
for(i in 1:length(degs)){
	countcellwall[rownames(countcellwall) == names(degs)[[i]],1] <- sum((degs[[i]]$V1 %in% goi) * (degs[[i]]$log2FoldChange > 0))
	countcellwall[rownames(countcellwall) == names(degs)[[i]],2] <- sum((degs[[i]]$V1 %in% goi) * (degs[[i]]$log2FoldChange < 0))
}
countcellwall$sample <- factor(rownames(countcellwall),levels = c("60min_vs_std","4hrs_vs_std","12hrsvs_std","24hrs_vs_std","72hrs_vs_std","120hrs_vs_std"))
countcellwall$Down <- (-1)*countcellwall$Down
dat <- melt(countcellwall)
dat$DEGs <- factor(dat$variable,levels=c("Down","Up"))
ggplot(dat,aes(y=value,x=sample,fill=variable)) +
	geom_col()

#measured thickness
cellwallthickness <- reshape(read.table(here("doc/Cell-wall-measurement.tsv"),header = TRUE), idvar = "time_in_hours", timevar = "treatment", direction = "wide")
cellwallthickness$minus <- cellwallthickness$width_in_um.cold - cellwallthickness$width_in_um.control
cellwallthickness$divide <- cellwallthickness$width_in_um.cold / cellwallthickness$width_in_um.control
cellwallthickness$time_in_hours <- factor(cellwallthickness$time_in_hours, levels = c(1,4,12,24,48,72,120,240))
ggplot(cellwallthickness,aes(y=minus,x=time_in_hours,group=1)) +
	geom_line()
ggplot(cellwallthickness,aes(y=divide,x=time_in_hours,group=1)) +
	geom_line()

dat2 <- cellwallthickness[(cellwallthickness$time_in_hours != 240)&(cellwallthickness$time_in_hours != 48),c(1,4)]
dat2$sample <- factor(c("60min_vs_std","4hrs_vs_std","12hrsvs_std","24hrs_vs_std","72hrs_vs_std","120hrs_vs_std"), levels = c("60min_vs_std","4hrs_vs_std","12hrsvs_std","24hrs_vs_std","72hrs_vs_std","120hrs_vs_std"))
ggplot(dat2,aes(y=minus,x=sample,group=1)) +
	geom_line()

ggplot() + geom_col(data=dat,mapping = aes(y=value,x=sample,fill=DEGs)) +
	geom_line(data=dat2,group=1,mapping=aes(x=sample,y=minus*1000)) +
	scale_y_continuous(name = "Number of cell wall-related differentially expressed genes",
										 sec.axis = sec_axis(trans = ~ . / 1000,
										 										name = "Differences between cell wall thickness of cold vs control (um)"))

########################################Cold##############################################

toi <- annotation[grepl("cold",annotation$`Annotation GO Term`),c("Annotation GO ID","Annotation GO Term")]

toi %<>% separate_longer_delim(everything(),delim="|") %>% distinct() %>% filter(grepl("cold",`Annotation GO Term`))

#' ## Gene of interest
goi <- unlist(annotation[Reduce("|",lapply(toi$`Annotation GO ID`,grepl,annotation$`Annotation GO ID`)),"Sequence Name"])

#' sanity
stopifnot(length(goi)==length(unique(goi)))

countcold <- as.data.frame(matrix(data = 0, nrow = length(degs), ncol = 2))
colnames(countcold) <- c("Up","Down")
rownames(countcold) <- c("60min_vs_std","4hrs_vs_std","12hrsvs_std","24hrs_vs_std","72hrs_vs_std","120hrs_vs_std")
for(i in 1:length(degs)){
	countcold[rownames(countcold) == names(degs)[[i]],1] <- sum((degs[[i]]$V1 %in% goi) * (degs[[i]]$log2FoldChange > 0))
	countcold[rownames(countcold) == names(degs)[[i]],2] <- sum((degs[[i]]$V1 %in% goi) * (degs[[i]]$log2FoldChange < 0))
}
countcold$sample <- factor(rownames(countcellwall),levels = c("60min_vs_std","4hrs_vs_std","12hrsvs_std","24hrs_vs_std","72hrs_vs_std","120hrs_vs_std"))
#countcold$Down <- (-1)*countcold$Down
countcold$Nonsig <- length(goi) - (countcold$Up + countcold$Down)
dat <- melt(countcold)
ggplot(dat,aes(y=value,x=sample,fill=variable)) +
	geom_col(position = "dodge")

###########################Getting the table##########################
goicellwall <- read.table(here("doc/GOI_cellwall.txt"))
colnames(goicellwall) <- "ID"
goicellwall$group <- "CellWall" 
goicoldres <- read.table(here("doc/GOI_coldres.txt"))
colnames(goicoldres) <- "ID"
goicoldres$group <- "ColdRes" 
goi <- bind_rows(goicellwall,goicoldres)

files <- list.files(here("data/analysis/DE"), pattern = "_std_genes.csv$", full.names = TRUE)
degs <- lapply(files,fread)
names(degs) <- paste0("de",sub("_genes","",sub(".csv","",basename(files))))

goi$de60min_vs_std <- "Nonsig"
goi$de4hrs_vs_std <- "Nonsig"
goi$de12hrsvs_std <- "Nonsig"
goi$de24hrs_vs_std <- "Nonsig"
goi$de72hrs_vs_std <- "Nonsig"
goi$de120hrs_vs_std <- "Nonsig"

load(here("data/analysis/salmon/dds-sample-swap-corrected.rda"))

for(i in 1:length(goi$ID)){
	for(j in 1:length(degs)){
		if(goi$ID[i] %in% degs[[j]]$V1){
			if(degs[[j]]$log2FoldChange[degs[[j]]$V1 == goi$ID[i]] > 0) {
				goi[i,names(degs)[j]] <- "Up"
			} else if(degs[[j]]$log2FoldChange[degs[[j]]$V1 == goi$ID[i]] < 0){
				goi[i,names(degs)[j]] <- "Down"
			}
		} else if(!goi$ID[i] %in% rownames(dds)){
			goi[i,names(degs)[j]] <- "NA"
		}
	}
}

goi[!goi$ID %in% rownames(dds),]
goi[goi$ID %in% rownames(dds),]
nrow(dds) #49477
sum(goi$ID %in% rownames(dds)) #212
length(goi$ID) #326
length(unique(goi$ID)) #324
goi[duplicated(goi$ID),]
goi[grep("TRINITY_DN48375",goi$ID),] #both cell wall and cold res but non sig in all

write.csv(goi,here("doc/GOI_cellwall_and_coldres.csv"), row.names = FALSE, quote = FALSE)

##############################The lipid and carb##############################
#From Olivia
selectedCarb <- c(
	"carbohydrate metabolic process",
	"carbohydrate biosynthetic process",
	"carbohydrate derivative biosynthetic process",
	"response to carbohydrate",
	"carbohydrate transport",
	"carbohydrate transmembrane transport",
	"cellular carbohydrate metabolic process",
	"carbohydrate utilization",
	"cellular carbohydrate biosynthetic process")

selectedLipid <- c(
	"lipid metabolic process",
	"glycolipid biosynthetic process",
	"phospholipid transport",
	"lipid transport",
	"lipid storage",
	"phospholipid binding",
	"lipid transporter activity",
	"positive regulation of lipid storage",
	"glycerophospholipid biosynthetic process",
	"galactolipid biosynthetic process",
	"phospholipid biosynthetic process",
	"membrane lipid biosynthetic process",
	"lipid droplet",
	"lipid biosynthetic process",
	"lipid localization",
	"phospholipid homeostasis",
	"cellular lipid catabolic process",
	"lipid homeostasis",
	"phospholipid catabolic process",
	"cellular lipid metabolic process",
	"regulation of lipid biosynthetic process",
	"regulation of lipid catabolic process",
	"regulation of membrane lipid distribution",
	"regulation of lipid storage",
	"lipid transport involved in lipid storage",
	"positive regulation of lipid biosynthetic process",
	"intermembrane lipid transfer",
	"lipid transfer activity",
	"negative regulation of lipid catabolic process",
	"positive regulation of lipid catabolic process")

goicarb <- unlist(annotation[grepl(paste(selectedCarb,collapse = "|"),annotation$`Sequence Description`) | 
														grepl(paste(selectedCarb,collapse = "|"),annotation$`Annotation GO Term`),"Sequence Name"],use.names=FALSE)
goilipid <- unlist(annotation[grepl(paste(selectedLipid,collapse = "|"),annotation$`Sequence Description`) | 
														 	grepl(paste(selectedLipid,collapse = "|"),annotation$`Annotation GO Term`),"Sequence Name"],use.names=FALSE)

goicellwall <- read.table(here("doc/GOI_cellwall.txt"))
colnames(goicellwall) <- "ID"
goicellwall$group <- "CellWall" 
goicoldres <- read.table(here("doc/GOI_coldres.txt"))
colnames(goicoldres) <- "ID"
goicoldres$group <- "ColdRes" 
goi <- bind_rows(goicellwall,goicoldres)

files <- list.files(here("data/analysis/DE"), pattern = "_std_genes.csv$", full.names = TRUE)
degs <- lapply(files,fread)
names(degs) <- paste0("de",sub("_genes","",sub(".csv","",basename(files))))

goi$de60min_vs_std <- "Nonsig"
goi$de4hrs_vs_std <- "Nonsig"
goi$de12hrsvs_std <- "Nonsig"
goi$de24hrs_vs_std <- "Nonsig"
goi$de72hrs_vs_std <- "Nonsig"
goi$de120hrs_vs_std <- "Nonsig"

load(here("data/analysis/salmon/dds-sample-swap-corrected.rda"))

for(i in 1:length(goi$ID)){
	for(j in 1:length(degs)){
		if(goi$ID[i] %in% degs[[j]]$V1){
			if(degs[[j]]$log2FoldChange[degs[[j]]$V1 == goi$ID[i]] > 0) {
				goi[i,names(degs)[j]] <- "Up"
			} else if(degs[[j]]$log2FoldChange[degs[[j]]$V1 == goi$ID[i]] < 0){
				goi[i,names(degs)[j]] <- "Down"
			}
		} else if(!goi$ID[i] %in% rownames(dds)){
			goi[i,names(degs)[j]] <- "NA"
		}
	}
}

goi[!goi$ID %in% rownames(dds),]
goi[goi$ID %in% rownames(dds),]
nrow(dds) #49477
sum(goi$ID %in% rownames(dds)) #212
length(goi$ID) #326
length(unique(goi$ID)) #324
goi[duplicated(goi$ID),]
goi[grep("TRINITY_DN48375",goi$ID),] #both cell wall and cold res but non sig in all

write.csv(goi,here("doc/GOI_cellwall_and_coldres.csv"), row.names = FALSE, quote = FALSE)


############### Usinf updated GO table
annot <- read_tsv(here("data/analysis/cold-response/algae_GO_annotation.tsv"),
									show_col_types=FALSE)

allannot <- read_tsv(here("data/b2g/blast2go_20190117_export.txt.gz"),show_col_types=FALSE) %>% 
	mutate(`Sequence Name`=sub("\\.p\\d+$","",`Sequence Name`))

#' ### Getting GO related to cold

toicold <- annot[grepl("cold",annot$Term),c("GOID","Term")]
toicold %<>% separate_longer_delim(everything(),delim="|") %>% distinct() %>% filter(grepl("cold",Term))

message(sprintf("There are %s GO terms related to cold in the assembly",length(toicold$GOID)))
knitr::kable(toicold)

#' ### Getting genes with GO related to cold

goicold <- unlist(annot[Reduce("|",lapply(toicold$GOID,grepl,annot$GOID)),"TxID"])
message(sprintf("There are %s candidates in the assembly",length(goicold)))
#85

#Try from allannot
toicold <- allannot[grepl("cold",allannot$`Annotation GO Term`),c("Annotation GO ID","Annotation GO Term")]
toicold %<>% separate_longer_delim(everything(),delim="|") %>% distinct() %>% filter(grepl("cold",`Annotation GO Term`))

message(sprintf("There are %s GO terms related to cold in the assembly",length(toicold$`Annotation GO ID`)))
knitr::kable(toicold)

#' ### Getting genes with GO related to cold

goicold2 <- unlist(allannot[Reduce("|",lapply(toicold$`Annotation GO ID`,grepl,allannot$`Annotation GO ID`)),"Sequence Name"])
message(sprintf("There are %s candidates in the assembly",length(goicold2)))
#123

#count if they are in dds
sum(goicold %in% rownames(dds))
sum(goicold2 %in% rownames(dds))
#all 85

#So, from allannot, not all of them are in dds
#from annot, all of them are in dds

#' 1. Getting GO related to cell wall
toicellwall <- annot[grepl("cell wall",annot$Term),c("GOID","Term")]
toicellwall %<>% separate_longer_delim(everything(),delim="|") %>% distinct() %>% filter(grepl("cell wall",Term))

message(sprintf("There are %s GO terms related to cold in the assembly",nrow(toicellwall)))
knitr::kable(toicellwall)
#8

#' 2. Getting genes with GO related to cold
goicellwall <- unlist(as.vector(annot[Reduce("|",lapply(toicellwall$GOID,grepl,annot$GOID)),"TxID"]))
message(sprintf("There are %s candidates in the assembly",length(goicellwall)))
#123

#Try from allannot
goicellwall2 <- unlist(allannot[grepl("cell wall",allannot$`Sequence Description`) | 
																 	grepl("cell wall",allannot$`Annotation GO Term`),"Sequence Name"],use.names=FALSE)

#' replace (\\.p\\d+$) .p followed by any number of digits. The $ the sign that indicates the end of the text
goicellwall2 <- sub("\\.p\\d+$","",goicellwall2)
#203

sum(goicellwall %in% rownames(dds))
sum(goicellwall2 %in% rownames(dds))
#all 127

#additional request from Olivia
#venn
res.list <- list("1 hr"=h1.vs.std,
								 "4 hr"=h4.vs.std,
								 "12 hr"=h12.vs.std,
								 "24 hr"=h24.vs.std,
								 "72 hr"=h72.vs.std,
								 "120 hr"=h120.vs.std)

ggVennDiagram(res.list, label_alpha = 0, label = "count") + 
	scale_fill_gradient(low="white",high = "dodgerblue") + 
	scale_x_continuous(expand = expansion(mult = .1))

upset(fromList(res.list),order.by="freq")

#
#So, from allannot, not all of them are in dds
#from annot, all of them are in dds

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```



