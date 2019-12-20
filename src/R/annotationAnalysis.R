#' ---
#' title: "Annotation analysis"
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
  library(here)
  library(igraph)
  library(RSQLite)
  library(tidyverse)
})

#' * Graphics
mar <- par("mar")

#' * B2G Annotation
annot <- read_tsv(here("data/b2g/blast2go_20190117_export.txt.gz")) %>% 
  filter(`Sequence Description` != "---NA---") %>% 
  mutate(Tax=factor(gsub(".*Tax=| TaxID=.*","",`Sequence Description`)),
         TaxID=factor(gsub(".*TaxID=| RepID=.*","",`Sequence Description`)),
         MappingTaxID=factor(as.integer(names(
           sapply(lapply(lapply(`Mapping Taxa ID`,table),sort,decreasing=TRUE),"[",1)))))



#' # Taxonomy
tab <- sort(table(annot$Tax),decreasing=TRUE)
par(mar=c(12.1,4.1,0.1,0.1))
barplot(tab,las=2)
barplot(tab[1:10],las=2)
par(mar=mar)

# ouch
sum(annot$MappingTaxID==annot$TaxID,na.rm = TRUE)
sum(grepl("^\\d+$",levels(annot$TaxID)))
length(grepl("\\d+",levels(annot$TaxID)))
# needs fixing too!
# tail(levels(annot$TaxID))
# [1] "UniRef90_Q9P0J0NADH dehydrogenase"                                      
# [2] "UniRef90_R7YVC9Saccharopine dehydrogenase"                              
# [3] "UniRef90_UPI00029B39775,10-methylenetetrahydrofolate reductase n=1 Tax="
# [4] "UniRef90_UPI0006EAFD0F glycogen"                                        
# [5] "UniRef90_UPI000D1B4B2DFAD-dependent oxidoreductase n=1 Tax="            
# [6] "UniRef90_UPI000F74087Balcohol dehydrogenase"

# load the database and check the mapping IDs
con <- dbConnect(dbDriver("SQLite"),
                 dbname=here("taxonomy/20190820/taxonomy.sqlite"))

ids <- unique(c(levels(annot$MappingTaxID),
                levels(annot$TaxID)[grepl("^\\d+$",levels(annot$TaxID))]))
taxid_divnames <- dbGetQuery(con, paste0("SELECT d.div_nam, n.tax_id FROM division AS d ",
                                         "LEFT JOIN node AS n ON d.div_id = n.div_id ",
                                         "WHERE n.tax_id IN (",
                                         paste(ids, collapse = ","),
                                         ")"))

# Diversity: we get a lot of Viridiplantae
barplot(table(taxid_divnames$div_nam),las=2)

# Abundance
tb <- table(annot$MappingTaxID)
barplot(by(as.vector(tb[match(names(tb),taxid_divnames$tax_id)]),
           taxid_divnames$div_nam,sum))

tb <- table(annot$TaxID)
tb <- tb[grepl("^\\d+$",names(tb))]
m.sel <- match(names(tb),taxid_divnames$tax_id)
barplot(by(as.vector(tb[m.sel]),
           taxid_divnames$div_nam[m.sel],sum),las=2)

# This is also fishy
#head(levels(annot$Tax))
# [1] ""                                   "50 kb inversion clade"             
# [3] "Abrus precatorius"                  "Absidia glauca"              
# Get the ones by name
nam.div <- dbGetQuery(con,paste("SELECT d.div_nam, t.nam from division d",
                                "LEFT JOIN node n ON d.div_id == n.div_id",
                                "LEFT JOIN taxonomy t ON t.tid == n.tax_id",
                                "WHERE t.nam IN (",paste("'",
                                                         unique(gsub(" \\(.*| var\\..*","",levels(annot$Tax))),
                                                         "'",sep="",collapse=","),");"))

barplot(table(nam.div$div_nam),las=2)

tb <- table(gsub(" \\(.*| var\\..*","",annot$Tax))
m.sel <- match(names(tb),nam.div$nam)
barplot(by(as.vector(tb[m.sel]),
           nam.div$div_nam[m.sel],sum),las=2)

# find the Embryophyta / ChloroPhyta
edges <- dbGetQuery(con,"SELECT tax_id, parent_id from node;")

## get tax id
taxids <- dbGetQuery(con,"SELECT tid, nam from taxonomy where nam in ('Chlorophyta','Embryophyta','Streptophyta','Viridiplantae')")

## create the tax graph, break it at Chlorophyta and get the membership
graph <- graph.edgelist(as.matrix(edges[!edges$parent_id %in% 
                                          edges[edges$tax_id == 
                                                  taxids[taxids$nam=="Embryophyta","tid"],
                                                "parent_id"],]))
mem <- clusters(graph)$membership

## get the taxon name part of the Chlorophyta division
ephylum <- dbGetQuery(con,paste("SELECT DISTINCT tid,nam from taxonomy WHERE tid in (",paste(which(mem==mem[taxids[taxids$nam=="Embryophyta","tid"]]),collapse=","),")"))


annot$Chlorophyta <- annot$UniRef90.taxon %in% cphylum
table(annot$Chlorophyta)

dbDisconnect(con)

# cross validate

# look at the GC vs expression of the unannotated sequences


# select only the algae
IDs <- taxa$ID[taxa$Tax == "Tetradesmus obliquus"]

# do we have unique IDs
sum(duplicated(sub("\\.p[0-9]+$","",IDs)))
IDs <- sub("\\.p[0-9]+$","",IDs)

dir.create(here("data/analysis/annotation"),recursive = TRUE)
write(as.character(IDs),here("data/analysis/annotation/algae-IDs.txt"))
