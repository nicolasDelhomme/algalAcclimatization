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
  library(S4Vectors)
  library(here)
  library(igraph)
  library(magrittr)
  library(pander)
  library(RSQLite)
  library(tidyverse)
})

#' * Graphics
mar <- par("mar")

#' # Process
#' ## B2G Annotation
annot <- read_tsv(here("data/b2g/blast2go_20190117_export.txt.gz"),
                  col_types = cols(
                    .default = col_character(),
                    `Sequence Length` = col_double(),
                    `Annotation GO Count` = col_double(),
                    `Blast Hits Count` = col_double(),
                    `Blast Min E-Value` = col_double(),
                    `Blast Top Hit E-Value` = col_double(),
                    `Blast Top Hit Length` = col_double(),
                    `Blast Top Hit Alignment Length` = col_double(),
                    `Blast Top Hit Positives` = col_double(),
                    `Blast Top Hit HSPs Number` = col_double(),
                    `Blast Top Hit Frame` = col_double()
                  )) %>% 
  filter(`Sequence Description` != "---NA---") %>% 
  mutate(Tax=ifelse(grepl("Tax=",`Sequence Description`),
                    gsub(".*Tax=| TaxID=.*","",`Sequence Description`),NA),
         TaxID=factor(ifelse(grepl("TaxID=",`Sequence Description`),
                             gsub(".*TaxID=| RepID=.*","",`Sequence Description`),NA)))

#' Update some of the Taxonomy names
df=data.frame(
  from=c("Lyngbya sp.","Cyanothece sp.","Scenedesmus fuscus","Pediculus humanus subsp. corporis"),
  to=c("Lyngbya","Cyanothece","Scenedesmus","Pediculus humanus"),stringsAsFactors = FALSE)

uq <- unique(annot$Tax)
mp <- lapply(df$from,grep,uq)
uq[unlist(mp)] <- rep(df$to,elementNROWS(mp))

annot %<>% mutate(Tax=factor(Tax,levels=unique(uq)))

#' Validation
stopifnot(all(grepl("^\\d+$",levels(annot$TaxID))))

#' ## Taxonomy
#' Establish the connection
con <- dbConnect(dbDriver("SQLite"),
                 dbname=here("taxonomy/20190820/taxonomy.sqlite"))

#' ### Taxonomy ID - division mapping
mids <- strsplit(annot$`Mapping Taxa ID`,"\\|")
ids <- unique(c(na.omit(unique(unlist(mids))),levels(annot$TaxID)))
taxid_divnames <- dbGetQuery(con, paste0("SELECT d.div_nam, n.tax_id FROM division AS d ",
                                         "LEFT JOIN node AS n ON d.div_id = n.div_id ",
                                         "WHERE n.tax_id IN (",
                                         paste(ids, collapse = ","),
                                         ")"))

#' ### Taxonomy name - division mapping
nam.div <- dbGetQuery(con,paste("SELECT d.div_nam, t.nam from division d",
                                "LEFT JOIN node n ON d.div_id == n.div_id",
                                "LEFT JOIN taxonomy t ON t.tid == n.tax_id",
                                "WHERE t.nam IN (",paste("'",
                                                         unique(gsub(" \\(.*| var\\..*","",levels(annot$Tax))),
                                                         "'",sep="",collapse=","),");"))


#' ### Embryophyta name and ID
edges <- dbGetQuery(con,"SELECT tax_id, parent_id from node;")

#' Get the Embryophyta taxonomy ID
ephyta <- dbGetQuery(con,"SELECT tid from taxonomy where nam == 'Embryophyta'")$tid

#' Create the taxonomy graph, break it at Embryophyta and extract 
#' the members of the Embryophyta clade
graph <- graph.edgelist(as.matrix(edges[!edges$parent_id %in% 
                                          edges[edges$tax_id == ephyta,"parent_id"],]))
mem <- clusters(graph)$membership

#' And get their taxon name and ID
ephylum <- dbGetQuery(con,paste("SELECT DISTINCT tid,nam from taxonomy WHERE tid in (",
                                paste(which(mem==mem[ephyta]),collapse=","),")"))

#' Done, disconnect
dbDisconnect(con)

#' # Report
#' ## Individual Taxon by name
tab <- sort(table(annot$Tax),decreasing=TRUE)
par(mar=c(12.1,4.1,0.1,0.1))
barplot(tab,las=2)
barplot(tab[1:10],las=2)

#' Add a column to the annot
annot %<>% mutate(Tax.Div=nam.div$div_nam[match(levels(annot$Tax),
                                                nam.div$nam)][as.integer(annot$Tax)])

barplot(table(annot$Tax.Div),las=2)

#' ## Individual Taxon by ID
annot %<>% mutate(TaxID.Div=taxid_divnames$div_nam[match(levels(annot$TaxID),
                                                         taxid_divnames$tax_id)][as.integer(annot$TaxID)])

barplot(table(annot$TaxID.Div),las=2)

#' We are in good agreement
table(annot$Tax.Div == annot$TaxID.Div,useNA="ifany")

#' ## Mapping Taxa by ID
res <- lapply(split(taxid_divnames$div_nam[match(unlist(mids),taxid_divnames$tax_id)],
                    rep(1:length(mids),elementNROWS(mids))),table)

#' Most of the Mapping Taxon agree
pander(table(elementNROWS(res)))
barplot(table(elementNROWS(res)))

#' Merge the names, merge the abundances
annot %<>% mutate(MappingTaxIdDiv=factor(sapply(lapply(res,names),paste,collapse="|")),
                  MappingTaxIdDivAbundance=sapply(lapply(res,as.integer),paste,collapse="|"))

pander(table(annot$MappingTaxIdDiv))
barplot(table(annot$MappingTaxIdDiv),las=2,cex.names=.8)

annot %<>% mutate(MappingTaxIdViridiplantae=as.integer(MappingTaxIdDiv) %in% 
                    which(grepl("Viridiplantae",levels(MappingTaxIdDiv))))

table(annot$MappingTaxIdViridiplantae)

#' Again we are in good agreement
sum(annot$MappingTaxIdViridiplantae & annot$Tax.Div == "Viridiplantae",na.rm = TRUE)
par(mar=mar)

#' # Select
#' We base the selection on the IDs retrieved from the TaxID 
#' (the uniref diamond annotation)
annot %<>% mutate(TaxIdEmbryophyta=annot$TaxID %in% ephylum$tid)
message(sprintf("There are %s sequences annotated as land plants (Embryophyta))",sum(annot$TaxIdEmbryophyta)))

#' # Export
IDs <- annot$`Sequence Name`[annot$Tax.Div == "Viridiplantae" & ! annot$TaxIdEmbryophyta & ! is.na(annot$Tax.Div)]
IDs <- sub("\\.p[0-9]+","",IDs)

dir.create(here("data/analysis/annotation"),recursive = TRUE,showWarnings = FALSE)

write(as.character(IDs),here("data/analysis/annotation/algae-IDs.txt"))

#' # Conclusion
#' We have 49,477 sequence that are identified as being of algal origin
#' There are many more sequences that we could try to identify based on 
#' a GC / expression grouping (possibly) - or run a ML approach.

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
