library(here)
library(RSQLite)

test <- read.delim(here("data/b2g/blast2go_20190117_export.txt.gz"),as.is = TRUE)

sel <- test$Mapping.Taxa.ID!=""

mar <- par("mar")

sum(sel)
length(sel)

sel2 <- test$Sequence.Description != "---NA---"
sum(sel2)
sum(sel | sel2)

taxa<-data.frame(
  ID=test[sel2,"Sequence.Name"],
  MappingTaxID=as.integer(names(
    sapply(lapply(lapply(test[sel2,"Mapping.Taxa.ID"],table),
                  sort,decreasing=TRUE),"[",1))),
  Tax=gsub(".*Tax=| TaxID=.*","",test[sel2,"Sequence.Description"]),
  TaxID=gsub(".*TaxID=| RepID=.*","",test[sel2,"Sequence.Description"]))

barplot(sort(table(taxa$Tax),decreasing=TRUE),las=2)

tab <- sort(table(taxa$Tax),decreasing=TRUE)

# ouch
sum(taxa$MappingTaxID==taxa$TaxID,na.rm = TRUE)

par(mar=c(12.1,4.1,0.1,0.1))
barplot(tab[1:10],las=2)
par(mar=mar)

names(tab)[1:20]

# load the database and check the mapping IDs
con <- dbConnect(dbDriver("SQLite"),
                 dbname="/mnt/picea/storage/reference/Taxonomy/20150407/taxonomy.sqlite")
nam.div <- dbGetQuery(con,paste("SELECT d.div_nam, t.nam from division d",
                                "LEFT JOIN node n ON d.div_id == n.div_id",
                                "LEFT JOIN taxonomy t ON t.tid == n.tax_id",
                                "WHERE t.nam IN (",paste("'",
                                                         unique(gsub(" \\(.*| var\\..*","",names(tab)))
                                                         ,"'",sep="",collapse=","),");"))
dbDisconnect(con)

# cross validate
taxid_divnames <- dbGetQuery(con, paste0("SELECT d.div_nam, n.tax_id FROM division AS d ",
                                         "LEFT JOIN node AS n ON d.div_id = n.div_id ",
                                         "WHERE n.tax_id IN (",
                                         paste(unique(taxID), collapse = ","),
                                         ")"))

# look at the GC vs expression of the unannotated sequences


# select only the algae
IDs <- taxa$ID[taxa$Tax == "Tetradesmus obliquus"]

# do we have unique IDs
sum(duplicated(sub("\\.p[0-9]+$","",IDs)))
IDs <- sub("\\.p[0-9]+$","",IDs)

dir.create(here("data/analysis/annotation"),recursive = TRUE)
write(as.character(IDs),here("data/analysis/annotation/algae-IDs.txt"))
