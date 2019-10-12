library(here)
library(tidyverse)

transIDs <- read_delim(here("data/TransDecoder/Trinity.fasta.transdecoder.pep.ID.txt"),
                       delim=" ",col_names=c("pID","tdID","ORF",
                                             "type","length","score","pos")) %>% 
    mutate(pID=sub(">","",pID),
           type=sub("type:","",type)) %>% 
    mutate(tID=sub("\\.p[0-9]+","",pID))

allIDs <- gsub("^>| .*","",scan(here("data/trinity/TrinityIDs.txt"),
               sep="\n",what="character"))

message(sprintf("Out of %s transcripts, %s are putatively protein-coding",
                length(allIDs),
                nrow(transIDs)))

barplot(table(transIDs$type))

metadata <- tibble(tID=allIDs) %>% left_join(transIDs) %>% select(-ORF)

write_rds(metadata,path=here("data/analysis/metadata.rds"))
