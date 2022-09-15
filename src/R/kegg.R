#' ---
#' title: "Retrieve lignin pathway"
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
suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(KEGGREST)
  library(magrittr)
  library(matrixStats)
  library(pathview)
  library(readr)
  library(tibble)
  library(tidyr)
})

#' * Annotation
annotation <- read_tsv(here("data/b2g/blast2go_20190117_export.txt.gz"),show_col_types=FALSE) %>% 
  mutate(`Sequence Name`=sub("\\.p\\d+$","",`Sequence Name`))

#' We extract the enzyme code for every transcript (that has one)
enzyme_map <- annotation %>% select(`Sequence Name`,`Enzyme Code`) %>% 
    rename(ID=`Sequence Name`,EC=`Enzyme Code`) %>% 
    filter(!is.na(EC)) %>% separate_rows(EC,sep="\\|")

#' * Expression
expression <- read_tsv(here("data/analysis/salmon/variance-stabilised_model-aware_gene-expression_data.tsv"),
                       show_col_types=FALSE)

#' Reverse sample swap
expression %<>% rename(B2_2_24hrs_D="B2_2_12hrs_D",B2_2_12hrs_D="B2_2_24hrs_D")

#' # KEGG data
#' 
#' Retrieve the pathway information. 
#' 
#' pathways in KEGG have multiple prefixes. `map` describe the pathway, while other prefixes
#' defines organisms or set thereof.
#' 
#' `ko` stands for kegg orthologs while `cre` is for clamydomonas
#' 
#' _e.g._ cre00010 is clamydomonas' glucolysis pathway
#' 
#' Here we retrieve the kjegg ortholog lignin pathway (phenylpropanoid)
ptw <- keggGet("ko00940")

#' We extract a map of its enzymes and corresponding codes
ec_ko_map <- tibble(EC=gsub(".*\\[|\\]","",ptw[[1]]$ORTHOLOGY),
                    KO=names(ptw[[1]]$ORTHOLOGY)) %>% 
  mutate(EC=gsub(" "," EC:",EC)) %>% 
  separate_rows(EC,sep=" ")

#' We look them up in our annotation
lignin <- enzyme_map %>% filter(EC %in% ec_ko_map$EC)

#' Finally we extract our gene of interest
goi_lignin <- lignin %>% select(ID) %>% unlist()

#' # Further analysis examples
#' 
#' ## Comparing set of genes of interest
#' Common between lignin and cell wall
# intersect(goi_cellwall,goi_lignin)
#' Unique in Cell wall
# setdiff(goi_cellwall,goi_lignin)
#' Unique in Lignin
# setdiff(goi_lignin, goi_cellwall) => unique in lignin
#' Common to both
# union(goi_cellwall,goi_lignin) => everything

#' ## Pathview visualisation example
#' We pick the std samples. This is done using a string that is unique to the column name.
#' You can change that variable content to check other time points.
timepoint="std"

#' We associate the expression data with the enzyme code we identified in the annotation
#' that exists in the selected pathway. We filter for those that are expressed
ko_expression = left_join(
  left_join(lignin, ec_ko_map,by=c("EC")),
  expression %>% filter(rowname %in% goi_lignin) %>% 
    select("rowname",contains(timepoint)) %>% rowwise() %>% 
    summarise(ID=rowname,m=median(c_across(contains(timepoint)))) %>% 
    filter(m>0),
  by=c("ID")) %>% filter(!is.na(m)) %>% 
  select(KO,m) %>% group_by(KO) %>% 
  summarise(exp=sum(m))

#' We transform the expression data to the required input format
gene.data <- unlist(ko_expression$exp)
names(gene.data) <- unlist(ko_expression$KO)

#' We visualise the pathway. Here I use the `ko` species and 
#' as the `id` of the pathway, only its numerical part.
#' 
#' The figure will be saved under the pathway name and the extension `.pathview.png`
#' and saved in the `data/analysis/cell-wall` folder.
pathview(gene.data=gene.data,pathway.id="00940",
         species="ko",kegg.dir=here("data/analysis/cell-wall"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
