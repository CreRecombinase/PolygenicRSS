library(readr)
library(dplyr)
library(tidyr)
Nnum <- as.integer(snakemake@params[["sample_n"]])
stopifnot(length(Nnum)==1)


read_delim(snakemake@input[["inf"]],
           delim="\t")  %>%
    select(-sample_id) %>%
    sample_n(size=as.integer(snakemake@params[["sample_n"]]),replace=F)%>%
    unite("vcfid",everything(),sep="_") %>% 
write_delim(snakemake@output[["outf"]],delim="\t",col_names=F)

