library(readr)
library(dplyr)

read_delim(snakemake@input[["inf"]],
           delim="_",
           col_names=c("grp_id","trait_id","trait"))  %>%
    mutate(sample_id=1:n()) %>%
    write_delim(snakemake@output[["outf"]],delim="\t")

