library(tidyverse)

inputf <-snakemake@input[["inputf"]]
chrom <- snakemake@params[["chrom"]]

all_gwas <- map_df(inputf,read_delim,delim="\t"
