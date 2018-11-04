library(tidyverse)

inputf <- snakemake@input[["inputf"]]
chrom <- snakemake@params[["chrom"]]
out_pref <-snakemake@params[["out_pref"]]
chunksize <- as.integer(snakemake@params[["chunk_size"]])
names(chunksize)  <- chrom
all_gwas <- map_df(inputf, read_delim, delim = "\t")

all_gwas %>% group_by( chrom ) %>%
    mutate(region_id = formatC(
               x = as.integer(gl(
                   n = ceiling(n() / chunksize),
                   k = chunksize,
                   length = n())) - 1L,
             , format = "d",
               digits = 1, flag = "0")) %>%
    ungroup()  %>%
    group_by(chrom, region_id) %>%
    do(write_delim(., path = paste0(
                          out_pref[.$chrom[1]],
                          .$region_id[1], ".txt")),delim="\t")
