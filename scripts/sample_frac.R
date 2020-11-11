library(dplyr)

tibble(snps=scan(snakemake@input[["snp_list"]],what=character())) %>%
  sample_frac(size=as.numeric(snakemake@params[["sample_frac"]]),replace=FALSE) %>% 
  readr::write_tsv(snakemake@output[["snp_list"]],col_names=FALSE)
