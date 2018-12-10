library(dplyr)
library(purrr)
library(feather)


out_tf <- snakemake@input[["outf_tf"]]
id_f <- snakemake@output[["id_f"]]
map_df(out_tf,~feather::read_feather(.x,columns=c("gwas_snp_id","chrom","region_id"))) %>% arrange(gwas_snp_id) %>% feather::write_feather(path=id_f)
