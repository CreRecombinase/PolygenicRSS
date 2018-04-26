library(tidyverse)

isnpf <- snakemake@input[["snp_bim"]]
bsnpf <- snakemake@output[["dup_snps"]]
snp_df <- read_delim(isnpf,delim="\t",col_names=c("chrom","SNP","b","pos","ref","alt"))
filter(snp_df,pos==0 | duplicated(paste0(chrom,":",pos))) %>% select(chrom,pos) %>% inner_join(snp_df) %>% select(SNP) %>% write_delim(bsnpf,col_names=F)

