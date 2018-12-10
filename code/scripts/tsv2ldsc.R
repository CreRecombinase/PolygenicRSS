library(dplyr)
library(fst)
library(readr)

fst::read_fst(snakemake@input[["gwasf"]])  %>% mutate(Z=beta_hat/se) %>% select(SNP=snp,A1=ref_allele,A2=alt_allele,N=sample_size,Z) %>% write_delim(snakemake@output[["ldscf"]],delim="\t")
