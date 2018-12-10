library(readr)
library(fst)
library(tidyverse)
library(dplyr)
#library(broom)

gwasif <- snakemake@input[["gwasf"]]
gwasof <- snakemake@output[["gwasf"]]
#gwas_ct <- snakemake@output[["gwas_ct"]]

gwas_cols <-readr::cols_only("variant"="c",
                        "rsid"="c",
                        "nCompleteSamples"="i",
                        "AC"="d",
                        "ytx"="c",
                        "beta"="d",
                        "se"="d",
                        "tstat"="d",
                        "pval"="d")


gwas_df <-read_delim(gwasif,delim="\t",col_types=gwas_cols) %>%
    separate(
    dplyr::distinct(chrom,pos,.keep_all=T)

## optimise_f <- function(N,beta_hat,se,p_value){
##     return(sum(abs(pt(beta_hat/se,df=N,lower.tail=FALSE)-p_value)))
## }

## sample_n(gwas_df,10000) %>%

write_fst(gwas_df,path=gwasof,compress=10)
