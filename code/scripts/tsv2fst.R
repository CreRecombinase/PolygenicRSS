library(readr)
library(fst)
library(dplyr)
#library(broom)

gwasif <- snakemake@input[["gwasf"]]
gwasof <- snakemake@output[["gwasf"]]
#gwas_ct <- snakemake@output[["gwas_ct"]]

gwas_cols <-readr::cols_only("chrom"="i",
                        "pos"="i",
                        "snp"="c",
                        "ref_allele"="c",
                        "alt_allele"="c",
                        "beta_hat"="d",
                        "se"="d",
                        "sample_size"="d",
                        "p_value"="d")


gwas_df <-read_delim(gwasif,delim="\t",col_types=gwas_cols) %>%
    dplyr::distinct(chrom,pos,.keep_all=T)

## optimise_f <- function(N,beta_hat,se,p_value){
##     return(sum(abs(pt(beta_hat/se,df=N,lower.tail=FALSE)-p_value)))
## }

## sample_n(gwas_df,10000) %>%

write_fst(gwas_df,path=gwasof,compress=10)
