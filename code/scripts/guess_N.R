library(readr)
library(fst)
library(tidyverse)

gwasfst <- snakemake@input[["gwasfst"]]
gwasf <- snakemake@input[["gwasf"]]
ogwasf <- snakemake@output[["ogwasf"]]
pl <- snakemake@params

param_df <- as_data_frame(pl[names(pl)!=""])%>% mutate(tNA=NA)

gwas_df <- read_fst(gwasfst)
if(is.null(gwas_df[["sample_size"]])){
    gwas_cols <-readr::cols_only("chrom"="i",
                        "pos"="i",
                        "snp"="c",
                        "ref_allele"="c",
                        "alt_allele"="c",
                        "beta_hat"="d",
                        "se"="d",
                        "sample_size"="d",
                        "p_value"="d")


    gwas_df <-read_delim(gwasif,delim="\t",col_types=gwas_cols)

}

guess_samp <- all(is.na(gwas_df[["sample_size"]])) || (any(gwas_df[["sample_size"]]<0))
if(guess_samp){

    optimise_f <- function(N,z,p_value){
        return(sum(abs(2*pt(abs(z),df=N,lower.tail=F)-p_value)))
    }

    calc_N <- function(z,p_value){
        return(optimise(optimise_f,interval=c(1,500000),z=z,p_value=p_value)$minimum)
    }


    gwas_df <- sample_n(gwas_df,50000,replace=F) %>% mutate(z=beta_hat/se,sample_size=map2_dbl(z,p_value,calc_N))
}

summ_df <- gwas_df %>%  summarise(
                            max_N=max(sample_size,na.rm=T),
                            median_N=median(sample_size,na.rm=T),
                            mean_N=mean(sample_size,na.rm=T),
                            var_N=var(sample_size,na.rm=T)
                        ) %>%
    mutate(tNA=NA) %>%
    inner_join(param_df) %>%
    select(-tNA)

write_delim(summ_df,ogwasf,delim="\t")
