library(EigenH5)
library(tidyverse)
pheno_h5f <- "../../../../Desktop/scratch/polyg_scratch/RSSp_sim_gwas_pheno/combined_wtcc_NoConfoundSmall_trait.h5"
indf <- "../../../../Desktop/scratch/polyg_scratch/plink/combined_19.grm.id"
pheno_of <- "../../../../Desktop/scratch/polyg_scratch/RSSp_sim_gwas_pheno/combined_wtcc_NoConfoundSmall_trait.pheno"

pheno_h5f <- snakemake@input[["phenof"]]
pheno_of <- snakemake@output[["phenof"]]
indf <- snakemake@input[["grmidbin"]]

info_df <- read_delim(indf,delim="\t",col_names=c("family","ID"))
stopifnot(!is.null(pheno_h5f),!is.null(pheno_of))
ymat <- read_matrix_h5(pheno_h5f,"trait","ymat")
#N <- nrow(ymat)
cbind(info_df,ymat) %>% as_data_frame() %>% write_delim(pheno_of,col_names=F,delim="\t")
