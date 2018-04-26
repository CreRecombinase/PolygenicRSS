
# save.image()
# stop()
# setwd("~/Dropbox/PolygenicRSS/code/snakemake_files/")
# load(".RData")
# save.image()
# stop()
library(SeqSupport)

library(tidyverse)
library(progress)
library(SeqArray)
library(RSSp)

# tym <- read_matrix_h5(ymf,"trait","ymat")

gdsf <- snakemake@input[["gdsf"]]
subsnpf <- snakemake@input[["subsnpf"]]
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
beta_h5file <- as.character(snakemake@output[["beta_hf"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))

out_h5f <- snakemake@output[["h5f"]]

stopifnot(!is.null(mfgeneid))
snp_df <- read_delim(subsnpf,delim="\t")

n <- calc_N_h5(c(gdsf,"/","dosage"))
p <- nrow(snp_df)
tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n = n,p = p,mfgeneid) %>% mutate(n=n,p=p)
                                        # save.image()

ymat <- gen_sim_phenotype_h5(snp_df,gdsf,beta_h5file,tparam_df)
stopifnot(all(!is.na(ymat)))
EigenH5::write_matrix_h5(out_h5f,"trait","ymat",ymat)
EigenH5::write_df_h5(tparam_df,groupname = "SimulationInfo",out_h5f)

