library(SeqSupport)
library(tidyverse)
library(progress)
library(RSSp)
#
# load("ty.RData")
# save.image("ty.RData")
# stop()
# tym <- read_matrix_h5(ymf,"trait","ymat")

gdsf <- snakemake@input[["gdsf"]]
subsnpf <- snakemake@input[["subsnpf"]]
subgwasf <- snakemake@input[["subgwasf"]]
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
beta_h5file <- as.character(snakemake@output[["beta_hf"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))

out_h5f <- snakemake@output[["h5f"]]

stopifnot(!is.null(mfgeneid))
snp_df <- read_delim(subsnpf,delim="\t")
ind_v <- readRDS(subgwasf)

n <- length(ind_v)
p <- nrow(snp_df)
tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n = n,p = p,mfgeneid) %>% mutate(n=n,p=p)
                                        # save.image()

stopifnot(!is.null(snp_df$AF),all(!is.na(snp_df$AF)))
g <- nrow(tparam_df)
ymat <- gen_sim_phenotype_h5(snp_df,gdsf,beta_h5file,tparam_df,AF=numeric(),ind_v)
if(ncol(ymat)==g && nrow(ymat)==n){
    ymat <- t(ymat)
}
stopifnot(nrow(ymat)==g,
          ncol(ymat)==n)
stopifnot(all(!is.na(ymat)))
EigenH5::write_matrix_h5(out_h5f,"trait","ymat",ymat)
EigenH5::write_df_h5(tparam_df,groupname = "SimulationInfo",out_h5f)

pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""])
EigenH5::write_df_h5(pl,groupname = "Wildcards",filename=out_h5f)

