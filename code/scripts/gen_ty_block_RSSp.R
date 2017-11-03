#Rprof(filename=snakemake@output[["proff"]],append=F)
library(SeqSupport)
library(RcppEigenH5)
library(tidyverse)
library(progress)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
beta_h5file <- as.character(snakemake@output[["beta_hf"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))

out_h5f <- snakemake@output[["h5f"]]



gds <- seqOpen(gdsf, readonly = T)
n <- length(seqGetData(gds, "sample.id"))
p <- length(seqGetData(gds, "variant.id"))

tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n = n,p = p,mfgeneid) %>% mutate(n=n,p=p)
# save.image()
ymat <- gen_sim_phenotype_gds(gds,tparam_df = tparam_df,fgeneid = mfgeneid,beta_h5file=beta_h5file)
write_mat_h5(out_h5f,"trait","ymat",ymat,deflate_level = 0L)
write_df_h5(tparam_df,groupname = "SimulationInfo",outfile = out_h5f,deflate_level = 0L)
#Rprof(NULL)
