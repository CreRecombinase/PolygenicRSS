                                        #Rprof(filename=snakemake@output[["proff"]],append=F)
# save.image()
# stop()
library(SeqSupport)

library(tidyverse)
library(progress)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
n_pgdsf <- snakemake@input[["n_pgdsf"]]
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
beta_infile <- snakemake@input[["beta_hf"]]
beta_h5file <- as.character(snakemake@output[["beta_hf"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))

out_h5f <- snakemake@output[["h5f"]]

stopifnot(!is.null(mfgeneid))

if(!is.null(n_pgdsf)){
  gds <- seqOpen(n_pgdsf, readonly = T)
  n <- length(seqGetData(gds, "sample.id"))
  p <- length(seqGetData(gds, "variant.id"))
  seqClose(gds)
  gds <- seqOpen(gdsf)
}else{

  n <- calc_N_h5(c(gdsf,"/","dosage"))
  p <- calc_p_h5(c(gdsf,"/","dosage"))
}
tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n = n,p = p,mfgeneid) %>% mutate(n=n,p=p)
                                        # save.image()
if(!is.null(beta_infile)){
  betamat <- read_2d_mat_h5(beta_infile,"/","beta")
  if(nrow(betamat)!=p){
    stopifnot(ncol(betamat)==p)
    betamat <- t(betamat)
  }
  stopifnot(ncol(betamat)==nrow(tparam_df))
}else{
  betamat <- NULL
}
ymat <- gen_sim_phenotype_h5(
    h5file=gdsf,
    tparam_df = tparam_df,
    fgeneid = mfgeneid,
    betamat=betamat,
    beta_h5file=beta_h5file)
EigenH5::write_matrix_h5(out_h5f,"trait","ymat",ymat)
EigenH5::write_df_h5(tparam_df,groupname = "SimulationInfo",outfile = out_h5f)
#Rprof(NULL)
