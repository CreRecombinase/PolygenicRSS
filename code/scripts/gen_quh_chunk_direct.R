                                        #Rprof(filename=snakemake@output[["proff"]],append=F)

library(SeqSupport)
library(RSSp)
library(dplyr)
library(rhdf5)
#


stopifnot(!is.null(snakemake@params[["N"]]))
stopifnot(!is.null(snakemake@params[["P"]]))
# stop()


evdf <- snakemake@input[["evdf"]]
LDchunk <- as.character(snakemake@params[["LDchunk"]])
quhf <- snakemake@output[["quhf"]]


n <- as.integer(as.numeric(snakemake@params[["N"]]))
p <- as.integer(as.numeric(snakemake@params[["P"]]))



mfgeneid <- as.character(snakemake@params[["fgeneid"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))



tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n = n,p = p,mfgeneid) %>% mutate(n=n,p=p)


D <- read_vec(evdf,"EVD/D")

gw_snpi <- read_vec(evdf,"LDinfo/snp_id")


resl <- sim_quh_dir_df(tparam_df,D=D,seed=NULL,snp_id=gw_snpi)

write_mat_h5(quhf,
             LDchunk,
             "quh",
             resl$quh,
             deflate_level=0,
             doTranspose=F)

write_vec(quhf,LDchunk,"D",resl$D)

write_vec(quhf,LDchunk,"snp_id",resl$snp_id)
write_df_h5(tparam_df,"SimulationInfo",quhf)
warnings()
#Rprof(NULL)
