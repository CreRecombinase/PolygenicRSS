# save.image("sim.RData")
# stop()
library(SeqSupport)
library(tidyverse)
library(progress)
library(EigenH5)
library(RSSp)
library(glue)
gdsf <- snakemake@input[["gdsf"]]
subsnpf <- snakemake@input[["subsnpf"]]
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
beta_h5file <- as.character(snakemake@output[["beta_hf"]])
evd_h5file  <- snakemake@input[["evdf"]]
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
chrom <- snakemake@params[["chrom"]]
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))

snpctf <- snakemake@input[["snpctf"]]
stopifnot(all(file.exists(snpctf)))
snpct <- sum(map_int(snpctf,~scan(.x,what=integer())))

out_h5f <- snakemake@output[["h5f"]]

stopifnot(!is.null(mfgeneid))
if(file.exists(beta_h5file)){
  file.remove(beta_h5file)
}
stopifnot(!file.exists(beta_h5file),file.exists(evd_h5file))
#snp_df <- read_delim(subsnpf,delim="\t")
snp_df <- read_df_h5(evd_h5file,"LDinfo")


n <- dim_h5(gdsf,"SampleInfo/sample_id")

p <- snpct
tparam_df <- tparamdf_norm(pve,bias,nreps,n = n,p = p,mfgeneid) %>% mutate(n=n,p=p)


#stopifnot(!is.null(snp_df$AF),all(!is.na(snp_df$AF)))
g <- nrow(tparam_df)
# save.image("sim.RData")
# stop()
ymat <- gen_ty_h5(snp_df = snp_df,
                  snp_h5file = gdsf,
                  beta_h5file = beta_h5file,
                  tparam_df = tparam_df,
                  sim_ind = 1:n)
if(ncol(ymat)==g && nrow(ymat)==n){
    ymat <- t(ymat)
}
stopifnot(nrow(ymat)==g,
          ncol(ymat)==n)
stopifnot(all(!is.na(ymat)))

EigenH5::write_matrix_h5(ymat,out_h5f,glue("genetic_trait/genetic_ymat"))
EigenH5::write_df_h5(tparam_df,out_h5f,"SimulationInfo")
EigenH5::write_df_h5(read_df_h5(gdsf,"SampleInfo"),out_h5f,"SampleInfo")

pl <- snakemake@wildcards
pl <- as_tibble(pl[names(pl)!=""])
EigenH5::write_df_h5(pl,filename=out_h5f, "Wildcards")
