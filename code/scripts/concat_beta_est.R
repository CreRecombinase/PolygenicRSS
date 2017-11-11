##Rproffilename=snakemake@output[["proff"]],append=F)
library(dplyr)
library(SeqSupport)
library(purrr)
library(progress)

inf <- snakemake@input[["beta_hf"]]
LDchunk <- snakemake@params[["LDchunk"]]
outf <- snakemake@output[["nbeta_hf"]]

nchunks <- length(inf)

beta_mat <- RcppEigenH5::concat_mat_chunks(
                             inf,
                             rep("/",length(LDchunk)),
                             rep("beta_est", length(LDchunk)))

write_mat_h5(h5file=outf,dataname="beta_est",data=beta_mat,deflate_level=1L,doTranspose=T)
rm(beta_mat)

beta_mat <- RcppEigenH5::concat_mat_chunks(
                             inf,
                             rep("/",length(LDchunk)),
                             rep("beta_oracle", length(LDchunk)))

write_mat_h5(h5file=outf,dataname="beta_oracle",data=beta_mat,deflate_level=1L,doTranspose=T)
rm(beta_mat)



gc()
ymat <- RcppEigenH5::sum_mats(inf,rep("/",length(LDchunk)),rep("y_est",nchunks))
write_mat_h5(h5file=outf,dataname="y_est",data=ymat,deflate_level=1L,doTranspose=T)
ymat <- RcppEigenH5::sum_mats(inf,rep("/",length(LDchunk)),rep("y_oracle",nchunks))
write_mat_h5(h5file=outf,dataname="y_oracle",data=ymat,deflate_level=1L,doTranspose=T)


