#Rprof(filename=snakemake@output[["proff"]],append=F)
library(dplyr)
library(SeqSupport)
library(purrr)
library(progress)

inf <- snakemake@input[["evdf"]]
LDchunk <- as.character(snakemake@params[["LDchunk"]])
outf <- snakemake@output[["rdsf"]]

nchunks <- length(inf)
iresl <- list()
snp_id<- map2(inf,paste0(LDchunk,"/snp_id"),SeqSupport::read_vec) %>% as_vector("integer")
write_vec(outf,dataname="snp_id",data=snp_id,deflate_level=1L)

tparam_df <- read_df_h5(inf[1],"SimulationInfo")
write_df_h5(tparam_df,"tparam_df",outf)
uh_mat <- RcppEigenH5::concat_mat_chunks(inf, LDchunk,
                                        rep("uh", length(LDchunk)))
write_mat_h5(h5file=outf,dataname="uh",data=uh_mat,deflate_level=1L,doTranspose=T)
#saveRDS(iresl, outf)
#Rprof(NULL)
