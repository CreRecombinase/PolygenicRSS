##Rproffilename=snakemake@output[["proff"]],append=F)

# save.image()
# stop()
library(dplyr)
#library(SeqSupport)
library(purrr)
library(progress)
library(EigenH5)
inf <- snakemake@input[["evdf"]]

outf <- snakemake@output[["rdsf"]]

nchunks <- length(inf)
tparam_df <- EigenH5::read_df_h5(inf,"SimulationInfo")
D <- map2(inf, LDchunk, EigenH5::read_vector_h5,dataname="D") %>% as_vector("numeric")
snp_info <- read_df_h5(inf,"SNPinfo")
u_reg_id <- as.character(sort(unique(snp_info$region_id)))
write_df_h5(snp_info,"snp_info",outf)
write_df_h5(tparam_df,"tparam_df",outf)


quh_mat <- RcppEigenH5::concat_mat_chunks(inf, LDchunk,
                                          rep("quh", length(LDchunk)))

write_matrix_h5(filename=outf,groupname="/",dataname="quh",data=quh_mat,doTranspose=F)
rm(quh_mat)
gc()
#se_mat <- RcppEigenH5::concat_mat_chunks(inf, LDchunk,
                                         #rep("se", length(LDchunk)))

#write_mat_h5(h5file=outf,dataname="se",data=se_mat,deflate_level=1L,doTranspose=T)

write_vector_h5(outf,"/",dataname="D",data=D)

