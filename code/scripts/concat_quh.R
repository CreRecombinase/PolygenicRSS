##Rproffilename=snakemake@output[["proff"]],append=F)
library(dplyr)
library(SeqSupport)
library(purrr)
library(progress)

inf <- snakemake@input[["evdf"]]
LDchunk <- as.character(snakemake@params[["LDchunk"]])
outf <- snakemake@output[["rdsf"]]

nchunks <- length(inf)
iresl <- list()
tparam_df <- read_df_h5(inf[1],"SimulationInfo")
D <- map2(inf, paste0(LDchunk,"/D"), read_vec) %>% as_vector("numeric")
snp_info<- map2_df(inf,LDchunk,function(file,LDc){
    data_frame(
        snp_id=SeqSupport::read_vec(file,paste0("/",LDc,"/snp_id")),
        region_id=LDc)})
write_df_h5(snp_info,"snp_info",outf,deflate_level=4L)
write_df_h5(tparam_df,"tparam_df",outf)
quh_mat <- RcppEigenH5::concat_mat_chunks(inf, LDchunk,
                             rep("quh", length(LDchunk)))
write_mat_h5(h5file=outf,dataname="quh",data=quh_mat,deflate_level=1L,doTranspose=T)

write_vec(outf,dataname="D",data=D,deflate_level=1L)
stopifnot(all.equal(sum(D),nrow(quh_mat)))
