#Rprof(filename=snakemake@output[["proff"]],append=F)
library(SeqSupport)

library(RSSp)

# save.image()
# stop()


evdf <- snakemake@input[["evdf"]]
uhf <- snakemake@input[["uhf"]]
LDchunk <- as.character(snakemake@params[["LDchunk"]])
quhf <- snakemake@output[["quhf"]]

tparam_df <- read_df_h5(uhf,"SimulationInfo")

gw_snpi <- read_vec(uhf,paste0("/",LDchunk,"/snp_id"))
uh_chunk <- read_2d_mat_h5(uhf, LDchunk,"uh")
stopifnot(!any(is.na(uh_chunk)))
resl <- gen_quh_chunk_mat(uh_chunk,evdf,gw_snpi)

write_mat_h5(quhf,
             LDchunk,
             "quh",
             resl$quh,
             deflate_level=0,
             doTranspose=F)

write_vec(quhf,LDchunk,"D",resl$D)

write_vec(quhf,LDchunk,"snp_id",resl$snp_id)
write_df_h5(tparam_df,"SimulationInfo",quhf)
#Rprof(NULL)
