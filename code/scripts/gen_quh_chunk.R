                                        #Rprof(filename=snakemake@output[["proff"]],append=F)

library(SeqSupport)

library(RSSp)

# save.image()
# stop()


evdf <- path.expand(snakemake@input[["evdf"]])
uhf <- path.expand(snakemake@input[["uhf"]])
LDchunk <- as.character(snakemake@params[["LDchunk"]])
perc_variance <- as.numeric(snakemake@params[["perc_variance"]])
if(length(perc_variance)==0){
  perc_variance <- 1.0
}
quhf <- snakemake@output[["quhf"]]

tparam_df <- read_df_h5(uhf,"SimulationInfo")

gw_snpi <- read_vec(uhf,paste0("/",LDchunk,"/snp_id"))
uh_chunk <- read_2d_mat_h5(uhf, LDchunk,"uh")
se_chunk <- read_2d_mat_h5(uhf,LDchunk,"se")
stopifnot(!any(is.na(uh_chunk)))
resl <- gen_quh_chunk_mat(uh_chunk,evdf,gw_snpi,perc_variance)

write_mat_h5(quhf,
             LDchunk,
             "quh",
             resl$quh,
             deflate_level=0,
             doTranspose=F)
write_mat_h5(quhf,
             LDchunk,
             "se",
             se_chunk,
             deflate_level=0,
             doTranspose=F)


write_vec(quhf,LDchunk,"D",resl$D)

write_vec(quhf,LDchunk,"snp_id",resl$snp_id)
write_df_h5(tparam_df,"SimulationInfo",quhf)
#Rprof(NULL)
