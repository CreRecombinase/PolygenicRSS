#                                         #Rprof(filename=snakemake@output[["proff"]],append=F)
# save.image()
# stop()

library(LDshrink)
library(tidyverse)
library(SeqSupport)
#library(RSSp)
library(EigenH5)

h5f <- snakemake@input[["h5f"]]
uhf <- snakemake@output[["uhf"]]
stopifnot(!is.null(h5f),!is.null(uhf))
ymatf <- normalizePath(snakemake@input[["ymatf"]])
stopifnot(!is.null(ymatf))
cores <- snakemake@threads

# data("break_df")
# break_dfl <- split(break_df,break_df$chr)
# break_df <- map_df(break_dfl,function(x){
# mutate(x,start=ifelse(start==min(start),0,start),stop=ifelse(stop==max(stop),.Machine$integer.max,stop))
# })

snp_df <- read_df_h5(h5f,"SNPinfo")


# ld_r <- LDshrink::set_ld_region(ld_regions = break_df,snp_info = snp_df)
#snp_id <- seqGetData(gds,var.name="variant.id")
tparam_df <- EigenH5::read_df_h5(ymatf, "SimulationInfo")


exp_dff <- chunk_df_h5(filename=ymatf,groupname="trait",dataname="ymat")
snp_dff <- chunk_df_h5(filename=h5f,groupname="/",dataname = "dosage",chunksize_row = 10000)
out_dff <- gen_outer_df(exp_dff,snp_dff)
uh_dff <- out_dff %>% mutate(filenames=uhf,groupnames="/",datanames="uh")
se_dff <- out_dff %>% mutate(filenames=uhf,groupnames="/",datanames="se")
EigenH5::create_matrix_h5(uhf,"/","uh",numeric(),dims=c(nrow(snp_df),nrow(tparam_df)),chunksizes=c(1000,nrow(tparam_df)))
EigenH5::create_matrix_h5(uhf,"/","se",numeric(),dims=c(nrow(snp_df),nrow(tparam_df)),chunksizes=c(1000,nrow(tparam_df)))

stopifnot(nrow(tparam_df)>0)
write_df_h5(snp_df,"SNPinfo",uhf)
write_df_h5(tparam_df,"SimulationInfo",uhf)


SeqSupport::map_eQTL_chunk_h5(snp_dff,exp_dff,uh_dff,se_dff)

