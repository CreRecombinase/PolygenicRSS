#                                         #Rprof(filename=snakemake@output[["proff"]],append=F)
# save.image()
# stop()

library(LDshrink)
library(tidyverse)
#library(RSSp)
library(EigenH5)

h5f <- snakemake@input[["h5f"]]
uhf <- snakemake@output[["uhf"]]
stopifnot(!is.null(h5f),!is.null(uhf))
ymatf <- normalizePath(snakemake@input[["ymatf"]])
stopifnot(!is.null(ymatf))
cores <- snakemake@threads

data("break_df")
break_dfl <- split(break_df,break_df$chr)
break_df <- map_df(break_dfl,function(x){
mutate(x,start=ifelse(start==min(start),0,start),stop=ifelse(stop==max(stop),.Machine$integer.max,stop))
})

snp_df <- read_df_h5(h5f,"SNPinfo")

ld_r <- LDshrink::set_ld_region(ld_regions = break_df,snp_info = snp_df)
#snp_id <- seqGetData(gds,var.name="variant.id")
tparam_df <- EigenH5::read_df_h5(ymatf, "SimulationInfo")
stopifnot(nrow(tparam_df)>0)
write_df_h5(snp_df,"SNPinfo",uhf)
write_df_h5(tparam_df,"SimulationInfo",uhf)
library(SeqSupport)
SeqSupport::map_eQTL_chunk_h5(h5f,ymatf,uhf,ld_r)

