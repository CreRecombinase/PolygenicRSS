

library(LDshrink)
library(tidyverse)
library(SeqSupport)

library(EigenH5)

h5f <- snakemake@input[["h5f"]]
uhf <- snakemake@output[["uhf"]]
subsnpf <- snakemake@input[["subsnpf"]]
stopifnot(!is.null(h5f),!is.null(uhf),!is.null(subsnpf))
ymatf <- normalizePath(snakemake@input[["ymatf"]])
stopifnot(!is.null(ymatf))
cores <- snakemake@threads


snp_df <- read_delim(subsnpf,delim="\t")

tparam_df <- EigenH5::read_df_h5(ymatf, "SimulationInfo")

p <- nrow(snp_df)
N <- as.integer(unique(tparam_df$n))
g <- nrow(tparam_df)
chunksize <- 10000
num_chunks <- ceiling(p/chunksize)

snp_dff <- dplyr::mutate(snp_df,nchunk_id=sort(as.integer(gl(n = num_chunks,k=chunksize,length = p)))) %>%
    EigenH5::split_chunk_df(pos_id=snp_id,group_id=nchunk_id) %>%
    dplyr::mutate(chunk_group=nchunk_id) %>% mutate(filenames=h5f,groupnames="/",datanames="dosage")


snp_dff  <- BBmisc::chunk(snp_df$snp_id,chunk.size = chunksize) %>%
  map_df(cont_reg) %>%
  mutate(filenames=h5f,groupnames="/",datanames="dosage",snp_chunk=1:n(),
         row_offsets=in_start,row_chunksizes=in_stop-in_start+1,col_offsets=0L,col_chunksizes=-1L)


exp_dff  <- data_frame(filenames=ymatf,groupnames="trait",datanames="ymat",row_offsets=0L,row_chunksizes=-1L,col_offsets=0L,col_chunksizes=-1L)


uh_dff <- BBmisc::chunk(1:p,chunk.size = chunksize) %>%
  map_df(cont_reg) %>%
  mutate(filenames=uhf,groupnames="/",datanames="uh",
         row_offsets=in_start,row_chunksizes=in_stop-in_start+1,
         col_offsets=0L,col_chunksizes=-1L,chunk_group=1:n())
se_dff <- uh_dff %>% mutate(filenames=uhf,groupnames="/",datanames="se")

EigenH5::create_matrix_h5(uhf,"/","uh",numeric(),dims=c(nrow(snp_df),nrow(tparam_df)),chunksizes=c(1000L,nrow(tparam_df)))
EigenH5::create_matrix_h5(uhf,"/","se",numeric(),dims=c(nrow(snp_df),nrow(tparam_df)),chunksizes=c(1000L,nrow(tparam_df)))

stopifnot(nrow(tparam_df)>0)

write_df_h5(snp_df,"SNPinfo",uhf)
write_df_h5(tparam_df,"SimulationInfo",uhf)


SeqSupport::map_eQTL_chunk_h5(snp_dff,exp_dff,uh_dff,se_dff)

