
#
# save.image()
# # stop()
# load(".RData")
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
g <- as.integer(nrow(tparam_df))
chunksize <- as.integer(10000)
num_chunks <- ceiling(p/chunksize)
cat("Chunking SNP data\n")
snp_lff  <- BBmisc::chunk(snp_df$snp_id,chunk.size=chunksize) %>% map(~list(subset_rows=.x,filename=h5f,datapath="dosage"))
exp_lff <- list(list(filename=ymatf,datapath="trait/ymat"))

## exp_dff  <- data_frame(filenames=ymatf,groupnames="trait",datanames="ymat",row_offsets=0L,row_chunksizes=N,col_offsets=0L,col_chunksizes=g)

uh_lff <-  BBmisc::chunk(1:p,chunk.size=chunksize) %>% map(~list(subset_rows=.x,filename=uhf,datapath="uh"))
se_lff <-  map(uh_lff,~update_list(.,datapath="se"))

cat("Creating output matrices\n")
EigenH5::create_matrix_h5(uhf,"/","uh",numeric(),dims=c(p,g),chunksizes=c(1000L,g))
EigenH5::create_matrix_h5(uhf,"/","se",numeric(),dims=c(p,g),chunksizes=c(1000L,g))

stopifnot(nrow(tparam_df)>0)

cat("Writing simulation/data info\n")
write_df_h5(snp_df,"SNPinfo",uhf)
write_df_h5(tparam_df,"SimulationInfo",uhf)

cat("Mapping traits\n")
SeqSupport::map_eQTL_chunk_h5(snp_lff,exp_lff,uh_lff,se_lff,EXP_first = F,SNP_first = T)

cat("Checking uh\n")
tuh <- EigenH5::read_matrix_h5(uhf,"/","uh",
                               subset_rows=sort(sample(1:p,min(p,100),replace=F)),
                               subset_cols=sort(sample(1:g,min(g,100),replace=F)))
stopifnot(all(!is.na(c(tuh))))
