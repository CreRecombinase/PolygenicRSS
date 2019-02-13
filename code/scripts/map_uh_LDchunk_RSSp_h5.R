library(LDshrink)
library(tidyverse)
library(SeqSupport)
library(EigenH5)

h5f <- snakemake@input[["h5f"]]
uhf <- snakemake@output[["uhf"]]
subsnpf <- snakemake@input[["subsnpf"]]
subgwasf <- snakemake@input[["subgwasf"]]

stopifnot(!is.null(h5f),!is.null(uhf),!is.null(subsnpf))
ymatf <- normalizePath(snakemake@input[["ymatf"]])
stopifnot(!is.null(ymatf))
cores <- snakemake@threads


snp_df <- read_delim(subsnpf,delim="\t")
ind_v <- readRDS(subgwasf)
tparam_df <- EigenH5::read_df_h5(ymatf, "SimulationInfo")

p <- nrow(snp_df)
N <- as.integer(unique(tparam_df$n))
stopifnot(length(ind_v)==N)
g <- as.integer(nrow(tparam_df))
expdims <- dim_h5(ymatf,"trait/ymat")
stopifnot(g==expdims[1],
          length(ind_v)==expdims[2])


chunksize <- as.integer(50000)
num_chunks <- ceiling(p/chunksize)
cat("Chunking SNP data\n")
snp_df <- mutate(snp_df,snp_chunk=gl(n = num_chunks,k=chunksize,length = p),snp_chunk_id=1:n())
snp_l <-  split(select(snp_df,snp_id,snp_chunk_id),snp_df$snp_chunk)

snp_lff  <- snp_l %>% map(~list(subset_rows=.x$snp_id,
                                filename=h5f,
                                subset_cols=ind_v,
                                datapath="dosage"))
exp_lff <- list(list(filename=ymatf,datapath="trait/ymat"))

uh_lff <-  snp_l %>% map(~list(subset_rows=.x$snp_chunk_id,filename=uhf,datapath="uh"))
se_lff <-  map(uh_lff,~update_list(.,datapath="se"))

cat("Creating output matrices\n")
EigenH5::create_matrix_h5(uhf,"/","uh",numeric(),dims=c(p,g),chunksizes=c(1000L,g))
EigenH5::create_matrix_h5(uhf,"/","se",numeric(),dims=c(p,g),chunksizes=c(1000L,g))

stopifnot(nrow(tparam_df)>0)
stopifnot(!is.null(snp_df$AF),
          all(!is.na(snp_df$AF)))

cat("Writing simulation/data info\n")
write_df_h5(snp_df,"SNPinfo",uhf)
write_df_h5(tparam_df,"SimulationInfo",uhf)
pl <- snakemake@wildcards
pl <- as_tibble(pl[names(pl)!=""])
write_df_h5(pl,groupname = "Wildcards",filename=uhf)

cat("Mapping traits\n")
SeqSupport::map_eQTL_chunk_h5(snp_lff,exp_lff,uh_lff,se_lff)


cat("Checking uh\n")
tuh <- EigenH5::read_matrix_h5(uhf,"/","uh",
                               subset_rows=sort(sample(1:p,min(p,100),replace=F)),
                               subset_cols=sort(sample(1:g,min(g,100),replace=F)))
stopifnot(all(!is.na(c(tuh))))
