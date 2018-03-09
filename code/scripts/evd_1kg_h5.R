library(EigenH5)
library(LDshrink)
library(tidyverse)




cutoff <- 1e-3
hdff <- snakemake@input[["hdff"]]
mapf <- snakemake@input[["mapf"]]

subsnpf <- snakemake@input[["subsnpf"]]
bdf <- snakemake@input[["bdf"]]
outf <- snakemake@output[["evdf"]]

cutoff <- as.numeric(snakemake@params[["cutoff"]])
if(length(cutoff)==0){
    cutoff <- formals(LDshrink::chunkwise_LDshrink_h5)[["cutoff"]]
}

stopifnot( !is.null(hdff), !is.null(outf),!is.null(mapf),!is.null(bdf))

stopifnot(file.exists(hdff), !file.exists(outf),file.exists(mapf),file.exists(bdf))
break_df <- read_delim(bdf,delim="\t")



snp_df <- read_delim(subsnpf,delim="\t")
snp_df <-assign_snp_block(snp_df,break_df,assign_all = F) %>% filter(!is.na(region_id))
p <- nrow(snp_df)

stopifnot(sorted_snp_df(snp_df))

map_df <- read_df_h5(mapf,"SNPinfo")
snp_df <- assign_map(snp_df,map_df)
rm(map_df)
stopifnot(group_by(snp_df,chr) %>% summarise(sorted=!is.unsorted(map)) %>% summarise(sorted=all(sorted)) %>% pull(1))

LDshrink::chunkwise_LDshrink_h5(input_file = hdff,output_file = outf,snp_df = snp_df,evd = T,svd = T,cutoff=cutoff,df=F)
