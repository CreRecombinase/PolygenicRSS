#Rprof(filename=snakemake@output[["proff"]],append=F)
library(EigenH5)
library(LDshrink)
library(tidyverse)


data("break_df")
hdff <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/h5/1KG19_seq_1kg_haplo.h5"


mapf <- "/run/media/nwknoblauch/Data/1kg/1000-genomes-genetic-maps/interpolated_hapmap.h5"
chunk_snps <- 0
outf <- paste0("/home/nwknoblauch/Desktop/scratch/polyg_scratch/EVD_H5/1KG19_1kg_hapmap_0.h5")
if(file.exists(outf)){
  file.remove(outf)
}
cutoff <- 1e-3
hdff <- snakemake@input[["hdff"]]
mapf <- snakemake@input[["mapf"]]
outf <- snakemake@output[["evdf"]]
chunk_snps <- as.integer(snakemake@params[["chunk_snps"]])
cutoff <- as.numeric(snakemake@params[["cutoff"]])
if(length(cutoff)==0){
    cutoff <- formals(LDshrink::chunkwise_LDshrink_h5)[["cutoff"]]
}
if(length(chunk_snps)==0){
  chunk_snps <- 0
}

# nbreak_df <- read_delim("https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/EUR/fourier_ls-chr1.bed",delim="\t")
stopifnot( !is.null(hdff), !is.null(outf),!is.null(mapf))

stopifnot(file.exists(hdff), !file.exists(outf),file.exists(mapf))
# break_df <-spread_ld_region(break_df)
snp_df <- read_df_h5(hdff,"SNPinfo")
snp_df <-assign_snp_block(snp_df,break_df,assign_all = F) %>% filter(!is.na(region_id))
p <- nrow(snp_df)
# snp_df <-collapse_ld_region(snp_df,break_df,min_block_size=chunk_snps)
stopifnot(sorted_snp_df(snp_df))

map_df <- read_df_h5(mapf,"SNPinfo")
snp_df <- assign_map(snp_df,map_df)
rm(map_df)


LDshrink::chunkwise_LDshrink_h5(input_file = hdff,output_file = outf,snp_df = snp_df,evd = F,svd = F,cutoff=cutoff,df=T)
