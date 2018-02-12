library(SeqArray)
library(EigenH5)
#library(LDshrink)
library(SeqSupport)
library(tidyverse)
#data("break_df")
#inf <- "~/Desktop/scratch/polyg_scratch/gds/scombined_19.gds"
inf <- snakemake@input[["input_gds"]]
outf <- snakemake@output[["outf"]]
#mapf <- snakemake@input[["mapf"]]

tgds <- seqOpen(inf)
si_df <- read_SNPinfo_gds(tgds) %>% mutate(chr=as.integer(chr))
#stopifnot(LDshrink::sorted_snp_df(si_df))
SeqSupport::gds2hdf5(tgds,outf)
#snp_df <- read_df_h5(outf,"SNPinfo")
#ld_r <- LDshrink::set_ld_region(ld_regions = break_df,snp_info = snp_df)
#map_df <- read_df_h5(mapf,"SNPinfo")

#snp_dfl <- split(snp_df,snp_df$chr)
#map_dfl <- split(map_df,map_df$chr)

#snp_df <- map2_df(map_dfl,snp_dfl,~mutate(.y,map=LDshrink::interpolate_map(.x$map,.x$pos,.y$pos)))
#write_vector_h5(outf,"SNPinfo","map",snp_df$map)
#write_vector_h5(outf,"SNPinfo","region_id",ld_r)
