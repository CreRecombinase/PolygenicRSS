library(EigenH5)
library(LDshrink)
library(tidyverse)
data("break_df")
inf <- snakemake@input[["inf"]]
mapf <- snakemake@input[["mapf"]]
outf <- snakemake@output[["outf"]]
snp_df <- read_df_h5(inf,"SNPinfo")
ld_r <- LDshrink::set_ld_region(ld_regions = break_df,snp_info = snp_df)
map_df <- read_df_h5(mapf,"SNPinfo")


snp_dfl <- split(snp_df,snp_df$chr)
map_dfl <- split(map_df,map_df$chr)
snp_df <- map2_df(map_dfl,snp_dfl,~mutate(.y,map=LDshrink::interpolate_map(.x$map,.x$pos,.y$pos)))
file.copy(inf,outf)
write_vector_h5(outf,"SNPinfo","map",snp_df$map)
write_vector_h5(outf,"SNPinfo","region_id",ld_r)


