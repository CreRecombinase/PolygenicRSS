library(tidyverse)

library(EigenH5)


map_files <- snakemake@input[["mapf"]]
outf <- snakemake@output[["mapf"]]
chr <- as.integer(snakemake@params[["chrom"]])

stopifnot(length(chr)==length(map_files), !is.null(chr))
map_file_df <- tibble(map_file = map_files,
                          chr = chr)

map_df <- map2_dfr(map_file_df$map_file, map_file_df$chr,
                   function(filename, chrom){
                       tdf <- read_delim(filename,
                                  col_names = c("pos", "map"),
                                  delim = " ", col_types = c("_id")) %>% mutate(chr = chrom)
                       stopifnot(!is.unsorted(tdf$map))
                       return(tdf)
                   }) %>% distinct(map,.keep_all=T)


write_df_h5(map_df,outf, "SNPinfo")
