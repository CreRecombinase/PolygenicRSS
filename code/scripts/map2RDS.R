library(tidyverse)
library(EigenH5)


map_files <- snakemake@input[["mapf"]]
outf <- snakemake@output[["mapf"]]
chr <- as.integer(snakemake@params[["chrom"]])
stopifnot(length(chr)==length(map_files),!is.null(chr))
map_file_df <- data_frame(map_file = map_files,
                          chr = chr)

map_df <- map2_dfr(map_file_df$map_file, map_file_df$chr,
                   function(filename, chrom){
                       read_delim(filename,
                                  col_names = c("SNP", "pos", "map"),
                                  delim = " ",col_types=c("_in")) %>% mutate(chr = chrom)
                   })


write_df_h5(map_df,"SNPinfo",outf)
