library(SeqSupport)
library(tidyverse)

inputf <- snakemake@input[["input_h5"]]
outputf <- snakemake@output[["output_df"]]
stopifnot(!is.null(inputf),!is.null(outputf))

read_df_h5(inputf,"SNPinfo") %>% write_delim(outputf,delim="\t")
