library(tidyverse)

input_f <- snakemake@input[["input_files"]]

output_f <- snakemake@output[["output_file"]]
delim <- snakemake@params[["delim"]] %||% "\t"

map_df(input_f,read_delim,delim=delim) %>% write_delim(path=output_f,delim=delim)
