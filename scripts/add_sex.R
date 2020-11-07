library(dplyr)
library(readr)
read_delim(snakemake@input[["onpf"]], delim=" ") %>% 
  inner_join(read_delim(snakemake@input[["samplef"]], delim=" ")) %>%
  write_delim(snakemake@output[["onpf"]], delim=" ")
