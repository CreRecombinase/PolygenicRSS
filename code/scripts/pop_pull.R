library(tidyverse)


panel_pop  <- snakemake@params[["pop"]]
stopifnot(!is.null(panel_pop))
panel_url <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

panel_ind <- read_delim(file = panel_url,delim = "\t",trim_ws = T) %>% filter(pop==panel_pop) %>% select(sample) %>% write_delim(snakemake@output[["pop_list"]],delim="\t",col_names=F)


