library(tidyverse)


panel_pop  <- snakemake@params[["pop"]]
stopifnot(!is.null(panel_pop))
panel_url <- snakemake@input[["panel_file"]]
output_f <- snakemake@output[["pop_list"]]
#panel_url <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"


panel_ind <- read_delim(file = panel_url,delim = "\t",trim_ws = T)
if(panel_pop %in% panel_ind$pop){
    sub_df <- panel_ind %>% filter(pop==panel_pop) %>% select(sample)
    stopifnot(nrow(sub_df)>0)
    write_delim(sub_df,,delim="\t",col_names=F)
}else{
    sub_df <- panel_ind %>% filter(panel_pop==super_pop) %>% select(sample)
    stopifnot(nrow(sub_df)>0)
    write_delim(sub_df,output_f,delim="\t",col_names=F)
}
