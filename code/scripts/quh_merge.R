library(dplyr)
library(purrr)

quhf <- snakemake@input[["quhf"]]
pvv <- snakemake@params[["pvv"]] %||% 0.0


samp_df <- read_delim(samp_dff,delim="\t")
quh_df <- imap_dfr(quhf,~read_df_h5(.x,"quh_df") %>% mutate(tc=.y))

wc_df <- imap_dfr(quh_dff,~read_df_h5(.x,"Wildcards") %>% select(-chrom))  %>% distinct() %>% mutate(ttca=NA) %>% inner_join(samp_df) %>% inner_join(pl) %>% select(-ttca)
