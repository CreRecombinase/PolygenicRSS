library(dplyr)
library(purrr)
af <- as.numeric(snakemake@params[["snp_freq"]])
max_p <- snakemake@params[["max_p"]] 
idf <- vroom::vroom(snakemake@input[["freqf"]]) %>% 
  filter(`#CHROM`==as.integer(snakemake@params[["chrom"]] %||% 19),
         between(ALT_FREQS,af,1-af))
idf %>% 
slice(seq_len(as.integer(max_p %||% nrow(idf)))) %>% 
select(ID) %>% 
readr::write_tsv(snakemake@output[["snp_list"]],col_names=FALSE)
