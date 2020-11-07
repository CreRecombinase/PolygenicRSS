library(dplyr)
vroom::vroom(snakemake@input[["gwasf"]]) %>% summarise(samplesize=mean(n_complete_samples,na.rm=TRUE)) %>% pull(samplesize) %>% saveRDS(snakemake@output[["outputf"]])
