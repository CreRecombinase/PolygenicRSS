library(dplyr)
af <- as.numeric(snakemake@params[["snp_freq"]])
ind_df <- vroom::vroom(snakemake@input[["freqf"]]) %>% 
filter(between(ALT_FREQS,af,1-af)) %>% 
select(ID) %>% 
dplyr::count(ID) %>% 
filter(n==1) %>% 
select(ID)

vroom::vroom(snakemake@input[["freqpf"]]) %>% 
filter(between(ALT_FREQS,af,1-af)) %>% 
select(ID)  %>% 
dplyr::count(ID) %>% 
filter(n==1) %>% 
select(ID) %>% 
semi_join(ind_df,by="ID") %>% 
readr::write_tsv(snakemake@output[["snp_list"]],col_names=FALSE)
