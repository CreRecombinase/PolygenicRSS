library(purrr)
library(dplyr)
library(readr)
library(tidyr)
library(SeqSupport)

inf <- snakemake@input[["logf"]]
fgeneid <- snakemake@params[["fgeneid"]]
outf <- snakemake@output[["logf"]]

#save.image()
res_df <- map2_dfr(inf,fgeneid,function(filen,genen){
    parse_ldsc_h2log(filen) %>%
        mutate(fgeneid = as.character(genen)) %>% filter(Variable!="Ratio")%>% select(-SD)
}
) %>% spread(Variable,Est)

write_delim(res_df, path = outf, delim = "\t")
