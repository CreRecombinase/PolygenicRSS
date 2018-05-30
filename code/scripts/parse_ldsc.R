# save.image()
# stop()

library(purrr)
library(dplyr)
library(readr)
library(tidyr)
library(EigenH5)
library(SeqSupport)

inf <- snakemake@input[["logf"]]
true_f <- snakemake@input[["true_f"]]
tparam_df <- read_df_h5(true_f,"SimulationInfo")
# rdsf <- snakemake@input[["rdsf"]]
fgeneid <- snakemake@params[["fgeneid"]]
outf <- snakemake@output[["logf"]]

pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""]) %>% mutate(ttca=NA)
# tparam_df <- read_df_h5(rdsf,"SimulationInfo") %>% mutate(ttca=NA) %>% inner_join(pl) %>% select(-ttca)
#save.image()
res_df <- map2_dfr(inf,fgeneid,function(filen,genen){
    parse_ldsc_h2log(filen) %>%
        mutate(fgeneid = as.character(genen)) %>% filter(Variable!="Ratio")%>% select(-SD)
}
) %>% spread(Variable,Est) %>% mutate(ttca=NA) %>% inner_join(pl) %>% select(-ttca) %>% rename(pve=Total_Observed_scale_h2) %>% inner_join(tparam_df)

write_delim(res_df, path = outf, delim = "\t")
#warnings()
