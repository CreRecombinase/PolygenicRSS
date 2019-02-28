# save.image()
# stop()

library(tidyverse)

ldsc_noi <- read_tsv(snakemake@input[["ldsc_noif"]])%>% mutate(sigu=RSSp:::calc_sigu(pve = pve,p_n = p/n))
ldsc_i <- read_tsv(snakemake@input[["ldsc_if"]]) %>% mutate(sigu=RSSp:::calc_sigu(pve = pve,p_n = p/n))
rssp_res_df <- readRDS(snakemake@input[["rsspf"]]) %>% mutate(ncovar=as.integer(ncovar))

fit_parameters <-c("pop_a","pop_b","scenario","useLDetect","geneticMap","useLDshrink","ncovar","n","p","scenario")

nldsc_noi <- select(ldsc_noi,-one_of(fit_parameters),-Lambda_GC,-Mean_Chi_2,-Total_time_elapsed,-pcovar,-fgeneid,-useIntercept) %>% rename(bias=Intercept) %>%  mutate(method="LDSC_o",convergence=0,bias=0.0)
nldsc_i <- select(ldsc_i,-one_of(fit_parameters),-Lambda_GC,-Mean_Chi_2,-Total_time_elapsed,-pcovar,-fgeneid,-useIntercept) %>%rename(bias=Intercept) %>%  mutate(method="LDSC_i",convergence=0)

nrssp_res_df <- select(rssp_res_df,-one_of(fit_parameters),-bias,-nterms,-pcovar,-fgeneid,-lnZ,-trait,-simulation) %>% mutate(method="RSSp",bias=0,trait_id=as.integer(trait_id))

param_df <-bind_rows(
select(ldsc_i,one_of(fit_parameters)) %>% distinct() %>% mutate(method="LDSC_i",useLDetect=as.character(useLDetect),useLDshrink=as.character(useLDshrink)),
select(ldsc_noi,one_of(fit_parameters)) %>% distinct() %>% mutate(method="LDSC_o",useLDetect=as.character(useLDetect),useLDshrink=as.character(useLDshrink)),
select(rssp_res_df,one_of(fit_parameters)) %>% distinct() %>% mutate(method="RSSp",useLDetect=as.character(useLDetect),useLDshrink=as.character(useLDshrink)))


all_res <- bind_rows(nldsc_noi,nldsc_i,nrssp_res_df)

saveRDS(all_res,snakemake@output[["data_f"]])
saveRDS(param_df,snakemake@output[["param_f"]])



