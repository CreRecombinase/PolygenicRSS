library(ldmap)
library(EigenH5)
library(rbgen)
library(dplyr)
input_f <- snakemake@input[["bgen"]]
output_f <- snakemake@output[["h5"]]
ldmr <- snakemake@params[["ldmr"]]
ldid <- ldetect_EUR[ldmr]
new_name <- magrittr::set_names(ldmap:::chromosome_levels(),stringr::str_pad(as.character(1:24),2,pad="0"))



range_df <- tibble(ldmr=ldid) %>% explode_ldmap_region() %>% transmute(chromosome=as.character(fct_recode(chrom,!!!new_name)),start,end)
data_r <- bgen.load(input_f,range_df)
