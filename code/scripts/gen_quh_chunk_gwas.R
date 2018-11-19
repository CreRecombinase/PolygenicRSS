library(tidyverse)
library(EigenH5)
library(RSSp)
library(glue)
library(progress)
input <- snakemake@input[["evdf"]]
output <- snakemake@output[["quhf"]]


snp_df <- read_df_h5(input,"LDinfo")

quh_fun <- function(x){
    read_matrix_h5(input,glue::glue("EVD/{x}/Q"))
}

snp_fun <- gensplit_uh_uf(snp_df$beta_hat/snp_df$se,snp_df$region_id)

regid <- unique(snp_df$region_id)
head(snp_df)
cat(regid)
pb <- dplyr::progress_estimated(length(regid))
Dvec <- unlist(map(glue("EVD/{regid}/D"),~{pb$tick()$print()
    read_vector_h5(input,.x)}))
rssp_df <- data_frame(quh=unlist(map_convert_quh(regid,snp_fun,quh_fun)),region_id=snp_df$region_id,D=Dvec)


pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""])
write_df_h5(pl, filename = output, datapath = "Wildcards")

"Writing"
write_df_h5(rssp_df,output,"quh_df")
