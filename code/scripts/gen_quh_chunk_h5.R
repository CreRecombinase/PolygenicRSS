                                        #Rprof(filename=snakemake@output[["proff"]],append=F)

#library(SeqSupport)
library(tidyverse)
library(RSSp)




evdf <- normalizePath(snakemake@input[["evdf"]])
uhf <- normalizePath(snakemake@input[["uhf"]])
quhf <- snakemake@output[["quhf"]]
stopifnot(!is.null(quhf))
quhf <- file.path(normalizePath(dirname(quhf)),basename(quhf))
perc_variance <- as.numeric(snakemake@params[["perc_variance"]])
if(length(perc_variance)==0){
  perc_variance <- 1.0
}

stopifnot(!is.null(perc_variance))

snp_df <- EigenH5::read_df_h5(evdf,"LDinfo")
ld_id <- as.character(sort(unique(snp_df$region_id)))
stopifnot(length(ld_id)>0)
D <- map(paste0("EVD/",ld_id),~EigenH5::read_vector_h5(evdf,.x,"D")) %>% as_vector()
stopifnot(is.numeric(D))
EigenH5::write_vector_h5(quhf,"/","D",D)

tparam_df <- EigenH5::read_df_h5(uhf,"SimulationInfo")

EigenH5::write_df_h5(tparam_df,"SimulationInfo",quhf)
EigenH5::write_df_h5(snp_df,"SNPinfo",quhf)
SeqSupport::crossprod_chunk_h5(evdf,uhf,ld_id,quhf,"Q","uh","/","quh")


#Rprof(NULL)
