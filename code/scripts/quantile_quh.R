library(EigenH5)
library(tidyverse)

inf <- snakemake@input[["quhf"]]
stopifnot(!is.null(inf))
outf <- snakemake@output[["quantf"]]

tparam_df <- read_df_h5(inf,"SimulationInfo")

write_df_h5(tparam_df,"SimulationInfo",outf)

quh_m  <- read_matrix_h5(inf,"/","quh")
D  <- read_vector_h5(inf,"/","D")

pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""])
write_df_h5(pl,groupname = "Wildcards",filename=outf)
tt_sd <- sqrt(rep(tparam_df$tsigu^2,length(D))*D^2+D)
quantile_m <- matrix(NA,nrow(quh_m),ncol(quh_m))

quantile_m[] <- pnorm(q = quh_m[],mean=0,sd=tt_sd)
# iwalk(tparam_df$tsigu,dbp)
stopifnot(!any(is.na(quantile_m)))
write_matrix_h5(outf,"/","quantile",data=quantile_m)
