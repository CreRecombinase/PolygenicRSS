library(SeqSupport)
library(EigenH5)
library(tidyverse)
library(LDshrink)


evdf <- snakemake@input[["evdf"]]
traitf <- snakemake@input[["traitf"]]


quhf <- snakemake@output[["quhf"]]


tparam_df <-read_df_h5(traitf,"SimulationInfo")

ld_gpn <- "EVD"
ld_dn <- "D"
ld_qn <- "Q"
ld_grp <-get_objs_h5(evdf,ld_gpn)
D_df <- map_df(ld_grp,~tibble(region_id=as.integer(.x),D=read_vector_h5(evdf,paste0(ld_gpn,"/",.x),ld_dn)))%>% arrange(region_id)


D <- zapsmall(D_df$D)
p <- length(D)

quh <- do.call("cbind",purrr::map(tparam_df$tsigu,~rnorm(n = p,mean = 0,sd = sqrt(.x^2*D^2+D))))
stopifnot(p==nrow(quh))
stopifnot(p==dim_h5(evdf,"LDinfo/pos"))

write_matrix_h5(quhf,
                "/",
                "quh",
                quh)
pl <- snakemake@wildcards
if(is.null(pl[["simulation"]])){
  pl[["simulation"]] <-"direct"
}
pl <- as_tibble(pl[names(pl)!=""])
write_df_h5(pl,groupname = "Wildcards",filename=quhf)

write_vector_h5(quhf,"/","D",D)

write_df_h5(tparam_df,"SimulationInfo",quhf)
write_df_h5(read_df_h5(evdf,"LDinfo"),"SNPinfo",quhf)

#warnings()
#Rprof(NULL)
warnings()
