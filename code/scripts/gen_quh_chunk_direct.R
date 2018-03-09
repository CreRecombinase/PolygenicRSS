                                        #Rprof(filename=snakemake@output[["proff"]],append=F)

library(SeqSupport)
library(EigenH5)
library(tidyverse)
library(LDshrink)


evdf <- snakemake@input[["evdf"]]
traitf <- snakemake@input[["traitf"]]
pvv <- as.numeric(snakemake@params[["pvv"]])
svd <- snakemake@params[["svd"]]
if(is.null(svd)){
    svd <- "F"
}
svd <- svd=="T"
quhf <- snakemake@output[["quhf"]]


tparam_df <-read_df_h5(traitf,"SimulationInfo")

# D_df <- map_df(ld_grp,~data_frame(region_id=.x,D=read_vector_h5(evdf,paste0("EVD/",.x),"D")))
#It doesn't actually matter what order D is in
                                        #D <- purrr::map(ld_grp,~read_vector_h5(evdf,paste0("EVD/",.x),"D")) %>% purrr::flatten_dbl()
ld_gpn <- ifelse(svd,"SVD","EVD")
ld_dn <- ifelse(svd,"d","D")
ld_grp <-get_objs_h5(evdf,ld_gpn)
D_df <- map_df(ld_grp,~data_frame(region_id=as.integer(.x),D=read_vector_h5(evdf,paste0(ld_gpn,"/",.x),ld_dn))) %>% arrange(as.integer(region_id))
if(svd){
    D_df <- mutate(D_df,D=D^2)
}

sub_D_pvv <- filter_pvv(D_df,pvv)
D <- sub_D_pvv$D
#gw_snpi <- read_vec(evdf,"LDinfo/snp_id")
p <- length(D)
quh <- do.call("cbind",purrr::map(tparam_df$tsigu,~rnorm(n = p,mean = 0,sd = sqrt(.x^2*D^2+D))))


write_matrix_h5(quhf,
                "/",
                "quh",
                quh)

write_vector_h5(quhf,"/","D",D)

write_df_h5(tparam_df,groupname = "SimulationInfo",quhf)
write_df_h5(read_df_h5(evdf,"LDinfo"),"SNPinfo",quhf)

#warnings()
#Rprof(NULL)
