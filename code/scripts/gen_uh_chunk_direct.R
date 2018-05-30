library(tidyverse)
library(RSSp)
library(LDshrink)
library(EigenH5)


evdf <- snakemake@input[["evdf"]]
quhf <- snakemake@input[["quhf"]]
uhf <- snakemake@output[["uhf"]]

stopifnot(!is.null(quhf),
          !is.null(evdf),
          file.exists(quhf),
          !file.exists(uhf))


snp_df <- EigenH5::read_df_h5(evdf,"LDinfo")
tparam_df <- EigenH5::read_df_h5(quhf,"SimulationInfo")

EigenH5::write_df_h5(tparam_df,"SimulationInfo",uhf)

pl <- snakemake@wildcards
if(is.null(pl[["simulation"]])){
  pl[["simulation"]] <-"direct"
}
pl <- as_data_frame(pl[names(pl)!=""])
write_df_h5(pl,groupname = "Wildcards",filename=uhf)

ld_grp <-ls_h5(evdf,"EVD")
D_df <- map_df(ld_grp,~data_frame(region_id=as.integer(.x),D=read_vector_h5(evdf,paste0("EVD/",.x),"D"))) %>%
    arrange(as.integer(region_id))


stopifnot(nrow(snp_df)==nrow(D_df))


Qf_dff <- EigenH5::gen_matslice_df(evdf,"EVD","Q")

p <- nrow(snp_df)


uh_dff <- mutate(snp_df,snp_id=1:n()) %>%
    select(snp_id,region_id) %>%
    split_chunk_df(snp_id,region_id,rowsel=T,colsel=F) %>%
    mutate(filenames=uhf,groupnames="/",datanames="uh")


quh_dff <-mutate(uh_dff,filenames=quhf,
                 groupnames="/",
                 datanames="quh")


#stopifnot(slice(quh_dff,n()) %>% mutate(lrow=row_offsets+row_chunksizes) %>% pull(lrow)==nrow(sub_D_pvv) )
#stopifnot(is.numeric(sub_D_pvv$D))
#EigenH5::write_vector_h5(quhf,"/","D",D_df$D)
create_matrix_h5(uhf,"/","uh",numeric(),dims=c(nrow(D_df),nrow(tparam_df)),chunksizes = c(pmin(1024,nrow(D_df)),1))




EigenH5::write_df_h5(snp_df,"SNPinfo",uhf)
SeqSupport::crossprod_quh_h5(Qf_dff,quh_dff,uh_dff,doTranspose=T)



#Rprof(NULL)
