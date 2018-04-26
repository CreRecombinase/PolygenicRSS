library(tidyverse)
library(RSSp)
library(LDshrink)
library(EigenH5)


evdf <- snakemake@input[["evdf"]]
quhf <- snakemake@output[["quhf"]]
uhf <- snakemake@input[["uhf"]]
stopifnot(!is.null(quhf),
          !is.null(evdf),
          file.exists(evdf),
          file.exists(uhf),
          !file.exists(quhf),
          file.exists(uhf))


snp_df <- EigenH5::read_df_h5(evdf,"LDinfo")
snp_df_u <- EigenH5::read_df_h5(uhf,"SNPinfo",subcols=c("chr","pos","allele"))

stopifnot(nrow(snp_df)==nrow(snp_df_u))
stopifnot(all(snp_df$chr==snp_df_u$chr),
              all(snp_df$pos==snp_df_u$pos),
              all(snp_df$allele==snp_df_u$allele))

tparam_df <- EigenH5::read_df_h5(uhf,"SimulationInfo")


p <- nrow(snp_df)
uh_d <- get_dims_h5(uhf,"uh")
stopifnot(uh_d[1]==p,
          uh_d[2]==nrow(tparam_df))



EigenH5::write_df_h5(tparam_df,"SimulationInfo",quhf)

ld_grp <-ls_h5(evdf,"EVD")
D_df <- map_df(ld_grp,~data_frame(region_id=as.integer(.x),D=read_vector_h5(evdf,paste0("EVD/",.x),"D"))) %>%
    arrange(as.integer(region_id))
LDsnp_df <- map_df(ld_grp,~data_frame(region_id=as.integer(.x),
                                      pos=read_vector_h5(evdf,paste0("LDi/",.x),"pos"),
                                      chr=read_vector_h5(evdf,paste0("LDi/",.x),"chr"),
                                      allele=read_vector_h5(evdf,paste0("LDi/",.x),"allele")
                                      )) %>%
    arrange(as.integer(region_id))

stopifnot(all(LDsnp_df$chr==snp_df_u$chr),
              all(LDsnp_df$pos==snp_df_u$pos),
              all(LDsnp_df$allele==snp_df_u$allele))


stopifnot(nrow(snp_df)==nrow(D_df))

Qf_dff <- EigenH5::gen_matslice_df(evdf,"EVD","Q")

uh_dff <- mutate(snp_df,
                snp_id=1:n()) %>%
    select(snp_id,region_id) %>%
    split_chunk_df(snp_id,region_id,rowsel=T,colsel=F) %>%
    mutate(filenames=uhf,
           groupnames="/",
           datanames="uh")


quh_dff <-mutate(uh_dff,filenames=quhf,
                 groupnames="/",
                 datanames="quh")

EigenH5::write_vector_h5(quhf,"/","D",D_df$D)
create_matrix_h5(quhf,"/","quh",numeric(),dims=c(nrow(D_df),nrow(tparam_df)),chunksizes = c(pmin(1024,nrow(D_df)),1))




EigenH5::write_df_h5(snp_df,"SNPinfo",quhf)
SeqSupport::crossprod_quh_h5(Qf_dff,uh_dff,quh_dff)



#Rprof(NULL)
