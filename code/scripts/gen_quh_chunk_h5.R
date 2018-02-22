                                        #Rprof(filename=snakemake@output[["proff"]],append=F)

#library(SeqSupport)
library(tidyverse)
library(RSSp)
library(LDshrink)
library(EigenH5)


evdf <- "../../../../Desktop/scratch/polyg_scratch/EVD_H5/RA-CAD_wtcc_hapmap.h5"
uhf <- "../../../../Desktop/scratch/polyg_scratch/RSSp_genome_gwas_uh/RA-CAD_wtcc_NoConfoundSmaller_sim.h5"
pvv <-0.875
quhf <- "../../../../Desktop/scratch/polyg_scratch/RSSp_genome_gwas_quh_chunk/RA-CAD_wtcc_NoConfoundSmaller_0.875.h5"
evdf <- normalizePath(snakemake@input[["evdf"]])
uhf <- normalizePath(snakemake@input[["uhf"]])
#pvv <-
quhf <- snakemake@output[["quhf"]]
stopifnot(!is.null(quhf))
quhf <- file.path(normalizePath(dirname(quhf)),basename(quhf))
pvv <- as.numeric(snakemake@params[["pvv"]])
if(length(pvv)==0){
  pvv <- 1.0
}

stopifnot(!is.null(pvv))

snp_df <- EigenH5::read_df_h5(evdf,"LDinfo")
tparam_df <- EigenH5::read_df_h5(uhf,"SimulationInfo")

EigenH5::write_df_h5(tparam_df,"SimulationInfo",quhf)
# ld_id <- as.character(sort(unique(snp_df$region_id)))
# stopifnot(length(ld_id)>0)

ld_grp <-get_objs_h5(evdf,"EVD")
D_df <- map_df(ld_grp,~data_frame(region_id=as.integer(.x),D=read_vector_h5(evdf,paste0("EVD/",.x),"D")))

sub_D_pvv <- filter_pvv(D_df,pvv)

Qf_df <- group_by(sub_D_pvv,region_id) %>%
  slice(n()) %>%
  ungroup() %>%
  select(region_id,col_chunksizes=chunk_id,row_chunksizes=chunksize) %>%
  mutate(filenames=evdf,groupnames=paste0("EVD/",region_id),datanames="Q",row_offsets=0,col_offsets=0)
stopifnot(nrow(Qf_df)==length(ld_grp))


uh_dff <- select(Qf_df,-region_id) %>% mutate(row_offsets=as.integer(c(0,head(cumsum(row_chunksizes),-1))),
                                              filenames=uhf,groupnames="/",
                                              datanames="uh",
                                              col_chunksizes=nrow(tparam_df))
stopifnot(slice(uh_dff,n()) %>% mutate(lrow=row_offsets+row_chunksizes) %>% pull(lrow)==nrow(snp_df) )


quh_dff <-select(Qf_df,-region_id) %>% mutate(row_offsets=as.integer(c(0,head(cumsum(col_chunksizes),-1))),
                                              row_chunksizes=col_chunksizes,
                                              filenames=quhf,groupnames="/",
                                              datanames="quh",
                                              col_chunksizes=nrow(tparam_df))
stopifnot(slice(quh_dff,n()) %>% mutate(lrow=row_offsets+row_chunksizes) %>% pull(lrow)==nrow(sub_D_pvv) )
stopifnot(is.numeric(sub_D_pvv$D))
EigenH5::write_vector_h5(quhf,"/","D",sub_D_pvv$D)
create_matrix_h5(quhf,"/","quh",numeric(),dims=c(nrow(sub_D_pvv),nrow(tparam_df)),chunksizes = c(1000,1))

#D <- map(paste0("EVD/",ld_id),~EigenH5::read_vector_h5(evdf,.x,"D")) %>% as_vector()

#
# Ql <- read_mat_l(Qf_df)
# uhl <- read_mat_l(uh_dff)
# quh_l <- purrr::map2(Ql,uhl,crossprod)
# quh_r <- do.call("rbind",quh_l)
# quh_mat <- read_matrix_h5(quhf,"/","quh")


EigenH5::write_df_h5(snp_df,"SNPinfo",quhf)
SeqSupport::crossprod_quh_h5(Qf_df,uh_dff,quh_dff)

#Rprof(NULL)
