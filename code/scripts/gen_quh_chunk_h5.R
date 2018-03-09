                                        #Rprof(filename=snakemake@output[["proff"]],append=F)

#library(SeqSupport)
library(tidyverse)
library(RSSp)
library(LDshrink)
library(EigenH5)


evdf <- normalizePath(snakemake@input[["evdf"]])
uhf <- normalizePath(snakemake@input[["uhf"]])
svd <- snakemake@params[["svd"]]
if(is.null(svd)){
    svd <- "F"
}
svd <- svd =="T"
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

ld_gpn <- ifelse(svd,"SVD","EVD")
ld_dn <- ifelse(svd,"d","D")
ld_qn <- ifelse(svd,"V","Q")
ld_grp <-get_objs_h5(evdf,ld_gpn)
D_df <- map_df(ld_grp,~data_frame(region_id=as.integer(.x),D=read_vector_h5(evdf,paste0(ld_gpn,"/",.x),ld_dn))) %>% arrange(as.integer(region_id))
if(svd){
    D_df <- mutate(D_df,D=D^2)
}
## ld_grp <-get_objs_h5(evdf,"EVD")
## D_df <- map_df(ld_grp,~data_frame(region_id=as.integer(.x),D=read_vector_h5(evdf,paste0("EVD/",.x),"D"))) %>% arrange(as.integer(region_id))
nD <- mutate(D_df,row_id=1:n()) %>% group_by(region_id) %>% mutate(sD=sort(D,decreasing=T))%>% filter(D!=sD) %>% ungroup()

stopifnot(nrow(snp_df)==nrow(D_df))
sub_D_pvv <- filter_pvv(D_df,pvv)

Qf_df <- group_by(sub_D_pvv,region_id) %>%
  slice(n()) %>%
  ungroup() %>%
  select(region_id,col_chunksizes=chunk_id,row_chunksizes=chunksize) %>%
  mutate(filenames=evdf,groupnames=paste0(ld_gpn,"/",region_id),datanames=ld_qn,row_offsets=0,col_offsets=0) %>% arrange(as.integer(region_id))
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
create_matrix_h5(quhf,"/","quh",numeric(),dims=c(nrow(sub_D_pvv),nrow(tparam_df)),chunksizes = c(pmin(1024,nrow(sub_D_pvv)),1))




EigenH5::write_df_h5(snp_df,"SNPinfo",quhf)
SeqSupport::crossprod_quh_h5(Qf_df,uh_dff,quh_dff)

#Rprof(NULL)
