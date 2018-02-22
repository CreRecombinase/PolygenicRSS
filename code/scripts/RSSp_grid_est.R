library(RSSp)
library(dplyr)
                                        #library(SeqSupport)
library(EigenH5)
library(purrr)
library(readr)

fix_path <- function(x){
    dirnames <- normalizePath(dirname(x))
    files <- basename(x)
    return(file.path(dirnames,files))
}
rdsf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/RSSp_genome_gwas_quh_chunk/sscombined_wtcc_NoConfoundSmall_sim.h5"

rdsf <- fix_path(snakemake@input[["rdsf"]])
outf <- fix_path(snakemake@output[["dff"]])
#proff <- fix_path(snakemake@output[["proff"]])

stopifnot(file.exists(rdsf),
          !is.null(rdsf))

tparam_df <- read_df_h5(rdsf,"SimulationInfo")

exec_fn <-  function(quhl,tpdfl,D){
  data_frame(fgeneid=tpdfl$fgeneid,
             quh=quhl,
             D=D,
             p_n=tpdfl$p/tpdfl$n,
             tbias=tpdfl$tbias,
             tpve=tpdfl$tpve,
             tsigu=tpdfl$tsigu) %>%
    RSSp_estimate_grid(bias_bounds=c(0,0)) %>%
    inner_join(as_data_frame(tpdfl),by="fgeneid")
}


#n_grp <- length(ldgrp)
#RcppEigenH5::concat_mat_chunks(rep(rdsf,n_grp),

rss_res <- map2_df(
  array_branch(read_mat_h5(rdsf,"/","quh"),margin=2),
  transpose(tparam_df),exec_fn,D=read_vector_h5(rdsf,"/","D"))


# group_by(rss_res,fgeneid) %>% mutate(slnZ=scale(lnZ,center=T,scale=T)) %>% ungroup() %>% ggplot(aes(x=pve,y=slnZ))+geom_point()+facet_wrap(~tpve)


write_delim(rss_res,outf,delim="\t")
