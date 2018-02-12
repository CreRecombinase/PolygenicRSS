

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

rdsf <- fix_path(snakemake@input[["rdsf"]])
outf <- fix_path(snakemake@output[["dff"]])
#proff <- fix_path(snakemake@output[["proff"]])

LDchunk <- as.character(snakemake@params[["LDchunk"]])

stopifnot(file.exists(rdsf),
          !is.null(rdsf))

tparam_df <- read_df_h5(rdsf,"SimulationInfo")
ldgrp <- ifelse(is.null(snakemake@params[["LDchunk"]]),"/",snakemake@params[["LDchunk"]])

exec_fn <-  function(quhl,tpdfl,D){
  data_frame(fgeneid=tpdfl$fgeneid,
             quh=quhl,
             D=D,
             p_n=tpdfl$p/tpdfl$n,
             tbias=tpdfl$tbias,
             tpve=tpdfl$tpve,
             tsigu=tpdfl$tsigu) %>%
    RSSp_estimate(doConfound=F) %>%
    inner_join(as_data_frame(tpdfl),by="fgeneid")
}


#n_grp <- length(ldgrp)
#RcppEigenH5::concat_mat_chunks(rep(rdsf,n_grp),

rss_res <- map2_df(
  array_branch(read_mat_h5(rdsf,ldgrp,"quh"),margin=2),
  transpose(tparam_df),exec_fn,D=read_vector_h5(rdsf,ldgrp,"D"))
# colnames(quh_mat) <- as.character(tparam_df$fgeneid)

# quh_mat <- quh_mat[,1,drop=F]
# tparam_df <- slice(tparam_df,1)

# n <-unique(tparam_df$n)
# stopifnot(length(n)==1)



# all_RSS_estimate <- function(data_df){
#     tp_df <- distinct(data_df,fgeneid,.keep_all=T) %>% select(-D,-quh)
#     return(cross_df(list(doConfound=c(T,F),log_params=c(F),useGradient=c(T))) %>%
#            pmap_dfr(RSSp_estimate,data_df=data_df) %>% inner_join(tp_df,by="fgeneid"))
#     }



# rss_res <- group_by(inputs,fgeneid) %>% do(all_RSS_estimate(.)) %>% ungroup()

write_delim(rss_res,outf,delim="\t")
