
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



stopifnot(file.exists(rdsf),
          !is.null(rdsf))

tparam_df <- read_df_h5(rdsf,"SimulationInfo")
ldgrp <- "/"

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




rss_res <- map2_df(
  array_branch(read_mat_h5(rdsf,ldgrp,"quh"),margin=2),
  transpose(tparam_df),exec_fn,D=read_vector_h5(rdsf,ldgrp,"D"))


write_delim(rss_res,outf,delim="\t")
