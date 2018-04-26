
library(RSSp)
library(dplyr)

                                        #library(SeqSupport)
library(EigenH5)
library(purrr)
library(readr)



# save.image()
# stop()
fix_path <- function(x){
    dirnames <- normalizePath(dirname(x))
    files <- basename(x)
    return(file.path(dirnames,files))
}


rdsf <- snakemake@input[["rdsf"]]
outf <- snakemake@output[["dff"]]
doConfound <- snakemake@params[["doConfound"]]
pvv <- snakemake@params[["pvv"]]
doMoM <- snakemake@params[["doMoM"]]

if(is.null(pvv)){
    pvv <- 1
}else{
    pvv <- as.numeric(pvv)
}
if(is.null(doConfound)){
    doConfound <- F
}else{
    doConfound <- doConfound =="T"
}
if(is.null(doMoM)){
  doMoM <- F
}else{
  doMoM <- doMoM =="T"
}





stopifnot(file.exists(rdsf),
          !is.null(rdsf))
pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""]) %>% mutate(ttca=NA)
tparam_df <- read_df_h5(rdsf,"SimulationInfo") %>% mutate(ttca=NA) %>% inner_join(pl) %>% select(-ttca)
ldgrp <- "/"
if(doMoM){
  exec_fn <-  function(quhl,tpdfl,D,region_id){
      data_frame(fgeneid=tpdfl$fgeneid,
                 region_id=region_id,
                 quh=quhl,
                 D=D,
                 p_n=tpdfl$p/tpdfl$n,
                 tbias=tpdfl$tbias,
                 tpve=tpdfl$tpve,
                 tsigu=tpdfl$tsigu) %>% filter_pvv(pvv) %>%
          RSSp_mom(doConfound=doConfound) %>%
          inner_join(as_data_frame(tpdfl),by="fgeneid")
  }

}else{
  exec_fn <-  function(quhl,tpdfl,D,region_id){
      data_frame(fgeneid=tpdfl$fgeneid,
                 region_id=region_id,
                 quh=quhl,
                 D=D,
                 p_n=tpdfl$p/tpdfl$n,
                 tbias=tpdfl$tbias,
                 tpve=tpdfl$tpve,
                 tsigu=tpdfl$tsigu) %>% filter_pvv(pvv) %>%
        RSSp_estimate(doConfound=doConfound) %>%
        inner_join(as_data_frame(tpdfl),by="fgeneid")
  }
}



rss_res <- map2_df(
  array_branch(read_matrix_h5(rdsf,ldgrp,"quh"),margin=2),
  transpose(tparam_df),exec_fn,D=read_vector_h5(rdsf,ldgrp,"D"),
  region_id=read_vector_h5(rdsf,"SNPinfo","region_id"))


write_delim(rss_res,outf,delim="\t")
