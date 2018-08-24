
library(RSSp)
library(dplyr)

                                        #library(SeqSupport)
library(EigenH5)
library(purrr)
library(readr)


#
# # save.image("trsp.RData")
# stop()
# setwd("~/Dropbox/PolygenicRSS/code/snakemake_files")
# load("trsp.RData")



rdsf <- snakemake@input[["rdsf"]]
outf <- snakemake@output[["dff"]]
doConfound <- snakemake@params[["doConfound"]]
pvv <- snakemake@params[["pvv"]]

y_grp <- snakemake@params[["y_grp"]]

if(is.null(y_grp)){
    y_grp  <- "SimulationInfo"
}

if(is.null(pvv)){
    pvv <- 0
}else{
    pvv <- as.numeric(pvv)
}

stopifnot(!is.na(pvv),
          !is.na(as.numeric(pvv)),
          length(pvv)==1)
if(is.null(doConfound)){
    doConfound <- F
}else{
    doConfound <- doConfound =="T"
}






stopifnot(file.exists(rdsf),
          !is.null(rdsf))
pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""]) %>% mutate(ttca=NA)
tparam_df <- read_df_h5(rdsf,y_grp) %>% mutate(ttca=NA) %>% inner_join(pl) %>% select(-ttca)


exec_fn <-  function(quhl,tpdfl,D,region_id){
    max_sigu <- calc_sigu(1,tpdfl$p/tpdfl$n)[1]
  data_frame(fgeneid=tpdfl$fgeneid,
             region_id=region_id,
             quh=quhl,
             D=D,
             p_n=tpdfl$p/tpdfl$n,
             tbias=tpdfl$tbias,
             tpve=tpdfl$tpve,
             tsigu=tpdfl$tsigu) %>%filter_normD(normD=pvv) %>%
    RSSp_estimate(pve_bounds=c(.Machine$double.eps, 4 - .Machine$double.eps),doConfound=doConfound,eigenvalue_cutoff=0) %>%
    inner_join(as_data_frame(tpdfl),by="fgeneid")
}




rss_res <- map2_df(
  array_branch(read_matrix_h5(rdsf,"quh"),margin=2),
  transpose(tparam_df),exec_fn,D=read_vector_h5(rdsf,datapath="/D"),
  region_id=read_vector_h5(rdsf,"SNPinfo/region_id"))


write_delim(rss_res,outf,delim="\t")


cat("Done!")
