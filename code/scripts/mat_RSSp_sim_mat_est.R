# save.image()
# stop()

library(RSSp)
library(dplyr)
library(SeqSupport)
library(purrr)
library(readr)
rdsf <- snakemake@input[["rdsf"]]
outf <- snakemake@output[["dff"]]
proff <- snakemake@output[["proff"]]

LDchunk <- as.character(snakemake@params[["LDchunk"]])

stopifnot(file.exists(rdsf),
          !is.null(rdsf))

if(is.null(snakemake@params[["LDchunk"]])){
    quh_mat <- read_2d_mat_h5(rdsf,"/","quh")
    D <- read_vec(rdsf,"D")
    tparam_df <- read_df_h5(rdsf,"tparam_df")
}else{
    quh_mat <- read_2d_mat_h5(rdsf,LDchunk,"quh")
    D <- read_vec(rdsf,paste0("/",LDchunk,"/D"))
    tparam_df <- read_df_h5(rdsf,"SimulationInfo")
}


colnames(quh_mat) <- as.character(tparam_df$fgeneid)


n <-unique(tparam_df$n)
stopifnot(length(n)==1)

inputs <- map2_df(
    array_branch(quh_mat,margin=2),
    transpose(tparam_df),
    function(quhl,tpdfl,D){
        return(data_frame(fgeneid=tpdfl$fgeneid,
                          quh=quhl,
                          D=D,
                          p_n=tpdfl$p/tpdfl$n,
                          tbias=tpdfl$tbias,
                          tpve=tpdfl$tpve,
                          tsigu=tpdfl$tsigu))
        },D=D)

all_RSS_estimate <- function(data_df){
    tp_df <- distinct(data_df,fgeneid,.keep_all=T) %>% select(-D,-quh)
    return(cross_df(list(doConfound=c(T,F),log_params=c(F),useGradient=c(T))) %>%
           pmap_dfr(RSSp_estimate,data_df=data_df) %>% inner_join(tp_df,by="fgeneid"))
    }



rss_res <- group_by(inputs,fgeneid) %>% do(all_RSS_estimate(.)) %>% ungroup()

write_delim(rss_res,outf,delim="\t")
