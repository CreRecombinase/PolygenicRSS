#Rprof(filename=snakemake@output[["proff"]],append=F)
library(RSSp)
library(SeqSupport)
library(dplyr)
library(purrr)
library(readr)
rss_rdsf <- snakemake@input[["rdsf"]]
outf <- snakemake@output[["dff"]]

quh_mat <- read_2d_mat_h5(rss_rdsf,"/","quh")
D <- read_vec(rss_rdsf,"D")
tparam_df <- transpose(read_df_h5(rss_rdsf,"tparam_df"))

colnames(quh_mat) <- as.character(tparam_df$fgeneid) 
rss_res <- map2_df(
    array_branch(quh_mat,margin=2),
    tparam_df,
    function(x,y,dvec){
        RSSp_oracle(
            quh=x,
            dvec=dvec,
            tsigu=y$tsigu,
            tbias=y$tbias,
            fgeneid=y$fgeneid
        )
    },dvec=D
)



write_delim(rss_res,outf,delim="\t")
