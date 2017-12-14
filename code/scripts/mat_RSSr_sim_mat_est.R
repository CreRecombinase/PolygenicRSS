#setwd("/home/nwknoblauch/Dropbox/snakemake_files")
## load(".RData")
library(RSSp)
library(dplyr)
library(SeqSupport)
library(rssr)
library(purrr)
library(readr)

uhf <- snakemake@input[["uhf"]]
ldf <- snakemake@input[["evdf"]]
outf <- snakemake@output[["dff"]]

LDchunk <- as.character(snakemake@params[["LDchunk"]])



tparam_df <- read_df_h5(uhf,"SimulationInfo")
ldgrp <- ifelse(is.null(snakemake@params[["LDchunk"]]),"",snakemake@params[["LDchunk"]])
bh <- read_2d_mat_h5(uhf,ldgrp,"uh")
se <- read_2d_mat_h5(uhf,ldgrp,"se")
bh <- bh*se
colnames(bh) <- tparam_df$fgeneid
ns <- ncol(bh)
R <- read_2d_mat_h5(ldf,"LD","R")

rssr_res <- map2(array_branch(bh,2),array_branch(se,2),RSSr_estimate,R=R,threads=1) %>% map2_df(as.character(1:ns),~mutate(.x,fgeneid=.y))%>% 
    inner_join(tparam_df,by="fgeneid") %>% mutate(region_id=LDchunk)
## rssr_res <- RSSr_estimate_mat(R=R,
##                               betahat=bh,
##                               se=se,grid_total=100) %>% mutate(fgeneid=as.character(fgeneid)) %>% 


    
write_delim(rssr_res,outf,delim="\t")
