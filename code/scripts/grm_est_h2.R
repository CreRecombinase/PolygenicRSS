library(SeqSupport)
library(dplyr)
library(purrr)
library(readr)

xf <- snakemake@input[["gdsf"]]
yf <-snakemake@input[["h5f"]]
of <- snakemake@output[["gremef"]]


X <- scale(t(RcppEigenH5::read_2d_mat_h5(xf,"/","dosage")),center=T,scale=F)
ym <- scale(RcppEigenH5::read_2d_mat_h5(yf,"trait","ymat"),center=T,scale=F)
si <- RcppEigenH5::read_df_h5(yf,"SimulationInfo")


optim_pg <- function(par,X,y){
    ret <-polygenic.model(X = X,y = y,h = par[1])
    gc()
    return(-ret$logw)
}

all_res_df <- apply(ym,2,function(yv,x){
    optimise(optim_pg,interval = c(1e-4,1),X=X,y=ym)
},x=X) %>%
    map_df(as_data_frame) %>%
    mutate(fgeneid=as.character(1:n())) %>%
    inner_join(si)

write_delim(all_res_df,path = of,delim="\t")







