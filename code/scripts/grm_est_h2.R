library(SeqSupport)
library(dplyr)
library(purrr)
library(readr)

xf <- snakemake@input[["gdsf"]]
yf <-snakemake@input[["h5f"]]
of <- snakemake@output[["grmef"]]
stopifnot(!is.null(of))

X <- scale(t(EigenH5::read_matrix_h5(xf,"/","dosage")),center=T,scale=F)
ym <- scale(EigenH5::read_matrix_h5(yf,"trait","ymat"),center=T,scale=F)
snp_si <- read_df_h5(xf,"SNPinfo") 
si <- EigenH5::read_df_h5(yf,"SimulationInfo")
sx <- sum(apply(X,2,var))
p <- ncol(X)
stopifnot(nrow(snp_si)==p)
svdX <- svd(X/sqrt(p))

U <- svdX$u
D <- svdX$d
rm(svdX,X)
gc()

optim_pg <- function(par,D,U,y,sx,p){
    ret <- polygenic.model_svd(U = U,D=D,y = y,h = par[1],sx=sx,p=p)
    return(-ret$logw)
}

## optim_cg <- function(par,x,y){
##     ret <- polygenic.model(X=x,y=y,h=par[1])
##     return(-ret$logw)
## }
    


## tyml <- yml[1:5]
## p_res <- tyml %>% imap_dfr(function(yv,fgeneid,X){
##     optimise(optim_cg,interval=c(1e-4,1),x=X,y=yv) %>%
##         as_data_frame() %>% mutate(fgeneid=as.character(fgeneid))
## },X=X) %>% inner_join(si)


all_res_df <- array_branch(ym,margin=2) %>% imap_dfr(function(yv,fgeneid,U,D,sx,p){
    optimise(optim_pg,interval = c(1e-4,1),U=U,D=D,sx=sx,p=p,y=yv) %>%
        as_data_frame() %>% mutate(fgeneid=as.character(fgeneid))
},U=U,D=D,sx=sx,p=p) %>%
    inner_join(si)



write_delim(all_res_df,path = of,delim="\t")







