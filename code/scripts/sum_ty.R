# save.image("sum.RData")
# stop()
#load("sum.RData")
library(SeqSupport)
library(EigenH5)
library(tidyverse)

input_f <- snakemake@input[["input_f"]]
covarf <- snakemake@input[["covarf"]]
outputf  <- snakemake@output[["h5f"]]
n_covars <- as.integer(snakemake@params[["nc"]] %||% 0L)
if(is.na(n_covars)){
    cat(n_covars,"\n")
    stop("nope!")
}

pl <- snakemake@wildcards
pl <- as_tibble(pl[names(pl)!=""])
EigenH5::write_df_h5(pl,filename=outputf, "Wildcards")

input_f_a <- input_f[1]

tp_df <- read_df_h5(input_f_a,"SimulationInfo") %>% mutate(trait_id=1:n())
write_df_h5(tp_df,outputf,"SimulationInfo")
ymat <- read_matrix_h5(input_f_a,"genetic_trait/genetic_ymat")
if(length(input_f)>1){

    rest_f <- input_f[-1]

    for( tif in rest_f){
        tymat <- read_matrix_h5(tif,"genetic_trait/genetic_ymat")
        stopifnot(all(dim(tymat)==dim(ymat)))
        ymat <- ymat+tymat
    }

}

stopifnot(!is.null(covarf))
if(!is.null(covarf) ){
    Q <- read_matrix_h5(covarf,"covariates")[,seq_len(n_covars)]
    stopifnot(nrow(Q)==ncol(ymat))

    real_ymat <- gen_sim_resid(t(ymat),tp_df,Q=Q)
}else{
    real_ymat <- gen_sim_resid(t(ymat),tp_df)
}



write_matrix_h5(real_ymat,outputf,"ymat/trait")
