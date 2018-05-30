library(tidyverse)
library(EigenH5)
inf <- snakemake@input[["gdsf"]]
outf_gwas <- snakemake@output[["gwaso"]]
outf_ld <- snakemake@output[["ldo"]]
doSplit <- snakemake@params[["doSplit"]]=="T"

inds  <- read_vector_h5(inf,"SampleInfo","sample_id")
N <- length(inds)
#stopifnot(all(inds==1:N))
if(doSplit){
    gwas_inds <- sort(sample(1:N,ceiling(N/2),replace = F))
    LD_inds  <- (1:N)[-gwas_inds]
}else{
    gwas_inds <- 1:N
    LD_inds  <- 1:N
}
saveRDS(gwas_inds,outf_gwas)
saveRDS(LD_inds,outf_ld)
