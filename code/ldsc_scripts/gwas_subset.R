#Rprof(filename=snakemake@output[["proff"]],append=F)
library(readr)
library(dplyr)
library(RcppEigenH5)

save.image()
inputf <- snakemake@input[["inputf"]]
LDchunkf <- snakemake@output[["outf"]]
LDchunk <- snakemake@params[["LDchunk"]]



input_df <- read_delim(inputf,delim="\t")
input_dfl<- split(input_df,input_df$region_id)

tparam_df <- data_frame(
    fgeneid = "1",
    tbias = NA_real_,
    tpve = NA_real_,
    tsigu = NA_real_,
    p= nrow(input_df),
    n=unique(input_df$N))

for(i in names(input_dfl)){
    tdfl <- input_dfl[[i]]
    stopifnot(sum(LDchunk==i)==1)
    j <- which(LDchunk==i)
    uhm <- matrix(tdfl$Z,length(tdfl$Z),1)
    snp_id <- tdfl$snp_id
    write_ivec_h5(LDchunkf[j],as.character(LDchunk)[j],"snp_id",snp_id)
    write_mat_h5(LDchunkf[j],as.character(LDchunk)[j],"uh",uhm,deflate_level=0,doTranspose=F)
    write_df_h5(tparam_df,"SimulationInfo",LDchunkf[j])
    }
#Rprof(NULL)
