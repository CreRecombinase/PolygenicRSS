

#library(profvis)
library(EigenH5)
library(LDshrink)
library(tidyverse)
# library(future)
library(progress)
# plan(sequential)
# mb <-profvis({


                                        # load("ssSNP.RData")
## save.image("evd.RData")
## stop()

#file.remove("evd.RData")
cutoff <- 1e-3
input_file <- snakemake@input[["input_file"]]
input_subset_file <- snakemake@input[["input_subset_file"]]
output_file <- snakemake@output[["evdf"]]
useLDshrink <- snakemake@params[["useLDshrink"]]=="T"
chunking <-snakemake@params[["useLDetect"]]
subldf <- snakemake@input[["subldf"]]
stopifnot(!is.null(chunking))
region_id <- snakemake@params[["region_id"]]

useLDetect <-chunking=="T"
useChunking <- !is.na(as.numeric(chunking))
if(useChunking){
    chunking <- as.integer(chunking)
}

cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]
m <- formals(LDshrink::LDshrink)[["m"]]
Ne <- formals(LDshrink::LDshrink)[["Ne"]]

dosage_dims <-EigenH5::dim_h5(input_file, "dosage")
tch <- EigenH5::dim_h5(input_file, "SNPinfo/chr")
SNPfirst <-  dosage_dims[1]==tch
if(!SNPfirst){
  stopifnot(dosage_dims[2]==tch)
}


stopifnot( !is.null(input_file),!is.null(input_subset_file), !is.null(output_file))

if(!is.null(subldf)){
    if(tools::file_ext(subldf)=="h5"){
        ind_v <- read_vector_h5(subldf,"SampleInfo/sample_id")
    }else{
        ind_v <- readRDS(subldf)
    }
    N <- length(ind_v)

}else{
    if(SNPfirst){
        N <-dosage_dims[2]
    }else{
        N <- dosage_dims[1]
    }
    ind_v <- 1:N
}


pl <- snakemake@wildcards
pl <- as_tibble(pl[names(pl)!=""])
write_df_h5(pl, filename = output_file, datapath = "Wildcards")


                                        #for(i in 1:num_b){

tdf <- read_df_h5(input_subset_file,paste0(region_id,"/SNPinfo"))
write_df_h5(tdf, output_file, "LDinfo")
if(SNPfirst){
    dosage <- EigenH5::read_matrix_h5(input_file,"dosage",subset_rows=tdf$ld_snp_id,subset_cols=ind_v,doTranspose=T)
}else{
    dosage <- EigenH5::read_matrix_h5(input_file,"dosage",subset_cols=tdf$ld_snp_id,subset_rows=ind_v)
}
mrid <- unique(tdf$region_id)
stopifnot(length(mrid)==1)
retl <- LDshrink_evd(dosage,
                     tdf$map,
                     m,
                     Ne,
                     cutoff,
                     useLDshrink=useLDshrink,
                     na.rm=T)
stopifnot(length(retl$D)==length(tdf$pos))
EigenH5::write_vector_h5(retl$D,output_file,"D")
EigenH5::write_matrix_h5(retl$Q,output_file,"Q")
EigenH5::write_vector_h5(retl$L2,output_file,"L2")

                                        #}
                                        #   },packages=c("EigenH5","LDshrink")))
                                        # },m=m,Ne=Ne,cutoff=cutoff,useLDshrink=useLDshrink,SNPfirst=SNPfirst,ind_l=ind_v)


# mmret <- SeqSupport::waitr(retl)
# mret <- purrr::map(retl,future::value)
# })
