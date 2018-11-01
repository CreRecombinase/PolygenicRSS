library(EigenH5)
library(LDshrink)
library(tidyverse)
library(progress)
writefn <- saveRDS


input_file <- snakemake@input[["input_file"]]
snp_dff <- snakemake@input[["subsnpf"]]
useLDshrink <- T




chunk_tot <- as.integer(snakemake@params[["chunk_tot"]])
chunk_ind <- as.integer(snakemake@params[["chunk_ind"]])
output_file <- snakemake@output[["evdf"]]

useLDshrink <- snakemake@params[["useLDshrink"]]=="T"
useLDT <- snakemake@params[["useLDetect"]]=="T"
pop <- snakemake@params[["pop"]]
if(pop=="omni"){
    pop <- "CEU"
}
r2_cutoff <- snakemake@params[["r2c"]]
stopifnot(!is.null(r2_cutoff))
r2_cutoff <- as.numeric(r2_cutoff)
stopifnot(!is.na(r2_cutoff))

data(map_parameters)
my_map <- filter(map_parameters, Pop==pop)

stopifnot(nrow(my_map)==1)
stopifnot(!is.null(input_file), !is.null(output_file), !useLDT)

cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]

m <- my_map$m
Ne <- my_map$Ne



if(!useLDT){
    tri_block_fn <- function(tot_size, chunk_ind){
        mat_size <- (1/2)*(sqrt(8*tot_size+1)-1)
        stopifnot(mat_size==round(mat_size))
        tmat <- matrix(0, mat_size, mat_size)
        tmat[upper.tri(tmat, diag=T)] <- 1:tot_size
        stopifnot(all(tmat[upper.tri(tmat, diag=T)]>0))
        trow <- row(tmat)[tmat==chunk_ind]
        tcol <- col(tmat)[tmat==chunk_ind]
        retdf <- data_frame(row=trow, col=tcol, chunk_tot=mat_size, chunk_ind=chunk_ind)
        stopifnot(nrow(retdf)==1)
        return(retdf)
    }
    my_chunk <- tri_block_fn(tot_size = chunk_tot, chunk_ind = chunk_ind)
    chunk_a <- my_chunk$row
    chunk_b <- my_chunk$col
}else{
    chunk_a <- chunk_ind
    chunk_b <- chunk_ind
}



snp_dfa <- read_df_h5(snp_dff,paste0(chunk_a,"/SNPinfo"))

dosage_dims <-EigenH5::dim_h5(input_file, "dosage")
tch <- EigenH5::dim_h5(input_file, "SNPinfo/chr")
SNPfirst <-  dosage_dims[1]==tch
if(!SNPfirst){
  stopifnot(dosage_dims[2]==tch)
}

if(SNPfirst){
    dosage_a <- scale(EigenH5::read_matrix_h5(input_file,"dosage",subset_rows=snp_dfa$snp_id,doTranspose=T),center=T,scale=F)
}else{
    dosage <- scale(EigenH5::read_matrix_h5(input_file,"dosage",subset_cols=snp_dfa$snp_id),center=T,scale=F)
}

if(chunk_a==chunk_b){
    writefn(ld2df(dosage_a,
                  snp_dfa$map,
                  snp_dfa$SNP,
                  m,
                  Ne,
                  cutoff,
                  r2_cutoff,
                  progress=T,
                  useLDshrink = useLDshrink),  output_file)
}else{
    snp_dfb <- read_df_h5(snp_dff,paste0(chunk_b,"/SNPinfo"))
    if(SNPfirst){
        dosage_b <- EigenH5::read_matrix_h5(input_file,"dosage",subset_rows=snp_dfb$snp_id,doTranspose=T)
    }else{
        dosage_b <- EigenH5::read_matrix_h5(input_file,"dosage",subset_cols=snp_dfb$snp_id)
    }
    writefn(ld2df_p(data_a = dosage_a,
                    data_b = dosage_b,
                    mapd_a = snp_dfa$map,
                    mapd_b = snp_dfb$map,
                    rsid_a = snp_dfa$SNP,
                    rsid_b = snp_dfb$SNP,
                    m=m,
                    Ne=Ne,
                    cutoff = cutoff,
                    r2cutoff = r2_cutoff,
                    progress = F,
                    useLDshrink = T), output_file)
}
