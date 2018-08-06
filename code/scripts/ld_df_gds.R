library(EigenH5)
library(LDshrink)
library(tidyverse)
library(SeqArray)
library(progress)


writefn <- saveRDS


input_file <- snakemake@input[["input_file"]]

useLDshrink <- T


snp_dfa <- readRDS(snakemake@input[["snpdff_a"]])
snp_dfb <- readRDS(snakemake@input[["snpdff_a"]])

subldf <- snakemake@input[["subldf"]]

output_file <- snakemake@output[["evdf"]]

useLDshrink <- snakemake@params[["useLDshrink"]]=="T"
useLDT <- snakemake@params[["useLDetect"]]=="T"
pop <- snakemake@params[["pop"]]
r2_cutoff <- snakemake@params[["r2c"]]
stopifnot(!is.null(r2_cutoff))
r2_cutoff <- as.numeric(r2_cutoff)
stopifnot(!is.na(r2_cutoff))
chunk_a <- as.integer(snakemake@params[["chunk_a"]])
chunk_b <- as.integer(snakemake@params[["chunk_b"]])
chunk_tot <- as.integer(snakemake@params[["chunk_tot"]])

data(map_parameters)
my_map <- filter(map_parameters, Pop==pop)

stopifnot(nrow(my_map)==1)

cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]

m <- my_map$m
Ne <- my_map$Ne



stopifnot(!is.null(input_file), !is.null(output_file), !useLDT)

stopifnot(file.exists(input_file), !file.exists(output_file))



ind_v <- readRDS(subldf)
## num_b <- length(snp_dfl)
tgds_a <- seqOpen(input_file,allow.duplicate=TRUE)

seqSetFilter(tgds_a,variant.id=snp_dfa$ld_snp_id,sample.sel=ind_v)
dosage_a <- scale(seqGetData(tgds_a,"$dosage"),center=T,scale=F)
if(chunk_a==chunk_b){
    seqClose(tgds_a)
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
    seqSetFilter(tgds_a, variant.id = snp_dfb$ld_snp_id, sample.sel=ind_v)
    dosage_b <- scale(seqGetData(tgds_a,"$dosage"), center=T, scale=F)
    seqClose(tgds_a)
    writefn(ld2df_p(scaled_data_a = dosage_a,
                    scaled_data_b = dosage_b,
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
      # write_fst(df_r, output_file, 100)

   # pb$tick()
}
