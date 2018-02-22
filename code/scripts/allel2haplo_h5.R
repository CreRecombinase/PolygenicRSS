library(EigenH5)
library(tidyverse)
library(SeqSupport)
library(progress)
library(LDshrink)
# inpf <- dir("/home/nwknoblauch/Desktop/scratch/polyg_scratch/impute",full.names=T)
#
# outf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/h5/EUR1_seq_wtcc_geno.h5"
# afc <- 0.01
tinpf <- snakemake@input[["legf"]]
outf <- snakemake@output[["h5f"]]

chunksize <- 20000


# tinpf <- inpf[1]
snp_df <- read_SNPinfo_allel(tinpf)
p<- nrow(snp_df)
good_snp_df <-filter(snp_df,numalt==1,VT=="SNP") %>% arrange(chr,pos) %>% select(-svlen) %>% mutate(nsnp_id=1:n())
np <- nrow(good_snp_df)
read_chunk <-gl(n=ceiling(np/chunksize),k=chunksize,length = np)
trl <-split(good_snp_df,read_chunk)
d_dims <- EigenH5::get_dims_h5(tinpf,"/","dosage")
N <- d_dims[2]
EigenH5::create_matrix_h5(filename = outf,groupname = "/",dataname = "dosage",data = integer(),doTranspose = F,dims = c(N,np),chunksizes = c(N,1000))
pb <- progress::progress_bar$new(total = length(trl))
for(x in trl){
  tD <- t(EigenH5::read_matrix_h5(filename = tinpf,groupname = "/",dataname = "dosage",subset_rows = x$snp_id))
  EigenH5::write_matrix_h5(filename = outf,groupname = "/",dataname = "dosage",data = tD,offsets = c(0,x$nsnp_id[1]-1))
  pb$tick()
}
good_snp_df <- mutate(good_snp_df,snp_id=nsnp_id) %>% select(-nsnp_id)
write_df_h5(df = good_snp_df,outfile = outf,groupname = "SNPinfo")




