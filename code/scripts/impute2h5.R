library(tidyverse)
library(readr)
library(dplyr)
library(SeqArray)
library(SeqSupport)
library(SNPRelate)

hap.fn <- snakemake@input[["hapf"]]
leg.fn <- snakemake@input[["legf"]]

chrom <- snakemake@params[["chrom"]]
out.hdf5 <- snakemake@output[["hdff"]]

snp_df <- read_delim(leg.fn,delim=" ",col_names=c("SNP","pos","a0","a1")) %>% unite(allele,a0,a1,sep=",")

p <- nrow(snp_df)
snp_df <- mutate(snp_df,SNP=as.integer(gsub("rs([0-9]+)","\\1",SNP)))
N <- scan(hap.fn,what = integer(),nlines = 1,sep=" ")

write_df_h5(snp_df,out.hdf5,"SNPinfo")
gz2hdf5(hap.fn,out.hdf5,"/","dosage",p,N,10000)
                           


