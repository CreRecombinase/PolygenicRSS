
library(SeqArray)
library(SeqSupport)
library(readr)
library(dplyr)
library(RSSp)
library(purrr)

inputf <- snakemake@input[["traitf"]]
geno_gdsf <- snakemake@input[["geno_gdsf"]]
haplo_gdsf <- snakemake@input[["haplo_gdsf"]]




train_geno_gdsf<- snakemake@output[["geno_gdsf_train"]]
test_geno_gdsf<- snakemake@output[["geno_gdsf_test"]]

train_geno_hdf5<- snakemake@output[["geno_hdf5_train"]]
test_geno_hdf5<- snakemake@output[["geno_hdf5_test"]]


train_haplo_gdsf<- snakemake@output[["haplo_gdsf_train"]]
test_haplo_gdsf<- snakemake@output[["haplo_gdsf_test"]]

train_haplo_hdf5<- snakemake@output[["haplo_hdf5_train"]]
test_haplo_hdf5<- snakemake@output[["haplo_hdf5_test"]]


osnpf <- snakemake@output[["snpif"]]
snptxtf <- snakemake@output[["snptxtf"]]

test_frac <-  as.numeric(snakemake@params[["test_frac"]])



input_df <- read_delim(inputf,delim="\t") %>% select(-one_of(c("chr","snp_id")))



gds <- seqOpen(geno_gdsf)
samples <- seqGetData(gds,"sample.id")
N <- length(samples)
test_num <- ceiling(N*test_frac)
test_id <- sample(1:N,test_num,replace=F)
train_id <- (1:N)[-test_id]

test_sample_id <- samples[test_id]
train_sample_id <- samples[train_id]
cat("test_num_geno:",test_num,"\n")
subset_export_gds(gds,test_sample_id,input_df,test_geno_gdsf,test_geno_hdf5)
gds <- seqOpen(geno_gdsf)
subset_export_gds(gds,train_sample_id,input_df,train_geno_gdsf,train_geno_hdf5)




gds <- seqOpen(haplo_gdsf)
samples <- seqGetData(gds,"sample.id")
N <- length(samples)
test_num <- ceiling(N*test_frac)
cat("test_num_haplo:",test_num,"\n")
test_id <- sample(1:N,test_num,replace=F)
train_id <- (1:N)[-test_id]

test_sample_id <- samples[test_id]
train_sample_id <- samples[train_id]
subset_export_gds(gds,test_sample_id,input_df,test_haplo_gdsf,test_haplo_hdf5)
gds <- seqOpen(haplo_gdsf)
output_df <- subset_export_gds(gds,train_sample_id,input_df,train_haplo_gdsf,train_haplo_hdf5)



saveRDS(output_df,osnpf)
write_delim(output_df,snptxtf,delim="\t")
