# save.image()
# stop()
library(SeqArray)
library(SeqSupport)
library(readr)
library(dplyr)
library(RSSp)
library(purrr)

inputf <- snakemake@input[["traitf"]]
geno_gdsf <- snakemake@input[["geno_gdsf"]]
haplo_gdsf <- snakemake@input[["haplo_gdsf"]]




train_geno_f<- snakemake@output[["geno_train"]]
test_geno_f<- snakemake@output[["geno_test"]]


train_haplo_f<- snakemake@output[["haplo_train"]]
test_haplo_f<- snakemake@output[["haplo_test"]]


test_frac <-  as.numeric(snakemake@params[["test_frac"]])



gds <- seqOpen(geno_gdsf)
samples <- seqGetData(gds,"sample.id")
N <- length(samples)
test_num <- ceiling(N*test_frac)
test_id_geno <- sample(1:N,test_num,replace=F)
test_samples_geno <- samples[test_id_geno]
test_samples_haplo <- c(sapply(test_samples_geno,function(x)c(paste0(x,"-1"),paste0(x,"-2"))))

train_id_geno <- (1:N)[-test_id_geno]
train_samples_geno <- samples[train_id_geno]
train_samples_haplo <- c(sapply(train_samples_geno,function(x)c(paste0(x,"-1"),paste0(x,"-2"))))


write(test_samples_geno,test_geno_f,sep="\n")
write(train_samples_geno,train_geno_f,sep="\n")

write(test_samples_haplo,test_haplo_f,sep="\n")
write(train_samples_haplo,train_haplo_f,sep="\n")
