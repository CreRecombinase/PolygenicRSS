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

ogeno_hdf5 <- snakemake@output[["geno_hdf5"]]
ohaplo_hdf5 <- snakemake@output[["haplo_hdf5"]]

ogeno_gdsf <- snakemake@output[["geno_gdsf"]]
ohaplo_gdsf <- snakemake@output[["haplo_gdsf"]]




osnpf <- snakemake@output[["snpif"]]
snptxtf <- snakemake@output[["snptxtf"]]


input_df <- read_delim(inputf,delim="\t") %>% select(-one_of(c("chr","snp_id")))



gds <- seqOpen(geno_gdsf)
samples <- seqGetData(gds,"sample.id")
tr <- subset_export_gds(gds,samples,input_df,ogeno_gdsf,ogeno_hdf5)


gds <- seqOpen(haplo_gdsf)
samples <- seqGetData(gds,"sample.id")
output_df <- subset_export_gds(gds,samples,input_df,ohaplo_gdsf,ohaplo_hdf5)



saveRDS(output_df,osnpf)
write_delim(output_df,snptxtf,delim="\t")
