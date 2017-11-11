# save.image()
# stop()
library(SeqArray)
library(SeqSupport)
library(readr)
library(dplyr)
library(RSSp)
library(purrr)

inputf <- snakemake@input[["traitf"]]
igdsf <- snakemake@input[["gdsf"]]
samplef <- snakemake@input[["tt_traitf"]]



gdsf <- snakemake@output[["gdsf"]]
hdf5f <- snakemake@output[["hdf5f"]]


samples <- scan(samplef,what=character(),sep="\n")


gds <- seqOpen(igdsf)
dsamples <- seqGetData(gds,"sample.id")
stopifnot(sum(samples %in% dsamples)>0)

output_df <- subset_export_gds(gds,samples,NULL,gdsf,hdf5f)

