library(SeqArray)
library(SeqSupport)
library(readr)
library(dplyr)
library(RSSp)
library(purrr)

inputf <- snakemake@input[["inputf"]]
haplo_gdsf <- snakemake@input[["haplo_gdsf"]]
ohaplo_gdsf <- snakemake@output[["haplo_gdsf"]]
ohaplo_hdf5 <- snakemake@output[["haplo_hdf5"]]


stopifnot(file.exists(c(haplo_gdsf)),
          !is.null(haplo_gdsf),
          !is.null(ohaplo_hdf5),
          !is.null(ohaplo_gdsf),
          !is.null(inputf))

input_df <- read_delim(inputf,delim="\t") %>% select(-one_of(c("chr",
                                                               "snp_id",
                                                               "mat_snp_id","pos")))

gds <- seqOpen(haplo_gdsf)
samples <- seqGetData(gds,"sample.id")
output_df <- subset_export_gds(gds,samples,input_df,ohaplo_gdsf,ohaplo_hdf5)
stopifnot(nrow(output_df)==nrow(input_df))
