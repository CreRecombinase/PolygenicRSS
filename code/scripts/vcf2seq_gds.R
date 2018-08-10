library(SeqArray)
library(tidyverse)


inpf <- snakemake@input[["vcff"]]
toutf <- snakemake@output[["temp_gds"]]

cores <- as.integer(snakemake@threads)
seqParallelSetup(cores)
seqVCF2GDS(vcf.fn = inpf, out.fn = toutf,digest=FALSE,
           storage.option = "LZ4_RA.fast")
