#Rprof(filename=snakemake@output[["proff"]],append=F)
library(SeqSupport)
library(tidyverse)
# save.image()
# stop()


gdsf <- snakemake@input[["gdsf"]]
chrom <-as.integer(as.numeric(snakemake@params[["chrom"]]))
out_dir <-snakemake@params[["outdir"]]
#save.image()
chunkwise_LDshrink_ldsc(gdsf,chrom,out_dir)

#Rprof(NULL)
