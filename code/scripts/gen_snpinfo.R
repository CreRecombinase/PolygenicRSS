#Rprof(filename=snakemake@output[["proff"]],append=F)
library(SeqSupport)
library(tidyverse)
library(SeqArray)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
osnpf <- snakemake@output[["snpif"]]
snptxtf <- snakemake@output[["snptxtf"]]
gds <- seqOpen(gdsf)
snpi <- read_SNPinfo_gds(gds)
saveRDS(snpi,osnpf)
write_delim(snpi,snptxtf,delim="\t")
#Rprof(NULL)
