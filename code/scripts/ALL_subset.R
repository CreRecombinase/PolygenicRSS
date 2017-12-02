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






snptxtf <- snakemake@output[["traitf"]]



gds <- seqOpen(geno_gdsf)
N <- calc_N(gds)

snpi <- read_SNPinfo_ldsc(gds) %>% mutate(N=N,Z=rnorm(n()))
seqClose(gds)


write_delim(snpi,snptxtf,delim="\t")
