library(SNPRelate)
library(tidyverse)
inpf <- snakemake@input[["vcff"]]
toutf <- snakemake@output[["temp_gds"]]
#cores <- as.integer(snakemake@threads)
#seqParallelSetup(cores)
snpgdsVCF2GDS(inpf, toutf, method="biallelic.only",compress.geno="LZ4_RA.fast")
## seqVCF2GDS(vcf.fn = inpf, out.fn = toutf,
##            storage.option = )
