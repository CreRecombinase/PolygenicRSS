library(SeqArray)
library(tidyverse)


## inpf <- "/run/media/nwknoblauch/Data/wtcc_input/vcf/t1d_19_dup.vcf.gz"
## toutf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/gds/scombined_19.gds"
inpf <- snakemake@input[["vcff"]]
toutf <- snakemake@output[["temp_gds"]]

cores <- as.integer(snakemake@threads)

seqParallelSetup(cores)
seqVCF2GDS(vcf.fn = inpf, out.fn = toutf,
           storage.option = "LZ4_RA.fast")
## tgds <- seqOpen(toutf)

## tm <- seqMissing(tgds,T,T)

## seqApply(tgds, "$dosage", duplicated, margin="by.variant", parallel=cores,.progress=TRUE,MARGIN=2)
## tx <- seqBlockApply(gdsfile = tgds,var.name = "$dosage",FUN = duplicated.array,margin = c("by.variant"),as.is="unlist",.progress=T,MARGIN=2)
