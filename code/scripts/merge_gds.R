library(SeqArray)
library(SeqSupport)
library(purrr)
library(magrittr)

input_f <- snakemake@input[["gdsf"]]

MAF_cutoff <- as.numeric(snakemake@params[["MAF"]] %||% 0.0)
MAF_cutoff <- sort(abs(c(MAF_cutoff,1-MAF_cutoff)))
output_f <- snakemake@output[["gdsf"]]


## gds <- seqOpen(x)
## seqSetFilterCond(gds,maf=MAF_cutoff)
## tf <- tempfile()
## seqExport(gds,out.fn=tf)
## seqClose(gds)

seqMerge(input_f,out.fn=output_f,storage.option="LZ4_RA.fast")
