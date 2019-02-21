library(SeqArray)
library(SeqSupport)
library(purrr)
library(EigenH5)
library(magrittr)

input_f <- snakemake@input[["gdsf"]]
samp_f <- snakemake@input[["hf"]]
MAF_cutoff <- as.numeric(snakemake@params[["MAF"]] %||% 0.0)
MAF_cutoff <- sort(abs(c(MAF_cutoff,1-MAF_cutoff)))
output_f <- snakemake@output[["gdsf"]]

samp_df <- read_df_h5(samp_f,"SampleInfo")
stopifnot(!is.null(output_f),!file.exists(output_f))
stopifnot(all(file.exists(input_f)),!is.null(input_f))

gds <- seqOpen(input_f)
seqSetFilterCond(gds,maf=MAF_cutoff)
seqSetFilter(gds,sample.id=samp_df$sample_id,action="intersect")
stopifnot(calc_N(gds)==nrow(samp_df))
seqExport(gds,out.fn=output_f)
seqClose(gds)

## seqMerge(seq_files,out.fn=output_f,storage.option="LZ4_RA.fast")
## unlink(seq_files)
