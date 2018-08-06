library(SeqArray)
seqGDS2SNP(snakemake@input[["input_file"]],snakemake@output[["output_file"]],compress.geno="LZ4_RA.fast")
