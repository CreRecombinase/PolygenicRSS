library(SeqSupport)


hdf_f <- snakemake@input[["hdff"]]
chrom <- as.character(snakemake@params[["chrom"]])
cutoff <- as.numeric(as.character(snakemake@params[["cutoff"]]))


