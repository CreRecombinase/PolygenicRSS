library(SeqArray)

SeqArray::seqBED2GDS(bed.fn = snakemake@input[["bedf"]],
                     bim.fn = snakemake@input[["bimf"]],
                     fam.fn = snakemake@input[["famf"]],
                     out.gdsfn = snakemake@output[["gds"]],
                     compress.geno = snakemake@params[["geno_compression"]],
                     compress.annotation = snakemake@params[["anno_compression"]]
                     )
cat("Success!")
