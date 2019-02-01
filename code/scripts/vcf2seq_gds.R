library(SeqArray)



inpf <- snakemake@input[["vcf"]]
toutf <- snakemake@output[["gds"]]

seqVCF2GDS(vcf.fn = inpf,
           out.fn = toutf,
           storage.option = snakemake@params[["geno_compression"]])

cat("Success!")
