library(SeqSupport)

evdf <- snakemake@input[["evdf"]]

quh_outf <- snakemake@output[["quh_hf"]]
uh_outf <- snakemake@output[["uh_hf"]]

LDchunks <- snakemake@params[["LDchunk"]]

mfgeneid <- as.character(snakemake@params[["fgeneid"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))



stopifnot(!is.null(LDchunks))
LDchunks <- as.character(LDchunks)




