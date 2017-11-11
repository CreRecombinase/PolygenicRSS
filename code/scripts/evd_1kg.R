#Rprof(filename=snakemake@output[["proff"]],append=F)


gdsf <- snakemake@input[["gdsf"]]
outf <- snakemake@output[["evdf"]]
region_id <- as.integer(as.numeric(snakemake@params[["region_id"]]))
cutoff <- as.numeric(snakemake@params[["cutoff"]])
if(length(cutoff)==0){
    cutoff <- formals(SeqSupport::chunkwise_LDshrink)[["cutoff"]]
}
library(rhdf5)


stopifnot(!is.null(region_id), !is.null(gdsf), !is.null(outf))

stopifnot(file.exists(gdsf), !file.exists(outf))

SeqSupport::chunkwise_LDshrink(gds_file = gdsf, region_id = region_id, outfile = outf,cutoff=cutoff)

#Rprof(NULL)
