library(SeqArray)

input <- snakemake@input[["gdsf"]]


AF <- as.numeric(snakemake@params[["AF"]])
output <- snakemake@output[["temp_gds"]]


in_gds <- seqOpen(input)
data_af <- seqAlleleFreq(in_gds, .progress=TRUE)

data_selection <- dplyr::between(data_af,AF,1-AF)
seqSetFilter(in_gds,variant.sel=data_selection)
seqExport(in_gds,output)
