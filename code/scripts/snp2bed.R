library(SeqArray)
library(SNPRelate)


inf <- snakemake@input[["gdsf"]]
out_pref <- trimws(snakemake@params[["pref"]])
tempf <- tempfile()

gds <- seqOpen(inf)

seqGDS2SNP(gds,tempf)
seqClose(gds)

gds <- snpgdsOpen(tempf)

snpgdsGDS2BED(gds,bed.fn=out_pref)
snpgdsClose(gds)
file.remove(tempf)

