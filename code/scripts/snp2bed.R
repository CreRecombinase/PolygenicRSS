library(SeqArray)
library(SNPRelate)
library(SeqSupport)
library(dplyr)
library(readr)

inf <- snakemake@input[["gdsf"]]
bim_outf <- snakemake@output[["bimf"]]
out_pref <- trimws(snakemake@params[["pref"]])

tempf <- tempfile()

gds <- seqOpen(inf)
snp_info <- read_SNPinfo_gds(gds)
seqGDS2SNP(gds,tempf)
seqClose(gds)

gds <- snpgdsOpen(tempf)

                          
snpgdsGDS2BED(gds,bed.fn=out_pref)
snpgdsClose(gds)
file.remove(tempf)
bimd <- readr::read_delim(bim_outf,delim="\t",col_names=c("chr",
                                                          "temp",
                                                          "map",
                                                          "pos",
                                                          "ref",
                                                          "alt"),
                          col_types=cols(chr="c",
                                         temp="c",
                                         map="i",
                                         pos="i",
                                         ref="c",
                                         alt="c"))
p <- nrow(bimd)
nbimd <- left_join(bimd,snp_info) %>% select(chr,SNP,map,pos,ref,alt)
stopifnot(all(!is.na(nbimd$SNP)),
          nrow(nbimd)==p,all.equal(bimd$pos,nbimd$pos))
write_delim(nbimd,bim_outf,delim="\t",col_names=F)

