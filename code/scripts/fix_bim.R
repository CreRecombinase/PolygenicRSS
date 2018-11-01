library(readr)
library(dplyr)

inf <- snakemake@input[["inf"]]
chrom <- snakemake@params[["chrom"]]
outf <- snakemake@output[["outf"]]


stopifnot(file.exists(inf),!is.null(inf),!is.null(chrom))

read_delim(inf,delim="\t",col_names=c("ch","id","map","pos","ref","alt")) %>% mutate(ch=as.integer(chrom),
                                                                                     pos=as.integer(gsub("[0-9]+:([0-9]+).*","\\1",id)))  %>% write_delim(outf,col_names=F)
