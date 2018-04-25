library(tidyverse)
inf <- snakemake@input[["ldetectf"]]
allchr <- snakemake@params[["allchr"]]=="T"
cat("HEY!",snakemake@params[["allchr"]],"\n")
chrom <- as.integer(snakemake@params[["chrom"]])
if(length(chrom)==0){
    chrom <- 1:22
}
outf <- snakemake@output[["bdf"]]
stopifnot(!is.null(inf),!is.null(outf))

break_df <- read_delim(inf,delim="\t",trim_ws = T) %>% mutate(chr=as.integer(gsub("chr","",chr))) %>% filter(chr %in% chrom) %>% mutate(region_id=1:n())
if(!allchr){
    break_df %>% write_delim(path=outf,delim="\t")
}else{
    break_df %>% group_by(chr) %>% summarise(start=min(start),stop=max(stop)) %>% mutate(region_id=1:n()) %>% write_delim(path=outf,delim="\t")
}