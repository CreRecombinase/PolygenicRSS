library(tidyverse)

ldetect <- snakemake@params[["use_ldetect"]]
use_ldetect <- ldetect=="T"

inf  <- snakemake@input[["inf"]]
if(!file.exists(inf)){
    inf <- "./fourier_ls-all.bed"
}
stopifnot(file.exists(inf))
cat("inf:",inf,"\n")
chrom <- as.integer(snakemake@params[["chrom"]] %||% 1:22)
outf <- snakemake@output[["bdf"]]
stopifnot(!is.null(inf),!is.null(outf))
if(!is.na(as.numeric(ldetect))){
    file.copy(inf,outf)
}else{
    break_df <- read_delim(inf,delim="\t",trim_ws = T) %>% mutate(chr=as.integer(gsub("chr","",chr))) %>% filter(chr %in% chrom) %>% mutate(region_id=1:n())
    if(use_ldetect){
        break_df %>% write_delim(path=outf,delim="\t")
    }else{
        break_df %>% group_by(chr) %>% summarise(start=0,stop=max(stop)+10) %>% mutate(region_id=1:n())  %>%  write_delim(path=outf,delim="\t")
    }
}
