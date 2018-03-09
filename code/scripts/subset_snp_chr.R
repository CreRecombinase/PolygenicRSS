library(tidyverse)
library(EigenH5)
inf <- snakemake@input[["gdsf"]]
outf <- snakemake@output[["outf"]]
chromosome <- snakemake@params[["chrom"]]
cat("chromosome:",chromosome)
if(grepl("-",chromosome)){
    chr_bounds <- as.integer(strsplit(chromosome,"-")[[1]])
    tchromosome <- chr_bounds[1]:chr_bounds[2]
}else{
    if(grepl(",",chromosome)){
        tchromosome <- as.integer(strsplit(chromosome,",")[[1]])
    }else{
        tchromosome <- as.integer(chromosome)
    }
}
stopifnot(all(!is.na(tchromosome)))

read_df_h5(inf,"SNPinfo") %>% filter(chr%in%tchromosome) %>% write_delim(outf,delim="\t")
