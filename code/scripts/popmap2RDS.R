library(tidyverse)
library(LDshrink)
library(EigenH5)


                                        #pop <- "CEU"
## save.image()
## stop()

outf <- snakemake@output[["mapf"]]
pop  <- snakemake@params[["pop"]]
if(pop=="EUR"){
    pop <- "CEU"
}

destination <- tempfile()
base_url <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/"
download.file(paste0(base_url, pop, "_omni_recombination_20130507.tar"), destfile = destination)
destination_dir  <- tempdir()
untar(destination, exdir=destination_dir)

## result_files <- , full.names = T)



map_file_df <- tibble(filename=dir(file.path(destination_dir, pop), full.names=T)) %>%
    mutate(chrom=as.integer(gsub(paste0(pop, "-([0-9]+)-final.txt.gz"), "\\1", basename(filename))))


rff <- function(filename, chrom){
    read_delim(filename, delim="\t", col_names=c("pos", "rate", "map", "filtered"), skip = 1, trim_ws = T) %>% mutate(chr=chrom) %>% return()}
map_df  <- group_by(map_file_df,filename) %>% do(rff(.$filename,.$chrom)) %>% ungroup() %>% select(chr,map,pos)  %>% arrange(chr,pos)

## %>% unnest()
## map_df <- pmap_dfr(map_file_df, ) %>% arrange(chr, pos) %>% select(chr, map, pos)
## data(map_parameters)
saveRDS(map_df, outf)
