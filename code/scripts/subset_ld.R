

library(EigenH5)
library(LDshrink)
library(tidyverse)
library(progress)


input_file <- snakemake@input[["input_file"]]
mapf <- snakemake@input[["mapf"]]

subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subldf"]]
bdf <- snakemake@input[["bdf"]]
output_file <- snakemake@output[["output_f"]]
useLDetect <-(snakemake@params[["useLDetect"]] %||% "T")=="T"
chr_df <- data_frame(chr = as.integer(snakemake@params[["chrom"]] %||% 1:22))
chunk_tot <- as.integer(snakemake@params[["chunk_tot"]] %||% 1)
chunk_ind <- as.integer(snakemake@params[["chunk_a"]] %||% chunk_tot)

stopifnot(!is.null(input_file),
          !is.null(output_file),
          !is.null(mapf),
          !is.null(bdf),
          file.exists(input_file),
          !file.exists(output_file),
          file.exists(mapf),
          file.exists(bdf))

break_df <- read_delim(bdf,delim="\t")  %>% mutate(chr=as.integer(chr),
                                                   start=as.integer(start),
                                                   stop=as.integer(stop),
                                                   region_id=as.integer(region_id))  %>% semi_join(chr_df)

if(!is.null(subsnpf)){
    if(tools::file_ext(subsnpf)=="h5"){
        snp_df <- read_df_h5(subsnpf,"SNPinfo")
    }else{
        snp_df <- read_delim(subsnpf,delim="\t")
    }
}else{
    snp_df <- read_df_h5(input_file,"SNPinfo")
    all_alleles  <- outer(c("A","C","T","G"),c("A","C","T","G"),function(x,y)paste(x,y,sep=","))
    all_alleles <- data_frame(allele=c(all_alleles[upper.tri(all_alleles)],all_alleles[lower.tri(all_alleles)]))
    snp_df <- semi_join(snp_df,all_alleles)
}
snp_df <- semi_join(snp_df,break_df) %>% group_by(chr) %>% distinct(pos,.keep_all=T)

if(!useLDetect){
    tri_block_fn <- function(tot_size, chunk_ind){
        mat_size <- (1/2)*(sqrt(8*tot_size+1)-1)
        stopifnot(mat_size==round(mat_size))
        tmat <- matrix(0, mat_size, mat_size)
        tmat[upper.tri(tmat, diag=T)] <- 1:tot_size
        stopifnot(all(tmat[upper.tri(tmat, diag=T)]>0))
        trow <- row(tmat)[tmat==chunk_ind]
        tcol <- col(tmat)[tmat==chunk_ind]
        retdf <- data_frame(row=trow, col=tcol, chunk_tot=mat_size, chunk_ind=chunk_ind)
        stopifnot(nrow(retdf)==1)
        return(retdf)
    }
    my_chunks <- map_df(1:chunk_tot, ~tri_block_fn(chunk_tot, .x))
    snp_df <- chunk_genome(snp_df, n_chunks = my_chunks$chunk_tot[1]) %>% mutate(chunk=region_id)
}else{
    my_chunks <- mutate(break_df, chunk=as.integer(gl(n = n_chunks,k = ceiling(n()/n_chunks), length=n())))
    snp_df <-assign_snp_block(snp_df, break_df, assign_all = T)
    semi_join(break_df, snp_df) %>% distinct(chr)
    snp_df <- inner_join(snp_df, break_df)
}

map_df <- read_df_h5(mapf, "SNPinfo")
cat("Assigning Map\n")
snp_df <- assign_map(snp_df, map_df)
cat("Removing Map\n")
rm(map_df)
stopifnot(!is.unsorted(snp_df$chr))
stopifnot(!is.unsorted(snp_df$snp_id, strictly = T))
stopifnot(group_by(snp_df, chr) %>%
          summarise(sorted=!is.unsorted(map)) %>%
          summarise(sorted=all(sorted)) %>%
          pull(1))

stopifnot(file.exists(input_file),
          !file.exists(output_file),
          !is.null(snp_df[["region_id"]]),
          !is.null(snp_df[["map"]]))


snp_dfl <- split(snp_df, snp_df$chunk)
cat("writing chunks")
num_b <- length(snp_dfl)
for(i in 1:num_b){
    tdf <- snp_dfl[[i]]
    mrid <- unique(tdf$chunk)
    write_df_h5(tdf, output_file, paste0(mrid, "/SNPinfo"))
}
