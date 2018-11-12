

library(tidyverse)
library(ldshrink)
library(glue)


inputf <- snakemake@input[["inputf"]]
output_file <-snakemake@output[["outf_tf"]]
chunks <- as.integer(snakemake@params[["chunks"]])
stopifnot(length(chunks)==length(output_file),all(!is.na(chunks)))
                                        #chunksize <- as.integer(snakemake@params[["chunk_size"]])
num_chunks <- length(chunks)

useLDetect <- snakemake@params[["useLDetect"]]
if(useLDetect %in% c("T","F")){
    useLDetect <- useLDetect=="T"
    chunksize <- NA_integer_
}else{
    useLDetect <- F
    chunksize <- as.integer(snakemake@params[["useLDetect"]])
}


cat("inputf:",inputf,"\n")

if(useLDetect==F){
    subsnp_df <- fst::read_fst(inputf) %>% dplyr::rename(chr=chrom)  %>%
        chunk_genome(chunk_size=chunksize) %>%
        dplyr::rename(chrom=chr)
}else{
    break_df <- readr::read_delim(snakemake@input[["ldetectf"]],delim="\t",trim_ws=T) %>%
        dplyr::mutate(region_id=1:n())  %>%
        dplyr::mutate(chr=as.integer(gsub("chr","",chr)),start=as.integer(start),stop=as.integer(stop)) %>%
        group_by(chr) %>%
        mutate(start=dplyr::if_else(start==min(start),0L,start),stop=if_else(stop==max(stop),as.integer(max(stop)*10L),stop)) %>%
        ungroup()
    subsnp_df <- fst::read_fst(inputf) %>%
        dplyr::rename(chr=chrom) %>%
        ldshrink::assign_snp_block(break_df=break_df,assign_all=F) %>% filter(!is.na(region_id)) %>%
        dplyr::rename(chrom=chr)
}

cat("Creating chunks\n")
chunk_df <- distinct(subsnp_df,region_id) %>% mutate(chunk_reg=as.integer(ggplot2::cut_number(as.integer(region_id),num_chunks,labels=F)))

filter(chunk_df,chunk_reg==chunks[1])
                                        #saveRDS(t_res,"ws.RDS")
saveret <- function(x,out_file){

    ## mchunk <-    formatC(x$region_id[1],format = "d",digits = 1,flag = "0")
    ## out_file <-    paste0(output_pref,".",mchunk,".chunked.fst")
    cat("Output_file is :",out_file,"\n")
    cat("nrow is :",nrow(x),"\n")
    feather::write_feather(x,out_file)
}
cat("chunking and writing\n")
for(i in chunks){
    outf <- output_file[as.integer(i)]
    semi_join(subsnp_df,filter(chunk_df,chunk_reg==as.integer(i))) %>% saveret(out_file=outf)
}
## subsnp_df  %>%    group_by(region_id) %>%
##     do(saveret(.))
