# save.image()
# stop()


library(EigenH5)
library(LDshrink)
library(tidyverse)
library(future)
library(progress)
plan(sequential)

cutoff <- 1e-3
input_file <- snakemake@input[["input_file"]]
mapf <- snakemake@input[["mapf"]]

subsnpf <- snakemake@input[["subsnpf"]]
bdf <- snakemake@input[["bdf"]]
output_file <- snakemake@output[["evdf"]]
useLDshrink <- snakemake@params[["useLDshrink"]]=="T"
## cutoff <- as.numeric(snakemake@params[["cutoff"]])
## if(length(cutoff)==0){
cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]
m <- formals(LDshrink::LDshrink)[["m"]]
Ne <- formals(LDshrink::LDshrink)[["Ne"]]
## }

stopifnot( !is.null(input_file), !is.null(output_file),!is.null(mapf),!is.null(bdf))

stopifnot(file.exists(input_file), !file.exists(output_file),file.exists(mapf),file.exists(bdf))
break_df <- read_delim(bdf,delim="\t") %>% group_by(chr) %>% mutate(start=ifelse(start==min(start),0,start),
                                                                    stop=ifelse(stop==max(stop),max(stop)*2,stop)) %>% ungroup()



snp_df <- read_delim(subsnpf,delim="\t")
snp_df <-assign_snp_block(snp_df,break_df,assign_all = T)
p <- nrow(snp_df)

stopifnot(sorted_snp_df(snp_df))

map_df <- read_df_h5(mapf,"SNPinfo")
cat("Assigning Map\n")
snp_df <- assign_map(snp_df,map_df)
cat("Removing Map\n")
rm(map_df)
stopifnot(group_by(snp_df,chr) %>% summarise(sorted=!is.unsorted(map)) %>% summarise(sorted=all(sorted)) %>% pull(1))

## LDshrink::chunkwise_LD_h5(input_file = input_file,
##                                 output_file = output_file,
##                                 snp_df = snp_df,
##                                 useLDshrink=useLDshrink)


stopifnot(file.exists(input_file),
          !file.exists(output_file),
          !is.null(snp_df[["region_id"]]),
          !is.null(snp_df[["map"]]))

p <- nrow(snp_df)
dosage_dims <-EigenH5::dim_h5(input_file,"dosage")
tch <- EigenH5::dim_h5(input_file,"SNPinfo/chr")


SNPfirst <-  dosage_dims[1]==tch
if(!SNPfirst){
  N <- dosage_dims[1]
  stopifnot(dosage_dims[2]==tch)
}else{
  N <- dosage_dims[2]
}

snp_df <- dplyr::mutate(snp_df,ld_snp_id=snp_id)
write_df_h5(snp_df,"LDinfo",output_file)
snp_dfl <- split(snp_df,snp_df$region_id)

cat("Creating LD futures\n")

retl <- snp_dfl %>% purrr::map(function(tdf,m,Ne,cutoff,useLDshrink,SNPfirst){
  return(future::future({
    if(SNPfirst){
      dosage <- t(EigenH5::read_matrix_h5(input_file,"/","dosage",subset_rows=tdf$ld_snp_id))
    }else{
      dosage <- EigenH5::read_matrix_h5(input_file,"/","dosage",subset_cols=tdf$ld_snp_id)
    }
    mrid <- unique(tdf$region_id)
    stopifnot(length(mrid)==1)
    retl <- LDshrink_evd(dosage,
                         tdf$map,
                         m,
                         Ne,
                         cutoff,
                         useLDshrink=useLDshrink,
                         na.rm=T)
    EigenH5::write_vector_h5(output_file,paste0("EVD/",mrid),"D",retl$D)
    EigenH5::write_matrix_h5(output_file,paste0("EVD/",mrid),"Q",retl$Q)
    EigenH5::write_vector_h5(output_file,paste0("L2/",mrid),"L2",retl$L2)
    EigenH5::write_matrix_h5(output_file,paste0("LD/",mrid),"R",retl$R)

    EigenH5::write_vector_h5(output_file,paste0("LDi/",mrid),"chr",tdf$chr)
    EigenH5::write_vector_h5(output_file,paste0("LDi/",mrid),"pos",tdf$pos)
    EigenH5::write_vector_h5(output_file,paste0("LDi/",mrid),"allele",tdf$allele)
    # EigenH5::write_vector_h5(output_file,paste0("LDi/",mrid),"allele",tdf$snp_id)


    TRUE
  },packages=c("EigenH5","LDshrink")))
},m=m,Ne=Ne,cutoff=cutoff,useLDshrink=useLDshrink,SNPfirst=SNPfirst)
cat("Estimating LD")

mmret <- SeqSupport::waitr(retl)
mret <- purrr::map(retl,future::value)
