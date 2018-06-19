

library(profvis)
library(EigenH5)
library(LDshrink)
library(tidyverse)
# library(future)
library(progress)
# plan(sequential)
# mb <-profvis({

# setwd("~/Dropbox/PolygenicRSS/code/snakemake_files/")
# load("ssSNP.RData")
# stop()


cutoff <- 1e-3
input_file <- snakemake@input[["input_file"]]
mapf <- snakemake@input[["mapf"]]

subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subldf"]]
bdf <- snakemake@input[["bdf"]]
output_file <- snakemake@output[["evdf"]]
useLDshrink <- snakemake@params[["useLDshrink"]]=="T"
#
# if(useLDshrink){
#   # save
  # load("evd_t.RData")
#   stop()
# }
## cutoff <- as.numeric(snakemake@params[["cutoff"]])
## if(length(cutoff)==0){
cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]
m <- formals(LDshrink::LDshrink)[["m"]]
Ne <- formals(LDshrink::LDshrink)[["Ne"]]
## }

stopifnot( !is.null(input_file), !is.null(output_file),!is.null(mapf),!is.null(bdf))

stopifnot(file.exists(input_file), !file.exists(output_file),file.exists(mapf),file.exists(bdf))
break_df <- read_delim(bdf,delim="\t") #%>% group_by(chr) %>% mutate(start=ifelse(start==min(start),0,start),
                                       #                             stop=ifelse(stop==max(stop),max(stop)*2,stop)) %>% ungroup()


ind_v <- readRDS(subldf)
snp_df <- read_delim(subsnpf,delim="\t")
op <- nrow(snp_df)
tsnp_df <- snp_df
snp_df <-assign_snp_block(snp_df,break_df,assign_all = T)
semi_join(break_df,snp_df) %>% distinct(chr)
stopifnot(nrow(snp_df)==op)
p <- nrow(snp_df)

stopifnot(sorted_snp_df(snp_df))

map_df <- read_df_h5(mapf,"SNPinfo")
cat("Assigning Map\n")
snp_df <- assign_map(snp_df,map_df)
cat("Removing Map\n")
rm(map_df)
stopifnot(!is.unsorted(snp_df$chr))
stopifnot(!is.unsorted(snp_df$snp_id,strictly = T))

stopifnot(group_by(snp_df,chr) %>%
            summarise(sorted=!is.unsorted(map)) %>%
            summarise(sorted=all(sorted)) %>%
            pull(1))

stopifnot(file.exists(input_file),
          !file.exists(output_file),
          !is.null(snp_df[["region_id"]]),
          !is.null(snp_df[["map"]]))




p <- nrow(snp_df)
dosage_dims <-EigenH5::dim_h5(input_file,"dosage")
tch <- EigenH5::dim_h5(input_file,"SNPinfo/chr")


SNPfirst <-  dosage_dims[1]==tch
N <- length(ind_v)
if(!SNPfirst){
  stopifnot(dosage_dims[2]==tch)
}

snp_df <- dplyr::mutate(snp_df,ld_snp_id=snp_id)
write_df_h5(snp_df,"LDinfo",output_file)
pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""])
write_df_h5(pl,groupname = "Wildcards",filename=output_file)
snp_dfl <- split(snp_df,snp_df$region_id)
cat("Estimating LD")
num_b <- length(snp_dfl)
pb <- progress::progress_bar$new(total = num_b)
for(i in 1:num_b){
  tdf <- snp_dfl[[i]]
  if(SNPfirst){
    dosage <- EigenH5::read_matrix_h5(input_file,"/","dosage",subset_rows=tdf$ld_snp_id,subset_cols=ind_v,doTranspose=T)
  }else{
    dosage <- EigenH5::read_matrix_h5(input_file,"/","dosage",subset_cols=tdf$ld_snp_id,subset_rows=ind_v)
  }
  mrid <- unique(tdf$region_id)
  stopifnot(length(mrid)==1)
  R <- cor(dosage)
  ld_df <- ld2df(R,tdf$SNP)
  EigenH5::write_df_h5(ld_df,paste0("LD_DF/",mrid),output_file)
  pb$tick()
}
#   },packages=c("EigenH5","LDshrink")))
# },m=m,Ne=Ne,cutoff=cutoff,useLDshrink=useLDshrink,SNPfirst=SNPfirst,ind_l=ind_v)


# mmret <- SeqSupport::waitr(retl)
# mret <- purrr::map(retl,future::value)
# })
