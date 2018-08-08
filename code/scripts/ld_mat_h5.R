library(EigenH5)
library(LDshrink)
library(tidyverse)
library(magrittr)
# library(SeqArray)
library(progress)

cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]
m <- formals(LDshrink::LDshrink)[["m"]]
Ne <- formals(LDshrink::LDshrink)[["Ne"]]


## input_file <- "/home/nwknoblauch/Dropbox/scratch/polyg_scratch/h5/ALL_bd_geno.h5"
## output_file <- "/home/nwknoblauch/Dropbox/scratch/polyg_scratch/EVD_H5/sparse_LD/bd_19.h5"
## file.remove(output_file)
## snp_dff <- "/home/nwknoblauch/Dropbox/PolygenicRSS/data/Snakemake_inputs/chr1-22AF0.01SNP0N0_cd_bd_LD.txt.gz"
## mapf <-"/run/media/nwknoblauch/Data/1kg/1000-genomes-genetic-maps/CEU_map.h5"
input_file <- snakemake@input[["input_file"]]
subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subldf"]]
bdf <- snakemake@input[["bdf"]]
output_file <- snakemake@output[["evdf"]]
useLDshrink <- snakemake@params[["useLDshrink"]]=="T"

ind_v <- readRDS(subldf)
#snp_df <- read_delim(subsnpf,delim="\t")
#op <- nrow(snp_df)
#tsnp_df <- snp_df
#snp_df <-assign_snp_block(snp_df,break_df,assign_all = T)
#semi_join(break_df,snp_df) %>% distinct(chr)
#stopifnot(nrow(snp_df)==op)
#p <- nrow(snp_df)

map_df <- read_df_h5(mapf, "SNPinfo")
snp_df <- assign_map(readr::read_delim(subsnpf,delim="\t"), map_df)

t_subsnp_df <- group_by(snp_df,chr) %>% do(chunk_genome(., n_chunks=10) %>% mutate(new_index=0:(n()-1))) %>% ungroup() %>% split(.$chr)


chn <- names(t_subsnp_df)

for(ch in chn){
  cat("chrom:",ch,"\n")
  subsnp_df <-t_subsnp_df[[ch]]
  snp_dfl <- split(subsnp_df,subsnp_df$region_id)
  total_p <- nrow(subsnp_df)
  ndfl <-length(snp_dfl)
  num_pb <-((10^2-2)/2)
  pb <- progress_bar$new(total = num_pb)
  i <- 1
  for(i in 1:ndfl){
    subsnp_dfa <- snp_dfl[[i]]
    mchr <-unique(subsnp_dfa$chr)
    dosage_a <- read_matrix_h5(input_file,"dosage",subset_rows=subsnp_dfa$snp_id)
    res <- ld2df(data = dosage_a,
                 mapd = subsnp_dfa$map,
                 rsid=subsnp_dfa$new_index,
                 m=m,Ne = Ne,cutoff = cutoff,r2cutoff = 0,
                 progress = F,useLDshrink = T)
    EigenH5::write_df_h5(df = res,filename = output_file,datapath = paste0("LD_DF/",mchr),max_dims=NA_integer_,append=T)
    pb$tick()
    if(i!=ndfl){
      for(j in (i+1):ndfl){
        subsnp_dfb <- snp_dfl[[j]]

        dosage_b <- read_matrix_h5(input_file,"dosage",subset_rows=subsnp_dfb$snp_id)

        res <- ld2df_p(data_a = dosage_a,data_b = dosage_b,
                       mapd_a = subsnp_dfa$map,mapd_b=subsnp_dfb$map,
                       rsid_a=subsnp_dfa$new_index,rsid_b = subsnp_dfb$new_index,
                       m=m,Ne = Ne,cutoff = cutoff,r2cutoff = 0,
                       progress = F,useLDshrink = T)
        EigenH5::write_df_h5(df = res,filename = output_file,datapath=paste0("LD_DF/",mchr),append=T)
        pb$tick()
      }
    }
  }
}
# big_df <- read_df_h5(output_file,paste0("LD_DF/",mchr))
# filter(big_df,rowsnp!=colsnp) %>% ggplot(aes(x=abs(r)))+geom_histogram(bins=100)+scale_x_log10()
# l2df <- filter(big_df,rowsnp!=colsnp) %>% group_by(rowsnp) %>% summarise(l2=sum(r^2)-1) %>% rename(new_index=rowsnp)
# nsubsnp_df <- inner_join(l2df,subsnp_df)
# anti_join(subsnp_df,l2df)
#
# ggplot()


