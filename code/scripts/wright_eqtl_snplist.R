library(tidyverse)
library(SeqSupport)
library(SeqArray)
library(EigenH5)
library(progress)
tgdsf <- "~/Dropbox/scratch/polyg_scratch/gds/ALL_EUR_geno.gds"
expinfof <-"/run/media/nwknoblauch/Data/eqtl_summary/gene.probeid.txt.gz"
eqtl_dir <- "/run/media/nwknoblauch/Data/eQTL_SumStats/Wrights_summary_hapmap"
eqtlf <-"/run/media/nwknoblauch/Data/eQTL_SumStats/WrighteQTL_zstd.h5"
beqtlf <- "/run/media/nwknoblauch/Data/eQTL_SumStats/WrighteQTL.h5"

system.time(tm <- read_matrix_h5(eqtlf,"/","uh",subset_cols=c(200L,300L,1000L)))
system.time(bm <- read_matrix_h5(beqtlf,"/","uh",subset_cols=c(200L,300L,1000L)))

snpif <-"/run/media/nwknoblauch/Data/eqtl_summary/allsnps_hapmap_common_positions_gmap.txt.gz"
exp_idf <- read_delim(expinfof,delim="\t") %>% mutate(chr=gsub("chr","",g0),filename=file.path(eqtl_dir,g0,paste0(probeID,".Rdata"))) %>% filter(file.exists(filename)) %>% mutate(exp_id=1:n())
subsnp_df <-read_delim(snpif,delim=" ",col_names = c("chr","SNP","pos","map")) %>% mutate(chr=as.character(chr)) %>% mutate(subset_id=1:n())
sub_p <- nrow(subsnp_df)
# all(1:sub_p==subsnp_df$snporder)
tgds <-seqOpen(tgdsf)
snp_df <- read_SNPinfo_gds(tgds,alleles = T)
bsnp_df <- inner_join(snp_df,subsnp_df,by=c("chr","pos")) %>% arrange(snp_id) %>% select(SNP=SNP.x,snp_id,chr,pos,allele,subset_id)
# basnp_df <- anti_join(subsnp_df,snp_df,by=c("chr","pos"))
# bbasnp_df <- anti_join(subsnp_df,snp_df,by=c("SNP"))

sid <- bsnp_df$subset_id
bsnp_df <- select(bsnp_df,-subset_id)

new_p <-nrow(bsnp_df)
# load(texp_df$filename[2])
# length(tstat.int)

g <-nrow(exp_idf)
create_matrix_h5(eqtlf,groupname = "/",dataname = "uh",data = numeric(),dims=c(new_p,g),chunksizes=c(as.integer(new_p/10),1L),filter_options=c(22L))
# load(exp_idf$filename[1])
# sub_t <- tstat.int[bsnp_df$subset_id]
ret_t <- function(fn,subset_ind){
  load(fn)
  tstat.int <- tstat.int[subset_ind]/100
  attr(tstat.int,"dim") <- c(new_p,1)
  return(tstat.int)
}
sret_t <- safely(ret_t)
rm(snp_df)
rm(subsnp_df)
gc()
error_i <- integer()
#error_i <-c(33375L,38001L)
pb <- progress_bar$new(total=g)
library(furrr)
plan(multicore)
chunk_g <-BBmisc::chunk(1:g,chunk.size = 100)

for(i in 1:length(chunk_g)){
  # load(exp_idf$filename[i])
  tcg <-chunk_g[[i]]
  tres <- future_map(exp_idf$filename[tcg],~sret_t(.x,sid)) %>%transpose()
  tmat <-do.call(cbind,compact(tres$result))
  # tsr <- sret_t(exp_idf$filename[i],sid)
  good_res <- map_lgl(tres$error,~is.null(.x))
  good_ind <- tcg[good_res]
  if(sum(good_res)>0){
    write_matrix_h5(filename =eqtlf,
                    groupname = "/",
                    dataname = "uh",
                    data = tmat,
                    subset_cols=tcg[good_res])
  }else{
    cat("Error!",i,"\n")
    error_i <- c(error_i,tcg[!good_res])
  }
  pb$tick(len = length(tcg))
  if(i %%10==0){
    gc()
  }
}
# error <- c(




