#Rprof(filename=snakemake@output[["proff"]],append=F)

library(SeqSupport)
library(tidyverse)
library(SeqArray)
library(EigenH5)
library(progress)
#gdsf <- snakemake@input[["gdsf"]]
evdf <- snakemake@input[["evdf"]]
outf <- snakemake@output[["outf"]]
soutf <- snakemake@output[["soutf"]]


file.create(outf)
file.create(soutf)

tdf <-data_frame(CHR=character(),
                 SNP=character(),
                 BP=integer(),
                 CM=numeric(),
                 MAF=numeric(),
                 L2=numeric())
walk(outf,write_delim,x=tdf,delim="\t")
walk(soutf,write,x=0L)

snp_df <- read_df_h5(evdf,"LDinfo") %>% rename(CHR=chr,BP=pos)
snp_df <- group_by(snp_df,region_id) %>% do(mutate(.,L2=read_vector_h5(evdf,paste0("LDSC/",as.character(.$region_id[1])),"L2"))) %>% ungro

#pb <- progress_bar$new(total=length(LDchunk))

## ldsc_read <-function(fn,LDc){
##   pb$tick()
## return(read_df_h5(fn,groupname="LDinfo") %>% mutate(region_id=as.integer(LDc),L2=colSums(read_2d_mat_h5(fn,"LD","R")^2)-1) %>% rename(CHR=chr,BP=pos))
## }

## evd_info <- map2_dfr(evdf,LDchunk,ldsc_read)




## gds <- seqOpen(gdsf)
## pb <- progress_bar$new(total=length(outf))
## snpi <- read_SNPinfo_ldsc_ld(gds) %>% inner_join(evd_info) %>% mutate(out_file=outf[as.integer(CHR)],sout_file=soutf[as.integer(CHR)])
split(snp_df,snpi$CHR) %>% walk(function(df){
  pb$tick()
  ldscoref <- outf[df$CHR[1]]
  countf <- soutf[df$CHR[1]]
  readr::write_delim(select(df,CHR,SNP,BP,CM,MAF,L2),path=ldscoref,delim="\t")
  nc <- dplyr::filter(df,MAF>0.05) %>% nrow()
  write(x=nc,file=countf)
  stopifnot(all(file.exists(c(ldscoref,countf))))
})


#Rprof(NULL)
