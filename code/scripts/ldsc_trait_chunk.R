#Rprof(filename=snakemake@output[["proff"]],append=F)

library(SeqSupport)
library(tidyverse)
library(SeqArray)
library(progress)
gdsf <- snakemake@input[["gdsf"]]
evdf <- snakemake@input[["evdf"]]
outf <- snakemake@output[["outf"]]
soutf <- snakemake@output[["soutf"]]


file.create(outf)
file.create(soutf)

tdf <-data_frame(CHR=character(),
                 SNP=character(),
                 BP=integer(),
                 CM=numeric(),MAF=numeric(),L2=numeric())
walk(outf,write_delim,x=tdf,delim="\t")
walk(soutf,write,x=0L)

LDchunk <- as.character(snakemake@params[["LDchunk"]])
pb <- progress_bar$new(total=length(LDchunk))

ldsc_read <-function(fn,LDc){
  pb$tick()
return(read_df_h5(fn,groupname="LDinfo") %>% mutate(region_id=as.integer(LDc),L2=colSums(read_2d_mat_h5(fn,"LD","R")^2)-1) %>% rename(CHR=chr,BP=pos))
}

evd_info <- map2_dfr(evdf,LDchunk,ldsc_read)

gds <- seqOpen(gdsf)
pb <- progress_bar$new(total=length(outf))
snpi <- read_SNPinfo_ldsc_ld(gds) %>% inner_join(evd_info) %>% mutate(out_file=outf[as.integer(CHR)],sout_file=soutf[as.integer(CHR)])
snpil <- split(snpi,snpi$CHR) %>% walk(function(df){
  pb$tick()
  ldscoref <- df$out_file[1]
  countf <- df$sout_file[1]
  readr::write_delim(select(df,CHR,SNP,BP,CM,MAF,L2),path=ldscoref,delim="\t")
  nc <- dplyr::filter(df,MAF>0.05) %>% nrow()
  write(x=nc,file=countf)
  stopifnot(all(file.exists(c(ldscoref,countf))))
})


#Rprof(NULL)
