library(tidyverse)
library(SeqArray)
library(EigenH5)
library(LDshrink)
library(progress)
# data("break_df")
# evdf <- "/run/media/nwknoblauch/Data/EVD_H5/chr19_bd_bd_F_omni_T_F.h5"
#
# outf <- paste0("/home/nwknoblauch/Desktop/scratch/polyg_scratch/eur_w_ld_chr_sscombined/",1:22,".l2.ldscore.gz")
# soutf <- paste0("/home/nwknoblauch/Desktop/scratch/polyg_scratch/eur_w_ld_chr_sscombined/",1:22,".l2.M_5_50")
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

snp_df <- read_df_h5(evdf,"LDinfo") %>% rename(CHR=chr,BP=pos,CM=map)

#select(break_df,chr,region_id) %>%
  # inner_join(read_df_h5(evdf,"LDinfo")) %>%
  #    %>% select(-snp_id)

if(is.null(snp_df[["MAF"]])){
    stopifnot(!is.null(snp_df[["AF"]]))
    snp_df <- mutate(snp_df,MAF=AF)
}
snp_df <- group_by(snp_df,region_id) %>% do(mutate(.,L2=read_vector_h5(evdf,paste0("L2/",as.character(.$region_id[1]),"/L2")))) %>% ungroup() %>% select(-region_id)


pb <- progress_bar$new(total=length(outf))
split(snp_df,snp_df$CHR) %>% walk(function(df){
  pb$tick()
  ldscoref <- outf[df$CHR[1]]
  countf <- soutf[df$CHR[1]]
  mutate(df,SNP=paste0("rs",SNP)) %>% select(CHR,SNP,BP,CM,MAF,L2) %>% readr::write_delim(path=ldscoref,delim="\t")
  nc <- dplyr::filter(df,MAF>0.05) %>% nrow()
  write(x=nc,file=countf)
  stopifnot(all(file.exists(c(ldscoref,countf))))
})


#Rprof(NULL)
