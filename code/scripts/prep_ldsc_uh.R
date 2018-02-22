library(tidyverse)
library(EigenH5)
# library(SeqSupport)
# library(readr)
# library(purrr)
# library(progress)
rss_rdsf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/RSSp_genome_gwas_uh/sscombined_wtcc_NoConfoundSmall_sim.h5"
tparam_df <- read_df_h5(rss_rdsf,"SimulationInfo")
outf <- paste0("/home/nwknoblauch/Desktop/scratch/polyg_scratch/ldsc_sim_gwas_genome/sscombined_sim_",tparam_df$fgeneid,"_uh_NoConfoundSmall.txt")
mfgeneid <- as.character(tparam_df$fgeneid)
names(outf) <- mfgeneid
rss_rdsf <- snakemake@input[["rdsf"]]
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
outf <- snakemake@output[["ldscf"]]
names(outf) <- mfgeneid



tparam_df <- read_df_h5(rss_rdsf,"SimulationInfo")
nfgeneid <- tparam_df$fgeneid
#tparam_df <- mutate(tparam_df,ldscf=outf[as.character(fgeneid)])

groupl <- as.integer(get_objs_h5(rss_rdsf))
groupl <- sort(groupl[!is.na(groupl)])

uh_df <- map_dfr(as.character(groupl),~as_data_frame(read_matrix_h5(rss_rdsf,.x,dataname="uh"))) %>%
  magrittr::set_colnames(as.character(tparam_df$fgeneid)) %>%
  mutate(snp_id=1:n()) %>%
  gather("fgeneid","uh",-snp_id)



N <- unique(tparam_df$n)
## snp_id <- read_vec(rss_rdsf,"snp_id")
## uh_l <- array_branch(read_2d_mat_h5(rss_rdsf,"/","uh"),margin=2)




## seqSetFilter(gds,variant.id=snp_id)
snpinfo <- read_df_h5(rss_rdsf,"SNPinfo") %>%
  select(SNP,allele) %>%
  separate(allele,c("A1","A2")) %>% mutate(N=N,SNP=paste0("rs",SNP))


split(uh_df,uh_df$fgeneid) %>% walk(function(x){
  mutate(snpinfo,Z=x$uh) %>% write_delim(path =outf[x$fgeneid[1]],delim="\t")
})
# group_by(uh_df,fgeneid) %>% do()
# walk2(uh_l,transpose(tparam_df),
#       function(Zd,tp_dfl,snpinfo,Nd){
#           cat(tp_dfl$ldscf,"\n")
#           mutate(snpinfo,Z=Zd,N=Nd) %>%
#               write_delim(path=tp_dfl$ldscf,delim="\t")
#           },snpinfo=snpinfo,Nd=N)
#Rprof(NULL)
