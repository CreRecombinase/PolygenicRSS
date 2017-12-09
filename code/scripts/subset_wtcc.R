library(SeqSupport)
library(tidyverse)
inf <- "/media/nwknoblauch/Data/combined.h5"
rhdf5::h5ls(inf)

gdsf <- "~/Desktop/scratch/polyg_scratch/all_gds/ALL_seq_hapmap_geno.gds"
gds <- SeqArray::seqOpen(gdsf)
si_map <- read_SNPinfo_gds(gds,alleles=T,map=T) %>% separate(allele,into=c("minor","major"))


snp_info_f <- function(file_name,grp_name,col_names=c("chr","labels","major","minor","pos")){
  read_df_h5(file_name,grp_name,subcols=col_names,filtervec=list(NULL,NULL)) %>%
    mutate(mat_snp_id=1:n(),grp=grp_name)
}
gwas_names <- c("bd","cad","cd","ht","ra","t1d","t2d")

snp_info <- map(gwas_names,snp_info_f,file_name=inf) %>% map(function(x){
  mutate(x,labels=paste0("rs",as.integer(labels)),major=R.oo::intToChar(major),minor=R.oo::intToChar(minor),pos=as.integer(pos),chr=as.character(chr)) %>%
    rename(SNP=labels)
}) %>% set_names(gwas_names)
yl <-  map(gwas_names,snp_info_f,file_name=inf,col_names="y") %>% set_names(gwas_names)
sum(map_dbl(yl,nrow))


inter_snpsl<- map(snp_info,semi_join,y=si_map,by=c("SNP","chr"))

inter_snpsl <- map(inter_snpsl,function(x){select(x,-mat_snp_id,-grp)})
inter_snps <- reduce(inter_snpsl,inner_join)

sub_snp_info <- map(snp_info,semi_join,y=inter_snps)

noutf <- "/media/nwknoblauch/Data/combined_harmonized.h5"
map_si <- select(si_map,SNP,chr,map)
subset_write <- function(df,input_file,output_file,map_si){
  grp_name <- unique(df$grp)
  stopifnot(length(grp_name)==1)
  subset_id <- df$mat_snp_id
  df <- inner_join(df,map_si)
  sub_X <- rhdf5::h5read(input_file,paste0(grp_name,"/X"))[,subset_id]
  gc()





}
iwalk(gwas_names)

