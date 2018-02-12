library(SeqSupport)
library(tidyverse)
library(EigenH5)

inf <- snakemake@input[["input_h5"]]
outf <- snakemake@output[["output_h5"]]
gn <-snakemake@params[["gn"]]


param_df <- map_dfc(col_names,~set_names(data_frame(c(read_matrix_h5(tbdf,groupname = "/",dataname = .x))),.x))
snp_info_f <- function(file_name,grp_name,col_names=c("chr","labels","major","minor","pos")){
  read_df_h5(file_name,grp_name,subcols=col_names,filtervec=list(NULL,NULL)) %>%
    mutate(mat_snp_id=1:n(),grp=grp_name)
}
gwas_names <- c("bd","cad","cd","ht","ra","t1d","t2d")
stopifnot(any(gn %in% gwas_names))

snp_info <- map(gwas_names,snp_info_f,file_name=inf) %>%
    map(function(x){
 mutate(x,labels=paste0("rs",as.integer(labels)),major=R.oo::intToChar(major),minor=R.oo::intToChar(minor),pos=as.integer(pos),chr=as.character(chr)) %>%
    rename(SNP=labels)
}) %>% set_names(gwas_names)
yl <-  map(gwas_names,snp_info_f,file_name=inf,col_names="y") %>% set_names(gwas_names)
sum(map_dbl(yl,nrow))

inter_snpsl<- map(snp_info,semi_join,y=si_map,by=c("SNP","chr"))

inter_snpsl <- map(inter_snpsl,function(x){select(x,-mat_snp_id,-grp)})
inter_snps <- reduce(inter_snpsl,inner_join)

sub_snp_info <- map(snp_info,semi_join,y=inter_snps) %>% set_names(gwas_names)
map_si <- select(si_map,SNP,chr,map,region_id)

input_file <- inf
output_file <- outf
df <- sub_snp_info[[gn]]

grp_name <- unique(df$grp)
stopifnot(length(grp_name)==1)
subset_id <- df$mat_snp_id
df <- inner_join(df,map_si)
cat("Reading rows\n")
sub_X <- EigenH5::read_mat_rows_h5(input_file,grp_name,"X",subset_id)
gc()
N <- ncol(sub_X)
cat("Writing dosage\n")
EigenH5::write_mat_h5(output_file,
                      groupname="/",
                      dataname="dosage",
                      data = sub_X,
                      doTranspose = F,chunksizes=as.integer(c(250,N)))
write_df_h5(df,groupname="SNPinfo",output_file)

