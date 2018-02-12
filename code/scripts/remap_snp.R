#library(SeqSupport)
# save.image()
# stop()
library(tidyverse)
library(EigenH5)
library(LDshrink)
data("break_df")

inf <- snakemake@input[["input_h5"]]
snp_dff <- snakemake@input[["snp_df"]]
mapf <- snakemake@input[["mapf"]]
marchf <- snakemake@input[["marchf"]]

outf <- snakemake@output[["outf"]]
gn <-snakemake@params[["gn"]]


stopifnot(!is.null(inf),!is.null(snp_dff),!is.null(mapf),!is.null(marchf),!is.null(outf),!is.null(gn))

gwas_names <- c("bd","cad","cd","ht","ra","t1d","t2d")
stopifnot(any(gn %in% gwas_names))
map_dat <- read_df_h5(mapf,"SNPinfo")


merge_df <- read_df_h5(marchf,"mergeInfo")


col_names <- c("labels","major","minor")
param_df <- map_dfc(col_names,~set_names(data_frame(c(read_matrix_h5(inf,groupname = "/",dataname = .x))),.x)) %>%
  mutate(labels=abs(as.integer(labels)),
         major=R.oo::intToChar(major),
         minor=R.oo::intToChar(minor),
         allele=paste0(major,",",minor)) %>% rename(SNP=labels) %>% select(-major,-minor) %>% mutate(snp_id=1:n())
snp_df <- read_df_h5(snp_dff,"SNPinfo")
param_df <-left_join(param_df,merge_df,by=c("SNP")) %>% mutate(SNP=ifelse(is.na(rsCurrent),SNP,rsCurrent)) %>% select(-rsCurrent)
param_df <- left_join(param_df,snp_df,by=c("SNP"))
rm(snp_df)
gc()
param_df <- filter(param_df,!is.na(pos)) %>% arrange(chr,pos)
ld_r <- LDshrink::set_ld_region(ld_regions = break_df,snp_info = param_df)
param_df <- mutate(param_df,region_id=ld_r)
param_df <- filter(param_df,region_id!=-1)


map_l <- split(map_dat,map_dat$chr)
param_l <- split(param_df,param_df$chr)
param_df <- map2_df(map_l,param_l,~mutate(.y,map=LDshrink::interpolate_map(.x$map,.x$pos,.y$pos)))
rm(map_dat,param_l)
param_df <- distinct(param_df,snp_id,.keep_all=T) %>% arrange(chr,pos)
gc()
sub_X <-EigenH5::read_matrix_h5(inf,groupname = "/",dataname = "X",subset_rows = param_df$snp_id)
param_df <- mutate(param_df,snp_id=1:n());
N <- ncol(sub_X)
cat("Writing dosage\n")
EigenH5::write_matrix_h5(outf,
                      groupname="/",
                      dataname="dosage",
                      data = sub_X,
                      doTranspose = F,
                      chunksizes=as.integer(c(250,N)))

EigenH5::write_df_h5(param_df,groupname="SNPinfo",outf)

